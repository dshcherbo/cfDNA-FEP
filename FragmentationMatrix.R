library(tidyverse)

#### Settings ####

#settings <- readRDS("settings.rds")

settings <- list()

settings$testRange <- c(0, 300) # Region to analyze downstream of each primer
settings$ratioLimits <- c(0, 100, 100, 300) # short/long ratio limits: low short, high short, low long, high long

settings$dir <- paste0("data/", settings$experimentID, "/") # alignment directory
settings$dirFastq <- paste0("fastq/", settings$experimentID, "/") # FASTQ directory
settings$metaDataFile <- "experiments.xlsx" # excel table with sample metadata

settings$allPrimersBed <- c("primer_design/00-all_primers.bed") # BED with coords for all designed primers
settings$targetPrimersList <- c("primer_design/cfDNA_FEP_0.4.tsv") # TSV with name and seq of primers used in the current experiment
settings$nReadsCutOff <- 200

# Loading data
dfList_batch1 <- readRDS(paste0("data/batch1/", "dfList.rds"))  # Fragment end positions in target regions for each sample
dfList_batch2 <- readRDS(paste0("data/batch2/", "dfList.rds"))  # Fragment end positions in target regions for each sample

metaDataTable_batch1 <- openxlsx::read.xlsx(settings$metaDataFile, sheet = "batch1")
metaDataTable_batch2 <- openxlsx::read.xlsx(settings$metaDataFile, sheet = "batch2")

metaDataTable <- bind_rows(batch1 = metaDataTable_batch1, batch2 = metaDataTable_batch2, .id = "batch") %>% as_tibble()

#### Functions ####
auc_dens <- function(coords, lowlim = 0, highlim = 100)
{
  if (length(coords) < 2) return(NA)
  else{
    coords_d <- stats::density(coords, from = lowlim, to = highlim)
    id <- order(coords_d$x)
    return(sum(diff(coords_d$x[id])*zoo::rollmean(coords_d$y[id],2)))
  }
}

find_top_peaks_for_primer <- function(coords, from_pos = 0, to_pos = 300, peaks_threshold = 100, return_peak = "both")
{
  # find top peaks in densities
  if (length(coords) < 2) return(NA)
  
  coords_d <- stats::density(coords, from = from_pos, to = to_pos, n = to_pos - from_pos)
  peak_coords <- which(diff(sign(diff(coords_d$y)))==-2)+1 # finding coordinates of the peaks

  top_peak_val_short <- max(coords_d$y[peak_coords[which(peak_coords<peaks_threshold)]])
  top_peak_val_long <- max(coords_d$y[peak_coords[which(peak_coords>=peaks_threshold)]])
  
  short_peak_coord <- round(coords_d$x[which(coords_d$y == top_peak_val_short)][1])
  long_peak_coord <- round(coords_d$x[which(coords_d$y == top_peak_val_long)][1])
  
  if(return_peak == "both") return(c(short_peak_coord, long_peak_coord))
  else if(return_peak == "short") return(c(short_peak_coord))
  else if(return_peak == "long") return(c(long_peak_coord))
  else stop("Unkonown reuturn_peak value, should be both, short or long")
}

peaks_ratio <- function(coords, from_pos = 0, to_pos = 300, short_peak = 75, long_peak = 150, n = 5)
{
  # compute ratio of the two peaks
  
  if (length(coords) < 2) return(NA)
  if (anyNA(short_peak, long_peak)) return(NA)
  
  coords_d <- stats::density(coords, from = from_pos, to = to_pos, n = to_pos - from_pos)
  
  # finding area under the peak +-n bp from the peak center
  short_peak_start <- short_peak - n
  short_peak_end <- short_peak + n
  
  long_peak_start <- long_peak - n
  long_peak_end <- long_peak + n
  
  if(short_peak_start <= from_pos) short_peak_start <- from_pos
  if(long_peak_end >= to_pos) long_peak_end <- to_pos
  
  auc_short_peak <- MESS::auc(x = coords_d$x, y = coords_d$y, from = short_peak_start, to = short_peak_end)
  auc_long_peak <- MESS::auc(x = coords_d$x, y = coords_d$y, from = long_peak_start, to = long_peak_end)
  
  return( c(auc_short_peak/auc_long_peak) )
}


#### Preparing df_raw from dfList ####

df_raw <- c(dfList_batch1, dfList_batch2) %>%
  bind_rows(.id = "sample") %>%
  left_join(., metaDataTable %>% select(SampleID, diagnosis = Diagnosis, stage = Stage), by = c("sample" = "SampleID")) %>% # Adding metadata
  group_by(sample, primer) %>%
  mutate(
    nreads = n()
  ) %>%
  ungroup() %>%
  mutate(across(where(is.character), as.factor))

glimpse(df_raw)

length(unique(df_raw$sample))

#### Pre-filtering: Removing low read samples ####

nreads_per_sample <- df_raw %>%
  group_by(sample) %>%
  select(sample, primer, nreads) %>%
  unique() %>% 
  summarize(nreads = sum(nreads))

low_read_cutoff <- 10000

low_read_samples <- nreads_per_sample %>% filter(nreads <= low_read_cutoff) %>% pull(sample)

low_read_samples

saveRDS(low_read_samples, paste0("data/low_read_samples.rds")) 

# removing low read samples now

length(unique(df_raw$sample))

df_raw <- df_raw %>% filter(!sample %in% low_read_samples)

length(unique(df_raw$sample))

#### Pre-filtering: Removing low count primers ####

low_read_primers <- df_raw %>%
  ungroup() %>%
  select(sample, primer, nreads) %>%
  unique() %>%
  group_by(primer) %>%
  summarise(q_reads = quantile(nreads)[2]) %>%
  filter(q_reads < settings$nReadsCutOff ) %>%
  pull(primer)

df_raw %>%
  ungroup() %>% 
  select(sample, primer, nreads) %>% 
  unique() %>%
  mutate(primer = reorder(primer, nreads)) %>% # set order
  ggplot(aes(x = primer, y = nreads, color = primer %in% low_read_primers)) +
  geom_boxplot(alpha = 0.3) +  
  geom_hline( yintercept = settings$nReadsCutOff, col = 'red') + # adding CutOff
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


## filtering out low read primers now
# low_read_primers <- readRDS("data/low_read_primers.rds") 

length(unique(df_raw$primer))

df_raw <- df_raw %>% filter(!primer %in% low_read_primers)

length(unique(df_raw$primer))

# saveRDS(low_read_primers, paste0("data/low_read_primers.rds")) 

#### Calculating values after pre-filetering ####

df_raw <- 
    df_raw %>%
    ungroup() %>% 
    group_by(primer) %>%
    mutate( #primer-level calculations
           short_peak_pos = find_top_peaks_for_primer(rpos, from_pos = settings$testRange[1], to_pos = settings$testRange[2], return_peak = "short"),
           long_peak_pos = find_top_peaks_for_primer(rpos,  from_pos = settings$testRange[1], to_pos = settings$testRange[2], return_peak = "long")
          ) %>%
    group_by(sample, primer) %>%
    mutate( #sample-level calculations
      dens_peaks_ratio = peaks_ratio(rpos, 
                                     from_pos = settings$testRange[1], 
                                     to_pos = settings$testRange[2], 
                                     short_peak = short_peak_pos, 
                                     long_peak = long_peak_pos, 
                                     n = 5)
    ) %>%
    ungroup()

# gathering motif frequences
motif_frequences_per_sample <- df_raw %>%
  group_by(sample, seq) %>%
  filter(!str_detect(seq, 'N')) %>% 
  summarise(count = n()) %>%
  group_by(sample) %>%
  mutate(freq = count/sum(count)) %>%
  select(-count) %>%
  pivot_wider(names_from = seq, names_glue = "{.value}_{seq}", values_from = freq)


#### Exporting fragmentation matrix & saving data ####

fragmetationScoreMatrixRatio <- df_raw %>%
  left_join(motif_frequences_per_sample, by = "sample") %>%
  select(sample, primer, diagnosis, stage, contains("freq_"), dens_peaks_ratio) %>%
  unique() %>% 
  pivot_wider(names_from = primer, values_from = dens_peaks_ratio)

# skimr::skim(fragmetationScoreMatrixRatio)


saveRDS(fragmetationScoreMatrixRatio, file = "data/fragmetationScoreMatrixRatio.rds")
saveRDS(df_raw, file = "data/df_raw.rds")

