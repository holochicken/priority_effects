# Cumulative mag counts for metabolic_networks - - updated 04/09/2023


# packages
library(GGally)
library(ggrepel)
library(gridExtra)
library(sjPlot)
library(tidyverse)
library(vegan)

# 1- Load data  -----------------
# Abundance table
mag_counts <- read_tsv(file = "data/mags/mag_counts.tsv") %>%
  arrange(., by_group = mag_id) %>%
  column_to_rownames(var = "mag_id")

# MAG stats
stats <- read_tsv("data/mags/stats.tsv") %>%
  mutate(correction_factor = median(mag_length) / mag_length)

# MAG taxonomy
taxonomy <- read_tsv("data/mags/taxonomy.tsv")

# 2- Standardization and correction  -----------------
# Correction
mag_weighted <- round(
  sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0
) %>%
  t() %>%
  as.data.frame()

# Standardization
total <- decostand(mag_weighted, "total")

# 3- Metadata and sample selection  -----------------
# set.seed(1)
dmm <- read_tsv("results/tables/dmm.tsv")

metadata <- read_tsv(file = "data/kpi/metadata.tsv")  %>%
  filter(animal_code %in% colnames(mag_counts)) %>%
  left_join(dmm, by = 'animal_code') %>%
  mutate(enterotype = case_when(dmm == '2' ~ "distinct",
                                dmm != '2' ~ "standard")) %>%
  filter(sampling_time == '7') %>%
  write_tsv(file = 'results/tables/metadata_d7_all_trials.tsv')

# 4- Calculate cumulative abundances
loop_values <- metadata$animal_code
animals_d7 <- list()

for (i in loop_values) {
  animals_d7[[i]] <-
    total %>%
    filter(rownames(.) == i) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename(rel_abu = i) %>%
    arrange(desc(rel_abu)) %>%
    mutate(cum_abu = cumsum(rel_abu)*100/sum(rel_abu))
}

save(animals_d7, file = 'results/tables/animals_d7_all_trials.RData')

