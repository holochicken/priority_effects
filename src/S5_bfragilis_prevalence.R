# 5- Bacteroides fragilis_A prevalence- updated 04/09/2023


# packages
library(sjPlot)
library(MuMIn)
library(nlme)
library(vegan)
library(microbiome)
library(tidyverse)

# 1- Load data  -----------------
# Abundance table
mag_counts <-
  read_tsv(file = "data/mags/mag_counts.tsv") %>%
  arrange(., by_group = mag_id) %>%
  pivot_longer(!mag_id, names_to = "sample_id", values_to = "value") %>%
  filter(grepl('CC', sample_id)) %>%
  pivot_wider(names_from = "sample_id", values_from = "value") %>%
  column_to_rownames(var = 'mag_id')

# Metadata
dmm <-
  read_tsv("results/tables/dmm.tsv")

metadata <-
  read_tsv(file = "results/clean/metadata.tsv") %>%
  filter(animal_code %in% colnames(mag_counts)) %>%
  mutate(sampling_time = factor(sampling_time,
                                levels = c("7",
                                           "21",
                                           "35"))) %>%
  left_join(dmm, by = 'animal_code') %>%
  column_to_rownames(var = "animal_code")

# MAG stats
stats <-
  read_tsv( "data/mags/stats.tsv") %>%
  mutate(correction_factor = median(mag_length) / mag_length)


# 2- Standardization and correction  -----------------
# Correction
mag_weighted <-
  round(sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0) %>%
  t() %>%
  as.data.frame()

# 3- B. fragilis abundance distribution on day 7 ------
# Transform data
clr_table <-
  microbiome::transform(phyloseq::otu_table(mag_weighted,
                                            taxa_are_rows = F),
                        transform = "clr") %>%
  data.frame

B_fragilis_A <-
  clr_table %>%
  select(cmag_558)

d7 <-
  metadata %>%
  filter(sampling_time == '7') %>%
  filter(trial == 'CC')

B_fragilis_A_filt <-
  B_fragilis_A %>%
  rownames_to_column('animal_code') %>%
  filter(animal_code %in% row.names(d7)) %>%
  column_to_rownames('animal_code')

mean(row.names(B_fragilis_A_filt) == row.names(d7))
B_fragilis_A <- B_fragilis_A_filt$cmag_558

# Linear models - adding pen as random effect
M <- lme(B_fragilis_A ~ sex + breed + treatment + age,
         random = ~1|pen,
         data = d7)

plot(M)
summary(M)

r.squaredGLMM(M)
# pen variation: 0.7881083
0.8308582 - 0.04335823

rm(list = ls())
