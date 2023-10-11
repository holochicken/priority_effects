#  Campylobacter relative abundances - updated 04/09/2023

# Packages
library(vegan)
library(tidyverse)

# 1- Load data ---------------------
# Abundance table
mag_counts <-
  read_tsv(file = "data/mags/mag_counts.tsv") %>%
  arrange(., by_group = mag_id) %>%
  pivot_longer(!mag_id, names_to = "sample_id", values_to = "value") %>%
  pivot_wider(names_from = "sample_id", values_from = "value") %>%
  column_to_rownames(var = 'mag_id')

# Metadata
metadata <-
  read_tsv(file = "data/kpi/metadata.tsv") %>%
  filter(animal_code %in% colnames(mag_counts)) %>%
  mutate(sampling_time = factor(sampling_time, levels = c("7", "21", "35"))) %>%
  dplyr::rename(pcr_campylo = has_campylobacter) %>%
  select(animal_code, trial, sampling_time, pcr_campylo)

stats <-
  read_tsv( "data/mags/stats.tsv") %>%
  mutate(correction_factor = median(mag_length) / mag_length)

# MAG taxonomy
taxonomy <-
  read_tsv("data/mags/taxonomy.tsv")


# 2- Standardization and correction  -----------------
# Correction
mag_weighted <-
  round(sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0) %>%
  t() %>%
  data.frame()

campylos <-
  decostand(mag_counts, "total") %>%
  t() %>%
  data.frame() %>%
  rownames_to_column("animal_code") %>%
  select(animal_code, cmag_519, cmag_395, cmag_332) %>%
  dplyr::rename(C_coli = cmag_519, C_jejuni = cmag_395, H_pullorum = cmag_332) %>%
  rowwise() %>%
  left_join(metadata) %>%
  gather(C_coli:H_pullorum,
         key = "mag",
         value = "counts") %>%
  mutate(mag = factor(mag, levels = c('C_jejuni',
                                      'C_coli',
                                      'H_pullorum'))) %>%
  mutate(pcr_campylo = if_else(is.na(pcr_campylo) & trial == 'CC', TRUE, pcr_campylo)) %>%
  mutate(pcr_campylo = if_else(is.na(pcr_campylo) & trial != 'CC', FALSE, pcr_campylo)) %>%
  mutate(counts_modified = ifelse(pcr_campylo == 'FALSE', '0', as.numeric(as.character(counts)))) %>%
  mutate(counts_modified = as.numeric(counts_modified) * 100)

# Plot
plot_camp <-
  ggplot(campylos) +
  geom_boxplot(aes(x = mag,
                   y = counts_modified,
                   color = trial),
               outlier.shape = NA,
               shape = 0.1) +
  geom_point(aes(x = mag,
                 y = counts_modified,
                 color = trial),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
             alpha = 0.5,
             size = 0.3) +
  scale_color_manual(values = c('#56B4E9', '#009E73', '#D55E00')) +
  facet_wrap(~ sampling_time) +
  ylab('Relative abundance') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6))


# Save plot
pdf(file = "results/figures/fig1a_campylo.pdf", width = 7, height = 5)

plot_camp + scale_y_continuous(trans = "log1p", breaks = c(0,2,10,15), labels = c('0', '0.02', '0.1', '0.15'))

dev.off()

rm(list = ls())
