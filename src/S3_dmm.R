# Community structure  using Dirichlet Multinomial Mixture models - updated 04/09/2023
# https://microbiome.github.io/tutorials/DMM.html

# packages
library(microbiome)
library(phylosmith)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(gridExtra)
library(tidyverse)

# Load data  -----------------
# Abundance table
mag_counts <-
  read_tsv(file = "data/mags/mag_counts.tsv") %>%
  arrange(by_group = mag_id) %>%
  column_to_rownames('mag_id')

metadata <-
  read_tsv(file = "data/kpi/metadata.tsv") %>%
  filter(animal_code %in% colnames(mag_counts)) %>%
  mutate(
    sampling_time = factor(sampling_time, levels = c('7', '21', '35'))) %>%
  column_to_rownames('animal_code')

# MAG stats
stats <-
  read_tsv("data/mags/stats.tsv") %>%
  mutate(correction_factor = median(mag_length) / mag_length)

# MAG taxonomy
taxonomy <-
  read_tsv("data/mags/taxonomy.tsv") %>%
  column_to_rownames(var = 'mag_id')


# 1- Standardization and correction  -----------------
# Correction
mag_weighted <-
  round(sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0)

rm(stats, mag_counts)


# 2- Mag community structure analysis -----------------
# Phyloseq objects
mag_counts_pseq <- otu_table(mag_weighted, taxa_are_rows = TRUE)
metadata_pseq <- sample_data(metadata)
taxonomy_pseq <- tax_table(as.matrix(taxonomy))
pseq <- phyloseq(mag_counts_pseq, metadata_pseq, taxonomy_pseq)
pseq <- set_sample_order(pseq, 'sampling_time')

rm(taxonomy, mag_weighted, mag_counts_pseq, taxonomy_pseq, metadata_pseq)


## Day 7 - DMM  -----------------
pseq_day7 <- subset_samples(pseq,sampling_time == "7")
day7_mags <- t(abundances(pseq_day7))

set.seed(1)
day7_fit <- lapply(1:5, dmn, count = day7_mags, verbose = TRUE)

day7_lplc <- base::sapply(day7_fit, DirichletMultinomial::laplace)
day7_aic <- base::sapply(day7_fit, DirichletMultinomial::AIC)
day7_bic <- base::sapply(day7_fit, DirichletMultinomial::BIC)

plot(day7_lplc, type = 'b', xlab = 'Number of Dirichlet Components', ylab = 'Model Fit')
lines(day7_aic, type = 'b', lty = 2)
lines(day7_bic, type = 'b', lty = 3)
(day7_best = day7_fit[[which.min(unlist(day7_lplc))]])

day7_ass <- apply(mixture(day7_best), 1, which.max)

rm(day7_fit, day7_lplc, day7_aic, day7_bic)


## Day 21 - DMM -----------------
set.seed(1)
pseq_day21 <- subset_samples(pseq,sampling_time == "21")
day21_mags <- t(abundances(pseq_day21))

set.seed(1)
day21_fit <- lapply(1:5, dmn, count = day21_mags, verbose = TRUE)

day21_lplc <- base::sapply(day21_fit, DirichletMultinomial::laplace)
day21_aic <- base::sapply(day21_fit, DirichletMultinomial::AIC)
day21_bic <- base::sapply(day21_fit, DirichletMultinomial::BIC)

plot(day21_lplc, type = 'b', xlab = 'Number of Dirichlet Components', ylab = 'Model Fit')
lines(day21_aic, type = 'b', lty = 2)
lines(day21_bic, type = 'b', lty = 3)
(day21_best = day21_fit[[which.min(unlist(day21_lplc))]])

day21_ass <- apply(mixture(day21_best), 1, which.max)

rm(day21_fit, day21_lplc, day21_aic, day21_bic)


## Day 35 - DMM  -----------------
set.seed(1)
pseq_day35 <- subset_samples(pseq,sampling_time == "35")
day35_mags <- t(abundances(pseq_day35))

set.seed(1)
day35_fit <- lapply(1:5, dmn, count = day35_mags, verbose = TRUE)

day35_lplc <- base::sapply(day35_fit, DirichletMultinomial::laplace)
day35_aic <- base::sapply(day35_fit, DirichletMultinomial::AIC)
day35_bic <- base::sapply(day35_fit, DirichletMultinomial::BIC)

plot(day35_lplc, type = 'b', xlab = 'Number of Dirichlet Components', ylab = 'Model Fit')
lines(day35_aic, type = 'b', lty = 2)
lines(day35_bic, type = 'b', lty = 3)
(day35_best = day35_fit[[which.min(unlist(day35_lplc))]])

day35_ass <- apply(mixture(day35_best), 1, which.max)

rm(day35_fit, day35_lplc, day35_aic, day35_bic)


# Save tables for lates analyses
d7_ass <-
  day7_ass %>%
  as.data.frame() %>%
  dplyr::rename(dmm = '.') %>%
  rownames_to_column('animal_code')

d21_ass <-
  day21_ass %>%
  as.data.frame() %>%
  dplyr::rename(dmm = '.') %>%
  rownames_to_column('animal_code')

d35_ass <-
  day35_ass %>%
  as.data.frame() %>%
  dplyr::rename(dmm = '.') %>%
  rownames_to_column('animal_code')

dmm <-
  bind_rows(d7_ass, d21_ass, d35_ass) %>%
  write_tsv("results/tables/dmm.tsv")


# 3- Community drivers - Top MAGs  -----------------
taxonomy <- read_tsv("data/mags/taxonomy.tsv")

order_colors <- c(Methanomicrobiales = "#f2f2f2",
                  Methanobacteriales = "#bfbfbf",
                  Saccharimonadales = "#dacce3",
                  Campylobacterales = "#cdadca",
                  Synergistales = "#c08eb3",
                  Methanomassiliicoccales = "#777dae",
                  Victivallales = "#0066ae",
                  Opitutales = "#1c7ebc",
                  Verrucomicrobiales = "#4a96cc",
                  Desulfovibrionales = "#68b0dc",
                  Enterobacterales = "#60cfde",
                  RF32 = "#60dfd2",
                  Burkholderiales = "#5cdfb5",
                  'Rs-D84' = "#40df91",
                  Bacteroidales = "#4dc87c",
                  Flavobacteriales = "#88c88b",
                  Gastranaerophilales = "#92e09f",
                  Coriobacteriales = "#c9e0af",
                  Actinomycetales = "#d8e093",
                  Deferribacterales = "#dff77e",
                  RF39 = "#ecf76d",
                  Lactobacillales = "#feef68",
                  'ML615J-28' = "#fde671",
                  Erysipelotrichales = "#ffd366",
                  Acholeplasmatales = "#fdc151",
                  Bacillales = "#fc953d",
                  RFN20 = "#fd8035",
                  Oscillospirales = "#fd8854",
                  TANB77 = "#f4814d",
                  Christensenellales = "#c26340",
                  Lachnospirales = "#c28c5c",
                  Clostridiales = "#ce9360",
                  Monoglobales_A = "#ce734f",
                  Peptostreptococcales = "#c65631",
                  Monoglobales = "#b93725",
                  UBA1212 = "#b10f19",
                  UBA4068 = "#8b1222",
                  Selenomonadales = "#7a1a12",
                  Acidaminococcales = "#64121f",
                  Veillonellales = "#5a0c37"
)

## Day 7  -----------------
## Day 7 - Enterotype 1
Day7_CT1 <- melt(fitted(day7_best))
colnames(Day7_CT1) <- c("MAG", "cluster", "value")

Day7_CT1_filtered <-
  subset(Day7_CT1, cluster == 1) %>%
  # Arrange MAGs by assignment strength
  arrange(value) %>%
  mutate(mag_id = factor(MAG, levels = unique(MAG))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.97)) %>%
  left_join(taxonomy, by = 'mag_id') %>%
  mutate(mag_id = reorder(mag_id, value))

p_Day7_CT1 <-
  ggplot(Day7_CT1_filtered, aes(x = mag_id, y = value, fill = order)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = order_colors) +
  # ylim(0,15) +
  labs(title = "a) d7 - community type 1") +
  theme_bw() +
  theme(legend.position = 'none') +
  coord_flip()

## Day 7 - Enterotype 2
Day7_CT2 <- melt(fitted(day7_best))
colnames(Day7_CT2) <- c("MAG", "cluster", "value")

Day7_CT2_filtered <-
  subset(Day7_CT2, cluster == 2) %>%
  arrange(value) %>%
  mutate(mag_id = factor(MAG, levels = unique(MAG))) %>%
  filter(abs(value) > quantile(abs(value), 0.97)) %>%
  left_join(taxonomy, by = 'mag_id') %>%
  mutate(mag_id = reorder(mag_id, value))

p_Day7_CT2 <-
  ggplot(Day7_CT2_filtered, aes(x = mag_id, y = value, fill = order)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = order_colors) +
  labs(title = "b) d7 - community type 2") +
  theme_bw() +
  theme(legend.position = 'none') +
  coord_flip()


# Day 21  -----------------
## Day 21 - Enterotype 1
Day21_CT1 <- melt(fitted(day21_best))
colnames(Day21_CT1) <- c("MAG", "cluster", "value")

Day21_CT1_filtered <-
  subset(Day21_CT1, cluster == 1) %>%
  # Arrange MAGs by assignment strength
  arrange(value) %>%
  mutate(mag_id = factor(MAG, levels = unique(MAG))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.97)) %>%
  left_join(taxonomy, by = 'mag_id') %>%
  mutate(mag_id = reorder(mag_id, value))

p_Day21_CT1 <-
  ggplot(Day21_CT1_filtered, aes(x = mag_id, y = value, fill = order)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = order_colors) +
  # ylim(0,15) +
  labs(title = "a) d21 - type 1") +
  theme_bw() +
  theme(legend.position = 'none') +
  coord_flip()

## Day 7 - Enterotype 2
Day21_CT2 <- melt(fitted(day21_best))
colnames(Day21_CT2) <- c("MAG", "cluster", "value")

Day21_CT2_filtered <-
  subset(Day21_CT2, cluster == 2) %>%
  arrange(value) %>%
  mutate(mag_id = factor(MAG, levels = unique(MAG))) %>%
  filter(abs(value) > quantile(abs(value), 0.97)) %>%
  left_join(taxonomy, by = 'mag_id') %>%
  mutate(mag_id = reorder(mag_id, value))

p_Day21_CT2 <-
  ggplot(Day21_CT2_filtered, aes(x = mag_id, y = value, fill = order)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = order_colors) +
  labs(title = "b) d21 - type 2") +
  theme_bw() +
  theme(legend.position = 'none') +
  coord_flip()


# Day 35  -----------------
## Day 35 - Enterotype 1
Day35_CT1 <- melt(fitted(day35_best))
colnames(Day35_CT1) <- c("MAG", "cluster", "value")

Day35_CT1_filtered <-
  subset(Day35_CT1, cluster == 1) %>%
  # Arrange MAGs by assignment strength
  arrange(value) %>%
  mutate(mag_id = factor(MAG, levels = unique(MAG))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.97)) %>%
  left_join(taxonomy, by = 'mag_id') %>%
  mutate(mag_id = reorder(mag_id, value))

p_Day35_CT1 <-
  ggplot(Day35_CT1_filtered, aes(x = mag_id, y = value, fill = order)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = order_colors) +
  # ylim(0,15) +
  labs(title = "a) d35 - type 1") +
  theme_bw() +
  theme(legend.position = 'none') +
  coord_flip()

## Day 35 - Enterotype 2
Day35_CT2 <- melt(fitted(day35_best))
colnames(Day35_CT2) <- c("MAG", "cluster", "value")

Day35_CT2_filtered <-
  subset(Day35_CT2, cluster == 2) %>%
  arrange(value) %>%
  mutate(mag_id = factor(MAG, levels = unique(MAG))) %>%
  filter(abs(value) > quantile(abs(value), 0.97)) %>%
  left_join(taxonomy, by = 'mag_id') %>%
  mutate(mag_id = reorder(mag_id, value))

p_Day35_CT2 <-
  ggplot(Day35_CT2_filtered, aes(x = mag_id, y = value, fill = order)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = order_colors) +
  labs(title = "b) d35 - type 2") +
  theme_bw() +
  theme(legend.position = 'none') +
  coord_flip()

## Day 35 - Enterotype 3
Day35_CT3 <- melt(fitted(day35_best))
colnames(Day35_CT3) <- c("MAG", "cluster", "value")

Day35_CT3_filtered <-
  subset(Day35_CT3, cluster == 3) %>%
  arrange(value) %>%
  mutate(mag_id = factor(MAG, levels = unique(MAG))) %>%
  filter(abs(value) > quantile(abs(value), 0.97)) %>%
  left_join(taxonomy, by = 'mag_id') %>%
  mutate(mag_id = reorder(mag_id, value))

p_Day35_CT3 <-
  ggplot(Day35_CT3_filtered, aes(x = mag_id, y = value, fill = order)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = order_colors) +
  labs(title = "b) d35 - type 3") +
  theme_bw() +
  theme(legend.position = 'none') +
  coord_flip()


p_blank <- ggplot() + theme_minimal()

pdf(file = "results/figures/supp_fig_s1.pdf", width = 7, height = 7)

grid.arrange(p_Day7_CT1, p_Day7_CT2, p_blank,
             p_Day21_CT1, p_Day21_CT2, p_blank,
             p_Day35_CT1, p_Day35_CT2, p_Day35_CT3,
             ncol = 3)

dev.off()

rm(list = ls())
