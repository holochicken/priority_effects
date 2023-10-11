#!/usr/bin/env Rscript

# 1 Libraries ----
source("src/S11_02_tidybulk_helpers.R")
library(tidyverse)
library(tidybulk)

# 2 Import data ----
sample_metadata <-
  read_tsv("results/tables/sample_metadata.tsv") %>%
  filter(tissue == "caecum")

counts <-
  read_tsv("data/transcriptomics/counts_caecum.tsv") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "counts") %>%
  left_join(sample_metadata) %>%
  filter(
    !is.na(sequencing_batch), ## remove CB17.08
    !is.na(has_campylobacter)
  ) %>%
  group_by(sample) %>%
  filter(sum(counts) > 10e6) %>%
  ungroup() %>%
  tidybulk(sample, gene, counts) %>%
  filter(day == "21") %>%
  mutate( # Overwrite what sample table says
    has_campylobacter = if_else(trial == "CC", TRUE, FALSE)
  )


# 3 Distribution of animals ----
chickens <-
  counts %>%
  pivot_sample() %>%
  select(sample, trial, day, has_campylobacter, sex, breed)

# Chickens by day and trial
chickens %>%
  select(has_campylobacter) %>%
  table(useNA = "ifany")

# 4 Overall distribution ----
## 4.1 Scale ----
counts_scaled <-
  counts %>%
  keep_abundant(factor_of_interest = has_campylobacter) %>%
  scale_abundance() %>%
  adjust_abundance(~ has_campylobacter + batch_effect)

counts_scaled %>%
  pivot_longer(
    cols = c("counts", "counts_scaled", "counts_scaled_adjusted"),
    names_to = "source",
    values_to = "abundance"
  ) %>%
  ggplot(aes(x = abundance + 1, color = sample)) +
  geom_density() +
  facet_wrap(~source) +
  scale_x_log10()

## 4.2 PCA ----
# counts_scaled_pca <-
counts_scaled %>%
  reduce_dimensions(
    top = 25000,
    .abundance = counts_scaled_adjusted,
    method = "PCA",
    .dims = 10
  ) %>%
  select(contains("PC"), everything()) %>%
  pivot_sample() %>%
  # filter(PC4 > 100)  # CC12.10
  # filter(PC2 > 150)  # CA04.12
  # filter(PC4 < -100 | PC4 > 50)  # CC21.07 CC04.08 CC18.12
  select(PC1:PC10, everything()) %>%
  ggpairs_pca(has_campylobacter)
# CA04.12 CC12.10 CC21.07 CC04.08 CC18.12

# 5 DE ----
de <-
  counts_scaled %>%
  test_differential_abundance(~ has_campylobacter + batch_effect) %>%
  annotate_differential_expression() %>%
  write_tsv("results/tables/de_21.tsv")

plot_ma(de)
plot_volcano(de)

de %>%
  filter(significant) %>%
  pull(direction) %>%
  table()

de %>%
  filter(gene %in% yada_abundant) %>%
  arrange(gene)


# 6 Enrichment
BiocManager::install("clusterProfiler")
source("src/S11_03_clusterprofiler_helpers.R")
de %>%
  mutate(direction = 1) %>%
  compute_go_terms("up") %>%
  clusterProfiler::cnetplot(showCategory = 100)
compute_go_terms(de, "up") %>% clusterProfiler::cnetplot(showCategory = 100)
compute_go_terms(de, "down") %>% clusterProfiler::cnetplot(showCategory = 100)

de %>%
  mutate(direction = 1) %>%
  compute_kegg_terms("up") %>%
  clusterProfiler::cnetplot(showCategory = 100)
compute_kegg_terms(de, "up") %>% clusterProfiler::cnetplot(showCategory = 100)
compute_kegg_terms(de, "down") %>% clusterProfiler::cnetplot(showCategory = 100)

remove.packages("clusterProfiler")



counts_scaled %>%
  filter(gene == "TNFR1") %>%
  ggplot(aes(x = sample, y = counts_scaled_adjusted, group = has_campylobacter, fill = has_campylobacter)) +
  geom_bar(stat = "identity") +
  scale_y_log10()
