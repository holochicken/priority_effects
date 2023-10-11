library(tidyverse)

de_07 <- read_tsv("results/tables/de_07.tsv") %>% mutate(day = "07")
de_21 <- read_tsv("results/tables/de_21.tsv") %>% mutate(day = "21")
de_35 <- read_tsv("results/tables/de_35.tsv") %>% mutate(day = "35")

descriptions <-
  read_tsv("data/tables/gene_annotation.tsv") %>%
  select(gene, description, gene_biotype)

st01 <-
  bind_rows(de_07, de_21, de_35) %>%
  filter(significant) %>%
  select(day, gene, logFC, FDR, ) %>%
  left_join(descriptions) %>%
  mutate(description = if_else(is.na(description), "", description)) %>%
  write_tsv("results/tables/st01.tsv")
