#!/usr/bin/env Rscript

library(tidyverse)

batch2_rrna <- c(
  "CA21.02",
  "CB20.02", "CB23.14",
  "CC01.10", "CC02.16", "CC03.02", "CC04.08", "CC04.09", "CC04.16", "CC05.03",
  "CC05.04", "CC06.05", "CC06.17", "CC07.04", "CC07.07",
  # "CC08.06",  # This one lies with the healthy samples in PC1-PC2
  "CC11.11", "CC11.18", "CC12.09", "CC12.10", "CC13.11",
  "CC15.16", "CC16.05", "CC16.10", "CC17.11", "CC22.17"
)

batches <-
  read_tsv("data/transcriptomics/mrna_seq_batches.txt") %>%
  mutate(
    tissue = case_when(
      tissue == "E1a" ~ "Caecum",
      tissue == "B1a" ~ "Ileum"
    ),
    lab_batch = case_when(
      sequencing_batch == 3 ~ "lab2",
      TRUE ~ "lab1"
    ),
    sequencing_batch = str_glue("seq{sequencing_batch}"),
    ribo_batch = case_when(
      animal_code %in% batch2_rrna ~ "ribo1",
      sequencing_batch == "seq3" ~ "ribo2",
      TRUE ~ "noribo"
    )
  )

# Load sex correction

sex_correction <-
  read_tsv("data/transcriptomics/sex_predictions.tsv") %>%
  rename(sample = animal_code)


sample_data <-
  read_tsv(file = "data/kpi/metadata.tsv") %>%
  filter(age > 0) %>%
  select(
    animal_code, trial, sampling_time, sex, age, breed, treatment, has_campylobacter
  ) %>%
  right_join(batches) %>% # CB17.08 disappears here
  rename(sample = animal_code) %>%
  mutate(
    sex = str_to_lower(sex),
    breed = str_to_lower(breed),
    tissue = str_to_lower(tissue),
    day = str_remove(sampling_time, "^Day ") %>% str_pad(width = 2, pad = "0")
  ) %>%
  select(-sampling_time) %>%
  left_join(sex_correction) %>%
  mutate(
    sex = case_when(
      !is.na(predicted_sex) ~ predicted_sex,
      is.na(predicted_sex) ~ sex
    ),
    batch_effect = str_glue("{sequencing_batch}_{ribo_batch}_{breed}_{sex}"),
  ) %>%
  select(-predicted_sex) %>%
  write_tsv("data/metatranscriptomics/sample_metadata.tsv")
