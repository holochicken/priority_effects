# MAG gifts and community level mci values - updated 22/03/2023

# Packages
library(distillR)
library(vegan)
library(tidyverse)
db <-  GIFT_db1

# 1- Load data ---------------------
# Create directories
dir.create(path = "results/figures/", showWarnings = FALSE, recursive = TRUE)
dir.create(path = "results/tables/", showWarnings = FALSE, recursive = TRUE)

# MAG counts and metadata
mag_counts <-
  read_tsv(file = "data/mags/mag_counts.tsv") %>%
  arrange(., by_group = mag_id) %>%
  column_to_rownames(var = 'mag_id')

metadata <-
  read_tsv("data/kpi/metadata.tsv") %>%
  mutate(sampling_time = factor(sampling_time, levels = c('7',
                                                          '21',
                                                          '35')))

# MAG characteristics
mag_annotations <-
  read_tsv("data/mags/gene_annotations.tsv.xz")

ena_to_mag_id <-
  read_tsv("data/mags/ena_assembly_to_mag_id.tsv")

stats <-
  read_tsv("data/mags/stats.tsv") %>%
  mutate(correction_factor = median(mag_length) / mag_length)


# 2- Standardisation ---------------------
# By MAG genome length
mag_weighted <-
  round(sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0) %>%
  t() %>%
  data.frame()

# By sequencing depth
hel <-
  decostand(mag_weighted, 'hellinger') %>%
  t() %>%
  as.data.frame()


# 3- Distillate annotations ---------------------
gifts <- distill(mag_annotations, GIFT_db1, genomecol = 2, annotcol = c(9,10,19))

gifts_raw <-
  gifts %>%
  data.frame() %>%
  rownames_to_column(var = 'mag_name') %>%
  left_join(ena_to_mag_id %>%
              select(mag_name, mag_id),
            by = 'mag_name') %>%
  select(-mag_name) %>%
  relocate(mag_id) %>%
  write_tsv("results/tables/gifts_raw.tsv") %>%
  column_to_rownames('mag_id') %>%
  as.matrix()

genome_completeness <-
  stats %>%
  select(mag_id, completeness_score) %>%
  as.matrix()


# 4- Define completeness corretion ---------------------
# Correction function defined in
# https://github.com/anttonalberdi/DAMMA/blob/main/R/damma_correction.R
# Based on https://www.nature.com/articles/s43705-023-00221-z
gifts_correction <- function(MCI_table, genome_completeness, stats = TRUE) {
  #Sort Genomes and test matching
  MCI_table <- MCI_table[order(rownames(MCI_table)),]
  genome_completeness <- genome_completeness[order(genome_completeness[,1]),]
  if (all(as.character(rownames(MCI_table)) !=
          as.character(genome_completeness[,1])))
    stop("Pathway table is missing")
  #Get completeness values
  Genome_completeness <- as.numeric(genome_completeness[,2])
  #Create corrected fullness matrix
  MCI_table_corrected <- matrix(0,nrow = nrow(MCI_table),ncol = ncol(MCI_table))
  colnames(MCI_table_corrected) = colnames(MCI_table)
  rownames(MCI_table_corrected) = rownames(MCI_table)
  #Iterate modelling and correction for each function
  suppressWarnings(
    for (i in 1:ncol(MCI_table)) {
      Model <- glm(MCI_table[,i]~Genome_completeness,family = "binomial")
      slope_coef <- coef(Model)[2]
      if (slope_coef > 0) {
        for (j in 1:nrow(MCI_table)) {
          # Model prediction of fullness if completeness was 100%
          pred_100 <- round(predict(Model,
                                    newdata = data.frame(
                                      Genome_completeness = 100),
                                    type = "response"),1)
          # Model prediction of fullness for actual completeness of the focal MAG
          pred_focal <- round(predict(Model,
                                      newdata = data.frame(
                                        Genome_completeness = Genome_completeness[j]),
                                      type = "response"), 1)
          # The expected change in function fullness if focal MAG was 100% complete
          pred_diff <- pred_100 - pred_focal
          MCI_table_corrected[j,i] = MCI_table[j,i] + pred_diff
        }
      } else if (slope_coef <= 0) {
        MCI_table_corrected[,i] <- MCI_table[,i]
      }
    }
  )
  # If corrected fullness >1, convert it to 1.
  MCI_table_corrected[MCI_table_corrected > 1] <- 1
  # Outout overall correction statistics on screen
  total <- nrow(MCI_table)*ncol(MCI_table)
  changes <- c(MCI_table == MCI_table_corrected)
  changes <- length(changes[!(changes)])
  percentage <- round(changes / total * 100,1)
  cat(paste0(changes," out of ",total," (",percentage,"%) fullness values were corrected\n"))
  if (stats == TRUE) {
    #Outout overall correction statistics on screen
    total <- ncol(MCI_table_corrected)
    for (r in rownames(MCI_table_corrected)) {
      completeness <- round(as.numeric(genome_completeness[genome_completeness[,1] == r,2]),1)
      changes <- MCI_table[r,] == MCI_table_corrected[r,]
      changes <- length(changes[!(changes)])
      percentage <- round(changes / total * 100,1)
      cat(paste0("\t",r," (",completeness,"%): ",changes,
                 "/",total," (",percentage,"%) fullness values were corrected\n"))
    }
  }
  #Output corrected table
  return(MCI_table_corrected)
}


# 5- Apply completeness corretion ---------------------
gifts_raw <- read_tsv("results/tables/gifts_raw.tsv") %>%
  column_to_rownames('mag_id')

gifts_corrected <- gifts_correction(gifts_raw, genome_completeness)

# Aggregate into compound level
gifts_elements <- to.elements(gifts_corrected,GIFT_db1)

gifts_elements <-
  gifts_elements %>%
  as.data.frame() %>%
  rownames_to_column(var = 'mag_id') %>%
  write_tsv("results/tables/gifts_elements.tsv") %>%
  column_to_rownames(var = 'mag_id')

# Aggregate into function level
gifts_elements <- read_tsv("results/tables/gifts_elements.tsv") %>%
  column_to_rownames(var = 'mag_id')

gifts_functions <- to.functions(gifts_elements, GIFT_db1)


rm(list = ls())
