# Temporal changes and sources of variation on the expression of the microbiota
# - updated 04/09/2023


# Packages
library(data.table)
library(devtools)
library(tidyverse)
library(reshape2)
library(compositions)
library(nlme)
library(ggplot2)
library(distillR)
library(vegan)
library(ggrepel)
library(grid)
library(gridExtra)
library(phyloseq)
library(microbiome)

# 1- Load distillates ---------------------
load("data/metatranscriptomics/distilled_caecum_1.RData")
load("data/metatranscriptomics/distilled_caecum_2.RData")
load("data/metatranscriptomics/distilled_caecum_3.RData")
load("data/metatranscriptomics/distilled_caecum_4.RData")
load("data/metatranscriptomics/distilled_caecum_5.RData")
load("data/metatranscriptomics/distilled_caecum_6.RData")
load("data/metatranscriptomics/distilled_caecum_7.RData")
load("data/metatranscriptomics/distilled_caecum_8.RData")

#Merge distillates
dist_exp <- c(distilled_expression_caecum_1,
              distilled_expression_caecum_2,
              distilled_expression_caecum_3,
              distilled_expression_caecum_4,
              distilled_expression_caecum_5,
              distilled_expression_caecum_6,
              distilled_expression_caecum_7,
              distilled_expression_caecum_8)

rm(distilled_expression_caecum_1,
   distilled_expression_caecum_2,
   distilled_expression_caecum_3,
   distilled_expression_caecum_4,
   distilled_expression_caecum_5,
   distilled_expression_caecum_6,
   distilled_expression_caecum_7,
   distilled_expression_caecum_8)


# 2- Element-based ---------------------
#Aggregate GIFTs to elements
dist_exp_elements_by_mag <- lapply(dist_exp,function(x) to.elements(x,GIFT_db))

#Convert per-MAG to per-sample list
dist_exp_elements_by_sample <- sweep_matrix_list(dist_exp_elements_by_mag)

#Convert correct tibble for applying CLR conversion
dist_exp_elements_mrgd <-
  dist_exp_elements_by_mag %>%
  lapply(function(x) t(x)) %>%
  lapply(function(x) as.data.frame(x)) %>%
  Map(cbind, ., MAG = names(.)) %>%
  lapply(function(x) rownames_to_column(x, "Element")) %>% # move row names to first column
  do.call(rbind, .) %>% # bind all tables
  as.data.frame() %>% # convert to data.frame
  relocate(MAG, .before = Element) %>%
  rename_with(~ gsub("F1a", "", .x, fixed = TRUE)) # remove F1a suffix from sample names

# 3- CLR normalisation at the Element level ---------------------
# apply CLR conversion (keeping MAG and element information)
mag_vector <- dist_exp_elements_mrgd$MAG
elements_vector <- unique(dist_exp_elements_mrgd$Element)

dist_exp_elements_mrgd_clr_mrgd <-
  dist_exp_elements_mrgd %>%
  group_by(Element) %>%
  summarise_at(vars(-MAG), sum) %>%
  select(-Element) %>%
  otu_table(.,taxa_are_rows = TRUE) %>%
  transform("clr") %>%
  as.data.frame() %>%
  mutate(Element = elements_vector) %>%
  relocate(Element)


# 4- Plot GIFT ordination ---------------------
mci_mt <- as.data.frame(dist_exp_elements_mrgd_clr_mrgd)
rownames(mci_mt) <- mci_mt[,1]
mci_mt <- mci_mt[,-1]
colnames(mci_mt) <- gsub("F1a","",colnames(mci_mt))
mci_mt <- t(mci_mt)


# Load metadata
dmm <-
  read_tsv("results/tables/dmm.tsv")

metadata_mt <-
  read_tsv(file = "data/kpi/metadata.tsv") %>%
  filter(animal_code %in% rownames(mci_mt)) %>%
  mutate(sampling_time = factor(sampling_time, levels = c("7",
                                                          "21",
                                                          "35"))) %>%
  left_join(dmm, by = 'animal_code') %>%
  mutate(enterotype = case_when(dmm == '2' ~ "distinct",
                                dmm != '2' ~ "standard")) %>%
  column_to_rownames(var = "animal_code")

design_mt <-
  metadata_mt %>%
  select(pen, trial, age, sampling_time, enterotype, breed, sex, treatment)

# Check both tables
mci_mt <- mci_mt[rownames(mci_mt) %in% rownames(design_mt),]
design_mt <- design_mt[match(rownames(mci_mt),rownames(design_mt)),]
mean(rownames(design_mt) == rownames(mci_mt))

# Anova
mci_mt_red <- mci_mt[!is.na(design_mt$enterotype),]
design_mt_red <- design_mt[!is.na(design_mt$enterotype),]
summary(rda(mci_mt_red ~ enterotype * age, data = design_mt_red))
anova(rda(mci_mt_red ~ enterotype * age, data = design_mt_red), by = "terms")
RsquareAdj(rda(mci_mt_red ~ enterotype * age, data = design_mt_red))

# Perform dbRDA
db_rda_scores <-
  rda(mci_mt_red ~ enterotype * age, data = design_mt_red) %>%
  scores(display = "wa") %>%
  as.data.frame()

db_rda_scores_toplot <- data.frame(db_rda_scores,design_mt_red)

db_rda_module_scores <-
  rda(mci_mt_red ~ enterotype * age, data = design_mt_red) %>%
  scores(display = "species") %>%
  as.data.frame()

db_rda_module_scores$Function <- substring(rownames(db_rda_module_scores),first = 1,last = 3)
db_rda_module_scores_func <- aggregate(db_rda_module_scores[,-3],by = list(db_rda_module_scores$Function),mean)

# Remove structural modules
db_rda_module_scores_red <- db_rda_module_scores[-c((nrow(db_rda_module_scores) - 5):(nrow(db_rda_module_scores))),]

# Remove structural functions
db_rda_module_scores_func <- db_rda_module_scores_func[-c(19:21),]
db_rda_constraining_scores <-
  rda(mci_mt_red ~ enterotype*age, data = design_mt_red) %>%
  scores(display = "bp") %>%
  as.data.frame()

# dbRDA arrows and centroids
db_rda_constraining_scores_enterotype <- db_rda_constraining_scores[1,]*-1
db_rda_constraining_scores_age <- db_rda_constraining_scores[2,]

Day7A_centroid <- as.data.frame(t(as.matrix(apply(db_rda_scores_toplot[db_rda_scores_toplot$trial == "CA" &
                                                                       db_rda_scores_toplot$sampling_time == 7,1:2],
                                                  MARGIN = 2,
                                                  FUN = mean))))

Day7B_centroid <- as.data.frame(t(as.matrix(apply(db_rda_scores_toplot[db_rda_scores_toplot$trial == "CB" &
                                                                       db_rda_scores_toplot$sampling_time == 7,1:2],
                                                  MARGIN = 2,
                                                  FUN = mean))))

Day7Cstandard_centroid <- as.data.frame(t(as.matrix(apply(db_rda_scores_toplot[db_rda_scores_toplot$trial == "CC" &
                                                                               db_rda_scores_toplot$enterotype == "standard" &
                                                                               db_rda_scores_toplot$sampling_time == 7,1:2],
                                                          MARGIN = 2,
                                                          FUN = mean))))

Day7Cdistinct_centroid <- as.data.frame(t(as.matrix(apply(db_rda_scores_toplot[db_rda_scores_toplot$trial == "CC" &
                                                                               db_rda_scores_toplot$enterotype == "distinct" &
                                                                               db_rda_scores_toplot$sampling_time == 7,1:2],
                                                          MARGIN = 2,
                                                          FUN = mean))))

Day21A_centroid <- as.data.frame(t(as.matrix(apply(db_rda_scores_toplot[db_rda_scores_toplot$trial == "CA" &
                                                                      db_rda_scores_toplot$sampling_time == 21,1:2],
                                                   MARGIN = 2,
                                                   FUN = mean))))

Day21B_centroid <- as.data.frame(t(as.matrix(apply(db_rda_scores_toplot[db_rda_scores_toplot$trial == "CB" &
                                                                      db_rda_scores_toplot$sampling_time == 21,1:2],
                                                   MARGIN = 2,
                                                   FUN = mean))))

Day21C_centroid <- as.data.frame(t(as.matrix(apply(db_rda_scores_toplot[db_rda_scores_toplot$trial == "CC" &
                                                                      db_rda_scores_toplot$sampling_time == 21,1:2],
                                                   MARGIN = 2,
                                                   FUN = mean))))

Day35A_centroid <- as.data.frame(t(as.matrix(apply(db_rda_scores_toplot[db_rda_scores_toplot$trial == "CA" &
                                                                      db_rda_scores_toplot$sampling_time == 35,1:2],
                                                   MARGIN = 2,
                                                   FUN = mean))))

Day35B_centroid <- as.data.frame(t(as.matrix(apply(db_rda_scores_toplot[db_rda_scores_toplot$trial == "CB" &
                                                                      db_rda_scores_toplot$sampling_time == 35,1:2],
                                                   MARGIN = 2,
                                                   FUN = mean))))

Day35C_centroid <- as.data.frame(t(as.matrix(apply(db_rda_scores_toplot[db_rda_scores_toplot$trial == "CC" &
                                                                      db_rda_scores_toplot$sampling_time == 35,1:2],
                                                   MARGIN = 2,
                                                   FUN = mean))))

windows(height = 7, width = 14)

grid.arrange(
  db_rda_scores_toplot %>%
    mutate(RDA1cent = case_when(trial == "CA" & sampling_time == 7 ~ as.numeric(Day7A_centroid[1]),
                                trial == "CB" & sampling_time == 7 ~ as.numeric(Day7B_centroid[1]),
                                trial == "CC" & sampling_time == 7 & enterotype == "standard" ~ as.numeric(Day7Cstandard_centroid[1]),
                                trial == "CC" & sampling_time == 7 & enterotype == "distinct" ~ as.numeric(Day7Cdistinct_centroid[1]),
                                trial == "CA" & sampling_time == 21 ~ as.numeric(Day21A_centroid[1]),
                                trial == "CB" & sampling_time == 21 ~ as.numeric(Day21B_centroid[1]),
                                trial == "CC" & sampling_time == 21 ~ as.numeric(Day21C_centroid[1]),
                                trial == "CA" & sampling_time == 35 ~ as.numeric(Day35A_centroid[1]),
                                trial == "CB" & sampling_time == 35 ~ as.numeric(Day35B_centroid[1]),
                                trial == "CC" & sampling_time == 35 ~ as.numeric(Day35C_centroid[1])),
           RDA2cent = case_when(trial == "CA" & sampling_time == 7 ~ as.numeric(Day7A_centroid[2]),
                                trial == "CB" & sampling_time == 7 ~ as.numeric(Day7B_centroid[2]),
                                trial == "CC" & sampling_time == 7 & enterotype == "standard" ~ as.numeric(Day7Cstandard_centroid[2]),
                                trial == "CC" & sampling_time == 7 & enterotype == "distinct" ~ as.numeric(Day7Cdistinct_centroid[2]),
                                trial == "CA" & sampling_time == 21 ~ as.numeric(Day21A_centroid[2]),
                                trial == "CB" & sampling_time == 21 ~ as.numeric(Day21B_centroid[2]),
                                trial == "CC" & sampling_time == 21 ~ as.numeric(Day21C_centroid[2]),
                                trial == "CA" & sampling_time == 35 ~ as.numeric(Day35A_centroid[2]),
                                trial == "CB" & sampling_time == 35 ~ as.numeric(Day35B_centroid[2]),
                                trial == "CC" & sampling_time == 35 ~ as.numeric(Day35C_centroid[2]))) %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0 , linetype = 2) +
    geom_point(mapping = aes(x = RDA1,
                             y = RDA2,
                             color = trial,
                             shape = sampling_time,
                             fill = enterotype),
               size = 1,
               stroke = 0.9,
               alpha = 0.7) +
    scale_shape_manual(values = c(21, 24, 22)) +
    scale_fill_manual(values = c('#D55E00', 'white')) +
    geom_segment(aes(x = RDA1,
                     y = RDA2,
                     xend = RDA1cent,
                     yend = RDA2cent,
                     linetype = enterotype,
                     color = trial)) +
    scale_color_manual(values = c("#56B4E9", "#009E73", "#D55E00")) +
    geom_segment(data = db_rda_constraining_scores_enterotype,
                 mapping = aes(x = 0,
                               y = 0,
                               xend = RDA1*2,
                               yend = RDA2*2),
                 arrow = arrow(length = unit(0.3,"cm"))) +
    geom_label_repel(data = db_rda_constraining_scores_enterotype,
                     mapping = aes(x = RDA1*2,
                                   y = RDA2*2,
                                   label = "Distinct enterotype")) +
    geom_segment(data = db_rda_constraining_scores_age,
                 mapping = aes(x = 0,
                               y = 0,
                               xend = RDA1*2,
                               yend = RDA2*2),
                 arrow = arrow(length = unit(0.3,"cm"))) +
    geom_label_repel(data = db_rda_constraining_scores_age,
                     mapping = aes(x = RDA1*2,
                                   y = RDA2*2,
                                   label = "Chicken age")) +
    ylab("RDA2 (17%)") +
    xlab("RDA1 (70%)") +
    theme_bw() +
    theme(plot.margin = unit(c(2,0.5,2,0.5),"cm"),
          panel.border = element_blank())
  ,
  ggplot() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_segment(data = db_rda_module_scores_func,
                 aes(x = 0,
                     xend = RDA1, y = 0, yend = RDA2),
                 arrow = arrow(length = unit(0.25, "cm")),
                 color = "black") +
    geom_label_repel(data = db_rda_module_scores_func,
                     aes(x = RDA1,
                         y = RDA2,
                         label = Group.1),
                     size = 3) +
    geom_segment(data = db_rda_constraining_scores_enterotype,
                 mapping = aes(x = 0,
                               y = 0,
                               xend = RDA1/4,
                               yend = RDA2/4),
                 arrow = arrow(length = unit(0.3,"cm"))) +
    geom_label_repel(data = db_rda_constraining_scores_enterotype,
                     mapping = aes(x = RDA1/4,
                                   y = RDA2/4,
                                   label = "Distinct enterotype")) +
    geom_segment(data = db_rda_constraining_scores_age,
                 mapping = aes(x = 0,
                               y = 0,
                               xend = RDA1/4,
                               yend = RDA2/4),
                 arrow = arrow(length = unit(0.3,"cm"))) +
    geom_label_repel(data = db_rda_constraining_scores_age,
                     mapping = aes(x = RDA1/4,
                                   y = RDA2/4,
                                   label = "Chicken age")) +
    xlab("RDA1 (70%)") +
    ylab("RDA2 (17%)") +
    theme_bw() +
    coord_cartesian(xlim = c(-0.4,0.4), ylim = c(-0.2,0.2)) +
    theme(plot.margin = unit(c(2,0.5,2,0.5),"cm"),
          panel.border = element_blank())
  ,nrow = 1
)

ggsave("results/figures/fig_3ab_db_rda_elements.pdf", height = 4, width = 12)
