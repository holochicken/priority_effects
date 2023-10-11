# Differential expression of microbiota - updated 04/09/2023

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

#Load distillates
load("src/mags/holofood_metatranscriptomics2/distillR/distilled_caecum_1.RData")
load("src/mags/holofood_metatranscriptomics2/distillR/distilled_caecum_2.RData")
load("src/mags/holofood_metatranscriptomics2/distillR/distilled_caecum_3.RData")
load("src/mags/holofood_metatranscriptomics2/distillR/distilled_caecum_4.RData")
load("src/mags/holofood_metatranscriptomics2/distillR/distilled_caecum_5.RData")
load("src/mags/holofood_metatranscriptomics2/distillR/distilled_caecum_6.RData")
load("src/mags/holofood_metatranscriptomics2/distillR/distilled_caecum_7.RData")
load("src/mags/holofood_metatranscriptomics2/distillR/distilled_caecum_8.RData")

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

#Convert correct tibble
dist_exp_elements_mrgd <-
  dist_exp_elements_by_mag %>%
  lapply(function(x) t(x)) %>%
  lapply(function(x) as.data.frame(x)) %>%
  Map(cbind, ., MAG = names(.)) %>%
  lapply(function(x) rownames_to_column(x, "Element")) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  relocate(MAG, .before = Element) %>%
  rename_with(~ gsub("F1a", "", .x, fixed = TRUE))

# 3- CLR normalisation at the Element level ---------------------
mag_vector <- dist_exp_elements_mrgd$MAG
elements_vector <- unique(dist_exp_elements_mrgd$Element) r

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


# 4- Differential expression analysis ---------------------
# Transform expression data
mci_mt <- as.data.frame(dist_exp_elements_mrgd_clr_mrgd)
rownames(mci_mt) <- mci_mt$Element
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

design_mt <- metadata_mt %>%
  select(pen, trial, age, sampling_time, enterotype, dmm, breed, sex, treatment)

## 4.1- Day 7 -----
# Filter day 7
mci_mt_7 <- mci_mt[rownames(mci_mt) %in% rownames(design_mt[which(design_mt$sampling_time == 7),]),]
design_mt_7 <- droplevels(design_mt[which(design_mt$sampling_time == 7),])

design_mt_7$Sample_ID <- rownames(design_mt_7)
design_mt_7 <- design_mt_7[match(rownames(mci_mt_7),rownames(design_mt_7)),]
mean(rownames(design_mt_7) == rownames(mci_mt_7))

# Run the GLM
p_val <- vector(mode = "numeric", length = ncol(mci_mt_7))
parameter <- vector(mode = "numeric", length = ncol(mci_mt_7))
p_val_adj <- vector(mode = "numeric", length = ncol(mci_mt_7))
results_mt <- data.frame(parameter, p_val, p_val_adj)
rownames(results_mt) <- colnames(mci_mt)

for (i in 1:ncol(mci_mt_7)) {
  comp <- mci_mt_7[,i]
  if (length(unique(mci_mt_7[,i])) > 9) {
    M <- lm(comp ~ dmm + trial + sex + breed + treatment + age, data = design_mt_7)
    results_mt[i,1] <- summary(M)$coefficients[2,1]
    results_mt[i,2] <- summary(M)$coefficients[2,4]
  }else{
    results_mt[i,1] <- NA
    results_mt[i,2] <- NA
  }
}

results_mt$p_val_adj <- p.adjust(results_mt$p_val, method = "fdr")
results_mt[order(results_mt$p_val),]
results_mt_day7 <- results_mt[results_mt$p_val < 0.05 & !is.na(results_mt$p_val),]
results_mt_day7$sampling_time <- "Day7"
results_mt_day7$element <- rownames(results_mt_day7)
results_mt_day7$funct <- substring(results_mt_day7$element, first = 1, last = 3)

# Drop structural element
results_mt_day7 <- results_mt_day7[!grepl("S",results_mt_day7$element),]


## 4.2- Day 21 ----
# Filter day 21
mci_mt_21 <- mci_mt[rownames(mci_mt) %in% rownames(design_mt[which(design_mt$sampling_time == 21),]),]
design_mt_21 <- droplevels(design_mt[which(design_mt$sampling_time == 21),])
design_mt_21[design_mt_21$trial == "CA", "trial"] <- "CAB"
design_mt_21[design_mt_21$trial == "CB", "trial"] <- "CAB"

design_mt_21$Sample_ID <- rownames(design_mt_21)
design_mt_21 <- design_mt_21[match(rownames(mci_mt_21),rownames(design_mt_21)),]
mean(rownames(design_mt_21) == rownames(mci_mt_21))

#Run the GLM
p_val <- vector(mode = "numeric",length = ncol(mci_mt_21))
parameter <- vector(mode = "numeric",length = ncol(mci_mt_21))
p_val_adj <- vector(mode = "numeric",length = ncol(mci_mt_21))
results_mt <- data.frame(parameter, p_val, p_val_adj)
rownames(results_mt) <- colnames(mci_mt)

 for (i in 1:ncol(mci_mt_21)) {
  if (length(unique(mci_mt_21[,i])) > 9) {
    comp <- mci_mt_21[,i]
    M <- lm(comp ~ trial + sex + breed + treatment + age, data = design_mt_21)
    results_mt[i,1] <- summary(M)$coefficients[2,1]
    results_mt[i,2] <- summary(M)$coefficients[2,4]
  }else{
    results_mt[i,1] <- NA
    results_mt[i,2] <- NA
  }
 }

results_mt$p_val_adj <- p.adjust(results_mt$p_val, method = "fdr")
results_mt[order(results_mt$p_val),]
results_mt_day21 <- results_mt[results_mt$p_val < 0.05 & !is.na(results_mt$p_val),]
results_mt_day21$sampling_time <- "Day21"
results_mt_day21$element <- rownames(results_mt_day21)
results_mt_day21$funct <- substring(results_mt_day21$element, first = 1, last = 3)

# Drop structural element
results_mt_day21 <- results_mt_day21[!grepl("S",results_mt_day21$element),]


## 4.3- Day 35 -----
#Filter day 35 only
mci_mt_35 <- mci_mt[rownames(mci_mt) %in% rownames(design_mt[which(design_mt$sampling_time == 35),]),]
design_mt_35 <- droplevels(design_mt[which(design_mt$sampling_time == 35),])
design_mt_35[design_mt_35$trial == "CA", "trial"] <- "CAB"
design_mt_35[design_mt_35$trial == "CB", "trial"] <- "CAB"

design_mt_35$Sample_ID <- rownames(design_mt_35)
design_mt_35 <- design_mt_35[match(rownames(mci_mt_35),rownames(design_mt_35)),]
mean(rownames(design_mt_35) == rownames(mci_mt_35))

#Run the GLM
p_val <- vector(mode = "numeric",length = ncol(mci_mt_35))
parameter <- vector(mode = "numeric",length = ncol(mci_mt_35))
p_val_adj <- vector(mode = "numeric",length = ncol(mci_mt_35))
results_mt <- data.frame(parameter, p_val, p_val_adj)
rownames(results_mt) <- colnames(mci_mt)

for (i in 1:ncol(mci_mt_35)) {
  if (length(unique(mci_mt_35[,i])) > 9) {
    comp <- mci_mt_35[,i]
    M <- lm(comp~ trial + sex + breed + treatment + age, data = design_mt_35)
    results_mt[i,1] <- summary(M)$coefficients[2,1]
    results_mt[i,2] <- summary(M)$coefficients[2,4]
  }else{
    results_mt[i,1] <- NA
    results_mt[i,2] <- NA
  }
}

results_mt$p_val_adj <- p.adjust(results_mt$p_val,method = "fdr")
results_mt[order(results_mt$p_val),]
results_mt_day35 <- results_mt[results_mt$p_val < 0.05 & !is.na(results_mt$p_val),]
results_mt_day35$sampling_time <- "Day35"
results_mt_day35$element <- rownames(results_mt_day35)
results_mt_day35$funct <- substring(results_mt_day35$element, first = 1, last = 3)

# Drop structural element
results_mt_day35 <- results_mt_day35[!grepl("S", results_mt_day35$element),]


## 5- Plot GLM results ---------------------
results_mt_toplot <- data.frame(rbind(results_mt_day7,results_mt_day21,results_mt_day35))
results_mt_toplot$toplot <- NA
results_mt_toplot$toplot[results_mt_toplot$parameter < 0 & results_mt_toplot$p_val < 0.05]  <- 1
results_mt_toplot$toplot[results_mt_toplot$parameter < 0 & results_mt_toplot$p_val_adj < 0.05]  <- 2
results_mt_toplot$toplot[results_mt_toplot$parameter > 0 & results_mt_toplot$p_val < 0.05] <- 1
results_mt_toplot$toplot[results_mt_toplot$parameter > 0 & results_mt_toplot$p_val_adj < 0.05] <- 2

results_mt_toplot$sampling_time <- factor(results_mt_toplot$sampling_time,
                                          levels = c("Day7", "Day21", "Day35"))

windows(height = 5, width = 12)

ggplot(data = results_mt_toplot,
       mapping = aes(x = sampling_time, y = element)) +
  geom_tile(mapping = aes(fill = toplot, width = 0.9, height = 0.9)) +
  scale_fill_gradientn(colours = c("darkblue", "lightblue", "red","darkred"),
                       guide = "none") +
  theme_bw() +
  facet_grid(. ~ funct, space = "free_x", scales = "free_x", switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  coord_flip()

ggsave("results/figures/fig_3c_heatmap_elements.pdf", height = 4, width = 18)
