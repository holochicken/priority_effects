# Analysis of source metabolites for Campylobacyer produced by enterotypes


library(lme4)
library(lmerTest)
library(multcomp)
library(MASS)
library(tidyverse)
library(ggplot2)

# 1- Load MAG-CMAG mapping table  -----------------
mapping_file <- read.table("data/mags/ena_assembly_to_mag_id.tsv", sep = "\t", header = T)

metadata <-
  read.table("results/tables/metadata_d7_all_trials.tsv", sep = "\t", header = T) %>%
  dplyr::select(names = animal_code, group = dmm, trial = trial,
                breed, sex, treatment, age, pen)

#Load relative abundance table
load("results/tables/animals_d7_all_trials.RData")


# 2- Analysis of source metabolites of Campylobacter  -----------------
# Load depencence table and rename to cmag
dependence <-
  read.table("data/mags/dependency_ERR4836918_bin_11.tsv", sep = "\t", header = F) %>%
  dplyr::select(mag_id = 1, value = 3, metabolites = 4) %>%
  left_join(mapping_file, by = c("mag_id" = "mag_id")) %>%
  dplyr::select(mag_name = mag_name, value = value, metabolites = metabolites)

dependence_c_coli <-
  read.table("data/mags/dependency_ERR4836965_bin_9.tsv", sep = "\t", header = F) %>%
  dplyr::select(mag_id = 1, value = 3, metabolites = 4) %>%
  left_join(mapping_file, by = c("mag_id" = "mag_id")) %>%
  dplyr::select(mag_name = mag_name, value = value, metabolites = metabolites)

#Calculate dependencies for all samples
dependence_results <-
  tibble(names = names(table)) %>%
  mutate(values = map_dbl(seq_along(table), ~ {
    table[[.x]] %>%
      as.data.frame() %>%
      rownames_to_column(var = "mag_name") %>%
      left_join(dependence, by = c("mag_name" = "mag_name")) %>%
      mutate(weight = rel_abu * value) %>%
      summarize(weighted_average = sum(weight)) %>%
      pull(weighted_average) %>%
      as.numeric()
  }))


# Overall calculations
dependence_results %>%
  left_join(metadata, by = c("names" = "names")) %>%
  group_by(group) %>%
  summarize(mean = mean(values))

# 3- Number of source metabolites each enterotype is capable of producing   -----------------
# Plot
pdf("results/figures/fig_2b_source_met_ent.pdf", width = 6, height = 4)

dependence_results %>%
  left_join(metadata, by = c("names" = "names")) %>%
  ggplot(aes(x = as.factor(group), y = values, shape = trial, color = trial)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.3) +
  scale_fill_manual(values = c('#56B4E9', '#009E73', '#D55E00')) +
  scale_color_manual(values = c('#56B4E9', '#009E73', '#D55E00')) +
  theme_classic()

dev.off()

# Linear models
formula <- as.formula(paste("values ~ group + trial + sex + breed + treatment + age + (1|pen)"))
metadata$group <- as.factor(metadata$group)
# Fit the LMM with random effect
dependence_results %>%
  left_join(metadata, by = c("names" = "names")) %>%
  lmer(formula, data = .) %>%
  summary()

# Results: group2       2.22542    0.40962  89.39776   5.433 4.74e-07 ***


# 4- Metabolite-specific calculations  -----------------
# List metabolits
source_mets <-
  unlist(strsplit(dependence$metabolites,",")) %>%
  gsub(" ","",.) %>%
  gsub("_c","",.) %>%
  unique()  # 124 source metabolites

source_mets_c_coli <-
  unlist(strsplit(dependence_c_coli$metabolites,",")) %>%
  gsub(" ","",.) %>%
  gsub("_c","",.) %>%
  unique()  # 129 source metabolites

# Get source metabolites relative abundance table
source_met_table <- tibble(sample = names(table))
n = 0
for (metabolite in source_mets) {
  n = n + 1
  print(n)

  source_mets2 <- dependence %>%
    mutate(value = if_else(str_detect(metabolites, metabolite), 1, 0)) %>%
    dplyr::select(mag_name, value)

  source_mets2_results <- tibble(sample = names(table)) %>%
    mutate(values = map_dbl(seq_along(table), ~ {
      table[[.x]] %>%
        as.data.frame() %>%
        rownames_to_column(var = "mag_name") %>%
        left_join(source_mets2, by = c("mag_name" = "mag_name")) %>%
        mutate(weight = round(rel_abu * value,2)) %>%
        summarize(weighted_average = sum(weight)) %>%
        pull(weighted_average) %>%
        as.numeric()
    })) %>%
    rename(., !!metabolite := values)

  source_met_table <- source_met_table %>%
    left_join(source_mets2_results, by = c("sample" = "sample"))
}
# write.table(source_met_table,
#             "source_ERR4836918_bin_11_metabolites_relabun.tsv",
#             col.names = T, quote = F)


# 5- Compute per-metabolite enterotype differences   -----------------
# source_met_table <-
#   read.table("source_ERR4836918_bin_11_metabolites_relabun.tsv")

#Add metadata to table
source_met_table_metadata <-
  source_met_table %>%
  left_join(metadata, by = c("sample" = "names"))

#rename metabolite names to avoid issues
#(remove any nonalphabetical charaacter)
source_met_table_metadata <-
  source_met_table_metadata %>%
  rename_all(~make.names(., unique = TRUE))

#rename enterotypes to 0-1
source_met_table_metadata$group <- ifelse(source_met_table_metadata$group == 1, 0, 1)


source_met_results <- list()
# Loop through each explanatory variable
n = 0
for (metabolite in names(source_met_table_metadata)[c(2:125)]) {
  n = n + 1
  print(n)
  tryCatch({
    # Create the formula for the GLM
    formula <- as.formula(paste(metabolite, "~ group + trial + sex +
                                               breed + treatment + age"))
    # Fit the GLMM with random effect
    model <- glmmPQL(formula,
                     data = source_met_table_metadata,
                     family = "quasibinomial",
                     random = ~1|pen)
    # Calculate p-values and adjusted p-values
    p_value <- summary(model)$tTable[, "p-value"][2]
    # adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
    #Get average values
    means <- source_met_table_metadata %>%
      group_by(group) %>%
      summarize(mean = mean(get(metabolite)))
    # Store the results in the list
    source_met_results[[metabolite]] <- tibble(Metabolite = metabolite,
                                                  p_value = p_value,
                                                  standard = means[1,2],
                                                  distinct = means[2,2],
                                                  difference = means[2,2] - means[1,2])
  }, error = function(e) {
    # Code to handle the error (optional)
    print(paste("Error occurred:", e))
  })
}

# Combine the results into a single tibble
source_met_results_table <- bind_rows(source_met_results)

source_met_results_table$adjusted_p_value <- p.adjust(source_met_results_table$p_value, method = "fdr")
summary(source_met_results_table$adjusted_p_value)

source_met_results_table %>%
  filter(adjusted_p_value < 0.05) %>%
  filter(difference > 0) %>%
  nrow()

# 31 source metabolites diff capacity for C. jejuni
# 35 source metabolites diff capacity for C. coli


# 6- Exploring top increased and decreased metabolites
## 6.1 - Number of increased metabolites
source_met_results_table %>%
  filter(adjusted_p_value < 0.05) %>%
  filter(difference > 0) %>%
  arrange(-difference) %>%
  nrow()

# Top 30 of increased metabolites
source_met_results_table %>%
  filter(adjusted_p_value < 0.05) %>%
  filter(difference > 0) %>%
  arrange(-difference) %>%
  print(n = 30)

## 6.2- Number of decreased metabolites
source_met_results_table %>%
  filter(adjusted_p_value < 0.05) %>%
  filter(difference < 0) %>%
  arrange(difference) %>%
  nrow()

# Top 30 of decreased metabolites
source_met_results_table %>%
  filter(adjusted_p_value < 0.05) %>%
  filter(difference < 0) %>%
  print(n = 30)


# 7- Metabolite-specific visualisation  -----------------
# Largest overall difference
### 1-
metabolite = "COPROPORPHYRIN_III"  # https://www.biorxiv.org/content/biorxiv/early/2023/03/30/2023.03.30.534706.full.pdf
### 2-
metabolite = "CPD.31"  # citramalate
### 3-
metabolite = "SULFATE"
### 4-
metabolite = "MOCS3-Cysteine"  # https://onlinelibrary.wiley.com/doi/full/10.1111/mmi.12732


# Cumulative relative abundance of bacteria
# capable of producing selected metabolites
dependence2 <-
  dependence %>%
  mutate(value = if_else(str_detect(metabolites, metabolite), 1, 0))

dependence2_results <-
  tibble(names = names(table)) %>%
  mutate(values = map_dbl(seq_along(table), ~ {
    table[[.x]] %>%
      as.data.frame() %>%
      rownames_to_column(var = "mag_name") %>%
      left_join(dependence2, by = c("mag_name" = "mag_name")) %>%
      mutate(weight = rel_abu * value) %>%
      summarize(weighted_average = sum(weight)) %>%
      pull(weighted_average) %>%
      as.numeric()
  }))


# Plot selected metabolites
pdf(
  "results/figures/fig_2d_COPROPORPHYRIN_III.pdf",
  # "results/figures/fig_2e_CPD_31_citramalate.pdf",
  # "results/figures/fig_2g_SULFATE.pdf",
  # "results/figures/fig_2f_MOCS3-Cysteine.pdf",
  width = 5, height = 3)

dependence2_results %>%
  left_join(metadata, by = c("names" = "names")) %>%
  ggplot(aes(x = as.factor(group),
             y = values,
             shape = trial,
             color = trial,
             fill = trial)) +
  geom_boxplot(lwd = 0.5,
               outlier.color = NA,
               alpha = 0.3,
               position = position_dodge(1.2)) +
  geom_jitter(alpha = 0.8) +
  ylim(c(0,1)) +
  scale_fill_manual(values = c('#56B4E9', '#009E73', '#D55E00')) +
  scale_color_manual(values = c('#56B4E9', '#009E73', '#D55E00')) +
  theme_classic() +
  ggtitle(
    "Coproporphyrin_III - adj_p-value = 1.02e-14"
    # "CPD-31 - adj_p-value = 9.63e-7"
    # "Sulfate adj_p-value = 1.34e-8"
    # "MOCS3-Cysteine adj_p-value = 5.06e-6"
    )

dev.off()
