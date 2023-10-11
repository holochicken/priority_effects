# Analysis of all metabolites produced by enterotypes
# Comparison of the production capacity between enterotypes


library(lme4)
library(lmerTest)
library(multcomp)
library(MASS)
library(tidyverse)


# 1- Load all metabolites -----------------
# Metadata
mapping_file <- read.table("data/mags/mag_trends.tsv", sep = "\t", header = T)

metadata <-
  read.table("results/tables/metadata_d7_all_trials.tsv", sep = "\t", header = T) %>%
  dplyr::select(names = animal_code, group = dmm, trial = trial,
                breed, sex, treatment, age, pen)

#Load relative abundance table
load("results/tables/animals_d7_all_trials.RData")


# 2- All metabolites -----------------
mag_mets <-
  read.table("data/mags/all_metabolites.tsv", sep = "\t", header = F) %>%
  dplyr::select(mag_id = 1, metabolites = 2)  %>%
  left_join(mapping_file, by = c("mag_id" = "mag_id")) %>%
  dplyr::select(mag_name = mag_name, metabolites = metabolites)

all_mets <-
  unlist(strsplit(mag_mets$metabolites,",")) %>%
  gsub(" ","",.) %>%
  gsub("_c","",.) %>%
  unique()

length(all_mets)  # 4,533 metabolites

# Iterate computation across metabolites
all_met_table <- tibble(sample = names(table))
n = 0
for (metabolite in all_mets) {
  n = n + 1
  print(n)

  mag_mets2 <-
    mag_mets %>%
    mutate(value = if_else(str_detect(metabolites, metabolite), 1, 0)) %>%
    dplyr::select(mag_name, value)

  mag_mets2_results <-
    tibble(sample = names(table)) %>%
    mutate(values = map_dbl(seq_along(table), ~ {
      table[[.x]] %>%
        as.data.frame() %>%
        rownames_to_column(var = "mag_name") %>%
        left_join(mag_mets2, by = c("mag_name" = "mag_name")) %>%
        mutate(weight = round(rel_abu * value,2)) %>%
        summarize(weighted_average = sum(weight)) %>%
        pull(weighted_average) %>%
        as.numeric()
    }))

  all_met_table <- all_met_table %>%
    left_join(mag_mets2_results, by = c("sample" = "sample"))
}

write.table(all_met_table, "results/tables/all_metabolites_relabun.tsv", col.names = T , quote = F)


# 3- Perform linear models for each metabolite -----------------
# Add metadata to table
all_met_table <- read.table("results/tables/all_metabolites_relabun.tsv")

all_met_table_metadata <-
  all_met_table %>%
  left_join(metadata, by = c("sample" = "names"))

# Rename metabolite names to avoid issues (remove any non-alphabetical character)
all_met_table_metadata <-
  all_met_table_metadata %>%
  rename_all(~make.names(., unique = TRUE))

# Rename enterotypes to 0-1
all_met_table_metadata$group <- ifelse(all_met_table_metadata$group == 1, 0, 1)


# Linear model
results <- list()
# Loop through each explanatory variable
n = 0
for (metabolite in names(all_met_table_metadata)[c(2:4534)]) {
  n = n + 1
  print(n)
  tryCatch({
    # Create the formula for the GLM
    formula <- as.formula(paste(metabolite, "~ group + trial + sex +
                                               breed + treatment + age"))
    # Fit the GLM with random effect
    model <- glmmPQL(formula,
                     data = all_met_table_metadata,
                     family = "quasibinomial",
                     random = ~1|pen)
    # Calculate p-values and adjusted p-values
    p_value <- summary(model)$tTable[, "p-value"][2]
    #Get average values
    means <- all_met_table_metadata %>%
      group_by(group) %>%
      summarize(mean = mean(get(metabolite)))
    # Store the results in the list
    results[[metabolite]] <- tibble(Metabolite = metabolite,
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
results_table <- bind_rows(results)

# Adjust p-values
results_table$adjusted_p_value <- p.adjust(results_table$p_value, method = "fdr")
summary(results_table$adjusted_p_value)

#Number of increased metabolites
results_table %>%
  filter(adjusted_p_value < 0.05) %>%
  filter(difference > 0) %>%
  nrow()

# 1,043 overexpressed metabolites in the distinct enterotype

#Show top metabolites showing differences between enterotypes
results_table %>%
  arrange(-difference) %>%
  print(n = 60)


# 4- Ploting specific metabolites -----------------
### 1-
metabolite = "D.glucopyranose.6.phosphate" # possible sugar examples

### 2-
metabolite = "propionate" # possible SCFA examples

# Cumulative relative abundance of bacteria
# capable of producing selected metabolites
mag_mets2 <-
  mag_mets %>%
  mutate(value = if_else(str_detect(metabolites, metabolite), 1, 0))

mag_mets2_results <-
  tibble(names = names(table)) %>%
  mutate(values = map_dbl(seq_along(table), ~ {
    table[[.x]] %>%
      as.data.frame() %>%
      rownames_to_column(var = "mag_name") %>%
      left_join(mag_mets2, by = c("mag_name" = "mag_name")) %>%
      mutate(weight = rel_abu * value) %>%
      summarize(weighted_average = sum(weight)) %>%
      pull(weighted_average) %>%
      as.numeric()
  }))


# Plot selected metabolites
pdf(
  "results/figures/fig_2h_D.glucopyranose.6.phosphate.pdf",
  # "results/figures/fig_2i_propionate.pdf",
  width = 5, height = 3)

mag_mets2_results %>%
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
    "D-Glucopyranose-6-Phosphate - adj_p-value = 0.1436277"
    # "Propionate - adj_p-value = 1.731827e-01"
    )

dev.off()
