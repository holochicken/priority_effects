# Community mci


library(distillR)
library(vegan)
library(tidyverse)
library(gridExtra)
library(sjPlot)

# Load data  -----------------
# Abundance table
mag_counts <-
  read_tsv(file = "data/mags/mag_counts.tsv") %>%
  arrange(., by_group = mag_id) %>%
  column_to_rownames(var = "mag_id")

# Metadata
dmm <- read_tsv("results/tables/dmm.tsv")

metadata <-
  read_tsv(file = "data/kpi/metadata.tsv") %>%
  filter(animal_code %in% colnames(mag_counts)) %>%
  mutate(
    sampling_time = factor(sampling_time, levels = c("7", "21", "35"))
  ) %>%
  left_join(dmm, by = 'animal_code') %>%
  mutate(enterotype = case_when(dmm == '2' ~ "distinct",
                                dmm != '2' ~ "standard"))

# MAG stats
stats <-
  read_tsv('data/mags/stats.tsv') %>%
  mutate(correction_factor = median(mag_length) / mag_length)

# Functional capacities
gifts <-
  read_tsv('results/tables/gifts_elements.tsv') %>%
  column_to_rownames(var = 'mag_id')

#Aggregate element-level GIFTs into the function level
functs <- to.functions(gifts,GIFT_db1)
domains <- to.domains(functs,GIFT_db1)


# 1- Standardization and correction  -----------------
# Correction
mag_weighted <-
  round(sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0) %>%
  t() %>%
  as.data.frame()

# Standardization
total <- decostand(mag_weighted, 'total') %>%
  t() %>%
  as.data.frame()

# Calculate community MCI
domains_community <- to.community(domains,sweep(total, 2, colSums(total), FUN = "/"), GIFT_db1)

domains_community <-
  domains_community %>%
  as.data.frame() %>%
  rownames_to_column(var = 'animal_code') %>%
  left_join(metadata, by = 'animal_code')


# -1 artean gaudenez, orduan eredu mixto generalizatua. Pisua aldiz numerikoa da
library(MASS)
library(emmeans)
library(sjPlot)


# Fit the GLMM with random effect -
model <- glmmPQL(Overall ~ log(age)*trial + sex + breed,
                 family = "quasibinomial",
                 random = ~1|pen,
                 data = domains_community)

summary(model)


# By day
model_emm_7  <- emmeans(model, "trial", by = "age", at = list(age = c(7))) # Tuckey-ren konfiantza tarteak
model_emm_21 <- emmeans(model, "trial", by = "age", at = list(age = c(21)))
model_emm_35 <- emmeans(model, "trial", by = "age", at = list(age = c(35)))

model_emm_7_to_plot <- data.frame(model_emm_7)
model_emm_21_to_plot <- data.frame(model_emm_21)
model_emm_35_to_plot <- data.frame(model_emm_35)

day7_p <-
  ggplot(model_emm_7_to_plot,
         aes(x = trial,
             y = exp(emmean)/(1 + exp(emmean)),
         color = trial)) +
  geom_point() +
  scale_color_manual(values = c('#56B4E9', '#009E73', '#D55E00')) +
  ggtitle("a) Day 7") +
  geom_errorbar(aes(ymin = exp(lower.CL)/(1 + exp(lower.CL)),
                    ymax = exp(upper.CL)/(1 + exp(upper.CL))),
                width = 0.1) +
  ylab("Overall community MCI") +
  scale_x_discrete(labels = c("Trial A","Trial B","Trial C")) +
  theme_minimal() +
  geom_text(aes(y = exp(upper.CL + 0.02)/(1 + exp(upper.CL))),
            label = c("A","B","C"),
            size = 5) +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 14))


day21_p <-
  ggplot(model_emm_21_to_plot,
         aes(x = trial,
             y = exp(emmean)/(1 + exp(emmean)),
         color = trial)) +
    geom_point() +
    scale_color_manual(values = c('#56B4E9', '#009E73', '#D55E00')) +
    ggtitle("b) Day 21") +
    geom_errorbar(aes(ymin = exp(lower.CL)/(1 + exp(lower.CL)),
                      ymax = exp(upper.CL)/(1 + exp(upper.CL)),
                      width = 0.1)) +
    ylab("Overall community MCI") +
    scale_x_discrete(labels = c("Trial A","Trial B","Trial C")) +
    theme_minimal() +
    geom_text(aes(y = exp(upper.CL + 0.02)/(1 + exp(upper.CL))),
                  label = c("A","B","B"),
                  size = 5) +
    theme(legend.position = "none",
          axis.title.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 14))


day35_p <-
  ggplot(model_emm_35_to_plot,
         aes(x = trial,
             y = exp(emmean)/(1 + exp(emmean)),
         color = trial)) +
    geom_point() +
    scale_color_manual(values = c('#56B4E9', '#009E73', '#D55E00')) +
    ggtitle("c) Day 35") +
    geom_errorbar(aes(ymin = exp(lower.CL)/(1 + exp(lower.CL)),
                      ymax = exp(upper.CL)/(1 + exp(upper.CL)),
                      width = 0.1)) +
    ylab("Overall community MCI") +
    scale_x_discrete(labels = c("Trial A","Trial B","Trial C")) +
    theme_minimal() +
    geom_text(aes(y = exp(upper.CL + 0.02)/(1 + exp(upper.CL))),
                  label = c("A","B","B"),
                  size = 5) +
    theme(legend.position = "none",
          axis.title.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 14))

grid.arrange(day7_p,day21_p, day35_p, ncol = 3, nrow = 1)
