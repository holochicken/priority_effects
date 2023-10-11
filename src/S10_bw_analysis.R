# Libraries ----


library(compositions)
library(nlme)
library(emmeans)
library(sjPlot)
library(grid)
library(gridExtra)
library(tidyverse)

# 1- Load data --------
mag_counts <-
  read_tsv(file = "data/mags/mag_counts.tsv") %>%
  arrange(., by_group = mag_id) %>%
  column_to_rownames(var = "mag_id")

# Metadata
dmm <-
  read_tsv("results/tables/dmm.tsv")

metadata <-
  read_tsv(file = "data/kpi/metadata.tsv") %>%
  filter(animal_code %in% colnames(mag_counts)) %>%
  mutate(sampling_time = factor(sampling_time, levels = c("7",
                                                          "21",
                                                          "35"))) %>%
  mutate(trial = factor(trial, levels = c("CA", "CB", "CC"))) %>%
  left_join(dmm, by = 'animal_code') %>%
  mutate(enterotype = case_when(dmm == '2' ~ "distinct",
                                dmm != '2' ~ "standard")) %>%
  column_to_rownames(var = "animal_code")

rm(dmm, mag_counts)


# 2- BW temporal evolution --------
M_BW <- lme(log(chicken_body_weight) ~ log(age) * trial + breed + sex,
            random = ~1|pen,
            data = metadata)

plot_grid(plot_model(M_BW, type = "diag"))
summary(M_BW)
anova(M_BW)

get_model_data(M_BW, type = "emm",terms = c("age","trial"))

M_BW_7_emm <- emmeans(M_BW, "trial", by = "age", at = list(age = c(7)))
M_BW_7_emm_toplot <- data.frame(M_BW_7_emm)

day7_p <-
  ggplot(M_BW_7_emm_toplot, aes(x = trial,
                                y = exp(emmean),
                                color = trial)) +
  geom_point() +
  scale_color_manual(values = c('#56B4E9', '#009E73', '#D55E00')) +
  ggtitle("a) Day 7") +
  geom_errorbar(aes(ymin = exp(lower.CL),
                    ymax = exp(upper.CL)),
                width = .1,
                position = position_dodge(0.9)) +
  ylab("Chicken body weight (g)") +
  # scale_x_discrete(labels = c("Trial A","Trial B","Trial C")) +
  theme_minimal() +
  geom_text(aes(y = exp(upper.CL + 0.02),
                label = c("A","A","B"),
                size = 5)) +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 14))

M_BW_21_emm <- emmeans(M_BW, "trial", by = "age", at = list(age = c(21)))
M_BW_21_emm_toplot <- data.frame(M_BW_21_emm)

day21_p <-
  ggplot(M_BW_21_emm_toplot, aes(x = trial,
                                 y = exp(emmean),
                                 color = trial)) +
  geom_point() +
  scale_color_manual(values = c('#56B4E9', '#009E73', '#D55E00')) +
  ggtitle("b) Day 21") +
  geom_errorbar(aes(ymin = exp(lower.CL),
                    ymax = exp(upper.CL)),
                width = .1,
                position = position_dodge(0.9)) +
  ylab("Chicken body weight (g)") +
  scale_x_discrete(labels = c("Trial A","Trial B","Trial C")) +
  theme_minimal() +
  geom_text(aes(y = exp(upper.CL + 0.02),
                label = c("A","A","A"),
                size = 5)) +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 14))

M_BW_35_emm <- emmeans(M_BW, "trial", by = "age", at = list(age = c(35)))
M_BW_35_emm_toplot <- data.frame(M_BW_35_emm)

day35_p <-
  ggplot(M_BW_35_emm_toplot, aes(x = trial,
                                 y = exp(emmean),
                                 color = trial)) +
  geom_point() +
  scale_color_manual( values = c('#56B4E9', '#009E73', '#D55E00')) +
  ggtitle("c) Day 35") +
  geom_errorbar(aes(ymin = exp(lower.CL),
                    ymax = exp(upper.CL)),
                width = .1,
                position = position_dodge(0.9)) +
  ylab("Chicken body weight (g)") +
  scale_x_discrete(labels = c("Trial A","Trial B","Trial C")) +
  theme_minimal() +
  geom_text(aes(y = exp(upper.CL + 0.02),
                label = c("A","A","B"),
                size = 5)) +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 14))

myfun <- function(x) exp(x)
overall_p <-
  plot_model(M_BW,
             type = "emm",
             terms = c("age", "trial"),
             title = "d)",
             colors = c('#56B4E9', '#009E73', '#D55E00')) +
  theme_minimal() +
  ylab("Chicken body weight (g)") +
  xlab("Age (days)") +
  theme(plot.title = element_text(size = 14))

windows(w = 10)

grid.arrange(arrangeGrob(day7_p,day21_p, day35_p, ncol = 3, nrow = 1),
             arrangeGrob(overall_p, ncol = 1, nrow = 1),
             heights = c(2,3))
