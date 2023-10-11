# Temporal changes and sources of variation on the microbiota - updated 04/09/2023

# packages
library(vegan)
library(microbiome)
library(ggrepel)
library(tidyverse)

# 1- Load data  -----------------
# Abundance table
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
  left_join(dmm, by = 'animal_code') %>%
  mutate(enterotype = case_when(dmm == '2' ~ "distinct",
                                dmm != '2' ~ "standard")) %>%
  column_to_rownames(var = "animal_code")

# MAG stats
stats <-
  read_tsv("data/mags/stats.tsv") %>%
  mutate(correction_factor = median(mag_length) / mag_length)

# MAG taxonomy
taxonomy <-
  read_tsv("data/mags/taxonomy.tsv")


# 2- Standardization and correction  -----------------
# Correction
mag_weighted <-
  round(sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0) %>%
  t() %>%
  as.data.frame()

# # metadata summary
design <-
  metadata %>%
  dplyr::select('trial', 'pen', 'age', 'sampling_time',
         'breed', 'sex', 'treatment', 'enterotype')

# Transformation
mag_clr <-
  microbiome::transform(phyloseq::otu_table(mag_weighted,
                                            taxa_are_rows = F),
                        transform = "clr")


## 3- Distance-based RDA  -----------------
set.seed(4)
rda_table <- rda(mag_clr ~ trial * age, data = design)

# Anova
summary(rda_table)
anova(rda_table, by = "term")
RsquareAdj(rda_table)

rda_scores <- data.frame(scores(rda_table, display = "wa"),
                           pen = design$pen,
                           trial = design$trial,
                           time = design$sampling_time,
                           enterotype = design$enterotype)

rda_scores$animal_code <- row.names(rda_scores)

rda_scores$enterotype <-  as.factor(rda_scores$enterotype)

# For age arrows
db_rda_constraining_scores <- data.frame(scores(rda_table,display = "bp"))
db_rda_constraining_scores_age <- db_rda_constraining_scores[4,]

# Detach packages causing error
detach(package:microbiome, unload = TRUE)
detach(package:vegan, unload = TRUE)

# Plot figure
pdf(file = "results/figures/fig_1c_rda.pdf", width = 7, height = 5)

ggplot(rda_scores) +
  geom_point(aes(x = RDA1,
                 y = -RDA2,
                 shape = time,
                 colour = trial,
                 fill = enterotype),
             size = 3,
             stroke = 0.9) +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_color_manual(values = c('#56B4E9', '#009E73', '#D55E00')) +
  scale_fill_manual(values = alpha(c("#D55E00", "white"), 0.7)) +
  ylab("RDA2 (36%)") +
  xlab("RDA1 (55%)") +
  theme_bw() +
  theme(panel.border = element_blank())

dev.off()


# 4- Filter d7 samples -----------------
metadata_d7 <-
  metadata %>%
  filter(sampling_time == '7') %>%
  rownames_to_column("animal_code")

# count table
mag_table_d7 <-
  mag_clr %>%
  data.frame() %>%
  filter(rownames(.) %in% metadata_d7$animal_code)

mean(row.names(mag_table_d7)) == metadata_d7$animal_code

# Calculate community MCI
library(distillR)

total <-
  decostand(mag_weighted, "total") %>%
  t() %>%
  as.data.frame()

gifts <-
  read_tsv('results/tables/gifts_elements.tsv') %>%
  column_to_rownames(var = 'mag_id')

gifts_community <-
  to.community(gifts,
               sweep(total, 2, colSums(total), FUN = "/"),
               GIFT_db1)

gifts_community_d7 <-
  gifts_community %>%
  as.data.frame() %>%
  filter(rownames(.) %in% metadata_d7$animal_code)


# 5- PERMANOVA ---------------------
library(vegan)
perm <- how(nperm = 9999)
setBlocks(perm) <- with(metadata_d7, pen)

adonis2(mag_table_d7 ~ trial + sex + treatment + breed + age + enterotype,
        method = "euclidean",
        data = metadata_d7)

adonis2(gifts_community_d7 ~ trial + sex + treatment + breed + age + enterotype,
        method = "euclidean",
        data = metadata_d7)

rm(list = ls())
