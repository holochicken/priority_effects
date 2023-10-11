# Phylogenetic tree and relative expression barplots  - updated 04/09/2023

# Packages
library(vegan)
library(ggtree)
library(ggtreeExtra)
library(phyloseq)
library(treeio)
library(ggstar)
library(scales)
library(ggnewscale)
library(gridExtra)
library(tidytree)
library(tidyverse)

# 1- Load data ---------------------
# MAG codes
ena_assembly_to_mag_id <-
  read_tsv("data/mags/ena_assembly_to_mag_id.tsv")

# Taxonomy
taxonomy <-
  read_tsv("data/mags/taxonomy.tsv") %>%
  select(mag_id, order)

# MT table
mt_table_raw <-
  read_tsv("data/metatranscriptomics/MT_per1000000.tsv") %>%
  rename(mag_id = MAG) %>%
  left_join(ena_assembly_to_mag_id, by = "mag_id") %>%
  select(-mag_id) %>%
  rename(mag_id = mag_name) %>%
  left_join(taxonomy) %>%
  relocate(mag_id)

# Tree
tree <-
  read.tree("data/mags/tree.nwk")

# Taxonomy
taxonomy <-
  read_tsv("data/mags/taxonomy.tsv") %>%
  rename(label = 'mag_id')

# Joining taxonomy to tree
tax_tree <-
  tree %>%
  as_tibble() %>%
  left_join(taxonomy, by = 'label') %>%
  as.treedata()

tax_tree_tibble <-
  tax_tree %>%
  as_tibble() %>%
  filter(grepl('mag', label))


# 2- ordering MAGs ---------------------
label_order <- tree$tip.label

tax_order_correct <- c('Oscillospirales',
                       'Monoglobales',
                       'Monoglobales_A',
                       'TANB77',
                       'UBA1212',
                       'Lachnospirales',
                       'Christensenellales',
                       'Peptostreptococcales',
                       'Clostridiales',
                       'Veillonellales',
                       'Acidaminococcales',
                       'Selenomonadales',
                       'UBA4068',
                       'RF39',
                       'Erysipelotrichales',
                       'RFN20',
                       'ML615J-28',
                       'Acholeplasmatales',
                       'Lactobacillales',
                       'Bacillales',
                       'Coriobacteriales',
                       'Actinomycetales',
                       'Gastranaerophilales',
                       'Saccharimonadales',
                       'Bacteroidales',
                       'Flavobacteriales',
                       'Burkholderiales',
                       'Enterobacterales',
                       'RF32',
                       'Rs-D84',
                       'Desulfovibrionales',
                       'Deferribacterales',
                       'Victivallales',
                       'Opitutales',
                       'Verrucomicrobiales',
                       'Campylobacterales',
                       'Synergistales')

order_colors <- c(Methanomicrobiales = "#f2f2f2",
                  Methanobacteriales = "#bfbfbf",
                  Saccharimonadales = "#dacce3",
                  Campylobacterales = "#cdadca",
                  Synergistales = "#c08eb3",
                  Methanomassiliicoccales = "#777dae",
                  Victivallales = "#0066ae",
                  Opitutales = "#1c7ebc",
                  Verrucomicrobiales = "#4a96cc",
                  Desulfovibrionales = "#68b0dc",
                  Enterobacterales = "#60cfde",
                  RF32 = "#60dfd2",
                  Burkholderiales = "#5cdfb5",
                  'Rs-D84' = "#40df91",
                  Bacteroidales = "#4dc87c",
                  Flavobacteriales = "#88c88b",
                  Gastranaerophilales = "#92e09f",
                  Coriobacteriales = "#c9e0af",
                  Actinomycetales = "#d8e093",
                  Deferribacterales = "#dff77e",
                  RF39 = "#ecf76d",
                  Lactobacillales = "#feef68",
                  'ML615J-28' = "#fde671",
                  Erysipelotrichales = "#ffd366",
                  Acholeplasmatales = "#fdc151",
                  Bacillales = "#fc953d",
                  RFN20 = "#fd8035",
                  Oscillospirales = "#fd8854",
                  TANB77 = "#f4814d",
                  Christensenellales = "#c26340",
                  Lachnospirales = "#c28c5c",
                  Clostridiales = "#ce9360",
                  Monoglobales_A = "#ce734f",
                  Peptostreptococcales = "#c65631",
                  Monoglobales = "#b93725",
                  UBA1212 = "#b10f19",
                  UBA4068 = "#8b1222",
                  Selenomonadales = "#7a1a12",
                  Acidaminococcales = "#64121f",
                  Veillonellales = "#5a0c37")


###############################################################################
###############################################################################
## B06 - Organic anion biosynthesis
mt_table <-
  mt_table_raw %>%
  filter(Function == 'B06')

## Day 35 standard enterotype
standard_sum_d35 <-
  mt_table %>%
  filter(enterotype == 'standard') %>%
  filter(sampling_time == '35') %>%
  mutate(order = factor(order, levels = tax_order_correct)) %>%
  mutate(mag_id = factor(mag_id, levels = label_order))

s_d35_b06 <-
  ggplot(standard_sum_d35,
         aes(x = mean_value,
             y = interaction(mag_id, order),
             fill = order)) +
  geom_bar(alpha = 0.7,
           stat = "identity",
           position = "dodge",
           width = 5,
           size = 0.2) +
  scale_fill_manual(values = order_colors) +
  theme(
    legend.position = 'none',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = unit(c(1, 0.25, 0, 0.5), "cm"))

s_d35_b06_edited <-
  s_d35_b06 +
  scale_x_continuous(
    expand = c(0, 0),
    trans = 'reverse',
    limits = c(350, 0),
    breaks = c(350, 250, 150, 50, 0))


## Day 35 distinct enterotype
distinct_sum_d35 <-
  mt_table %>%
  filter(enterotype == 'distinct') %>%
  filter(sampling_time == '35') %>%
  mutate(order = factor(order, levels = tax_order_correct)) %>%
  mutate(mag_id = factor(mag_id, levels = label_order))

d_d35_b06 <-
  distinct_sum_d35 %>%
  ggplot() +
  geom_bar(aes(y = interaction(mag_id, order),
               x = mean_value,
               fill = order),
           alpha = 0.7,
           orientation = "y",
           stat = "identity",
           width = 5,
           size = 0.2) +
  scale_fill_manual(values = order_colors) +
  theme(
    legend.position = 'none',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = unit(c(1, 0.5, 0, 0.25), "cm"))

d_d35_b06_edited <-
  d_d35_b06 +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 350),
    breaks = c(0, 50, 150, 250, 350))



###############################################################################
###############################################################################
## D06- Nitrogen compound degradation
mt_table <-
  mt_table_raw %>%
  filter(Function == 'D06')

# standard d35 - Need to filter function**
standard_sum_d35 <-
  mt_table %>%
  filter(enterotype == 'standard') %>%
  filter(sampling_time == '35') %>%
  mutate(order = factor(order, levels = tax_order_correct)) %>%
  mutate(mag_id = factor(mag_id, levels = label_order))

s_d35_d06 <-
  ggplot(standard_sum_d35,
         aes(x = mean_value,
             y = interaction(mag_id, order),
             fill = order)) +
  geom_bar(alpha = 0.7,
           stat = "identity",
           position = "dodge",
           width = 7,
           size = 0.4) +
  scale_fill_manual(values = order_colors) +
  theme(
    legend.position = 'none',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = unit(c(1, 0.25, 0, 0.5), "cm"))

s_d35_d06_edited <-
  s_d35_d06 +
  scale_x_continuous(
    expand = c(0, 0),
    trans = 'reverse',
    limits = c(350, 0),
    breaks = c(350, 250, 150, 50, 0))

# Distinct enterotype
distinct_sum_d35 <-
  mt_table %>%
  filter(enterotype == 'distinct') %>%
  filter(sampling_time == '35') %>%
  mutate(order = factor(order, levels = tax_order_correct)) %>%
  mutate(mag_id = factor(mag_id, levels = label_order))

d_d35_d06 <-
  distinct_sum_d35 %>%
  ggplot() +
  geom_bar(aes(y = interaction(mag_id, order),
               x = mean_value,
               fill = order),
           alpha = 0.7,
           orientation = "y",
           stat = "identity",
           width = 7,
           size = 0.4) +
  scale_fill_manual(values = order_colors) +
  scale_x_continuous(trans = "reverse") +
  theme(
    legend.position = 'none',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = unit(c(1, 0.5, 0, 0.25), "cm"))

d_d35_d06_edited <-
  d_d35_d06 +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 350),
    breaks = c(0, 50, 150, 250, 350))


###############################################################################
###############################################################################
## B03 - Amino acid derivative biosynthesis
mt_table <-
  mt_table_raw %>%
  filter(Function == 'B03')

# standard d35 - Need to filter function**
standard_sum_d35 <-
  mt_table %>%
  filter(enterotype == 'standard') %>%
  filter(sampling_time == '35') %>%
  mutate(order = factor(order, levels = tax_order_correct)) %>%
  mutate(mag_id = factor(mag_id, levels = label_order))

s_d35_b03 <-
  ggplot(standard_sum_d35,
         aes(x = mean_value,
             y = interaction(mag_id, order),
             fill = order)
  ) +
  geom_bar(alpha = 0.7,
           stat = "identity",
           position = "dodge",
           width = 7,
           size = 0.4) +
  scale_fill_manual(values = order_colors) +
  theme(
    legend.position = 'none',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = unit(c(1, 0.25, 0, 0.5), "cm"))

s_d35_b03_edited <-
  s_d35_b03 +
  scale_x_continuous(
    expand = c(0, 0),
    trans = 'reverse',
    limits = c(350, 0),
    breaks = c(350, 250, 150, 50, 0))

# Distinct enterotype
distinct_sum_d35 <-
  mt_table %>%
  filter(enterotype == 'distinct') %>%
  filter(sampling_time == '35') %>%
  mutate(order = factor(order, levels = tax_order_correct)) %>%
  mutate(mag_id = factor(mag_id, levels = label_order))

d_d35_b03 <-
  distinct_sum_d35 %>%
  ggplot() +
  geom_bar(aes(y = interaction(mag_id, order),
               x = mean_value,
               fill = order),
           alpha = 0.7,
           orientation = "y",
           stat = "identity",
           width = 7,
           size = 0.4) +
  scale_fill_manual(values = order_colors) +
  scale_x_continuous(trans = "reverse") +
  theme(
    legend.position = 'none',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = unit(c(1, 0.5, 0, 0.25), "cm"))

d_d35_b03_edited <-
  d_d35_b03 +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 350),
    breaks = c(0, 50, 150, 250, 350))


###############################################################################
###############################################################################
## D01- Lipid degradation
mt_table <-
  mt_table_raw %>%
  filter(Function == 'D01')

# standard d35 - Need to filter function**
standard_sum_d35 <-
  mt_table %>%
  filter(enterotype == 'standard') %>%
  filter(sampling_time == '35') %>%
  mutate(order = factor(order, levels = tax_order_correct)) %>%
  mutate(mag_id = factor(mag_id, levels = label_order))

s_d35_d01 <-
  ggplot(standard_sum_d35,
         aes(x = mean_value,
             y = interaction(mag_id, order),
             fill = order)
  ) +
  geom_bar(alpha = 0.7,
           stat = "identity",
           position = "dodge",
           width = 7,
           size = 0.4) +
  scale_fill_manual(values = order_colors) +
  xlab('Relative expression') +
  theme(
    legend.position = 'none',
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = unit(c(1, 0.25, 0, 0.5), "cm"))

s_d35_d01_edited <-
  s_d35_d01 +
  scale_x_continuous(
    expand = c(0, 0),
    trans = 'reverse',
    limits = c(350, 0),
    breaks = c(350, 250, 150, 50, 0))

# Distinct enterotype
distinct_sum_d35 <-
  mt_table %>%
  filter(enterotype == 'distinct') %>%
  filter(sampling_time == '35') %>%
  mutate(order = factor(order, levels = tax_order_correct)) %>%
  mutate(mag_id = factor(mag_id, levels = label_order))

d_d35_d01 <-
  distinct_sum_d35 %>%
  ggplot() +
  geom_bar(aes(y = interaction(mag_id, order),
               x = mean_value,
               fill = order),
           alpha = 0.7,
           orientation = "y",
           stat = "identity",
           width = 7,
           size = 0.4) +
  scale_fill_manual(values = order_colors) +
  xlab('Relative expression') +
  theme(
    legend.position = 'none',
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = unit(c(1, 0.5, 0, 0.25), "cm"))

d_d35_d01_edited <-
  d_d35_d01 +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 350),
    breaks = c(0, 50, 150, 250, 350))


grid.arrange(s_d35_b06_edited, d_d35_b06_edited,
             s_d35_d06_edited, d_d35_d06_edited,
             s_d35_b03_edited, d_d35_b03_edited,
             s_d35_d01_edited, d_d35_d01_edited,
             ncol = 2)

dev.off()
