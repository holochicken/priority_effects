# Phylogenetic tree and relative abundance barplots  - updated 04/09/2023


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

# 1- Load and tidy data ---------------------
# 1.1- Relative abundance table
mag_counts <-
  read_tsv("data/mags/mag_counts.tsv") %>%
  arrange(., by_group = mag_id) %>%
  column_to_rownames(var = 'mag_id')

# # 1.2- Metadata
dmm <- read_tsv("results/tables/dmm.tsv")

metadata <-
  read_tsv(file = "data/kpi/metadata.tsv")  %>%
  left_join(dmm, by = 'animal_code') %>%
  mutate(enterotype = case_when(dmm == '2' ~ "distinct",
                                dmm != '2' ~ "standard"))

# # 1.3- MAG stats
stats <-
  read_tsv("data/mags/stats.tsv") %>%
  mutate(correction_factor = median(mag_length) / mag_length)

# # 1.4- Standardisation
# # By MAG genome length
mag_weighted <-
  round(sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0) %>%
  t() %>%
  data.frame()

# # By sequencing depth
total <- decostand(mag_weighted, 'total')
rowSums(total) # check if all values are 1

# 1.6- Taxonomy
taxonomy <-
  read_tsv("data/mags/taxonomy.tsv")

# Calculate rel counts per mag and per sampling day
mag_counts_day <-
  total %>%
  rownames_to_column('animal_code') %>%
  left_join(metadata %>% select(animal_code,
                                sampling_time,
                                enterotype),
            by = 'animal_code') %>%
  pivot_longer(!c(animal_code, sampling_time, enterotype),
               names_to = 'mag_id', values_to = 'values') %>%
  group_by(sampling_time, enterotype, mag_id) %>%
  summarise(total = colMeans(across(where(is.numeric)))) %>%
  arrange(., by_group = mag_id) %>%
  left_join(taxonomy %>% select(mag_id, order), by = 'mag_id') %>%
  mutate(
    sampling_time = factor(sampling_time, levels = c('7',
                                                     '21',
                                                     '35'))) %>%
  mutate(
    enterotype = factor(enterotype, levels = c('standard',
                                               'distinct'))) %>%
  relocate(mag_id)

# 1.5- Phylogenetic tree
tree <-
  read.tree("data/mags/tree.nwk")

taxonomy <-
  taxonomy %>%
  rename(label = 'mag_id')

# 1.6- Joining taxonomy to tree
tax_tree <-
  tree %>%
  as_tibble() %>%
  left_join(taxonomy, by = 'label') %>%
  as.treedata()

rm(dmm, mag_counts, mag_weighted, stats)

# 2- Vertical plot ---------------------
# Define colors
phylum_colors <- c(Halobacteriota = "#f2f2f2",
                   Methanobacteriota = "#bfbfbf",
                   Patescibacteria = "#dacce3",
                   Campylobacterota = "#cdadca",
                   Synergistota = "#c08eb3",
                   Thermoplasmatota = "#777dae",
                   Verrucomicrobiota = "#0066ae",
                   Desulfobacterota = "#68b0dc",
                   Proteobacteria = "#60cfde",
                   Bacteroidota = "#4dc87c",
                   Cyanobacteria = "#92e09f",
                   Actinobacteriota = "#c9e0af",
                   Deferribacterota = "#dff77e",
                   Firmicutes = "#ffd366",
                   Firmicutes_A = "#fd8854",
                   Firmicutes_B = "#8b1222",
                   Firmicutes_C = "#5a0c37")

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


# a) Draw vertical tree
p <-
  ggtree(tax_tree, branch.length = 'none', color = '#36454F') +
  theme(legend.position = 'none')
p

# b) caculate MRCA for phylum and illustrate them
phylum_ids <-
  taxonomy %>%
  group_by(phylum) %>%
  reframe(phylum_node_id = MRCA(tree, label)) %>%
  rename(phylum_id = 'phylum')

p1 <-
  p +
  geom_highlight(data = phylum_ids,
                 mapping = aes(node = phylum_node_id,
                               fill = phylum_id),
                 show.legend = FALSE,
                 align = 'both',
                 alpha = 0.4,
                 to.bottom = TRUE) +
  scale_fill_manual(values = phylum_colors)
p1

# c) Add order column
order_ids <-
  taxonomy %>%
  group_by(order) %>%
  reframe(order_node_id = MRCA(tree, label)) %>%
  rename(order_id = 'order')

p2 <-
  p1 +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = mag_id, fill = order),
    offset = 0.1,
    width = 5,
    grid.params = list(color = NA)) +
  scale_fill_manual(values = order_colors) +
  geom_cladelab(
    data = order_ids,
    mapping = aes(node = order_node_id, label = order_id),
    barcolour = NA,
    align = TRUE,
    barsize = 1,
    offset = 6,
    fontsize = 1.5)
p2

# plot tree
pdf(file = "results/figures/fig_1c_phylo_tree.pdf", width = 5, height = 10)
p2
dev.off()


# 3- Relative abundance barplots
# Defining log transformations for axis
reverse_log_trans <-
  function() {
    trans_new(
    "reverse_log",
    transform = function(x) -log1p(x),
    inverse = function(x) exp(-x) - 1
  )
}

#
log_trans <-
  function() {
  trans_new(
    "reverse_log",
    transform = function(x) log1p(x),
    inverse = function(x) exp(x) - 1
  )
}

# Organising taxonomy orders
tax_tree_tibble <-
  tax_tree %>%
  as_tibble() %>%
  filter(grepl('mag', label))

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

###############################################################################
###############################################################################
## Day 7
standard_sum_d7 <-
  mag_counts_day %>%
  filter(enterotype == 'standard') %>%
  filter(sampling_time == '7') %>%
  mutate(order = factor(order, levels = tax_order_correct)) %>%
  mutate(mag_id = factor(mag_id, levels = label_order))

standard_sum_d7$total <- standard_sum_d7$total * 100

s_d7 <-
  ggplot(standard_sum_d7,
         aes(x = total,
             y = interaction(mag_id, desc(order)),
             fill = order)) +
  geom_bar(alpha = 0.7,
           # orientation = "y",
           stat = "identity",
           position = "dodge",
           width = 5,
           size = 0.2) +
  scale_fill_manual(values = order_colors) +
  # xlim(0,0.3) +
  # scale_x_log10() +
  # coord_flip() +
  # scale_x_reverse() +
  # theme_bw() +
  xlab('Relative abundance') +
  ggtitle('Standard - d7') +
  theme(
    legend.position = 'none',
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank())
s_d7

# editing x axis
s_d7_edited_manually <- s_d7 +
  scale_x_continuous(trans = reverse_log_trans(),
                     limits = c(30, 0),
                     breaks = c(30, 10, 2, 0),
                     labels = c('0.3', '0.1', '0.02', '0'))
s_d7_edited_manually

# Distinct enterotype
distinct_sum_d7 <-
  mag_counts_day %>%
  filter(enterotype == 'distinct') %>%
  filter(sampling_time == '7') %>%
  mutate(order = factor(order, levels = tax_order_correct)) %>%
  mutate(mag_id = factor(mag_id, levels = label_order))

distinct_sum_d7$total <- distinct_sum_d7$total * 100

d_d7 <-
  distinct_sum_d7 %>%
  ggplot() +
  geom_bar(aes(y = interaction(mag_id, desc(order)), x = total, fill = order),
           alpha = 0.7,
           orientation = "y",
           stat = "identity",
           width = 5,
           size = 0.2) +
  scale_fill_manual(values = order_colors) +
  scale_x_continuous(trans = "reverse") +
  xlab('Relative abundance') +
  ggtitle('Distinct - d7') +
  theme(
    legend.position = 'none',
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank())

d_d7_edited_manually <-
  d_d7 +
  scale_x_continuous(trans = log_trans(),
                     limits = c(0, 30),
                     breaks = c(0,2,10,30),
                     labels = c('0', '0.02', '0.1', '0.3')
  )
d_d7_edited_manually

pdf(file = "results/figures/fig_1c_phylo_rel_abu.pdf", width = 7, height = 10)

grid.arrange(s_d7_edited_manually, d_d7_edited_manually, ncol = 2)

dev.off()

###############################################################################
###############################################################################
## Day 21 standard enterotype
standard_sum_d21 <-
  mag_counts_day %>%
  filter(enterotype == 'standard') %>%
  filter(sampling_time == '21') %>%
  mutate(order = factor(order, levels = tax_order_correct)) %>%
  mutate(mag_id = factor(mag_id, levels = label_order))

standard_sum_d21$total <- standard_sum_d21$total * 100

s_d21 <-
  ggplot(standard_sum_d21, aes(x = total, y = interaction(mag_id, desc(order)), fill = order)) +
  geom_bar(alpha = 0.7,
           # orientation = "y",
           stat = "identity",
           position = "dodge",
           width = 5,
           size = 0.2) +
  scale_fill_manual(values = order_colors) +
  xlab('Relative abundance') +
  ggtitle('Standard - d21') +
  theme(
    legend.position = 'none',
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank())
s_d21

s_d21_edited_manually <- s_d21 +
  scale_x_continuous(trans = reverse_log_trans(),
                     limits = c(30, 0),
                     breaks = c(30, 10, 2, 0),
                     labels = c('0.3', '0.1', '0.02', '0'))
s_d21_edited_manually


## Day 21 distinct enterotype
distinct_sum_d21 <-
  mag_counts_day %>%
  filter(enterotype == 'distinct') %>%
  filter(sampling_time == '21') %>%
  mutate(order = factor(order, levels = tax_order_correct)) %>%
  mutate(mag_id = factor(mag_id, levels = label_order))

distinct_sum_d21$total <- distinct_sum_d21$total * 100

d_d21 <-
  distinct_sum_d21 %>%
  ggplot() +
  geom_bar(aes(y = interaction(mag_id, desc(order)), x = total, fill = order),
           alpha = 0.7,
           orientation = "y",
           stat = "identity",
           width = 5,
           size = 0.2) +
  scale_fill_manual(values = order_colors) +
  scale_x_continuous(trans = "reverse") +
  xlab('Relative abundance') +
  ggtitle('Distinct - d21') +
  theme(
    legend.position = 'none',
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank())

d_d21_edited_manually <-
  d_d21 +
  scale_x_continuous(trans = log_trans(),
                     limits = c(0, 30),
                     breaks = c(0,2,10,30),
                     labels = c('0', '0.02', '0.1', '0.3')
  )

d_d21_edited_manually

###############################################################################
###############################################################################
## Day 35 standard enterotype
standard_sum_d35 <-
  mag_counts_day %>%
  filter(enterotype == 'standard') %>%
  filter(sampling_time == '35') %>%
  mutate(order = factor(order, levels = tax_order_correct)) %>%
  mutate(mag_id = factor(mag_id, levels = label_order))

standard_sum_d35$total <- standard_sum_d35$total * 100

s_d35 <-
  ggplot(standard_sum_d35,
         aes(x = total, y = interaction(mag_id, desc(order)),
             fill = order)) +
  geom_bar(alpha = 0.7,
           stat = "identity",
           position = "dodge",
           width = 5,
           size = 0.2) +
  scale_fill_manual(values = order_colors) +
  xlab('Relative abundance') +
  ggtitle('Standard - d35') +
  theme(
    legend.position = 'none',
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank())
s_d35

s_d35_edited_manually <-
  s_d35 +
  scale_x_continuous(trans = reverse_log_trans(),
                     limits = c(30, 0),
                     breaks = c(30, 10, 2, 0),
                     labels = c('0.3', '0.1', '0.02', '0'))
s_d35_edited_manually


## Day 35 distinct enterotype
distinct_sum_d35 <-
  mag_counts_day %>%
  filter(enterotype == 'distinct') %>%
  filter(sampling_time == '35') %>%
  mutate(order = factor(order, levels = tax_order_correct)) %>%
  mutate(mag_id = factor(mag_id, levels = label_order))

distinct_sum_d35$total <- distinct_sum_d35$total * 100

d_d35 <-
  distinct_sum_d35 %>%
  ggplot() +
  geom_bar(aes(y = interaction(mag_id, desc(order)), x = total, fill = order),
           alpha = 0.7,
           orientation = "y",
           stat = "identity",
           width = 5,
           size = 0.2) +
  scale_fill_manual(values = order_colors) +
  scale_x_continuous(trans = "reverse") +
  xlab('Relative abundance') +
  ggtitle('Distinct - d35') +
  theme(
    legend.position = 'none',
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank())

d_d35_edited_manually <-
  d_d35 +
  scale_x_continuous(trans = log_trans(),
                     limits = c(0, 30),
                     breaks = c(0,2,10,30),
                     labels = c('0', '0.02', '0.1', '0.3')
  )

d_d35_edited_manually


# Plot figures
pdf(file = "results/figures/fig_supp_s2_phylo_rel_abu.pdf", width = 12, height = 10)

grid.arrange(s_d7_edited_manually, d_d7_edited_manually,
             s_d21_edited_manually, d_d21_edited_manually,
             s_d35_edited_manually, d_d35_edited_manually,
             ncol = 6)

dev.off()


rm(list = ls())
