# Plotting code for "Alternative MAPK drivers in KRAS WT pancreatic cancer - Letter"
# Visualize cancer effect sizes and epistatic effects in PAAD
# Last updated: April 26, 2024

# load libraries
library(cancereffectsizeR) # using version 2.8.0
library(tidyverse)
library(ces.refset.hg19) # using version 1.2.2
library(TCGAbiolinks)
library(ggplot2)
library(data.table)
library(dplyr)
options(scipen = 100, digits = 4) # for plot axis scaling

# load the complete CESA as generated in "paad_analysis.R"
paad_cesa <- load_cesa('paad_cesa.rds')

# extract selection results from CESAnalysis and take top variants for visualization
results <- paad_cesa$selection$kras
results$selection_intensity = log10(res$selection_intensity)
gene_names <- c('KRAS', 'BRAF', 'NRAS', 'IDH1', 'SMAD4', 'TP53', 'CTNNB1', 'CDKN2A.p16INK4a', 'GNAS', 'ARID2', 'PTEN', 
                'NF2', 'WRN', 'PMS2', 'NTRK1', 'ROS1', 'PIK3CA', 'MYC') # frequently mutated genes in PAAD described in the original Singh et al. article

# create and save plot
paad_plot <- plot_effects(results[gene %in% gene_names],
                          group_by = "gene",
                          color_by = "#df8f44",
                          legend_size_name = "Variant Prevalence") +
  theme_classic() + 
  theme(legend.position = c(0.8, 0.2), legend.direction = "horizontal", legend.title = element_text(hjust = 0.5)) +
  labs(x = "Scaled selection coefficient", y = "Gene") + 
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000))+ theme(text = element_text(family = "Verdana"),
                                                                             axis.text = element_text(size = 14), 
                                                                             axis.title = element_text(size = 16))

paad_si_plot <- ggarrange(paad_plot, nrow = 1, ncol = 1, labels = c('A'), font.label = list(size = 16, face = "bold", family = "Verdana"))
ggsave("paad_si_plot.jpeg", paad_si_plot, width = 11.5, height = 6)

# perform pairwise epistasis analysis
# use baseline KRAS selection from default model

# create list of all other genes mutated in combined paad dataset besides KRAS
other_genes = paad_cesa$gene_rates$gene[paad_cesa$gene_rates$gene != 'KRAS']
gene_pairs = lapply(other_genes, c, 'KRAS')

# calculate pairwise epistatic selection intensitites for each gene paired with KRAS
paad_cesa <- ces_gene_epistasis(paad_cesa, genes = gene_pairs, variants = "recurrent", run_name = "kras_epi")

# subset results for know driver genes and IDH1, which was found to have particularly significant suppressive effects on the selection of KRAS mutations
gene_ep_results <- subset(paad_cesa_complete$epistasis$kras_epi, 
                          variant_A %in% c('BRAF', 'NRAS', 'IDH1', 'SMAD4', 'TP53', 'CTNNB1', 'CDKN2A.p16INK4a'))

# calculate ratios of selection, the mean fold-change of selection on mutant background over wildtype background
gene_ep_results$change_in_vB <-  gene_ep_results$ces_B_on_A/gene_ep_results$ces_B0
gene_ep_results$change_in_vA <- gene_ep_results$ces_A_on_B/gene_ep_results$ces_A0

# generate a plot of epistasis results where PAAD driver genes are mutated in the background of mutated and wildtype KRAS
# create breaks based on quintiles for epistatic pair prevalence
breaks <- round(quantile(gene_ep_results$nAB, probs = seq(0, 1, by = 0.25)))
gene_ep_results_kras_background <- gene_ep_results

# change labels to just mutated gene
gene_ep_results_kras_background = 
  gene_ep_results_kras_background %>%
  mutate(label = paste0(variant_A))

axis_labels = gene_ep_results_kras_background %>%
  arrange((change_in_vA)) %>%
  pull(label)

# generate plot
kras_background_plot <- gene_ep_results_kras_background %>%
  mutate(label = factor(label, levels = axis_labels)) %>%
  ggplot(aes(x = label, y = ces_A_on_B)) + 
  geom_errorbar(aes(ymin = ci_low_95_ces_A0, ymax = ci_high_95_ces_A0), width = 0.1, color = '#80796B') + 
  geom_point(aes(y = ces_A0, size = nA0), fill = "#df8f44", color = "black", shape = 21) +
  geom_errorbar(aes(ymin = ci_low_95_ces_A_on_B, ymax = ci_high_95_ces_A_on_B), width = 0.1, color = '#80796B') + 
  geom_point(aes(size = nAB), fill = "#373e55", color = "black", shape = 21) +
  coord_flip() + 
  theme_classic() +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = 3)) +
  labs(x = "Gene", y = "Scaled selection coefficient") + 
  scale_size_continuous(breaks = breaks, labels = breaks) + 
  scale_y_log10() +
  guides(size = guide_legend(title.position = "top", title.hjust = 0.5)) +
  theme(legend.position = c(0.2, 0.7), legend.direction = "horizontal", legend.box.background = element_rect(color = "black"), legend.text = element_text(size = 10), legend.title = element_text(size = 11)) +
  labs(size = "Mutated gene on mutated KRAS") + 
  theme(text = element_text(family = "Verdana"),
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16))

kras_background_plot_save <- ggarrange(kras_background_plot, nrow = 1, ncol = 1, labels = c('B'), font.label = list(size = 16, face = "bold", family = "Verdana"))
ggsave("kras_background_plot.jpeg", kras_background_plot_save, width = 11.5, height = 6)

# generate a plot of epistasis results where KRAS is mutated in the background of mutated and wildtype PAAD driver genes
gene_ep_results_kras_mut <- gene_ep_results

# change labels to just mutated gene
gene_ep_results_kras_mut = 
  gene_ep_results_kras_mut %>%
  mutate(label2 = paste0(variant_A))

axis_labels2 = gene_ep_results_kras_mut %>%
  arrange((change_in_vB)) %>%
  pull(label2)

# fix NA CIs by flooring them
gene_ep_results_kras_wt$ces_B_on_A[1] = 1
gene_ep_results_kras_wt$ces_B_on_A[7] = 1

# generate plot
kras_mut_plot <- gene_ep_results_kras_mut %>% 
  mutate(label2 = factor(label2, levels = axis_labels2)) %>%
  ggplot(aes(x=label2, y=ces_B_on_A)) + 
  geom_errorbar(aes(ymin = ci_low_95_ces_B0, ymax = ci_high_95_ces_B0), width = 0.1, color = '#80796B') + #, position = position_nudge(gene_ep_results$nudge_dist)) + 
  geom_point(aes(y=ces_B0, size = nB0), fill = "#df8f44", color = "black", shape = 21) + #, position = position_nudge(gene_ep_results$nudge_dist)) +
  geom_errorbar(aes(ymin = ci_low_95_ces_B_on_A, ymax = ci_high_95_ces_B_on_A), width = 0.1, color = '#80796B') + #, position = position_nudge(gene_ep_results$nudge_dist)) + 
  geom_point(aes(y=ces_B_on_A, size = nAB), fill = "#373e55", color = "black", shape = 21) + #, position = position_nudge(gene_ep_results$nudge_dist)) +
  coord_flip() + 
  scale_size_continuous(breaks = breaks, labels = breaks) + 
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000)) + 
  theme_classic() +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = 3)) +
  labs(x = "Gene", y = "Scaled selection coefficient") +
  guides(size = guide_legend(title.position = "top", title.hjust = 0.5)) +
  theme(legend.direction = "horizontal", legend.box.background = element_rect(color = "black"), legend.text = element_text(size = 10), legend.title = element_text(size = 11)) +
  labs(size = "Mutated KRAS on mutated gene") + 
  theme(text = element_text(family = "Verdana"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

kras_mut_plot_save <- ggarrange(kras_mut_plot, nrow = 1, ncol = 1, labels = c('C'), font.label = list(size = 16, face = "bold", family = "Verdana"))
ggsave("kras_mut_plot.jpeg", kras_mut_plot_save, width = 11.5, height = 6)
