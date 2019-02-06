#Written for the analysis of a range of inductions in the episomal library
#DNA_BC: DNA
#R0A_BC and R0B_BC: RNA at 0 µM Forsk replicate A and B
#R2#A_BC and R2#B_BC: RNA at 2^# µM Forsk replicate A and B, includes negative 
#and positive #'s

#tested concentrations: 0, 2^-5, 2^-4, 2^-3, 2^-2, 2^-1, 2^0, 2^2 µM Forsk

#All figures require sections "Load index and bcmap files" to 'Median BC 
#expression"

#Establish workspace------------------------------------------------------------

#import necessary libraries

library(tidyverse)
library(stringr)
library(viridis)
library(cowplot)
library(ggExtra)
library(modelr)
library(lazyeval)
library(splines)
library(broom)
library(GGally)
library(lemon)
library(devtools)
library(updateR)
library(ggsignif)

#General figure customizations

cbPalette7 <- c('#440154FF', '#39568CFF', '#287D8EFF', '#20A387FF', '#73D055FF',
                '#B8DE29FF', '#FDE725FF')

cbPalette7_grad_light <- c('white', '#FDE725FF', '#B8DE29FF', '#55C667FF', 
                           '#1F968BFF', '#39568CFF', '#482677FF')

spacing_5_20_palette <- c('gray20', 'dodgerblue3', 'indianred2', '#55C667FF')

figurefont_theme <- theme(text = element_text(size = 8)) +
  theme(axis.title = element_text(size = 8)) +
  theme(legend.title = element_text(size = 9)) +
  theme(axis.text = element_text(size = 8))

#General functions

var_log10 <- function(df) {
  log_ratio_df <- df %>% 
    mutate_if(is.double, funs(log10(.)))
  return(log_ratio_df)
}

#Load index and bcmap files-----------------------------------------------------

bc_DNA <- read_tsv('BCreads_txts/DNA_BC.txt')
bc_R0A <- read_tsv('BCreads_txts/R0A_BC.txt')
bc_R0B <- read_tsv('BCreads_txts/R0B_BC.txt')
bc_R2_5A <- read_tsv('BCreads_txts/R2-5A_BC.txt')
bc_R2_5B <- read_tsv('BCreads_txts/R2-5B_BC.txt')
bc_R2_4A <- read_tsv('BCreads_txts/R2-4A_BC.txt')
bc_R2_4B <- read_tsv('BCreads_txts/R2-4B_BC.txt')
bc_R2_3A <- read_tsv('BCreads_txts/R2-3A_BC.txt')
bc_R2_3B <- read_tsv('BCreads_txts/R2-3B_BC.txt')
bc_R2_2A <- read_tsv('BCreads_txts/R2-2A_BC.txt')
bc_R2_2B <- read_tsv('BCreads_txts/R2-2B_BC.txt')
bc_R2_1A <- read_tsv('BCreads_txts/R2-1A_BC.txt')
bc_R2_1B <- read_tsv('BCreads_txts/R2-1B_BC.txt')
bc_R20A <- read_tsv('BCreads_txts/R20A_BC.txt')
bc_R20B <- read_tsv('BCreads_txts/R20B_BC.txt')
bc_R22A <- read_tsv('BCreads_txts/R22A_BC.txt')
bc_R22B <- read_tsv('BCreads_txts/R22B_BC.txt')


#Load barcode mapping table, sequences are rcomp due to sequencing format

barcode_map <- read_tsv('../../BCMap/uniqueSP2345.txt', 
                        col_names = c(
                          'fluff', 'barcode', 'name', 'most_common'), 
                        skip = 1) %>%
  select(-fluff) %>%
  mutate(subpool = ifelse(startsWith(name, 'subpool'), 
                          substr(name, 1, 8), 
                          'control')) 


#pick out controls, SP3, and SP5 in the bcmap that were used in this assay

SP3_SP5_map <- barcode_map %>%
  mutate(subpool = ifelse(
    startsWith(name, 'subpool'), 
    substr(name, 1, 8), 
    'control')) %>%
  filter(subpool != 'subpool2') %>%
  filter(subpool != 'subpool4')


#Join reads to bcmap------------------------------------------------------------

#Join BC reads to BC mapping, keeping the reads only appearing in barcode 
#mapping and replacing na with 0 reads. Determine normalized reads per million

bc_map_join_bc <- function(df1, df2) {
  df2 <- df2 %>%
    mutate(normalized = as.numeric((num_reads * 1000000) / (sum(num_reads))))
  keep_bc <- left_join(df1, df2, by = 'barcode') %>%
    mutate(normalized = if_else(is.na(normalized), 
                                0, 
                                normalized)) %>%
    mutate(num_reads = if_else(is.na(num_reads), 
                               as.integer(0), 
                               num_reads))
  return(keep_bc)
}

bc_join_DNA <- bc_map_join_bc(SP3_SP5_map, bc_DNA)
bc_join_R0A <- bc_map_join_bc(SP3_SP5_map, bc_R0A)
bc_join_R0B <- bc_map_join_bc(SP3_SP5_map, bc_R0B)
bc_join_R2_5A <- bc_map_join_bc(SP3_SP5_map, bc_R2_5A)
bc_join_R2_5B <- bc_map_join_bc(SP3_SP5_map, bc_R2_5B)
bc_join_R2_4A <- bc_map_join_bc(SP3_SP5_map, bc_R2_4A)
bc_join_R2_4B <- bc_map_join_bc(SP3_SP5_map, bc_R2_4B)
bc_join_R2_3A <- bc_map_join_bc(SP3_SP5_map, bc_R2_3A)
bc_join_R2_3B <- bc_map_join_bc(SP3_SP5_map, bc_R2_3B)
bc_join_R2_2A <- bc_map_join_bc(SP3_SP5_map, bc_R2_2A)
bc_join_R2_2B <- bc_map_join_bc(SP3_SP5_map, bc_R2_2B)
bc_join_R2_1A <- bc_map_join_bc(SP3_SP5_map, bc_R2_1A)
bc_join_R2_1B <- bc_map_join_bc(SP3_SP5_map, bc_R2_1B)
bc_join_R20A <- bc_map_join_bc(SP3_SP5_map, bc_R20A)
bc_join_R20B <- bc_map_join_bc(SP3_SP5_map, bc_R20B)
bc_join_R22A <- bc_map_join_bc(SP3_SP5_map, bc_R22A)
bc_join_R22B <- bc_map_join_bc(SP3_SP5_map, bc_R22B)

#Median BC expression-----------------------------------------------------------

#Join DNA BC reads > 6 reads and join to RNA. Take ratio of RNA/DNA normalized
#reads per million

dna7_join_rna_rep <- function(df1, df2) {
  filter_DNA <- filter(df1, num_reads > 6)
  DNA_RNA_join <- left_join(filter_DNA, df2,
                            by = c("barcode", "name", "subpool", 
                                   "most_common"), 
                            suffix = c('_DNA', '_RNA')) %>%
    mutate(ratio = normalized_RNA/normalized_DNA)
  return(DNA_RNA_join)
}

bc_DNA_RNA_0A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R0A)
bc_DNA_RNA_0B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R0B)
bc_DNA_RNA_2_5A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_5A)
bc_DNA_RNA_2_5B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_5B)
bc_DNA_RNA_2_4A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_4A)
bc_DNA_RNA_2_4B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_4B)
bc_DNA_RNA_2_3A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_3A)
bc_DNA_RNA_2_3B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_3B)
bc_DNA_RNA_2_2A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_2A)
bc_DNA_RNA_2_2B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_2B)
bc_DNA_RNA_2_1A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_1A)
bc_DNA_RNA_2_1B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_1B)
bc_DNA_RNA_20A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R20A)
bc_DNA_RNA_20B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R20B)
bc_DNA_RNA_22A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R22A)
bc_DNA_RNA_22B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R22B)


#Count barcodes per variant per DNA and RNA sample, set minimum of 8 BC's per 
#variant in DNA sample, take median BC RNA/DNA per variant, then per variant 
#determine the median absolute deviation of all barcode ratios. 

ratio_bc_med_var <- function(df) {
  bc_count_DNA <- df %>%
    group_by(subpool, name, most_common) %>%
    summarize(barcodes_DNA = n()) %>%
    filter(barcodes_DNA > 7)
  bc_count_RNA <- df %>%
    group_by(subpool, name, most_common) %>%
    filter(num_reads_RNA != 0) %>%
    summarize(barcodes_RNA = n())
  bc_DNA_RNA <- left_join(bc_count_DNA, bc_count_RNA, 
                          by = c('subpool', 'name', 'most_common')) %>%
    ungroup()
  bc_min_8_df <- left_join(bc_DNA_RNA, df, 
                           by = c('subpool', 'name', 'most_common')) %>%
    ungroup()
  med_ratio <- bc_min_8_df %>%
    group_by(subpool, name, most_common) %>%
    summarize(med_ratio = median(ratio))
  mad_ratio <- bc_min_8_df %>%
    group_by(subpool, name, most_common) %>%
    summarize(mad = mad(ratio))
  med_mad <- inner_join(med_ratio, mad_ratio, 
                        by = c('subpool', 'name', 'most_common')) %>%
    mutate(mad_over_med = as.double(mad/med_ratio)) %>%
    mutate(mad_over_med = if_else(
      is.na(mad_over_med),
      as.double(0), 
      mad_over_med))
  bc_med <- inner_join(med_mad, bc_DNA_RNA, 
                       by = c('subpool', 'name', 'most_common')) %>%
    ungroup() %>%
    filter(med_ratio > 0)
  return(bc_med)
}

med_ratio_R0A <- ratio_bc_med_var(bc_DNA_RNA_0A)
med_ratio_R0B <- ratio_bc_med_var(bc_DNA_RNA_0B)
med_ratio_R2_5A <- ratio_bc_med_var(bc_DNA_RNA_2_5A)
med_ratio_R2_5B <- ratio_bc_med_var(bc_DNA_RNA_2_5B)
med_ratio_R2_4A <- ratio_bc_med_var(bc_DNA_RNA_2_4A)
med_ratio_R2_4B <- ratio_bc_med_var(bc_DNA_RNA_2_4B)
med_ratio_R2_3A <- ratio_bc_med_var(bc_DNA_RNA_2_3A)
med_ratio_R2_3B <- ratio_bc_med_var(bc_DNA_RNA_2_3B)
med_ratio_R2_2A <- ratio_bc_med_var(bc_DNA_RNA_2_2A)
med_ratio_R2_2B <- ratio_bc_med_var(bc_DNA_RNA_2_2B)
med_ratio_R2_1A <- ratio_bc_med_var(bc_DNA_RNA_2_1A)
med_ratio_R2_1B <- ratio_bc_med_var(bc_DNA_RNA_2_1B)
med_ratio_R20A <- ratio_bc_med_var(bc_DNA_RNA_20A)
med_ratio_R20B <- ratio_bc_med_var(bc_DNA_RNA_20B)
med_ratio_R22A <- ratio_bc_med_var(bc_DNA_RNA_22A)
med_ratio_R22B <- ratio_bc_med_var(bc_DNA_RNA_22B)

#Combine biological replicates

var_conc_rep_med <- function(df0A, df0B, df2_5A, df2_5B, df2_4A, df2_4B, df2_3A, 
                             df2_3B, df2_2A, df2_2B, df2_1A, df2_1B, df20A, df20B, 
                             df22A, df22B) {
  join_0 <- inner_join(df0A, df0B, 
                       by = c("name", "subpool", "most_common", "barcodes_DNA"), 
                       suffix = c("_0A", "_0B"))
  join_2_5 <- inner_join(df2_5A, df2_5B, 
                         by = c("name", "subpool", "most_common", 
                                "barcodes_DNA"), 
                         suffix = c("_2_5A", "_2_5B"))
  join_2_4 <- inner_join(df2_4A, df2_4B, 
                         by = c("name", "subpool", "most_common", 
                                "barcodes_DNA"), 
                         suffix = c("_2_4A", "_2_4B"))
  join_2_3 <- inner_join(df2_3A, df2_3B, 
                         by = c("name", "subpool", "most_common", 
                                "barcodes_DNA"), 
                         suffix = c("_2_3A", "_2_3B"))
  join_2_2 <- inner_join(df2_2A, df2_2B, 
                         by = c("name", "subpool", "most_common", 
                                "barcodes_DNA"), 
                         suffix = c("_2_2A", "_2_2B"))
  join_2_1 <- inner_join(df2_1A, df2_1B, 
                         by = c("name", "subpool", "most_common", 
                                "barcodes_DNA"), 
                         suffix = c("_2_1A", "_2_1B"))
  join_20 <- inner_join(df20A, df20B, 
                        by = c("name", "subpool", "most_common", 
                               "barcodes_DNA"), 
                        suffix = c("_20A", "_20B"))
  join_22 <- inner_join(df22A, df22B, 
                        by = c("name", "subpool", "most_common", 
                               "barcodes_DNA"), 
                        suffix = c("_22A", "_22B"))
  join_0_2_5 <- inner_join(join_0, join_2_5, 
                           by = c("name", "subpool", "most_common", 
                                  "barcodes_DNA"))
  join_0_2_4 <- inner_join(join_0_2_5, join_2_4, 
                           by = c("name", "subpool", "most_common", 
                                  "barcodes_DNA"))
  join_0_2_3 <- inner_join(join_0_2_4, join_2_3, 
                           by = c("name", "subpool", "most_common", 
                                  "barcodes_DNA"))
  join_0_2_2 <- inner_join(join_0_2_3, join_2_2, 
                           by = c("name", "subpool", "most_common", 
                                  "barcodes_DNA"))
  join_0_2_1 <- inner_join(join_0_2_2, join_2_1, 
                           by = c("name", "subpool", "most_common", 
                                  "barcodes_DNA"))
  join_0_20 <- inner_join(join_0_2_1, join_20, 
                          by = c("name", "subpool", "most_common", 
                                 "barcodes_DNA"))
  join_0_22 <- inner_join(join_0_20, join_22, 
                          by = c("name", "subpool", "most_common", 
                                 "barcodes_DNA")) %>%
    ungroup()
  print('processed dfs in order of samples: 0A, 0B, 2_5A, 2_5B, 2_4A, 2_4B, 
        2_3A, 2_3B, 2_2A, 2_2B, 2_1A, 2_1B, 20A, 20B, 22A, 22B')
  return(join_0_22)
}

med_rep_0_22_A_B <- var_conc_rep_med(med_ratio_R0A, med_ratio_R0B, 
                                     med_ratio_R2_5A, med_ratio_R2_5B,
                                     med_ratio_R2_4A, med_ratio_R2_4B,
                                     med_ratio_R2_3A, med_ratio_R2_3B,
                                     med_ratio_R2_2A, med_ratio_R2_2B,
                                     med_ratio_R2_1A, med_ratio_R2_1B,
                                     med_ratio_R20A, med_ratio_R20B,
                                     med_ratio_R22A, med_ratio_R22B)

#Figure 1D, episomal------------------------------------------------------------

#Plot replicability with backgrounds in orange and positive control in red

p_fig1_epi_med_rep <- med_rep_0_22_A_B %>%
  ggplot(aes(med_ratio_22A, med_ratio_22B)) +
  geom_point(alpha = 0.1, size = 0.75) +
  geom_point(data = filter(med_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             fill = 'orange', shape = 21, size = 1.75) + 
  geom_point(data = filter(med_rep_0_22_A_B, 
                           name == 'pGL4.29 Promega 1-63 + 1-87'), 
             fill = 'red', shape = 21, size = 1.75) +
  annotation_logticks(scaled = TRUE) +
  xlab("Expression (a.u.) replicate 1") +
  ylab("Expression (a.u.) replicate 2") +
  scale_x_log10(limits = c(0.03, 20), breaks = c(0.1, 1, 10)) + 
  scale_y_log10(limits = c(0.03, 20), breaks = c(0.1, 1, 10)) +
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.line.y = element_line(), panel.spacing.x=unit(1, "lines")) +
  figurefont_theme

ggsave('final_figs/p_fig1_epi_med_rep.png', p_fig1_epi_med_rep,
       width = 2, height = 2, units = 'in')

log10_med_rep_0_22_A_B <- var_log10(med_rep_0_22_A_B)

pearsons <- tibble(
  sample = c('all', 'subpool3', 'subpool5'),
  pearsons = c(cor(log10_med_rep_0_22_A_B$med_ratio_22A, 
                   log10_med_rep_0_22_A_B$med_ratio_22B, 
                   use = "pairwise.complete.obs", method = "pearson"),
               cor(filter(log10_med_rep_0_22_A_B, 
                          subpool == 'subpool3')$med_ratio_22A,
                   filter(log10_med_rep_0_22_A_B, 
                          subpool == 'subpool3')$med_ratio_22B,
                   use = "pairwise.complete.obs", method = "pearson"),
               cor(filter(log10_med_rep_0_22_A_B, 
                          subpool == 'subpool5')$med_ratio_22A,
                   filter(log10_med_rep_0_22_A_B, 
                          subpool == 'subpool5')$med_ratio_22B,
                   use = "pairwise.complete.obs", method = "pearson")))

#


