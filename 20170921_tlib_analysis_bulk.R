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

cbPalette7 <- c('#440154FF', '#39568CFF', '#287D8EFF', '#20A387FF', '#73D055FF',
                '#B8DE29FF', '#FDE725FF')
cbPalette7_grad_light <- c('white', '#FDE725FF', '#B8DE29FF', '#55C667FF', 
                           '#1F968BFF', '#39568CFF', '#482677FF')
spacing_5_20_palette <- c('gray20', 'dodgerblue3', 'indianred2', '#55C667FF')

#Written for the analysis of a range of inductions in the transient library
#DNA_BC: DNA
#R0A_BC and R0B_BC: RNA at 0 µM Forsk replicate A and B
#R2#A_BC and R2#B_BC: RNA at 2^# µM Forsk replicate A and B, includes negative 
#and positive #'s

#tested concentrations: 0, 2^-5, 2^-4, 2^-3, 2^-2, 2^-1, 2^0, 2^2 µM Forsk


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


#Load barcode mapping table, remember sequences are rcomp due to sequencing 
#format

barcode_map <- read_tsv('../../BCMap/uniqueSP2345.txt', 
                        col_names = c(
                          'fluff', 'barcode', 'name', 'most_common'), 
                        skip = 1) %>%
  select(-fluff) %>%
  mutate(subpool = ifelse(startsWith(name, 'subpool'), 
                          substr(name, 1, 8), 
                          'control')) 


#pick out SP3 and SP5 in the bcmap that were used in this assay

SP3_SP5_map <- barcode_map %>%
  mutate(subpool = ifelse(
    startsWith(name, 'subpool'), 
    substr(name, 1, 8), 
    'control')) %>%
  filter(subpool != 'subpool2') %>%
  filter(subpool != 'subpool4')


#Join reads to bcmap------------------------------------------------------------

#Join BC reads to BC mapping, keeping the reads only appearing in barcode 
#mapping and replacing na with 0 reads.

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


#Determine variant counts by summing--------------------------------------------

#sum unique barcodes and normalized bc reads per variant. Set-up minimum BC's 
#per variant in DNA at 8 for reliability

var_sum_bc_num <- function(df1) {
  bc_count <- df1 %>%
    filter(df1$num_reads > 0) %>%
    group_by(subpool, name, most_common) %>%
    summarise(barcodes = n())
  variant_sum <- df1 %>%
    group_by(subpool, name, most_common) %>%
    count(name, wt = normalized) %>%
    rename(sum = n)
  bc_sum <- right_join(variant_sum, bc_count, 
                       by = c("name", "subpool", "most_common")) %>%
    ungroup()
  return(bc_sum)
}

variant_counts_DNA <- var_sum_bc_num(bc_join_DNA) %>%
  filter(barcodes > 7)
variant_counts_R0A <- var_sum_bc_num(bc_join_R0A)
variant_counts_R0B <- var_sum_bc_num(bc_join_R0B)
variant_counts_R2_5A <- var_sum_bc_num(bc_join_R2_5A)
variant_counts_R2_5B <- var_sum_bc_num(bc_join_R2_5B)
variant_counts_R2_4A <- var_sum_bc_num(bc_join_R2_4A)
variant_counts_R2_4B <- var_sum_bc_num(bc_join_R2_4B)
variant_counts_R2_3A <- var_sum_bc_num(bc_join_R2_3A)
variant_counts_R2_3B <- var_sum_bc_num(bc_join_R2_3B)
variant_counts_R2_2A <- var_sum_bc_num(bc_join_R2_2A)
variant_counts_R2_2B <- var_sum_bc_num(bc_join_R2_2B)
variant_counts_R2_1A <- var_sum_bc_num(bc_join_R2_1A)
variant_counts_R2_1B <- var_sum_bc_num(bc_join_R2_1B)
variant_counts_R20A <- var_sum_bc_num(bc_join_R20A)
variant_counts_R20B <- var_sum_bc_num(bc_join_R20B)
variant_counts_R22A <- var_sum_bc_num(bc_join_R22A)
variant_counts_R22B <- var_sum_bc_num(bc_join_R22B)


#Join RNA to DNA and determine expression from summing--------------------------

#Left_join DNA with RNA and determine expression ratio as sum nrpm RNA/sum nrpm
#DNA. Filter out 0 sum expression.

var_expression <- function(df1, df2) {
  RNA_DNA <- left_join(df1, df2, 
                        by = c("name", "subpool", "most_common"), 
                        suffix = c("_DNA", "_RNA")) %>%
    filter(sum_RNA > 0) %>%
    mutate(ratio = sum_RNA / sum_DNA)
  print('x defined as DNA, y defined as RNA in var_expression(x,y)')
  return(RNA_DNA)
}

RNA_DNA_0A <- var_expression(variant_counts_DNA, variant_counts_R0A)
RNA_DNA_0B <- var_expression(variant_counts_DNA, variant_counts_R0B)
RNA_DNA_2_5A <- var_expression(variant_counts_DNA, variant_counts_R2_5A)
RNA_DNA_2_5B <- var_expression(variant_counts_DNA, variant_counts_R2_5B)
RNA_DNA_2_4A <- var_expression(variant_counts_DNA, variant_counts_R2_4A)
RNA_DNA_2_4B <- var_expression(variant_counts_DNA, variant_counts_R2_4B)
RNA_DNA_2_3A <- var_expression(variant_counts_DNA, variant_counts_R2_3A)
RNA_DNA_2_3B <- var_expression(variant_counts_DNA, variant_counts_R2_3B)
RNA_DNA_2_2A <- var_expression(variant_counts_DNA, variant_counts_R2_2A)
RNA_DNA_2_2B <- var_expression(variant_counts_DNA, variant_counts_R2_2B)
RNA_DNA_2_1A <- var_expression(variant_counts_DNA, variant_counts_R2_1A)
RNA_DNA_2_1B <- var_expression(variant_counts_DNA, variant_counts_R2_1B)
RNA_DNA_20A <- var_expression(variant_counts_DNA, variant_counts_R20A)
RNA_DNA_20B <- var_expression(variant_counts_DNA, variant_counts_R20B)
RNA_DNA_22A <- var_expression(variant_counts_DNA, variant_counts_R22A)
RNA_DNA_22B <- var_expression(variant_counts_DNA, variant_counts_R22B)


#combine biological replicates--------------------------------------------------

var_conc_rep <- function(df0A, df0B, df2_5A, df2_5B, df2_4A, df2_4B, df2_3A, 
                         df2_3B, df2_2A, df2_2B, df2_1A, df2_1B, df20A, df20B, 
                         df22A, df22B) {
  join_0 <- inner_join(df0A, df0B, 
                       by = c("name", "subpool", "most_common", "sum_DNA", 
                         "barcodes_DNA"), 
                       suffix = c("_0A", "_0B"))
  join_2_5 <- inner_join(df2_5A, df2_5B, 
                         by = c("name", "subpool", "most_common", "sum_DNA", 
                                "barcodes_DNA"), 
                         suffix = c("_2_5A", "_2_5B"))
  join_2_4 <- inner_join(df2_4A, df2_4B, 
                         by = c("name", "subpool", "most_common", "sum_DNA", 
                                "barcodes_DNA"), 
                         suffix = c("_2_4A", "_2_4B"))
  join_2_3 <- inner_join(df2_3A, df2_3B, 
                         by = c("name", "subpool", "most_common", "sum_DNA", 
                                "barcodes_DNA"), 
                         suffix = c("_2_3A", "_2_3B"))
  join_2_2 <- inner_join(df2_2A, df2_2B, 
                         by = c("name", "subpool", "most_common", "sum_DNA", 
                                "barcodes_DNA"), 
                         suffix = c("_2_2A", "_2_2B"))
  join_2_1 <- inner_join(df2_1A, df2_1B, 
                         by = c("name", "subpool", "most_common", "sum_DNA", 
                                "barcodes_DNA"), 
                         suffix = c("_2_1A", "_2_1B"))
  join_20 <- inner_join(df20A, df20B, 
                        by = c("name", "subpool", "most_common", "sum_DNA", 
                               "barcodes_DNA"), 
                        suffix = c("_20A", "_20B"))
  join_22 <- inner_join(df22A, df22B, 
                        by = c("name", "subpool", "most_common", "sum_DNA", 
                               "barcodes_DNA"), 
                        suffix = c("_22A", "_22B"))
  join_0_2_5 <- inner_join(join_0, join_2_5, 
                           by = c("name", "subpool", "most_common", "sum_DNA", 
                                  "barcodes_DNA"))
  join_0_2_4 <- inner_join(join_0_2_5, join_2_4, 
                           by = c("name", "subpool", "most_common", "sum_DNA", 
                                  "barcodes_DNA"))
  join_0_2_3 <- inner_join(join_0_2_4, join_2_3, 
                           by = c("name", "subpool", "most_common", "sum_DNA", 
                                  "barcodes_DNA"))
  join_0_2_2 <- inner_join(join_0_2_3, join_2_2, 
                           by = c("name", "subpool", "most_common", "sum_DNA", 
                                  "barcodes_DNA"))
  join_0_2_1 <- inner_join(join_0_2_2, join_2_1, 
                           by = c("name", "subpool", "most_common", "sum_DNA", 
                                  "barcodes_DNA"))
  join_0_20 <- inner_join(join_0_2_1, join_20, 
                          by = c("name", "subpool", "most_common", "sum_DNA", 
                                 "barcodes_DNA"))
  join_0_22 <- inner_join(join_0_20, join_22, 
                          by = c("name", "subpool", "most_common", "sum_DNA", 
                                 "barcodes_DNA")) %>%
    ungroup()
  print('processed dfs in order of samples: 0A, 0B, 2_5A, 2_5B, 2_4A, 2_4B, 
        2_3A, 2_3B, 2_2A, 2_2B, 2_1A, 2_1B, 20A, 20B, 22A, 22B')
  return(join_0_22)
}

rep_0_22_A_B <- var_conc_rep(RNA_DNA_0A, RNA_DNA_0B, RNA_DNA_2_5A, RNA_DNA_2_5B,
                             RNA_DNA_2_4A, RNA_DNA_2_4B, RNA_DNA_2_3A, 
                             RNA_DNA_2_3B, RNA_DNA_2_2A, RNA_DNA_2_2B, 
                             RNA_DNA_2_1A, RNA_DNA_2_1B, RNA_DNA_20A, 
                             RNA_DNA_20B, RNA_DNA_22A, RNA_DNA_22B)

#determine the log(RNA/DNA) for each sample (this takes the log of sum_RNA and 
#sum_DNA as well). Log2 is useful for replicate plots for expression and log10 
#is useful for barcode read analysis. This is useful for replicate plots, but 
#further manipulations should use trans_back_norm_rep_0_22_A_B

var_log2 <- function(df) {
  log_ratio_df <- df %>% 
    mutate_if(is.double, funs(log2(.)))
  return(log_ratio_df)
}

var_log10 <- function(df) {
  log_ratio_df <- df %>% 
    mutate_if(is.double, funs(log10(.)))
  return(log_ratio_df)
}

log10_rep_0_22_A_B <- var_log10(rep_0_22_A_B)


#Median analysis----------------------------------------------------------------

#Join DNA BC reads > 6 and join to RNA. Take ratio of RNA/DNA norm reads

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


#Count barcodes per variant per DNA and RNA, set minimum of 8 BC's per variant 
#in DNA, take median RNA/DNA per variant, then per variant determine 
#the median absolute deviation of all barcode ratios. Then filter out variants 
#with 0 median expression

ratio_bc_med_var <- function(df) {
  bc_count_DNA <- df %>%
    group_by(subpool, name, most_common) %>%
    summarize(barcodes_DNA = n()) %>%
    filter(barcodes_DNA > 7)
  bc_count_RNA <- df %>%
    group_by(subpool, name, most_common) %>%
    filter(num_reads_RNA != 0) %>%
    summarize(barcodes_RNA = n())
  bc_DNA_RNA <- inner_join(bc_count_DNA, bc_count_RNA, 
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

log10_med_rep_0_22_A_B <- var_log10(med_rep_0_22_A_B)

ggplot(med_rep_0_22_A_B, aes(med_ratio_22A, med_ratio_22B)) +
  geom_point(alpha = 0.1) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks(sides = 'bl')


#Plot median vs. sum for uninduced and induced

med_sum_join <- function(sum, med) {
  sum_select <- sum %>%
    select(name, ratio_22A, ratio_22B)
  med_select <- med %>%
    select(name, med_ratio_22A, med_ratio_22B)
  sum_med <- inner_join(sum_select, med_select, by = 'name') %>%
    select(-name)
  return(sum_med)
}

med_vs_sum <- med_sum_join(rep_0_22_A_B, med_rep_0_22_A_B) %>%
  var_log10()

p_ave_med_vs_sum <- med_vs_sum %>%
  mutate(ave_sum = (ratio_22A + ratio_22B)/2) %>%
  mutate(ave_med = (med_ratio_22A + med_ratio_22B)/2) %>%
  ggplot(aes(ave_sum, ave_med)) +
  geom_point(alpha = 0.1)

my_points <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.1, size = 0.75) +
    scale_x_continuous(limits = c(-1.5, 1.5), breaks = c(-1:1)) + 
    scale_y_continuous(limits = c(-1.5, 1.5), breaks = c(-1:1)) +
    annotation_logticks(sides = 'bl')
}

my_density <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_density(kernel = 'gaussian') +
    scale_x_continuous(limits = c(-1.5, 1.5), breaks = c(-3:1)) +
    scale_y_continuous(limits = c(-0.1, 1.5)) +
    annotation_logticks(sides = 'b')
}

p_med_vs_sum <- ggpairs(med_vs_sum, 
                        columnLabels = c('Sum Exp.\nRep. 1', 'Sum Exp.\nRep. 2',
                                         'Med Exp.\nRep. 1', 'Med Exp.\nRep. 2'),
                        lower = list(continuous = my_points),
                        diag = list(continuous = my_density)) +
  panel_border() + 
  theme(panel.grid.major = element_blank())

save_plot('plots/p_med_vs_sum_22.pdf', p_med_vs_sum, scale = 1.5)
  


#Normalize to background--------------------------------------------------------

back_norm <- function(df1) {
  gsub_0_22 <- df1 %>%
    ungroup () %>%
    filter(subpool != 'control') %>%
    mutate(
      name = gsub('Smith R. Vista chr9:83712599-83712766', 'v chr9', name),
      name = gsub('Vista Chr5:88673410-88674494', 'v chr5', name),
      name = gsub('scramble pGL4.29 Promega 1-63 \\+ 1-87', 's pGl4', name)
    ) %>%
    mutate(background = name) %>%
    mutate(background = str_sub(background, 
                                nchar(background)-5, 
                                nchar(background)))
  backgrounds <- gsub_0_22 %>%
    filter(startsWith(name, 
                      'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')) %>%
    select(background, ratio_0A, ratio_0B, ratio_2_5A, ratio_2_5B, ratio_2_4A, 
           ratio_2_4B, ratio_2_3A, ratio_2_3B, ratio_2_2A, ratio_2_2B, 
           ratio_2_1A, ratio_2_1B, ratio_20A, ratio_20B, ratio_22A, 
           ratio_22B) %>%
    rename(ratio_0A_back = ratio_0A) %>%
    rename(ratio_0B_back = ratio_0B) %>%
    rename(ratio_2_5A_back = ratio_2_5A) %>%
    rename(ratio_2_5B_back = ratio_2_5B) %>%
    rename(ratio_2_4A_back = ratio_2_4A) %>%
    rename(ratio_2_4B_back = ratio_2_4B) %>%
    rename(ratio_2_3A_back = ratio_2_3A) %>%
    rename(ratio_2_3B_back = ratio_2_3B) %>%
    rename(ratio_2_2A_back = ratio_2_2A) %>%
    rename(ratio_2_2B_back = ratio_2_2B) %>%
    rename(ratio_2_1A_back = ratio_2_1A) %>%
    rename(ratio_2_1B_back = ratio_2_1B) %>%
    rename(ratio_20A_back = ratio_20A) %>%
    rename(ratio_20B_back = ratio_20B) %>%
    rename(ratio_22A_back = ratio_22A) %>%
    rename(ratio_22B_back = ratio_22B)
  back_join_norm <- left_join(gsub_0_22, backgrounds, by = 'background') %>%
    mutate(ratio_0A_norm = ratio_0A/ratio_0A_back) %>%
    mutate(ratio_0B_norm = ratio_0B/ratio_0B_back) %>%
    mutate(ratio_2_5A_norm = ratio_2_5A/ratio_2_5A_back) %>%
    mutate(ratio_2_5B_norm = ratio_2_5B/ratio_2_5B_back) %>%
    mutate(ratio_2_4A_norm = ratio_2_4A/ratio_2_4A_back) %>%
    mutate(ratio_2_4B_norm = ratio_2_4B/ratio_2_4B_back) %>%
    mutate(ratio_2_3A_norm = ratio_2_3A/ratio_2_3A_back) %>%
    mutate(ratio_2_3B_norm = ratio_2_3B/ratio_2_3B_back) %>%
    mutate(ratio_2_2A_norm = ratio_2_2A/ratio_2_2A_back) %>%
    mutate(ratio_2_2B_norm = ratio_2_2B/ratio_2_2B_back) %>%
    mutate(ratio_2_1A_norm = ratio_2_1A/ratio_2_1A_back) %>%
    mutate(ratio_2_1B_norm = ratio_2_1B/ratio_2_1B_back) %>%
    mutate(ratio_20A_norm = ratio_20A/ratio_20A_back) %>%
    mutate(ratio_20B_norm = ratio_20B/ratio_20B_back) %>%
    mutate(ratio_22A_norm = ratio_22A/ratio_22A_back) %>%
    mutate(ratio_22B_norm = ratio_22B/ratio_22B_back) %>%
    mutate(ave_ratio_0_norm = (ratio_0A_norm + ratio_0B_norm)/2) %>%
    mutate(ave_ratio_2_5_norm = (ratio_2_5A_norm + ratio_2_5B_norm)/2) %>%
    mutate(ave_ratio_2_4_norm = (ratio_2_4A_norm + ratio_2_4B_norm)/2) %>%
    mutate(ave_ratio_2_3_norm = (ratio_2_3A_norm + ratio_2_3B_norm)/2) %>%
    mutate(ave_ratio_2_2_norm = (ratio_2_2A_norm + ratio_2_2B_norm)/2) %>%
    mutate(ave_ratio_2_1_norm = (ratio_2_1A_norm + ratio_2_1B_norm)/2) %>%
    mutate(ave_ratio_20_norm = (ratio_20A_norm + ratio_20B_norm)/2) %>%
    mutate(ave_ratio_22_norm = (ratio_22A_norm + ratio_22B_norm)/2) %>%
    mutate(induction = ave_ratio_22_norm/ave_ratio_0_norm)
  return(back_join_norm)
}

trans_back_norm_rep_0_22 <- back_norm(rep_0_22_A_B)

#Add pGl4.29 pc's into background-normalization

trans_back_norm_pc_spGl4 <- rep_0_22_A_B %>%
  filter(name == 'pGL4.29 Promega 1-63 + 1-87' | name == 'pGL4.29 Promega 1-87') %>%
  mutate(name = str_c(name, '_scramble pGL4.29 Promega 1-63 + 1-87')) %>%
  mutate(subpool = 'subpool3') %>%
  rbind(rep_0_22_A_B) %>%
  back_norm()

log10_trans_back_norm_rep_0_22 <- var_log10(trans_back_norm_rep_0_22)


#Make df with concentration and expression as a variable------------------------

#Make untidy data with expression and concentration as variables

var_conc_exp_nonnorm <- function(df) {
  df <- df %>%
    ungroup() %>%
    filter(subpool != 'control') %>%
    mutate(
      name = gsub('Smith R. Vista chr9:83712599-83712766', 'v chr9', name),
      name = gsub('Vista Chr5:88673410-88674494', 'v chr5', name), 
      name = gsub('scramble pGL4.29 Promega 1-63 \\+ 1-87', 's pGl4', name)
      ) %>%
      mutate(background = name) %>%
      mutate(background = str_sub(background, 
                                  nchar(background)-5, 
                                  nchar(background)))
  df_0 <- df %>%
    mutate(ave_barcode = (barcodes_RNA_0A + barcodes_RNA_0B)/2) %>%
    mutate(ave_ratio = (ratio_0A + ratio_0B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode, ave_ratio) %>%
    mutate(conc = 0)
  df_2_5 <- df %>%
    mutate(ave_barcode = (barcodes_RNA_2_5A + barcodes_RNA_2_5B)/2) %>%
    mutate(ave_ratio = (ratio_2_5A + ratio_2_5B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode, ave_ratio) %>%
    mutate(conc = 2^-5)
  df_2_4 <- df %>%
    mutate(ave_barcode = (barcodes_RNA_2_4A + barcodes_RNA_2_4B)/2) %>%
    mutate(ave_ratio = (ratio_2_4A + ratio_2_4B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode, ave_ratio) %>%
    mutate(conc = 2^-4)
  df_2_3 <- df %>%
    mutate(ave_barcode = (barcodes_RNA_2_3A + barcodes_RNA_2_3B)/2) %>%
    mutate(ave_ratio = (ratio_2_3A + ratio_2_3B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode, ave_ratio) %>%
    mutate(conc = 2^-3)
  df_2_2 <- df %>%
    mutate(ave_barcode = (barcodes_RNA_2_2A + barcodes_RNA_2_2B)/2) %>%
    mutate(ave_ratio = (ratio_2_2A + ratio_2_2B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode, ave_ratio) %>%
    mutate(conc = 2^-2)
  df_2_1 <- df %>%
    mutate(ave_barcode = (barcodes_RNA_2_1A + barcodes_RNA_2_1B)/2) %>%
    mutate(ave_ratio = (ratio_2_1A + ratio_2_1B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode, ave_ratio) %>%
    mutate(conc = 2^-1)
  df_20 <- df %>%
    mutate(ave_barcode = (barcodes_RNA_20A + barcodes_RNA_20B)/2) %>%
    mutate(ave_ratio = (ratio_20A + ratio_20B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode, ave_ratio) %>%
    mutate(conc = 2^0)
  df_22 <- df %>%
    mutate(ave_barcode = (barcodes_RNA_22A + barcodes_RNA_22B)/2) %>%
    mutate(ave_ratio = (ratio_22A + ratio_22B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode, ave_ratio) %>%
    mutate(conc = 2^2)
  df_0_22 <- rbind(df_0, df_2_5, df_2_4, df_2_3, df_2_2, df_2_1, df_20, df_22)
  return(df_0_22)
}

trans_conc <- var_conc_exp_nonnorm(rep_0_22_A_B)


#Make untidy data with background-normalized data

var_conc_exp <- function(df) {
  df_0 <- df %>%
    mutate(ave_barcode_0 = (barcodes_RNA_0A + barcodes_RNA_0B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_0, 
           ave_ratio_0_norm) %>%
    mutate(conc = 0) %>%
    rename(ave_ratio_norm = ave_ratio_0_norm) %>%
    rename(ave_barcode = ave_barcode_0)
  df_2_5 <- df %>%
    mutate(ave_barcode_2_5 = (barcodes_RNA_2_5A + barcodes_RNA_2_5B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_2_5, 
           ave_ratio_2_5_norm) %>%
    mutate(conc = 2^-5) %>%
    rename(ave_ratio_norm = ave_ratio_2_5_norm) %>%
    rename(ave_barcode = ave_barcode_2_5)
  df_2_4 <- df %>%
    mutate(ave_barcode_2_4 = (barcodes_RNA_2_4A + barcodes_RNA_2_4B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_2_4, 
           ave_ratio_2_4_norm) %>%
    mutate(conc = 2^-4) %>%
    rename(ave_ratio_norm = ave_ratio_2_4_norm) %>%
    rename(ave_barcode = ave_barcode_2_4)
  df_2_3 <- df %>%
    mutate(ave_barcode_2_3 = (barcodes_RNA_2_3A + barcodes_RNA_2_3B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_2_3, 
           ave_ratio_2_3_norm) %>%
    mutate(conc = 2^-3) %>%
    rename(ave_ratio_norm = ave_ratio_2_3_norm) %>%
    rename(ave_barcode = ave_barcode_2_3)
  df_2_2 <- df %>%
    mutate(ave_barcode_2_2 = (barcodes_RNA_2_2A + barcodes_RNA_2_2B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_2_2, 
           ave_ratio_2_2_norm) %>%
    mutate(conc = 2^-2) %>%
    rename(ave_ratio_norm = ave_ratio_2_2_norm) %>%
    rename(ave_barcode = ave_barcode_2_2)
  df_2_1 <- df %>%
    mutate(ave_barcode_2_1 = (barcodes_RNA_2_1A + barcodes_RNA_2_1B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_2_1, 
           ave_ratio_2_1_norm) %>%
    mutate(conc = 2^-1) %>%
    rename(ave_ratio_norm = ave_ratio_2_1_norm) %>%
    rename(ave_barcode = ave_barcode_2_1)
  df_20 <- df %>%
    mutate(ave_barcode_20 = (barcodes_RNA_20A + barcodes_RNA_20B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_20, 
           ave_ratio_20_norm) %>%
    mutate(conc = 2^0) %>%
    rename(ave_ratio_norm = ave_ratio_20_norm) %>%
    rename(ave_barcode = ave_barcode_20)
  df_22 <- df %>%
    mutate(ave_barcode_22 = (barcodes_RNA_22A + barcodes_RNA_22B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_22, 
           ave_ratio_22_norm) %>%
    mutate(conc = 2^2) %>%
    rename(ave_ratio_norm = ave_ratio_22_norm) %>%
    rename(ave_barcode = ave_barcode_22)
  df_0_22 <- rbind(df_0, df_2_5, df_2_4, df_2_3, df_2_2, df_2_1, df_20, df_22)
  return(df_0_22)
}

trans_back_norm_conc <- var_conc_exp(trans_back_norm_rep_0_22)

trans_back_norm_pc_spGl4_conc <- var_conc_exp(trans_back_norm_pc_spGl4)


#Separate into subpools---------------------------------------------------------

#Subpool 3 contains 2 consensus binding sites with flanks (ATTGACGTCAGC) that 
#vary in distance from one another by 0 (no inner flanks), 5, 10, 15, 20 and 
#70 bp (all but 0 appear as -4 bp spacing). Each site distance combination is 
#then moved along the backgrounds at 1 bp increments starting from closest to 
#the minP. Separation lists the spacing between sites and distance (start of 
#consensus and flanks). Added 2 to all distances to measure end of background to
#start of BS and then added 64 to measure to minimal promoter. Added 4 to all 
#spacing but 0 to measure difference between start of sites

subpool3 <- function(df) {
  df <- df %>%
    filter(subpool == "subpool3") %>%
    ungroup() %>%
    select(-subpool) %>%
    mutate(name = gsub('2BS ', '', name), 
           name = gsub(' bp spacing ', '_', name)) %>%
    separate(name, 
             into = c("subpool", "spacing", "fluff2", "fluff3", "dist", "fluff4"),
             sep = "_", convert = TRUE) %>%
    select(-subpool, -fluff2, -fluff3, -fluff4) %>%
    mutate(dist = as.integer(dist + 2 + 64)) %>%
    mutate(spacing = 
             ifelse(spacing != as.integer(0), 
                    as.integer(spacing + 4), as.integer(spacing)))
}

s3_tidy <- rep_0_22_A_B %>%
  ungroup() %>%
  mutate(
    name = gsub('Smith R. Vista chr9:83712599-83712766', 'v chr9', name),
    name = gsub('Vista Chr5:88673410-88674494', 'v chr5', name),
    name = gsub('scramble pGL4.29 Promega 1-63 \\+ 1-87', 's pGl4', name)
  ) %>%
  mutate(background = name) %>%
  mutate(background = str_sub(background, 
                              nchar(background)-5,
                              nchar(background))) %>%
  subpool3()

s3_untidy <- subpool3(trans_back_norm_conc)
  

#Subpool 5 contains 6 equally spaced sites spaced 13 bp apart and starting from 
#furthest to the minP. These sites are filled with sites of either the consensus
#site, a weak site or no site. Both the weak and consensus sites are flanked by 
#the same flanking sequence. It's too complicated now to rename these sites 
#instead of the 1 -> 6 system and instead with actual distance to the minimal 
#promoter (site 1, 2, 3, 4, 5, 6 equate to -191, -166, -141, -116, -91 and -66 
#respectively)

subpool5 <- function(df) {
  df <- df %>%
    filter(subpool == "subpool5") %>%
    ungroup() %>%
    select(-subpool) %>%
    mutate(name = gsub('no_site', 'nosite', name)) %>%
    separate(name, into = c("subpool", "site1", "site2", "site3", "site4", 
                            "site5", "site6", "fluff"), sep = "_") %>%
    select(-subpool, -fluff) %>%
    mutate(consensus = str_detect(site1, "consensus") + 
             str_detect(site2, "consensus") + 
             str_detect(site3, "consensus") + 
             str_detect(site4, "consensus") + 
             str_detect(site5, "consensus") + 
             str_detect(site6, "consensus")) %>%
    mutate(weak = str_detect(site1, "weak") +
             str_detect(site2, "weak") +
             str_detect(site3, "weak") +
             str_detect(site4, "weak") +
             str_detect(site5, "weak") +
             str_detect(site6, "weak")) %>%
    mutate(nosite = str_detect(site1, "nosite") +
             str_detect(site2, "nosite") +
             str_detect(site3, "nosite") +
             str_detect(site4, "nosite") +
             str_detect(site5, "nosite") +
             str_detect(site6, "nosite")) %>%
    mutate(total_sites = consensus + weak) %>%
    mutate(site_combo = 
             ifelse(weak == 0 & consensus > 0, 
                    'consensus', 'mixed')) %>%
    mutate(site_combo = 
             ifelse(consensus == 0 & weak > 0, 
                    'weak', site_combo)) %>%
    mutate(site_combo = 
             ifelse(consensus == 0 & weak == 0, 
                    'none', site_combo)) 
}

s5_tidy <- rep_0_22_A_B %>%
  ungroup() %>%
  mutate(
    name = gsub('Smith R. Vista chr9:83712599-83712766', 'v chr9', name),
    name = gsub('Vista Chr5:88673410-88674494', 'v chr5', name),
    name = gsub('scramble pGL4.29 Promega 1-63 \\+ 1-87', 's pGl4', name)
  ) %>%
  mutate(background = name) %>%
  mutate(background = str_sub(background, 
                              nchar(background)-5,
                              nchar(background))) %>%
  subpool5()

s5_untidy <- subpool5(trans_back_norm_conc)
  
#Do I still look at controls?
controls <- 
  filter(rep_0_22_A_B, subpool == "control") %>%
  ungroup() %>%
  mutate(ave_ratio_0 = (ratio_0A + ratio_0B)/2) %>%
  mutate(ave_ratio_2_5 = (ratio_2_5A + ratio_2_5B)/2) %>%
  mutate(ave_ratio_2_4 = (ratio_2_4A + ratio_2_4B)/2) %>%
  mutate(ave_ratio_2_3 = (ratio_2_3A + ratio_2_3B)/2) %>%
  mutate(ave_ratio_2_2 = (ratio_2_2A + ratio_2_2B)/2) %>%
  mutate(ave_ratio_2_1 = (ratio_2_1A + ratio_2_1B)/2) %>%
  mutate(ave_ratio_20 = (ratio_20A + ratio_20B)/2) %>%
  mutate(ave_ratio_22 = (ratio_22A + ratio_22B)/2)
  

#Plot subpool3 expression features----------------------------------------------

#Rotation about a dna with 5 and 10 bp spacing

df <- tibble(bp = seq(from = 0, to = 10, by = 1), 
             main = seq(from = 0, to = 343, by = 34.3),
             five = seq(from = 85.9, to = 428.9, by = 34.3),
             ten = seq(from = 257.4, to = 600.4, by = 34.3)) %>%
  mutate(five = if_else(five >= 360, five-360, five)) %>%
  mutate(ten = if_else(ten >= 360, ten-360, ten)) %>%
  mutate(main = if_else((main > 90 & main < 180), 90-(main-90), main)) %>%
  mutate(main = if_else((main > 180 & main < 270), -(main-180), main)) %>%
  mutate(main = if_else((main > 270 & main < 360), -(90-(main-270)), main)) %>%
  mutate(five = if_else((five > 90 & five < 180), 90-(five-90), five)) %>%
  mutate(five = if_else((five > 180 & five < 270), -(five-180), five)) %>%
  mutate(five = if_else((five > 270 & five < 360), -(90-(five-270)), five)) %>%
  mutate(ten = if_else((ten > 90 & ten < 180), 90-(ten-90), ten)) %>%
  mutate(ten = if_else((ten > 180 & ten < 270), -(ten-180), ten)) %>%
  mutate(ten = if_else((ten > 270 & ten < 360), -(90-(ten-270)), ten))

ggplot(df, aes(bp)) +
  geom_point(aes(y = main)) +
  geom_line(aes(y = main)) +
  geom_point(aes(y = five), color = 'red') +
  geom_line(aes(y = five), color = 'red') +
  geom_point(aes(y = ten), color = 'blue') +
  geom_line(aes(y = ten), color = 'blue')

bp <- seq(0, 10, 1)
main <- tibble(main = sin((1/10.5*(2*pi))*bp)) 
five <- tibble(five = sin((1/10.5*(2*pi))*bp + 1.5))
ten <- tibble(ten = sin((1/10.5*(2*pi))*bp + 4.49))
all <- cbind(bp, main, five, ten) %>%
  gather(main, five, ten, key = 'CREB', value = 'y')

ggplot(all, aes(bp, y, color = CREB)) +
  geom_line()

main_five <- tibble(main_five = 2*(cos(1.5/2)*sin((1/10.5*(2*pi))*bp + 1.5/2)))
main_ten <- tibble(main_ten = 2*(cos(4.49/2)*sin((1/10.5*(2*pi))*bp + 4.49/2)))
main_all <- cbind(bp, main_five, main_ten)  %>%
  gather(main_five, main_ten, key = 'spacing', value = 'y')

ggplot(main_all, aes(bp, y, color = spacing)) +
  geom_line() + 
  geom_point()


#v chr9 spacing overlay plots using 3 bp moving average and plot as line over 
#original data in overlays

library(caTools)

moveavg_dist3 <- function(df) {
  df <- df %>%
    mutate(ave_3 = runmean(ave_ratio_22, 3, alg = 'R', endrule = 'NA'))
}

s3_tidy_moveavg3 <- s3_tidy %>%
  mutate(ave_ratio_22 = (ratio_22A + ratio_22B)/2) %>%
  select(background, spacing, dist, ave_ratio_22) %>%
  group_by(background, spacing) %>%
  arrange(dist, .by_group = TRUE) %>%
  nest() %>%
  mutate(ave_3 = map(.$data, moveavg_dist3)) %>%
  unnest() %>%
  select(-dist1, -ave_ratio_221)

test <- s3_tidy_moveavg3 %>%
  filter(background == 'v chr9' & (spacing == 5 | spacing == 10))
  
p_subpool3_spa_4_vchr9_5_10 <- s3_tidy_moveavg3 %>%
  filter(background == 'v chr9' & dist < 127 & (spacing == 5 | spacing == 10)) %>%
  ggplot(aes(x = dist, y = ave_ratio_22, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c('gray20', 'dodgerblue3'), name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  geom_vline(xintercept = c(78, 88.5), color = 'gray20', linetype = 2, 
             alpha = 0.5) +
  geom_vline(xintercept = c(82.5, 92.5), color = 'dodgerblue3', linetype = 2, 
             alpha = 0.5) +
  scale_y_log10(limits = c(0.1, 2)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 66, to = 126, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

p_subpool3_spa_4_vchr9_5_15 <- s3_tidy_moveavg3 %>%
  filter(background == 'v chr9' & dist < 127 & (spacing == 5 | spacing == 15)) %>%
  ggplot(aes(x = dist, y = ave_ratio_22, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c('gray20', 'firebrick3'), 
                     name = 'spacing (bp)') +
  scale_fill_manual(values = c('gray20', 'firebrick3'), 
                    name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  geom_vline(xintercept = c(78, 88.5), color = 'gray20', linetype = 2, 
             alpha = 0.5) +
  scale_y_log10(limits = c(0.1, 2)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 66, to = 126, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

p_subpool3_spa_4_vchr9_10_20 <- s3_tidy_moveavg3 %>%
  filter(background == 'v chr9' & dist < 127 & (spacing == 10 | spacing == 20)) %>%
  ggplot(aes(x = dist, y = ave_ratio_22, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c('dodgerblue3', '#55C667FF'),
                     name = 'spacing (bp)') +
  scale_fill_manual(values = c('dodgerblue3', '#55C667FF'),
                    name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  geom_vline(xintercept = c(82.5, 92.5), color = 'dodgerblue3', linetype = 2, 
             alpha = 0.5) +
  scale_y_log10(limits = c(0.1, 2)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 66, to = 126, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

ggsave('plots/p_subpool3_spa_4_vchr9_5_10.pdf', p_subpool3_spa_4_vchr9_5_10, 
       scale = 1.3, height = 1.65, width = 4.8, units = 'in')

ggsave('plots/p_subpool3_spa_4_vchr9_5_15.pdf', p_subpool3_spa_4_vchr9_5_15, 
       scale = 1.3, height = 1.65, width = 4.8, units = 'in')

ggsave('plots/p_subpool3_spa_4_vchr9_10_20.pdf', p_subpool3_spa_4_vchr9_10_20, 
       scale = 1.3, height = 1.65, width = 4.8, units = 'in')


#Plot subpool5 expression features----------------------------------------------

#Concensus expression uninduced and induced

p_s5_consnum_0_4 <- s5_untidy %>%
  filter(weak == 0) %>%
  filter(conc == 0 | conc == 4) %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  ggplot(aes(as.factor(consensus), ave_ratio_norm, fill = as.factor(conc))) +
  geom_boxplot(outlier.size = 1, size = 0.3, outlier.shape = 21,
               outlier.alpha = 1, position = position_dodge(0.75),
               show.legend = TRUE) + 
  scale_fill_manual(values = forskolin2, name = 'forskolin') +
  facet_grid(~ background) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  geom_vline(xintercept = c(1.5:6.5), alpha = 0.25) +
  panel_border(colour = 'black') +
  scale_y_log10() +
  annotation_logticks(sides = 'l') +
  ylab('Average normalized\nexpression (a.u.)') +
  xlab('Consensus sites')

save_plot('plots/p_s5_consnum_0_4.pdf', p_s5_consnum_0_4, scale = 1.3,
          base_width = 5, base_height = 2.5)

#Weak expression uninduced and induced

p_s5_weaknum_0_4 <- s5_untidy %>%
  filter(consensus == 0) %>%
  filter(conc == 0 | conc == 4) %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  ggplot(aes(as.factor(weak), ave_ratio_norm, fill = as.factor(conc))) +
  geom_boxplot(outlier.size = 1, size = 0.3, outlier.shape = 21,
               outlier.alpha = 1, position = position_dodge(0.75),
               show.legend = TRUE) + 
  scale_fill_manual(values = forskolin2, name = 'forskolin') +
  facet_grid(~ background) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  geom_vline(xintercept = c(1.5:6.5), alpha = 0.25) +
  panel_border(colour = 'black') +
  xlab('Weak sites') +
  ylab('Average normalized\nexpression (a.u.)') +
  scale_y_continuous(limits = c(0.75, 2), breaks = c(0.75, 1, 1.25, 1.5))

save_plot('plots/p_s5_weaknum_0_4.pdf', p_s5_weaknum_0_4, scale = 1.3,
          base_width = 5, base_height = 2.5)

#Consensus with weak induced expression

p_s5_num_cons_num_weak <- s5_tidy %>%
  filter(background == 's pGl4') %>%
  ggplot(aes(as.factor(consensus), ave_ratio_22_norm)) +
  geom_boxplot(aes(fill = as.factor(weak)), outlier.size = 1, size = 0.3, 
               outlier.shape = 21, outlier.alpha = 1, 
               position = position_dodge(0.75), show.legend = TRUE) +
  scale_y_log10() + 
  panel_border(colour = 'black') +
  annotation_logticks(sides = 'l') +
  scale_fill_manual(name = 'number of\nweak sites', 
                    values = cbPalette7_grad_light) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  background_grid(major = 'y', minor = 'none') + 
  geom_vline(xintercept = c(1.5:6.5), alpha = 0.25) +
  ylab('Average normalized\n expression (a.u.)') +
  xlab('Number of consensus sites')

save_plot('plots/p_s5_num_cons_num_weak.pdf', p_s5_num_cons_num_weak,
          scale = 1.3, base_height = 2.25, base_width = 5.25)

p_s5_num_cons_num_weak_allback_4 <- s5_tidy %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  ggplot(aes(as.factor(consensus), ave_ratio_22_norm)) +
  facet_grid(background ~ .) +
  geom_boxplot(aes(fill = as.factor(weak)), outlier.size = 1, size = 0.3, 
               outlier.shape = 21, outlier.alpha = 1, 
               position = position_dodge(0.75), show.legend = TRUE) +
  scale_y_log10() + 
  panel_border(colour = 'black') +
  annotation_logticks(sides = 'l') +
  scale_fill_manual(name = 'number of\nweak sites', values = cbPalette7_grad_light)  +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  background_grid(major = 'y', minor = 'none') + 
  geom_vline(xintercept = c(1.5:6.5), alpha = 0.25) +
  ylab('Average normalized\n expression (a.u.)') +
  xlab('Number of consensus sites')

save_plot('plots/p_s5_num_cons_num_weak_allback_4.pdf', 
          p_s5_num_cons_num_weak_allback_4,
          scale = 1.3, base_height = 4.75, base_width = 5.25)

p_s5_num_cons_num_weak_allback_0 <- s5_tidy %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  ggplot(aes(as.factor(consensus), ave_ratio_0_norm)) +
  facet_grid(background ~ .) +
  geom_boxplot(aes(fill = as.factor(weak)), outlier.size = 1, size = 0.3, 
               outlier.shape = 21, outlier.alpha = 1, 
               position = position_dodge(0.75), show.legend = TRUE) +
  panel_border(colour = 'black') +
  scale_fill_manual(name = 'number of\nweak sites', values = cbPalette7_grad_light) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  background_grid(major = 'y', minor = 'none') + 
  geom_vline(xintercept = c(1.5:6.5), alpha = 0.25) +
  ylab('Average normalized\n expression (a.u.)') +
  xlab('Number of consensus sites')

save_plot('plots/p_s5_num_cons_num_weak_allback_0.pdf', 
          p_s5_num_cons_num_weak_allback_0,
          scale = 1.3, base_height = 4.75, base_width = 5.25)

#Looking at individual site expression

test <- s5_untidy %>%
  filter(conc == 4 & consensus >=5 & background == 'v chr5' & weak == 0) %>%
  write_csv('s5_sitedesign.csv')

s5_single_site_exp <- s5_untidy %>%
  filter(site_combo != 'mixed' & total_sites == 1) %>%
  mutate(site1 = str_detect(site1, 
                            paste(c("consensus", 'weak'), 
                                  collapse = '|')) * -127) %>%
  mutate(site2 = str_detect(site2, 
                            paste(c("consensus", 'weak'), 
                                  collapse = '|')) * -102) %>%
  mutate(site3 = str_detect(site3, 
                            paste(c("consensus", 'weak'), 
                                  collapse = '|')) * -77) %>%
  mutate(site4 = str_detect(site4, 
                            paste(c("consensus", 'weak'), 
                                  collapse = '|')) * -52) %>%
  mutate(site5 = str_detect(site5, 
                            paste(c("consensus", 'weak'), 
                                  collapse = '|')) * -27) %>%
  mutate(site6 = str_detect(site6, 
                            paste(c("consensus", 'weak'), 
                                  collapse = '|')) * -2) %>%
  mutate(site = site1 + site2 + site3 + site4 + site5 + site6)

p_s5_single_site_exp <- s5_single_site_exp %>%
  filter(conc == 4, consensus == 1) %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  ggplot(aes(site, ave_ratio_norm)) +
  geom_point(aes(color = background), shape = 21) +
  geom_line(aes(color = background)) +
  geom_hline(yintercept = 1, linetype = 2, alpha = 0.5) +
  scale_color_manual(name = 'background', values = cbPalette3) +
  ylab('Average normalized\nexpression (a.u.)') +
  xlab('Distance along background (bp)')

save_plot('plots/p_s5_single_site_exp.pdf', p_s5_single_site_exp,
          scale = 1.3, base_height = 2, base_width = 4)

#Look at effect of adding weak sites to response curves

s5_binom <- s5_untidy %>%
  mutate(ave_ratio_norm = log10(ave_ratio_norm)) %>%
  mutate(ave_ratio_norm_max = max(ave_ratio_norm)) %>%
  mutate(ave_ratio_norm = ave_ratio_norm/ave_ratio_norm_max)

test <- s5_binom %>%
  mutate(conc = conc + 0.005) %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  filter(consensus > 1 & weak == 0) %>%
  ggplot(aes(conc, ave_ratio_norm, color = as.factor(consensus))) +
  facet_grid(background ~ .) +
  scale_x_log10() +
  annotation_logticks(sides = 'b') +
  geom_point(alpha = 0.75) +
  geom_smooth(alpha = 0.5, method = 'loess') +
  scale_color_viridis(discrete = TRUE)

#Induction analysis-------------------------------------------------------------

s5_typecount_siteloc <- function(df) {
  df <- df %>%
    mutate(ave_ratio_22 = (ratio_22A + ratio_22B)/2)
  csite1 <- df %>%
    group_by(most_common, background, consensus, weak, 
             site1, site2, site3, site4, site5, site6, ave_ratio_22) %>%
    count(site1, site1 == 'consensus') %>%
    filter(site1 == 'consensus') %>%
    mutate(csite = '-191') %>%
    select(most_common, background, weak, consensus, site1, site2, site3, site4, 
           site5, site6, csite, ave_ratio_22) %>%
    ungroup()
  csite2 <- df %>%
    group_by(most_common, background, consensus, weak, 
             site1, site2, site3, site4, site5, site6, ave_ratio_22) %>%
    count(site2, site2 == 'consensus') %>%
    filter(site2 == 'consensus') %>%
    mutate(csite = '-166') %>%
    select(most_common, background, weak, consensus, site1, site2, site3, site4, 
           site5, site6, csite, ave_ratio_22) %>%
    ungroup()
  csite3 <- df %>%
    group_by(most_common, background, consensus, weak, 
             site1, site2, site3, site4, site5, site6, ave_ratio_22) %>%
    count(site3, site3 == 'consensus') %>%
    filter(site3 == 'consensus') %>%
    mutate(csite = '-141') %>%
    select(most_common, background, weak, consensus, site1, site2, site3, site4, 
           site5, site6, csite, ave_ratio_22) %>%
    ungroup()
  csite4 <- df %>%
    group_by(most_common, background, consensus, weak, 
             site1, site2, site3, site4, site5, site6, ave_ratio_22) %>%
    count(site4, site4 == 'consensus') %>%
    filter(site4 == 'consensus') %>%
    mutate(csite = '-116') %>%
    select(most_common, background, weak, consensus, site1, site2, site3, site4, 
           site5, site6, csite, ave_ratio_22) %>%
    ungroup()
  csite5 <- df %>%
    group_by(most_common, background, consensus, weak, 
             site1, site2, site3, site4, site5, site6, ave_ratio_22) %>%
    count(site5, site5 == 'consensus') %>%
    filter(site5 == 'consensus') %>%
    mutate(csite = '-91') %>%
    select(most_common, background, weak, consensus, site1, site2, site3, site4, 
           site5, site6, csite, ave_ratio_22) %>%
    ungroup()
  csite6 <- df %>%
    group_by(most_common, background, consensus, weak, 
             site1, site2, site3, site4, site5, site6, ave_ratio_22) %>%
    count(site6, site6 == 'consensus') %>%
    filter(site6 == 'consensus') %>%
    mutate(csite = '-66') %>%
    select(most_common, background, weak, consensus, site1, site2, site3, site4, 
           site5, site6, csite, ave_ratio_22) %>%
    ungroup()
  nocsites <- df %>%
    filter(consensus == 0) %>%
    mutate(csite = 'none') %>%
    select(most_common, background, weak, consensus, site1, site2, site3, site4, 
           site5, site6, csite, ave_ratio_22)
  csites <- rbind(csite1, csite2, csite3, csite4, csite5, csite6, nocsites)
  wsite1 <- df %>%
    group_by(most_common, background, consensus, weak, 
             site1, site2, site3, site4, site5, site6, ave_ratio_22) %>%
    count(site1, site1 == 'weak') %>%
    filter(site1 == 'weak') %>%
    mutate(wsite = -191) %>%
    select(most_common, background, weak, consensus, site1, site2, site3, site4, 
           site5, site6, wsite, ave_ratio_22) %>%
    ungroup()
  wsite2 <- df %>%
    group_by(most_common, background, consensus, weak, 
             site1, site2, site3, site4, site5, site6, ave_ratio_22) %>%
    count(site2, site2 == 'weak') %>%
    filter(site2 == 'weak') %>%
    mutate(wsite = -166) %>%
    select(most_common, background, weak, consensus, site1, site2, site3, site4, 
           site5, site6, wsite, ave_ratio_22) %>%
    ungroup()
  wsite3 <- df %>%
    group_by(most_common, background, consensus, weak, 
             site1, site2, site3, site4, site5, site6, ave_ratio_22) %>%
    count(site3, site3 == 'weak') %>%
    filter(site3 == 'weak') %>%
    mutate(wsite = -141) %>%
    select(most_common, background, weak, consensus, site1, site2, site3, site4, 
           site5, site6, wsite, ave_ratio_22) %>%
    ungroup()
  wsite4 <- df %>%
    group_by(most_common, background, consensus, weak, 
             site1, site2, site3, site4, site5, site6, ave_ratio_22) %>%
    count(site4, site4 == 'weak') %>%
    filter(site4 == 'weak') %>%
    mutate(wsite = -116) %>%
    select(most_common, background, weak, consensus, site1, site2, site3, site4, 
           site5, site6, wsite, ave_ratio_22) %>%
    ungroup()
  wsite5 <- df %>%
    group_by(most_common, background, consensus, weak, 
             site1, site2, site3, site4, site5, site6, ave_ratio_22) %>%
    count(site5, site5 == 'weak') %>%
    filter(site5 == 'weak') %>%
    mutate(wsite = -91) %>%
    select(most_common, background, weak, consensus, site1, site2, site3, site4, 
           site5, site6, wsite, ave_ratio_22) %>%
    ungroup()
  wsite6 <- df %>%
    group_by(most_common, background, consensus, weak, 
             site1, site2, site3, site4, site5, site6, ave_ratio_22) %>%
    count(site6, site6 == 'weak') %>%
    filter(site6 == 'weak') %>%
    mutate(wsite = -66) %>%
    select(most_common, background, weak, consensus, site1, site2, site3, site4, 
           site5, site6, wsite, ave_ratio_22) %>%
    ungroup()
  nowsites <- df %>%
    filter(weak == 0) %>%
    mutate(wsite = 'none') %>%
    select(most_common, background, weak, consensus, site1, site2, site3, site4, 
           site5, site6, wsite, ave_ratio_22)
  wsites <- rbind(wsite1, wsite2, wsite3, wsite4, wsite5, wsite6, nowsites)
  c_w_sites <- inner_join(csites, wsites,
                          by = c('most_common', 'background', 'weak', 
                                 'consensus', 'site1', 'site2', 'site3', 
                                 'site4', 'site5', 'site6', 'ave_ratio_22'))
  return(c_w_sites)
}

s5_cons_weak_csite_wsite_ind <- s5_typecount_siteloc(s5_tidy) %>%
  mutate(csite = factor(csite, levels = c('none', '-191', '-166', '-141', '-116',
                                          '-91', '-66'))) %>%
  filter(weak == 0, consensus == 1) %>%
  ggplot(aes(csite, ave_ratio_22)) +
  geom_boxplot() +
  scale_y_log10() +
  annotation_logticks(sides = 'l')

test <- s5_cons_weak_csite_wsite_ind %>%
  filter(background == 'v chr9') %>%
  mutate(csite = factor(csite, levels = c('none', '-127', '-102', '-77', '-52',
                                          '-27', '-2'))) %>%
  ggplot(aes(csite, induction, color = as.factor(consensus))) +
  geom_boxplot()

p_s5_cons_weak_csite_wsite_med_ind_vchr9 <- s5_cons_weak_csite_wsite_med_ind %>%
  filter(background == 'v chr9') %>%
  ungroup() %>%
  mutate(csite = factor(csite, levels = c('none', '-127', '-102', '-77', '-52',
                                          '-27', '-2'))) %>%
  mutate(wsite = factor(wsite, levels = c('none', '-127', '-102', '-77', '-52',
                                          '-27', '-2'))) %>%
  mutate(weak = factor(weak, levels = c(6, 5, 4, 3, 2, 1, 0))) %>%
  ggplot(aes(csite, wsite)) +
  geom_tile(aes(fill = log10(med_induction))) +
  facet_grid(weak ~ consensus) +
  scale_fill_viridis(name = 'log10(med(induction))') +
  panel_border(colour = 'black') +
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        panel.spacing = unit(0.10, "lines"))

save_plot('plots/p_s5_cons_weak_csite_wsite_med_ind_vchr9.pdf', 
          p_s5_cons_weak_csite_wsite_med_ind_vchr9, scale = 1.3,
          base_width = 5, base_height = 3)

p_s5_cons_weak_csite_wsite_med_ind_spGl4 <- s5_cons_weak_csite_wsite_med_ind %>%
  filter(background == 's pGl4') %>%
  ungroup() %>%
  mutate(csite = factor(csite, levels = c('none', '-127', '-102', '-77', '-52',
                                          '-27', '-2'))) %>%
  mutate(wsite = factor(wsite, levels = c('none', '-127', '-102', '-77', '-52',
                                          '-27', '-2'))) %>%
  mutate(weak = factor(weak, levels = c(6, 5, 4, 3, 2, 1, 0))) %>%
  ggplot(aes(csite, wsite)) +
  geom_tile(aes(fill = log10(med_induction))) +
  facet_grid(weak ~ consensus) +
  scale_fill_viridis(name = 'log10(med(induction))', limits = c(-0.05, 1.5)) +
  panel_border(colour = 'black') +
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        panel.spacing = unit(0.10, "lines"))

save_plot('plots/p_s5_cons_weak_csite_wsite_med_ind_spGl4.pdf', 
          p_s5_cons_weak_csite_wsite_med_ind_spGl4, scale = 1.3,
          base_width = 5, base_height = 3)

p_s5_cons_weak_csite_wsite_med_ind_vchr5 <- s5_cons_weak_csite_wsite_med_ind %>%
  filter(background == 'v chr5') %>%
  ungroup() %>%
  mutate(csite = factor(csite, levels = c('none', '-127', '-102', '-77', '-52',
                                          '-27', '-2'))) %>%
  mutate(wsite = factor(wsite, levels = c('none', '-127', '-102', '-77', '-52',
                                          '-27', '-2'))) %>%
  mutate(weak = factor(weak, levels = c(6, 5, 4, 3, 2, 1, 0))) %>%
  ggplot(aes(csite, wsite)) +
  geom_tile(aes(fill = log10(med_induction))) +
  facet_grid(weak ~ consensus) +
  scale_fill_viridis(name = 'log10(med(induction))', limits = c(-0.05, 1.5)) +
  panel_border(colour = 'black') +
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        panel.spacing = unit(0.10, "lines"))

save_plot('plots/p_s5_cons_weak_csite_wsite_med_ind_vchr5.pdf', 
          p_s5_cons_weak_csite_wsite_med_ind_vchr5, scale = 1.3,
          base_width = 5, base_height = 3)


ind_top50 <- function(df) {
  vchr9 <- df %>%
    filter(background == 'v chr9') %>%
    arrange(desc(induction)) %>%
    slice(1:50)
  spgl4 <- df %>%
    filter(background == 's pGl4') %>%
    arrange(desc(induction)) %>%
    slice(1:50)
  vchr5 <- df %>%
    filter(background == 'v chr5') %>%
    arrange(desc(induction)) %>%
    slice(1:50)
  back_bind <- rbind(vchr9, spgl4, vchr5)
  site1 <- back_bind %>%
    group_by(background) %>%
    count(site1, site1 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site1) %>%
    select(-3) %>%
    mutate(site = 1)
  site2 <- back_bind %>%
    group_by(background) %>%
    count(site2, site2 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site2) %>%
    select(-3) %>%
    mutate(site = 2)
  site3 <- back_bind %>%
    group_by(background) %>%
    count(site3, site3 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site3) %>%
    select(-3) %>%
    mutate(site = 3)
  site4 <- back_bind %>%
    group_by(background) %>%
    count(site4, site4 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site4) %>%
    select(-3) %>%
    mutate(site = 4)
  site5 <- back_bind %>%
    group_by(background) %>%
    count(site5, site5 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site5) %>%
    select(-3) %>%
    mutate(site = 5)
  site6 <- back_bind %>%
    group_by(background) %>%
    count(site6, site6 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site6) %>%
    select(-3) %>%
    mutate(site = 6)
  site_join <- bind_rows(site1, site2, site3, site4, site5, site6) %>%
    ungroup()
  return(site_join)
}
  
s5_ind50_site_count <- ind_top50(s5_tidy) %>%
  mutate(type = factor(type, levels = c('nosite', 'weak', 'consensus'))) %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5')))
  
p_s5_ind50_site_count <- ggplot(s5_ind50_site_count,
                                aes(as.factor(site), counts, fill = type)) +
  facet_grid(~ background) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_viridis(discrete = TRUE) + 
  xlab('Site position') + 
  panel_border()


#BC analysis--------------------------------------------------------------------

#Determining subpool proportions in DNA

DNA_nonzero <- bc_join_DNA %>%
  filter(num_reads > 0) %>%
  filter(subpool == 'subpool5')


#Number of BC's/variant in each sample

p_BC_per_variant <- ggplot(rep_0_22_A_B, aes(x = "", y = NULL)) +
  geom_boxplot(aes(x = "DNA", y = barcodes_DNA, color = subpool)) +
  geom_boxplot(aes(x = "R0.00A", y = barcodes_RNA_0A, color = subpool)) + 
  geom_boxplot(aes(x = "R0.00B", y = barcodes_RNA_0B, color = subpool)) + 
  geom_boxplot(aes(x = "R0.03A", y = barcodes_RNA_2_5A, color = subpool)) + 
  geom_boxplot(aes(x = "R0.03B", y = barcodes_RNA_2_5B, color = subpool)) + 
  geom_boxplot(aes(x = "R0.06A", y = barcodes_RNA_2_4A, color = subpool)) +
  geom_boxplot(aes(x = "R0.06B", y = barcodes_RNA_2_4B, color = subpool)) + 
  geom_boxplot(aes(x = "R0.13A", y = barcodes_RNA_2_3A, color = subpool)) +
  geom_boxplot(aes(x = "R0.13B", y = barcodes_RNA_2_3B, color = subpool)) + 
  geom_boxplot(aes(x = "R0.25A", y = barcodes_RNA_2_2A, color = subpool)) +
  geom_boxplot(aes(x = "R0.25B", y = barcodes_RNA_2_2B, color = subpool)) + 
  geom_boxplot(aes(x = "R0.50A", y = barcodes_RNA_2_1A, color = subpool)) +
  geom_boxplot(aes(x = "R0.50B", y = barcodes_RNA_2_1B, color = subpool)) + 
  geom_boxplot(aes(x = "R1.00A", y = barcodes_RNA_20A, color = subpool)) +
  geom_boxplot(aes(x = "R1.00B", y = barcodes_RNA_20B, color = subpool)) + 
  geom_boxplot(aes(x = "R4.00A", y = barcodes_RNA_22A, color = subpool)) +
  geom_boxplot(aes(x = "R4.00B", y = barcodes_RNA_22B, color = subpool)) + 
  xlab("") +
  ylab("Barcodes per variant") +
  background_grid(major = c('y'), minor = c('y')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot('plots/BC_per_variant.png',
          p_BC_per_variant, scale = 2.8)


#replicate plots----------------------------------------------------------------

#plot replicates for summed variant expression at 4 µM for fig. 1

p_fig1_trans_rep <- rep_0_22_A_B %>%
  ggplot(aes(ratio_22A, ratio_22B)) +
  geom_point(alpha = 0.1, size = 0.75) +
  geom_point(data = filter(rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             fill = 'orange', shape = 21, size = 2.25) + 
  geom_point(data = filter(rep_0_22_A_B, 
                           name == 'pGL4.29 Promega 1-63 + 1-87'), 
             fill = 'red', shape = 21, size = 2.25) +
  annotation_logticks(scaled = TRUE) +
  xlab("Expression (a.u.) replicate 1") +
  ylab("Expression (a.u.) replicate 2") +
  scale_x_log10(limits = c(0.01, 20), breaks = c(0.01, 0.1, 1, 10)) + 
  scale_y_log10(limits = c(0.01, 20), breaks = c(0.01, 0.1, 1, 10)) +
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.line.y = element_line(), panel.spacing.x=unit(1, "lines")) 

ggsave('plots/p_fig1_trans_rep.pdf', p_fig1_trans_rep, scale = 1.3,
       width = 2.1, height = 1.8, units = 'in')

trans_4_pearsons <- tibble(
  sample = c('all', 'subpool3', 'subpool5'),
  pearsons = c(round(cor(log10_rep_0_22_A_B$ratio_22A, 
                         log10_rep_0_22_A_B$ratio_22B, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(filter(log10_rep_0_22_A_B, 
                                subpool == 'subpool3')$ratio_22A,
                         filter(log10_rep_0_22_A_B, 
                                subpool == 'subpool3')$ratio_22B,
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(filter(log10_rep_0_22_A_B, 
                                subpool == 'subpool5')$ratio_22A,
                         filter(log10_rep_0_22_A_B, 
                                subpool == 'subpool5')$ratio_22B,
                         use = "pairwise.complete.obs", method = "pearson"), 
                     3)))

write_csv(trans_4_pearsons, 'trans_4_pearsons.csv')


#plot all concentrations faceted for titration figure

var_conc_ratio <- function(df) {
  df_0 <- df %>%
    mutate(ave_barcode_0 = (barcodes_RNA_0A + barcodes_RNA_0B)/2) %>%
    select(subpool, name, most_common, background, ratio_0A_norm, 
           ratio_0B_norm) %>%
    mutate(conc = 0) %>%
    rename(ratio_A_norm = ratio_0A_norm) %>%
    rename(ratio_B_norm = ratio_0B_norm)
  df_2_5 <- df %>%
    mutate(ave_barcode_2_5 = (barcodes_RNA_2_5A + barcodes_RNA_2_5B)/2) %>%
    select(subpool, name, most_common, background, ratio_2_5A_norm, 
           ratio_2_5B_norm) %>%
    mutate(conc = 2^-5) %>%
    rename(ratio_A_norm = ratio_2_5A_norm) %>%
    rename(ratio_B_norm = ratio_2_5B_norm)
  df_2_4 <- df %>%
    mutate(ave_barcode_2_4 = (barcodes_RNA_2_4A + barcodes_RNA_2_4B)/2) %>%
    select(subpool, name, most_common, background, ratio_2_4A_norm, ratio_2_4B_norm) %>%
    mutate(conc = 2^-4) %>%
    rename(ratio_A_norm = ratio_2_4A_norm) %>%
    rename(ratio_B_norm = ratio_2_4B_norm)
  df_2_3 <- df %>%
    mutate(ave_barcode_2_3 = (barcodes_RNA_2_3A + barcodes_RNA_2_3B)/2) %>%
    select(subpool, name, most_common, background, ratio_2_3A_norm, ratio_2_3B_norm) %>%
    mutate(conc = 2^-3) %>%
    rename(ratio_A_norm = ratio_2_3A_norm) %>%
    rename(ratio_B_norm = ratio_2_3B_norm)
  df_2_2 <- df %>%
    mutate(ave_barcode_2_2 = (barcodes_RNA_2_2A + barcodes_RNA_2_2B)/2) %>%
    select(subpool, name, most_common, background, ratio_2_2A_norm, ratio_2_2B_norm) %>%
    mutate(conc = 2^-2) %>%
    rename(ratio_A_norm = ratio_2_2A_norm) %>%
    rename(ratio_B_norm = ratio_2_2B_norm)
  df_2_1 <- df %>%
    mutate(ave_barcode_2_1 = (barcodes_RNA_2_1A + barcodes_RNA_2_1B)/2) %>%
    select(subpool, name, most_common, background, ratio_2_1A_norm, ratio_2_1B_norm) %>%
    mutate(conc = 2^-1) %>%
    rename(ratio_A_norm = ratio_2_1A_norm) %>%
    rename(ratio_B_norm = ratio_2_1B_norm)
  df_20 <- df %>%
    mutate(ave_barcode_20 = (barcodes_RNA_20A + barcodes_RNA_20B)/2) %>%
    select(subpool, name, most_common, background, ratio_20A_norm, ratio_20B_norm) %>%
    mutate(conc = 2^0) %>%
    rename(ratio_A_norm = ratio_20A_norm) %>%
    rename(ratio_B_norm = ratio_20B_norm)
  df_22 <- df %>%
    mutate(ave_barcode_22 = (barcodes_RNA_22A + barcodes_RNA_22B)/2) %>%
    select(subpool, name, most_common, background, ratio_22A_norm, ratio_22B_norm) %>%
    mutate(conc = 2^2) %>%
    rename(ratio_A_norm = ratio_22A_norm) %>%
    rename(ratio_B_norm = ratio_22B_norm)
  df_0_22 <- rbind(df_0, df_2_5, df_2_4, df_2_3, df_2_2, df_2_1, df_20, df_22)
  return(df_0_22)
}

rep_norm_ratio_conc <- var_conc_ratio(trans_back_norm_pc_spGl4)

p_var_rep_all <- ggplot(rep_norm_ratio_conc, aes(ratio_A_norm, ratio_B_norm)) +
  facet_rep_wrap(~ conc, nrow = 4, ncol = 2) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_point(data = filter(rep_norm_ratio_conc, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             fill = 'orange', shape = 21, size = 1.75) + 
  geom_point(data = filter(rep_norm_ratio_conc, name == 'pGL4.29 Promega 1-63 + 1-87_s pGl4'), 
             fill = 'red', shape = 21, size = 1.75) +
  annotation_logticks(scaled = TRUE) +
  xlab("Normalized expression (a.u.)\nreplicate 1") +
  ylab("Normalized expression (a.u.) replicate 2") +
  scale_x_log10(limits = c(0.5, 175), breaks = c(1, 10, 100)) + 
  scale_y_log10(limits = c(0.5, 175), breaks = c(1, 10, 100)) +
  background_grid(major = 'xy', minor = 'none') + 
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.line.y = element_line(), panel.spacing.x=unit(1, "lines"))

save_plot('plots/var_rep_all.pdf', p_var_rep_all, scale = 1.3, 
          base_width = 2.5, base_height = 5)

log10_trans_back_norm_pc_spGl4 <- var_log10(trans_back_norm_pc_spGl4)

pearsons_conc <- tibble(
  conc = c(0, 2^-5, 2^-4, 2^-3, 2^-2, 2^-1, 2^0, 2^2),
  pearsons = c(round(cor(log10_trans_back_norm_pc_spGl4$ratio_0A_norm, 
                         log10_trans_back_norm_pc_spGl4$ratio_0B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_trans_back_norm_pc_spGl4$ratio_2_5A_norm, 
                         log10_trans_back_norm_pc_spGl4$ratio_2_5B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_trans_back_norm_pc_spGl4$ratio_2_4A_norm, 
                         log10_trans_back_norm_pc_spGl4$ratio_2_4B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_trans_back_norm_pc_spGl4$ratio_2_3A_norm, 
                         log10_trans_back_norm_pc_spGl4$ratio_2_3B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_trans_back_norm_pc_spGl4$ratio_2_2A_norm, 
                         log10_trans_back_norm_pc_spGl4$ratio_2_2B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_trans_back_norm_pc_spGl4$ratio_2_1A_norm, 
                         log10_trans_back_norm_pc_spGl4$ratio_2_1B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_trans_back_norm_pc_spGl4$ratio_20A_norm, 
                         log10_trans_back_norm_pc_spGl4$ratio_20B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_trans_back_norm_pc_spGl4$ratio_22A_norm, 
                         log10_trans_back_norm_pc_spGl4$ratio_22B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3))
  )

write_csv(pearsons_conc, 'pearsons_conc.csv')


#Comparison to integrated-------------------------------------------------------

int_rep_1_2 <- read_tsv('../20171129_intLib/rep_1_2.txt') %>%
  mutate(ave_med_ratio = (med_ratio_br1 + med_ratio_br2)/2)

int_trans <- rep_0_22_A_B %>%
  select(subpool, name, most_common, barcodes_DNA, ratio_22A, barcodes_RNA_22A, 
         ratio_22B, barcodes_RNA_22B) %>%
  mutate(ave_ratio_22 = (ratio_22A + ratio_22B)/2) %>%
  inner_join(int_rep_1_2, by = c('subpool', 'name', 'most_common'))

MPRA_ave <- int_trans %>%
  ungroup() %>%
  filter(subpool != 'control') %>%
  mutate(
        name = gsub('Smith R. Vista chr9:83712599-83712766', 'v chr9', name),
        name = gsub('Vista Chr5:88673410-88674494', 'v chr5', name),
        name = gsub('scramble pGL4.29 Promega 1-63 \\+ 1-87', 's pGl4', name)
        ) %>%
  mutate(background = name) %>%
  mutate(background = str_sub(background, 
                              nchar(background)-5,
                              nchar(background))) %>%
  select(subpool, name, background, barcodes_RNA_br1, barcodes_RNA_br2, 
         med_ratio_br1, med_ratio_br2, ave_med_ratio, barcodes_RNA_22A, 
         barcodes_RNA_22B, ratio_22A, ratio_22B, ave_ratio_22) %>%
  mutate(integrated = (barcodes_RNA_br1 + barcodes_RNA_br2)/2) %>%
  mutate(episomal = (barcodes_RNA_22A + barcodes_RNA_22B)/2) %>%
  select(-barcodes_RNA_br1, -barcodes_RNA_br2, -barcodes_RNA_22A, 
         -barcodes_RNA_22B) %>%
  gather(integrated, episomal, key = 'MPRA', value = 'barcodes') %>%
  mutate(integrated = ave_med_ratio) %>%
  mutate(episomal = ave_ratio_22) %>%
  gather(integrated, episomal, key = 'MPRA2', value = 'ave_ratio') %>%
  filter((MPRA == 'integrated' & MPRA2 == 'integrated') | (MPRA == 'episomal' & MPRA2 == 'episomal')) %>%
  select(-MPRA2)

int_trans_log10 <- var_log10(int_trans)

int_trans_pearsons <- tibble(
  sample = c('all', 'subpool3', 'subpool5'),
  pearsons = c(round(cor(int_trans_log10$ave_med_ratio, 
                         int_trans_log10$ave_ratio_22, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(filter(int_trans_log10, 
                                subpool == 'subpool3')$ave_med_ratio,
                         filter(int_trans_log10, 
                                subpool == 'subpool3')$ave_ratio_22,
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(filter(int_trans_log10, 
                                subpool == 'subpool5')$ave_med_ratio,
                         filter(int_trans_log10, 
                                subpool == 'subpool5')$ave_ratio_22,
                         use = "pairwise.complete.obs", method = "pearson"), 
                     3)))

write_csv(int_trans_pearsons, 'int_trans_pearsons.csv')

#Comparison of subpool3 between integrated and episomal-------------------------

p_int_trans_ave_med_rep_sp3 <- int_trans %>%
  filter(subpool == 'subpool3') %>%
  ggplot(aes(ave_med_ratio, ave_ratio_22)) +
  geom_point(alpha = 0.2, size = 1) +
  xlab("Average integrated expression (a.u.)") +
  ylab("Average episomal expression (a.u.)") +
  scale_x_log10(limits = c(0.01, 5)) + 
  scale_y_log10(limits = c(0.1, 5)) +
  annotation_logticks(scaled = TRUE)

ggsave('plots/p_int_trans_ave_med_rep_sp3.pdf', p_int_trans_ave_med_rep_sp3, 
       scale = 1.3, width = 3, height = 2.7, units = 'in')


#Compare subpool3 features between integrated and transient

s3_int_trans <- MPRA_ave %>%
  filter(subpool == 'subpool3') %>%
  subpool3()


library(caTools)

int_epi_moveavg_dist3 <- function(df) {
  df <- df %>%
    mutate(ave_3 = runmean(ave_ratio, 3, alg = 'R', endrule = 'NA'))
}

s3_int_trans_moveavg3 <- s3_int_trans %>%
  select(spacing, dist, background, MPRA, ave_ratio) %>%
  group_by(background, spacing, MPRA) %>%
  arrange(dist, .by_group = TRUE) %>%
  nest() %>%
  mutate(ave_3 = map(.$data, int_epi_moveavg_dist3)) %>%
  unnest() %>%
  select(-dist1, -ave_ratio1)

test <- s3_int_trans_moveavg3 %>%
  filter(background == 's pGl4' & MPRA == 'episomal' & spacing == 20)


#Subpool3 despite noise seems to have same distance effects across background

p_space_dist_int_trans <- s3_int_trans_moveavg3 %>%
  filter(background != 'v chr9') %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  ggplot(aes(x = dist, y = ave_ratio, color = MPRA)) + 
  geom_point(alpha = 0.5, size = 1.2) +
  geom_line(aes(y = ave_3), size = 0.4) +
  facet_grid(spacing ~ background) +
  scale_color_manual(values = c('#3CBB75FF', 'gray35')) +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10() +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'x', colour.major = 'grey90',
                  colour.minor = 'grey95') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 70, to = 190, by = 20)) +
  theme(legend.position = 'right',
        strip.background = element_rect(colour="black", fill="white"))

ggsave('plots/p_space_dist_int_trans.pdf', p_space_dist_int_trans, 
          scale = 1.3, width = 8.5, height = 6.5, units = 'in')

p_space15_dist_int_trans <- s3_int_trans_moveavg3 %>%
  filter(background != 'v chr9' & spacing == 15) %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  ggplot(aes(x = dist, y = ave_ratio, color = MPRA)) + 
  geom_point(alpha = 0.5, size = 1.2) +
  geom_line(aes(y = ave_3), size = 0.4) +
  facet_grid(background ~ .) +
  scale_color_manual(values = c('#3CBB75FF', 'gray35')) +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10() +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'x', colour.major = 'grey90',
                  colour.minor = 'grey95') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 70, to = 190, by = 20)) +
  theme(legend.position = 'right',
        strip.background = element_rect(colour="black", fill="white"))


#distance effects in episomal and integrated

s3_int_trans_bin20bp <- s3_int_trans %>%
  filter(dist <= 176) %>%
  mutate(bin = cut(dist, seq(from = 66, to = 176, by = 22),
                   labels = c('67-88', '89-110', '111-132', '133-154',
                              '155-176')))

p_s3_dist_int_bin20bp <- s3_int_trans_bin20bp %>%
  mutate(background = factor(background, levels = c('s pGl4', 'v chr5'))) %>%
  filter(spacing != 0 & spacing != 70 & background != 'v chr9' & MPRA == 'integrated') %>%
  ggplot(aes(bin, ave_ratio)) +
  geom_signif(comparisons = list(c('67-88', '89-110')), y_position = log10(5)) +
  geom_signif(comparisons = list(c('67-88', '111-132')), y_position = log10(15)) +
  facet_grid(. ~ background) +
  geom_jitter(aes(color = as.factor(spacing)), 
              position=position_jitter(width=0.3, height=0), alpha = 0.75,
              size = 0.75) +
  geom_boxplot(outlier.shape=NA, size = 0.3, position = position_dodge(1),
               show.legend = FALSE, alpha = 0) +
  scale_color_manual(values = spacing_5_20_palette, name = 'spacing (bp)') +
  theme(legend.position = 'right', axis.ticks.x = element_blank(), 
        strip.background = element_rect(colour="black", fill="white"),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_log10(limits = c(0.015, 30)) +
  annotation_logticks(sides = 'l') +
  panel_border(colour = 'black') +
  ylab('Average expression (a.u.)') +
  xlab('Distance to minimal promoter (bp)')

p_s3_dist_trans_bin20bp <- s3_int_trans_bin20bp %>%
  mutate(background = factor(background, levels = c('s pGl4', 'v chr5'))) %>%
  filter(spacing != 0 & spacing != 70 & background != 'v chr9' & MPRA == 'episomal') %>%
  ggplot(aes(bin, ave_ratio)) +
  geom_signif(comparisons = list(c('67-88', '89-110')), y_position = log10(6)) +
  geom_signif(comparisons = list(c('67-88', '111-132')), y_position = log10(15)) +
  facet_grid(. ~ background) +
  geom_jitter(aes(color = as.factor(spacing)), 
              position=position_jitter(width=0.3, height=0), alpha = 0.75,
              size = 0.75) +
  geom_boxplot(outlier.shape=NA, size = 0.3, position = position_dodge(1),
               show.legend = FALSE, alpha = 0) +
  scale_color_manual(values = spacing_5_20_palette, name = 'spacing (bp)') +
  theme(legend.position = 'right', axis.ticks.x = element_blank(), 
        strip.background = element_rect(colour="black", fill="white"),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_log10(limits = c(0.1, 30)) +
  annotation_logticks(sides = 'l') +
  panel_border(colour = 'black') +
  ylab('Average expression (a.u.)') +
  xlab('Distance to minimal promoter (bp)')

ggsave('plots/p_s3_dist_int_bin20bp.pdf', p_s3_dist_int_bin20bp,
       scale = 1.3, width = 5, height = 2.5, unit = 'in')

ggsave('plots/p_s3_dist_trans_bin20bp.pdf', p_s3_dist_trans_bin20bp,
       scale = 1.3, width = 5, height = 2.5, unit = 'in')


#spacing effects in integrated and episomal backgrounds

#fix axis limits between the graphs above and the graphs below

p_s3_space_trans <- s3_int_trans %>%
  filter(dist <= 110 & background != 'v chr9' & MPRA == 'episomal') %>%
  mutate(background = factor(background, levels = c('s pGl4', 'v chr5'))) %>%
  ggplot(aes(as.factor(spacing), ave_ratio)) +
  facet_grid(. ~ background) +
  geom_signif(comparisons = list(c('5', '10')), y_position = log10(3)) +
  geom_signif(comparisons = list(c('5', '15')), y_position = log10(7)) +
  geom_signif(comparisons = list(c('5', '20')), y_position = log10(15)) +
  geom_boxplot(size = 0.3, position = position_dodge(1)) +
  theme(axis.ticks.x = element_blank(), 
        strip.background = element_rect(colour="black", fill="white")) + 
  scale_y_log10(limits = c(0.1, 30)) +
  annotation_logticks(sides = 'l') +
  panel_border(colour = 'black') +
  ylab('Average expression (a.u.)') +
  xlab('Spacing (bp)')

p_s3_space_int <- s3_int_trans %>%
  filter(dist <= 110 & background != 'v chr9' & MPRA == 'integrated') %>%
  mutate(background = factor(background, levels = c('s pGl4', 'v chr5'))) %>%
  ggplot(aes(as.factor(spacing), ave_ratio)) +
  facet_grid(. ~ background) +
  geom_signif(comparisons = list(c('5', '10')), y_position = log10(2)) +
  geom_signif(comparisons = list(c('5', '15')), y_position = log10(6)) +
  geom_signif(comparisons = list(c('5', '20')), y_position = log10(15)) +
  geom_boxplot(size = 0.3, position = position_dodge(1)) +
  theme(axis.ticks.x = element_blank(), 
        strip.background = element_rect(colour="black", fill="white")) + 
  scale_y_log10(limits = c(0.010, 30)) +
  annotation_logticks(sides = 'l') +
  panel_border(colour = 'black') +
  ylab('Average expression (a.u.)') +
  xlab('Spacing (bp)')

ggsave('plots/p_s3_space_trans.pdf', p_s3_space_trans, scale = 1.3,
       width = 4, height = 2, units = 'in')

ggsave('plots/p_s3_space_int.pdf', p_s3_space_int, scale = 1.3,
       width = 4, height = 2, units = 'in')
  
  
#spacing offsets in integrated and episomal backgrounds

p_subpool3_spa_spgl4_trans_5_10 <- s3_int_trans_moveavg3 %>%
  filter(background == 's pGl4' & (spacing == 5 | spacing == 10) & MPRA == 'episomal') %>%
  ggplot(aes(x = dist, y = ave_ratio, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_vline(xintercept = c(94, 104), color = 'gray20', linetype = 2, 
             alpha = 0.5) +
  geom_vline(xintercept = c(99, 109.5), color = 'dodgerblue3', linetype = 2, 
             alpha = 0.5) +
  scale_color_manual(values = c('gray20', 'dodgerblue3'), name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10(limits = c(0.1, 0.7)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

p_subpool3_spa_spgl4_trans_5_15 <- s3_int_trans_moveavg3 %>%
  filter(background == 's pGl4' & (spacing == 5 | spacing == 15) & MPRA == 'episomal') %>%
  ggplot(aes(x = dist, y = ave_ratio, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_vline(xintercept = c(94, 104), color = 'gray20', linetype = 2, 
             alpha = 0.5) +
  scale_color_manual(values = c('gray20', 'firebrick3'), 
                     name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10(limits = c(0.1, 0.7)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

p_subpool3_spa_spgl4_trans_10_20 <- s3_int_trans_moveavg3 %>%
  filter(background == 's pGl4' & (spacing == 10 | spacing == 20) & MPRA == 'episomal') %>%
  ggplot(aes(x = dist, y = ave_ratio, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_vline(xintercept = c(99, 109.5), color = 'dodgerblue3', linetype = 2, 
             alpha = 0.5) +
  scale_color_manual(values = c('dodgerblue3', '#55C667FF'),
                     name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10(limits = c(0.1, 0.7)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10),
                     limits = c(67, 191)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

ggsave('plots/p_subpool3_spa_spgl4_trans_5_10.pdf', p_subpool3_spa_spgl4_trans_5_10, 
       scale = 1.3, height = 1.4, width = 6.5, units = 'in')

ggsave('plots/p_subpool3_spa_spgl4_trans_5_15.pdf', p_subpool3_spa_spgl4_trans_5_15, 
       scale = 1.3, height = 1.4, width = 6.5, units = 'in')

ggsave('plots/p_subpool3_spa_spgl4_trans_10_20.pdf', p_subpool3_spa_spgl4_trans_10_20, 
       scale = 1.3, height = 1.4, width = 6.5, units = 'in')
  

p_subpool3_spa_spgl4_int_5_10 <- s3_int_trans_moveavg3 %>%
  filter(background == 's pGl4' & (spacing == 5 | spacing == 10) & MPRA == 'integrated') %>%
  ggplot(aes(x = dist, y = ave_ratio, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c('gray20', 'dodgerblue3'), name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10(limits = c(0.01, 0.4)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

p_subpool3_spa_spgl4_int_5_15 <- s3_int_trans_moveavg3 %>%
  filter(background == 's pGl4' & (spacing == 5 | spacing == 15) & MPRA == 'integrated') %>%
  ggplot(aes(x = dist, y = ave_ratio, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c('gray20', 'firebrick3'), 
                     name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10(limits = c(0.01, 0.4)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

p_subpool3_spa_spgl4_int_10_20 <- s3_int_trans_moveavg3 %>%
  filter(background == 's pGl4' & (spacing == 10 | spacing == 20) & MPRA == 'integrated') %>%
  ggplot(aes(x = dist, y = ave_ratio, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c('dodgerblue3', '#55C667FF'),
                     name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10(limits = c(0.01, 0.4)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10),
                     limits = c(67, 191)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

ggsave('plots/p_subpool3_spa_spgl4_int_5_10.pdf', p_subpool3_spa_spgl4_int_5_10, 
       scale = 1.3, height = 1.4, width = 6.5, units = 'in')

ggsave('plots/p_subpool3_spa_spgl4_int_5_15.pdf', p_subpool3_spa_spgl4_int_5_15, 
       scale = 1.3, height = 1.4, width = 6.5, units = 'in')

ggsave('plots/p_subpool3_spa_spgl4_int_10_20.pdf', p_subpool3_spa_spgl4_int_10_20, 
       scale = 1.3, height = 1.4, width = 6.5, units = 'in')


p_subpool3_spa_vchr5_trans_5_10 <- s3_int_trans_moveavg3 %>%
  filter(background == 'v chr5' & (spacing == 5 | spacing == 10) & MPRA == 'episomal') %>%
  ggplot(aes(x = dist, y = ave_ratio, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c('gray20', 'dodgerblue3'), name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10(limits = c(0.2, 4)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

p_subpool3_spa_vchr5_trans_5_15 <- s3_int_trans_moveavg3 %>%
  filter(background == 'v chr5' & (spacing == 5 | spacing == 15) & MPRA == 'episomal') %>%
  ggplot(aes(x = dist, y = ave_ratio, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c('gray20', 'firebrick3'), 
                     name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10(limits = c(0.2, 4)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

p_subpool3_spa_vchr5_trans_10_20 <- s3_int_trans_moveavg3 %>%
  filter(background == 'v chr5' & (spacing == 10 | spacing == 20) & MPRA == 'episomal') %>%
  ggplot(aes(x = dist, y = ave_ratio, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c('dodgerblue3', '#55C667FF'),
                     name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10(limits = c(0.2, 4)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10),
                     limits = c(67, 191)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

ggsave('plots/p_subpool3_spa_vchr5_trans_5_10.pdf', p_subpool3_spa_vchr5_trans_5_10, 
       scale = 1.3, height = 1.4, width = 6.5, units = 'in')

ggsave('plots/p_subpool3_spa_vchr5_trans_5_15.pdf', p_subpool3_spa_vchr5_trans_5_15, 
       scale = 1.3, height = 1.4, width = 6.5, units = 'in')

ggsave('plots/p_subpool3_spa_vchr5_trans_10_20.pdf', p_subpool3_spa_vchr5_trans_10_20, 
       scale = 1.3, height = 1.4, width = 6.5, units = 'in')


p_subpool3_spa_vchr5_int_5_10 <- s3_int_trans_moveavg3 %>%
  filter(background == 'v chr5' & (spacing == 5 | spacing == 10) & MPRA == 'integrated') %>%
  ggplot(aes(x = dist, y = ave_ratio, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c('gray20', 'dodgerblue3'), name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10(limits = c(0.03, 4)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

p_subpool3_spa_vchr5_int_5_15 <- s3_int_trans_moveavg3 %>%
  filter(background == 'v chr5' & (spacing == 5 | spacing == 15) & MPRA == 'integrated') %>%
  ggplot(aes(x = dist, y = ave_ratio, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c('gray20', 'firebrick3'), 
                     name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10(limits = c(0.03, 4)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

p_subpool3_spa_vchr5_int_10_20 <- s3_int_trans_moveavg3 %>%
  filter(background == 'v chr5' & (spacing == 10 | spacing == 20) & MPRA == 'integrated') %>%
  ggplot(aes(x = dist, y = ave_ratio, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c('dodgerblue3', '#55C667FF'),
                     name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10(limits = c(0.03, 4)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10),
                     limits = c(67, 191)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

ggsave('plots/p_subpool3_spa_vchr5_int_5_10.pdf', p_subpool3_spa_vchr5_int_5_10, 
       scale = 1.3, height = 1.4, width = 6.5, units = 'in')

ggsave('plots/p_subpool3_spa_vchr5_int_5_15.pdf', p_subpool3_spa_vchr5_int_5_15, 
       scale = 1.3, height = 1.4, width = 6.5, units = 'in')

ggsave('plots/p_subpool3_spa_vchr5_int_10_20.pdf', p_subpool3_spa_vchr5_int_10_20, 
       scale = 1.3, height = 1.4, width = 6.5, units = 'in')


#Comparison of subpool5 between integrated and episomal-------------------------

p_int_trans_ave_med_rep_sp5 <- int_trans %>%
  filter(subpool == 'subpool5') %>%
  ggplot(aes(ave_med_ratio, ave_ratio_22)) +
  geom_point(alpha = 0.2, size = 1) + 
  annotation_logticks() +
  xlab("Average integrated expression (a.u.)") +
  ylab("Average episomal expression (a.u.)") +
  scale_x_log10(limits = c(0.01, 20)) + 
  scale_y_log10(limits = c(0.09, 20))

ggsave('plots/p_int_trans_ave_med_rep_sp5.pdf', p_int_trans_ave_med_rep_sp5, 
       scale = 1.3, width = 3, height = 2.7, units = 'in')


#Compare subpool5 features between integrated and transient

s5_int_trans <- MPRA_ave %>%
  filter(subpool == 'subpool5') %>%
  subpool5()

pred_resid <- function(df1, x) {
  df2 <- df1 %>%
    add_predictions(x)
  df3 <- df2 %>%
    add_residuals(x)
  return(df3)
  print('processed pre_res_trans_int(df1, df2) in order of (data, model)')
}


#Consensus and weak vs. expression

p_s5_num_cons_num_weak_epi_int_spgl4 <- s5_int_trans %>%
  filter(background == 's pGl4') %>%
  ggplot(aes(as.factor(consensus), ave_ratio)) +
  facet_grid(MPRA ~ ., scale = 'free') +
  geom_boxplot(aes(fill = as.factor(weak)), outlier.size = 1, size = 0.3, 
               outlier.shape = 21, outlier.alpha = 1, 
               position = position_dodge(0.75), show.legend = FALSE) +
  scale_y_log10() + 
  panel_border(colour = 'black') +
  annotation_logticks(sides = 'l') +
  scale_fill_manual(name = 'number of\nweak sites', 
                    values = cbPalette7_grad_light) +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  ylab('Average expression (a.u.)') +
  xlab('Consensus CREs')

ggsave('plots/p_s5_num_cons_num_weak_epi_int_spgl4.pdf', 
       p_s5_num_cons_num_weak_epi_int_spgl4,
       scale = 1.3, width = 5, height = 3.25, units = 'in')

p_s5_num_cons_num_weak_epi_int_vchr9 <- s5_int_trans %>%
  filter(background == 'v chr9') %>%
  ggplot(aes(as.factor(consensus), ave_ratio)) +
  facet_grid(MPRA ~ ., scale = 'free') +
  geom_boxplot(aes(fill = as.factor(weak)), outlier.size = 1, size = 0.3, 
               outlier.shape = 21, outlier.alpha = 1, 
               position = position_dodge(0.75), show.legend = FALSE) +
  scale_y_log10() + 
  panel_border(colour = 'black') +
  annotation_logticks(sides = 'l') +
  scale_fill_manual(name = 'number of\nweak sites', 
                    values = cbPalette7_grad_light) +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  ylab('Average expression (a.u.)') +
  xlab('Consensus CREs')

ggsave('plots/p_s5_num_cons_num_weak_epi_int_vchr9.pdf', 
       p_s5_num_cons_num_weak_epi_int_vchr9,
       scale = 1.3, width = 5, height = 3.25, units = 'in')

p_s5_num_cons_num_weak_epi_int_vchr5 <- s5_int_trans %>%
  filter(background == 'v chr5') %>%
  ggplot(aes(as.factor(consensus), ave_ratio)) +
  facet_grid(MPRA ~ ., scale = 'free') +
  geom_boxplot(aes(fill = as.factor(weak)), outlier.size = 1, size = 0.3, 
               outlier.shape = 21, outlier.alpha = 1, 
               position = position_dodge(0.75), show.legend = FALSE) +
  scale_y_log10() + 
  panel_border(colour = 'black') +
  annotation_logticks(sides = 'l') +
  scale_fill_manual(name = 'number of\nweak sites', 
                    values = cbPalette7_grad_light) +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  ylab('Average expression (a.u.)') +
  xlab('Consensus CREs')

ggsave('plots/p_s5_num_cons_num_weak_epi_int_vchr5.pdf', 
       p_s5_num_cons_num_weak_epi_int_vchr5,
       scale = 1.3, width = 5, height = 3.25, units = 'in')


#Linear models

#Including weak sites, all independent, independent background 

ind_site_ind_back <- function(df) {
  model <- lm(ave_ratio ~ background + site1 + site2 + site3 + site4 + site5 + site6, 
              data = df)
}

ind_site_ind_back_sumtidy <- function(df) {
  df <- tidy(df)
  sites <- df %>%
    filter(str_detect(term, '^site')) %>%
    mutate(term = gsub('consensus', '_consensus', term)) %>%
    mutate(term = gsub('weak', '_weak', term)) %>%
    separate(term, into = c('variable', 'type'), sep = "_")
  background <- df %>%
    mutate(term = gsub('\\(Intercept\\)', 'backgroundv chr9', term)) %>%
    filter(str_detect(term, '^background')) %>%
    mutate(term = gsub('background', 'background_', term)) %>%
    separate(term, into = c('variable', 'type'), sep = '_')
  sum <- rbind(sites, background)
  return(sum)
}

subpool5_ncw <- s5_int_trans %>%
  mutate(site1 = gsub('nosite', 'anosite', site1)) %>%
  mutate(site2 = gsub('nosite', 'anosite', site2)) %>%
  mutate(site3 = gsub('nosite', 'anosite', site3)) %>%
  mutate(site4 = gsub('nosite', 'anosite', site4)) %>%
  mutate(site5 = gsub('nosite', 'anosite', site5)) %>%
  mutate(site6 = gsub('nosite', 'anosite', site6)) %>%
  mutate(background = gsub('v chr9', 'av chr9', background)) %>%
  var_log10()

ind_site_ind_back_epi <- subpool5_ncw %>%
  filter(MPRA == 'episomal') %>%
  ind_site_ind_back()

ind_site_ind_back_sum_epi <- ind_site_ind_back_sumtidy(ind_site_ind_back_epi)

ind_site_ind_back_anova_epi <- tidy(anova(ind_site_ind_back_epi)) %>%
  mutate(term_fctr = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq)

ind_site_ind_back_p_r_epi <- pred_resid(filter(subpool5_ncw, MPRA == 'episomal'), 
                                    ind_site_ind_back_epi)

lessthan1_2color <- c('red', 'black', 'black', 'black', 'black', 'black', 'black')

p_ind_site_ind_back_epi <- ggplot(ind_site_ind_back_p_r_epi, 
                              aes(ave_ratio, pred, 
                                  color = as.factor(consensus))) +
  geom_point(alpha = 0.2, size = 1, show.legend = FALSE) +
  scale_color_manual(values = lessthan1_2color) +
  scale_x_continuous(name = 'Measured expression', breaks = c(-1:1),
                     limits = c(-1.5, 1.8)) + 
  scale_y_continuous(name = 'Predicted expression', breaks = c(-1:1),
                     limits = c(-1.5, 1.8)) +
  annotation_logticks(sides = 'bl') +
  annotate("text", x = -0.5, y = 1, 
           label = paste('r =', 
                         round(cor(ind_site_ind_back_p_r_epi$pred,
                                   ind_site_ind_back_p_r_epi$ave_ratio,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))

p_ind_site_ind_back_sum_epi <- ind_site_ind_back_sum_epi %>%
  mutate(type = factor(type, 
                       levels = c('v chr9', 's pGl4', 'v chr5', 'consensus', 
                                  'weak'))) %>%
  ggplot(aes(variable, estimate, fill = type)) + 
  geom_bar(stat = 'identity', position = 'dodge', color = 'gray60', 
           size = 0.3) + 
  geom_hline(yintercept = 0, size = 0.25) +
  scale_x_discrete(position = 'bottom') + 
  scale_fill_viridis(discrete = TRUE) + 
  theme(axis.ticks.x = element_blank(), legend.position = 'top',
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank()) +
  ylab('Weight')

p_ind_site_ind_back_anova_epi <- ind_site_ind_back_anova_epi %>%
  ggplot(aes(term_fctr, per_sumsq)) + 
  geom_bar(stat = 'identity') + 
  ylab('Proportion of\nvariance explained') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank())

ggsave('plots/p_ind_site_ind_back_epi.pdf', p_ind_site_ind_back_epi, 
       scale = 1.3, width = 2.3, height = 2.3, units = 'in')

ggsave('plots/p_ind_site_ind_back_sum_epi.pdf', p_ind_site_ind_back_sum_epi,
       scale = 1.3, width = 3, height = 2.5, units = 'in')

ggsave('plots/p_ind_site_ind_back_anova_epi.pdf', p_ind_site_ind_back_anova_epi,
       scale = 1.3, width = 2.5, height = 2.5)



ind_site_ind_back_int <- subpool5_ncw %>%
  filter(MPRA == 'integrated') %>%
  ind_site_ind_back()

ind_site_ind_back_sum_int <- ind_site_ind_back_sumtidy(ind_site_ind_back_int)

ind_site_ind_back_anova_int <- tidy(anova(ind_site_ind_back_int)) %>%
  mutate(term_fctr = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq)

ind_site_ind_back_p_r_int <- pred_resid(filter(subpool5_ncw, MPRA == 'integrated'), 
                                        ind_site_ind_back_int)

p_ind_site_ind_back_int <- ggplot(ind_site_ind_back_p_r_int, 
                                  aes(ave_ratio, pred, 
                                      color = as.factor(consensus))) +
  geom_point(alpha = 0.2, size = 1, show.legend = FALSE) +
  scale_color_manual(values = lessthan1_2color) +
  scale_x_continuous(name = 'Measured expression', breaks = c(-2:1),
                     limits = c(-2.1, 1.8)) + 
  scale_y_continuous(name = 'Predicted expression', breaks = c(-2:1),
                     limits = c(-2.1, 1.8)) +
  annotation_logticks(sides = 'bl') +
  annotate("text", x = -0.5, y = 1, 
           label = paste('r =', 
                         round(cor(ind_site_ind_back_p_r_int$pred,
                                   ind_site_ind_back_p_r_int$ave_ratio,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))

p_ind_site_ind_back_sum_int <- ind_site_ind_back_sum_int %>%
  mutate(type = factor(type, 
                       levels = c('v chr9', 's pGl4', 'v chr5', 'consensus', 
                                  'weak'))) %>%
  ggplot(aes(variable, estimate, fill = type)) + 
  geom_bar(stat = 'identity', position = 'dodge', color = 'gray60', 
           size = 0.3) + 
  geom_hline(yintercept = 0, size = 0.25) +
  scale_x_discrete(position = 'bottom') + 
  scale_fill_viridis(discrete = TRUE) + 
  theme(axis.ticks.x = element_blank(), legend.position = 'top',
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank()) +
  scale_y_continuous(breaks = c(-2, -1, 0, 0.5)) +
  ylab('Weight')

p_ind_site_ind_back_anova_int <- ind_site_ind_back_anova_int %>%
  ggplot(aes(term_fctr, per_sumsq)) + 
  geom_bar(stat = 'identity') + 
  ylab('Proportion of\nvariance explained') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank())

ggsave('plots/p_ind_site_ind_back_int.pdf', p_ind_site_ind_back_int, 
       scale = 1.3, width = 2.3, height = 2.3, units = 'in')

ggsave('plots/p_ind_site_ind_back_sum_int.pdf', p_ind_site_ind_back_sum_int,
       scale = 1.3, width = 3, height = 2.5, units = 'in')

ggsave('plots/p_ind_site_ind_back_anova_int.pdf', p_ind_site_ind_back_anova_int,
       scale = 1.3, width = 2.5, height = 2.5)


#break down integrated vs. transient by site combos

cons_int_epi_lm <- lm(ave_ratio_22 ~ ave_med_ratio, 
                      data = var_log10(filter(s5_int_trans, 
                                              site_combo == 'consensus')))

cons_pearsons_log10 <- s5_int_trans %>%
  var_log10() %>%
  filter(site_combo == 'consensus')

cons_pearsons <- round(cor(cons_pearsons_log10$ave_med_ratio, 
                           cons_pearsons_log10$ave_ratio_22, 
                           use = "pairwise.complete.obs", 
                           method = "pearson"), 2)

s5_int_trans_cons_lm <- pred_resid(var_log10(s5_int_trans), cons_int_epi_lm)

p_s5_int_trans_site_combo <- s5_int_trans_cons_lm %>%
  filter(site_combo != 'none') %>%
  mutate(site_combo = factor(site_combo, 
                             levels = c('consensus', 'weak', 'mixed'))) %>%
  ggplot(aes(ave_med_ratio, ave_ratio_22)) +
  facet_grid(. ~ site_combo) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_line(aes(ave_med_ratio, pred), color = 'red', size = 0.5) +
  annotation_logticks() +
  scale_y_continuous(breaks = seq(from = -1, to = 1, by =1)) +
  xlab('Average integrated expression (a.u.)') +
  ylab('Average episomal\nexpression (a.u.)') +
  panel_border(colour = 'black') +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

ggsave('plots/p_s5_int_trans_site_combo.pdf', p_s5_int_trans_site_combo,
       scale = 1.3, width = 4, height = 2.25, units = 'in')

p_s5_int_trans_site_combo_resid <- s5_int_trans_cons_lm %>%
  ggplot(aes(as.factor(consensus), resid, fill = as.factor(weak))) +
  geom_boxplot(outlier.size = 1, size = 0.3, 
               outlier.shape = 21, outlier.alpha = 1, 
               position = position_dodge(0.75), show.legend = FALSE) +
  scale_fill_manual(name = 'number of\nweak sites', 
                    values = cbPalette7_grad_light) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5, color = 'red') +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

ggsave('plots/p_s5_int_trans_site_combo_resid.pdf', 
       p_s5_int_trans_site_combo_resid, scale = 1.3, width = 4, height = 2.25,
       units = 'in')


#Look at competition between sites as forskolin decreases

MPRA_ave_allconc <- int_rep_1_2 %>%
  select(-med_ratio_br1, -med_ratio_br2, -mad_br1, -mad_br2, -mad_over_med_br1,
         -mad_over_med_br2, -barcodes_DNA_br1, -barcodes_DNA_br2) %>%
  ungroup() %>%
  filter(subpool != 'control') %>%
  mutate(
    name = gsub('Smith R. Vista chr9:83712599-83712766', 'v chr9', name),
    name = gsub('Vista Chr5:88673410-88674494', 'v chr5', name),
    name = gsub('scramble pGL4.29 Promega 1-63 \\+ 1-87', 's pGl4', name)
  ) %>%
  mutate(background = name) %>%
  mutate(background = str_sub(background, 
                              nchar(background)-5,
                              nchar(background))) %>%
  inner_join(trans_back_norm_conc, 
             by = c('subpool', 'name', 'most_common', 'background')) %>%
  mutate(integrated = (barcodes_RNA_br1 + barcodes_RNA_br2)/2) %>%
  rename(episomal = ave_barcode) %>%
  select(-barcodes_RNA_br1, -barcodes_RNA_br2) %>%
  gather(integrated, episomal, key = 'MPRA', value = 'barcodes') %>%
  mutate(integrated = ave_med_ratio) %>%
  mutate(episomal = ave_ratio_norm) %>%
  rename(ave_sum_ratio = ave_ratio_norm) %>%
  gather(integrated, episomal, key = 'MPRA2', value = 'ave_ratio') %>%
  filter((MPRA == 'integrated' & MPRA2 == 'integrated') | (MPRA == 'episomal' & MPRA2 == 'episomal')) %>%
  select(-MPRA2)

subpool5(MPRA_ave_allconc) %>%
  filter(site_combo == 'consensus') %>%
  ggplot(aes(ave_med_ratio, ave_sum_ratio, color = consensus)) +
  scale_color_viridis() +
  geom_point(alpha = 0.2) +
  facet_wrap(~ conc) +
  scale_x_log10(name = 'integrated expression') +
  scale_y_log10(name = 'transient expression') + 
  panel_border(colour = 'black') +
  annotation_logticks(sides = 'bl')

cons_int_epi_lm_conc <- function(df) {
  df <- df %>%
    filter(subpool == 'subpool5') %>%
    subpool5() %>%
    filter(site_combo == 'consensus')
  fit <- lm(ave_sum_ratio ~ ave_med_ratio, data = df)
  return(fit)
}


#fit a linear model to other concentrations in episomal

subpool5_ncw_allconc <- subpool5(MPRA_ave_allconc) %>%
  mutate(site1 = gsub('nosite', 'anosite', site1)) %>%
  mutate(site2 = gsub('nosite', 'anosite', site2)) %>%
  mutate(site3 = gsub('nosite', 'anosite', site3)) %>%
  mutate(site4 = gsub('nosite', 'anosite', site4)) %>%
  mutate(site5 = gsub('nosite', 'anosite', site5)) %>%
  mutate(site6 = gsub('nosite', 'anosite', site6)) %>%
  mutate(background = gsub('v chr9', 'av chr9', background)) %>%
  mutate(ave_ratio = log10(ave_ratio)) %>%
  filter(MPRA == 'episomal' & conc == 2^-4)
  
ind_site_ind_back_epi_2_1 <- subpool5_ncw_allconc %>%
  ind_site_ind_back()

ind_site_ind_back_sum_epi_2_1 <- ind_site_ind_back_sumtidy(ind_site_ind_back_epi_2_1)

ind_site_ind_back_anova_epi_2_1 <- tidy(anova(ind_site_ind_back_epi_2_1)) %>%
  mutate(term_fctr = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq)

ind_site_ind_back_p_r_epi_2_1 <- pred_resid(subpool5_ncw_allconc, 
                                            ind_site_ind_back_epi_2_1)

ggplot(ind_site_ind_back_p_r_epi_2_1, 
       aes(ave_ratio, pred, color = as.factor(consensus))) +
  geom_point(alpha = 0.2, size = 1, show.legend = FALSE) +
  scale_color_manual(values = lessthan1_2color) +
  scale_x_continuous(name = 'Measured expression') + 
  scale_y_continuous(name = 'Predicted expression') +
  annotation_logticks(sides = 'bl') +
  annotate("text", x = -0.5, y = 1, 
           label = paste('r =', 
                         round(cor(ind_site_ind_back_p_r_epi_2_1$pred,
                                   ind_site_ind_back_p_r_epi_2_1$ave_ratio,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))

ind_site_ind_back_sum_epi_2_1 %>%
  mutate(type = factor(type, 
                       levels = c('v chr9', 's pGl4', 'v chr5', 'consensus', 
                                  'weak'))) %>%
  ggplot(aes(variable, estimate, fill = type)) + 
  geom_bar(stat = 'identity', position = 'dodge', color = 'gray60', 
           size = 0.3) + 
  geom_hline(yintercept = 0, size = 0.25) +
  scale_x_discrete(position = 'bottom') + 
  scale_fill_viridis(discrete = TRUE) + 
  theme(axis.ticks.x = element_blank(), legend.position = 'top',
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank()) +
  ylab('Weight')

ind_site_ind_back_anova_epi_2_1 %>%
  ggplot(aes(term_fctr, per_sumsq)) + 
  geom_bar(stat = 'identity') + 
  ylab('Proportion of\nvariance explained') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank())




#Comparison to integrated with median episomal----------------------------------

#First look at periodicity offset in episomal for main figure 

s3_tidy <- med_rep_0_22_A_B %>%
  ungroup() %>%
  mutate(
    name = gsub('Smith R. Vista chr9:83712599-83712766', 'v chr9', name),
    name = gsub('Vista Chr5:88673410-88674494', 'v chr5', name),
    name = gsub('scramble pGL4.29 Promega 1-63 \\+ 1-87', 's pGl4', name)
  ) %>%
  mutate(background = name) %>%
  mutate(background = str_sub(background, 
                              nchar(background)-5,
                              nchar(background))) %>%
  subpool3()


s3_tidy_moveavg3 <- s3_tidy %>%
  mutate(ave_ratio_22 = (med_ratio_22A + med_ratio_22B)/2) %>%
  select(background, spacing, dist, ave_ratio_22) %>%
  group_by(background, spacing) %>%
  arrange(dist, .by_group = TRUE) %>%
  nest() %>%
  mutate(ave_3 = map(.$data, moveavg_dist3)) %>%
  unnest() %>%
  select(-dist1, -ave_ratio_221)

p_subpool3_spa_4_vchr9_5_10 <- s3_tidy_moveavg3 %>%
  filter(background == 'v chr9' & dist < 124 & (spacing == 5 | spacing == 10)) %>%
  ggplot(aes(x = dist, y = ave_ratio_22, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c('gray20', 'dodgerblue3'), name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  geom_vline(xintercept = c(78, 88.5), color = 'gray20', linetype = 2, 
             alpha = 0.5) +
  geom_vline(xintercept = c(82.5, 92.5), color = 'dodgerblue3', linetype = 2, 
             alpha = 0.5) +
  scale_y_log10(limits = c(0.1, 2)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 124, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

p_subpool3_spa_4_vchr9_5_15 <- s3_tidy_moveavg3 %>%
  filter(background == 'v chr9' & dist < 124 & (spacing == 5 | spacing == 15)) %>%
  ggplot(aes(x = dist, y = ave_ratio_22, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c('gray20', 'firebrick3'), 
                     name = 'spacing (bp)') +
  scale_fill_manual(values = c('gray20', 'firebrick3'), 
                    name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  geom_vline(xintercept = c(78, 88.5), color = 'gray20', linetype = 2, 
             alpha = 0.5) +
  scale_y_log10(limits = c(0.1, 2)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 124, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

p_subpool3_spa_4_vchr9_10_20 <- s3_tidy_moveavg3 %>%
  filter(background == 'v chr9' & dist < 124 & (spacing == 10 | spacing == 20)) %>%
  ggplot(aes(x = dist, y = ave_ratio_22, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c('dodgerblue3', '#55C667FF'),
                     name = 'spacing (bp)') +
  scale_fill_manual(values = c('dodgerblue3', '#55C667FF'),
                    name = 'spacing (bp)') +
  ylab('Average expression (a.u.)') + 
  panel_border(colour = 'black') +
  geom_vline(xintercept = c(82.5, 92.5), color = 'dodgerblue3', linetype = 2, 
             alpha = 0.5) +
  scale_y_log10(limits = c(0.1, 2)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 124, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

ggsave('plots/p_subpool3_med_spa_4_vchr9_5_10.pdf', p_subpool3_spa_4_vchr9_5_10, 
       scale = 1.3, height = 1.65, width = 4.8, units = 'in')

ggsave('plots/p_subpool3_med_spa_4_vchr9_5_15.pdf', p_subpool3_spa_4_vchr9_5_15, 
       scale = 1.3, height = 1.65, width = 4.8, units = 'in')

ggsave('plots/p_subpool3_med_spa_4_vchr9_10_20.pdf', p_subpool3_spa_4_vchr9_10_20, 
       scale = 1.3, height = 1.65, width = 4.8, units = 'in')


#Compare features to integrated

int_rep_1_2 <- read_tsv('../20171129_intLib/rep_1_2.txt') %>%
  mutate(ave_med_ratio = (med_ratio_br1 + med_ratio_br2)/2)

int_trans <- med_rep_0_22_A_B %>%
  select(subpool, name, most_common, barcodes_DNA, med_ratio_22A, 
         barcodes_RNA_22A, med_ratio_22B, barcodes_RNA_22B) %>%
  mutate(ave_ratio_22 = (med_ratio_22A + med_ratio_22B)/2) %>%
  inner_join(int_rep_1_2, by = c('subpool', 'name', 'most_common'))

MPRA_ave <- int_trans %>%
  ungroup() %>%
  filter(subpool != 'control') %>%
  mutate(
    name = gsub('Smith R. Vista chr9:83712599-83712766', 'v chr9', name),
    name = gsub('Vista Chr5:88673410-88674494', 'v chr5', name),
    name = gsub('scramble pGL4.29 Promega 1-63 \\+ 1-87', 's pGl4', name)
  ) %>%
  mutate(background = name) %>%
  mutate(background = str_sub(background, 
                              nchar(background)-5,
                              nchar(background))) %>%
  select(subpool, name, background, barcodes_RNA_br1, barcodes_RNA_br2, 
         med_ratio_br1, med_ratio_br2, ave_med_ratio, barcodes_RNA_22A, 
         barcodes_RNA_22B, med_ratio_22A, med_ratio_22B, ave_ratio_22) %>%
  mutate(integrated = (barcodes_RNA_br1 + barcodes_RNA_br2)/2) %>%
  mutate(episomal = (barcodes_RNA_22A + barcodes_RNA_22B)/2) %>%
  select(-barcodes_RNA_br1, -barcodes_RNA_br2, -barcodes_RNA_22A, 
         -barcodes_RNA_22B) %>%
  gather(integrated, episomal, key = 'MPRA', value = 'barcodes') %>%
  mutate(integrated = ave_med_ratio) %>%
  mutate(episomal = ave_ratio_22) %>%
  gather(integrated, episomal, key = 'MPRA2', value = 'ave_ratio') %>%
  filter((MPRA == 'integrated' & MPRA2 == 'integrated') | (MPRA == 'episomal' & MPRA2 == 'episomal')) %>%
  select(-MPRA2)

int_trans_log10 <- var_log10(int_trans)

int_trans_pearsons <- tibble(
  sample = c('all', 'subpool3', 'subpool5'),
  pearsons = c(round(cor(int_trans_log10$ave_med_ratio, 
                         int_trans_log10$ave_ratio_22, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(filter(int_trans_log10, 
                                subpool == 'subpool3')$ave_med_ratio,
                         filter(int_trans_log10, 
                                subpool == 'subpool3')$ave_ratio_22,
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(filter(int_trans_log10, 
                                subpool == 'subpool5')$ave_med_ratio,
                         filter(int_trans_log10, 
                                subpool == 'subpool5')$ave_ratio_22,
                         use = "pairwise.complete.obs", method = "pearson"), 
                     3)))

write_csv(int_trans_pearsons, 'int_trans_med_pearsons.csv')

p_int_trans_ave_med_rep_sp3 <- int_trans %>%
  filter(subpool == 'subpool3') %>%
  ggplot(aes(ave_med_ratio, ave_ratio_22)) +
  geom_point(alpha = 0.2, size = 1) +
  xlab("Average integrated expression (a.u.)") +
  ylab("Average episomal expression (a.u.)") +
  scale_x_log10() + 
  scale_y_log10() +
  annotation_logticks(scaled = TRUE)

ggsave('plots/p_int_med_trans_med_ave_rep_sp3.pdf', p_int_trans_ave_med_rep_sp3, 
       scale = 1.3, width = 3, height = 2.7, units = 'in')


#Compare subpool3 features between integrated and transient

s3_int_trans <- MPRA_ave %>%
  filter(subpool == 'subpool3') %>%
  subpool3()

s3_int_trans_bin20bp <- s3_int_trans %>%
  mutate(bin = cut(dist, seq(from = 64, to = 196, by = 22),
                   labels = c('64-86', '86-108', '108-130', '130-152',
                              '152-174', '174-196')))

p_s3_dist_trans_bin20bp <- s3_int_trans_bin20bp %>%
  mutate(background = factor(background, levels = c('s pGl4', 'v chr5'))) %>%
  filter(spacing != 0 & spacing != 70 & background != 'v chr9' & bin != '174-196' & MPRA == 'episomal') %>%
  ggplot(aes(bin, ave_ratio)) +
  geom_signif(comparisons = list(c('64-86', '86-108')), y_position = log10(5)) +
  geom_signif(comparisons = list(c('64-86', '108-130')), y_position = log10(15)) +
  facet_grid(. ~ background) +
  geom_jitter(aes(color = as.factor(spacing)), 
              position=position_jitter(width=0.3, height=0), alpha = 0.75,
              size = 0.75) +
  geom_boxplot(outlier.shape=NA, size = 0.3, position = position_dodge(1),
               show.legend = FALSE, alpha = 0) +
  scale_color_manual(values = spacing_5_20_palette, name = 'spacing (bp)') +
  theme(legend.position = 'right', axis.ticks.x = element_blank(), 
        strip.background = element_rect(colour="black", fill="white"),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_log10(limits = c(0.065, 30)) +
  annotation_logticks(sides = 'l') +
  panel_border(colour = 'black') +
  ylab('Average expression (a.u.)') +
  xlab('Distance to minimal promoter (bp)')

ggsave('plots/p_s3_dist_trans_med_bin20bp.pdf', p_s3_dist_trans_bin20bp,
       scale = 1.3, width = 5, height = 2.5, unit = 'in')

p_s3_space_trans <- s3_int_trans %>%
  filter(dist <= 108 & background != 'v chr9' & MPRA == 'episomal') %>%
  mutate(background = factor(background, levels = c('s pGl4', 'v chr5'))) %>%
  ggplot(aes(as.factor(spacing), ave_ratio)) +
  facet_grid(. ~ background) +
  geom_signif(comparisons = list(c('5', '10')), y_position = log10(3)) +
  geom_signif(comparisons = list(c('5', '15')), y_position = log10(7)) +
  geom_signif(comparisons = list(c('5', '20')), y_position = log10(15)) +
  geom_boxplot(size = 0.3, position = position_dodge(1)) +
  theme(axis.ticks.x = element_blank(), 
        strip.background = element_rect(colour="black", fill="white")) + 
  scale_y_log10(limits = c(0.060, 20)) +
  annotation_logticks(sides = 'l') +
  panel_border(colour = 'black') +
  ylab('Average expression (a.u.)') +
  xlab('Spacing (bp)')

ggsave('plots/p_s3_space_trans_med.pdf', p_s3_space_trans, scale = 1.3,
       width = 4, height = 2, units = 'in')


#Compare subpool 5 between integrated and transient.....fix this

p_int_trans_ave_med_rep_sp5 <- int_trans %>%
  filter(subpool == 'subpool5') %>%
  ggplot(aes(ave_med_ratio, ave_ratio_22)) +
  geom_point(alpha = 0.2, size = 1) + 
  annotation_logticks() +
  xlab("Average integrated expression (a.u.)") +
  ylab("Average episomal expression (a.u.)") +
  scale_x_log10() + 
  scale_y_log10()

ggsave('plots/p_int_med_trans_med_ave_rep_sp5.pdf', p_int_trans_ave_med_rep_sp5, 
       scale = 1.3, width = 3, height = 2.7, units = 'in')


#Compare subpool5 features between integrated and transient

s5_int_trans <- MPRA_ave %>%
  filter(subpool == 'subpool5') %>%
  subpool5()

subpool5_ncw <- s5_int_trans %>%
  mutate(site1 = gsub('nosite', 'anosite', site1)) %>%
  mutate(site2 = gsub('nosite', 'anosite', site2)) %>%
  mutate(site3 = gsub('nosite', 'anosite', site3)) %>%
  mutate(site4 = gsub('nosite', 'anosite', site4)) %>%
  mutate(site5 = gsub('nosite', 'anosite', site5)) %>%
  mutate(site6 = gsub('nosite', 'anosite', site6)) %>%
  mutate(background = gsub('v chr9', 'av chr9', background)) %>%
  var_log10()

ind_site_ind_back_epi <- subpool5_ncw %>%
  filter(MPRA == 'episomal') %>%
  ind_site_ind_back()

ind_site_ind_back_sum_epi <- ind_site_ind_back_sumtidy(ind_site_ind_back_epi)

ind_site_ind_back_anova_epi <- tidy(anova(ind_site_ind_back_epi)) %>%
  mutate(term_fctr = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq)

ind_site_ind_back_p_r_epi <- pred_resid(filter(subpool5_ncw, MPRA == 'episomal'), 
                                        ind_site_ind_back_epi)

p_ind_site_ind_back_epi <- ggplot(ind_site_ind_back_p_r_epi, 
                                  aes(ave_ratio, pred)) +
  geom_point(alpha = 0.2, size = 1) +
  scale_x_continuous(name = 'Measured expression', breaks = c(-1:1),
                     limits = c(-1.6, 2)) + 
  scale_y_continuous(name = 'Predicted expression', breaks = c(-1:1),
                     limits = c(-1.6, 2)) +
  annotation_logticks(sides = 'bl') +
  annotate("text", x = -0.5, y = 1, 
           label = paste('r =', 
                         round(cor(ind_site_ind_back_p_r_epi$pred,
                                   ind_site_ind_back_p_r_epi$ave_ratio,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))

p_ind_site_ind_back_sum_epi <- ind_site_ind_back_sum_epi %>%
  mutate(type = factor(type, 
                       levels = c('v chr9', 's pGl4', 'v chr5', 'consensus', 
                                  'weak'))) %>%
  ggplot(aes(variable, estimate, fill = type)) + 
  geom_bar(stat = 'identity', position = 'dodge', color = 'gray60', 
           size = 0.3) + 
  geom_hline(yintercept = 0, size = 0.25) +
  scale_x_discrete(position = 'bottom') + 
  scale_fill_viridis(discrete = TRUE) + 
  theme(axis.ticks.x = element_blank(), legend.position = 'top',
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank()) +
  ylab('Weight')

p_ind_site_ind_back_anova_epi <- ind_site_ind_back_anova_epi %>%
  ggplot(aes(term_fctr, per_sumsq)) + 
  geom_bar(stat = 'identity') + 
  ylab('Proportion of\nvariance explained') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank())

ggsave('plots/p_ind_site_ind_back_epi_med.pdf', p_ind_site_ind_back_epi, 
       scale = 1.3, width = 2.3, height = 2.3, units = 'in')

ggsave('plots/p_ind_site_ind_back_sum_epi_med.pdf', p_ind_site_ind_back_sum_epi,
       scale = 1.3, width = 3, height = 2.5, units = 'in')

ggsave('plots/p_ind_site_ind_back_anova_epi_med.pdf', p_ind_site_ind_back_anova_epi,
       scale = 1.3, width = 2.5, height = 2.5)


cons_int_epi_lm <- lm(ave_ratio_22 ~ ave_med_ratio, 
                      data = var_log10(filter(s5_int_trans, 
                                              site_combo == 'consensus')))

s5_int_trans_cons_lm <- pred_resid(var_log10(s5_int_trans), cons_int_epi_lm)

p_s5_int_trans_site_combo <- s5_int_trans_cons_lm %>%
  filter(site_combo != 'none') %>%
  mutate(site_combo = factor(site_combo, 
                             levels = c('consensus', 'weak', 'mixed'))) %>%
  ggplot(aes(ave_med_ratio, ave_ratio_22)) +
  facet_grid(. ~ site_combo) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_line(aes(ave_med_ratio, pred), color = 'red', size = 0.5) +
  annotation_logticks() +
  scale_y_continuous(breaks = seq(from = -1, to = 1, by =1)) +
  xlab('Average integrated expression (a.u.)') +
  ylab('Average episomal\nexpression (a.u.)') +
  panel_border(colour = 'black') +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

ggsave('plots/p_s5_int_med_trans_med_site_combo.pdf', p_s5_int_trans_site_combo,
       scale = 1.3, width = 4, height = 2.25, units = 'in')

p_s5_int_trans_site_combo_resid <- s5_int_trans_cons_lm %>%
  ggplot(aes(as.factor(consensus), resid, fill = as.factor(weak))) +
  geom_boxplot(outlier.size = 1, size = 0.3, 
               outlier.shape = 21, outlier.alpha = 1, 
               position = position_dodge(0.75), show.legend = FALSE) +
  scale_fill_manual(name = 'number of\nweak sites', 
                    values = cbPalette7_grad_light) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5, color = 'red') +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

ggsave('plots/p_s5_int_med_trans_med_site_combo_resid.pdf', 
       p_s5_int_trans_site_combo_resid, scale = 1.3, width = 4, height = 2.25,
       units = 'in')


#Titration range between genomic and episomal using luciferase------------------

#This csv file was transformed to make R plotting more straightforward, 
#replicates between genomic and episomal assays are not the same

titration_luc <- read_csv("../../../Plate Reader/170726_trans_int_R.csv") %>%
  mutate(RLU_epi = luciferase_epi/renilla_epi) %>%
  mutate(forskolin = log2(forskolin))

p_gen_titration_luc <- titration_luc %>%
  ggplot(aes(forskolin, luciferase_gen)) +
  geom_point() +
  geom_smooth() +
  ylab('bulk luminescence (a.u.)') +
  scale_x_continuous(breaks = (-2:5), 'log2 forskolin (µM)') +
  annotation_logticks(sides = 'b')

p_epi_titration_luc <- titration_luc %>%
  ggplot(aes(forskolin, RLU_epi)) +
  geom_point() +
  geom_smooth() +
  ylab('bulk relative\nluminescence (a.u.)') +
  scale_x_continuous(breaks = (-2:5), 'log2 forskolin (µM)') +
  annotation_logticks(sides = 'b')

ggsave('plots/p_gen_titration_luc.pdf', p_gen_titration_luc, scale = 1.3,
       width = 4.5, height = 2.5, units = 'in')

ggsave('plots/p_epi_titration_luc.pdf', p_epi_titration_luc, scale = 1.3,
       width = 4.5, height = 2.5, units = 'in')


#Comparison to Luciferase assays------------------------------------------------

#Since using background-normalized reads to account for read stealing. Also, 
#since using background-normalized reads, will not compare negative controls as 
#these values are all 1 in sequencing data

#Want to add pc control pGL4.29 Promega 1-63 + 1-87 to background-normalized 
#reads

trans_back_norm_pc_spGl4 <- rep_0_22_A_B %>%
  filter(name == 'pGL4.29 Promega 1-63 + 1-87') %>%
  mutate(name = 'pGL4.29 Promega 1-63 + 1-87_scramble pGL4.29 Promega 1-63 + 1-87') %>%
  mutate(subpool = 'subpool3') %>%
  rbind(rep_0_22_A_B) %>%
  back_norm()


#Isolating variants used in both assays from the transient data, make df untidy
#with conc, Expave/Expave0, assay type and name
trans_untidy_conc_expto0 <- function(df) {
  sequences <- df %>%
    filter(name %in% c('subpool3_2BS 16 bp spacing consensus+flank x2_dist_79_s pGl4',
                       'subpool3_2BS 6 bp spacing consensus+flank x2_dist_9_s pGl4',
                       'subpool3_2BS 11 bp spacing consensus+flank x2_dist_14_s pGl4',
                       'subpool5_no_site_consensus_weak_consensus_no_site_weak_v chr9',
                       'pGL4.29 Promega 1-63 + 1-87_s pGl4')) %>%
    select(subpool, name, ave_ratio_0_norm, ave_ratio_2_5_norm, 
           ave_ratio_2_4_norm, ave_ratio_2_3_norm, ave_ratio_2_2_norm, 
           ave_ratio_2_1_norm, ave_ratio_20_norm, ave_ratio_22_norm) %>%
    mutate(is0 = ave_ratio_0_norm)
  df_0 <- sequences %>%
    select(subpool, name, is0, ave_ratio_0_norm) %>%
    mutate(conc = 0) %>%
    rename(ave_ratio_norm = ave_ratio_0_norm)
  df_2_5 <- sequences %>%
    select(subpool, name, is0, ave_ratio_2_5_norm) %>%
    mutate(conc = 2^-5) %>%
    rename(ave_ratio_norm = ave_ratio_2_5_norm)
  df_2_4 <- sequences %>%
    select(subpool, name, is0, ave_ratio_2_4_norm) %>%
    mutate(conc = 2^-4) %>%
    rename(ave_ratio_norm = ave_ratio_2_4_norm)
  df_2_3 <- sequences %>%
    select(subpool, name, is0, ave_ratio_2_3_norm) %>%
    mutate(conc = 2^-3) %>%
    rename(ave_ratio_norm = ave_ratio_2_3_norm)
  df_2_2 <- sequences %>%
    select(subpool, name, is0, ave_ratio_2_2_norm) %>%
    mutate(conc = 2^-2) %>%
    rename(ave_ratio_norm = ave_ratio_2_2_norm)
  df_2_1 <- sequences %>%
    select(subpool, name, is0, ave_ratio_2_1_norm) %>%
    mutate(conc = 2^-1) %>%
    rename(ave_ratio_norm = ave_ratio_2_1_norm)
  df_20 <- sequences %>%
    select(subpool, name, is0, ave_ratio_20_norm) %>%
    mutate(conc = 2^0) %>%
    rename(ave_ratio_norm = ave_ratio_20_norm)
  df_22 <- sequences %>%
    select(subpool, name, is0, ave_ratio_22_norm) %>%
    mutate(conc = 2^2) %>%
    rename(ave_ratio_norm = ave_ratio_22_norm)
  df_0_22 <- rbind(df_0, df_2_5, df_2_4, df_2_3, df_2_2, df_2_1, df_20, 
                   df_22) %>%
    mutate(ave_exp_to0 = ave_ratio_norm/is0) %>%
    select(-is0, -ave_ratio_norm) %>%
    mutate(assay = 'trans')
  return(df_0_22)
}
  
trans_comp <- trans_untidy_conc_expto0(trans_back_norm_pc_spGl4)

#Import luciferase assay df

int_luc <- read_csv('Int_luc_4var_pncontrols.csv') %>%
  rename(conc = Conc)

#Will leave manipulation of all variants from luciferase assay in function, the
#rbind function at the end determines which are combined and compared to other
#assays

int_untidy_conc_expto0 <- function(df) {
  nc_Chr9 <- df %>%
    select(conc, nc_Chr9a, nc_Chr9b, nc_Chr9c) %>%
    mutate(ave_lum = (nc_Chr9a + nc_Chr9b + nc_Chr9c)/3) %>%
    select(conc, ave_lum) %>%
    mutate(ave_exp_to0 = ave_lum/ave_lum[1]) %>%
    select(-ave_lum) %>%
    mutate(name = 'subpool5_no_site_no_site_no_site_no_site_no_site_no_site_v chr9') %>%
    mutate(subpool = 'subpool5')
  pc_CRE <- df %>%
    select(conc, pc_CREa, pc_CREb, pc_CREc) %>%
    mutate(ave_lum = (pc_CREa + pc_CREb + pc_CREc)/3) %>%
    select(conc, ave_lum) %>%
    mutate(ave_exp_to0 = ave_lum/ave_lum[1])  %>%
    select(-ave_lum) %>%
    mutate(name = 'pGL4.29 Promega 1-63 + 1-87') %>%
    mutate(subpool = 'control')
  nc_CREmut <- df %>%
    select(conc, nc_CREmuta, nc_CREmutb, nc_CREmutc) %>%
    mutate(ave_lum = (nc_CREmuta + nc_CREmutb + nc_CREmutc)/3) %>%
    select(conc, ave_lum) %>%
    mutate(ave_exp_to0 = ave_lum/ave_lum[1]) %>%
    select(-ave_lum) %>%
    mutate(name = 'subpool5_no_site_no_site_no_site_no_site_no_site_no_site_s pGl4') %>%
    mutate(subpool = 'subpool5')
  var_11 <- df %>%
    select(conc, var_11a, var_11b, var_11c) %>%
    mutate(ave_lum = (var_11a + var_11b + var_11c)/3) %>%
    select(conc, ave_lum) %>%
    mutate(ave_exp_to0 = ave_lum/ave_lum[1])  %>%
    select(-ave_lum) %>%
    mutate(name = 'subpool3_2BS 16 bp spacing consensus+flank x2_dist_79_s pGl4') %>%
    mutate(subpool = 'subpool3')
  var_14 <- df %>%
    select(conc, var_14a, var_14b, var_14c) %>%
    mutate(ave_lum = (var_14a + var_14b + var_14c)/3) %>%
    select(conc, ave_lum) %>%
    mutate(ave_exp_to0 = ave_lum/ave_lum[1])  %>%
    select(-ave_lum) %>%
    mutate(name = 'subpool3_2BS 6 bp spacing consensus+flank x2_dist_9_s pGl4') %>%
    mutate(subpool = 'subpool3')
  var_17 <- df %>%
    select(conc, var_17a, var_17b, var_17c) %>%
    mutate(ave_lum = (var_17a + var_17b + var_17c)/3) %>%
    select(conc, ave_lum) %>%
    mutate(ave_exp_to0 = ave_lum/ave_lum[1])  %>%
    select(-ave_lum) %>%
    mutate(name = 'subpool3_2BS 11 bp spacing consensus+flank x2_dist_14_s pGl4') %>%
    mutate(subpool = 'subpool3')
  var_18 <- df %>%
    select(conc, var_18a, var_18b, var_18c) %>%
    mutate(ave_lum = (var_18a + var_18b + var_18c)/3) %>%
    select(conc, ave_lum) %>%
    mutate(ave_exp_to0 = ave_lum/ave_lum[1])  %>%
    select(-ave_lum) %>%
    mutate(name = 'subpool5_no_site_consensus_weak_consensus_no_site_weak_v chr9') %>%
    mutate(subpool = 'subpool3')
  all_var <- rbind(var_11, var_14, var_17, var_18, pc_CRE) %>%
    mutate(assay = 'int_luc')
  return(all_var)
}

luc_comp <- int_untidy_conc_expto0(int_luc)

comp_luc_trans <- rbind(trans_comp, luc_comp)

p_luc_comp <- ggplot(luc_comp, aes(conc, ave_exp_to0, color = name)) +
  geom_point() +
  geom_line() +
  scale_color_viridis(discrete = TRUE) + 
  xlab('Forskolin (µM)') +
  ylab('Normalized luminescence\nto 0 µM Forsk (AU)') +
  background_grid() + 
  ggtitle('Integrated Luminescence')

save_plot('plots/p_luc_comp.png', p_luc_comp, scale = 1.1,
          base_width = 10, base_height = 3)

p_trans_comp <- ggplot(trans_comp, aes(conc, ave_exp_to0, color = name)) +
  geom_point() +
  geom_line() +
  scale_color_viridis(discrete = TRUE) + 
  xlab('Forskolin (µM)') +
  ylab('Normalized expression\nto 0 µM Forsk (AU)') +
  background_grid() + 
  ggtitle('Transient MPRA         ')

save_plot('plots/p_trans_comp.png', p_trans_comp, scale = 1.1,
          base_width = 10, base_height = 3)


#K-means clustering on data-----------------------------------------------------

#Can I bin out noisy data manually just from induction level?

ind_less2 <- trans_back_norm_rep_0_22 %>%
  filter(induction_of_ave < 2) %>%
  var_conc_exp()

ind_more2 <- trans_back_norm_rep_0_22 %>%
  filter(induction_of_ave >= 2) %>%
  var_conc_exp()

ggplot(ind_less2, 
       aes(conc, ave_ratio_norm, color = name)) +
  geom_point(show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  xlab('Forskolin µM') +
  ylab('Average background-normalized\nsum RNA/DNA') +
  panel_border() +
  annotation_logticks(sides = 'b') +
  background_grid(major = 'xy', minor = 'none')


#Can I use k-means to bin out noisy data?

min_log10_back_norm_rep_0_22 <- log10_trans_back_norm_rep_0_22 %>%
  filter(subpool == 'subpool5') %>%
  select(ave_ratio_0_norm, ave_ratio_2_5_norm, ave_ratio_2_4_norm,
         ave_ratio_2_3_norm, ave_ratio_2_2_norm, ave_ratio_2_1_norm,
         ave_ratio_20_norm, ave_ratio_22_norm) %>%
  scale()

wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data) - 1) * sum(apply(data, 2, var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot(min_log10_back_norm_rep_0_22)

save_plot('plots/wsskmeans.png', wssplot(min_log10_back_norm_rep_0_22),
          base_height = 4.4, base_width = 5)

set.seed(1234)
min_km <- kmeans(min_log10_back_norm_rep_0_22, 7, nstart = 25)

min_km

s5_log10_bn_rep_0_22 <- log10_trans_back_norm_rep_0_22 %>%
  filter(subpool == 'subpool5')

clust_s5_bn_conc_rep_0_22 <- augment(min_km, 
                                     min_back_norm_rep_0_22) %>%
  rename(cluster = .cluster) %>%
  select(cluster) %>%
  cbind(s5_log10_bn_rep_0_22) %>%
  select(subpool, name, most_common, cluster) %>%
  left_join(trans_back_norm_conc, by = c('subpool', 'name', 'most_common'))

ggplot(filter(clust_s5_bn_conc_rep_0_22, cluster == 5 | cluster == 7), 
       aes(conc, ave_ratio_norm, color = name)) +
  geom_point(show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ cluster) +
  xlab('Forskolin µM') +
  ylab('Average background-normalized\nsum RNA/DNA') +
  panel_border() +
  annotation_logticks(sides = 'b') +
  background_grid(major = 'xy', minor = 'none')


#Titration Plots----------------------------------------------------------------

#Subtract ave expression at 0 µM from each variant's ave expression

trans_back_0_norm_conc <- trans_back_norm_pc_spGl4 %>%
  mutate(ave_ratio_2_5_norm = ave_ratio_2_5_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_2_4_norm = ave_ratio_2_4_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_2_3_norm = ave_ratio_2_3_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_2_2_norm = ave_ratio_2_2_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_2_1_norm = ave_ratio_2_1_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_20_norm = ave_ratio_20_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_22_norm = ave_ratio_22_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_0_norm = ave_ratio_0_norm - ave_ratio_0_norm) %>%
  var_conc_exp()

trans_back_norm_pc_spGl4_conc <- trans_back_norm_pc_spGl4_conc %>%
  mutate(conc = if_else(conc == 0,
                        2^-7,
                        conc)) %>%
  mutate(conc = log2(conc))

#plots of titration overall with 1 control (duplicated pGl4, original does not
#have high expression)

p_titr_pc_back <- trans_back_norm_pc_spGl4_conc %>%
  ggplot(aes(conc, ave_ratio_norm)) +
  geom_line(aes(group = name), alpha = 0.1) +
  geom_point(data = filter(trans_back_norm_pc_spGl4_conc, 
                           startsWith(name, 
                                      'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')),
             color = 'darkgoldenrod1', shape = 19, stroke = 1.25) +
  geom_line(data = filter(trans_back_norm_pc_spGl4_conc, 
                          startsWith(name, 
                                     'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')),
            color = 'darkgoldenrod1', size = 1.25) +
  geom_point(data = filter(trans_back_norm_pc_spGl4_conc, 
                           startsWith(name, 
                                      'pGL4.29 Promega 1-63 + 1-87')),
             color = 'firebrick2', shape = 19, stroke = 1.25) +
  geom_line(data = filter(trans_back_norm_pc_spGl4_conc, 
                          startsWith(name, 
                                     'pGL4.29 Promega 1-63 + 1-87')),
            color = 'firebrick2', size = 1.25) +
  ylab('Average normalized\nexpression (a.u.)') +
  annotation_logticks(sides = 'b') +
  scale_x_continuous(breaks = (-7:2), 'log2 forskolin (µM)')

save_plot('plots/p_titr_pc_back.pdf', p_titr_pc_back, scale = 1.3, 
          base_width = 4.5, base_height = 2.5)


#Fitting nested data

library(minpack.lm)

m_m_model_nlslm <- function(df) {
  m_m_nlslm <- nlsLM(
    ave_ratio_norm ~ (max_ave_ratio_norm * conc^n)/(conc_half_max^n + conc^n),
    data = df, 
    start = list(conc_half_max = (2^-3), max_ave_ratio_norm = 2, n = 1))
  return(m_m_nlslm)
}

m_m_nest_coef <- function(df1) {
  add_coef_unnest <- df1 %>%
    mutate(results = map(m_m_fit, tidy)) %>%
    select(-m_m_fit, -data) %>%
    unnest()
}

m_m_nest_pred_resid <- function(df1) {
  pred <- df1 %>%
    mutate(m_m_pred = map2(data, m_m_fit, add_predictions)) %>%
    select(-data, -m_m_fit) %>%
    unnest()
  resid <- df1 %>%
    mutate(m_m_resids = map2(data, m_m_fit, add_residuals)) %>%
    select(-data, -m_m_fit) %>%
    unnest()
  data <- df1 %>%
    select(-m_m_fit) %>%
    unnest()
  data_pred <- left_join(data, pred, 
                         by = c('subpool', 'name', 'most_common', 'background',
                                'conc', 'ave_ratio_norm'))
  data_pred_resid <- left_join(data_pred, resid,
                               by = c('subpool', 'name', 'most_common', 
                                      'background', 'conc', 'ave_ratio_norm')) %>%
    ungroup()
  return(data_pred_resid)
}


#Select higher expressing variants to fit curves to, had hard time fitting
#nested data to 2500-4172 variants, but these drove lower than 0.75 expression 
#above 0, would not expect a well-fitting curve for these lower-expressing 
#variants

trans_back_0_norm_conc_nest <- trans_back_0_norm_conc %>%
  select(-ave_barcode) %>%
  filter(!grepl('^subpool5_no_site_no_site_no_site_no_site_no_site_no_site', 
                name)) %>%
  filter(subpool == 'subpool5') %>%
  arrange(desc(ave_ratio_norm)) %>%
  group_by(subpool, name, most_common, background) %>%
  nest() %>%
  slice(1:1500)

m_m_nest_fit_nlslm <- trans_back_0_norm_conc_nest %>%
  mutate(m_m_fit = map(trans_back_0_norm_conc_nest$data, m_m_model_nlslm))

m_m_coef_nlslm <- m_m_nest_coef(m_m_nest_fit_nlslm)

#Test fits

hillcoef <- m_m_coef_nlslm %>%
  filter(term == 'n') %>%
  mutate(n_rse = std.error/estimate) %>%
  filter(n_rse <= 0.25) %>%
  rename(n = estimate) %>%
  select(-term, -std.error, -statistic, -p.value) %>%
  left_join(s5_untidy,
            by = c('most_common', 'background')) %>%
  select(-subpool) %>%
  ungroup()

hill_EC50 <- m_m_coef_nlslm %>%
  filter(term == 'conc_half_max') %>%
  mutate(EC50_rse = std.error/estimate) %>%
  filter(EC50_rse <= 0.25) %>%
  rename(EC50 = estimate) %>%
  select(-term, -std.error, -statistic, -p.value) %>%
  inner_join(hillcoef, by = c('name', 'most_common', 'background')) %>%
  ungroup()

m_m_p_r <- m_m_nest_pred_resid(m_m_nest_fit_nlslm) %>%
  rename(ave_ratio_norm_0 = ave_ratio_norm)

m_m_EC50_n_p_r <- left_join(hill_EC50, m_m_p_r, 
                            by = c('subpool', 'name', 'most_common', 
                                   'background', 'conc'))

p_resid_dens <- m_m_EC50_n_p_r %>%
  ggplot(aes(resid)) +
  geom_density(kernel = 'gaussian') +
  xlab('expression residual')

save_plot('plots/p_resid_dens.pdf', p_resid_dens, scale = 1.3, base_width = 2.5,
          base_height = 2)

p_resid_distr <- ggplot(m_m_EC50_n_p_r, aes(ave_ratio_norm_0, resid)) +
  geom_point(alpha = 0.1, size = 0.75) +
  xlab('Average normalized\nexpression (a.u.)\n-exp. at 0 µM')

save_plot('plots/p_resid_distr.pdf', p_resid_distr, scale = 1.3, 
          base_width = 3.5, base_height = 2)

p_EC50_rse_hist <- m_m_EC50_n_p_r %>%
  filter(conc == 4) %>%
  ggplot(aes(EC50_rse)) +
  geom_density(kernel = 'gaussian') +
  scale_x_continuous(breaks = c(0, 0.1, 0.2))

save_plot('plots/p_EC50_rse_hist.pdf', p_EC50_rse_hist, scale = 1.3, base_height = 2,
          base_width = 2.5)

p_n_rse_hist <- m_m_EC50_n_p_r %>%
  filter(conc == 4) %>%
  ggplot(aes(n_rse)) +
  geom_density(kernel = 'gaussian') +
  scale_x_continuous(breaks = c(0, 0.1, 0.2))

save_plot('plots/p_n_rse_hist.pdf', p_n_rse_hist, scale = 1.3, base_height = 2,
          base_width = 2.5)


#Figure 2-----------------------------------------------------------------------

#plot concentration vs. expression of variants that were fit for site number
#analysis in figure 2

p_response_cons_back <- m_m_EC50_n_p_r %>%
  mutate(conc = if_else(conc == 0,
                        0.005,
                        conc)) %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  filter(consensus > 1 & weak == 0) %>%
  ggplot(aes(conc, ave_ratio_norm, color = as.factor(consensus))) +
  facet_grid(background ~ .) +
  scale_x_log10(name = 'µM forskolin') +
  scale_y_log10(name = 'Average normalized expression (a.u.)') +
  annotation_logticks(sides = 'bl') +
  geom_point(alpha = 0.75, size = 1) +
  geom_smooth(alpha = 0.3, method = 'loess') +
  scale_color_viridis(discrete = TRUE, name = 'Consensus CREs') +
  panel_border(colour = 'black') +
  theme(legend.position = 'top',
    strip.background = element_rect(colour="black", fill="white"))

save_plot('plots/p_response_cons_back.pdf', p_response_cons_back, scale = 1.3,
          base_width = 4, base_height = 4)

p_n_cons_back <- hill_EC50 %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  filter(weak == 0 & consensus > 1 & conc == 4) %>%
  ggplot(aes(as.factor(consensus), n)) +
  facet_grid(background ~ .) +
  geom_hline(yintercept = 1.6, color = 'grey 60', linetype = 2) +
  geom_boxplot(outlier.size = 0.5, size = 0.3) +
  scale_color_viridis(discrete = TRUE) +
  xlab('consensus sites') +
  panel_border(colour = 'black') +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

save_plot('plots/p_n_cons_back.pdf', p_n_cons_back, scale = 1.3, 
          base_width = 2.25, base_height = 2.25)

p_EC50_cons_back <- hill_EC50 %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  filter(weak == 0 & consensus > 1 & conc == 4) %>%
  ggplot(aes(as.factor(consensus), EC50)) +
  facet_grid(background ~ .) +
  geom_boxplot(outlier.size = 0.5, size = 0.3) +
  scale_color_viridis(discrete = TRUE) +
  xlab('consensus sites') +
  panel_border(colour = 'black') +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

save_plot('plots/p_EC50_cons_back.pdf', p_EC50_cons_back, scale = 1.3, 
          base_width = 2.25, base_height = 2.25)

p_induct_cons_back <- hill_EC50 %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  filter(weak == 0 & consensus > 1 & conc == 4) %>%
  ggplot(aes(as.factor(consensus), induction)) +
  facet_grid(background ~ .) +
  geom_boxplot(outlier.size = 0.5, size = 0.3) +
  scale_color_viridis(discrete = TRUE) +
  xlab('consensus sites') +
  panel_border(colour = 'black') +
  background_grid(major = 'y') +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

save_plot('plots/p_induct_cons_back.pdf', p_induct_cons_back, scale = 1.3, 
          base_width = 2.25, base_height = 2.25)


#plot coefficients (n and EC50) and their relationship to architectures

p_n_consensus <- hill_EC50 %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  filter(weak == 0 & consensus > 1 & conc == 4) %>%
  ggplot(aes(as.factor(consensus), n)) +
  geom_jitter(aes(color = background), size = 1,
              position=position_jitter(width=0.3, height=0), alpha = 0.5) +
  geom_boxplot(outlier.shape=NA, alpha = 0, outlier.size = 0.5, size = 0.5, 
               outlier.alpha = 0.5) +
  scale_color_manual(values = cbPalette3, name = 'background') +
  xlab('consensus sites') +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

save_plot('plots/p_n_consensus.pdf', p_n_consensus, scale = 1.3, 
          base_width = 3.25, base_height = 2)

p_n_expression <- hill_EC50 %>%
  filter(conc == 4) %>%
  ggplot(aes(ave_ratio_norm, n)) +
  geom_point(alpha = 0.3, size = 0.75) +
  scale_x_log10() +
  annotation_logticks(sides = 'b') +
  xlab('induced expression') +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

save_plot('plots/p_n_expression.pdf', p_n_expression, scale = 1.3, 
          base_width = 2.5, base_height = 1.5)

p_n_rse <- hill_EC50 %>%
  filter(conc == 4) %>%
  ggplot(aes(n_rse, n)) +
  geom_point(alpha = 0.3, size = 0.75) +
  xlab('n rse') +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

save_plot('plots/p_n_rse.pdf', p_n_rse, scale = 1.3, 
          base_width = 2.5, base_height = 1.5)

p_n_consensus_weak <- hill_EC50 %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  filter(consensus > 1 & conc == 4) %>%
  ggplot(aes(as.factor(consensus), n)) +
  geom_boxplot(aes(fill = as.factor(weak)), outlier.size = 1, size = 0.3, 
               outlier.shape = 21, outlier.alpha = 0.75, 
               position = position_dodge(0.75)) +
  scale_fill_manual(values = cbPalette7_grad_light, name = 'weak sites') +
  xlab('consensus sites') +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

save_plot('plots/p_n_consensus_weak.pdf', p_n_consensus_weak, scale = 1.3, 
          base_width = 3.25, base_height = 2)

test <- hill_EC50 %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  filter(consensus == 3 & conc == 4 & site3 == 'consensus' & site5 == 'consensus') %>%
  ggplot(aes(as.factor(weak), n)) +
  facet_wrap(~ background) +
  geom_jitter(size = 1, color = 'dodgerblue2',
              position=position_jitter(width=0.3, height=0)) +
  geom_boxplot(outlier.shape=NA, alpha = 0, outlier.size = 0.5, size = 0.5, 
               outlier.alpha = 0.5) +
  scale_fill_manual(values = cbPalette7_grad_light, name = 'weak sites') +
  xlab('weak sites') +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))



p_s5_num_cons_back_EC50 <- hill_EC50 %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  filter(consensus > 1 & weak == 0 & conc == 4) %>%
  ggplot(aes(as.factor(consensus), EC50)) +
  facet_grid(~ background) +
  geom_boxplot(outlier.size = 0.75, size = 0.3, outlier.alpha = 1) +
  panel_border(colour = 'black') +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  background_grid(major = 'y', minor = 'none') + 
  ylab('EC50') +
  xlab('Number of consensus sites')

save_plot('plots/p_s5_num_cons_back_EC50.pdf', 
          p_s5_num_cons_back_EC50, scale = 1.3,
          base_width = 5, base_height = 2)
  
p_EC50_cons_weak <- hill_EC50 %>%
  filter(consensus > 1 & conc == 4) %>%
  ggplot(aes(as.factor(consensus), EC50)) +
  geom_boxplot(aes(fill = as.factor(weak)), outlier.size = 1, size = 0.3, 
               outlier.shape = 21, outlier.alpha = 0.75, 
               position = position_dodge(0.75)) +
  scale_fill_manual(values = cbPalette7_grad_light, name = 'weak sites') +
  xlab('consensus sites') +
  background_grid(major = 'y') +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

save_plot('plots/p_EC50_cons_weak.pdf', p_EC50_cons_weak, scale = 1.3,
          base_width = 4, base_height = 2)

p_EC50_expression <- hill_EC50 %>%
  filter(conc == 4) %>%
  ggplot(aes(ave_ratio_norm, EC50)) +
  geom_point(alpha = 0.3, size = 0.75) +
  scale_x_log10() +
  annotation_logticks(sides = 'b') +
  xlab('induced expression') +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

save_plot('plots/p_EC50_expression.pdf', p_EC50_expression, scale = 1.3,
          base_width = 2.5, base_height = 1.5)

p_EC50_n <- hill_EC50 %>%
  filter(conc == 4) %>%
  ggplot(aes(n, EC50)) +
  geom_point(alpha = 0.3, size = 0.75) +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

save_plot('plots/p_EC50_n.pdf', p_EC50_n, scale = 1.3,
          base_width = 2.5, base_height = 1.5)

p_EC50_rse <- hill_EC50 %>%
  filter(conc == 4) %>%
  ggplot(aes(EC50_rse, EC50)) +
  geom_point(alpha = 0.3, size = 0.75) +
  xlab('EC50 rse') +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

save_plot('plots/p_EC50_rse.pdf', p_EC50_rse, scale = 1.3,
          base_width = 2.5, base_height = 1.5)

#Fit sp3 induction curves-------------------------------------------------------

trans_back_0_norm_conc_nest_sp3 <- trans_back_0_norm_conc %>%
  filter(startsWith(name, 'subpool3')) %>%
  select(-ave_barcode) %>%
  arrange(desc(ave_ratio_norm)) %>%
  group_by(subpool, name, most_common, background) %>%
  nest() %>%
  slice(1:800)

m_m_nest_fit_nlslm_sp3 <- trans_back_0_norm_conc_nest_sp3 %>%
  mutate(m_m_fit = map(trans_back_0_norm_conc_nest_sp3$data, m_m_model_nlslm))

m_m_coef_nlslm_sp3 <- m_m_nest_coef(m_m_nest_fit_nlslm_sp3)

#Test fits

hillcoef_sp3 <- m_m_coef_nlslm_sp3 %>%
  filter(term == 'n') %>%
  mutate(n_rse = std.error/estimate) %>%
  filter(n_rse <= 0.25) %>%
  rename(n = estimate) %>%
  select(-term, -std.error, -statistic, -p.value) %>%
  left_join(s3_untidy,
            by = c('most_common', 'background')) %>%
  select(-subpool) %>%
  ungroup()

hill_EC50_sp3 <- m_m_coef_nlslm_sp3 %>%
  filter(term == 'conc_half_max') %>%
  mutate(EC50_rse = std.error/estimate) %>%
  filter(EC50_rse <= 0.25) %>%
  rename(EC50 = estimate) %>%
  select(-term, -std.error, -statistic, -p.value) %>%
  inner_join(hillcoef_sp3, by = c('name', 'most_common', 'background')) %>%
  ungroup()

m_m_p_r_sp3 <- m_m_nest_pred_resid(m_m_nest_fit_nlslm_sp3) %>%
  rename(ave_ratio_norm_0 = ave_ratio_norm)

m_m_EC50_n_p_r_sp3 <- left_join(hill_EC50_sp3, m_m_p_r_sp3, 
                            by = c('subpool', 'name', 'most_common', 
                                   'background', 'conc'))

#residuals and error plots

p_resid_distr_sp3 <- m_m_EC50_n_p_r_sp3 %>%
  ggplot(aes(ave_ratio_norm_0, resid)) +
  geom_point(alpha = 0.1, size = 0.75) +
  xlab('Average normalized\nexpression (a.u.)\n-exp. at 0 µM')

save_plot('plots/p_resid_distr_sp3.pdf', p_resid_distr_sp3, scale = 1.3, 
          base_width = 3.5, base_height = 2)

p_EC50_rse_hist_sp3 <- m_m_EC50_n_p_r_sp3 %>%
  filter(conc == 4) %>%
  ggplot(aes(EC50_rse)) +
  geom_density(kernel = 'gaussian') +
  scale_x_continuous(breaks = c(0, 0.1, 0.2))

save_plot('plots/p_EC50_rse_hist_sp3.pdf', p_EC50_rse_hist_sp3, scale = 1.3, 
          base_height = 2, base_width = 2.5)

p_n_rse_hist_sp3 <- m_m_EC50_n_p_r_sp3 %>%
  filter(conc == 4) %>%
  ggplot(aes(n_rse)) +
  geom_density(kernel = 'gaussian') +
  scale_x_continuous(breaks = c(0, 0.1, 0.2))

save_plot('plots/p_n_rse_hist_sp3.pdf', p_n_rse_hist_sp3, scale = 1.3, 
          base_height = 2, base_width = 2.5)

test <- m_m_EC50_n_p_r_sp3 %>%
  group_by(subpool, name, most_common, background) %>%
  nest() %>%
  sample_n(8) %>%
  unnest() %>%
  ggplot(aes(x = conc, color = name)) +
  geom_point(aes(y = ave_ratio_norm_0)) +
  geom_line(aes(y = ave_ratio_norm_0)) +
  geom_point(aes(y = pred), shape = 21) +
  geom_line(aes(y = pred), linetype = 2) +
  scale_color_viridis(discrete = TRUE)

#plot distance and spacing features

p_s3_n_dist_vchr5 <- m_m_EC50_n_p_r_sp3 %>%
  mutate(bin = cut(dist, seq(from = 0, to = 140, by = 20),
                   labels = c('0-20', '20-40', '40-60', '60-80',
                              '80-100', '100-120', '120-140'))) %>%
  filter(bin != '120-140' & conc == 4 & spacing != 70 & background == 'v chr5') %>%
  ggplot(aes(bin, n)) +
  geom_jitter(aes(color = as.factor(spacing)), 
              position=position_jitter(width=0.3, height=0), alpha = 0.75,
              size = 0.5) +
  geom_boxplot(outlier.shape=NA, size = 0.3, position = position_dodge(1),
               show.legend = FALSE, alpha = 0) +
  scale_color_manual(values = cbPalette4, name = 'spacing (bp)') +
  theme(legend.position = 'top', axis.ticks.x = element_blank(), 
        strip.background = element_rect(colour="black", fill="white"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  panel_border(colour = 'black') +
  ylab('Hill n') +
  xlab('Distance along background (bp)')

save_plot('plots/p_s3_n_dist_vchr5.pdf', p_s3_n_dist_vchr5, scale = 1.3,
          base_width = 3, base_height = 2.25)

p_s3_EC50_dist_vchr5 <- m_m_EC50_n_p_r_sp3 %>%
  mutate(bin = cut(dist, seq(from = 0, to = 140, by = 20),
                   labels = c('0-20', '20-40', '40-60', '60-80',
                              '80-100', '100-120', '120-140'))) %>%
  filter(bin != '120-140' & conc == 4 & spacing != 70 & background == 'v chr5') %>%
  ggplot(aes(bin, EC50)) +
  geom_jitter(aes(color = as.factor(spacing)), 
              position=position_jitter(width=0.3, height=0), alpha = 0.75,
              size = 0.5) +
  geom_boxplot(outlier.shape=NA, size = 0.3, position = position_dodge(1),
               show.legend = FALSE, alpha = 0) +
  scale_color_manual(values = cbPalette4, name = 'spacing (bp)') +
  theme(legend.position = 'top', axis.ticks.x = element_blank(), 
        strip.background = element_rect(colour="black", fill="white"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  panel_border(colour = 'black') +
  ylab('EC50 (µM)') +
  xlab('Distance along background (bp)')

save_plot('plots/p_s3_EC50_dist_vchr5.pdf', p_s3_EC50_dist_vchr5, scale = 1.3,
          base_width = 3, base_height = 2.25)

p_s3_n_space_vchr5 <- m_m_EC50_n_p_r_sp3 %>%
  filter(conc == 4 & background == 'v chr5') %>%
  ggplot(aes(as.factor(spacing), n)) +
  geom_boxplot(outlier.size = 1, size = 0.3, position = position_dodge(1)) +
  theme(axis.ticks.x = element_blank(), 
        strip.background = element_rect(colour="black", fill="white")) + 
  panel_border(colour = 'black') +
  ylab('Hill n') +
  xlab('Spacing (bp)')

save_plot('plots/p_s3_n_space_vchr5.pdf', p_s3_n_space_vchr5, scale = 1.3,
          base_width = 2.5, base_height = 1.75)

p_s3_EC50_space_vchr5 <- m_m_EC50_n_p_r_sp3 %>%
  filter(conc == 4 & background == 'v chr5' & dist < 40) %>%
  ggplot(aes(as.factor(spacing), EC50)) +
  geom_boxplot(outlier.size = 1, size = 0.3, position = position_dodge(1)) +
  theme(axis.ticks.x = element_blank(), 
        strip.background = element_rect(colour="black", fill="white")) + 
  panel_border(colour = 'black') +
  ylab('EC50 (µM)') +
  xlab('Spacing (bp)')

save_plot('plots/p_s3_EC50_space_vchr5.pdf', p_s3_EC50_space_vchr5, scale = 1.3,
          base_width = 2.5, base_height = 1.75)

test <- m_m_EC50_n_p_r_sp3 %>%
  filter(conc == 4 & background == 'v chr5') %>%
  ggplot(aes(ave_ratio_norm_0, resid)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  annotation_logticks(sides = 'b')


#Combine and compare expression across a different set of concentrations--------

trans_back_norm_rep_0_64 <- read_tsv(
  '../20170630_tlib/trans_back_norm_rep_0_64.txt')

#Compare background-normalized expression at conc 0, 1 and 4 µM that were tested
#in both experiments

conc_0_1_4_comp <- function(df1, df2) {
  two_df <- inner_join(df1, df2, 
                       by = c('subpool', 'name', 'most_common', 'background'),
                       suffix = c('_0921', '_0631')) %>%
    select(subpool, name, most_common, background, ave_ratio_0_norm_0921, 
           ave_ratio_0_norm_0631, ave_ratio_20_norm, ave_ratio_1_norm, 
           ave_ratio_22_norm, ave_ratio_4_norm) %>%
    rename(ave_ratio_1_norm_0921 = ave_ratio_20_norm) %>%
    rename(ave_ratio_1_norm_0631 = ave_ratio_1_norm) %>%
    rename(ave_ratio_4_norm_0921 = ave_ratio_22_norm) %>%
    rename(ave_ratio_4_norm_0631 = ave_ratio_4_norm)
}

trans_back_norm_comp_0921_0631 <- conc_0_1_4_comp(
  trans_back_norm_rep_0_22, trans_back_norm_rep_0_64) %>%
  var_log10()

trans_comp_lm_0 <- lm(ave_ratio_0_norm_0921 ~ ave_ratio_0_norm_0631, 
                      data = trans_back_norm_comp_0921_0631)

trans_comp_lm_1 <- lm(ave_ratio_1_norm_0921 ~ ave_ratio_1_norm_0631, 
                      data = trans_back_norm_comp_0921_0631)

trans_comp_lm_4 <- lm(ave_ratio_4_norm_0921 ~ ave_ratio_4_norm_0631, 
                      data = trans_back_norm_comp_0921_0631)

pre_res_trans_int <- function(df1, x) {
  df2 <- df1 %>%
    add_predictions(x)
  df3 <- df2 %>%
    add_residuals(x)
  return(df3)
  print('processed pre_res_trans_int(df1, df2) in order of (data, model)')
}

trans_comp_p_r_0 <- pre_res_trans_int(trans_back_norm_comp_0921_0631, trans_comp_lm_0)

trans_comp_p_r_1 <- pre_res_trans_int(trans_back_norm_comp_0921_0631, trans_comp_lm_1)

trans_comp_p_r_4 <- pre_res_trans_int(trans_back_norm_comp_0921_0631, trans_comp_lm_4)

p_corr_0921_0631_0 <- ggplot(data = NULL, aes(ave_ratio_0_norm_0921, ave_ratio_0_norm_0631,
                                              color = subpool)) +
  geom_point(data = filter(trans_back_norm_comp_0921_0631, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(trans_back_norm_comp_0921_0631, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = trans_back_norm_comp_0921_0631,
                 alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  annotation_logticks(scaled = TRUE) +
  xlab("Ave log10 variant norm.\nsum RNA/DNA 170921") +
  ylab("Ave log10 variant norm.\nsum RNA/DNA 170631") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7), 
                     limits = c(-0.5, 7.5)) + 
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7), 
                     limits = c(-0.5, 7.5)) + 
  annotate("text", x = 1, y = 6,
           label = paste('r =', round(
             cor(
               trans_back_norm_comp_0921_0631$ave_ratio_0_norm_0921, 
               trans_back_norm_comp_0921_0631$ave_ratio_0_norm_0631,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) +
  annotate("text", x = 1, y = 5, color = '#2D708EFF',
           label = paste('r =', round(
             cor(
               filter(trans_back_norm_comp_0921_0631, subpool == 'subpool3')$ave_ratio_0_norm_0921, 
               filter(trans_back_norm_comp_0921_0631, subpool == 'subpool3')$ave_ratio_0_norm_0631,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) + 
  annotate("text", x = 1, y = 4, color = '#482677FF',
           label = paste('r =', round(
             cor(
               filter(trans_back_norm_comp_0921_0631, subpool == 'subpool5')$ave_ratio_0_norm_0921, 
               filter(trans_back_norm_comp_0921_0631, subpool == 'subpool5')$ave_ratio_0_norm_0631,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  )


#Negatives plots----------------------------------------------------------------------------

log_neg_cont <- function(df1) {
  neg <- filter(df1, 
         grepl(
           'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
           name)) %>%
    mutate(name = str_sub(name, 58,)) %>%
    mutate(name = gsub('Smith R. Vista chr9:83712599-83712766','vista chr9', name
    ),
    name = gsub('Vista Chr5:88673410-88674494', 'vista chr5', name
    ),
    name = str_sub(name, 1,13)
    )
  cont <- filter(df1, grepl('control', name))
  neg_cont <- bind_rows(neg, cont)
  return(neg_cont)
}

log_rep_0_22_A_B_neg_cont <- log_neg_cont(log_rep_0_22_A_B)

p_nc_var_rep_0 <- ggplot(NULL, aes(log_ratio_0A, log_ratio_0B)) +
  geom_point(data = filter(log_rep_0_22_A_B_neg_cont, subpool == 'control'),
             alpha = 0.6) +
  geom_point(data = filter(log_rep_0_22_A_B_neg_cont, 
                           subpool == 'subpool5'), 
             color = 'red', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expr. Rep. 1") +
  ylab("Variant log10 Expr. Rep. 2") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, .5)) + 
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, .5)) + 
  annotate("text", x = -1, y = 0,
           label = paste(
             'r =', round(
               cor(
                 log_rep_0_22_A_B_neg_cont$log_ratio_0A,
                 log_rep_0_22_A_B_neg_cont$log_ratio_0B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_nc_var_rep_2_5 <- ggplot(NULL, aes(log_ratio_2_5A, log_ratio_2_5B)) +
  geom_point(data = filter(log_rep_0_22_A_B_neg_cont, subpool == 'control'),
             alpha = 0.6) +
  geom_point(data = filter(log_rep_0_22_A_B_neg_cont, 
                           subpool == 'subpool5'), 
             color = 'red', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expr. Rep. 1") +
  ylab("Variant log10 Expr. Rep. 2") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, .5)) + 
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, .5)) + 
  annotate("text", x = -1, y = 0,
           label = paste(
             'r =', round(
               cor(
                 log_rep_0_22_A_B_neg_cont$log_ratio_2_5A,
                 log_rep_0_22_A_B_neg_cont$log_ratio_2_5B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_nc_var_rep_2_4 <- ggplot(NULL, aes(log_ratio_2_4A, log_ratio_2_4B)) +
  geom_point(data = filter(log_rep_0_22_A_B_neg_cont, subpool == 'control'),
             alpha = 0.6) +
  geom_point(data = filter(log_rep_0_22_A_B_neg_cont, 
                           subpool == 'subpool5'), 
             color = 'red', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expr. Rep. 1") +
  ylab("Variant log10 Expr. Rep. 2") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, .5)) + 
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, .5)) + 
  annotate("text", x = -1, y = 0,
           label = paste(
             'r =', round(
               cor(
                 log_rep_0_22_A_B_neg_cont$log_ratio_2_4A,
                 log_rep_0_22_A_B_neg_cont$log_ratio_2_4B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_nc_var_rep_2_3 <- ggplot(NULL, aes(log_ratio_2_3A, log_ratio_2_3B)) +
  geom_point(data = filter(log_rep_0_22_A_B_neg_cont, subpool == 'control'),
             alpha = 0.6) +
  geom_point(data = filter(log_rep_0_22_A_B_neg_cont, 
                           subpool == 'subpool5'), 
             color = 'red', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expr. Rep. 1") +
  ylab("Variant log10 Expr. Rep. 2") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, .5)) + 
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, .5)) + 
  annotate("text", x = -1, y = 0,
           label = paste(
             'r =', round(
               cor(
                 log_rep_0_22_A_B_neg_cont$log_ratio_2_3A,
                 log_rep_0_22_A_B_neg_cont$log_ratio_2_3B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_nc_var_rep_2_2 <- ggplot(NULL, aes(log_ratio_2_2A, log_ratio_2_2B)) +
  geom_point(data = filter(log_rep_0_22_A_B_neg_cont, subpool == 'control'),
             alpha = 0.6) +
  geom_point(data = filter(log_rep_0_22_A_B_neg_cont, 
                           subpool == 'subpool5'), 
             color = 'red', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expr. Rep. 1") +
  ylab("Variant log10 Expr. Rep. 2") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, .5)) + 
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, .5)) + 
  annotate("text", x = -1, y = 0,
           label = paste(
             'r =', round(
               cor(
                 log_rep_0_22_A_B_neg_cont$log_ratio_2_2A,
                 log_rep_0_22_A_B_neg_cont$log_ratio_2_2B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_nc_var_rep_2_1 <- ggplot(NULL, aes(log_ratio_2_1A, log_ratio_2_1B)) +
  geom_point(data = filter(log_rep_0_22_A_B_neg_cont, subpool == 'control'),
             alpha = 0.6) +
  geom_point(data = filter(log_rep_0_22_A_B_neg_cont, 
                           subpool == 'subpool5'), 
             color = 'red', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expr. Rep. 1") +
  ylab("Variant log10 Expr. Rep. 2") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, .5)) + 
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, .5)) + 
  annotate("text", x = -1, y = 0,
           label = paste(
             'r =', round(
               cor(
                 log_rep_0_22_A_B_neg_cont$log_ratio_2_1A,
                 log_rep_0_22_A_B_neg_cont$log_ratio_2_1B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_nc_var_rep_20 <- ggplot(NULL, aes(log_ratio_20A, log_ratio_20B)) +
  geom_point(data = filter(log_rep_0_22_A_B_neg_cont, subpool == 'control'),
             alpha = 0.6) +
  geom_point(data = filter(log_rep_0_22_A_B_neg_cont, 
                           subpool == 'subpool5'), 
             color = 'red', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expr. Rep. 1") +
  ylab("Variant log10 Expr. Rep. 2") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, .5)) + 
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, .5)) + 
  annotate("text", x = -1, y = 0,
           label = paste(
             'r =', round(
               cor(
                 log_rep_0_22_A_B_neg_cont$log_ratio_20A,
                 log_rep_0_22_A_B_neg_cont$log_ratio_20B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_nc_var_rep_22 <- ggplot(NULL, aes(log_ratio_22A, log_ratio_22B)) +
  geom_point(data = filter(log_rep_0_22_A_B_neg_cont, subpool == 'control'),
             alpha = 0.6) +
  geom_point(data = filter(log_rep_0_22_A_B_neg_cont, 
                           subpool == 'subpool5'), 
             color = 'red', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expr. Rep. 1") +
  ylab("Variant log10 Expr. Rep. 2") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, .5)) + 
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, .5)) + 
  annotate("text", x = -1, y = 0,
           label = paste(
             'r =', round(
               cor(
                 log_rep_0_22_A_B_neg_cont$log_ratio_22A,
                 log_rep_0_22_A_B_neg_cont$log_ratio_22B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_nc_log_rep_grid <- plot_grid(
  p_nc_var_rep_0, p_nc_var_rep_2_5, p_nc_var_rep_2_4, p_nc_var_rep_2_3, 
  p_nc_var_rep_2_2, p_nc_var_rep_2_1, p_nc_var_rep_20, p_nc_var_rep_22,
  labels = c(
    "0 µM", "2^-5 µM", "2^-4 µM","2^-3 µM", "2^-2 µM", "2^-1 µM","2^0 µM", "2^2 µM"),
  nrow = 3, ncol = 3, align = 'hv', hjust = -3, vjust = 0.5, scale = 0.9)

save_plot('plots/p_nc_log_rep_grid.png', 
          p_nc_log_rep_grid, base_height = 10, base_width = 10)


