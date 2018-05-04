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

cbPalette7 <- c('#440154FF', '#39568CFF', '#287D8EFF', '#20A387FF', '#73D055FF',
                '#B8DE29FF', '#FDE725FF')

cbPalette3 <- c('#39568CFF', '#1F968BFF', '#73D055FF')

forskolin2 <- c('white', 'aquamarine3')

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
#mapping and replacing na with 0 reads. I mixed the transfected and sequenced 
#plasmid DNA separately from one another. Should be the same ratio so decided to 
#normalize reads in bulk, direct comparisons between subpools might not be as 
#precise because of this but it allows us to normaize reads to the background 
#sequence in SP5.

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

#sum unique barcodes and normalized bc reads per variant. Output is total 
#barcodes and sum normalized reads per variant. There is only 2 controls that 
#are not represented by a BC read in some of the samples, thus it is not present
#in joined tables later

var_sum_bc_num <- function(df1) {
  bc_count <- df1 %>%
    filter(df1$num_reads > 0) %>%
    group_by(subpool, name, most_common) %>%
    summarise(barcodes = n())
  variant_sum <- df1 %>%
    group_by(subpool, name, most_common) %>%
    count(name, wt = normalized) %>%
    rename(sum = n)
  bc_sum <- left_join(variant_sum, bc_count, 
                       by = c("name", "subpool", "most_common")) %>%
    ungroup()
  return(bc_sum)
}

variant_counts_DNA <- var_sum_bc_num(bc_join_DNA) %>%
  filter(barcodes > 2)

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

#combine DNA and RNA cumm. BC counts, only keeping instances in both sets (as 
#only 1 control drops out in some samples) and determining RNA/DNA per variant. 
#Ratio is summed normalized reads of RNA over DNA

var_expression <- function(df1, df2) {
  RNA_DNA <- inner_join(df1, df2, 
                        by = c("name", "subpool", "most_common"), 
                        suffix = c("_DNA", "_RNA")) %>%
    mutate(ratio = sum_RNA / sum_DNA) %>%
    filter(ratio > 0)
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

#After combining, rename backgrounds to simplified names, make background column 
#(excluding controls), separate out background values in each dataset and left 
#join to original dataset. Normalize expression of each variant to its 
#background in that biological replicate. Determine average normalized 
#expression across biological replicates.

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

#Normalize to background--------------------------------------------------------

#Ideally would have put this earlier as there are so many columns to modify, but
#normalizing to backgrounds excludes the controls, so I would have to repeat the
#code for the controls...

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
}

trans_back_norm_rep_0_22 <- back_norm(rep_0_22_A_B)

#Add pGl4.29 pc's into background-normalization

trans_back_norm_pc_spGl4 <- rep_0_22_A_B %>%
  filter(name == 'pGL4.29 Promega 1-63 + 1-87' | name == 'pGL4.29 Promega 1-87') %>%
  mutate(name = str_c(name, '_scramble pGL4.29 Promega 1-63 + 1-87')) %>%
  mutate(subpool = 'subpool3') %>%
  rbind(rep_0_22_A_B) %>%
  back_norm()


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
log10_trans_back_norm_rep_0_22 <- var_log10(trans_back_norm_rep_0_22)


#Make df with concentration and expression as a variable------------------------

#Make untidy data with expression and concentration as variables

var_conc_exp <- function(df) {
  df_0 <- df %>%
    mutate(ave_barcode_0 = (barcodes_RNA_0A + barcodes_RNA_0B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode_0, 
           ave_ratio_0_norm) %>%
    mutate(conc = 0) %>%
    rename(ave_ratio_norm = ave_ratio_0_norm) %>%
    rename(ave_barcode = ave_barcode_0)
  df_2_5 <- df %>%
    mutate(ave_barcode_2_5 = (barcodes_RNA_2_5A + barcodes_RNA_2_5B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode_2_5, 
           ave_ratio_2_5_norm) %>%
    mutate(conc = 2^-5) %>%
    rename(ave_ratio_norm = ave_ratio_2_5_norm) %>%
    rename(ave_barcode = ave_barcode_2_5)
  df_2_4 <- df %>%
    mutate(ave_barcode_2_4 = (barcodes_RNA_2_4A + barcodes_RNA_2_4B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode_2_4, 
           ave_ratio_2_4_norm) %>%
    mutate(conc = 2^-4) %>%
    rename(ave_ratio_norm = ave_ratio_2_4_norm) %>%
    rename(ave_barcode = ave_barcode_2_4)
  df_2_3 <- df %>%
    mutate(ave_barcode_2_3 = (barcodes_RNA_2_3A + barcodes_RNA_2_3B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode_2_3, 
           ave_ratio_2_3_norm) %>%
    mutate(conc = 2^-3) %>%
    rename(ave_ratio_norm = ave_ratio_2_3_norm) %>%
    rename(ave_barcode = ave_barcode_2_3)
  df_2_2 <- df %>%
    mutate(ave_barcode_2_2 = (barcodes_RNA_2_2A + barcodes_RNA_2_2B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode_2_2, 
           ave_ratio_2_2_norm) %>%
    mutate(conc = 2^-2) %>%
    rename(ave_ratio_norm = ave_ratio_2_2_norm) %>%
    rename(ave_barcode = ave_barcode_2_2)
  df_2_1 <- df %>%
    mutate(ave_barcode_2_1 = (barcodes_RNA_2_1A + barcodes_RNA_2_1B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode_2_1, 
           ave_ratio_2_1_norm) %>%
    mutate(conc = 2^-1) %>%
    rename(ave_ratio_norm = ave_ratio_2_1_norm) %>%
    rename(ave_barcode = ave_barcode_2_1)
  df_20 <- df %>%
    mutate(ave_barcode_20 = (barcodes_RNA_20A + barcodes_RNA_20B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode_20, 
           ave_ratio_20_norm) %>%
    mutate(conc = 2^0) %>%
    rename(ave_ratio_norm = ave_ratio_20_norm) %>%
    rename(ave_barcode = ave_barcode_20)
  df_22 <- df %>%
    mutate(ave_barcode_22 = (barcodes_RNA_22A + barcodes_RNA_22B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode_22, 
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
#consensus and flanks). Added 2 to all distances to measure to start of BS and 
#not to flank. Added 4 to all spacing but 0 to measure difference between start 
#of sites. Also took average of log2 med BC expression between biological 
#replicates for plotting

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
    mutate(dist = as.integer(dist + 2)) %>%
    mutate(spacing = 
             ifelse(spacing != as.integer(0), 
                    as.integer(spacing + 4), as.integer(spacing)))
}

s3_tidy <- subpool3(trans_back_norm_rep_0_22)
s3_untidy <- subpool3(trans_back_norm_conc)
  

#Subpool 5 contains 6 equally spaced sites spaced 13 bp apart and starting from 
#furthest to the minP. These sites are filled with sites of either the consensus
#site, a weak site or no site. Both the weak and consensus sites are flanked by 
#the same flanking sequence. 

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
    mutate(site_type = 
             ifelse(consensus == 0 & weak > 0, 
                    'weak', site_combo))
}

s5_tidy <- subpool5(trans_back_norm_rep_0_22)
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
  

#Plot subpool expression features-----------------------------------------------

#Subpool 3

p_subpool3_spa_back_norm <- ggplot(subpool3_log2_norm, aes(x = dist)) + 
  geom_point(aes(y = ave_ratio_0_norm), alpha = 0.5, size = 1.5, 
             color = '#440154FF') +
  geom_point(aes(y = ave_ratio_2_5_norm), alpha = 0.5, size = 1.5, 
             color = '#482677FF') +
  geom_point(aes(y = ave_ratio_2_4_norm), alpha = 0.5, size = 1.5, 
             color = '#39568CFF') +
  geom_point(aes(y = ave_ratio_2_3_norm), alpha = 0.5, size = 1.5, 
             color = '#2D708EFF') +
  geom_point(aes(y = ave_ratio_2_2_norm), alpha = 0.5, size = 1.5, 
             color = '#1F968BFF') +
  geom_point(aes(y = ave_ratio_2_1_norm), alpha = 0.5, size = 1.5, 
             color = '#29AF7FFF') +
  geom_point(aes(y = ave_ratio_20_norm), alpha = 0.5, size = 1.5, 
             color = '#73D055FF') +
  geom_point(aes(y = ave_ratio_22_norm), alpha = 0.5, size = 1.5, 
             color = '#B8DE29FF') +
  facet_grid(spacing ~ background) + 
  ylab('Log2 normalized average sum BC expression') + 
  panel_border() +
  background_grid(major = 'xy', minor = 'none') +
  scale_x_continuous(
    "Distance from First Site to Proximal Promoter End (bp)", 
    breaks = seq(from = 0, to = 150, by = 10))

save_plot('plots/p_subpool3_spa_back_norm.png', p_subpool3_spa_back_norm, 
          base_width = 46, base_height = 17, scale = 0.35)

p_subpool3_spa_back_norm_smooth <- ggplot(subpool3_log2_norm, aes(x = dist)) + 
  geom_point(aes(y = ave_ratio_0_norm), alpha = 0, size = 1.5, 
             color = '#440154FF') +
  geom_smooth(aes(y = ave_ratio_0_norm), 
              span = 0.1, size = 0.7, color = '#440154FF', se = FALSE) +
  geom_point(aes(y = ave_ratio_2_5_norm), alpha = 0, size = 1.5, 
             color = '#482677FF') +
  geom_smooth(aes(y = ave_ratio_2_5_norm), 
              span = 0.1, size = 0.7, color = '#482677FF', se = FALSE) +
  geom_point(aes(y = ave_ratio_2_4_norm), alpha = 0, size = 1.5, 
             color = '#39568CFF') +
  geom_smooth(aes(y = ave_ratio_2_4_norm), 
              span = 0.1, size = 0.7, color = '#39568CFF', se = FALSE) +
  geom_point(aes(y = ave_ratio_2_3_norm), alpha = 0, size = 1.5, 
             color = '#2D708EFF') +
  geom_smooth(aes(y = ave_ratio_2_3_norm), 
              span = 0.1, size = 0.7, color = '#2D708EFF', se = FALSE) +
  geom_point(aes(y = ave_ratio_2_2_norm), alpha = 0, size = 1.5, 
             color = '#1F968BFF') +
  geom_smooth(aes(y = ave_ratio_2_2_norm), 
              span = 0.1, size = 0.7, color = '#1F968BFF', se = FALSE) +
  geom_point(aes(y = ave_ratio_2_1_norm), alpha = 0, size = 1.5, 
             color = '#29AF7FFF') +
  geom_smooth(aes(y = ave_ratio_2_1_norm), 
              span = 0.1, size = 0.7, color = '#29AF7FFF', se = FALSE) +
  geom_point(aes(y = ave_ratio_20_norm), alpha = 0, size = 1.5, 
             color = '#73D055FF') +
  geom_smooth(aes(y = ave_ratio_20_norm), 
              span = 0.1, size = 0.7, color = '#73D055FF', se = FALSE) +
  geom_point(aes(y = ave_ratio_22_norm), alpha = 0, size = 1.5, 
             color = '#B8DE29FF') +
  geom_smooth(aes(y = ave_ratio_22_norm), 
              span = 0.1, size = 0.7, color = '#B8DE29FF', se = FALSE) +
  facet_grid(spacing ~ background) + 
  ylab('Log2 normalized average sum BC expression') + 
  panel_border() +
  background_grid(major = 'xy', minor = 'none') +
  scale_x_continuous(
    "Distance from First Site to Proximal Promoter End (bp)", 
    breaks = seq(from = 0, to = 150, by = 10))

save_plot('plots/p_subpool3_spa_back_norm_smooth.png', 
          p_subpool3_spa_back_norm_smooth, 
          base_width = 46, base_height = 17, scale = 0.35)

#Subpool 5

p_s5_consnum_0_4 <- s5_untidy %>%
  filter(weak == 0) %>%
  filter(conc == 0 | conc == 4) %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  group_by(background, consensus, conc) %>%
  mutate(med_back_cons_conc = median(ave_ratio_norm)) %>%
  ungroup() %>%
  ggplot(aes(as.factor(consensus), ave_ratio_norm, fill = as.factor(conc))) +
  geom_boxplot(outlier.size = 0.7, size = 0.3, 
               outlier.alpha = 0.5, position = position_dodge(0.75),
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

p_s5_weaknum_0_4 <- s5_untidy %>%
  filter(consensus == 0) %>%
  filter(conc == 0 | conc == 4) %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  group_by(background, consensus, conc) %>%
  mutate(med_back_cons_conc = median(ave_ratio_norm)) %>%
  ungroup() %>%
  ggplot(aes(as.factor(weak), ave_ratio_norm, fill = as.factor(conc))) +
  geom_boxplot(outlier.size = 0.7, size = 0.3, 
               outlier.alpha = 0.5, position = position_dodge(0.75),
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

#in order to facet by site_type, need to duplicate the 6 nosite baseline in each
#category

s5_dup_6no <- function(df) {
  no_cons <- df %>%
    filter(total_sites == 0) %>%
    mutate(site_type = 'consensus')
  no_weak <- df %>%
    filter(total_sites == 0) %>%
    mutate(site_type = 'weak')
  no_join_rm_mixed <- rbind(no_cons, no_weak, df) %>%
    filter(site_type != 'mixed')
  return(no_join_rm_mixed)
}

p_s5_sitenum_ind <- s5_tidy %>%
  s5_dup_6no() %>%
  mutate(background = factor(background, 
                             levels = c('v chr9', 's pGl4', 'v chr5'))) %>%
  ggplot(aes(as.factor(total_sites), induction)) +
  geom_boxplot(aes(color = background), outlier.size = 0.7, size = 0.3, 
               outlier.alpha = 0.5, position = position_dodge(1),
               show.legend = TRUE) +
  scale_y_log10() +
  facet_grid(site_type ~ .) +
  panel_border() +
  annotation_logticks(sides = 'l') +
  scale_color_manual(values = cbPalette3) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  background_grid(major = 'y', minor = 'none') + 
  geom_vline(xintercept = c(1.5:6.5), alpha = 0.5) +
  ylab('Induction (a.u.)') +
  xlab('Number of total sites')


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

#plot replicates for summed variant expression

#plot all concentrations faceted for main figure

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


#Comparison to integrated-------------------------------------------------------------------

#Join datasets, log10 transform 

int_back_norm_rep_1_2 <- read_tsv('../20171129_intLib/int_back_norm_rep_1_2.txt')

int_rep_1_2 <- read_tsv('../20171129_intLib/rep_1_2.txt')

int_trans <- left_join(int_rep_1_2, rep_0_22_A_B, 
                       by = c('subpool', 'name', 'most_common'))

log10_int_trans <- var_log10(int_trans) %>%
  mutate(ave_med_ratio = (med_ratio_br1 + med_ratio_br2)/2) %>%
  mutate(ave_ratio_22 = (ratio_22A + ratio_22B)/2)

log10_int_trans_sp5 <- log10_int_trans %>%
  filter(subpool == 'subpool5')

p_var_log10_int_trans_4_sp5 <- ggplot(log10_int_trans_sp5, 
                                  aes(ave_med_ratio, ave_ratio_22)) +
  geom_point(alpha = 0.2) +
  geom_point(data = filter(log10_int_trans_sp5, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Integrated\naverage log10\nmedian expression") +
  ylab("Transient\naverage log10\nsum expression") +
  annotate("text", x = -1.5, y = 1,
           label = paste('r =', 
                         round(cor(log10_int_trans_sp5$ave_med_ratio, 
                                   log10_int_trans_sp5$ave_ratio_22,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))

save_plot('plots/var_log10_int_trans_4_sp5.png', p_var_log10_int_trans_4_sp5)


#Fitting a simple linear model integrated_exp ~ transient_exp(4 µM)

log10_int_trans_22 <- log10_int_trans %>%
  select(subpool, name, background, barcodes_br1, barcodes_br2, ave_med_ratio_norm, 
         barcodes_RNA_22A, barcodes_RNA_22B, ave_ratio_22_norm) %>%
  mutate(integrated = (barcodes_br1 + barcodes_br2)/2) %>%
  mutate(transient = (barcodes_RNA_22A + barcodes_RNA_22B)/2) %>%
  select(-barcodes_br1, -barcodes_br2, -barcodes_RNA_22A, -barcodes_RNA_22B)

log10_MPRA <- log10_int_trans %>%
  select(subpool, name, background, barcodes_br1, barcodes_br2, ave_med_ratio_norm, 
         barcodes_RNA_22A, barcodes_RNA_22B, ave_ratio_22_norm) %>%
  mutate(integrated = (barcodes_br1 + barcodes_br2)/2) %>%
  mutate(transient = (barcodes_RNA_22A + barcodes_RNA_22B)/2) %>%
  select(-barcodes_br1, -barcodes_br2, -barcodes_RNA_22A, -barcodes_RNA_22B) %>%
  gather(integrated, transient, key = 'MPRA', value = 'barcodes') %>%
  rename(integrated = ave_med_ratio_norm) %>%
  rename(transient = ave_ratio_22_norm) %>%
  gather(integrated, transient, key = 'MPRA2', value = 'ave_ratio_norm') %>%
  filter((MPRA == 'integrated' & MPRA2 == 'integrated') | (MPRA == 'transient' & MPRA2 == 'transient')) %>%
  select(-MPRA2)
  

lm_trans_int <- function(df) {
  model <- lm(ave_med_ratio_norm ~ ave_ratio_22_norm, data = df)
}

pre_res_trans_int <- function(df1, x) {
  df2 <- df1 %>%
    add_predictions(x)
  df3 <- df2 %>%
    add_residuals(x)
  return(df3)
  print('processed pre_res_trans_int(df1, df2) in order of (data, model)')
}

#testing a simple linear model

lm_all <- lm_trans_int(log10_int_trans_22)
p_r_all <- pre_res_trans_int(log10_int_trans_22, lm_all)

p_lm_int_trans <- ggplot(p_r_all, aes(x = ave_ratio_22_norm)) +
  facet_grid(subpool ~ background) +
  geom_point(aes(y = ave_med_ratio_norm, color = ave_int_barcode), alpha = 0.5) +
  scale_color_viridis() +
  geom_line(aes(y = pred), color = 'red') +
  annotation_logticks(scaled = TRUE) +
  ylab("Ave log10 variant median\nnorm. RNA/DNA integrated") +
  xlab("Ave log10 variant norm.\nsum RNA/DNA transient") +
  panel_border()

save_plot('plots/p_lm_int_trans.png', p_lm_int_trans, base_width = 10, base_height = 7)

p_res_int_trans <- ggplot(p_r_all, aes(ave_ratio_22_norm, resid)) +
  facet_grid(subpool ~ background) +
  geom_point(aes(color = ave_int_barcode), alpha = 0.5) +
  scale_color_viridis() +
  geom_ref_line(h = 0, colour = 'black', size = 1) +
  annotation_logticks(scaled = TRUE, sides = 'b') +
  ylab("Residuals") +
  xlab("Ave log10 variant norm.\nsum RNA/DNA transient") +
  panel_border()

save_plot('plots/p_res_int_trans.png', p_res_int_trans, base_width = 13, base_height = 7)

p_res_int_trans_sp <- ggplot(p_r_all, aes(x = subpool, y = resid)) +
  facet_grid(. ~ background) +
  geom_boxplot(alpha = 0.5) +
  geom_ref_line(h = 0, colour = 'black', size = 1) +
  ylab("Residuals") +
  xlab("Subpool") +
  panel_border()

#Analyze differences between transient and integrated using features within each subpool

subpool3_int_trans <- 
  filter(log10_MPRA, subpool == "subpool3") %>%
  ungroup () %>%
  select(-subpool) %>%
  mutate(name = gsub('2BS ', '', name), 
         name = gsub(' bp spacing ', '_', name)) %>%
  separate(name, 
           into = c("subpool", "spacing", "fluff2", "fluff3", "dist", "fluff4"),
           sep = "_", convert = TRUE
  ) %>%
  select(-subpool, -fluff2, -fluff3, -fluff4) %>%
  mutate(dist = as.integer(dist + 2)) %>%
  mutate(spacing = 
           ifelse(spacing != as.integer(0), 
                  as.integer(spacing + 4), as.integer(spacing))) 

subpool5_int_trans <- 
  filter(log10_MPRA, subpool == "subpool5") %>%
  ungroup () %>%
  select(-subpool) %>%
  mutate(name = gsub('no_site', 'nosite', name)) %>%
  mutate(
    name = gsub('_v chr9', '', name),
    name = gsub('_v chr5', '', name),
    name = gsub('_s pGl4', '', name)
    ) %>%
  mutate(design = name) %>%
  separate(name, into = c(
    "subpool", "site1", "site2", "site3", "site4", "site5", "site6"), sep = "_"
    ) %>%
  select(-subpool) %>%
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
  mutate(total_sites = consensus + weak)

subpool5_int_trans_cons <- subpool5_int_trans %>%
  filter(weak == 0) %>%
  mutate(site1 = as.integer(str_detect(site1, 'consensus'))) %>%
  mutate(site2 = as.integer(str_detect(site2, 'consensus'))) %>%
  mutate(site3 = as.integer(str_detect(site3, 'consensus'))) %>%
  mutate(site4 = as.integer(str_detect(site4, 'consensus'))) %>%
  mutate(site5 = as.integer(str_detect(site5, 'consensus'))) %>%
  mutate(site6 = as.integer(str_detect(site6, 'consensus')))
  
#If want to have site as one variable figure out a better way than 
  #gather(site1, site2, site3, site4, site5, site6, key = 'site', value = 'filled')
  
#plot subpool features vs. residuals between transient to integrated model

p_subpool5_int_trans_consw <- ggplot(subpool5_int_trans, 
                                    aes(x = as.factor(consensus), y = ave_ratio_norm,
                                               color = MPRA)) +
  geom_boxplot() +
  facet_grid(. ~ background) +
  xlab('Number of consensus binding sites') +
  scale_y_continuous('Average\nnormalized expression') +
  panel_border()

save_plot('plots/p_subpool5_int_trans_consw.png', p_subpool5_int_trans_consw,
          base_height = 4, base_width = 8)

p_subpool5_int_trans_cons <- ggplot(filter(subpool5_int_trans, weak == 0), 
                                     aes(x = as.factor(consensus), y = ave_ratio_norm,
                                         color = MPRA)) +
  geom_boxplot() +
  facet_grid(. ~ background) +
  xlab('Number of consensus binding sites') +
  scale_y_continuous('Average\nnormalized expression') +
  panel_border()

save_plot('plots/p_subpool5_int_trans_cons.png', p_subpool5_int_trans_cons,
          base_height = 4, base_width = 8)

ggplot(NULL, aes(x = '', y = ave_ratio_norm, color = MPRA)) +
  facet_grid(background ~ .) +
  geom_point(data = filter(subpool5_int_trans_cons, site1 == 1 & total_sites == 1), 
             aes(x = 'site1'), size = 2, show.legend = FALSE) +
  geom_point(data = filter(subpool5_int_trans_cons, site2 == 1 & total_sites == 1), 
             aes(x = 'site2'), size = 2, show.legend = FALSE) +
  geom_point(data = filter(subpool5_int_trans_cons, site3 == 1 & total_sites == 1), 
             aes(x = 'site3'), size = 2, show.legend = FALSE) +
  geom_point(data = filter(subpool5_int_trans_cons, site4 == 1 & total_sites == 1), 
             aes(x = 'site4'), size = 2, show.legend = FALSE) +
  geom_point(data = filter(subpool5_int_trans_cons, site5 == 1 & total_sites == 1), 
             aes(x = 'site5'), size = 2, show.legend = FALSE) +
  geom_point(data = filter(subpool5_int_trans_cons, site6 == 1 & total_sites == 1), 
             aes(x = 'site6'), size = 2, show.legend = FALSE) +
  panel_border() +
  annotation_logticks(sides = 'l') +
  ylab('log10 average normalized ratio') +
  xlab('Binding site architecture')

#Plot comparisons between trans concentrations and integrated

p_var_log10_int_trans_0 <- ggplot(NULL, aes(ave_med_ratio_norm, ave_ratio_0_norm)) +
  geom_point(data = filter(log10_int_trans, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(log10_int_trans, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = log10_int_trans, alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log10_int_trans, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.5) +
  geom_point(data = filter(log10_int_trans, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Ave log10 variant median\nnorm. RNA/DNA integrated") +
  ylab("Ave log10 variant norm.\nsum RNA/DNA transient") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
                     limits = c(-2, 8.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
                     limits = c(-2, 8.5)) + 
  annotate("text", x = 0, y = 7,
           label = paste('r =', round(
             cor(
               log10_int_trans$ave_med_ratio_norm, 
               log10_int_trans$ave_ratio_0_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) +
  annotate("text", x = 0, y = 6, color = '#2D708EFF',
           label = paste('r =', round(
             cor(
               filter(log10_int_trans, subpool == 'subpool3')$ave_med_ratio_norm, 
               filter(log10_int_trans, subpool == 'subpool3')$ave_ratio_0_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) + 
  annotate("text", x = 0, y = 5, color = '#482677FF',
           label = paste('r =', round(
             cor(
               filter(log10_int_trans, subpool == 'subpool5')$ave_med_ratio_norm, 
               filter(log10_int_trans, subpool == 'subpool5')$ave_ratio_0_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  )

p_var_log10_int_trans_2_5 <- ggplot(NULL, aes(ave_med_ratio_norm, ave_ratio_2_5_norm)) +
  geom_point(data = filter(log10_int_trans, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(log10_int_trans, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = log10_int_trans, alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log10_int_trans, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.5) +
  geom_point(data = filter(log10_int_trans, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Ave log10 variant median\nnorm. RNA/DNA integrated") +
  ylab("Ave log10 variant norm.\nsum RNA/DNA transient") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
                     limits = c(-2, 8.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
                     limits = c(-2, 8.5)) + 
  annotate("text", x = 0, y = 7,
           label = paste('r =', round(
             cor(
               log10_int_trans$ave_med_ratio_norm, 
               log10_int_trans$ave_ratio_2_5_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) +
  annotate("text", x = 0, y = 6, color = '#2D708EFF',
           label = paste('r =', round(
             cor(
               filter(log10_int_trans, subpool == 'subpool3')$ave_med_ratio_norm, 
               filter(log10_int_trans, subpool == 'subpool3')$ave_ratio_2_5_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) + 
  annotate("text", x = 0, y = 5, color = '#482677FF',
           label = paste('r =', round(
             cor(
               filter(log10_int_trans, subpool == 'subpool5')$ave_med_ratio_norm, 
               filter(log10_int_trans, subpool == 'subpool5')$ave_ratio_2_5_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  )

p_var_log10_int_trans_2_4 <- ggplot(NULL, aes(ave_med_ratio_norm, ave_ratio_2_4_norm)) +
  geom_point(data = filter(log10_int_trans, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(log10_int_trans, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = log10_int_trans, alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log10_int_trans, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.5) +
  geom_point(data = filter(log10_int_trans, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Ave log10 variant median\nnorm. RNA/DNA integrated") +
  ylab("Ave log10 variant norm.\nsum RNA/DNA transient") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
                     limits = c(-2, 8.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
                     limits = c(-2, 8.5)) + 
  annotate("text", x = 0, y = 7,
           label = paste('r =', round(
             cor(
               log10_int_trans$ave_med_ratio_norm, 
               log10_int_trans$ave_ratio_2_4_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) +
  annotate("text", x = 0, y = 6, color = '#2D708EFF',
           label = paste('r =', round(
             cor(
               filter(log10_int_trans, subpool == 'subpool3')$ave_med_ratio_norm, 
               filter(log10_int_trans, subpool == 'subpool3')$ave_ratio_2_4_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) + 
  annotate("text", x = 0, y = 5, color = '#482677FF',
           label = paste('r =', round(
             cor(
               filter(log10_int_trans, subpool == 'subpool5')$ave_med_ratio_norm, 
               filter(log10_int_trans, subpool == 'subpool5')$ave_ratio_2_4_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  )

p_var_log10_int_trans_2_3 <- ggplot(NULL, aes(ave_med_ratio_norm, ave_ratio_2_3_norm)) +
  geom_point(data = filter(log10_int_trans, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(log10_int_trans, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = log10_int_trans, alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log10_int_trans, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.5) +
  geom_point(data = filter(log10_int_trans, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Ave log10 variant median\nnorm. RNA/DNA integrated") +
  ylab("Ave log10 variant norm.\nsum RNA/DNA transient") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
                     limits = c(-2, 8.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
                     limits = c(-2, 8.5)) + 
  annotate("text", x = 0, y = 7,
           label = paste('r =', round(
             cor(
               log10_int_trans$ave_med_ratio_norm, 
               log10_int_trans$ave_ratio_2_3_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) +
  annotate("text", x = 0, y = 6, color = '#2D708EFF',
           label = paste('r =', round(
             cor(
               filter(log10_int_trans, subpool == 'subpool3')$ave_med_ratio_norm, 
               filter(log10_int_trans, subpool == 'subpool3')$ave_ratio_2_3_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) + 
  annotate("text", x = 0, y = 5, color = '#482677FF',
           label = paste('r =', round(
             cor(
               filter(log10_int_trans, subpool == 'subpool5')$ave_med_ratio_norm, 
               filter(log10_int_trans, subpool == 'subpool5')$ave_ratio_2_3_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  )

p_var_log10_int_trans_2_2 <- ggplot(NULL, aes(ave_med_ratio_norm, ave_ratio_2_2_norm)) +
  geom_point(data = filter(log10_int_trans, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(log10_int_trans, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = log10_int_trans, alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log10_int_trans, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.5) +
  geom_point(data = filter(log10_int_trans, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Ave log10 variant median\nnorm. RNA/DNA integrated") +
  ylab("Ave log10 variant norm.\nsum RNA/DNA transient") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
                     limits = c(-2, 8.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
                     limits = c(-2, 8.5)) + 
  annotate("text", x = 0, y = 7,
           label = paste('r =', round(
             cor(
               log10_int_trans$ave_med_ratio_norm, 
               log10_int_trans$ave_ratio_2_2_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) +
  annotate("text", x = 0, y = 6, color = '#2D708EFF',
           label = paste('r =', round(
             cor(
               filter(log10_int_trans, subpool == 'subpool3')$ave_med_ratio_norm, 
               filter(log10_int_trans, subpool == 'subpool3')$ave_ratio_2_2_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) + 
  annotate("text", x = 0, y = 5, color = '#482677FF',
           label = paste('r =', round(
             cor(
               filter(log10_int_trans, subpool == 'subpool5')$ave_med_ratio_norm, 
               filter(log10_int_trans, subpool == 'subpool5')$ave_ratio_2_2_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  )

p_var_log10_int_trans_2_1 <- ggplot(NULL, aes(ave_med_ratio_norm, ave_ratio_2_1_norm)) +
  geom_point(data = filter(log10_int_trans, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(log10_int_trans, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = log10_int_trans, alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log10_int_trans, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.5) +
  geom_point(data = filter(log10_int_trans, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Ave log10 variant median\nnorm. RNA/DNA integrated") +
  ylab("Ave log10 variant norm.\nsum RNA/DNA transient") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
                     limits = c(-2, 8.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
                     limits = c(-2, 8.5)) + 
  annotate("text", x = 0, y = 7,
           label = paste('r =', round(
             cor(
               log10_int_trans$ave_med_ratio_norm, 
               log10_int_trans$ave_ratio_2_1_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) +
  annotate("text", x = 0, y = 6, color = '#2D708EFF',
           label = paste('r =', round(
             cor(
               filter(log10_int_trans, subpool == 'subpool3')$ave_med_ratio_norm, 
               filter(log10_int_trans, subpool == 'subpool3')$ave_ratio_2_1_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) + 
  annotate("text", x = 0, y = 5, color = '#482677FF',
           label = paste('r =', round(
             cor(
               filter(log10_int_trans, subpool == 'subpool5')$ave_med_ratio_norm, 
               filter(log10_int_trans, subpool == 'subpool5')$ave_ratio_2_1_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  )

p_var_log10_int_trans_20 <- ggplot(NULL, aes(ave_med_ratio_norm, ave_ratio_20_norm)) +
  geom_point(data = filter(log10_int_trans, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(log10_int_trans, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = log10_int_trans, alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log10_int_trans, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.5) +
  geom_point(data = filter(log10_int_trans, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Ave log10 variant median\nnorm. RNA/DNA integrated") +
  ylab("Ave log10 variant norm.\nsum RNA/DNA transient") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
                     limits = c(-2, 8.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
                     limits = c(-2, 8.5)) + 
  annotate("text", x = 0, y = 7,
           label = paste('r =', round(
             cor(
               log10_int_trans$ave_med_ratio_norm, 
               log10_int_trans$ave_ratio_20_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) +
  annotate("text", x = 0, y = 6, color = '#2D708EFF',
           label = paste('r =', round(
             cor(
               filter(log10_int_trans, subpool == 'subpool3')$ave_med_ratio_norm, 
               filter(log10_int_trans, subpool == 'subpool3')$ave_ratio_20_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) + 
  annotate("text", x = 0, y = 5, color = '#482677FF',
           label = paste('r =', round(
             cor(
               filter(log10_int_trans, subpool == 'subpool5')$ave_med_ratio_norm, 
               filter(log10_int_trans, subpool == 'subpool5')$ave_ratio_20_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  )

p_var_log10_int_trans_22 <- ggplot(NULL, aes(ave_med_ratio_norm, ave_ratio_22_norm)) +
  geom_point(data = filter(log10_int_trans, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(log10_int_trans, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = log10_int_trans, alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log10_int_trans, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.5) +
  geom_point(data = filter(log10_int_trans, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Ave log10 variant median\nnorm. RNA/DNA integrated") +
  ylab("Ave log10 variant norm.\nsum RNA/DNA transient") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
                     limits = c(-2, 8.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
                     limits = c(-2, 8.5)) + 
  annotate("text", x = 0, y = 7,
           label = paste('r =', round(
             cor(
               log10_int_trans$ave_med_ratio_norm, 
               log10_int_trans$ave_ratio_22_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) +
  annotate("text", x = 0, y = 6, color = '#2D708EFF',
           label = paste('r =', round(
             cor(
               filter(log10_int_trans, subpool == 'subpool3')$ave_med_ratio_norm, 
               filter(log10_int_trans, subpool == 'subpool3')$ave_ratio_22_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  ) + 
  annotate("text", x = 0, y = 5, color = '#482677FF',
           label = paste('r =', round(
             cor(
               filter(log10_int_trans, subpool == 'subpool5')$ave_med_ratio_norm, 
               filter(log10_int_trans, subpool == 'subpool5')$ave_ratio_22_norm,
               use = "pairwise.complete.obs", method = "pearson"
             ), 2
           )
           )
  )

p_var_log10_int_trans_grid <- plot_grid(
  p_var_log10_int_trans_0, p_var_log10_int_trans_2_5, p_var_log10_int_trans_2_4, 
  p_var_log10_int_trans_2_3, p_var_log10_int_trans_2_2, p_var_log10_int_trans_2_1, 
  p_var_log10_int_trans_20, p_var_log10_int_trans_22, 
  labels = c(
    "    0 µM", "2^-5 µM", "2^-4 µM","2^-3 µM", "2^-2 µM", "2^-1 µM"," 2^0 µM", " 2^2 µM"),
  nrow = 3, ncol = 3, align = 'hv', hjust = -2.5, vjust = 0.5, scale = 0.9)

save_plot('plots/p_var_log10_int_trans_grid.png', 
          p_var_log10_int_trans_grid, base_height = 12, base_width = 12)

save_plot('plots/p_var_log10_int_trans_22.png', p_var_log10_int_trans_22)


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

#Performing this on averaged data first, can go back and use both replicates to 
#fit later

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

#plots of titration overall with 1 control (duplicated pGl4, original does not
#have high expression)

p_titr_pc_back <- ggplot(trans_back_0_norm_conc, aes(conc, ave_ratio_norm)) +
  geom_line(aes(group = name), alpha = 0.1) +
  geom_point(data = filter(trans_back_0_norm_conc, 
                           startsWith(name, 
                                      'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')),
             color = 'darkgoldenrod1', shape = 19, stroke = 1.25) +
  geom_line(data = filter(trans_back_0_norm_conc, 
                          startsWith(name, 
                                     'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')),
            color = 'darkgoldenrod1', size = 1.25) +
  geom_point(data = filter(trans_back_0_norm_conc, 
                           startsWith(name, 
                                      'pGL4.29 Promega 1-63 + 1-87')),
             color = 'firebrick2', shape = 19, stroke = 1.25) +
  geom_line(data = filter(trans_back_0_norm_conc, 
                          startsWith(name, 
                                     'pGL4.29 Promega 1-63 + 1-87')),
            color = 'firebrick2', size = 1.25) +
  ylab('Average normalized\nexpression (a.u.)') +
  xlab('Forskolin concentration (µM)')

save_plot('plots/p_titr_pc_back.pdf', p_titr_pc_back, scale = 1.3, 
          base_width = 3.75, base_height = 2.5)


#michaelis model

m_m_model <- function(df) {
  conc_half_max_init <- (2^-3)
  max_ave_ratio_norm_init <- 0.5
  m_m_nls <- nls(
    ave_ratio_norm ~ (max_ave_ratio_norm * conc)/(conc_half_max + conc),
    data = df, start = c(conc_half_max = conc_half_max_init, 
                         max_ave_ratio_norm = max_ave_ratio_norm_init))
  return(m_m_nls)
}

#Trying to avoid the single gradient error message with nlslm........

library(minpack.lm)

m_m_model_nlslm <- function(df) {
  m_m_nlslm <- nlsLM(
    ave_ratio_norm ~ (max_ave_ratio_norm * conc)/(conc_half_max + conc),
    data = df, 
    start = list(conc_half_max = (2^-3), max_ave_ratio_norm = 2))
  return(m_m_nlslm)
}


#Fitting single variant

trans_back_0_norm_conc_sample1 <- trans_back_0_norm_conc %>%
  filter(subpool == 'subpool5') %>%
  select(-ave_barcode) %>%
  group_by(subpool, name, most_common, background) %>%
  nest() %>%
  sample_n(1) %>%
  unnest()

ggplot(trans_back_0_norm_conc_sample1, 
       aes(conc, ave_ratio_norm, color = name)) +
  geom_point(show.legend = FALSE) +
  geom_line() +
  xlab('Forskolin µM') +
  ylab('Average background-normalized\nsum RNA/DNA')

m_m_fit <- m_m_model(trans_back_0_norm_conc_sample1)
summary(m_m_fit)

m_m_fit_nlslm <- m_m_model_nlslm(trans_back_0_norm_conc_sample1)
summary(m_m_fit_nlslm)

pred_resid <- function(df1, x) {
  df2 <- df1 %>%
    add_predictions(x)
  df3 <- df2 %>%
    add_residuals(x)
  return(df3)
  print('processed pre_res_trans_int(df1, df2) in order of (data, model)')
}

m_m_p_r <- pred_resid(trans_back_0_norm_conc_sample1, m_m_fit)

ggplot(m_m_p_r, aes(x = conc)) +
  geom_point(aes(y = ave_ratio_norm), color = 'black') +
  geom_line(aes(y = ave_ratio_norm), color = 'black') +
  geom_point(aes(y = pred), color = 'red') +
  geom_line(aes(y = pred), color = 'red', linetype = 2) +
  xlab('Forskolin µM') +
  ylab('Average background-normalized\nsum RNA/DNA')


#Fitting nested data

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
                                      'background', 'conc', 'ave_ratio_norm'))
  return(data_pred_resid)
}

#Select higher expressing variants to fit curves to, had hard time fitting
#nested data to 2500-4172 variants, but these drove lower than 0.75 expression 
#above 0, would not expect a well-fitting curve for these lower-expressing 
#variants

trans_back_0_norm_conc_nest <- trans_back_0_norm_conc %>%
  select(-ave_barcode) %>%
  filter(!grepl('^subpool5_no_site_no_site_no_site_no_site_no_site_no_site', name)) %>%
  filter(subpool == 'subpool5') %>%
  arrange(desc(ave_ratio_norm)) %>%
  group_by(subpool, name, most_common, background) %>%
  nest() %>%
  slice(1:1500)

ggplot(unnest(trans_back_0_norm_conc_nest), 
       aes(conc, ave_ratio_norm, color = name)) +
  geom_point(show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  xlab('Forskolin µM') +
  scale_y_continuous() +
  ylab('Average background-normalized\nsum RNA/DNA')

m_m_nest_fit_nlslm <- trans_back_0_norm_conc_nest %>%
  mutate(m_m_fit = map(trans_back_0_norm_conc_nest$data, m_m_model_nlslm))

m_m_coef_nlslm <- m_m_nest_coef(m_m_nest_fit_nlslm)

m_m_p_r <- m_m_nest_pred_resid(m_m_nest_fit_nlslm)

test <- m_m_p_r %>%
  group_by(subpool, name, most_common, background) %>%
  nest() %>%
  sample_n(50) %>%
  unnest()

ggplot(test, aes(x = conc)) +
  facet_wrap(~ name) +
  panel_border() +
  geom_point(aes(y = ave_ratio_norm), color = 'black') +
  geom_line(aes(y = ave_ratio_norm), color = 'black') +
  geom_point(aes(y = pred), color = 'red') +
  geom_line(aes(y = pred), color = 'red', linetype = 2) +
  xlab('Forskolin µM') +
  scale_y_continuous() +
  ylab('Average background-normalized\nsum RNA/DNA')

EC50 <- m_m_coef_nlslm %>%
  filter(term == 'conc_half_max') %>%
  select(-term) %>%
  rename(EC50 = estimate) %>%
  mutate(rel_std_error = std.error/EC50) %>%
  filter(rel_std_error <= 0.25)

s5_sep <- function(df1) {
  df1 %>%
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
  mutate(site_type = 
           ifelse(consensus == 0 & weak > 0, 
                  'weak', site_combo))
}

s5_EC50 <- s5_sep(EC50)

ggplot(filter(s5_EC50, weak == 0), aes(factor(consensus), EC50)) +
  geom_boxplot(aes(color = background))

s5_hill_frac <- EC50 %>%
  left_join(trans_back_0_norm_conc,
            by = c('name', 'subpool', 'most_common', 'background')) %>%
  mutate(frac = conc/(conc + EC50)) %>%
  filter(conc > 0) %>%
  mutate(conc_over_EC50 = conc/EC50) %>%
  s5_sep()

ggplot(filter(s5_hill_frac, weak == 0), aes(log2(conc), log2(frac/(1-frac)))) +
  facet_grid(~ consensus) +
  geom_line(aes(color = most_common), show.legend = FALSE) + 
  background_grid(major = 'xy')


#Hill plots---------------------------------------------------------------------

#Sample data and plot to visualize trend

trans_back_norm_conc_50 <- trans_back_norm_conc_log2 %>%
  group_by(subpool, name, most_common, background) %>%
  nest() %>%
  sample_n(50) %>%
  unnest()

ggplot(trans_back_norm_conc_50, aes(conc, ave_ratio_norm, color = name)) +
  geom_point(show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ subpool) +
  scale_x_continuous(breaks = c(-7:2)) +
  xlab('log2 forskolin µM') +
  ylab('Average background-normalized\nsum RNA/DNA') +
  panel_border() +
  annotation_logticks(sides = 'b') +
  background_grid(major = 'xy', minor = 'none')


#Try fitting hill/log curves

#Pull out 1 variant to test

trans_back_norm_conc_sp5 <- trans_back_norm_conc_log2 %>%
  filter(subpool == 'subpool5') %>%
  group_by(subpool, name, most_common, background) %>%
  nest() %>%
  sample_n(1)

ggplot(unnest(trans_back_norm_conc_sp5), 
       aes(conc, ave_ratio_norm, color = name)) +
  geom_point(show.legend = FALSE) +
  geom_line() +
  scale_x_continuous(breaks = c(-7:2)) +
  xlab('log2 forskolin µM') +
  ylab('Average background-normalized\nsum RNA/DNA')

#Log curve on log2 conc

trans_back_norm_conc8_sp5 <- trans_back_norm_conc_log2 %>%
  unnest() %>%
  mutate(conc = log2(conc) + 8) %>%
  group_by(subpool, name, most_common, background) %>%
  nest()

log_curve_model <- function(df) {
  n_init <- 1
  conc_half_max_init <- -4
  max_ave_ratio_norm_init <- 10
  log_curve_nls <- nls(
    ave_ratio_25_norm ~ max_ave_ratio_norm/(1 + exp(-n * (conc - conc_half_max))),
    data = df, start = c(n = n_init, conc_half_max = conc_half_max_init, 
                         max_ave_ratio_norm = max_ave_ratio_norm_init))
  return(log_curve_nls)
}

pred_resid <- function(df1, x) {
  df2 <- df1 %>%
    add_predictions(x)
  df3 <- df2 %>%
    add_residuals(x)
  return(df3)
  print('processed pre_res_trans_int(df1, df2) in order of (data, model)')
}

log_curve_fit <- log_curve_model(bin_site_s5)
summary(log_curve_fit)
log_curve_p_r <- pred_resid(bin_site_s5, log_curve_fit)


#Try making fraction occupied hill plots

#Let's test this on subpool5 and facet by total # sites

sp5_conc_exp <- function(df) {
  df_0 <- df %>%
    mutate(ave_barcode_0 = (barcodes_RNA_0A + barcodes_RNA_0B)/2) %>%
    select(consensus, weak, nosite, most_common, background, ave_barcode_0, 
           ave_ratio_0_norm) %>%
    mutate(conc = 2^-7) %>%
    rename(ave_ratio_norm = ave_ratio_0_norm) %>%
    rename(ave_barcode = ave_barcode_0)
  df_2_5 <- df %>%
    mutate(ave_barcode_2_5 = (barcodes_RNA_2_5A + barcodes_RNA_2_5B)/2) %>%
    select(consensus, weak, nosite, most_common, background, ave_barcode_2_5, 
           ave_ratio_2_5_norm) %>%
    mutate(conc = 2^-5) %>%
    rename(ave_ratio_norm = ave_ratio_2_5_norm) %>%
    rename(ave_barcode = ave_barcode_2_5)
  df_2_4 <- df %>%
    mutate(ave_barcode_2_4 = (barcodes_RNA_2_4A + barcodes_RNA_2_4B)/2) %>%
    select(consensus, weak, nosite, most_common, background, ave_barcode_2_4, 
           ave_ratio_2_4_norm) %>%
    mutate(conc = 2^-4) %>%
    rename(ave_ratio_norm = ave_ratio_2_4_norm) %>%
    rename(ave_barcode = ave_barcode_2_4)
  df_2_3 <- df %>%
    mutate(ave_barcode_2_3 = (barcodes_RNA_2_3A + barcodes_RNA_2_3B)/2) %>%
    select(consensus, weak, nosite, most_common, background, ave_barcode_2_3, 
           ave_ratio_2_3_norm) %>%
    mutate(conc = 2^-3) %>%
    rename(ave_ratio_norm = ave_ratio_2_3_norm) %>%
    rename(ave_barcode = ave_barcode_2_3)
  df_2_2 <- df %>%
    mutate(ave_barcode_2_2 = (barcodes_RNA_2_2A + barcodes_RNA_2_2B)/2) %>%
    select(consensus, weak, nosite, most_common, background, ave_barcode_2_2, 
           ave_ratio_2_2_norm) %>%
    mutate(conc = 2^-2) %>%
    rename(ave_ratio_norm = ave_ratio_2_2_norm) %>%
    rename(ave_barcode = ave_barcode_2_2)
  df_2_1 <- df %>%
    mutate(ave_barcode_2_1 = (barcodes_RNA_2_1A + barcodes_RNA_2_1B)/2) %>%
    select(consensus, weak, nosite, most_common, background, ave_barcode_2_1, 
           ave_ratio_2_1_norm) %>%
    mutate(conc = 2^-1) %>%
    rename(ave_ratio_norm = ave_ratio_2_1_norm) %>%
    rename(ave_barcode = ave_barcode_2_1)
  df_20 <- df %>%
    mutate(ave_barcode_20 = (barcodes_RNA_20A + barcodes_RNA_20B)/2) %>%
    select(consensus, weak, nosite, most_common, background, ave_barcode_20, 
           ave_ratio_20_norm) %>%
    mutate(conc = 2^0) %>%
    rename(ave_ratio_norm = ave_ratio_20_norm) %>%
    rename(ave_barcode = ave_barcode_20)
  df_22 <- df %>%
    mutate(ave_barcode_22 = (barcodes_RNA_22A + barcodes_RNA_22B)/2) %>%
    select(consensus, weak, nosite, most_common, background, ave_barcode_22, 
           ave_ratio_22_norm) %>%
    mutate(conc = 2^2) %>%
    rename(ave_ratio_norm = ave_ratio_22_norm) %>%
    rename(ave_barcode = ave_barcode_22)
  df_0_22 <- rbind(df_0, df_2_5, df_2_4, df_2_3, df_2_2, df_2_1, df_20, df_22)
  return(df_0_22)
}

sp5_back_norm_conc_log2 <- sp5_conc_exp(subpool5) %>%
  mutate(conc = log2(conc)) %>%
  filter(weak == 0)

hill_fraction_sp5 <- function(df) {
  exp_2_2 <- df %>%
    filter(conc == 2) %>%
    select(most_common, consensus, background, ave_ratio_norm) %>%
    rename(ave_ratio_norm_2_2 = ave_ratio_norm) %>%
    mutate(ave_ratio_norm_2_2 = 1.1 * ave_ratio_norm_2_2) %>%
    filter(ave_ratio_norm_2_2 > 5)
  fraction <- right_join(df, exp_2_2, by = c('consensus', 'most_common','background')) %>%
    filter(conc != -7) %>%
    mutate(frac_active = log2((ave_ratio_norm/ave_ratio_norm_2_2)/(1 - ave_ratio_norm/ave_ratio_norm_2_2)))
  return(fraction)
}

trans_back_norm_conc_hillfrac <- hill_fraction_sp5(sp5_back_norm_conc_log2)

ggplot(trans_back_norm_conc_hillfrac, 
       aes(conc, frac_active, color = most_common)) +
  geom_point(show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ consensus) +
  scale_x_continuous(limits = c(-5, 0)) +
  geom_abline(slope = 1.5) +
  xlab('log2 forskolin µM') +
  ylab('log2 (y/1-y)') + 
  panel_border() +
  annotation_logticks(sides = 'bl') +
  background_grid(major = 'xy', minor = 'none')
  

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


