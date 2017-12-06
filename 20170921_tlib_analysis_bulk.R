library(tidyverse)
library(stringr)
library(viridis)
library(cowplot)
library(ggExtra)

#Written for the analysis of a range of inductions in the transient library
#DNA_BC: DNA
#R0A_BC and R0B_BC: RNA at 0 µM Forsk replicate A and B
#R2#A_BC and R2#B_BC: RNA at 2^# µM Forsk replicate A and B, includes negative and positive 
##'s

#tested concentrations: 0, 2^-5, 2^-4, 2^-3, 2^-2, 2^-1, 2^0, 2^2 µM Forsk


#Load index and bcmap files------------------------------------------------------------------

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


#Load barcode mapping table, remember sequences are rcomp due to sequencing format

barcode_map <- read_tsv('../../BCMap/uniqueSP2345.txt', 
                 col_names = c(
                   'fluff', 'barcode', 'name', 'most_common'
                   ), 
                 skip = 1) %>%
  select(-fluff)


#pick out SP3 and SP5 in the bcmap that were used in this assay

SP3_SP5_map <- barcode_map %>%
  mutate(subpool = ifelse(
    startsWith(name, 'subpool'), 
    substr(name, 1, 8), 
    'control'
    )
  ) %>%
  filter(subpool != 'subpool2') %>%
  filter(subpool != 'subpool4')


#Join reads to bcmap------------------------------------------------------------------------

#Join BC reads to BC mapping, keeping the reads only appearing in barcode mapping and 
#replacing na with 0 reads. I mixed the transfected and sequenced plasmid DNA separately 
#from one another. Should be the same ratio so decided to normalize reads in bulk, direct
#comparisons between subpools might not be as precise because of this but it allows us to 
#normaize reads to the background sequence in SP5.

bc_map_join_bc <- function(df1, df2) {
  keep_bc <- left_join(df1, df2, by = 'barcode') %>%
    mutate(num_reads = if_else(
      is.na(num_reads), 
      as.integer(0), 
      num_reads
      )
    ) %>%
    mutate(normalized = as.numeric((num_reads * 1000000) / (sum(num_reads))))
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


#Determine variant counts by summing--------------------------------------------------------

#sum unique barcodes and normalized bc reads per variant. Output is total barcodes and sum 
#normalized reads per variant. There is only 2 controls that are not represented by a BC read 
#in some of the samples, thus it is not present in joined tables later

var_sum_bc_num <- function(df1) {
  bc_count <- df1 %>%
    filter(df1$num_reads > 0) %>%
    group_by(subpool, name, most_common) %>%
    summarise(barcodes = n())
  variant_sum <- df1 %>%
    group_by(subpool, name, most_common) %>%
    count(name, wt = normalized) %>%
    rename(sum = n)
  bc_sum <- inner_join(variant_sum, bc_count, 
                       by = c("name", "subpool", "most_common")
                       ) %>%
    ungroup()
  return(bc_sum)
}

variant_counts_DNA <- var_sum_bc_num(bc_join_DNA)
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


#Join RNA to DNA and determine expression from summing--------------------------------------

#combine DNA and RNA cumm. BC counts, only keeping instances in both sets (as only 1 control 
#drops out in some samples) and determining RNA/DNA per variant. Ratio is summed normalized
#reads of RNA over DNA

var_expression <- function(df1, df2) {
  RNA_DNA <- inner_join(df1, df2, 
                        by = c("name", "subpool", "most_common"), 
                        suffix = c("_RNA", "_DNA")
  ) %>%
    mutate(ratio = sum_RNA / sum_DNA)
  print('x defined as RNA, y defined as DNA in var_expression(x,y)')
  return(RNA_DNA)
}

RNA_DNA_0A <- var_expression(variant_counts_R0A, variant_counts_DNA)
RNA_DNA_0B <- var_expression(variant_counts_R0B, variant_counts_DNA)
RNA_DNA_2_5A <- var_expression(variant_counts_R2_5A, variant_counts_DNA)
RNA_DNA_2_5B <- var_expression(variant_counts_R2_5B, variant_counts_DNA)
RNA_DNA_2_4A <- var_expression(variant_counts_R2_4A, variant_counts_DNA)
RNA_DNA_2_4B <- var_expression(variant_counts_R2_4B, variant_counts_DNA)
RNA_DNA_2_3A <- var_expression(variant_counts_R2_3A, variant_counts_DNA)
RNA_DNA_2_3B <- var_expression(variant_counts_R2_3B, variant_counts_DNA)
RNA_DNA_2_2A <- var_expression(variant_counts_R2_2A, variant_counts_DNA)
RNA_DNA_2_2B <- var_expression(variant_counts_R2_2B, variant_counts_DNA)
RNA_DNA_2_1A <- var_expression(variant_counts_R2_1A, variant_counts_DNA)
RNA_DNA_2_1B <- var_expression(variant_counts_R2_1B, variant_counts_DNA)
RNA_DNA_20A <- var_expression(variant_counts_R20A, variant_counts_DNA)
RNA_DNA_20B <- var_expression(variant_counts_R20B, variant_counts_DNA)
RNA_DNA_22A <- var_expression(variant_counts_R22A, variant_counts_DNA)
RNA_DNA_22B <- var_expression(variant_counts_R22B, variant_counts_DNA)


#Tidy data of combined conc. and rep--------------------------------------------------------

var_conc_rep <- function(
  df0A, df0B, df2_5A, df2_5B, df2_4A, df2_4B, df2_3A, df2_3B, df2_2A, df2_2B, df2_1A, df2_1B, 
  df20A, df20B, df22A, df22B
  ) {
  join_0 <- inner_join(df0A, df0B, 
                       by = c(
                         "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                         ), suffix = c("_0A", "_0B")
                       )
  join_2_5 <- inner_join(df2_5A, df2_5B, 
                         by = c(
                           "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                           ), suffix = c("_2_5A", "_2_5B")
                         )
  join_2_4 <- inner_join(df2_4A, df2_4B, 
                         by = c(
                           "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                           ), suffix = c("_2_4A", "_2_4B")
                         )
  join_2_3 <- inner_join(df2_3A, df2_3B, 
                         by = c(
                           "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                           ), suffix = c("_2_3A", "_2_3B")
                         )
  join_2_2 <- inner_join(df2_2A, df2_2B, 
                         by = c(
                           "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                           ), suffix = c("_2_2A", "_2_2B")
                         )
  join_2_1 <- inner_join(df2_1A, df2_1B, 
                         by = c(
                           "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                           ), suffix = c("_2_1A", "_2_1B")
                         )
  join_20 <- inner_join(df20A, df20B, 
                        by = c(
                          "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                          ), suffix = c("_20A", "_20B")
                        )
  join_22 <- inner_join(df22A, df22B, 
                        by = c(
                          "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                          ), suffix = c("_22A", "_22B")
                        )
  join_0_2_5 <- inner_join(join_0, join_2_5, 
                           by = c(
                             "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                             )
                           )
  join_0_2_4 <- inner_join(join_0_2_5, join_2_4, 
                           by = c(
                             "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                             )
                           )
  join_0_2_3 <- inner_join(join_0_2_4, join_2_3, 
                           by = c(
                             "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                             )
                           )
  join_0_2_2 <- inner_join(join_0_2_3, join_2_2, 
                           by = c(
                             "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                             )
                           )
  join_0_2_1 <- inner_join(join_0_2_2, join_2_1, 
                           by = c(
                             "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                             )
                           )
  join_0_20 <- inner_join(join_0_2_1, join_20, 
                          by = c(
                            "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                            )
                          )
  join_0_22 <- inner_join(join_0_20, join_22, 
                          by = c(
                            "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                            )
                          )
  print(
    'processed dfs in order of samples: 0A, 0B, 2_5A, 2_5B, 2_4A, 2_4B, 2_3A, 2_3B, 2_2A, 
    2_2B, 2_1A, 2_1B, 20A, 20B, 22A, 22B'
    )
  return(join_0_22)
}

rep_0_22_A_B <- var_conc_rep(RNA_DNA_0A, RNA_DNA_0B, RNA_DNA_2_5A, RNA_DNA_2_5B,
                             RNA_DNA_2_4A, RNA_DNA_2_4B, RNA_DNA_2_3A, RNA_DNA_2_3B,
                             RNA_DNA_2_2A, RNA_DNA_2_2B, RNA_DNA_2_1A, RNA_DNA_2_1B,
                             RNA_DNA_20A, RNA_DNA_20B, RNA_DNA_22A, RNA_DNA_22B)


#determine the log(RNA/DNA) for each sample

var_log <- function(df) {
  log_ratio_df <- rep_0_22_A_B %>% 
    mutate_if(is.double, 
              funs(log10(.))
              )
  return(log_ratio_df)
}

log_rep_0_22_A_B <- var_log(rep_0_22_A_B)


#BC analysis---------------------------------------------------------------------------------

#Determining subpool proportions in DNA

DNA_nonzero <- bc_join_DNA %>%
  filter(num_reads > 0) %>%
  filter(subpool == 'subpool5')


#Plot reads per BC

p_BC_num_reads_viol_full <- ggplot(NULL, aes(x = "", y = num_reads)) +
  geom_violin(data = bc_join_DNA, aes(x = "DNA", color = subpool)) +
  geom_violin(data = bc_join_R0A, aes(x = "R0A", color = subpool)) + 
  geom_violin(data = bc_join_R0B, aes(x = "R0B", color = subpool)) + 
  geom_violin(data = bc_join_R2_5A, aes(x = "R2^-5A", color = subpool)) + 
  geom_violin(data = bc_join_R2_5B, aes(x = "R2^-5B", color = subpool)) + 
  geom_violin(data = bc_join_R2_4A, aes(x = "R2^-4A", color = subpool)) +
  geom_violin(data = bc_join_R2_4B, aes(x = "R2^-4B", color = subpool)) + 
  geom_violin(data = bc_join_R2_3A, aes(x = "R2^-3A", color = subpool)) + 
  geom_violin(data = bc_join_R2_3B, aes(x = "R2^-3B", color = subpool)) + 
  geom_violin(data = bc_join_R2_2A, aes(x = "R2^-2A", color = subpool)) + 
  geom_violin(data = bc_join_R2_2B, aes(x = "R2^-2B", color = subpool)) + 
  geom_violin(data = bc_join_R2_1A, aes(x = "R2^-1A", color = subpool)) + 
  geom_violin(data = bc_join_R2_1B, aes(x = "R2^-1B", color = subpool)) + 
  geom_violin(data = bc_join_R20A, aes(x = "R2^0A", color = subpool)) + 
  geom_violin(data = bc_join_R20B, aes(x = "R2^0B", color = subpool)) + 
  geom_violin(data = bc_join_R22A, aes(x = "R2^2A", color = subpool)) + 
  geom_violin(data = bc_join_R22B, aes(x = "R2^2B", color = subpool)) +
  xlab("") +
  ylab("Reads per BC")

save_plot('plots/BC_num_reads_viol_full.png',
          p_BC_num_reads_viol_full, scale = 2.8)


p_BC_num_reads_box_zoom <- ggplot(NULL, aes(x = "", y = num_reads)) +
  geom_boxplot(data = bc_join_DNA, aes(x = "DNA", color = subpool)) +
  geom_boxplot(data = bc_join_R0A, aes(x = "R0A", color = subpool)) + 
  geom_boxplot(data = bc_join_R0B, aes(x = "R0B", color = subpool)) + 
  geom_boxplot(data = bc_join_R2_5A, aes(x = "R2^-5A", color = subpool)) + 
  geom_boxplot(data = bc_join_R2_5B, aes(x = "R2^-5B", color = subpool)) + 
  geom_boxplot(data = bc_join_R2_4A, aes(x = "R2^-4A", color = subpool)) +
  geom_boxplot(data = bc_join_R2_4B, aes(x = "R2^-4B", color = subpool)) + 
  geom_boxplot(data = bc_join_R2_3A, aes(x = "R2^-3A", color = subpool)) + 
  geom_boxplot(data = bc_join_R2_3B, aes(x = "R2^-3B", color = subpool)) + 
  geom_boxplot(data = bc_join_R2_2A, aes(x = "R2^-2A", color = subpool)) + 
  geom_boxplot(data = bc_join_R2_2B, aes(x = "R2^-2B", color = subpool)) + 
  geom_boxplot(data = bc_join_R2_1A, aes(x = "R2^-1A", color = subpool)) + 
  geom_boxplot(data = bc_join_R2_1B, aes(x = "R2^-1B", color = subpool)) + 
  geom_boxplot(data = bc_join_R20A, aes(x = "R2^0A", color = subpool)) + 
  geom_boxplot(data = bc_join_R20B, aes(x = "R2^0B", color = subpool)) + 
  geom_boxplot(data = bc_join_R22A, aes(x = "R2^2A", color = subpool)) + 
  geom_boxplot(data = bc_join_R22B, aes(x = "R2^2B", color = subpool)) +
  xlab("") +
  ylab("Reads") + ylim(0, 25)

save_plot('plots/BC_num_reads_box_zoom.png',
          p_BC_num_reads_box_zoom, scale = 2.8)


#replicate plots----------------------------------------------------------------------------

#plot replicates for summed variant expression

p_var_rep_0 <- ggplot(NULL, aes(ratio_0A, ratio_0B)) +
  geom_point(data = log_rep_0_22_A_B, alpha = 0.3) +
  geom_point(data = filter(log_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expr. Rep. 1") +
  ylab("Variant log10 Expr. Rep. 2") +
  scale_x_continuous(breaks = c(0, 1, 2), limits = c(-2, 2)) + 
  scale_y_continuous(breaks = c(0, 1, 2), limits = c(-2, 2)) + 
  annotate("text", x = -1, y = 1,
           label = paste(
             'r =', round(
               cor(
                 log_rep_0_22_A_B$ratio_0A,log_rep_0_22_A_B$ratio_0B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_rep_2_5 <- ggplot(NULL, aes(ratio_2_5A, ratio_2_5B)) +
  geom_point(data = log_rep_0_22_A_B, alpha = 0.3) +
  geom_point(data = filter(log_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expr. Rep. 1") +
  ylab("Variant log10 Expr. Rep. 2") +
  scale_x_continuous(breaks = c(0, 1, 2), limits = c(-2, 2)) + 
  scale_y_continuous(breaks = c(0, 1, 2), limits = c(-2, 2)) + 
  annotate("text", x = -1, y = 1,
           label = paste(
             'r =', round(
               cor(
                 log_rep_0_22_A_B$ratio_2_5A,log_rep_0_22_A_B$ratio_2_5B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_rep_2_4 <- ggplot(NULL, aes(ratio_2_4A, ratio_2_4B)) +
  geom_point(data = log_rep_0_22_A_B, alpha = 0.3) +
  geom_point(data = filter(log_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expr. Rep. 1") +
  ylab("Variant log10 Expr. Rep. 2") +
  scale_x_continuous(breaks = c(0, 1, 2), limits = c(-2, 2)) + 
  scale_y_continuous(breaks = c(0, 1, 2), limits = c(-2, 2)) + 
  annotate("text", x = -1, y = 1,
           label = paste(
             'r =', round(
               cor(
                 log_rep_0_22_A_B$ratio_2_4A,log_rep_0_22_A_B$ratio_2_4B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_rep_2_3 <- ggplot(NULL, aes(ratio_2_3A, ratio_2_3B)) +
  geom_point(data = log_rep_0_22_A_B, alpha = 0.3) +
  geom_point(data = filter(log_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expr. Rep. 1") +
  ylab("Variant log10 Expr. Rep. 2") +
  scale_x_continuous(breaks = c(0, 1, 2), limits = c(-2, 2)) + 
  scale_y_continuous(breaks = c(0, 1, 2), limits = c(-2, 2)) + 
  annotate("text", x = -1, y = 1,
           label = paste(
             'r =', round(
               cor(
                 log_rep_0_22_A_B$ratio_2_3A,log_rep_0_22_A_B$ratio_2_3B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_rep_2_2 <- ggplot(NULL, aes(ratio_2_2A, ratio_2_2B)) +
  geom_point(data = log_rep_0_22_A_B, alpha = 0.3) +
  geom_point(data = filter(log_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expr. Rep. 1") +
  ylab("Variant log10 Expr. Rep. 2") +
  scale_x_continuous(breaks = c(0, 1, 2), limits = c(-2, 2)) + 
  scale_y_continuous(breaks = c(0, 1, 2), limits = c(-2, 2)) + 
  annotate("text", x = -1, y = 1,
           label = paste(
             'r =', round(
               cor(
                 log_rep_0_22_A_B$ratio_2_2A,log_rep_0_22_A_B$ratio_2_2B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_rep_2_1 <- ggplot(NULL, aes(ratio_2_1A, ratio_2_1B)) +
  geom_point(data = log_rep_0_22_A_B, alpha = 0.3) +
  geom_point(data = filter(log_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expr. Rep. 1") +
  ylab("Variant log10 Expr. Rep. 2") +
  scale_x_continuous(breaks = c(0, 1, 2), limits = c(-2, 2)) + 
  scale_y_continuous(breaks = c(0, 1, 2), limits = c(-2, 2)) + 
  annotate("text", x = -1, y = 1,
           label = paste(
             'r =', round(
               cor(
                 log_rep_0_22_A_B$ratio_2_1A,log_rep_0_22_A_B$ratio_2_1B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_rep_20 <- ggplot(NULL, aes(ratio_20A, ratio_20B)) +
  geom_point(data = log_rep_0_22_A_B, alpha = 0.3) +
  geom_point(data = filter(log_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expr. Rep. 1") +
  ylab("Variant log10 Expr. Rep. 2") +
  scale_x_continuous(breaks = c(0, 1, 2), limits = c(-2, 2)) + 
  scale_y_continuous(breaks = c(0, 1, 2), limits = c(-2, 2)) + 
  annotate("text", x = -1, y = 1,
           label = paste(
             'r =', round(
               cor(
                 log_rep_0_22_A_B$ratio_20A,log_rep_0_22_A_B$ratio_20B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_rep_22 <- ggplot(NULL, aes(ratio_22A, ratio_22B)) +
  geom_point(data = log_rep_0_22_A_B, alpha = 0.3) +
  geom_point(data = filter(log_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expr. Rep. 1") +
  ylab("Variant log10 Expr. Rep. 2") +
  scale_x_continuous(breaks = c(0, 1, 2), limits = c(-2, 2)) + 
  scale_y_continuous(breaks = c(0, 1, 2), limits = c(-2, 2)) + 
  annotate("text", x = -1, y = 1,
           label = paste(
             'r =', round(
               cor(
                 log_rep_0_22_A_B$ratio_22A,log_rep_0_22_A_B$ratio_22B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_log_var_rep_grid <- plot_grid(
  p_var_rep_0, p_var_rep_2_5, p_var_rep_2_4, p_var_rep_2_3, p_var_rep_2_2, p_var_rep_2_1, 
  p_var_rep_20, p_var_rep_22, 
  labels = c(
    "0 µM", "2^-5 µM", "2^-4 µM","2^-3 µM", "2^-2 µM", "2^-1 µM","2^0 µM", "2^2 µM"),
  nrow = 3, ncol = 3, align = 'hv', hjust = -3, vjust = 0.5, scale = 0.9)

save_plot('plots/p_log_var_rep_grid.png', 
          p_log_var_rep_grid, base_height = 10, base_width = 10)


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


#Comparing sub-pool to bulk normalization
p_log_ratio_norm_0A <- ggplot(log_rep_0_22_A_B, 
                              aes(subpool, log_ratio_ratio_0A)) +
  geom_boxplot() +
  geom_hline(yintercept = 1) +
  ylab("Ratio of log10 Expr bulk norm. to\nlog10 Expr SP norm. Rep. 2") +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-3, 3))

p_log_ratio_norm_0B <- ggplot(log_rep_0_22_A_B, 
                              aes(subpool, log_ratio_ratio_0B)) +
  geom_boxplot() +
  geom_hline(yintercept = 1) +
  ylab("Ratio of log10 Expr bulk norm. to\nlog10 Expr SP norm. Rep. 2") +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-3, 3))

p_log_ratio_norm_0_grid <- plot_grid(p_log_ratio_norm_0A, p_log_ratio_norm_0B,
  labels = c("0 µM rep. A", "0 µM rep. B"),
  nrow = 1, ncol = 2, align = 'hv', hjust = -1, vjust = 0.5, scale = 0.9)

save_plot('plots/p_log_ratio_norm_0_grid.png', 
          p_log_ratio_norm_0_grid, base_height = 4, base_width = 9)

p_log_ratio_norm_2_5A <- ggplot(log_rep_0_22_A_B, aes(
  subpool, log_ratio_ratio_2_5A)) +
  geom_boxplot() +
  ylab("Ratio of log10 Expr bulk norm. to\nlog10 Expr SP norm. Rep. 2") +
  scale_y_continuous(limits = c(-3, 3))

p_log_ratio_norm_2_5B <- ggplot(log_rep_0_22_A_B, aes(
  subpool, log_ratio_ratio_2_5B)) +
  geom_boxplot() +
  ylab("Ratio of log10 Expr bulk norm. to\nlog10 Expr SP norm. Rep. 2") +
  scale_y_continuous(limits = c(-3, 3))

p_log_ratio_norm_2_4A <- ggplot(log_rep_0_22_A_B, aes(
  subpool, log_ratio_ratio_2_4A)) +
  geom_boxplot() +
  ylab("Ratio of log10 Expr bulk norm. to\nlog10 Expr SP norm. Rep. 2") +
  scale_y_continuous(limits = c(-3, 3))

p_log_ratio_norm_2_4B <- ggplot(log_rep_0_22_A_B, aes(
  subpool, log_ratio_ratio_2_4B)) +
  geom_boxplot() +
  ylab("Ratio of log10 Expr bulk norm. to\nlog10 Expr SP norm. Rep. 2") +
  scale_y_continuous(limits = c(-3, 3))

p_log_ratio_norm_2_3A <- ggplot(log_rep_0_22_A_B, aes(
  subpool, log_ratio_ratio_2_3A)) +
  geom_boxplot() +
  ylab("Ratio of log10 Expr bulk norm. to\nlog10 Expr SP norm. Rep. 2") +
  scale_y_continuous(limits = c(-3, 3))

p_log_ratio_norm_2_3B <- ggplot(log_rep_0_22_A_B, aes(
  subpool, log_ratio_ratio_2_3B)) +
  geom_boxplot() +
  ylab("Ratio of log10 Expr bulk norm. to\nlog10 Expr SP norm. Rep. 2") +
  scale_y_continuous(limits = c(-3, 3))

p_log_ratio_norm_2_2A <- ggplot(log_rep_0_22_A_B, aes(
  subpool, log_ratio_ratio_2_2A)) +
  geom_boxplot() +
  ylab("Ratio of log10 Expr bulk norm. to\nlog10 Expr SP norm. Rep. 2") +
  scale_y_continuous(limits = c(-3, 3))

p_log_ratio_norm_2_2B <- ggplot(log_rep_0_22_A_B, aes(
  subpool, log_ratio_ratio_2_2B)) +
  geom_boxplot() +
  ylab("Ratio of log10 Expr bulk norm. to\nlog10 Expr SP norm. Rep. 2") +
  scale_y_continuous(limits = c(-3, 3))

p_log_ratio_norm_2_1A <- ggplot(log_rep_0_22_A_B, aes(
  subpool, log_ratio_ratio_2_1A)) +
  geom_boxplot() +
  ylab("Ratio of log10 Expr bulk norm. to\nlog10 Expr SP norm. Rep. 2") +
  scale_y_continuous(limits = c(-3, 3))

p_log_ratio_norm_2_1B <- ggplot(log_rep_0_22_A_B, aes(
  subpool, log_ratio_ratio_2_1B)) +
  geom_boxplot() +
  ylab("Ratio of log10 Expr bulk norm. to\nlog10 Expr SP norm. Rep. 2") +
  scale_y_continuous(limits = c(-3, 3))

p_log_ratio_norm_20A <- ggplot(log_rep_0_22_A_B, aes(
  subpool, log_ratio_ratio_20A)) +
  geom_boxplot() +
  ylab("Ratio of log10 Expr bulk norm. to\nlog10 Expr SP norm. Rep. 2") +
  scale_y_continuous(limits = c(-3, 3))

p_log_ratio_norm_20B <- ggplot(log_rep_0_22_A_B, aes(
  subpool, log_ratio_ratio_20B)) +
  geom_boxplot() +
  ylab("Ratio of log10 Expr bulk norm. to\nlog10 Expr SP norm. Rep. 2") +
  scale_y_continuous(limits = c(-3, 3))

p_log_ratio_norm_22A <- ggplot(log_rep_0_22_A_B, aes(
  subpool, log_ratio_ratio_22A)) +
  geom_boxplot() +
  ylab("Ratio of log10 Expr bulk norm. to\nlog10 Expr SP norm. Rep. 2") +
  scale_y_continuous(limits = c(-3, 3))

p_log_ratio_norm_22B <- ggplot(log_rep_0_22_A_B, aes(
  subpool, log_ratio_ratio_22B)) +
  geom_boxplot() +
  ylab("Ratio of log10 Expr bulk norm. to\nlog10 Expr SP norm. Rep. 2") +
  scale_y_continuous(limits = c(-3, 3))

p_log_ratio_norm_grid <- plot_grid(
  p_log_ratio_norm_0A, p_log_ratio_norm_0B, p_log_ratio_norm_2_5A,
  p_log_ratio_norm_2_5B, p_log_ratio_norm_2_4A, p_log_ratio_norm_2_4B,
  p_log_ratio_norm_2_3A, p_log_ratio_norm_2_3B, p_log_ratio_norm_2_2A,
  p_log_ratio_norm_2_2B, p_log_ratio_norm_2_1A, p_log_ratio_norm_2_1B,
  p_log_ratio_norm_20A, p_log_ratio_norm_20B, p_log_ratio_norm_22A,
  p_log_ratio_norm_22B,
  labels = c(
    "0 µM rep. A", "0 µM rep. B", "2^-5 µM rep. A", "2^-5 µM rep. B", "2^-4 µM rep. A", 
    "2^-4 µM rep. B", "2^-3 µM rep. A", "2^-3 µM rep. B", "2^-2 µM rep. A", 
    "2^-2 µM rep. B", "2^-1 µM rep. A", "2^-1 µM rep. B", "2^0 µM rep. A", 
    "2^0 µM rep. B", "2^2 µM rep. A", "2^2 µM rep. B"),
  nrow = 4, ncol = 4, align = 'hv', hjust = -3, vjust = 0.5, scale = 0.9)

save_plot('plots/p_log_ratio_norm_grid.png', 
          p_log_ratio_norm_grid, base_height = 14, base_width = 14)


#High uninduced high induced variant for Rishi----------------------
reportcomp <- filter(log_rep_0_64_A_B, 
                     name == 'subpool5_weak_consensus_consensus_consensus_consensus_consensus_Vista Chr5:88673410-88674494' | name == 'pGL4.29 Promega 1-63 + 1-87' | name == 'pGL4.29 Promega 1-87')



#Untidy joining-----------------------------------------------------

var_expression <- function(df1, df2) {
  RNA_DNA <- inner_join(
    df1, df2, by = c("name", "subpool", "most_common"), 
    suffix = c("_RNA", "_DNA")
  ) %>%
    mutate(ratio = n_RNA / n_DNA)
  print('x defined as RNA, y defined as DNA in var_expression(x,y)')
  return(RNA_DNA)
}

RNA_DNA_0A <- var_expression(variant_counts_R0A, variant_counts_DNA) %>%
  mutate(conc = 0) %>%
  mutate(rep = 'A')
RNA_DNA_0B <- var_expression(variant_counts_R0B, variant_counts_DNA) %>%
  mutate(conc = 0) %>%
  mutate(rep = 'B')
RNA_DNA_2_5A <- var_expression(variant_counts_R2_5A, variant_counts_DNA) %>%
  mutate(conc = as.numeric(2^-5)) %>%
  mutate(rep = 'A')
RNA_DNA_2_5B <- var_expression(variant_counts_R2_5B, variant_counts_DNA) %>%
  mutate(conc = as.numeric(2^-5)) %>%
  mutate(rep = 'B')
RNA_DNA_2_4A <- var_expression(variant_counts_R2_4A, variant_counts_DNA) %>%
  mutate(conc = as.numeric(2^-4)) %>%
  mutate(rep = 'A')
RNA_DNA_2_4B <- var_expression(variant_counts_R2_4B, variant_counts_DNA) %>%
  mutate(conc = as.numeric(2^-4)) %>%
  mutate(rep = 'B')
RNA_DNA_2_3A <- var_expression(variant_counts_R2_3A, variant_counts_DNA) %>%
  mutate(conc = as.numeric(2^-3)) %>%
  mutate(rep = 'A')
RNA_DNA_2_3B <- var_expression(variant_counts_R2_3B, variant_counts_DNA) %>%
  mutate(conc = as.numeric(2^-3)) %>%
  mutate(rep = 'B')
RNA_DNA_2_2A <- var_expression(variant_counts_R2_2A, variant_counts_DNA) %>%
  mutate(conc = as.numeric(2^-2)) %>%
  mutate(rep = 'A')
RNA_DNA_2_2B <- var_expression(variant_counts_R2_2B, variant_counts_DNA) %>%
  mutate(conc = as.numeric(2^-2)) %>%
  mutate(rep = 'B')
RNA_DNA_2_1A <- var_expression(variant_counts_R2_1A, variant_counts_DNA) %>%
  mutate(conc = as.numeric(2^-1)) %>%
  mutate(rep = 'A')
RNA_DNA_2_1B <- var_expression(variant_counts_R2_1B, variant_counts_DNA) %>%
  mutate(conc = as.numeric(2^-1)) %>%
  mutate(rep = 'B')
RNA_DNA_20A <- var_expression(variant_counts_R20A, variant_counts_DNA) %>%
  mutate(conc = as.numeric(2^0)) %>%
  mutate(rep = 'A')
RNA_DNA_20B <- var_expression(variant_counts_R20B, variant_counts_DNA) %>%
  mutate(conc = as.numeric(2^0)) %>%
  mutate(rep = 'B')
RNA_DNA_22A <- var_expression(variant_counts_R22A, variant_counts_DNA) %>%
  mutate(conc = as.numeric(2^2)) %>%
  mutate(rep = 'A')
RNA_DNA_22B <- var_expression(variant_counts_R22B, variant_counts_DNA) %>%
  mutate(conc = as.numeric(2^2)) %>%
  mutate(rep = 'B')


#Untidy data of combined conc. and rep------------------------------------------------------
var_conc_rep <- function(
  df0A, df0B, df2_5A, df2_5B, df2_4A, df2_4B, df2_3A, df2_3B, df2_2A, df2_2B, df2_1A, df2_1B, 
  df20A, df20B, df22A, df22B) {
  join_0 <- inner_join(
    df0A, df0B, by = c(
      "name", "subpool", "most_common", "conc", "n_DNA", "barcodes_DNA"
    ), suffix = c("_0A", "_0B")
  )
  join_2_5 <- inner_join(
    df2_5A, df2_5B, by = c(
      "name", "subpool", "most_common", "conc", "n_DNA", "barcodes_DNA"
    ), suffix = c("_2_5A", "_2_5B")
  )
  join_2_4 <- inner_join(
    df2_4A, df2_4B, by = c(
      "name", "subpool", "most_common", "conc", "n_DNA", "barcodes_DNA"
    ), suffix = c("_2_4A", "_2_4B")
  )
  join_2_3 <- inner_join(
    df2_3A, df2_3B, by = c(
      "name", "subpool", "most_common", "conc", "n_DNA", "barcodes_DNA"
    ), suffix = c("_2_3A", "_2_3B")
  )
  join_2_2 <- inner_join(
    df2_2A, df2_2B, by = c(
      "name", "subpool", "most_common", "conc", "n_DNA", "barcodes_DNA"
    ), suffix = c("_2_2A", "_2_2B")
  )
  join_2_1 <- inner_join(
    df2_1A, df2_1B, by = c(
      "name", "subpool", "most_common", "conc", "n_DNA", "barcodes_DNA"
    ), suffix = c("_2_1A", "_2_1B")
  )
  join_20 <- inner_join(
    df20A, df20B, by = c(
      "name", "subpool", "most_common", "conc", "n_DNA", "barcodes_DNA"
    ), suffix = c("_20A", "_20B")
  )
  join_22 <- inner_join(
    df22A, df22B, by = c(
      "name", "subpool", "most_common", "conc", "n_DNA", "barcodes_DNA"
    ), suffix = c("_22A", "_22B")
  )
  join_0_2_5 <- inner_join(
    join_0, join_2_5, by = c(
      "name", "subpool", "most_common", "n_DNA", "barcodes_DNA"
    ), suffix = c('')
  )
  join_0_2_4 <- inner_join(
    join_0_2_5, join_2_4, by = c(
      "name", "subpool", "most_common", "n_DNA", "barcodes_DNA"
    )
  )
  join_0_2_3 <- inner_join(
    join_0_2_4, join_2_3, by = c(
      "name", "subpool", "most_common", "n_DNA", "barcodes_DNA"
    )
  )
  join_0_2_2 <- inner_join(
    join_0_2_3, join_2_2, by = c(
      "name", "subpool", "most_common", "n_DNA", "barcodes_DNA"
    )
  )
  join_0_2_1 <- inner_join(
    join_0_2_2, join_2_1, by = c(
      "name", "subpool", "most_common", "n_DNA", "barcodes_DNA"
    )
  )
  join_0_20 <- inner_join(
    join_0_2_1, join_20, by = c(
      "name", "subpool", "most_common", "n_DNA", "barcodes_DNA"
    )
  )
  join_0_22 <- inner_join(
    join_0_20, join_22, by = c(
      "name", "subpool", "most_common", "n_DNA", "barcodes_DNA"
    )
  )
  print(
    'processed dfs in format: 0A, 0B, 2_5A, 2_5B, 2_4A, 2_4B, 2_3A, 2_3B, 2_2A, 2_2B, 2_1A, 
    2_1B, 20A, 20B, 22A, 22B'
  )
  return(join_0_22)
}

rep_0_22_A_B <- var_conc_rep(RNA_DNA_0A, RNA_DNA_0B, RNA_DNA_2_5A, RNA_DNA_2_5B,
                             RNA_DNA_2_4A, RNA_DNA_2_4B, RNA_DNA_2_3A, RNA_DNA_2_3B,
                             RNA_DNA_2_2A, RNA_DNA_2_2B, RNA_DNA_2_1A, RNA_DNA_2_1B,
                             RNA_DNA_20A, RNA_DNA_20B, RNA_DNA_22A, RNA_DNA_22B)

