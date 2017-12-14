library(tidyverse)
library(stringr)
library(viridis)
library(cowplot)
library(ggExtra)

#Written for the bc analysis of a range of inductions in the transient library
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
#replacing na with pseudocount of 1 read. I remove barcodes with 1 pseudocount read in DNA
#as I can't determine their expression. I mixed the transfected and sequenced plasmid DNA 
#separately from one another. Should be the same ratio so decided to normalize reads in bulk, 
#direct comparisons between subpools might not be as precise because of this but it allows 
#us to normalize reads to the background sequence in SP5.

bc_map_join_bc <- function(df1, df2) {
  keep_bc <- left_join(df1, df2, by = 'barcode') %>%
    mutate(num_reads = if_else(
      is.na(num_reads), 
      as.integer(0), 
      num_reads
      )
    ) %>%
    mutate(norm = as.numeric((num_reads * 1000000) / (sum(num_reads))))
  return(keep_bc)
}

bc_join_DNA <- bc_map_join_bc(SP3_SP5_map, bc_DNA)
bc_join_DNA <- filter(bc_join_DNA, num_reads > 0)

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


#Join RNA to DNA and determine expression----------------------------------------------------

#combine DNA and RNA normalized reads, only keeping instances in both sets (as only 1 control 
#drops out in some samples) and determining RNA/DNA per barcode. Ratio is normalized reads of
#RNA over DNA per BC

bc_expression <- function(df1, df2) {
  RNA_DNA <- right_join(df1, df2, 
                        by = c("barcode", "name", "subpool", "most_common"), 
                        suffix = c("_RNA", "_DNA")
                        ) %>%
    mutate(ratio = norm_RNA / norm_DNA)
  print('x defined as RNA, y defined as DNA in bc_expression(x,y)')
  return(RNA_DNA)
}

RNA_DNA_0A <- bc_expression(bc_join_R0A, bc_join_DNA)
RNA_DNA_0B <- bc_expression(bc_join_R0B, bc_join_DNA)
RNA_DNA_2_5A <- bc_expression(bc_join_R2_5A, bc_join_DNA)
RNA_DNA_2_5B <- bc_expression(bc_join_R2_5B, bc_join_DNA)
RNA_DNA_2_4A <- bc_expression(bc_join_R2_4A, bc_join_DNA)
RNA_DNA_2_4B <- bc_expression(bc_join_R2_4B, bc_join_DNA)
RNA_DNA_2_3A <- bc_expression(bc_join_R2_3A, bc_join_DNA)
RNA_DNA_2_3B <- bc_expression(bc_join_R2_3B, bc_join_DNA)
RNA_DNA_2_2A <- bc_expression(bc_join_R2_2A, bc_join_DNA)
RNA_DNA_2_2B <- bc_expression(bc_join_R2_2B, bc_join_DNA)
RNA_DNA_2_1A <- bc_expression(bc_join_R2_1A, bc_join_DNA)
RNA_DNA_2_1B <- bc_expression(bc_join_R2_1B, bc_join_DNA)
RNA_DNA_20A <- bc_expression(bc_join_R20A, bc_join_DNA)
RNA_DNA_20B <- bc_expression(bc_join_R20B, bc_join_DNA)
RNA_DNA_22A <- bc_expression(bc_join_R22A, bc_join_DNA)
RNA_DNA_22B <- bc_expression(bc_join_R22B, bc_join_DNA)


#Tidy data of combined conc. and rep--------------------------------------------------------

bc_conc_rep <- function(
  df0A, df0B, df2_5A, df2_5B, df2_4A, df2_4B, df2_3A, df2_3B, df2_2A, df2_2B, df2_1A, df2_1B, 
  df20A, df20B, df22A, df22B
  ) {
  join_0 <- inner_join(df0A, df0B, 
                       by = c(
                         "barcode", "name", "subpool", "most_common", "num_reads_DNA", 
                         "norm_DNA"
                         ), suffix = c("_0A", "_0B")
                       )
  join_2_5 <- inner_join(df2_5A, df2_5B, 
                         by = c(
                           "barcode", "name", "subpool", "most_common", "num_reads_DNA", 
                           "norm_DNA"
                           ), suffix = c("_2_5A", "_2_5B")
                         )
  join_2_4 <- inner_join(df2_4A, df2_4B, 
                         by = c(
                           "barcode", "name", "subpool", "most_common", "num_reads_DNA", 
                           "norm_DNA"
                           ), suffix = c("_2_4A", "_2_4B")
                         )
  join_2_3 <- inner_join(df2_3A, df2_3B, 
                         by = c(
                           "barcode", "name", "subpool", "most_common", "num_reads_DNA", 
                           "norm_DNA"
                           ), suffix = c("_2_3A", "_2_3B")
                         )
  join_2_2 <- inner_join(df2_2A, df2_2B, 
                         by = c(
                           "barcode", "name", "subpool", "most_common", "num_reads_DNA", 
                           "norm_DNA"
                           ), suffix = c("_2_2A", "_2_2B")
                         )
  join_2_1 <- inner_join(df2_1A, df2_1B, 
                         by = c(
                           "barcode", "name", "subpool", "most_common", "num_reads_DNA", 
                           "norm_DNA"
                           ), suffix = c("_2_1A", "_2_1B")
                         )
  join_20 <- inner_join(df20A, df20B, 
                        by = c(
                          "barcode", "name", "subpool", "most_common", "num_reads_DNA", 
                          "norm_DNA"
                          ), suffix = c("_20A", "_20B")
                        )
  join_22 <- inner_join(df22A, df22B, 
                        by = c(
                          "barcode", "name", "subpool", "most_common", "num_reads_DNA", 
                          "norm_DNA"
                          ), suffix = c("_22A", "_22B")
                        )
  join_0_2_5 <- inner_join(join_0, join_2_5, 
                           by = c(
                             "barcode", "name", "subpool", "most_common", "num_reads_DNA", 
                             "norm_DNA"
                             )
                           )
  join_0_2_4 <- inner_join(join_0_2_5, join_2_4, 
                           by = c(
                             "barcode", "name", "subpool", "most_common", "num_reads_DNA", 
                             "norm_DNA"
                             )
                           )
  join_0_2_3 <- inner_join(join_0_2_4, join_2_3, 
                           by = c(
                             "barcode", "name", "subpool", "most_common", "num_reads_DNA", 
                             "norm_DNA"
                             )
                           )
  join_0_2_2 <- inner_join(join_0_2_3, join_2_2, 
                           by = c(
                             "barcode", "name", "subpool", "most_common", "num_reads_DNA", 
                             "norm_DNA"
                             )
                           )
  join_0_2_1 <- inner_join(join_0_2_2, join_2_1, 
                           by = c(
                             "barcode", "name", "subpool", "most_common", "num_reads_DNA", 
                             "norm_DNA"
                             )
                           )
  join_0_20 <- inner_join(join_0_2_1, join_20, 
                          by = c(
                            "barcode", "name", "subpool", "most_common", "num_reads_DNA", 
                            "norm_DNA"
                            )
                          )
  join_0_22 <- inner_join(join_0_20, join_22, 
                          by = c(
                            "barcode", "name", "subpool", "most_common", "num_reads_DNA", 
                            "norm_DNA"
                            )
                          )
  print(
    'processed dfs in order of samples: 0A, 0B, 2_5A, 2_5B, 2_4A, 2_4B, 2_3A, 2_3B, 2_2A, 
    2_2B, 2_1A, 2_1B, 20A, 20B, 22A, 22B'
    )
  return(join_0_22)
}

rep_0_22_A_B <- bc_conc_rep(RNA_DNA_0A, RNA_DNA_0B, RNA_DNA_2_5A, RNA_DNA_2_5B,
                             RNA_DNA_2_4A, RNA_DNA_2_4B, RNA_DNA_2_3A, RNA_DNA_2_3B,
                             RNA_DNA_2_2A, RNA_DNA_2_2B, RNA_DNA_2_1A, RNA_DNA_2_1B,
                             RNA_DNA_20A, RNA_DNA_20B, RNA_DNA_22A, RNA_DNA_22B)


#determine the log(RNA/DNA) for each sample (this takes the log of sum_RNA and sum_DNA 
#as well). Only perform this function if 0 reads have been replaced with 1

var_log <- function(df) {
  log_ratio_df <- df %>% 
    mutate_if(is.double, 
              funs(log10(.))
    )
  return(log_ratio_df)
}

log_rep_0_22_A_B <- var_log(rep_0_22_A_B)


#BC analysis---------------------------------------------------------------------------------

#Reads/BC of BC's with > 1 read in DNA, can have 0 reads in RNA

p_BC_num_reads_viol_full <- ggplot(rep_0_22_A_B, aes(x = "", y = NULL)) +
  geom_violin(aes(x = "DNA", y = num_reads_DNA, color = subpool)) +
  geom_violin(aes(x = "R0.00A", y = num_reads_RNA_0A, color = subpool)) + 
  geom_violin(aes(x = "R0.00B", y = num_reads_RNA_0B, color = subpool)) + 
  geom_violin(aes(x = "R0.03A", y = num_reads_RNA_2_5A, color = subpool)) + 
  geom_violin(aes(x = "R0.03B", y = num_reads_RNA_2_5B, color = subpool)) + 
  geom_violin(aes(x = "R0.06A", y = num_reads_RNA_2_4A, color = subpool)) +
  geom_violin(aes(x = "R0.06B", y = num_reads_RNA_2_4B, color = subpool)) + 
  geom_violin(aes(x = "R0.13A", y = num_reads_RNA_2_3A, color = subpool)) +
  geom_violin(aes(x = "R0.13B", y = num_reads_RNA_2_3B, color = subpool)) + 
  geom_violin(aes(x = "R0.25A", y = num_reads_RNA_2_2A, color = subpool)) +
  geom_violin(aes(x = "R0.25B", y = num_reads_RNA_2_2B, color = subpool)) + 
  geom_violin(aes(x = "R0.50A", y = num_reads_RNA_2_1A, color = subpool)) +
  geom_violin(aes(x = "R0.50B", y = num_reads_RNA_2_1B, color = subpool)) + 
  geom_violin(aes(x = "R1.00A", y = num_reads_RNA_20A, color = subpool)) +
  geom_violin(aes(x = "R1.00B", y = num_reads_RNA_20B, color = subpool)) + 
  geom_violin(aes(x = "R4.00A", y = num_reads_RNA_22A, color = subpool)) +
  geom_violin(aes(x = "R4.00B", y = num_reads_RNA_22B, color = subpool)) + 
  xlab("") +
  ylab("Reads per BC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot('plots/BC_num_reads_viol_full.png',
          p_BC_num_reads_viol_full, scale = 2.8)

p_BC_num_reads_box_zoom <- ggplot(rep_0_22_A_B, aes(x = "", y = NULL)) +
  geom_boxplot(aes(x = "DNA", y = num_reads_DNA, color = subpool)) +
  geom_boxplot(aes(x = "R0.00A", y = num_reads_RNA_0A, color = subpool)) + 
  geom_boxplot(aes(x = "R0.00B", y = num_reads_RNA_0B, color = subpool)) + 
  geom_boxplot(aes(x = "R0.03A", y = num_reads_RNA_2_5A, color = subpool)) + 
  geom_boxplot(aes(x = "R0.03B", y = num_reads_RNA_2_5B, color = subpool)) + 
  geom_boxplot(aes(x = "R0.06A", y = num_reads_RNA_2_4A, color = subpool)) +
  geom_boxplot(aes(x = "R0.06B", y = num_reads_RNA_2_4B, color = subpool)) + 
  geom_boxplot(aes(x = "R0.13A", y = num_reads_RNA_2_3A, color = subpool)) +
  geom_boxplot(aes(x = "R0.13B", y = num_reads_RNA_2_3B, color = subpool)) + 
  geom_boxplot(aes(x = "R0.25A", y = num_reads_RNA_2_2A, color = subpool)) +
  geom_boxplot(aes(x = "R0.25B", y = num_reads_RNA_2_2B, color = subpool)) + 
  geom_boxplot(aes(x = "R0.50A", y = num_reads_RNA_2_1A, color = subpool)) +
  geom_boxplot(aes(x = "R0.50B", y = num_reads_RNA_2_1B, color = subpool)) + 
  geom_boxplot(aes(x = "R1.00A", y = num_reads_RNA_20A, color = subpool)) +
  geom_boxplot(aes(x = "R1.00B", y = num_reads_RNA_20B, color = subpool)) + 
  geom_boxplot(aes(x = "R4.00A", y = num_reads_RNA_22A, color = subpool)) +
  geom_boxplot(aes(x = "R4.00B", y = num_reads_RNA_22B, color = subpool)) + 
  xlab("") +
  ylab("Reads per BC") + ylim(0, 25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot('plots/BC_num_reads_box_zoom.png',
          p_BC_num_reads_box_zoom, scale = 2.8)


#replicate plots

#plot replicates for BC's similarly

p_bc_rep_0 <- ggplot(NULL, aes(ratio_0A, ratio_0B)) +
  geom_point(data = log_rep_0_22_A_B, alpha = 0.1) +
  geom_point(data = filter(log_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red', alpha = 0.3) +
  annotation_logticks(scaled = TRUE) +
  xlab("Log10 BC RNA/DNA Rep. 1") +
  ylab("Log10 BC RNA/DNA Rep. 2") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.5, 2.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.5, 2.5)) + 
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

p_bc_rep_2_5 <- ggplot(NULL, aes(ratio_2_5A, ratio_2_5B)) +
  geom_point(data = log_rep_0_22_A_B, alpha = 0.1) +
  geom_point(data = filter(log_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red', alpha = 0.3) +
  annotation_logticks(scaled = TRUE) +
  xlab("Log10 BC RNA/DNA Rep. 1") +
  ylab("Log10 BC RNA/DNA Rep. 2") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.5, 2.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.5, 2.5)) + 
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

p_bc_rep_2_4 <- ggplot(NULL, aes(ratio_2_4A, ratio_2_4B)) +
  geom_point(data = log_rep_0_22_A_B, alpha = 0.1) +
  geom_point(data = filter(log_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red', alpha = 0.3) +
  annotation_logticks(scaled = TRUE) +
  xlab("Log10 BC RNA/DNA Rep. 1") +
  ylab("Log10 BC RNA/DNA Rep. 2") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.5, 2.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.5, 2.5)) + 
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

p_bc_rep_2_3 <- ggplot(NULL, aes(ratio_2_3A, ratio_2_3B)) +
  geom_point(data = log_rep_0_22_A_B, alpha = 0.1) +
  geom_point(data = filter(log_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red', alpha = 0.3) +
  annotation_logticks(scaled = TRUE) +
  xlab("Log10 BC RNA/DNA Rep. 1") +
  ylab("Log10 BC RNA/DNA Rep. 2") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.5, 2.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.5, 2.5)) + 
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

p_bc_rep_2_2 <- ggplot(NULL, aes(ratio_2_2A, ratio_2_2B)) +
  geom_point(data = log_rep_0_22_A_B, alpha = 0.1) +
  geom_point(data = filter(log_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red', alpha = 0.3) +
  annotation_logticks(scaled = TRUE) +
  xlab("Log10 BC RNA/DNA Rep. 1") +
  ylab("Log10 BC RNA/DNA Rep. 2") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.5, 2.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.5, 2.5)) + 
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

p_bc_rep_2_1 <- ggplot(NULL, aes(ratio_2_1A, ratio_2_1B)) +
  geom_point(data = log_rep_0_22_A_B, alpha = 0.1) +
  geom_point(data = filter(log_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red', alpha = 0.3) +
  annotation_logticks(scaled = TRUE) +
  xlab("Log10 BC RNA/DNA Rep. 1") +
  ylab("Log10 BC RNA/DNA Rep. 2") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.5, 2.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.5, 2.5)) + 
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

p_bc_rep_20 <- ggplot(NULL, aes(ratio_20A, ratio_20B)) +
  geom_point(data = log_rep_0_22_A_B, alpha = 0.1) +
  geom_point(data = filter(log_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red', alpha = 0.3) +
  annotation_logticks(scaled = TRUE) +
  xlab("Log10 BC RNA/DNA Rep. 1") +
  ylab("Log10 BC RNA/DNA Rep. 2") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.5, 2.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.5, 2.5)) + 
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

p_bc_rep_22 <- ggplot(NULL, aes(ratio_22A, ratio_22B)) +
  geom_point(data = log_rep_0_22_A_B, alpha = 0.1) +
  geom_point(data = filter(log_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red', alpha = 0.3) +
  annotation_logticks(scaled = TRUE) +
  xlab("Log10 BC RNA/DNA Rep. 1") +
  ylab("Log10 BC RNA/DNA Rep. 2") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.5, 2.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.5, 2.5)) + 
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

p_bc_rep_grid <- plot_grid(
  p_bc_rep_0, p_bc_rep_2_5, p_bc_rep_2_4, p_bc_rep_2_3, p_bc_rep_2_2, p_bc_rep_2_1, 
  p_bc_rep_20, p_bc_rep_22, 
  labels = c(
    "0 µM", "2^-5 µM", "2^-4 µM","2^-3 µM", "2^-2 µM", "2^-1 µM","2^0 µM", "2^2 µM"),
  nrow = 3, ncol = 3, align = 'hv', hjust = -3, vjust = 0.5, scale = 0.9)

save_plot('plots/p_bc_rep_grid.png', 
          p_bc_rep_grid, base_height = 10, base_width = 10)




