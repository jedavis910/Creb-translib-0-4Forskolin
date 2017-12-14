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
variant_counts_DNA <- filter(variant_counts_DNA, barcodes > 7)

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


#determine the log(RNA/DNA) for each sample (this takes the log of sum_RNA and sum_DNA 
#as well) This is useful for replicate plots, but manipulations should use rep_0_22_A_B

var_log2 <- function(df) {
  log_ratio_df <- df %>% 
    mutate_if(is.double, 
              funs(log2(.))
              )
  return(log_ratio_df)
}

var_log10 <- function(df) {
  log_ratio_df <- df %>% 
    mutate_if(is.double, 
              funs(log2(.))
    )
  return(log_ratio_df)
}

log2_rep_0_22_A_B <- var_log2(rep_0_22_A_B)


#Separate into subpools----------------------------------------------------------------------

#Perform subpool analyses without log2 data, figure out better way to condense operations







#Subpool 3 contains 2 consensus binding sites with flanks (ATTGACGTCAGC) that vary in 
#distance from one another by 0 (no inner flanks), 5, 10, 15, 20 and 70 bp (all but 0 appear 
#as -4 bp spacing). Each site distance combination is then moved along the backgrounds at 1 
#bp increments starting from closest to the minP. Separation lists the spacing between sites,
#distance (start of consensus and flanks) and the background. Added 2 to all distances to 
#measure to start of BS and not to flank. Added 4 to all spacing but 0 to measure difference 
#between start of sites. Also took average of log2 med BC expression between biological 
#replicates for plotting

subpool3 <- 
  filter(log2_rep_0_22_A_B, subpool == "subpool3") %>%
  ungroup() %>%
  select(-subpool) %>%
  mutate(name = gsub('2BS ', '', name), 
         name = gsub(' bp spacing ', '_', name)) %>%
  separate(name, 
           into = c("subpool", "spacing", "fluff2", "fluff3", "dist", "background"),
           sep = "_", convert = TRUE
  ) %>%
  select(-subpool, -fluff2, -fluff3) %>%
  mutate(background = gsub('Smith R. Vista chr9:83712599-83712766',
                           'vista chr9', background),
         background = gsub('Vista Chr5:88673410-88674494', 'vista chr5',
                           background),
         background = str_sub(background, 1, 13)
  ) %>%
  mutate(dist = dist + 2) %>%
  mutate(spacing = 
           ifelse(spacing != as.integer(0), 
                  spacing + 4, spacing)) %>%
  mutate(ave_ratio_0 = (ratio_0A + ratio_0B)/2) %>%
  mutate(ave_ratio_2_5 = (ratio_2_5A + ratio_2_5B)/2) %>%
  mutate(ave_ratio_2_4 = (ratio_2_4A + ratio_2_4B)/2) %>%
  mutate(ave_ratio_2_3 = (ratio_2_3A + ratio_2_3B)/2) %>%
  mutate(ave_ratio_2_2 = (ratio_2_2A + ratio_2_2B)/2) %>%
  mutate(ave_ratio_2_1 = (ratio_2_1A + ratio_2_1B)/2) %>%
  mutate(ave_ratio_20 = (ratio_20A + ratio_20B)/2) %>%
  mutate(ave_ratio_22 = (ratio_22A + ratio_22B)/2)

#Subpool 5 contains 6 equally spaced sites spaced 13 bp apart and starting from furthest to 
#the minP. These sites are filled with sites of either the consensus site, a weak site or no 
#site. Both the weak and consensus sites are flanked by the same flanking sequence. Also took
#average of log2 med BC expression between biological replicates for plotting

subpool5 <- 
  filter(log2_rep_0_22_A_B, subpool == "subpool5") %>%
  ungroup() %>%
  select(-subpool) %>%
  mutate(name = gsub('no_site', 'nosite', name)) %>%
  separate(name, into = c(
    "subpool", "site1", "site2", "site3", "site4", "site5", "site6", 
    "background"), sep = "_"
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
  mutate(total_sites = consensus + weak) %>%
  mutate(site_combo = 
           ifelse(weak == 0 & consensus > 0, 
                  'consensus', 'mixed')) %>%
  mutate(site_type = 
           ifelse(consensus == 0 & weak > 0, 
                  'weak', site_combo)) %>%
  mutate(background = gsub('Smith R. Vista chr9:83712599-83712766',
                           'vista chr9', background),
         background = gsub('Vista Chr5:88673410-88674494', 'vista chr5',
                           background),
         background = str_sub(background, 1, 13)) %>%
  mutate(ave_ratio_0 = (ratio_0A + ratio_0B)/2) %>%
  mutate(ave_ratio_2_5 = (ratio_2_5A + ratio_2_5B)/2) %>%
  mutate(ave_ratio_2_4 = (ratio_2_4A + ratio_2_4B)/2) %>%
  mutate(ave_ratio_2_3 = (ratio_2_3A + ratio_2_3B)/2) %>%
  mutate(ave_ratio_2_2 = (ratio_2_2A + ratio_2_2B)/2) %>%
  mutate(ave_ratio_2_1 = (ratio_2_1A + ratio_2_1B)/2) %>%
  mutate(ave_ratio_20 = (ratio_20A + ratio_20B)/2) %>%
  mutate(ave_ratio_22 = (ratio_22A + ratio_22B)/2)

controls <- 
  filter(log2_rep_0_22_A_B, subpool == "control") %>%
  ungroup()

#Normalize each reads within each subpool to background

backgrounds <- subpool5 %>%
  filter(total_sites == 0) %>%
  select(background, ave_ratio_0, ave_ratio_2_5, ave_ratio_2_4, ave_ratio_2_3, 
         ave_ratio_2_2, ave_ratio_2_1, ave_ratio_20, ave_ratio_22) %>%
  rename(ave_ratio_0_back = ave_ratio_0) %>%
  rename(ave_ratio_2_5_back = ave_ratio_2_5) %>%
  rename(ave_ratio_2_4_back = ave_ratio_2_4) %>%
  rename(ave_ratio_2_3_back = ave_ratio_2_3) %>%
  rename(ave_ratio_2_2_back = ave_ratio_2_2) %>%
  rename(ave_ratio_2_1_back = ave_ratio_2_1) %>%
  rename(ave_ratio_20_back = ave_ratio_20) %>%
  rename(ave_ratio_22_back = ave_ratio_22)

subpool3_norm <- left_join(subpool3, backgrounds, by = 'background') %>%
  mutate(ave_ratio_0_norm = ave_ratio_0/ave_ratio_0_back) %>%
  mutate(ave_ratio_2_5_norm = ave_ratio_2_5/ave_ratio_2_5_back) %>%
  mutate(ave_ratio_2_4_norm = ave_ratio_2_4/ave_ratio_2_4_back) %>%
  mutate(ave_ratio_2_3_norm = ave_ratio_2_3/ave_ratio_2_3_back) %>%
  mutate(ave_ratio_2_2_norm = ave_ratio_2_2/ave_ratio_2_2_back) %>%
  mutate(ave_ratio_2_1_norm = ave_ratio_2_1/ave_ratio_2_1_back) %>%
  mutate(ave_ratio_20_norm = ave_ratio_20/ave_ratio_20_back) %>%
  mutate(ave_ratio_22_norm = ave_ratio_22/ave_ratio_22_back)
  

#Plot subpool expression features-----------------------------------------------------------

#Subpool 3

p_subpool3_spa_back_norm <- ggplot(subpool3_norm, aes(x = dist)) + 
  geom_point(aes(y = ave_ratio_0_norm), alpha = 0.5, size = 1.5, color = '#440154FF') +
  geom_point(aes(y = ave_ratio_2_5_norm), alpha = 0.5, size = 1.5, color = '#482677FF') +
  geom_point(aes(y = ave_ratio_2_4_norm), alpha = 0.5, size = 1.5, color = '#39568CFF') +
  geom_point(aes(y = ave_ratio_2_3_norm), alpha = 0.5, size = 1.5, color = '#2D708EFF') +
  geom_point(aes(y = ave_ratio_2_2_norm), alpha = 0.5, size = 1.5, color = '#1F968BFF') +
  geom_point(aes(y = ave_ratio_2_1_norm), alpha = 0.5, size = 1.5, color = '#29AF7FFF') +
  geom_point(aes(y = ave_ratio_20_norm), alpha = 0.5, size = 1.5, color = '#73D055FF') +
  geom_point(aes(y = ave_ratio_22_norm), alpha = 0.5, size = 1.5, color = '#B8DE29FF') +
  facet_grid(spacing ~ background) + 
  ylab('Normalized average log2 sum BC expression') + 
  panel_border() +
  background_grid(major = 'xy', minor = 'none') +
  scale_x_continuous(
    "Distance from First Site to Proximal Promoter End (bp)", 
    breaks = seq(from = 0, to = 150, by = 10))

save_plot('plots/p_subpool3_spa_back_norm.png', p_subpool3_spa_back_norm, 
          base_width = 46, base_height = 17, scale = 0.35)


#BC analysis---------------------------------------------------------------------------------

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


#replicate plots----------------------------------------------------------------------------

#plot replicates for summed variant expression

p_var_rep_0 <- ggplot(NULL, aes(ratio_0A, ratio_0B)) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = log2_rep_0_22_A_B, alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.5) +
  geom_point(data = filter(log2_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log2(sum(RNA/DNA)) TR 1") +
  ylab("Variant log2(sum(RNA/DNA)) TR 2") +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4.5)) + 
  scale_y_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4.5)) + 
  annotate("text", x = -2, y = 2,
           label = paste(
             'r =', round(
               cor(
                 log2_rep_0_22_A_B$ratio_0A,log2_rep_0_22_A_B$ratio_0B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_rep_2_5 <- ggplot(NULL, aes(ratio_2_5A, ratio_2_5B)) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = log2_rep_0_22_A_B, alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.5) +
  geom_point(data = filter(log2_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log2(sum(RNA/DNA)) TR 1") +
  ylab("Variant log2(sum(RNA/DNA)) TR 2") +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4.5)) + 
  scale_y_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4.5)) + 
  annotate("text", x = -2, y = 2,
           label = paste(
             'r =', round(
               cor(
                 log2_rep_0_22_A_B$ratio_2_5A,log2_rep_0_22_A_B$ratio_2_5B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_rep_2_4 <- ggplot(NULL, aes(ratio_2_4A, ratio_2_4B)) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = log2_rep_0_22_A_B, alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.5) +
  geom_point(data = filter(log2_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log2(sum(RNA/DNA)) TR 1") +
  ylab("Variant log2(sum(RNA/DNA)) TR 2") +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4.5)) + 
  scale_y_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4.5)) + 
  annotate("text", x = -2, y = 2,
           label = paste(
             'r =', round(
               cor(
                 log2_rep_0_22_A_B$ratio_2_4A,log2_rep_0_22_A_B$ratio_2_4B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_rep_2_3 <- ggplot(NULL, aes(ratio_2_3A, ratio_2_3B)) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = log2_rep_0_22_A_B, alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.5) +
  geom_point(data = filter(log2_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log2(sum(RNA/DNA)) TR 1") +
  ylab("Variant log2(sum(RNA/DNA)) TR 2") +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4.5)) + 
  scale_y_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4.5)) + 
  annotate("text", x = -2, y = 2,
           label = paste(
             'r =', round(
               cor(
                 log2_rep_0_22_A_B$ratio_2_3A,log2_rep_0_22_A_B$ratio_2_3B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_rep_2_2 <- ggplot(NULL, aes(ratio_2_2A, ratio_2_2B)) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = log2_rep_0_22_A_B, alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.5) +
  geom_point(data = filter(log2_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log2(sum(RNA/DNA)) TR 1") +
  ylab("Variant log2(sum(RNA/DNA)) TR 2") +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4.5)) + 
  scale_y_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4.5)) + 
  annotate("text", x = -2, y = 2,
           label = paste(
             'r =', round(
               cor(
                 log2_rep_0_22_A_B$ratio_2_2A,log2_rep_0_22_A_B$ratio_2_2B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_rep_2_1 <- ggplot(NULL, aes(ratio_2_1A, ratio_2_1B)) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = log2_rep_0_22_A_B, alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.5) +
  geom_point(data = filter(log2_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log2(sum(RNA/DNA)) TR 1") +
  ylab("Variant log2(sum(RNA/DNA)) TR 2") +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4.5)) + 
  scale_y_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4.5)) + 
  annotate("text", x = -2, y = 2,
           label = paste(
             'r =', round(
               cor(
                 log2_rep_0_22_A_B$ratio_2_1A,log2_rep_0_22_A_B$ratio_2_1B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_rep_20 <- ggplot(NULL, aes(ratio_20A, ratio_20B)) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = log2_rep_0_22_A_B, alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.5) +
  geom_point(data = filter(log2_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log2(sum(RNA/DNA)) TR 1") +
  ylab("Variant log2(sum(RNA/DNA)) TR 2") +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4.5)) + 
  scale_y_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4.5)) + 
  annotate("text", x = -2, y = 2,
           label = paste(
             'r =', round(
               cor(
                 log2_rep_0_22_A_B$ratio_20A,log2_rep_0_22_A_B$ratio_20B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_rep_22 <- ggplot(NULL, aes(ratio_22A, ratio_22B)) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = log2_rep_0_22_A_B, alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log2_rep_0_22_A_B, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.5) +
  geom_point(data = filter(log2_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log2(sum(RNA/DNA)) TR 1") +
  ylab("Variant log2(sum(RNA/DNA)) TR 2") +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4.5)) + 
  scale_y_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4.5)) + 
  annotate("text", x = -2, y = 2,
           label = paste(
             'r =', round(
               cor(
                 log2_rep_0_22_A_B$ratio_22A,log2_rep_0_22_A_B$ratio_22B,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_log2_var_rep_grid <- plot_grid(
  p_var_rep_0, p_var_rep_2_5, p_var_rep_2_4, p_var_rep_2_3, p_var_rep_2_2, p_var_rep_2_1, 
  p_var_rep_20, p_var_rep_22, 
  labels = c(
    "    0 µM", "2^-5 µM", "2^-4 µM","2^-3 µM", "2^-2 µM", "2^-1 µM"," 2^0 µM", " 2^2 µM"),
  nrow = 3, ncol = 3, align = 'hv', hjust = -2.5, vjust = 0.5, scale = 0.9)

save_plot('plots/p_log2_var_rep_grid.png', 
          p_log2_var_rep_grid, base_height = 12, base_width = 12)


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









#Untidy joining---------------------------------------------------------------------------

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

