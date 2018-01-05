library(tidyverse)
library(stringr)
library(viridis)
library(cowplot)
library(ggExtra)
library(modelr)
library(lazyeval)
library(splines)

#Written for the analysis of a range of inductions in the transient library
#D3_BC: DNA
#R#A_BC: RNA at # µM Forsk replicate A
#R#B_BC: RNA at # µM Forsk replicate B

#tested concentrations: 0, 1, 2, 4, 8, 16, 25, 32, 64


#Load index and bcmap files------------------------------------------------------------------

bc_DNA <- read_tsv('D3_BC.txt')
bc_0A <- read_tsv('R0A_BC.txt')
bc_0B <- read_tsv('R0B_BC.txt')
bc_1A <- read_tsv('R1A_BC.txt')
bc_1B <- read_tsv('R1B_BC.txt')
bc_2A <- read_tsv('R2A_BC.txt')
bc_2B <- read_tsv('R2B_BC.txt')
bc_4A <- read_tsv('R4A_BC.txt')
bc_4B <- read_tsv('R4B_BC.txt')
bc_8A <- read_tsv('R8A_BC.txt')
bc_8B <- read_tsv('R8B_BC.txt')
bc_16A <- read_tsv('R16A_BC.txt')
bc_16B <- read_tsv('R16B_BC.txt')
bc_25A <- read_tsv('R25A_BC.txt')
bc_25B <- read_tsv('R25B_BC.txt')
bc_32A <- read_tsv('R32A_BC.txt')
bc_32B <- read_tsv('R32B_BC.txt')
bc_64A <- read_tsv('R64A_BC.txt')
bc_64B <- read_tsv('R64B_BC.txt')


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
bc_join_0A <- bc_map_join_bc(SP3_SP5_map, bc_0A)
bc_join_0B <- bc_map_join_bc(SP3_SP5_map, bc_0B)
bc_join_1A <- bc_map_join_bc(SP3_SP5_map, bc_1A)
bc_join_1B <- bc_map_join_bc(SP3_SP5_map, bc_1B)
bc_join_2A <- bc_map_join_bc(SP3_SP5_map, bc_2A)
bc_join_2B <- bc_map_join_bc(SP3_SP5_map, bc_2B)
bc_join_4A <- bc_map_join_bc(SP3_SP5_map, bc_4A)
bc_join_4B <- bc_map_join_bc(SP3_SP5_map, bc_4B)
bc_join_8A <- bc_map_join_bc(SP3_SP5_map, bc_8A)
bc_join_8B <- bc_map_join_bc(SP3_SP5_map, bc_8B)
bc_join_16A <- bc_map_join_bc(SP3_SP5_map, bc_16A)
bc_join_16B <- bc_map_join_bc(SP3_SP5_map, bc_16B)
bc_join_25A <- bc_map_join_bc(SP3_SP5_map, bc_25A)
bc_join_25B <- bc_map_join_bc(SP3_SP5_map, bc_25B)
bc_join_32A <- bc_map_join_bc(SP3_SP5_map, bc_32A)
bc_join_32B <- bc_map_join_bc(SP3_SP5_map, bc_32B)
bc_join_64A <- bc_map_join_bc(SP3_SP5_map, bc_64A)
bc_join_64B <- bc_map_join_bc(SP3_SP5_map, bc_64B)


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

variant_counts_0A <- var_sum_bc_num(bc_join_0A)
variant_counts_0B <- var_sum_bc_num(bc_join_0B)
variant_counts_1A <- var_sum_bc_num(bc_join_1A)
variant_counts_1B <- var_sum_bc_num(bc_join_1B)
variant_counts_2A <- var_sum_bc_num(bc_join_2A)
variant_counts_2B <- var_sum_bc_num(bc_join_2B)
variant_counts_4A <- var_sum_bc_num(bc_join_4A)
variant_counts_4B <- var_sum_bc_num(bc_join_4B)
variant_counts_8A <- var_sum_bc_num(bc_join_8A)
variant_counts_8B <- var_sum_bc_num(bc_join_8B)
variant_counts_16A <- var_sum_bc_num(bc_join_16A)
variant_counts_16B <- var_sum_bc_num(bc_join_16B)
variant_counts_25A <- var_sum_bc_num(bc_join_25A)
variant_counts_25B <- var_sum_bc_num(bc_join_25B)
variant_counts_32A <- var_sum_bc_num(bc_join_32A)
variant_counts_32B <- var_sum_bc_num(bc_join_32B)
variant_counts_64A <- var_sum_bc_num(bc_join_64A)
variant_counts_64B <- var_sum_bc_num(bc_join_64B)


#Join RNA to DNA and determine expression from summing--------------------------------------

#combine DNA and RNA cumm. BC counts, only keeping instances in both sets (as only 1 control 
#drops out in some samples) and determining RNA/DNA per variant. Ratio is summed normalized
#reads of RNA over DNA

var_expression <- function(df1, df2) {
  RNA_DNA <- inner_join(df1, df2, 
                        by = c("name", "subpool", "most_common"), 
                        suffix = c("_DNA", "_RNA")
  ) %>%
    mutate(ratio = sum_RNA / sum_DNA)
  print('x defined as DNA, y defined as RNA in var_expression(x,y)')
  return(RNA_DNA)
}

RNA_DNA_0A <- var_expression(variant_counts_DNA, variant_counts_0A)
RNA_DNA_0B <- var_expression(variant_counts_DNA, variant_counts_0B)
RNA_DNA_1A <- var_expression(variant_counts_DNA, variant_counts_1A)
RNA_DNA_1B <- var_expression(variant_counts_DNA, variant_counts_1B)
RNA_DNA_2A <- var_expression(variant_counts_DNA, variant_counts_2A)
RNA_DNA_2B <- var_expression(variant_counts_DNA, variant_counts_2B)
RNA_DNA_4A <- var_expression(variant_counts_DNA, variant_counts_4A)
RNA_DNA_4B <- var_expression(variant_counts_DNA, variant_counts_4B)
RNA_DNA_8A <- var_expression(variant_counts_DNA, variant_counts_8A)
RNA_DNA_8B <- var_expression(variant_counts_DNA, variant_counts_8B)
RNA_DNA_16A <- var_expression(variant_counts_DNA, variant_counts_16A)
RNA_DNA_16B <- var_expression(variant_counts_DNA, variant_counts_16B)
RNA_DNA_25A <- var_expression(variant_counts_DNA, variant_counts_25A)
RNA_DNA_25B <- var_expression(variant_counts_DNA, variant_counts_25B)
RNA_DNA_32A <- var_expression(variant_counts_DNA, variant_counts_32A)
RNA_DNA_32B <- var_expression(variant_counts_DNA, variant_counts_32B)
RNA_DNA_64A <- var_expression(variant_counts_DNA, variant_counts_64A)
RNA_DNA_64B <- var_expression(variant_counts_DNA, variant_counts_64B)


#combine biological replicates---------------------------------------------------------------

#After combining, rename backgrounds to simplified names, make background column (excluding 
#controls), separate out background values in each dataset and left join to original dataset.
#Normalize expression of each variant to its background in that biological replicate.
#Determine average normalized expression across biological replicates.

var_conc_rep <- function(
  df0A, df0B, df1A, df1B, df2A, df2B, df4A, df4B, df8A, df8B, df16A, df16B, 
  df25A, df25B, df32A, df32B, df64A, df64B
) {
  join_0 <- inner_join(df0A, df0B, 
                       by = c(
                         "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                       ), suffix = c("_0A", "_0B")
  )
  join_1 <- inner_join(df1A, df1B, 
                         by = c(
                           "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                         ), suffix = c("_1A", "_1B")
  )
  join_2 <- inner_join(df2A, df2B, 
                         by = c(
                           "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                         ), suffix = c("_2A", "_2B")
  )
  join_4 <- inner_join(df4A, df4B, 
                         by = c(
                           "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                         ), suffix = c("_4A", "_4B")
  )
  join_8 <- inner_join(df8A, df8B, 
                         by = c(
                           "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                         ), suffix = c("_8A", "_8B")
  )
  join_16 <- inner_join(df16A, df16B, 
                         by = c(
                           "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                         ), suffix = c("_16A", "_16B")
  )
  join_25 <- inner_join(df25A, df25B, 
                        by = c(
                          "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                        ), suffix = c("_25A", "_25B")
  )
  join_32 <- inner_join(df32A, df32B, 
                        by = c(
                          "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                        ), suffix = c("_32A", "_32B")
  )
  join_64 <- inner_join(df64A, df64B, 
                        by = c(
                          "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                        ), suffix = c("_64A", "_64B")
  )
  join_0_1 <- inner_join(join_0, join_1, 
                           by = c(
                             "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                           )
  )
  join_0_2 <- inner_join(join_0_1, join_2, 
                           by = c(
                             "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                           )
  )
  join_0_4 <- inner_join(join_0_2, join_4, 
                           by = c(
                             "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                           )
  )
  join_0_8 <- inner_join(join_0_4, join_8, 
                           by = c(
                             "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                           )
  )
  join_0_16 <- inner_join(join_0_8, join_16, 
                           by = c(
                             "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                           )
  )
  join_0_25 <- inner_join(join_0_16, join_25, 
                          by = c(
                            "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                          )
  )
  join_0_32 <- inner_join(join_0_25, join_32, 
                          by = c(
                            "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                          )
  )
  join_0_64 <- inner_join(join_0_32, join_64, 
                          by = c(
                            "name", "subpool", "most_common", "sum_DNA", "barcodes_DNA"
                          )
  )
  print(
    'processed dfs in order of samples: 0A, 0B, 1A, 1B, 2A, 2B, 4A, 4B, 8A, 
    8B, 16A, 16B, 25A, 25B, 32A, 32B, 64A, 64B'
  )
  return(join_0_64)
}

rep_0_64_A_B <- var_conc_rep(RNA_DNA_0A, RNA_DNA_0B, RNA_DNA_1A, RNA_DNA_1B,
                             RNA_DNA_2A, RNA_DNA_2B, RNA_DNA_4A, RNA_DNA_4B,
                             RNA_DNA_8A, RNA_DNA_8B, RNA_DNA_16A, RNA_DNA_16B,
                             RNA_DNA_25A, RNA_DNA_25B, RNA_DNA_32A, RNA_DNA_32B,
                             RNA_DNA_64A, RNA_DNA_64B)

#Normalize to background--------------------------------------------------------------------

#Ideally would have put this earlier as there are so many columns to modify, but normalizing
#to backgrounds excludes the controls, so I would have to repeat the code for the controls...

back_norm <- function(df1) {
  gsub_0_64 <- df1 %>%
    ungroup () %>%
    filter(subpool != 'control') %>%
    mutate(
      name = gsub('Smith R. Vista chr9:83712599-83712766', 'v chr9', name),
      name = gsub('Vista Chr5:88673410-88674494', 'v chr5', name),
      name = gsub('scramble pGL4.29 Promega 1-63 \\+ 1-87', 's pGl4', name)
    ) %>%
    mutate(background = name) %>%
    mutate(background = str_sub(background, nchar(background)-5, nchar(background)))
  backgrounds <- gsub_0_64 %>%
    filter(startsWith(name, 'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')) %>%
    select(background, ratio_0A, ratio_0B, ratio_1A, ratio_1B, ratio_2A, ratio_2B,
           ratio_4A, ratio_4B, ratio_8A, ratio_8B, ratio_16A, ratio_16B, ratio_25A,
           ratio_25B, ratio_32A, ratio_32B, ratio_64A, ratio_64B) %>%
    rename(ratio_0A_back = ratio_0A) %>%
    rename(ratio_0B_back = ratio_0B) %>%
    rename(ratio_1A_back = ratio_1A) %>%
    rename(ratio_1B_back = ratio_1B) %>%
    rename(ratio_2A_back = ratio_2A) %>%
    rename(ratio_2B_back = ratio_2B) %>%
    rename(ratio_4A_back = ratio_4A) %>%
    rename(ratio_4B_back = ratio_4B) %>%
    rename(ratio_8A_back = ratio_8A) %>%
    rename(ratio_8B_back = ratio_8B) %>%
    rename(ratio_16A_back = ratio_16A) %>%
    rename(ratio_16B_back = ratio_16B) %>%
    rename(ratio_25A_back = ratio_25A) %>%
    rename(ratio_25B_back = ratio_25B) %>%
    rename(ratio_32A_back = ratio_32A) %>%
    rename(ratio_32B_back = ratio_32B) %>%
    rename(ratio_64A_back = ratio_64A) %>%
    rename(ratio_64B_back = ratio_64B)
  back_join_norm <- left_join(gsub_0_64, backgrounds, by = 'background') %>%
    mutate(ratio_0A_norm = ratio_0A/ratio_0A_back) %>%
    mutate(ratio_0B_norm = ratio_0B/ratio_0B_back) %>%
    mutate(ratio_1A_norm = ratio_1A/ratio_1A_back) %>%
    mutate(ratio_1B_norm = ratio_1B/ratio_1B_back) %>%
    mutate(ratio_2A_norm = ratio_2A/ratio_2A_back) %>%
    mutate(ratio_2B_norm = ratio_2B/ratio_2B_back) %>%
    mutate(ratio_4A_norm = ratio_4A/ratio_4A_back) %>%
    mutate(ratio_4B_norm = ratio_4B/ratio_4B_back) %>%
    mutate(ratio_8A_norm = ratio_8A/ratio_8A_back) %>%
    mutate(ratio_8B_norm = ratio_8B/ratio_8B_back) %>%
    mutate(ratio_16A_norm = ratio_16A/ratio_16A_back) %>%
    mutate(ratio_16B_norm = ratio_16B/ratio_16B_back) %>%
    mutate(ratio_25A_norm = ratio_25A/ratio_25A_back) %>%
    mutate(ratio_25B_norm = ratio_25B/ratio_25B_back) %>%
    mutate(ratio_32A_norm = ratio_32A/ratio_32A_back) %>%
    mutate(ratio_32B_norm = ratio_32B/ratio_32B_back) %>%
    mutate(ratio_64A_norm = ratio_64A/ratio_64A_back) %>%
    mutate(ratio_64B_norm = ratio_64B/ratio_64B_back) %>%
    mutate(ave_ratio_0_norm = (ratio_0A_norm + ratio_0B_norm)/2) %>%
    mutate(ave_ratio_1_norm = (ratio_1A_norm + ratio_1B_norm)/2) %>%
    mutate(ave_ratio_2_norm = (ratio_2A_norm + ratio_2B_norm)/2) %>%
    mutate(ave_ratio_4_norm = (ratio_4A_norm + ratio_4B_norm)/2) %>%
    mutate(ave_ratio_8_norm = (ratio_8A_norm + ratio_8B_norm)/2) %>%
    mutate(ave_ratio_16_norm = (ratio_16A_norm + ratio_16B_norm)/2) %>%
    mutate(ave_ratio_25_norm = (ratio_25A_norm + ratio_25B_norm)/2) %>%
    mutate(ave_ratio_32_norm = (ratio_32A_norm + ratio_32B_norm)/2) %>%
    mutate(ave_ratio_64_norm = (ratio_64A_norm + ratio_64B_norm)/2)
}

trans_back_norm_rep_0_64 <- back_norm(rep_0_64_A_B)


#Write table to compare expression to other analyses-----------------------------------------

output_int <- trans_back_norm_rep_0_64 %>%
  write.table(
    "trans_back_norm_rep_0_64.txt", 
    sep = '\t', row.names = FALSE)


