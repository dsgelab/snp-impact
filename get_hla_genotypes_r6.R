setwd("/home/jsjukara/hla/files")
library(data.table)
library(tidyverse)

### read different HLA files to a list of data frames
df_list <- list(NULL)
df_list[[1]] <- read_tsv("/home/jsjukara/hla/allele_dosages/R6_HLA_A_imputed_allele_dosages.tsv")
df_list[[2]] <- read_tsv("/home/jsjukara/hla/allele_dosages/R6_HLA_B_imputed_allele_dosages.tsv")
df_list[[3]] <- read_tsv("/home/jsjukara/hla/allele_dosages/R6_HLA_C_imputed_allele_dosages.tsv")
df_list[[4]] <- read_tsv("/home/jsjukara/hla/allele_dosages/R6_HLA_DPB1_imputed_allele_dosages.tsv")
df_list[[5]] <- read_tsv("/home/jsjukara/hla/allele_dosages/R6_HLA_DQA1_imputed_allele_dosages.tsv")
df_list[[6]] <- read_tsv("/home/jsjukara/hla/allele_dosages/R6_HLA_DQB1_imputed_allele_dosages.tsv")
df_list[[7]] <- read_tsv("/home/jsjukara/hla/allele_dosages/R6_HLA_DRB1_imputed_allele_dosages.tsv")

#process names of the data frames and get ids
ids <- c(NULL)
for (i in 1:length(df_list)) {
    data <- df_list[[i]]
    names(data) <- paste("HLA_", names(data), sep="")
    names(data) <- gsub("\\*", "_", names(data))
    names(data) <- gsub("\\.", "_", names(data))
    names(data)[1] <- "FINNGENID"
    df_list[[i]] <- data
    ids <- unique(c(ids, data$FINNGENID))
}

#initialize df with ids to left_join the data into
hla_df <- data.frame(FINNGENID = ids)

for (i in 1:length(df_list)) {
    hla_df <- left_join(hla_df, df_list[[i]], by="FINNGENID")
}

#export data
saveRDS(hla_df, "/home/jsjukara/hla/hla_variant_df.rda")