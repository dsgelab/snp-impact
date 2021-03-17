setwd("/home/jsjukara/hla/files")
library(data.table)
library(tidyverse)

### read different HLA files to a list of data frames
df_list <- list(NULL)
df_list[[1]] <- read_tsv("R5_A_imputed.tsv", col_types="cccdd")
df_list[[2]] <- read_tsv("R5_B_imputed.tsv", col_types="cccdd")
df_list[[3]] <- read_tsv("R5_C_imputed.tsv", col_types="cccdd")
df_list[[4]] <- read_tsv("R5_DPB1_imputed.tsv", col_types="cccdd")
df_list[[5]] <- read_tsv("R5_DQA1_imputed.tsv", col_types="cccdd")
df_list[[6]] <- read_tsv("R5_DQB1_imputed.tsv", col_types="cccdd")
df_list[[7]] <- read_tsv("R5_DRB1_imputed.tsv", col_types="cccdd")

genes <- c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1")
for (i in 1:length(df_list)) {
    data <- df_list[[i]]
    data$allele1 <- paste("HLA-", genes[i], "*", data$allele1, sep="")
    data$allele2 <- paste("HLA-",genes[i], "*", data$allele2, sep="")
    data$gene <- genes[i]
    df_list[[i]] <- data
}

# join the list to a single data frame
data_all <- bind_rows(df_list)

names(data_all)[1] <- "FINNGENID"

### QC for R5 (might not be needed in R6)
# remove duplicated rows
dim(data_all)
data_all <- data_all[!duplicated(data_all),]
dim(data_all)

# remove rows with duplicated imputed alleles
data_all <- data_all[!duplicated(paste(data_all$FINNGENID, data_all$allele1, data_all$allele2)),]
dim(data_all)

## from duplicated rows for each FINNGENID and gene, only keep one of the rows
# shuffle rows to ensure that random imputed variant is kept from duplicated rows
data_all <- data_all[sample(1:nrow(data_all), replace=FALSE),]
# keep first row
data_all <- data_all %>%
    distinct(FINNGENID, gene, .keep_all=TRUE)
dim(temp)

### Transform the "long" format data into a "wide" format (each column is one HLA allele, each row one id)
# get unique allele names and ids
uniquenames <- unique(c(data_all$allele1, data_all$allele2))
ids <- unique(data_all$FINNGENID)

# make data frame where each row is one FINNGENID, each column is one variant
m1 <- as.data.frame(matrix(0,ncol=length(uniquenames), nrow=length(ids)))
colnames(m1) <- uniquenames
# assign FINNGENID as first column
m1$FINNGENID <- ids
m1 <- m1[,c(ncol(m1), 1:(ncol(m1)-1))]

# loop over rows in data_all, assign alleles to m1 in 0, 1, 2 coding
# takes little less than 2h in R5
for (i in 1:nrow(data_all)) {
    # get row of m1 for the current FINNGENID
    row <- which(m1$FINNGENID==data_all$FINNGENID[i])
    # assign both alleles to their places
    m1[row, data_all$allele1[i]] <- m1[row, data_all$allele1[i]] + 1
    m1[row, data_all$allele2[i]] <- m1[row, data_all$allele2[i]] + 1

}

# save
saveRDS(m1, "/home/jsjukara/hla/hla_variant_df.rda")