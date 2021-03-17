setwd("/home/jsjukara/ukbb/")
library(tidyverse)
library(data.table)

# read ukbb genotypes
geno <- fread("zcat skari_snps_source_ukbb.tsv.bgz")
geno <- transpose(geno)
geno <- as.data.frame(geno)

# wrangle ukbb genotypes
rsid_row <- geno[1,2:ncol(geno)]
geno <- geno[-1,]

locus <- geno[1,2:ncol(geno)]
locus <- gsub(":", "_", locus)
locus <- paste("chr", locus, sep="")
geno[1,] <- c("locus", locus)

alleles <- gsub('"', "", geno[2,2:ncol(geno)])
alleles <- gsub("\\[", "", alleles)
alleles <- gsub("\\]", "", alleles)
alleles <- gsub(",", "_", alleles)

# create data frame of loci and alleles from ukbb
ukbb_locus_allele <- data.frame(ukbb_locus=locus, alleles_ukbb=alleles, stringsAsFactors=FALSE)


setwd("/home/jsjukara/ukbb/")
library(tidyverse)
library(data.table)

# save ids and remove id column
geno_ids <- geno[3:nrow(geno),1]
geno <- geno[,-1]

# read variant mapping and create locus columns
variant_map <- fread("/home/jsjukara/finemap/finemap_variants_18102020.tsv")
variant_map$locus_hg19 <- paste("chr", gsub(":", "_", variant_map$locus_hg19), sep="")
variant_map$locus_hg38 <- paste("chr", gsub(":", "_", variant_map$locus_hg38), sep="")

# crate column that contains the finngen alleles
variant_map$alleles_finngen <- NA
for (i in 1:nrow(variant_map)) {
    vec <- unlist(strsplit(variant_map$variant_hg19[i], split="_"))
    allele <- paste(vec[3], vec[4], sep="_")
    variant_map$alleles_finngen[i] <- allele
}

# check number of multiallelic loci
table(variant_map$locus_hg19)[table(variant_map$locus_hg19)>1]

# join extracted ukbb alleles, with possible duplicated loci
variant_map$alleles_ukbb_extracted <- NA
for (i in 1:nrow(variant_map)) {
    # get current locus
    current_locus <- variant_map$locus_hg19[i]
    # get vector of loci and alleles extracted from ukbb
    extracted_loci <- locus[locus == current_locus]
    extracted_alleles <- alleles[locus == current_locus]
    # if no match, skip row
    if (length(extracted_loci) == 0) {
        next
    # if 1 match, assign
    } else if (length(extracted_loci) == 1) {
        variant_map$alleles_ukbb_extracted[i] <- extracted_alleles[1]
    # else assign first and duplicate rows and assign rest
    } else {
        variant_map$alleles_ukbb_extracted[i] <- extracted_alleles[1]
        counter = 2
        while (counter <= length(extracted_alleles)) {
            row <- variant_map[i,]
            row$alleles_ukbb_extracted <- extracted_alleles[counter]
            #row$duplicated <- 1
            variant_map <- rbind(variant_map, row)
            counter <- counter + 1
        }
        #print(variant_map[variant_map$locus_hg19 == current_locus,])
        
    }
}

# create reversed allele column
variant_map$alleles_ukbb_extracted_rev <- NA
for (i in 1:nrow(variant_map)) {
    if (!is.na(variant_map$alleles_ukbb_extracted[i])) {
        vec <- unlist(strsplit(variant_map$alleles_ukbb_extracted[i], split="_"))
        variant_map$alleles_ukbb_extracted_rev[i] <- paste(vec[2], vec[1], sep="_")
    }
}

# assign rows as "missing" if loci not extracted in ukbb, "nomatch" if alleles don't match
# "match" if reversed or unreversed allele matches
variant_map$alleles_match <- "missing"
variant_map$ukbb_reversed <- 0
for (i in 1:nrow(variant_map)) {
    if (!is.na(variant_map$alleles_ukbb_extracted[i])) {
            if (variant_map$alleles_finngen[i] != variant_map$alleles_ukbb_extracted[i] &
                variant_map$alleles_finngen[i] != variant_map$alleles_ukbb_extracted_rev[i]) {
                variant_map$alleles_match[i] <- "nomatch"
        }  else {
                variant_map$alleles_match[i] <- "match"
                if (variant_map$alleles_finngen[i] == variant_map$alleles_ukbb_extracted_rev[i]) {
                    variant_map$ukbb_reversed[i] <- 1
                }
            }
    }

}

# check reversed alleles
variant_map %>% filter(ukbb_reversed==1)


variant_map <- variant_map[!duplicated(paste(variant_map$variant_hg19,
                             variant_map$variant_hg38,
                            variant_map$alleles_finngen,
                            variant_map$alleles_ukbb_extracted)),]

# create variant_map that contains only matching columns for ukbb
variant_map_match <- variant_map %>%
    filter(alleles_match=="match")


### loop over geno columns, set allele ordering as same as in variant_map, store indices of geno columns to keep
### and store indices where the alleles are reversed
# initialize vector to contain geno column indices that will be kept 
col_indices_keep <- c(NULL)
# initialize vector to contain geno column indices that will have alleles reversed
col_indices_rev <- c(NULL)
colnames_ukbb <- c(NULL)
# loop over variants
for (j in 1:ncol(geno)) {
    # get current locus
    current_locus <- geno[1,j]
    # get current allele in "C_T" format from extracted UKBB data
    current_extracted_allele<-paste(unlist(strsplit(gsub('\\[|\\]|\"', "", geno[2,j]), ",")), collapse="_")
    
    # get rows of variant_map_match that match current loci
    df_loci_match <- variant_map_match %>% filter(locus_hg19 == current_locus)
    
    if (nrow(df_loci_match)>0) {
        # loop over the rows and assign matching alleles in same order/format as the variant_map one
        for (i in 1:nrow(df_loci_match)) {
            # define reversed version of FinnGen allele for checking match
            alleles_rev <- unlist(strsplit(df_loci_match$alleles_finngen[i], split="_"))
            alleles_rev <- paste(alleles_rev[2], alleles_rev[1], sep="_")
            # if allele matches with unreversed allele
            if (df_loci_match$alleles_finngen[i] == current_extracted_allele) {
                # set allele to unreversed allele
                geno[2,j] <- current_extracted_allele
                # save column indices in genotype data from UKBB to keep (those with matching locus + allele)
                col_indices_keep <- c(col_indices_keep, j)
            # if allele matches to reversed allele
            } else if (alleles_rev == current_extracted_allele) {
                # set allele to unreversed allele
                geno[2,j] <- df_loci_match$alleles_finngen[i]
                # save reversed column indices
                col_indices_rev <- c(col_indices_rev, j)
                # save column indices in genotype data from UKBB to keep (those with matching locus + allele)
                col_indices_keep <- c(col_indices_keep, j)
            }
        }
    }
    
}

# remove duplicated column indices
col_indices_keep <- unique(col_indices_keep)

# set column names
colnames(geno) <- paste(geno[1,], geno[2,], sep="_")
geno <- geno[-c(1,2),]
# reverse coding of reversed alleles
for (j in col_indices_rev) {
    geno[,j] <- as.integer(geno[,j])
    geno[,j] <- (geno[,j]-2)*-1
}

# leave out unmatched columns
geno_ukbb <- geno[,col_indices_keep]

# add id column
geno_ukbb <- add_column(geno_ukbb, eid=geno_ids, .before=1)


# rename UKBB variant ids from hg19 to hg38
for (j in 2:ncol(geno_ukbb)) {
    colname <- colnames(geno_ukbb)[j]
    fg_colname <- variant_map_match$variant_hg38[variant_map_match$variant_hg19==colname]
    names(geno_ukbb)[j] <- fg_colname
}

# convert ukbb genotype columns as integers from character
for (i in 2:ncol(geno_ukbb)) {
    geno_ukbb[,i] <- as.integer(geno_ukbb[,i])
}

### FINNGEN GENOTYPES
geno <- fread("/home/jsjukara/finngen_data/GT_21102020.tsv", data.table=FALSE)
dim(geno)
coln <- as.character(geno[,3])
loci <- as.character(geno[,1])
alleles <- as.character(geno[,2])
geno <- geno[,-c(1:3)]
ids <- colnames(geno)
geno <- t(geno)

geno <- as.data.frame(geno)
colnames(geno) <- coln
geno$FINNGENID <- ids
row.names(geno) <- NULL

head(geno)

#names(geno) <- c(coln, "FINNGENID")

#rsids1 <- coln

# go from 2,3,4 coding to 0 1 2 coding(..?)
#for (i in 1:(ncol(geno)-1)) {
#    geno[,i] <- as.integer(geno[,i])-2
#}

# convert to integers
#for (i in 1:(ncol(geno)-1)) {
#    geno[,i] <- as.integer(geno[,i])
#}

geno_fg <- geno
rm(geno)

# make sure all alleles are the right way around in geno_fg_fg
loci <- gsub(":", "_", loci)

fg_loci_alleles <- data.frame(locus_hg38 = loci, variant_finngen = names(geno_fg)[1:(ncol(geno_fg)-1)], stringsAsFactors=FALSE)

fg_loci_alleles <- left_join(fg_loci_alleles, variant_map[,c("locus_hg38", "variant_hg38")], by="locus_hg38")
head(fg_loci_alleles)
sum(is.na(fg_loci_alleles$variant_hg38))
#fg_loci_alleles[fg_loci_alleles$variant_hg38 != fg_loci_alleles$variant_finngen,]

fg_loci_alleles$alleles_extracted <- NA
fg_loci_alleles$alleles_extracted_rev <- NA
fg_loci_alleles$alleles_map <- NA
fg_loci_alleles$match <- 0
fg_loci_alleles$rev_match <- 0
for (i in 1:nrow(fg_loci_alleles)) {
    vec <- unlist(strsplit(fg_loci_alleles$variant_finngen[i], split="_"))
    alleles <- paste(vec[3], vec[4], sep="_")
    fg_loci_alleles$alleles_extracted[i] <- alleles
    
    alleles <- paste(vec[4], vec[3], sep="_")
    fg_loci_alleles$alleles_extracted_rev[i] <- alleles
    
    vec <- unlist(strsplit(fg_loci_alleles$variant_hg38[i], split="_"))
    alleles <- paste(vec[3], vec[4], sep="_")
    fg_loci_alleles$alleles_map[i] <- alleles
    
    if(fg_loci_alleles$alleles_extracted[i] == alleles) {
        fg_loci_alleles$match[i] <- 1
    } else if (fg_loci_alleles$alleles_extracted_rev[i] == alleles) {
        fg_loci_alleles$rev_match[i] <- 1
    }
}

rev_variants <- fg_loci_alleles$variant_finngen[fg_loci_alleles$rev_match==1]
rev_indices <- which(names(geno_fg) %in% rev_variants)
#head(geno_fg[,rev_indices])

# reverse coding of reversed alleles
for (j in 1:ncol(geno_fg)) {
    if (names(geno_fg)[j] %in% rev_variants) {
        vec <- unlist(strsplit(names(geno_fg)[j], split="_"))
        locus <- paste(vec[1], vec[2], sep="_")
        alleles <- paste(vec[4], vec[3], sep="_")
        names(geno_fg)[j] <- paste(locus, alleles, sep="_")
        geno_fg[,j] <- (geno_fg[,j]-2)*-1
    }
}

fg_keep_variants <- fg_loci_alleles$variant_hg38[fg_loci_alleles$match_any==1]
length(fg_keep_variants)
sum(fg_keep_variants %in% names(geno_fg))

variants <- intersect(names(geno_fg), names(geno_ukbb))

# function to generate data frame of genotype frequencies
get_genotype_freqs <- function(data, variants) {
    genotype_freqs <- data.frame(rsid = variants, p0=NA, p1=NA, p2=NA, stringsAsFactors=FALSE)
    for (i in 1:length(variants)) {
        freqs <- as.vector(table(data[,variants[i]])/nrow(data))
        if (length(freqs)==3) {
            genotype_freqs[i,2:4] <- freqs

        } else if (length(freqs)==2) {
            genotype_freqs[i,2:4] <- c(freqs,0)
        } else {
            print("error")
        }    
    }
    
    genotype_freqs$maf <- (1*genotype_freqs$p1+2*genotype_freqs$p2)/2
    genotype_freqs$mafover1 <- ifelse(genotype_freqs$maf>0.01, 1, 0)
    return(genotype_freqs)
}

genotype_freqs_ukbb <- get_genotype_freqs(geno_ukbb, variants=variants)
genotype_freqs_fg <- get_genotype_freqs(geno_fg, variants=variants)

# make FINNGENID as first column
n <- ncol(geno_fg)
geno_fg <- geno_fg[,c(n,2:(n-1))]

# temporary ad hoc fix to one discrepant allele
#geno_ukbb[,"chr1_145804971_G_C"] <- -1*(geno_ukbb[,"chr1_145804971_G_C"])+2

# plot maf to ensure coding is same way
plot(genotype_freqs_ukbb$maf, genotype_freqs_fg$maf)

# plot maf to ensure coding is same way
plot(genotype_freqs_ukbb$maf, genotype_freqs_fg$maf)

saveRDS(geno_fg, "/home/jsjukara/finngen_data/preprocessed_genotypes_fg.rds", compress=FALSE)
saveRDS(geno_ukbb, "/home/jsjukara/ukbb/preprocessed_genotypes_ukbb.rds", compress=FALSE)

setwd("/home/jsjukara/ukbb/")
library(tidyverse)
library(data.table)

# read ukbb genotypes
geno <- fread("zcat skari_snps_source_ukbb.tsv.bgz")
geno <- transpose(geno)
geno <- as.data.frame(geno)

# wrangle ukbb genotypes
rsid_row <- geno[1,2:ncol(geno)]
geno <- geno[-1,]

locus <- geno[1,2:ncol(geno)]
locus <- gsub(":", "_", locus)
locus <- paste("chr", locus, sep="")
geno[1,] <- c("locus", locus)

alleles <- gsub('"', "", geno[2,2:ncol(geno)])
alleles <- gsub("\\[", "", alleles)
alleles <- gsub("\\]", "", alleles)
alleles <- gsub(",", "_", alleles)

# create data frame of loci and alleles from ukbb
ukbb_locus_allele <- data.frame(ukbb_locus=locus, alleles_ukbb=alleles, stringsAsFactors=FALSE)

# save ids and remove id column
geno_ids <- geno[3:nrow(geno),1]
geno <- geno[,-1]

# read variant mapping and create locus columns
variant_map <- fread("/home/jsjukara/finemap/finemap_variants_18102020.tsv")
variant_map$locus_hg19 <- paste("chr", gsub(":", "_", variant_map$locus_hg19), sep="")
variant_map$locus_hg38 <- paste("chr", gsub(":", "_", variant_map$locus_hg38), sep="")

# crate column that contains the finngen alleles
variant_map$alleles_finngen <- NA
for (i in 1:nrow(variant_map)) {
    vec <- unlist(strsplit(variant_map$variant_hg19[i], split="_"))
    allele <- paste(vec[3], vec[4], sep="_")
    variant_map$alleles_finngen[i] <- allele
}

# check number of multiallelic loci
table(variant_map$locus_hg19)[table(variant_map$locus_hg19)>1]

# join extracted ukbb alleles, with possible duplicated loci
variant_map$alleles_ukbb_extracted <- NA
for (i in 1:nrow(variant_map)) {
    # get current locus
    current_locus <- variant_map$locus_hg19[i]
    # get vector of loci and alleles extracted from ukbb
    extracted_loci <- locus[locus == current_locus]
    extracted_alleles <- alleles[locus == current_locus]
    # if no match, skip row
    if (length(extracted_loci) == 0) {
        next
    # if 1 match, assign
    } else if (length(extracted_loci) == 1) {
        variant_map$alleles_ukbb_extracted[i] <- extracted_alleles[1]
    # else assign first and duplicate rows and assign rest
    } else {
        variant_map$alleles_ukbb_extracted[i] <- extracted_alleles[1]
        counter = 2
        while (counter <= length(extracted_alleles)) {
            row <- variant_map[i,]
            row$alleles_ukbb_extracted <- extracted_alleles[counter]
            #row$duplicated <- 1
            variant_map <- rbind(variant_map, row)
            counter <- counter + 1
        }
        #print(variant_map[variant_map$locus_hg19 == current_locus,])
        
    }
}

# create reversed allele column
variant_map$alleles_ukbb_extracted_rev <- NA
for (i in 1:nrow(variant_map)) {
    if (!is.na(variant_map$alleles_ukbb_extracted[i])) {
        vec <- unlist(strsplit(variant_map$alleles_ukbb_extracted[i], split="_"))
        variant_map$alleles_ukbb_extracted_rev[i] <- paste(vec[2], vec[1], sep="_")
    }
}

# assign rows as "missing" if loci not extracted in ukbb, "nomatch" if alleles don't match
# "match" if reversed or unreversed allele matches
variant_map$alleles_match <- "missing"
variant_map$ukbb_reversed <- 0
for (i in 1:nrow(variant_map)) {
    if (!is.na(variant_map$alleles_ukbb_extracted[i])) {
            if (variant_map$alleles_finngen[i] != variant_map$alleles_ukbb_extracted[i] &
                variant_map$alleles_finngen[i] != variant_map$alleles_ukbb_extracted_rev[i]) {
                variant_map$alleles_match[i] <- "nomatch"
        }  else {
                variant_map$alleles_match[i] <- "match"
                if (variant_map$alleles_finngen[i] == variant_map$alleles_ukbb_extracted_rev[i]) {
                    variant_map$ukbb_reversed[i] <- 1
                }
            }
    }

}

# check reversed alleles
#variant_map %>% filter(ukbb_reversed==1)

# remove possible duplicated rows
variant_map <- variant_map[!duplicated(paste(variant_map$variant_hg19,
                             variant_map$variant_hg38,
                            variant_map$alleles_finngen,
                            variant_map$alleles_ukbb_extracted)),]

# create variant_map that contains only matching columns for ukbb
variant_map_match <- variant_map %>%
    filter(alleles_match=="match")


### loop over geno columns, set allele ordering as same as in variant_map, store indices of geno columns to keep
### and store indices where the alleles are reversed
# initialize vector to contain geno column indices that will be kept 
col_indices_keep <- c(NULL)
# initialize vector to contain geno column indices that will have alleles reversed
col_indices_rev <- c(NULL)
colnames_ukbb <- c(NULL)
# loop over variants
for (j in 1:ncol(geno)) {
    # get current locus
    current_locus <- geno[1,j]
    # get current allele in "C_T" format from extracted UKBB data
    current_extracted_allele<-paste(unlist(strsplit(gsub('\\[|\\]|\"', "", geno[2,j]), ",")), collapse="_")
    
    # get rows of variant_map_match that match current loci
    df_loci_match <- variant_map_match %>% filter(locus_hg19 == current_locus)
    
    if (nrow(df_loci_match)>0) {
        # loop over the rows and assign matching alleles in same order/format as the variant_map one
        for (i in 1:nrow(df_loci_match)) {
            # define reversed version of FinnGen allele for checking match
            alleles_rev <- unlist(strsplit(df_loci_match$alleles_finngen[i], split="_"))
            alleles_rev <- paste(alleles_rev[2], alleles_rev[1], sep="_")
            # if allele matches with unreversed allele
            if (df_loci_match$alleles_finngen[i] == current_extracted_allele) {
                # set allele to unreversed allele
                geno[2,j] <- current_extracted_allele
                # save column indices in genotype data from UKBB to keep (those with matching locus + allele)
                col_indices_keep <- c(col_indices_keep, j)
            # if allele matches to reversed allele
            } else if (alleles_rev == current_extracted_allele) {
                # set allele to unreversed allele
                geno[2,j] <- df_loci_match$alleles_finngen[i]
                # save reversed column indices
                col_indices_rev <- c(col_indices_rev, j)
                # save column indices in genotype data from UKBB to keep (those with matching locus + allele)
                col_indices_keep <- c(col_indices_keep, j)
            }
        }
    }
    
}

# remove duplicated column indices
col_indices_keep <- unique(col_indices_keep)

# set column names
colnames(geno) <- paste(geno[1,], geno[2,], sep="_")
geno <- geno[-c(1,2),]
# reverse coding of reversed alleles
for (j in col_indices_rev) {
    geno[,j] <- as.integer(geno[,j])
    geno[,j] <- (geno[,j]-2)*-1
}

# leave out unmatched columns
geno_ukbb <- geno[,col_indices_keep]

# add id column
geno_ukbb <- add_column(geno_ukbb, eid=geno_ids, .before=1)


# rename UKBB variant ids from hg19 to hg38
for (j in 2:ncol(geno_ukbb)) {
    colname <- colnames(geno_ukbb)[j]
    fg_colname <- variant_map_match$variant_hg38[variant_map_match$variant_hg19==colname]
    names(geno_ukbb)[j] <- fg_colname
}

# convert ukbb genotype columns as integers from character
for (i in 2:ncol(geno_ukbb)) {
    geno_ukbb[,i] <- as.integer(geno_ukbb[,i])
}

### FINNGEN GENOTYPES
geno <- fread("/home/jsjukara/finngen_data/GT_21102020.tsv", data.table=FALSE)
dim(geno)
coln <- as.character(geno[,3])
loci <- as.character(geno[,1])
alleles <- as.character(geno[,2])
geno <- geno[,-c(1:3)]
head(geno)
ids <- colnames(geno)
geno <- t(geno)

geno <- as.data.frame(geno)
colnames(geno) <- coln
geno$FINNGENID <- ids
row.names(geno) <- NULL


#names(geno) <- c(coln, "FINNGENID")

#rsids1 <- coln

# go from 2,3,4 coding to 0 1 2 coding(..?)
#for (i in 1:(ncol(geno)-1)) {
#    geno[,i] <- as.integer(geno[,i])-2
#}

# convert to integers
#for (i in 1:(ncol(geno)-1)) {
#    geno[,i] <- as.integer(geno[,i])
#}

# make sure all alleles are the right way around in geno_fg_fg

geno_fg <- geno
rm(geno)

loci <- gsub(":", "_", loci)

fg_loci_alleles <- data.frame(locus_hg38 = loci, variant_finngen = names(geno_fg)[1:(ncol(geno_fg)-1)], stringsAsFactors=FALSE)

fg_loci_alleles <- left_join(fg_loci_alleles, variant_map[,c("locus_hg38", "variant_hg38")], by="locus_hg38")
head(fg_loci_alleles)
sum(is.na(fg_loci_alleles$variant_hg38))
#fg_loci_alleles[fg_loci_alleles$variant_hg38 != fg_loci_alleles$variant_finngen,]

rev_variants <- fg_loci_alleles$variant_finngen[fg_loci_alleles$rev_match==1]

fg_loci_alleles$alleles_extracted <- NA
fg_loci_alleles$alleles_extracted_rev <- NA
fg_loci_alleles$alleles_map <- NA
fg_loci_alleles$match <- 0
fg_loci_alleles$rev_match <- 0
for (i in 1:nrow(fg_loci_alleles)) {
    vec <- unlist(strsplit(fg_loci_alleles$variant_finngen[i], split="_"))
    alleles <- paste(vec[3], vec[4], sep="_")
    fg_loci_alleles$alleles_extracted[i] <- alleles
    
    alleles <- paste(vec[4], vec[3], sep="_")
    fg_loci_alleles$alleles_extracted_rev[i] <- alleles
    
    vec <- unlist(strsplit(fg_loci_alleles$variant_hg38[i], split="_"))
    alleles <- paste(vec[3], vec[4], sep="_")
    fg_loci_alleles$alleles_map[i] <- alleles
    
    if(fg_loci_alleles$alleles_extracted[i] == alleles) {
        fg_loci_alleles$match[i] <- 1
    } else if (fg_loci_alleles$alleles_extracted_rev[i] == alleles) {
        fg_loci_alleles$rev_match[i] <- 1
    }
}

rev_indices <- which(names(geno_fg) %in% rev_variants)
#head(geno_fg[,rev_indices])

# reverse coding of reversed alleles
for (j in 1:ncol(geno_fg)) {
    if (names(geno_fg)[j] %in% rev_variants) {
        vec <- unlist(strsplit(names(geno_fg)[j], split="_"))
        locus <- paste(vec[1], vec[2], sep="_")
        alleles <- paste(vec[4], vec[3], sep="_")
        names(geno_fg)[j] <- paste(locus, alleles, sep="_")
        geno_fg[,j] <- (geno_fg[,j]-2)*-1
    }
}

fg_keep_variants <- fg_loci_alleles$variant_hg38[fg_loci_alleles$match_any==1]
length(fg_keep_variants)
sum(fg_keep_variants %in% names(geno_fg))

variants <- intersect(names(geno_fg), names(geno_ukbb))

# function to generate data frame of genotype frequencies
get_genotype_freqs <- function(data, variants) {
    genotype_freqs <- data.frame(rsid = variants, p0=NA, p1=NA, p2=NA, stringsAsFactors=FALSE)
    for (i in 1:length(variants)) {
        freqs <- as.vector(table(data[,variants[i]])/nrow(data))
        if (length(freqs)==3) {
            genotype_freqs[i,2:4] <- freqs

        } else if (length(freqs)==2) {
            genotype_freqs[i,2:4] <- c(freqs,0)
        } else {
            print("error")
        }    
    }
    
    genotype_freqs$maf <- (1*genotype_freqs$p1+2*genotype_freqs$p2)/2
    genotype_freqs$mafover1 <- ifelse(genotype_freqs$maf>0.01, 1, 0)
    return(genotype_freqs)
}

genotype_freqs_ukbb <- get_genotype_freqs(geno_ukbb, variants=variants)
genotype_freqs_fg <- get_genotype_freqs(geno_fg, variants=variants)

# make FINNGENID as first column
n <- ncol(geno_fg)
geno_fg <- geno_fg[,c(n,2:(n-1))]


# plot maf to ensure coding is same way
plot(genotype_freqs_ukbb$maf, genotype_freqs_fg$maf)

saveRDS(geno_fg, "/home/jsjukara/finngen_data/preprocessed_genotypes_fg.rds", compress=FALSE)
saveRDS(geno_ukbb, "/home/jsjukara/ukbb/preprocessed_genotypes_ukbb.rds", compress=FALSE)

# flip variant coding when MAF over 0.5
#genotype_freqs$flipped <- 0
#for (j in 1:nrow(genotype_freqs)) {
#    if(genotype_freqs$maf[j]>0.5) {
#        p0 <- genotype_freqs$p0[j]
#        genotype_freqs$p0[j] <- genotype_freqs$p2[j]
#        genotype_freqs$p2[j] <- p0
#        genotype_freqs$maf[j] <- 1-genotype_freqs$maf[j]
#        genotype_freqs$flipped[j] <- 1
#    }
#}

# flip coding in fg_df
#for (i in 1:nrow(genotype_freqs)) {
#    if(genotype_freqs$flipped[i] == 1) {
#        rs <- as.character(genotype_freqs$rsid[i])
#        fg_df[,rs] <- (fg_df[,rs] - 2)*-1
#    }
#}


## Form FinnGen endpoints in UKBB

library(readxl)
library(tidyverse)
library(data.table)
library(lubridate)
library(tictoc)
setwd("/home/jsjukara/ukbb/")

remove_special <- function(x) {
    string_vec <- x
    # substitute underscores into spaces
    string_vec <- gsub(" ", "_", string_vec)
    # remove apostrophes
    string_vec <- gsub("\\'", "", string_vec)
    # remove commas
    string_vec <- gsub(",", "", string_vec)
    # remove dashes
    string_vec <- gsub("\\-", "_", string_vec)
    # remove \
    string_vec <- gsub("/", "_", string_vec)
    return(string_vec)
}

# define the span of follow-up years
fu_startyear <- 1998
fu_endyear <- 2020.33

# get ids of european ancestry individuals
pops <- fread("ukb_diverse_pops_pruned.tsv", data.table=FALSE)
ids_ukbb <- pops$s[pops$pop=="EUR"]
length(ids_ukbb)
nrow(pops)


# read FinnGen ontology excel
ontology <- suppressWarnings(read_excel("FINNGEN_ENDPOINTS_DF5_V2_2020-02-11_public.xlsx"))

# read the table of endpoints that contain the sets of FinnGen endpoints that define the new endpoint
# e.g. "C3_COLON|C3_RECTUM"
map <- fread("/home/jsjukara/Auxiliary data/gbd_to_finngen_ontology_manual_9_2020.txt")
map <- map[map$fg_endpoints!="",]
map <- map[map$include==1,]
names(map) <- tolower(names(map))
map$cause <- remove_special(map$cause)

# initialize column that will contain the regex for the custom endpoint ICD10 codes
map$icd10_regex <- NA

# define recursive function to get regex of ICD-10 codes for a set of endpoints
#e.g. endpoints = "C3_COLON|C3_RECTUM" will output regular expression for all diagnoses
# under the endpoints
get_icd10_regex <- function(endpoints, ontology=ontology) {
  endpoints <- unlist(strsplit(endpoints, split="\\|"))
  # if current endpoint not a composite endpoint
  if (length(endpoints) == 1) {
    # get corresponding from from finngen ontology
    row <- ontology[ontology$NAME == endpoints[1],]
    # if endpoint does not have children, get the regular expression
    if (is.na(row$INCLUDE[1])) {
      # get vector of ICD10 regex's from diagnoses, causes of death, and cancer topography
      # for cancer HD_ICD_10 might be empty, whereas the codes live in COD or TOPO
      icd10_regex <- c(row$HD_ICD_10[1], row$COD_ICD_10[1], row$CANC_TOPO[1])
      icd10_regex <- unique(icd10_regex[!is.na(icd10_regex)])
      # if empty, return "empty"
      if (length(icd10_regex) == 0) {
          return("empty")
      }
      # paste them with | together
      icd10_regex <- paste(icd10_regex, collapse="|")
      return(icd10_regex)
    } else { # else try to get diagnoses as above, and recurse into the composite endpoints
      icd10_regex <- c(row$HD_ICD_10[1], row$COD_ICD_10[1], row$CANC_TOPO[1])
      icd10_regex <- unique(icd10_regex[!is.na(icd10_regex)])
      # get composite endpoints
      composite_endpoints <- row$INCLUDE[1]
      # if empty, just recurse to composite endpoints
      if (length(icd10_regex) == 0) {
          return(get_icd10_regex(endpoints = composite_endpoints, ontology=ontology))
      }
      # else get recurse to composite endpoints and paste them with current regex
      icd10_regex_composite <- get_icd10_regex(endpoints = composite_endpoints, ontology=ontology)
      icd10_regex <- c(icd10_regex, icd10_regex_composite)
      icd10_regex <- paste(icd10_regex, collapse="|")
      return(icd10_regex)
    }
  } else { # else if current endpoint is a composite endpoint
    # initialize helper index for while loop
    i=1
    # loop along endpoints until first non empty regex is found
    while (is.na(get_icd10_regex(endpoints = endpoints[i], ontology=ontology))) {
      i <- i + 1
      if (i > length(endpoints)) {
        break
      }
    }
    # get get the first endpoint's regex by recursion
    icd10_regex_new <- get_icd10_regex(endpoints = endpoints[i], ontology=ontology)
    icd10_regex <- icd10_regex_new
    # loop over the rest of the endpoints, get regex by recursion
    for (endpoint in endpoints[-1]) {
      icd10_regex_new <- get_icd10_regex(endpoints = endpoint, ontology=ontology)
      #icd10_regex_new <- ontology$row$HD_ICD_10[ontology$NAME == endpoint]
      if (!is.na(icd10_regex_new)) {
        icd10_regex <- paste(icd10_regex, icd10_regex_new, sep="|")
      }
    }
    return(icd10_regex)
  }
}


for (i in 1:nrow(map)) {
  map$icd10_regex[i] <- get_icd10_regex(map$fg_endpoints[i], ontology=ontology)
}

# manually set AMD endpoint
map$icd10_regex[map$fg_endpoints=="H7_AMD"] <- "H353"

#remove $!$ that comes from empty endpoint definitions
map$icd10_regex <- gsub("\\|\\$!\\$", "", map$icd10_regex)
#remove "empty" from the regex, that result as byproduct of translation loop
map$icd10_regex <- gsub("empty", "", map$icd10_regex)
#replace || with |
map$icd10_regex <- gsub("\\|\\|", "|", map$icd10_regex)
#remove | at end
map$icd10_regex <- gsub("\\|$", "", map$icd10_regex)
#remove | at beginning
map$icd10_regex <- gsub("^\\|", "", map$icd10_regex)

#export mapping
write_tsv(map, "/home/jsjukara/Auxiliary data/gbd_to_finngen_ontology_with_ukbb_9_2020.tsv")

# causes of death
death <- fread("ukb31063.death.20200706.txt", data.table=FALSE)
death_cause <- fread("ukb31063.death_cause.20200706.txt", data.table=FALSE)

death$year_of_death <- decimal_date(as.Date(death$date_of_death, format="%d/%m/%Y"))
death$DEATH <- 1

death_cause <- left_join(death_cause, death[,c("eid", "year_of_death")])

death_cause <- death_cause[,-4]

## demographics
demographics <- fread("export_for_map_analysis.tsv", data.table=FALSE)
demographics <- demographics[,c("f.eid", "f.31.0.0", "f.53.0.0", "f.21022.0.0")]
names(demographics) <- c("eid", "FEMALE", "date_assessment", "age_at_assessment")

demographics$FEMALE <- ifelse(demographics$FEMALE == 1, 0, 1)
demographics$date_assessment <- as.Date(demographics$date_assessment, format="%Y-%m-%d")
demographics$year_assessment <- decimal_date(demographics$date_assessment)
demographics$YEAR_OF_BIRTH <- demographics$year_assessment - demographics$age_at_assessment
demographics$FU_START_AGE <- fu_startyear - demographics$YEAR_OF_BIRTH


# define end of follow-up as year of death or end of follow-up, whichever comes first
#join year of death and death indicator
demographics <- left_join(demographics, death[c("eid","year_of_death", "DEATH")], by="eid")
# if eid died, set as age of death
demographics$DEATH_AGE <- ifelse(is.na(demographics$year_of_death), NA, demographics$year_of_death - demographics$YEAR_OF_BIRTH)
# if did not die, set as age at end of follow-up
demographics$AGE_AT_DEATH_OR_NOW <- ifelse(is.na(demographics$DEATH_AGE), fu_endyear - demographics$YEAR_OF_BIRTH, demographics$DEATH_AGE)
#make values past 2020 be 2020
demographics$AGE_AT_DEATH_OR_NOW <- ifelse(demographics$AGE_AT_DEATH_OR_NOW>fu_endyear, fu_endyear, demographics$AGE_AT_DEATH_OR_NOW)

demographics$FU_END_AGE <- demographics$AGE_AT_DEATH_OR_NOW

# set NA DEATH as 0
demographics$DEATH <- ifelse(!is.na(demographics$DEATH), demographics$DEATH, 0) 

demographics$eid <- as.character(demographics$eid)

# remove duplicated rows
demographics <- demographics[!duplicated(demographics),]

# only keep needed columns
selected_colnames <- c("eid", "YEAR_OF_BIRTH", "FU_START_AGE",
                       "FU_END_AGE", "AGE_AT_DEATH_OR_NOW", "DEATH",
                       "DEATH_AGE", "FEMALE", paste("PC", 1:20, sep=""))
demographics <- demographics[,names(demographics) %in% selected_colnames]

# read QC file for PCs and join them
qc <- fread("ukb31063.sample_qc.tsv", data.table=FALSE)
names(qc)[1] <- "eid"
qc$eid <- as.character(qc$eid)
head(qc)
demographics <- left_join(demographics, qc[,c("eid", paste("PC", 1:20, sep=""))], by="eid")
head(demographics)

## read hospital data
hesin <- fread("ukb31063.hesin.20200810.txt", data.table=FALSE)
hesin_diag <- fread("ukb31063.hesin_diag.20200810.txt", data.table=FALSE)

# make date epistart as posixct, and create year column designating year of epistart
hesin$epistart <- as.Date(hesin$epistart, format="%d/%m/%Y")
hesin$year <- decimal_date(hesin$epistart)

# assume missing epidur are 0
hesin$epidur <- as.integer(ifelse(is.na(hesin$epidur), 0, hesin$epidur))

# if year is missing, try to get it from admidate - epidur (used for dsource == 10)
hesin$year <- ifelse(!is.na(hesin$year), hesin$year, decimal_date(as.Date(hesin$admidate, format="%d/%m/%Y") - hesin$epidur))

# remove episodes without year
hesin <- hesin %>% filter(!is.na(year))
# filter by year
hesin <- hesin %>% filter(year >= fu_startyear & year < fu_endyear)

# only keep columns to be used
hesin <- hesin[,c("eid", "ins_index", "year")]
hesin_diag <- hesin_diag[,c("eid", "ins_index", "arr_index", "diag_icd10")]

# get unique diagnoses and see which diagnoses the regular expression for endpoints matches to
#diagnoses <- unique(hesin_diag$diag_icd10)
#map$matches_ukbb <- NA
#for (i in 1:nrow(map)) {
#    map$matches_ukbb[i] <- paste(diagnoses[grepl(map$icd10_regex[i], diagnoses)], collapse="|")
#}

hesin_diag <- left_join(hesin_diag, hesin[,c("eid", "ins_index", "year")], by=c("eid", "ins_index"))

hesin_diag <- hesin_diag %>% filter(!is.na(year))

# join causes of death to hesin_diag to get causes of death as diagnoses
names(death_cause)[4:5] <- c("diag_icd10", "year")
hesin_diag <- rbind(hesin_diag, death_cause)



## Convert UKBB diagnoses into FinnGen ontology by matching them to the regular expressions
# takes around 5 hours
# preallocate vectors, assuming there wont be more than 0.6 times number of rows in hesin_diag as causes
# preallocating is necessary as appending results to a vector or to a data frame will be horribly slow
# as the size of the objects grows
matches <- rep(NA, floor(nrow(hesin_diag)*0.6))
ids <- matches
years <- matches
counter <- 1
tic()
for (i in 1:nrow(hesin_diag)) {
    # loop over map, for each icd10 match append cause, id and year
    # allows for an icd10 code to be matched with multiple causes
    for (j in 1:nrow(map)) {
        if (grepl(map$icd10_regex[j], hesin_diag$diag_icd10[i])) {
            matches[counter] <- map$cause[j]
            ids[counter] <- hesin_diag$eid[i]
            years[counter] <- hesin_diag$year[i]
            counter <- counter + 1
        }

    }
    
}
# formulate long_ukbb from the vectors
long_ukbb <- data.frame(cause=matches,
                        eid=as.character(ids),
                        year=years,
                        stringsAsFactors=FALSE)
toc()

long_ukbb <- long_ukbb[!is.na(long_ukbb$cause),]

long_ukbb <- long_ukbb[,c(2,3,1)]

long_ukbb <- long_ukbb %>%
    arrange(eid, year, cause)

# formulate first occurrence data frame in long format
long_ukbb_first <- long_ukbb %>%
    mutate(eid_cause = paste(eid, cause)) %>%
    distinct(eid_cause, .keep_all = TRUE) %>%
    select(eid, year, cause)

long_ukbb_first <- left_join(long_ukbb_first, demographics[,c("eid", "YEAR_OF_BIRTH")], by="eid")

long_ukbb_first$eid <- as.character(long_ukbb_first$eid)

# remove possible duplicated rows
long_ukbb_first <- long_ukbb_first[!duplicated(long_ukbb_first),]

# define wide format data for causes and cause_AGE's
# start with matrix that is converted to data frame with appropriate column names
n_causes <- length(unique(long_ukbb_first$cause))
n_individuals <- length(unique(ids_ukbb))
wide_ukbb_first <- matrix(NA, ncol=n_causes*2 + 1, nrow=n_individuals)
wide_ukbb_first <- as.data.frame(wide_ukbb_first)
colnames(wide_ukbb_first) <- c("eid", unique(long_ukbb_first$cause), paste(unique(long_ukbb_first$cause), "_AGE", sep=""))
wide_ukbb_first$eid <- as.character(unique(ids_ukbb))

ids <- as.character(wide_ukbb_first$eid)
# formulate wide format cause file, with binary cause as 0/1
# and cause_AGE as the year of the endpoint first happening


tic()
for (j in 2:((ncol(wide_ukbb_first)+1)/2)) {
    #get current cause
    current_cause <- names(wide_ukbb_first)[j]
    events <- long_ukbb_first %>%
        filter(cause == current_cause)
    joindf <- data.frame(eid = ids, stringsAsFactors=FALSE)
    joindf <- left_join(joindf, events, by="eid")
    joindf$age_at_event <- joindf$year - joindf$YEAR_OF_BIRTH
    wide_ukbb_first[,current_cause] <- ifelse(!is.na(joindf$cause), 1, 0)
    wide_ukbb_first[,paste(current_cause, "_AGE", sep="")] <- joindf$age_at_event
}

toc()

# join cause data to demographics
ukbb_df <- left_join(demographics, wide_ukbb_first, by="eid")

# for each cause_AGE column, for NA values (cause did not happen) substitute NA as the age at end of follow-up
for (j in 1:ncol(ukbb_df)) {
    colname <- names(ukbb_df)[j]
    if (grepl("_AGE", colname)) {
        ukbb_df[,colname] <- ifelse(!is.na(ukbb_df[,colname]),
                                    ukbb_df[,colname],
                                   ukbb_df$FU_END_AGE)
    }
}

#harmonize names
names(ukbb_df) <- remove_special(names(ukbb_df))


# join cause data to demographics
ukbb_df <- left_join(demographics, wide_ukbb_first, by="eid")

# for each cause_AGE column, for NA values (cause did not happen) substitute NA as the age at end of follow-up
for (j in 1:ncol(ukbb_df)) {
    colname <- names(ukbb_df)[j]
    if (grepl("_AGE", colname)) {
        ukbb_df[,colname] <- ifelse(!is.na(ukbb_df[,colname]),
                                    ukbb_df[,colname],
                                   ukbb_df$FU_END_AGE)
    }
}

#harmonize names
names(ukbb_df) <- remove_special(names(ukbb_df))

save.image(file = "ukbb_workspace.RData", compress=FALSE)

load("ukbb_workspace.RData")

ukbb_df <- ukbb_df %>% filter(eid %in% ids_ukbb)
ukbb_df <- ukbb_df %>% filter(eid %in% ids_ukbb)

# join genotypes to ukbb_df
#ukbb_df <- left_join(ukbb_df, geno_ukbb, by="eid")

# remove ids without genotypes
#ids_ukbb <- ids_ukbb[ids_ukbb %in% geno_ukbb$eid]
#ukbb_df <- ukbb_df %>% filter(eid %in% ids_ukbb)


saveRDS(ukbb_df, "/home/jsjukara/ukbb/first_endpoints_ukbb.rds", compress=FALSE)
saveRDS(long_ukbb, "/home/jsjukara/ukbb/longitudinal_endpoints_ukbb.rds", compress=FALSE)


#save.image(file = "ukbb_workspace.RData", compress=FALSE)
## Compare incidences in FinnGen vs UKBB
fg_df <- readRDS("/home/jsjukara/ukbb/fg_df.rds")
endpoint_map <- readRDS("/home/jsjukara/ukbb/endpoint_map.rds")

names(ukbb_df) <- gsub(" ", "_", names(ukbb_df))
map$cause <- gsub(" ", "_", map$cause)

for (j in 1:ncol(fg_df)) {
    colname <- names(fg_df)[j]
    for (i in 1:nrow(map)) {
        if (colname == map$fg_endpoints[i]) {
            print(paste(names(fg_df)[j], map$cause[i]))
            names(fg_df)[j] <- map$cause[i]
            
        } else if (colname == paste(map$fg_endpoints[i], "_AGE", sep="")) {
            print(paste(names(fg_df)[j], map$cause[i]))
            names(fg_df)[j] <- paste(map$cause[i], "_AGE", sep="")
        }
    }
}



get_incidences <- function(cause, cohort, data) {
  # define data frame by age groups that is to contain probability estimates in age groups for outcome and death
  age_df <- data.frame(age_name = c("Under 5", "5 to 9", "10 to 14", "15 to 19", "20 to 24", "25 to 29", "30 to 34",
                                    "35 to 39", "40 to 44", "45 to 49", "50 to 54", "55 to 59", "60 to 64", "65 to 69",
                                    "70 to 74", "75 to 79", "80 to 84"),
                       age_low = seq(0,80,5),
                       age_high = seq(5,85,5),
                       n=NA,
                       n_outcome=NA,
                       n_death=NA,
                       person_time=NA,
                       outcome_incidence=NA)
    age_df$age_name <- as.character(age_df$age_name)
    endpoint_age <- paste(endpoint, "_AGE", sep="")
        
    # loop over age groups
    for (i in 1:nrow(age_df)) {

        temp <- fg_df %>% filter(FU_START_AGE <= age_df$age_low[i], # excl. if FU start after start of interval
                           FU_END_AGE >= age_df$age_low[i], #excl if FU end before start of interval
                           #get(endpoint) == 0 | get(endpoint_age) >= age_df$age_low[i], # excl. if outcome before interval
                           get(endpoint_age) >= age_df$age_low[i])#, # excl. if outcome before interval
                           #DEATH == 0 | (DEATH == 1 & AGE_AT_DEATH_OR_NOW >= age_df$age_low[i]))#, # excl. if died before start of interval
                           #DEATH == 0 & FU_END_AGE < age_df$age_high[i]) #excl if did not die and FU ended before high

        # set outcome positive if outcome in interval
        temp <- temp %>% mutate(outcome = case_when(get(endpoint) == 1 &
                                                   get(endpoint_age) > age_df$age_low[i] &
                                                   get(endpoint_age) <= age_df$age_high[i] ~ 1,
                                                   TRUE ~ 0)
                                                   ) # hack to get else condition working

        # set futime

        temp <- temp %>% mutate(futime = case_when(outcome == 1 ~ get(endpoint_age),
                                                  TRUE ~ age_df$age_high[i])
        )

        # set death as positive in interval
        temp <- temp %>% mutate(death = case_when(DEATH == 1 &
                                                  DEATH_AGE >= age_df$age_low[i] &
                                                   DEATH_AGE < age_df$age_high[i] ~ 1,
                                                   TRUE ~ 0)
                                                   ) # hack to get else condition working

        temp$person_time <- temp$futime - age_df$age_low[i]

        age_df$person_time[i] <- sum(temp$person_time)
        age_df$n_death[i] = sum(temp$death)
        age_df$n[i] = nrow(temp)
        age_df$n_outcome[i] = sum(temp$outcome)

        }

    age_df$outcome_incidence <- age_df$n_outcome/age_df$person_time

    return(age_df)
}


plot_incidences <- function(cause) {
    coh <- "all"
    met <- "coxph"


    age_df1 <- get_incidences(endpoint=cause,
                                 method=met,
                                 cohort=coh,
                                 fg_df=ukbb_df)
    age_df1$source <- "ukbb"

    age_df2 <- get_incidences(endpoint=cause,
                                 method=met,
                                 cohort=coh,
                                 fg_df=fg_df)

    age_df2$source <- "finngen"

    age_df <- rbind(age_df1, age_df2)

    return(age_df %>%
        filter(age_low >=25, age_low<=75) %>% 
        ggplot(aes(x=age_low, y=outcome_incidence, colour=source)) +
        geom_point() +
        geom_line() +
        labs(title=cause))
}

# get temporary list of causes that are shared between the two data sets
causes <- intersect(names(ukbb_df), names(fg_df))
causes <- unique(gsub("_AGE", "", causes))
causes <- causes[-c(1:6)]
causes

plot(plot_incidences("Alzheimer\'s_disease_and_other_dementias"))
plot(plot_incidences("Ischemic_heart_disease"))
plot(plot_incidences("Asthma"))
plot(plot_incidences("Migraine"))
plot(plot_incidences("Chronic_obstructive_pulmonary_disease"))
plot(plot_incidences("Ischemic_stroke"))
plot(plot_incidences("Colon_and_rectum_cancer"))

for (cause in causes) {
    plot(plot_incidences(cause))
}

fg_df <- readRDS("/home/jsjukara/ukbb/fg_df_incidences.rds")
head(fg_df)

names(ukbb_df) <- remove_special(names(ukbb_df))

get_incidences <- function(cause, data) {
  # define data frame by age groups that is to contain probability estimates in age groups for outcome and death
  age_df <- data.frame(age_name = c("Under 5", "5 to 9", "10 to 14", "15 to 19", "20 to 24", "25 to 29", "30 to 34",
                                    "35 to 39", "40 to 44", "45 to 49", "50 to 54", "55 to 59", "60 to 64", "65 to 69",
                                    "70 to 74", "75 to 79", "80 to 84"),
                       age_low = seq(0,80,5),
                       age_high = seq(5,85,5),
                       person_time=NA,
                       n_outcome=NA,
                       incidence=NA,
                       incidence_loci=NA,
                       incidence_hici=NA,
                       cause_name=cause,
                      stringsAsFactors=FALSE)
    cause_age <- paste(cause, "_AGE", sep="")
    

    for (i in 1:nrow(age_df)) {
        temp <- data %>% filter(FU_START_AGE <= age_df$age_low[i], # excl. if FU start after start of interval
                           FU_END_AGE >= age_df$age_low[i], #excl if FU end before start of interval
                           #get(cause) == 0 | get(cause_age) >= age_df$age_low[i], # excl. if outcome before interval
                           get(cause_age) >= age_df$age_low[i])#, # excl. if outcome before interval
                           #DEATH == 0 | (DEATH == 1 & AGE_AT_DEATH_OR_NOW >= age_df$age_low[i]))#, # excl. if died before start of interval
                           #DEATH == 0 & FU_END_AGE < age_df$age_high[i]) #excl if did not die and FU ended before high

        # set outcome positive if outcome in interval
        temp <- temp %>% mutate(outcome = case_when(get(cause) == 1 &
                                                   get(cause_age) > age_df$age_low[i] &
                                                   get(cause_age) <= age_df$age_high[i] ~ 1,
                                                   TRUE ~ 0)
                                                   ) # hack to get else condition working

        # set futime
        temp <- temp %>% mutate(futime = case_when(outcome == 1 ~ get(cause_age),
                                                  TRUE ~ age_df$age_high[i])
        )

        temp$person_time <- temp$futime - age_df$age_low[i]

        age_df$person_time[i] <- sum(temp$person_time)
        age_df$n_outcome[i] = sum(temp$outcome)

        }

    age_df$incidence <- age_df$n_outcome/age_df$person_time
    
    se <- sqrt((1-age_df$incidence)/age_df$n_outcome)
    
    age_df$incidence_loci <- exp(log(age_df$incidence) - 1.96*se)
    age_df$incidence_hici <- exp(log(age_df$incidence) + 1.96*se)

    return(age_df)
}

# get causes both in fg_df and ukbb_df
common_causes <- map$cause[(map$cause %in% names(fg_df))]
common_causes[!(common_causes %in% names(ukbb_df))]

col_names <- c("location", "cause", "age_name", "age_low", "age_high", "incidence", "incidence_loci", "incidence_hici")
templist <- list(NULL)
for (cause in common_causes) {
    colnames <- c("FINNGENID", "FU_START_AGE","FU_END_AGE", cause, paste(cause, "_AGE", sep=""))
    templist[[cause]] <- get_incidences(cause, data=fg_df[,colnames])
}
fg_incidences <- bind_rows(templist, .id="cause")
fg_incidences$location <- "finngen"
fg_incidences <- fg_incidences[,col_names]

templist <- list(NULL)
for (cause in common_causes) {
    colnames <- c("eid", "FU_START_AGE","FU_END_AGE", cause, paste(cause, "_AGE", sep=""))
    templist[[cause]] <- get_incidences(cause, data=ukbb_df[,colnames])
}
ukbb_incidences <- bind_rows(templist, .id="cause")
ukbb_incidences$location <- "ukbb"
ukbb_incidences <- ukbb_incidences[,col_names]

#read GBD cause data
cause_data <- read.csv(stringsAsFactors=FALSE, "/home/jsjukara/Auxiliary data/IHME-GBD_2019_DATA_UK_FINLAND_16102020.csv")
cause_data$cause <- remove_special(cause_data$cause)

ages_vec <- c("Under 5", "5 to 9", "10 to 14", "15 to 19", "20 to 24", "25 to 29", "30 to 34",
                                    "35 to 39", "40 to 44", "45 to 49", "50 to 54", "55 to 59", "60 to 64", "65 to 69",
                                    "70 to 74", "75 to 79", "80 to 84")

gbd_incidences_cause <- cause_data %>% filter(year == 2019,
          metric == "Rate",
          sex=="Both",
          measure == "Incidence",
          age %in% ages_vec) %>%
    select(c("location", "measure", "sex", "age", "cause", "val", "upper", "lower")) %>%
    mutate(incidence=val / 100000,
          incidence_loci=lower / 100000,
          incidence_hici=upper / 100000,
          age_name=age)

gbd_incidences_finland <- gbd_incidences_cause %>% filter(location=="Finland")
gbd_incidences_uk <- gbd_incidences_cause %>% filter(location=="United Kingdom")

head(ukbb_incidences)
head(fg_incidences)

gbd_incidences <- ukbb_incidences[,1:5]
gbd_incidences$location <- "United Kingdom"
temp <- fg_incidences[,1:5]
temp$location <- "Finland"
gbd_incidences <- rbind(gbd_incidences, temp)

dim(gbd_incidences)
gbd_incidences <- left_join(gbd_incidences,
                            gbd_incidences_cause[c("location", "cause", "age_name", "incidence", "incidence_loci", "incidence_hici")])
dim(gbd_incidences)
head(gbd_incidences)

incidences_df <- rbind(fg_incidences, ukbb_incidences)
incidences_df <- rbind(incidences_df, gbd_incidences)
dim(incidences_df)
head(incidences_df)

length(unique(incidences_df$cause))

options(repr.plot.width=16, repr.plot.height=50)
incidences_df %>%
    #filter(cause == "Stomach_cancer") %>%
    ggplot(aes(x=age_low+2.5, y=incidence*10^5, colour=location)) +
    geom_point(size=2) +
    geom_line(size=1) +
    geom_errorbar(aes(ymin=incidence_loci*10^5, ymax=incidence_hici*10^5)) +
    facet_wrap(~cause, scales="free", ncol=4) +
    labs(x="Age", y="No. cases / 100 000") +
    theme_bw() +
    scale_color_manual(values=c("#8080ff", "#0000ff", "#cc0000", "#ff6666"),
                      name = "",
                      labels=c("GBD Fin", "FinnGen", "UKBB", "GBD UK")) +
    #scale_color_brewer(#palette = "Dark2",
    #                  values = c("#b3b3ff", "#0000ff", "#cc0000", "#ff9999"),
    #                  name = "",
    #                      labels=c("GBD Fin", "FinnGen", "UKBB", "GBD UK")) +
    theme(legend.position = "top",
         panel.grid.minor = element_blank(),
         legend.key.width = unit(1.5,"cm")) +
    guides(colour = guide_legend(override.aes = list(size = 2)))

options(repr.plot.width=16, repr.plot.height=50)
incidences_df %>%
    #filter(cause == "Stomach_cancer") %>%
    ggplot(aes(x=age_low+2.5, y=incidence*10^5, colour=location)) +
    geom_point(size=2) +
    geom_line(size=1) +
    #geom_errorbar(aes(ymin=incidence_loci*10^5, ymax=incidence_hici*10^5)) +
    facet_wrap(~cause, scales="free", ncol=4) +
    labs(x="Age", y="No. cases / 100 000") +
    theme_bw() +
    scale_color_manual(values=c("#8080ff", "#0000ff", "#cc0000", "#ff6666"),
                      name = "",
                      labels=c("GBD Fin", "FinnGen", "UKBB", "GBD UK")) +
    #scale_color_brewer(#palette = "Dark2",
    #                  values = c("#b3b3ff", "#0000ff", "#cc0000", "#ff9999"),
    #                  name = "",
    #                      labels=c("GBD Fin", "FinnGen", "UKBB", "GBD UK")) +
    theme(legend.position = "top",
         panel.grid.minor = element_blank(),
         legend.key.width = unit(1.5,"cm")) +
    guides(colour = guide_legend(override.aes = list(size = 2)))

endpoints <- c("MIGRAINE_TRIPTAN", "C3_PROSTATE", "Schizophrenia", "J10_COPD", "T1D_WIDE")
templist <- list(NULL)
for (ep in endpoints) {
    colnames <- c("FINNGENID", "FU_START_AGE","FU_END_AGE", ep, paste(ep, "_AGE", sep=""))
    templist[[ep]] <- get_incidences(ep, data=fg_df[,colnames])
}

tempdf <- bind_rows(templist, .id="endpoint_name")

endpoint_map$endpoint_name <- endpoint_map$endpoint

tempdf <- left_join(tempdf, endpoint_map[,c("endpoint_name", "cause_name")])

tempdf <- left_join(tempdf, gbd_subset[,c("cause_name", "age_name", "outcome_incidence")])

translate_df <- data.frame(endpoint_name=c("MIGRAINE_TRIPTAN", "C3_PROSTATE", "Schizophrenia", "J10_COPD", "T1D_WIDE"),
                          label_name=c("Migraine", "Prostate cancer", "Schizophrenia", "COPD", "T1DM"))

tempdf <- left_join(tempdf, translate_df)

#keep_endpoints <- c("C3_BREAST", "C3_PROSTATE", "I9_AF", "J10_ASTHMA",
#                    "I9_ISCHHEART", "I9_STR_SAH", "J10_COPD", "KRA_PSY_DEMENTIA",
#                   "MIGRAINE_TRIPTAN", "N14_PROSTHYPERPLA", "Schizophrenia", "T2D")
temp <- tempdf %>%
    pivot_longer(names_to="incidence_source",
                 values_to="incidence",
                 c("outcome_incidence", "incidence")) %>%
    mutate(incidence_source = case_when(incidence_source == "outcome_incidence" ~ "GBD incidence",
                                        TRUE ~ "FinnGen incidence")) %>%
    mutate(incidence_source = as.factor(incidence_source)) %>%
    mutate(incidence = incidence*10^5)
options(repr.plot.width=10, repr.plot.height=7.5)

temp %>%
    filter(endpoint_name != "T1D_WIDE") %>%
    ggplot(aes(x=age_low+2.5, y=incidence, colour=incidence_source)) +
    geom_point(size=3) +
    geom_line(size=1.5) +
    facet_wrap(~label_name, scales="free", ncol=2) +
    labs(x="Age", y="No. cases / 100 000") +
    theme_bw(base_size=20) +
    scale_color_brewer(palette = "Dark2",
                      name = "",
                          labels=c("FinnGen", "GBD")) +
    theme(legend.position = "top",
         panel.grid.minor = element_blank(),
         legend.key.width = unit(1.5,"cm")) +
    guides(colour = guide_legend(override.aes = list(size = 2)))

names(ukbb_df) <- remove_special(names(ukbb_df))

endpoints <- c("Migraine", "Prostate_cancer", "Schizophrenia", "Chronic_obstructive_pulmonary_disease")
templist <- list(NULL)
for (ep in endpoints) {
    colnames <- c("eid", "FU_START_AGE","FU_END_AGE", ep, paste(ep, "_AGE", sep=""))
    templist[[ep]] <- get_incidences(ep, data=ukbb_df[,colnames])
}

tempdf_ukbb <- bind_rows(templist, .id="endpoint_name")

tempdf_ukbb$cause_name <- tempdf_ukbb$endpoint_name

tempdf_ukbb$incidence_ukbb <- tempdf_ukbb$incidence

tempdf1 <- left_join(tempdf, tempdf_ukbb[c("cause_name", "age_name", "incidence_ukbb")])

temp <- tempdf1 %>%
    pivot_longer(names_to="incidence_source",
                 values_to="incidence",
                 c("outcome_incidence", "incidence", "incidence_ukbb")) %>%
    mutate(incidence_source = case_when(incidence_source == "outcome_incidence" ~ "GBD incidence",
                                        incidence_source == "incidence_ukbb" ~ "UKBB incidence",
                                        TRUE ~ "FinnGen incidence")) %>%
    mutate(incidence_source = as.factor(incidence_source)) %>%
    mutate(incidence = incidence*10^5)

options(repr.plot.width=10, repr.plot.height=7.5)
temp %>%
    filter(endpoint_name != "T1D_WIDE") %>%
    ggplot(aes(x=age_low+2.5, y=incidence, colour=incidence_source)) +
    geom_point(size=3) +
    geom_line(size=1.5) +
    facet_wrap(~label_name, scales="free", ncol=2) +
    labs(x="Age", y="No. cases / 100 000") +
    theme_bw(base_size=20) +
    scale_color_brewer(palette = "Dark2",
                      name = "",
                          labels=c("FinnGen", "GBD", "UKBB")) +
    theme(legend.position = "top",
         panel.grid.minor = element_blank(),
         legend.key.width = unit(1.5,"cm")) +
    guides(colour = guide_legend(override.aes = list(size = 2)))