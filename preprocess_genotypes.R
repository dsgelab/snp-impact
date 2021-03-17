setwd("/home/jsjukara/ukbb/")
library(tidyverse)
library(data.table)

# read ukbb genotypes
geno <- fread("zcat skari_snps_28072020.tsv.gz")
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
variant_map <- read.csv("/home/jsjukara/finemap/finemap_rsids_and_variant_ids_28072020.csv", stringsAsFactors=FALSE)
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
#table(variant_map$locus_hg19)[table(variant_map$locus_hg19)>1]

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

# remove duplicated rows
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
geno <- fread("/home/jsjukara/finngen_data/doc_jukarainen_SNPs_R6_25082020.tsv", data.table=FALSE)

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

geno_fg <- geno
rm(geno)

# make sure all alleles are the right way around in geno_fg_fg
loci <- gsub(":", "_", loci)

fg_loci_alleles <- data.frame(locus_hg38 = loci, variant_finngen = names(geno_fg)[1:(ncol(geno_fg)-1)], stringsAsFactors=FALSE)

fg_loci_alleles <- left_join(fg_loci_alleles, variant_map[,c("locus_hg38", "variant_hg38")], by="locus_hg38")
head(fg_loci_alleles)
sum(is.na(fg_loci_alleles$variant_hg38))
#fg_loci_alleles[fg_loci_alleles$variant_hg38 != fg_loci_alleles$variant_finngen,]

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
geno_ukbb[,"chr1_145804971_G_C"] <- -1*(geno_ukbb[,"chr1_145804971_G_C"])+2

# plot maf to ensure coding is same way
plot(genotype_freqs_ukbb$maf, genotype_freqs_fg$maf)

saveRDS(geno_fg, "/home/jsjukara/finngen_data/preprocessed_genotypes_fg.rds", compress=FALSE)
saveRDS(geno_ukbb, "/home/jsjukara/ukbb/preprocessed_genotypes_ukbb.rds", compress=FALSE)