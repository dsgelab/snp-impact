setwd("/home/jsjukara/")
library(data.table)
library(tidyverse)
library(survival)

#load genotypes
geno <- readRDS("/home/jsjukara/finngen_data/preprocessed_genotypes_fg.rds")
#geno_ukbb <- readRDS("/home/jsjukara/ukbb/preprocessed_genotypes_ukbb.rds")

#load small and big endpoint file
end <- fread("zcat /home/jsjukara/finngen_data/finngen_R6_cov_pheno_V1.txt.gz", data.table=FALSE)
endbig <- fread("zcat /home/jsjukara/finngen_data/finngen_R6_v2_endpoint.gz", data.table=FALSE)

age_name <- c("Under 5", "5 to 9", "10 to 14", "15 to 19", "20 to 24", "25 to 29", "30 to 34", "35 to 39", "40 to 44",
"45 to 49", "50 to 54", "55 to 59", "60 to 64", "65 to 69", "70 to 74", "75 to 79", "80 to 84")

endpoint_map <- as.data.frame(fread("Auxiliary data/gbd_to_finngen_ontology_manual_7_2020.txt"))
    
endpoint_map$fg_endpoints <- ifelse(endpoint_map$fg_endpoints=="", NA, endpoint_map$fg_endpoints)
endpoint_map <- endpoint_map %>%
    filter(!is.na(fg_endpoints) & include==1) %>%
    select(Cause, fg_endpoints)

modify_strings <- function(varname, data) {
    string_vec <- data[,varname]
    # substitute underscores into spaces
    string_vec <- gsub(" ", "_", string_vec)
    # remove apostrophes
    string_vec <- gsub("\\'", "", string_vec)
    # remove commas
    string_vec <- gsub(",", "", string_vec)
    data[,varname] <- string_vec
    return(data)
}

endpoint_map <- modify_strings("Cause", endpoint_map)

# get vector of unique FinnGen endpoints in endpoint_map
fg_endpoints_vec <- c(NULL)
for (i in 1:nrow(endpoint_map)) {
    vec <- unlist(strsplit(endpoint_map$fg_endpoints[i], split="\\|"))
    fg_endpoints_vec <- c(fg_endpoints_vec,vec)
}


# get endpoints without a match in FinnGen data, print them and remove them from endpoint_map
fg_endpoints_vec <- fg_endpoints_vec[(fg_endpoints_vec %in% names(endbig))]
endpoint_map$remove <- 0
for (i in 1:nrow(endpoint_map)) {
    vec <- unlist(strsplit(endpoint_map$fg_endpoints[i], split="\\|"))
    if (sum(vec %in% names(endbig)) == 0) {
        endpoint_map$remove[i] <- 1
    }
}
print(endpoint_map[endpoint_map$remove==1,])
endpoint_map <- endpoint_map[endpoint_map$remove==0,]

#set names for columns to keep from endpoint data frames, keep only them.
pcs <- paste("PC",1:20, sep="")
keepcols_end <- c("FINNGENID", "BL_AGE", "BL_YEAR", "SEX", "AGE_AT_DEATH_OR_NOW", "regionofbirth", "cohort", pcs)

fg_endpoints_vec <- fg_endpoints_vec[(fg_endpoints_vec %in% names(endbig))]
keepcols_endbig <- c("FINNGENID","FU_END_AGE", fg_endpoints_vec, paste(fg_endpoints_vec, "_AGE", sep=""),
                     'DEATH','DEATH_AGE','DEATH_YEAR')

#keepcols_endbig <- keepcols_endbig[(keepcols_endbig %in% names(endbig))]
fg_df <- endbig[,keepcols_endbig]
end <- end[,keepcols_end]
fg_df <- left_join(fg_df, end, by="FINNGENID")
fg_df <- fg_df[!is.na(fg_df$AGE_AT_DEATH_OR_NOW),]

fg_df$YEAR_OF_BIRTH <- fg_df$BL_YEAR - fg_df$BL_AGE
fg_df$FU_START_AGE <- 1972 - fg_df$YEAR_OF_BIRTH
fg_df$FEMALE <- ifelse(fg_df$SEX=="female", 1, 0)

fg_df$FU_START_AGE <- ifelse(fg_df$FU_START_AGE<0, 0, fg_df$FU_START_AGE)

# make NA endpoints 0 (they are controls with comorbidities)
for (endpoint in fg_endpoints_vec) {
    fg_df[,endpoint][is.na(fg_df[,endpoint])] <- 0
}

# join genotypes
fg_df <- left_join(fg_df, geno, by="FINNGENID")
rsids <- names(geno)[2:length(geno)]

# get genotype frequencies from nationally representative population samples in FinnGen
popsample_cohortnames <- c("THL BIOBANK FINRISK 1992",
                          "THL BIOBANK FINRISK 1997",
                          "THL BIOBANK FINRISK 2002",
                          "THL BIOBANK FINRISK 2007",
                          "THL BIOBANK FINRISK 2012",
                          "THL BIOBANK FINHEALTH 2017",
                          "THL BIOBANK HEALTH 2000",
                          "THL BIOBANK HEALTH 2011")
popsample_fg_df <- fg_df[fg_df$cohort %in% popsample_cohortnames,]
genotype_freqs <- data.frame(rsid = rsids, p0=NA, p1=NA, p2=NA, stringsAsFactors=FALSE)
for (i in 1:length(rsids)) {
    freqs <- as.vector(table(popsample_fg_df[,rsids[i]])/nrow(popsample_fg_df))
    if (length(freqs)==3) {
        genotype_freqs[i,2:4] <- freqs
        
    } else if (length(freqs)==2) {
        genotype_freqs[i,2:4] <- c(freqs,0)
    } else {
        print("error")
    }    
}


# Load Finnish DALYs
dalys <- read.csv(stringsAsFactors=FALSE, "Auxiliary data/IHME-GBD_2017_DATA-OVERALL_DALY.csv")
dalys <- modify_strings("cause_name", dalys)

# Load and filter GBD life tables -----------------------------------------
lt <- read.csv(stringsAsFactors=FALSE, "Auxiliary data/IHME_GBD_2017_ABRIDGED_LIFE_TABLES_2017_Y2018M11D08.CSV")
lt <- lt %>%
  filter(location_name == "Finland" &
         measure_name == "Probability of death" &
         sex_name == "Both")

names(lt)[c(6,12)] <- c("age_name", "prob_death")

# Load population data ----------------------------------------------------
pop <- read.csv(stringsAsFactors=FALSE, "Auxiliary data/IHME_GBD_2017_POP_2015_2017_Y2018M11D08.CSV")

pop <- pop %>%
  filter(location_name == "Finland" &
           year_id == 2017 &
         sex_name == "Both")

names(pop)[c(6,11)] <- c("age_name", "population")

# Get incidences of outcome & all cause mortality from GBD -----------------------------------------------
gbdfin <- read.csv(stringsAsFactors=FALSE, "Auxiliary data/IHME-GBD_2017_DATA_FINLAND/IHME-GBD_2017_DATA_FINLAND.csv")
gbdfin <- modify_strings("cause_name", gbdfin)

#df for all cause mortality
acm <- gbdfin %>% filter(year == 2017 &
                           cause_name == "All_causes" &
                           metric_name == "Number" &
                           sex_name == "Both" &
                           measure_name %in% c("Deaths"))

acm <- left_join(acm, pop[,c("age_name", "sex_name", "population")], by=c("age_name", "sex_name"))
acm$acm <- acm$val/acm$population # calculate all-cause mortality rate


#load population dalys
dalys_both <- dalys %>% filter(sex_name=="Both")

rm(end)
rm(endbig)

#gbd <- read.csv(stringsAsFactors=FALSE, "Auxiliary data/IHME-GBD_FINLAND_UK_2016_2017.csv")

# join dalys to endpoint_map
temp <- dalys_both
names(temp)[10] <- "Cause"
endpoint_map <- left_join(endpoint_map, temp[,c("Cause", "val")])
names(endpoint_map)[grep("val", names(endpoint_map))] <- "population_dalys"

# remove endpoints with no DALYs
endpoint_map <- endpoint_map[!is.na(endpoint_map$population_dalys),]
# remove endpoints with no incidences
cause_names_with_incidence <- unique(gbdfin$cause_name[gbdfin$measure_name=="Incidence"])
# get cause names without incidence in GBD
causes_wo_incidence <- endpoint_map$Cause[!(endpoint_map$Cause %in% cause_names_with_incidence)]

## go from "ENDPOINT1|ENDPOINT2" notation to generate custom composite endpoints in FinnGen
# initialize vectors to contain endpoint names
endpoints_gbd <- rep(NA,nrow(endpoint_map))
endpoints <- rep(NA,nrow(endpoint_map))
# loop through rows in endpoint_map
for (i in 1:nrow(endpoint_map)) {
    # if current row has a composite outcome of FG endpoints, generate new endpoint
    if(grepl("\\|", endpoint_map$fg_endpoints[i])) {
        # initialize new composite endpoint names
        ep_name <- endpoint_map$Cause[i]
        ep_age_name <- paste(endpoint_map$Cause[i], "_AGE", sep="")
        fg_df[,ep_name] <- NA
        fg_df[,ep_age_name] <- NA
        # split multiple finngen endpoints into a vector of endpoints
        vec <- unlist(strsplit(endpoint_map$fg_endpoints[i], "\\|"))
        vec_age <- paste(vec, "_AGE", sep="")
        # if all listed endpoints not in FG, remove those that aren't
        if (sum(!(vec %in% names(fg_df)))!=0) {
            print("Following missing from FinnGen endpoints, deleting from set:")
            print(vec[!(vec %in% names(fg_df))])
            
            vec_age <- vec_age[(vec %in% names(fg_df))]
            vec <- vec[(vec %in% names(fg_df))]
            # if no endpoints left, skip to next endpoint
            if (length(vec)==0) {
                print(paste("Skipping endpoint:",ep_name))
                next
            }
        }
        if (length(vec)==1) {
            fg_df[,ep_name] <- fg_df[,vec]
            fg_df[,ep_age_name] <- fg_df[,vec_age]
        } else {
            fg_df[,ep_name] <- ifelse(rowSums(fg_df[,vec])==0, 0, 1)
            fg_df[,ep_age_name] <- ifelse(fg_df[,ep_name]==1, apply(fg_df[,vec_age], 1, FUN=min), fg_df$AGE_AT_DEATH_OR_NOW)
        }
        endpoints_gbd[i] <- ep_name
        endpoints[i] <- ep_name
    } else { # else if not composite outcome, just set the endpoint names
        endpoints_gbd[i] <- endpoint_map$Cause[i]
        endpoints[i] <- endpoint_map$fg_endpoints[i]
    }
}

# remove gbd endpoints that are not in gbdfin
endpoints <- endpoints[(endpoints_gbd %in% unique(gbdfin$cause_name))]
endpoints_gbd <- endpoints_gbd[(endpoints_gbd %in% unique(gbdfin$cause_name))]
endpoint_map <- endpoint_map[(endpoint_map$Cause %in% endpoints_gbd),]
endpoint_map$endpoint <- endpoints
endpoint_map$endpoint_gbd <- endpoints_gbd

endpoint_map <- endpoint_map[endpoint_map$Cause %in% cause_names_with_incidence,]

# define number of outcomes in FinnGen to endpoint_map for each endpoint
endpoint_map$n_outcomes <- NA
for (i in 1:nrow(endpoint_map)) {
    endpoint_map$n_outcomes[i] <- sum(fg_df[,endpoint_map$endpoint[i]])
}

# remove outcomes that have too few endpoints or account for too few population DALYs
min_endpoints <- 500
min_dalys <- 0
endpoint_map <- endpoint_map[endpoint_map$n_outcomes > min_endpoints,]
endpoint_map <- endpoint_map[endpoint_map$population_dalys > min_dalys,]

# write RDS files of endpoints and endpoint map
#saveRDS(fg_df, "/home/jsjukara/ukbb/fg_df.rds")
#saveRDS(endpoint_map, "/home/jsjukara/ukbb/endpoint_map.rds")

#removed_endpoints <- c("C3_CANCER", "L12_ACNE", "M13_LOWBACKPAINORANDSCIATICA",
#                       "Q17_CONGEN_MALFO_DEFORMAT_CHROMOSOMAL_ABNORMALITI", "CONGEN_HEART_ARTER",
#                      "G6_HEADACHE", "KRA_PSY_ALCOH", "KRA_PSY_SUBSTANCE", "F5_ANOREX", "F5_EATING", "G6_OTHNEU",
#                      "DIABETES_FG", "Pancreatitis", "I9_CVD", "I9_STR_EXH", "C3_OTHER_SKIN", "Cardiomyopathy_and_myocarditis",
#                      "N14_FEMGENPROL", "Appendicitis", "Peptic_ulcer_disease", "Gastritis_and_duodenitis",
#                      "K11_REFLUX", "N14_ENDOMETRIOSIS", "I9_VHD", "C3_CORPUS_UTERI", "K11_ILEUS", "G6_EPLEPSY",
#                      "F5_BIPO")
#endpoint_map <- endpoint_map[!(endpoint_map$endpoint %in% removed_endpoints),]

# calculate minor allele frequencies
genotype_freqs$maf <- (1*genotype_freqs$p1+2*genotype_freqs$p2)/2
genotype_freqs$mafover1 <- ifelse(genotype_freqs$maf>0.01, 1, 0)
table(genotype_freqs$mafover1)

# filter maf over 1 variants
genotype_freqs <- genotype_freqs[genotype_freqs$mafover1==1,]

rsids <- genotype_freqs$rsid
# get number of variant-endpoint pairs
nrow(endpoint_map)*length(rsids)

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


# get vector of unique FinnGen endpoints in endpoint_map to use in getting finemapped snps
eps_for_snps <- c(NULL)
for (i in 1:nrow(endpoint_map)) {
    vec <- unlist(strsplit(endpoint_map$fg_endpoints[i], split="\\|"))
    eps_for_snps <- c(eps_for_snps,vec)
}

# eps_for_snps <- fg_endpoints_vec[fg_endpoints_vec %in% names(endbig)]
# get command for listing finemapping results for all endpoints for shell
paste("gsutil ls -r gs://finngen-production-library-green/finngen_R5/finngen_R5_analysis_data/finemap/ | grep -E \"", paste(unique(eps_for_snps), collapse="|"), "\" | grep \"SUSIE.snp\" > Directorylisting.txt", sep="")