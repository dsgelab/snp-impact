setwd("/home/jsjukara/")
library(data.table)
library(tidyverse)
library(survival)
library(tictoc)

setwd("/home/jsjukara/")
library(data.table)

library(tidyverse)
library(survival)g
library(tictoc)

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

get_age_df_se <- function(endpoint, rsid, method, cohort, fg_df) {
  # define data frame by age groups that is to contain probability estimates in age groups for outcome and death
  age_df <- data.frame(age_name = c("Under 5", "5 to 9", "10 to 14", "15 to 19", "20 to 24", "25 to 29", "30 to 34",
                                    "35 to 39", "40 to 44", "45 to 49", "50 to 54", "55 to 59", "60 to 64", "65 to 69",
                                    "70 to 74", "75 to 79", "80 to 84"),
                       age_low = seq(0,80,5),
                       age_high = seq(5,85,5),
                       loghr=NA,
                       loghr_se=NA,
                       hr01=NA,
                       hr01loci=NA,
                       hr01hici=NA,
                       hr02=NA,
                       n=NA,
                       n_outcome=NA,
                       n_death=NA,
                       person_time=NA,
                       outcome_incidence_fg=NA)
    age_df$age_name <- as.character(age_df$age_name)
    endpoint_age <- paste(endpoint, "_AGE", sep="")

    if (method == "coxph" | method == "crude_incidence") {
        
        # loop over age groups
        for (i in 1:nrow(age_df)) {
            if (cohort == "post1998") {
                temp <- fg_df %>% filter(1998 - YEAR_OF_BIRTH <=  age_df$age_low[i], # excl. if age at 1998 higher than start of interval
                                   #FU_START_AGE < age_df$age_low[i], 
                                   FU_END_AGE >= age_df$age_low[i], #excl if FU end before start of interval
                                   #get(endpoint) == 0 | get(endpoint_age) >= age_df$age_low[i], # excl. if outcome before interval,
                                   get(endpoint_age) >= age_df$age_low[i])#, # excl. if outcome before interval
                                   #DEATH == 0 | (DEATH == 1 & AGE_AT_DEATH_OR_NOW >= age_df$age_low[i]))#, # excl. if died before start of interval
                                   #DEATH == 0 & FU_END_AGE < age_df$age_high[i]) #excl if did not die and FU ended before high
            } else {
                temp <- fg_df %>% filter(FU_START_AGE <= age_df$age_low[i], # excl. if FU start after start of interval
                                   FU_END_AGE >= age_df$age_low[i], #excl if FU end before start of interval
                                   #get(endpoint) == 0 | get(endpoint_age) >= age_df$age_low[i], # excl. if outcome before interval
                                   get(endpoint_age) >= age_df$age_low[i])#, # excl. if outcome before interval
                                   #DEATH == 0 | (DEATH == 1 & AGE_AT_DEATH_OR_NOW >= age_df$age_low[i]))#, # excl. if died before start of interval
                                   #DEATH == 0 & FU_END_AGE < age_df$age_high[i]) #excl if did not die and FU ended before high
            }
            
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
            
            #temp %>% mutate(person_time = case_when(person_time > 5 ~ 5,
            #                                        !is.na(FINNGENID) ~ person_time))

            age_df$person_time[i] <- sum(temp$person_time)
            age_df$n_death[i] = sum(temp$death)
            age_df$n[i] = nrow(temp)
            age_df$n_outcome[i] = sum(temp$outcome)

            if (method == "coxph") {
              # predict outcome
              form <- formula(paste("Surv(futime, outcome) ~ ", rsid, " + ", covariates_formula, sep=""))
              fit_outcome <- tryCatch(suppressWarnings(coxph(formula = form,
                           data=temp)), error=function(err) NA)
              if (is.na(fit_outcome)[1]) {
                  next
              }
                
              #only if model is fit, store HRs
              if (!is.na(fit_outcome$coefficients[1])) {
                  # save HR
                  age_df$loghr[i] <- fit_outcome$coefficients[1]
                  age_df$loghr_se[i] <- summary(fit_outcome)$coefficients[1,3]

                  age_df$hr01[i] <- exp(age_df$loghr[i])
                  age_df$hr02[i] <- exp(age_df$loghr[i]*2)

                  # estimate 95% CI
                  #if (!is.na(age_df$hr01[i])) {
                  #    age_df$hr01loci[i] <- summary(fit_outcome)$conf.int[1,3]
                  #    age_df$hr01hici[i] <- summary(fit_outcome)$conf.int[1,4]
                  #}
                  age_df$loghr_se[i] <- summary(fit_outcome)$coefficients[1,3]
              }
              
            } else if (method == "crude_incidence") {
              # no. outcomes by genotype
              n0 <- sum(temp$outcome[temp[,rsid]==0])
              n1 <- sum(temp$outcome[temp[,rsid]==1])
              # no. person-time by genotype
              pt0 <- sum(temp$person_time[temp[,rsid]==0])
              pt1 <- sum(temp$person_time[temp[,rsid]==1])
              if (n0>0 & n1>0 & pt0>0 & pt1>0) {
                age_df$hr01[i] <- (n1/pt1)/(n0/pt0) # estimate incidence rate ratio
                # impute hr02 through additive model
                age_df$hr02[i] <- age_df$hr01[i]^2
                
              }

              # estimate 95% CI https://journals.plos.org/plosmedicine/article/file?type=supplementary&id=info:doi/10.1371/journal.pmed.1001270.s005
              #if (!is.na(age_df$hr01[i])) {
              #  se <- sqrt(1/n0+1/n1)
              #
              #  age_df$hr01loci[i] <- exp(log(age_df$hr01[i]) - 1.96*se)
              #  age_df$hr01hici[i] <- exp(log(age_df$hr01[i]) + 1.96*se)
              #}

            }

            # if missing HRs, impute 1
            #if (sum(is.na(age_df[,c("hr01", "hr02")]))!=0) {
                # print("Had to impute HRs to 1")
            #   age_df$hr01 <- ifelse(is.na(age_df$hr01), 1, age_df$hr01)
            #    age_df$hr02 <- ifelse(is.na(age_df$hr02), 1, age_df$hr02)
            #}

        }
    } else if (method == "coxph_single" & cohort == "all") {
        temp <- fg_df
        
        # set outcome positive if outcome in interval
        temp <- temp %>% mutate(outcome = get(endpoint))
        # set death as positive in interval
        temp <- temp %>% mutate(death = DEATH)

        #temp[,"FU_START_AGE"] <- ifelse(temp[,"FU_START_AGE"]<0, 0, temp[,"FU_START_AGE"])

        temp$person_time <- temp[,endpoint_age]

        age_df$person_time <- sum(temp$person_time)
        age_df$n_death = sum(temp$death)
        age_df$n = nrow(temp)
        age_df$n_outcome = sum(temp$outcome)


        # predict outcome
        form <- formula(paste("Surv(", endpoint_age, ", ", endpoint,") ~ ", rsid, " + ", covariates_formula, sep=""))
        fit_outcome <- tryCatch(suppressWarnings(coxph(formula = form,
                           data=temp)), error=function(err) NA)
        if (is.na(fit_outcome)[1]) {
          next
        }

        # save HR
        age_df$loghr <- fit_outcome$coefficients[1]
        age_df$loghr_se <- summary(fit_outcome)$coefficients[1,3]

        age_df$hr01 <- exp(age_df$loghr)
        age_df$hr02 <- exp(age_df$loghr*2)

        # save HRs
        #age_df$hr01 <- exp(fit_outcome$coefficients[1])
        #age_df$hr01loci <- summary(fit_outcome)$conf.int[1,3]
        #age_df$hr01hici <- summary(fit_outcome)$conf.int[1,4]

        #age_df$hr02 <- exp(fit_outcome$coefficients[1]*2)
    }
    age_df$raw_prob_death <- age_df$n_death / age_df$n
    age_df$raw_prob_outcome <- age_df$n_outcome / age_df$n
    age_df$outcome_incidence_fg <- age_df$n_outcome/age_df$person_time

    return(age_df)
}

#define function to output lifetime risk given age_df, version for incidence method
get_ltr_incidence <- function(outcome_name,
                              age_df,
                              gbdfin,
                              rsid,
                              genotype_freqs) {

    # get incidences by multiplying GBD incidence by HR

    # get proportions of genotypes
    props <- as.vector(unlist(genotype_freqs[genotype_freqs$rsid==rsid,2:4]))
    
    age_df$sex_name <- "Both"
    
    # get GBD data for outcome, incidence, outcome deaths and prevalence
    gbd_outcome <- gbdfin %>% filter(cause_name == outcome_name)

    names_join_from_gbd <- c("age_name",
      "sex_name",
      "population",
      #"outcome_deaths", 
      "outcome_incidence",
      "outcome_mortality",
      "prevalence",
      "acm")

    age_df <- left_join(age_df, gbd_outcome[,names_join_from_gbd], by=c("age_name", "sex_name"))
    
    # estimate incidence attributable to 0, 1, 2 alleles
    age_df$i0 <- (age_df$outcome_incidence*age_df$population) / (props[1] * age_df$population + age_df$hr01 * props[2] * age_df$population + age_df$hr02 * props[3] * age_df$population)
    age_df$i1 <- age_df$i0 * age_df$hr01
    age_df$i2 <- age_df$i0 * age_df$hr02
    
    # Estimate lifetime risk of outcome using incidence method -------------------------------------
    ltrvec <- rep(NA,3)
    for (dosage in 0:2) {
      tab <- data.frame(age_name = age_df$age_name, # age group name
                      l = NA, # population of hypothetical cohort
                      m = age_df$acm, # death rate (per resident)
                      mo = age_df$outcome_mortality, #death rate from outcome (per resident)
                      rp = age_df[,paste("i", dosage, sep="")], # incidence rate of outcome in population
                      a = NA, # no. new cases within interval
                      d = NA, # no. non-outcome related deaths among outcome free individuals
                      prevalence = age_df$prevalence) # prevalence proportion in population

      #tab$gp = 1 - exp(-5*tab$rp) #prob of outcome in total population
      #tab$g = tab$gp *) #prob of outcome in outcome-free population
      #tab$r = (-1/5) * log(1 - tab$g) # incidence of outcome in outcome free population

      tab$r = tab$rp * (1 / (1 - tab$prevalence)) # incidence of outcome in outcome free population, new method from stroke paper

      tab$l[1] <- 10*10^6 #set 10 million as population of hypothetical cohort
      
      # "simulate" population at each timepoint by removing outcome-free deaths and incident outcomes
      for (i in 1:nrow(tab)) {
          if (i >= 2) {
            tab$l[i] <- tab$l[i-1] * exp(-5*( (tab$m[i-1] - tab$mo[i-1]) + tab$r[i-1]))
          }
      }
      # calculate number of new cases in interval
      tab$a <- tab$l * (1 - exp(-5*(tab$m + tab$r - tab$mo))) * (tab$r / (tab$m + tab$r - tab$mo))
      # calculate number deaths in interval
      tab$d <- tab$l * (1 - exp(-5*(tab$m + tab$r - tab$mo))) * ((tab$m - tab$mo) / (tab$m + tab$r - tab$mo))

      # Calculate a and d for last open ended age interval
      ind <- nrow(tab)
      tab$a[ind] <- tab$l[ind] * (tab$r[ind]/(tab$r[ind] + tab$m[ind] - tab$mo[ind]))
      tab$d[ind] <- tab$l[ind] * ((tab$m[ind] - tab$mo[ind])/(tab$r[ind] + tab$m[ind] - tab$mo[ind]))

      # Calculate lifetime risk
      ltrvec[dosage+1] <- sum(tab$a)/tab$l[1]
    }
    return(ltrvec)
}

covariates_formula <- paste(" + FEMALE + ", paste(pcs, collapse=" + "), sep="")
methods <- c("coxph", "coxph_single", "coxph_bayes6")# "crude_incidence")
cohorts <- c("all")
ltr_df <- expand.grid(rsid=rsids,
            endpoint_name=endpoint_map$endpoint,
            cohort=cohorts,
            method=methods,
            stringsAsFactors=FALSE)
temp <- endpoint_map
temp$endpoint_name <- temp$endpoint
temp$cause_name <- temp$endpoint_gbd
ltr_df <- left_join(ltr_df, temp[,c("endpoint_name", "cause_name")], by="endpoint_name")
#ltr_df$cause_name <- endpoint_map$endpoint_gbd
ltr_df$ltr0 <- NA
ltr_df$ltr1 <- NA
ltr_df$ltr2 <- NA
ltr_df <- left_join(ltr_df, dalys_both[,c("cause_name", "val")], by="cause_name")
names(ltr_df)[grep("val", names(ltr_df))] <- "population_dalys"

ltr_df <- ltr_df[!(ltr_df$cohort!="all" & ltr_df$method=="coxph_single"),]

# gsub names to harmonize
ltr_df$cause_name <- gsub("\\-", "_", ltr_df$cause_name)
ltr_df$cause_name <- gsub("/", "_", ltr_df$cause_name)
ltr_df$endpoint_name <- gsub("\\-", "_", ltr_df$endpoint_name)
ltr_df$endpoint_name <- gsub("/", "_", ltr_df$endpoint_name)

names(fg_df) <- gsub("\\-", "_", names(fg_df))
names(fg_df) <- gsub("/", "_", names(fg_df))

# shuffle rows to ensure even distribution of parallelized computations later on
set.seed(1)
ltr_df <- ltr_df[sample(1:nrow(ltr_df)),]

# initialize list to contain all age_dfs
age_df_list <- list()

#make negative ages 0
names_age <- paste(unique(ltr_df$endpoint_name), "_AGE", sep="")
for (col in names_age) {
    fg_df[,col] <- ifelse(fg_df[,col]<0, 0, fg_df[,col])
}

sex_name="Both"

gbd_subset <- gbdfin %>%
    filter(year == 2017,
          metric_name == "Number",
          measure_name %in% c("Incidence", "Deaths", "Prevalence"),
          !(age_id %in% c(2,3,4,5,148,45,33,28,199,44, 21, 22, 32)),
          #age_name %in% age_name,
          sex_name=="Both",
          cause_name %in% ltr_df$cause_name) %>%
    select(c("measure_name", "sex_name", "age_name", "cause_name", "val")) %>%
    mutate(measure_name = case_when(measure_name == "Incidence" ~ "outcome_incidence",
                                    measure_name == "Deaths" ~ "outcome_deaths",
                                    measure_name == "Prevalence" ~ "prevalence")) %>%
    pivot_wider(names_from = measure_name, values_from = val)

gbd_subset$outcome_deaths[is.na(gbd_subset$outcome_deaths)] <- 0
gbd_subset$prevalence[is.na(gbd_subset$prevalence)] <- 0
gbd_subset$outcome_incidence[is.na(gbd_subset$outcome_incidence)] <- 0

gbd_subset <- left_join(gbd_subset, pop[,c("age_name", "sex_name", "population")], by=c("age_name", "sex_name"))
gbd_subset$prevalence <- gbd_subset$prevalence/gbd_subset$population
gbd_subset$outcome_incidence <- gbd_subset$outcome_incidence/gbd_subset$population
gbd_subset$outcome_mortality <- gbd_subset$outcome_deaths/gbd_subset$population
gbd_subset <- left_join(gbd_subset, acm[,c("age_name", "sex_name", "acm")], by=c("age_name", "sex_name"))
gbd_subset <- left_join(gbd_subset, lt[,c("age_name", "sex_name", "prob_death")], by=c("age_name", "sex_name"))

# function for taking ltr_df through future_pmap to parallelize calculations
get_age_df_list <- function(rsid,
                            endpoint_name,
                            cohort,
                            method,
                            cause_name,
                           fg_df,
                           ltr0,
                           ltr1,
                           ltr2,
                           population_dalys) {
    if (method %in% c("coxph", "coxph_single")) {
        endpoint <- endpoint_name
        gbd_endpoint <- cause_name
        endpoint_age <- paste(endpoint, "_AGE", sep="")
        #rsid <- ltr_df$rsid[i]
        temp_selected_colnames <- c(selected_colnames, endpoint, endpoint_age, rsid)
        temp_fg_df <- fg_df[,temp_selected_colnames]
        age_df <- get_age_df_se(endpoint=endpoint,
                             rsid=rsid,
                             method=method,
                             cohort=cohort,
                             fg_df=temp_fg_df)

        age_df$rsid <- rsid
        age_df$endpoint_name <- endpoint
        age_df$cause_name <- gbd_endpoint
        age_df$cohort <- cohort
        age_df$method <- method
        return(age_df)
    }
}


#save.image(file = "gbd_workspace.RData", compress=FALSE)

load("./gbd_workspace.RData")

#200 GB RAM 32 vCPU settings: approx 16h, 31 workers, mult = 7, did not return results
library(furrr)
#plan(multiprocess, workers=31)
mult=9
options(future.globals.maxSize= 891289600*mult) #850mb * x
# maximum number of workers by RAM
400/(mult*0.850)

# Run parallelized version
selected_colnames <- c("FINNGENID", "YEAR_OF_BIRTH", "FU_START_AGE",
                       "FU_END_AGE", "AGE_AT_DEATH_OR_NOW", "DEATH",
                       "DEATH_AGE", "FEMALE", pcs)
library(tictoc)
tic()
plan(multiprocess, workers=50)
results <- future_pmap(.l=ltr_df, .f=get_age_df_list, fg_df=fg_df)
toc()


# create partition indices
vec1 <- 1+(0:21)*10000
vec1
vec2 <- (1:22)*10000
vec2[22] <- nrow(ltr_df)
vec2

paste(vec1, ":", vec2, sep="")

selected_colnames <- c("FINNGENID", "YEAR_OF_BIRTH", "FU_START_AGE",
                       "FU_END_AGE", "AGE_AT_DEATH_OR_NOW", "DEATH",
                       "DEATH_AGE", "FEMALE", pcs)
results_list <- list()
i=1
partition <- ltr_df[vec1[i]:vec2[i],]
tic()
plan(multiprocess, workers=50)
results_list[[1]] <- future_pmap(.l=partition, .f=get_age_df_list, fg_df=fg_df)
toc()

56742.575/(60*60)

#15.76h runtime, 50 workers

selected_colnames <- c("FINNGENID", "YEAR_OF_BIRTH", "FU_START_AGE",
                       "FU_END_AGE", "AGE_AT_DEATH_OR_NOW", "DEATH",
                       "DEATH_AGE", "FEMALE", pcs)

tic()
results_list <- list()
for (i in 1:length(vec1)) {
    plan(multiprocess, workers=50)
    partition <- ltr_df[vec1[i]:vec2[i],]
    print(paste(vec1[i], ":", vec2[i], sep=""))
    results_list[[i]] <- future_pmap(.l=partition, .f=get_age_df_list, fg_df=fg_df)
}
toc()

#saveRDS(results_list, "/home/jsjukara/results.rds")

results_list <- readRDS("/home/jsjukara/results.rds")

# bind results from results_list to results
results <- bind_rows(results_list[[1]])

for (i in 2:length(results_list)) {
    results <- rbind(results, bind_rows(results_list[[i]]))
}

# Loop for SE version of age_df
selected_colnames <- c("FINNGENID", "YEAR_OF_BIRTH", "FU_START_AGE",
                       "FU_END_AGE", "AGE_AT_DEATH_OR_NOW", "DEATH",
                       "DEATH_AGE", "FEMALE", "PC1", "PC2", "PC3", "PC4")
library(tictoc)
tic()
for (i in 75750:75968) {
    if (ltr_df$method[i] %in% c("coxph", "coxph_single")) {
        endpoint <- ltr_df$endpoint_name[i]
        gbd_endpoint <- ltr_df$cause_name[i]
        endpoint_age <- paste(endpoint, "_AGE", sep="")
        rsid <- ltr_df$rsid[i]
        temp_selected_colnames <- c(selected_colnames, endpoint, endpoint_age, rsid)
        temp_fg_df <- fg_df[,temp_selected_colnames]
        age_df <- get_age_df_se(endpoint=endpoint,
                             rsid=rsid,
                             method=ltr_df$method[i],
                             cohort=ltr_df$cohort[i],
                             fg_df=temp_fg_df)

        age_df$rsid <- rsid
        age_df$endpoint_name <- endpoint
        age_df$cause_name <- gbd_endpoint
        age_df$cohort <- ltr_df$cohort[i]
        age_df$method <- ltr_df$method[i]
        age_df_list[[paste(ltr_df$cohort[i], ltr_df$method[i], rsid, endpoint)]] <- age_df
    }

}
toc()


workspace.size <- function() {
  ws <- sum(sapply(ls(envir=globalenv()), function(x)object.size(get(x))))
  class(ws) <- "object_size"
  ws
}

print(workspace.size(),units="Mb")
                   

adjust_hr_bayes <- function(age_df_all, mult) {
    # Bayesian adjustment of age group HR's
    bayes <- age_df_all %>% filter(method=="coxph")

    #get normally distributed prior distribution from coxph_single model
    temp <- age_df_all %>% filter(method=="coxph_single")

    #assign priors
    m.pr <- temp$loghr
    s.pr <- temp$loghr_se * mult
    #assign observed values
    m.obs <- bayes$loghr
    s.obs <- bayes$loghr_se
    
    #calculate posteriors
    s.po <- 1/sqrt(1/s.pr^2 + 1/s.obs^2)
    m.po <- s.po^2 * (m.pr/s.pr^2 + m.obs/s.obs^2)

    #if posteriors are NA (from missing observed values), use priors
    m.po <- ifelse(is.na(m.po), m.pr, m.po)
    s.po <- ifelse(is.na(s.po), s.pr, s.po)
    
    bayes$loghr <- m.po
    bayes$loghr_se <- s.po

    #calculate HR and 95% CI
    bayes$hr01 <- exp(m.po)
    bayes$hr01loci <- exp(log(bayes$hr01) - 1.96*s.po)
    bayes$hr01hici <- exp(log(bayes$hr01) + 1.96*s.po)
    
    #name the method
    bayes$method <- paste("coxph_bayes", mult, sep="")
    
    return(bayes)
}

# fix gbdfin cause_name formatting
gbdfin$cause_name <- gsub("\\-", "_", gbdfin$cause_name)
gbdfin$cause_name <- gsub("/", "_", gbdfin$cause_name)

m_vec <- c(6)
age_df_all <- results
age_df_all$mult <- NA
age_df_all$hr01loci <- exp(age_df_all$loghr - 1.96 * age_df_all$loghr_se)
age_df_all$hr01hici <- exp(age_df_all$loghr + 1.96 * age_df_all$loghr_se)

# loop over different values for mult to get different SD for prior distribution
# extend analyses
for (mult in m_vec) {
    temp <- adjust_hr_bayes(age_df_all, mult)
    temp$mult <- mult
    age_df_all <- rbind(age_df_all, temp)
}
age_df_all$hr02 <- exp(age_df_all$loghr*2)

age_df_all <- left_join(age_df_all, genotype_freqs[,c("rsid", "maf", "mafover1")])

bayes_include <- c("coxph_bayes6")
age_df_all <- age_df_all %>%
    filter(method %in% c("coxph", "coxph_single", bayes_include))

#temp fix to get ltr for na
temp_indices <- which(is.na(ltr_df$ltr0))
tic()
for (i in temp_indices) {
    endpoint <- ltr_df$endpoint_name[i]
    gbd_endpoint <- ltr_df$cause_name[i]
    endpoint_age <- paste(endpoint, "_AGE", sep="")
    rsid <- ltr_df$rsid[i]
    
    if (ltr_df$method[i] %in% c("coxph", "coxph_single")) {
        age_df <- age_df_all %>%
            filter(cohort==ltr_df$cohort[i],
                  method==ltr_df$method[i],
                  rsid==ltr_df$rsid[i],
                  endpoint_name==ltr_df$endpoint_name[i])
        ltr <- get_ltr_incidence(outcome_name = gbd_endpoint,
                                 age_df = age_df,
                                 gbdfin=gbd_subset,
                                 rsid=rsid,
                                 genotype_freqs=genotype_freqs)

        ltr_df[i, 6:8] <- ltr
    } else {
        age_df <- age_df_all %>%
            filter(cohort==ltr_df$cohort[i],
                  method==ltr_df$method[i],
                  rsid==ltr_df$rsid[i],
                  endpoint_name==ltr_df$endpoint_name[i])
        
        ltr <- get_ltr_incidence(outcome_name = gbd_endpoint,
                                 age_df = age_df,
                                 gbdfin=gbd_subset,
                                 rsid=rsid,
                                 genotype_freqs=genotype_freqs)

        ltr_df[i, 6:8] <- ltr
        
    }
}
toc()


tic()
for (i in 1:nrow(ltr_df)) {
    endpoint <- ltr_df$endpoint_name[i]
    gbd_endpoint <- ltr_df$cause_name[i]
    endpoint_age <- paste(endpoint, "_AGE", sep="")
    rsid <- ltr_df$rsid[i]
    
    if (ltr_df$method[i] %in% c("coxph", "coxph_single")) {
        age_df <- age_df_all %>%
            filter(cohort==ltr_df$cohort[i],
                  method==ltr_df$method[i],
                  rsid==ltr_df$rsid[i],
                  endpoint_name==ltr_df$endpoint_name[i])
        ltr <- get_ltr_incidence(outcome_name = gbd_endpoint,
                                 age_df = age_df,
                                 gbdfin=gbd_subset,
                                 rsid=rsid,
                                 genotype_freqs=genotype_freqs)

        ltr_df[i, 6:8] <- ltr
    } else {
        age_df <- age_df_all %>%
            filter(cohort==ltr_df$cohort[i],
                  method==ltr_df$method[i],
                  rsid==ltr_df$rsid[i],
                  endpoint_name==ltr_df$endpoint_name[i])
        
        ltr <- get_ltr_incidence(outcome_name = gbd_endpoint,
                                 age_df = age_df,
                                 gbdfin=gbd_subset,
                                 rsid=rsid,
                                 genotype_freqs=genotype_freqs)

        ltr_df[i, 6:8] <- ltr
        
    }
}
toc()


table(is.na(ltr_df$ltr1), ltr_df$method)

workspace.size <- function() {
  ws <- sum(sapply(ls(envir=globalenv()), function(x)object.size(get(x))))
  class(ws) <- "object_size"
  ws
}

print(workspace.size(),units="Mb")
                   
#save.image(file = "gbd_workspace1.RData", compress=FALSE)
#load("./gbd_workspace1.RData")

ltr_df <- left_join(ltr_df, genotype_freqs, by="rsid")

n_population <- 5511371
ltr_df$n0 <- n_population*ltr_df$p0
ltr_df$n1 <- n_population*ltr_df$p1
ltr_df$n2 <- n_population*ltr_df$p2
ltr_df$denom <- (ltr_df$ltr0 * ltr_df$p0 + ltr_df$ltr1 * ltr_df$p1 + ltr_df$ltr2 * ltr_df$p2)
ltr_df$pop_dalys0 <- ltr_df$ltr0*ltr_df$p0 / ltr_df$denom * ltr_df$population_dalys
ltr_df$pop_dalys1 <- ltr_df$ltr1*ltr_df$p1 / ltr_df$denom * ltr_df$population_dalys
ltr_df$pop_dalys2 <- ltr_df$ltr2*ltr_df$p2 / ltr_df$denom * ltr_df$population_dalys
ltr_df$person_dalys0 <- ltr_df$pop_dalys0 / ltr_df$n0
ltr_df$person_dalys1 <- ltr_df$pop_dalys1 / ltr_df$n1
ltr_df$person_dalys2 <- ltr_df$pop_dalys2 / ltr_df$n2

ltr_df$attr_person_dalys1 <- ltr_df$person_dalys1 - ltr_df$person_dalys0
ltr_df$attr_person_dalys2 <- ltr_df$person_dalys2 - ltr_df$person_dalys0

ltr_df$attr_population_dalys1 <- ltr_df$attr_person_dalys1 * ltr_df$n1
ltr_df$attr_population_dalys2 <- ltr_df$attr_person_dalys2 * ltr_df$n2

ltr_df$person_attr1_lifetime <- ltr_df$attr_person_dalys1*82
ltr_df$person_attr2_lifetime <- ltr_df$attr_person_dalys2*82
ltr_df$attr_population_dalys12 <- ltr_df$attr_population_dalys1 + ltr_df$attr_population_dalys2

ltr_df_b <- ltr_df %>% filter(method=="coxph_bayes6")


# function for getting n variants with top absolute lifetime DALYS
get_top_lifetime_variants <- function(method, ltr_df, n) {
    temp <- ltr_df %>%
        filter(method %in% c("coxph_bayes6")) %>%
        group_by(rsid) %>%
        summarize(dalys_lifetime = sum(person_attr1_lifetime, na.rm=T)) %>%
        #filter(sum_dalys_lifetime >0.5) %>%
        arrange(desc(abs(dalys_lifetime)))
    print(temp[1:n,])
    top_n <- temp$rsid[1:n]
}

# function for getting n variants with top absolute population DALYS
get_top_population_variants <- function(method, ltr_df, n) {
    temp <- ltr_df %>%
        filter(method %in% c("coxph_bayes6")) %>%
        group_by(rsid) %>%
        summarize(dalys_population = sum(attr_population_dalys12, na.rm=T)) %>%
        #filter(sum_dalys_lifetime >0.5) %>%
        arrange(desc(abs(dalys_population)))
    print(temp[1:n,])
    top_n <- temp$rsid[1:n]
}

top10_lifetime <- get_top_lifetime_variants(method="coxph_bayes6", ltr_df=ltr_df_b, n=10)
top10_population <- get_top_population_variants(method="coxph_bayes6", ltr_df=ltr_df_b, n=10)

get_info_rsid <- function(rs, ltr_df, met, population=TRUE) {
    data <- ltr_df[ltr_df$rsid==rs & ltr_df$method==met,]
    print(rs)
    print(paste("Sum of attributable lifetime DALYs for heterozygotes: ",
                round(sum(data$person_attr1_lifetime, na.rm=T),3), sep=""))
    print(paste("Sum of yearly population DALYs attributable to variant: ",
                round(sum(data$attr_population_dalys12, na.rm=T),0), sep=""))
    print(paste("MAF: ", round(data$maf[1]*100, 2), "%", sep=""))
    if (population==FALSE) {
        temp <- data %>%
            mutate(risk_difference1=ltr1-ltr0) %>%
            arrange(desc(abs(person_attr1_lifetime)))
        endpoints <- temp$endpoint_name[1:9]
        #temp1 <- temp %>%
        #    filter(endpoint_name %in% endpoints)
    
    } else {
        temp <- data %>%
            mutate(risk_difference1=ltr1-ltr0) %>%
            arrange(desc(abs(attr_population_dalys12)))
        endpoints <- temp$endpoint_name[1:9]
    }
    
    temp1 <- temp[1:9,c("endpoint_name", "person_attr1_lifetime", "attr_population_dalys12")]
    temp1$person_attr1_lifetime <- round(temp1$person_attr1_lifetime,3)
    temp1$attr_population_dalys12 <- round(temp1$attr_population_dalys12,0)
    print(temp1)
    
    plot(temp %>%
        filter(endpoint_name %in% endpoints) %>%
        gather(key=var, value=val, person_attr1_lifetime, risk_difference1, attr_population_dalys12) %>%
        ggplot(aes(x=endpoint_name, y=val)) +
        geom_point() +
        facet_wrap(~var, scales="free") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45)) +
        geom_hline(yintercept=0) +
        labs(title="Top 8 endpoints with strongest association to attributable lifetime DALYs",
            x="FinnGen Endpoint Name",
            y="Lifetime DALY's attributable to carrying 1 vs 0 variants"))
    data8 <- data %>%
        filter(endpoint_name %in% endpoints) %>%
        mutate(pop_dalys1m = pop_dalys1 - attr_population_dalys1,
              pop_dalys2m = pop_dalys2 - attr_population_dalys2) %>%    
        select(endpoint_name, pop_dalys0, pop_dalys1m, pop_dalys2m, attr_population_dalys1,
               attr_population_dalys2, population_dalys, attr_population_dalys2) %>%
        pivot_longer(-endpoint_name, names_to="var", values_to="val") %>%
        mutate(genotype = case_when(var %in% c("pop_dalys1m", "attr_population_dalys1") ~ "snp=1",
                                   var %in% c("pop_dalys2m", "attr_population_dalys2") ~ "snp=2",
                                   TRUE ~ "snp=0"),
              attr = case_when(var %in% c("attr_population_dalys1", "attr_population_dalys2") ~ "Attributable",
                              TRUE ~ "Not attributable"))
    
    plot(data8 %>%
        filter(var %in% c("pop_dalys0", "pop_dalys1m", "pop_dalys2m",
                          "attr_population_dalys1", "attr_population_dalys2")) %>%
        ggplot(aes(fill=attr, y=val, x=genotype)) +
            geom_bar(position="stack", stat="identity") +
            facet_wrap(~endpoint_name, scales="free"))
        
    
}
options(repr.plot.width=14, repr.plot.height=6)

get_info_rsid(rs=top10_lifetime[1], met="coxph_bayes6", ltr_df=ltr_df, population=FALSE)

get_info_rsid(rs=top10_lifetime[2], met="coxph_bayes6", ltr_df=ltr_df, population=FALSE)

get_info_rsid(rs=top10_lifetime[3], met="coxph_bayes6", ltr_df=ltr_df, population=FALSE)

get_info_rsid(rs=top10_lifetime[4], met="coxph_bayes6", ltr_df=ltr_df, population=FALSE)

genotype_freqs %>% filter(rsid == "chr1_145804971_G_C")

for (rsid in top10_population) {
    get_info_rsid(rs=rsid, met="coxph_bayes6", ltr_df=ltr_df, population=TRUE)
    Sys.sleep(0.2)
}


round(cor(fg_df[,top10_lifetime]),2)

round(cor(fg_df[,top10_population]),2)

ltr_df_nona <- ltr_df[!is.na(ltr_df$ltr0),]

ltr_df_aggregated <- ltr_df_nona %>%
    filter(method=="coxph_bayes6") %>%
    group_by(rsid) %>%
    summarize(sum_person_attr1_lifetime=sum(person_attr1_lifetime),
             sum_attr_population_dalys12=sum(attr_population_dalys12)) %>%
    arrange(desc(abs(sum_attr_population_dalys12)))

ltr_df_aggregated_top100i <- ltr_df_aggregated %>%
    arrange(desc(abs(sum_attr_population_dalys12)))

ltr_df_aggregated_top100i <- ltr_df_aggregated_top100i[1:100,]

head(ltr_df_aggregated_top100i)

ltr_df_aggregated_top100p <- ltr_df_aggregated %>%
    arrange(desc(abs(sum_person_attr1_lifetime)))

ltr_df_aggregated_top100p <- ltr_df_aggregated_top100p[1:100,]

head(ltr_df_aggregated_top100p)


#ltr_df_top100i <- ltr_df %>%
#    filter(method=="coxph_bayes6") %>%
#    arrange(desc(abs(person_attr1_lifetime)))


summary(ltr_df_aggregated)

library(lattice)
densityplot(ltr_df_aggregated$sum_person_attr1_lifetime)

qqnorm(ltr_df_aggregated$sum_person_attr1_lifetime)
qqline(ltr_df_aggregated$sum_person_attr1_lifetime)

options(repr.plot.width=10, repr.plot.height=8)
ltr_df_aggregated_top100i %>%
    ggplot(aes(x=reorder(rsid, sum_person_attr1_lifetime), y=sum_person_attr1_lifetime)) +
    geom_point() +
    geom_hline(yintercept=0) +
    labs(y="Attributable individual lifetime DALYs") +
    theme_bw() +
    coord_flip() +
    labs(x="", title="")

options(repr.plot.width=10, repr.plot.height=8)
ltr_df_aggregated_top100p %>%
    ggplot(aes(x=reorder(rsid, sum_attr_population_dalys12), y=sum_attr_population_dalys12)) +
    geom_point() +
    geom_hline(yintercept=0) +
    labs(y="Attributable individual lifetime DALYs") +
    theme_bw() +
    coord_flip() +
    labs(x="", title="")

ltr_df %>%
    filter(method=="coxph_bayes6")) %>%
    ggplot(aes(x=reorder(measure, estimate), y=estimate)) +
    geom_errorbar(aes(x=reorder(measure, estimate), ymin=loci, ymax=hici)) +
    geom_point() +
    geom_hline(yintercept=0) +
    labs(y="Attributable individual lifetime DALYs") +
    theme_bw() +
    coord_flip() +
    labs(x="", title=="")
    

options(repr.plot.width=5, repr.plot.height=6)
#temp <- get_age_df()
single <- age_df_all %>%
    filter(endpoint_name=="KRA_PSY_DEMENTIA" & rsid == "chr19_44908684_T_C" & method == "coxph_single") %>%
    pull(hr01)
age_df_all %>%
    filter(endpoint_name=="KRA_PSY_DEMENTIA",
           rsid == "chr19_44908684_T_C",
           method != "coxph_single",
           method %in% c("coxph", "coxph_bayes6")) %>%
    ggplot(aes(x=age_low, y=hr01, colour=method)) +
    geom_errorbar(aes(ymin=hr01loci, ymax=hr01hici), alpha=1, position=position_dodge(width=3)) +
    geom_point(alpha=1, position=position_dodge(width=3), size=1) +
    scale_y_continuous(limits=c(0,8), breaks=seq(0,8,1)) +
    geom_hline(yintercept=c(1), linetype=1) +
    geom_hline(yintercept=c(single), linetype=2) +
    labs(y="Hazard ratio, 0 vs 1 copies", x="Age group", title="Dementia and rs429358 APOE variant") +
    theme_bw() +
    scale_colour_discrete(name = "Method", labels=c("Stratified Cox model", "Posterior estimates, m=6")) +
    theme(legend.position = "top")

# look how much HR changes before 40 vs after 40
# pick m that in greater than 40 gives similar HR than original, before 40 maximizes difference

#plot differences between non prior vs prior HR in bins with more than 50-100 outcomes, and in bins with smaller