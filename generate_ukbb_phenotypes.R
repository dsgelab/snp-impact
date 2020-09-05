library(readxl)
library(tidyverse)
library(data.table)
library(lubridate)
library(tictoc)
setwd("/home/jsjukara/ukbb/")

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
map <- fread("gb_to_finngen_endpoints.csv")
map <- map[map$fg_endpoints!="",]
names(map) <- tolower(names(map))

# initialize column that will contain the regex for the custom endpoint ICD10 codes
map$icd10_regex<- NA

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

ukbb_df <- ukbb_df %>% filter(eid %in% ids_ukbb)

# join genotypes to ukbb_df
ukbb_df <- left_join(ukbb_df, geno_ukbb, by="eid")

# remove ids without genotypes
ids_ukbb <- ids_ukbb[ids_ukbb %in% geno_ukbb$eid]
ukbb_df <- ukbb_df %>% filter(eid %in% ids_ukbb)


saveRDS(ukbb_df, "/home/jsjukara/ukbb/first_endpoints_ukbb.rds", compress=FALSE)
saveRDS(long_ukbb, "/home/jsjukara/ukbb/longitudinal_endpoints_ukbb.rds", compress=FALSE)
