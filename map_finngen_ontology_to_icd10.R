library(readxl)
library(tidyverse)
library(data.table)
# read FinnGen ontology excel
ontology <- read_excel("/Users/jsjukara/Dropbox/FIMM/GWAS\ Impact/ukbb/FINNGEN_ENDPOINTS_DF5_V2_2020-02-11_public.xlsx")

# read the table of endpoints that contain the sets of FinnGen endpoints that define the new endpoint
# e.g. "C3_COLON|C3_RECTUM"
map <- fread("/Users/jsjukara/Dropbox/FIMM/GWAS\ Impact/ukbb/gb_to_finngen_endpoints.csv")

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
      # get the unique set of regex's
      icd10_regex <- unique(icd10_regex[!is.na(icd10_regex)])
      # paste them with | together
      icd10_regex <- paste(icd10_regex, collapse="\\|")
      return(icd10_regex)
    } else { # else recurse into the composite endpoints
      composite_endpoints <- row$INCLUDE[1]
      return(get_icd10_regex(endpoints = composite_endpoints, ontology=ontology))
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


for (endpoint in endpoints) {
  icd10_regex_new <- get_icd10_regex(endpoints = endpoint, ontology=ontology)
  #icd10_regex_new <- ontology$row$HD_ICD_10[ontology$NAME == endpoint]
  icd10_regex <- paste(icd10_regex, icd10_regex_new, sep="|")
  return(icd10_regex)
}

get_icd10_regex("D3_ANAEMIA_IRONDEF", ontology=ontology)
get_icd10_regex("C3_CANCER|C3_COLORECTAL", ontology=ontology)
get_icd10_regex("C3_COLORECTAL", ontology=ontology)
get_icd10_regex("C3_COLON", ontology=ontology)



for (i in 1:nrow(map)) {
  map$icd10_regex[i] <- get_icd10_regex(map$fg_endpoints[i], ontology=ontology)
  
}

grepl("H40[8-9]||H42", "H72")
