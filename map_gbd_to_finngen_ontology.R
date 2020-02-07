library(tidyverse)

gbdfin <- read_csv("/Users/sakarijukarainen/Library/Mobile Documents/com~apple~CloudDocs/Documents/FIMM/GWAS Impact/Auxiliary data/IHME-GBD_2017_DATA_FINLAND/IHME-GBD_2017_DATA_FINLAND.csv")

gbdfin %>%
  filter(measure_name=="Prevalence")
  

map_cod <- read_csv("/Users/sakarijukarainen/Library/Mobile Documents/com~apple~CloudDocs/Documents/FIMM/GWAS Impact/Auxiliary data/IHME_GBD_2017_ICD_CAUSE_MAPS/IHME_GBD_2017_ICD_CAUSE_MAP_CAUSES_OF_DEATH_Y2018M11D08.csv")

map_nonfatal <- read_csv("/Users/sakarijukarainen/Library/Mobile Documents/com~apple~CloudDocs/Documents/FIMM/GWAS Impact/Auxiliary data/IHME_GBD_2017_ICD_CAUSE_MAPS/IHME_GBD_2017_ICD_CAUSE_MAP_NONFATAL_CAUSES_Y2018M11D08.csv")


test <- as.character(map_nonfatal[1,2])

test
test <- gsub(" ", "", test)
test <- unlist(strsplit(test, ","))

test[1]
temp <- unlist(strsplit(test[1], "-"))
letter <- substr(temp[1], 1, 1)
nums <- substr(temp,2,20)
nums[1] <- as.numeric(nums[])

nums[1] <- as.integer(unlist(strsplit(nums[1], "\\.")))
dgs <- c()
for (i in nums[2] - nums[1])