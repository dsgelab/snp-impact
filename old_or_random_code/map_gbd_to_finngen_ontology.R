library(tidyverse)

# Load GBD data and causse maps to ICD-10 ---------------------------------
gbdfin <- read_csv("/Users/jsjukara/Dropbox/FIMM/GWAS\ Impact/Auxiliary data/IHME-GBD_2017_DATA_FINLAND/IHME-GBD_2017_DATA_FINLAND.csv")

# Causes of death
map_cod <- read_csv("/Users/jsjukara/Dropbox/FIMM/GWAS\ Impact/Auxiliary data/IHME_GBD_2017_ICD_CAUSE_MAPS/IHME_GBD_2017_ICD_CAUSE_MAP_CAUSES_OF_DEATH_Y2018M11D08.csv")
map_cod <- as.data.frame(map_cod)
map_cod <- map_cod[-c(311:314),]
# fix erroneous diagnosis code C47-C4A -> C47-C47.9
map_cod[149,"ICD10"] <- "C17-C17.9, C30-C31.9, C37-C38.8, C40-C41.9, C47-C47.9, C51-C52.9, C57-C57.8, C58-C58.0, C60-C60.9, C63-C63.8, C66-C66.9, C68.0-C68.8, C69-C69.9, C74-C75.8, D07.4, D09.2, D13.2-D13.3, D14.0, D15-D16.9, D28.0-D28.1, D28.7, D29.0, D30.2, D30.4-D30.8, D31-D31.9, D35-D35.2, D35.5-D36, D36.1-D36.7, D37.2, D38.2-D38.5, D39.2, D39.8, D41.2-D41.3, D44.1-D44.8, D48.0-D48.4"

# Nonfatal causes
map_nonfatal <- read_csv("/Users/jsjukara/Dropbox/FIMM/GWAS\ Impact/Auxiliary data/IHME_GBD_2017_ICD_CAUSE_MAPS/IHME_GBD_2017_ICD_CAUSE_MAP_NONFATAL_CAUSES_Y2018M11D08.csv")
map_nonfatal <- map_nonfatal[-302,]
map_nonfatal <- as.data.frame(map_nonfatal)
# fix erroneous "E90-E998" -> "E90-E99.8"
map_nonfatal[258,"ICD10"] <- "D66-D69.49, D69.6-D70.0, D70.2-D77, D80-D84.9, D86.8, D86.82-D86.84, D86.86-D86.87, D89-D89.9, E03-E03.1, E03.3-E06.3, E06.5-E07.9, E15-E16, E16.1-E16.9, E20-E23.0, E23.2-E24.1, E24.3-E27.2, E27.4-E28.1, E28.3-E32.9, E34-E35.8, E65-E66.09, E66.2-E68, E70-E80.09, E80.3-E88.9, E90-E99.8, Z83.4-Z83.49"

# Load GBD cause names by level
level3 <- read_csv("/Users/jsjukara/Dropbox/FIMM/GWAS\ Impact/Auxiliary data/Level 3 causes GBD.csv")
level3 <- filter(level3, metric_name=="Number", measure_name=="Deaths")
level3_id <- level3$cause_id
level3_name <- level3$cause_name

level4 <- read_csv("/Users/jsjukara/Dropbox/FIMM/GWAS\ Impact/Auxiliary data/Level 4 causes GBD.csv")
level4 <- filter(level4, metric_name=="Number", measure_name=="Deaths")
level4_id <- level4$cause_id
level4_name <- level4$cause_name

# Load ICD codes ----------------------------------------------------------
icd10_dgs <- read_csv("/Users/jsjukara/Dropbox/FIMM/GWAS\ Impact/Auxiliary data/icd10 tautiluokitus thl diagnoses vector 08012020.csv")
icd10_dgs <- as.data.frame(icd10_dgs)[,1]
icd10_dgs[nchar(icd10_dgs)==6]


# Add leading zeroes to three character codes to duplicate them with .0 at the end
# for example, when at "B24", will add "B24.0" at end of vector
# this is needed for checking whether a diagnosis from get_vector_of_codes()
# is actually an ICD-10 code, since get_vector_of_codes() might output "B24.0" whereas ICD-10 only has "B24"
for (i in 1:length(icd10_dgs)) {
  if (nchar(icd10_dgs[i])==3) {
    icd10_dgs <- c(icd10_dgs, paste(icd10_dgs[i], ".0", sep=""))
  }
}


# Define functions to extract individual diagnoses ------------------------

# Function that returns a vector of diagnoses between given "start" and "end" diagnoses
get_intervening_diagnoses <- function(start, end) {
  letter <- substr(start, 1, 1)
  endletter <- substr(end, 1, 1)
  start <- substr(start,2,20)
  end <- substr(end,2,20)
  if (grepl("\\.", start)) {
    temp <- unlist(strsplit(start, "\\."))
    start_major <- temp[1]
    start_minor <- temp[2]
  } else {
    start_major <- start
    start_minor <- "0"
  }
  if (grepl("\\.", end)) {
    temp <- unlist(strsplit(end, "\\."))
    end_major <- temp[1]
    end_minor <- temp[2]
  } else {
    end_major <- end
    end_minor <- "9"
  }
  current <- paste(start_major, start_minor, sep=".")
  current_major <- as.integer(start_major)
  current_minor <- as.integer(start_minor)
  current_letter <- letter
  letter_vec <- NULL
  diagnoses_vec <- NULL
  while (current_minor != as.integer(end_minor) | current_major != as.integer(end_major)) {
    # if has to change letter to new one:
    if (current_minor == 100 & current_major == 99) {
    current_minor <- 0
    current_major <- 0
    current_letter <- endletter
    letter_vec <- c(letter_vec, current_letter)
    } else if (current_minor < 100) { #grow minor until 99
      diagnoses_vec <- c(diagnoses_vec, paste(current_major, current_minor, sep="."))
      current_minor <- current_minor + 1
      letter_vec <- c(letter_vec, current_letter)
    } else { # increment major, set minor to 0
      current_major <- current_major + 1
      current_minor <- 0
    }
  }
  diagnoses_vec <- c(diagnoses_vec, paste(current_major, current_minor, sep="."))
  diagnoses_vec <- paste(letter_vec, diagnoses_vec, sep="")
  return(diagnoses_vec)
}
# Examples:
#get_intervening_diagnoses("B20","B23.8")
#get_intervening_diagnoses("B20.1","B23")
#get_intervening_diagnoses("B20.9","B23.2")

# Function that takes a string containing GBD ICD-10 codes (e.g. "B20-B23.8, B24-B24.0,...")
# and outputs a vector of individual diagnoses spanning the intervals
get_vector_of_codes <- function(x) {
  # initialize vector to contain separate diagnoses
  diagnoses_vec <- NULL
  x <- as.character(x)
  # remove spaces
  x <- gsub(" ", "", x)
  # split into vector at commas
  xvec  <- unlist(strsplit(x, ","))
  # loop over vector
  for (i in 1:length(xvec)) {
    # get ith diagnosis or range of diagnoses
    dg <- xvec[i]
    # if range, store a vector of diagnoses
    if (grepl("-", dg)) {
      temp <- unlist(strsplit(dg, "-"))
      start <- temp[1]
      start <- substr(start,0,6)
      end <- temp[2]
      end <- substr(end,0,6)
      new_dgs <- get_intervening_diagnoses(start, end)
      diagnoses_vec <- append(diagnoses_vec, new_dgs) 
    } else { # else just store one
      diagnoses_vec <- append(diagnoses_vec, dg)
    }
  }
  return(diagnoses_vec)
}

vec_to_char <- function(x) {
  return(paste(x, collapse=","))
}


# Extract individual diagnoses --------------------------------------------

map_nonfatal$ICD10_parsed <- NA
for (i in 1:nrow(map_nonfatal)) {
  vec <- get_vector_of_codes(map_nonfatal[i,2])
  vec <- vec[vec %in% icd10_dgs]
  map_nonfatal$ICD10_parsed[i] <- vec_to_char(vec)
}

map_nonfatal <- map_nonfatal[,c(1,2,6)]


map_cod$ICD10_parsed <- NA
for (i in 1:nrow(map_cod)) {
  vec <- get_vector_of_codes(map_cod[i,2])
  vec <- vec[vec %in% icd10_dgs]
  map_cod$ICD10_parsed[i] <- vec_to_char(vec)
}

map_cod <- map_cod[,c(1,2,4)]

write.csv(map_cod, "map_cod.csv", row.names = FALSE)
write.csv(map_nonfatal, "map_nonfatal.csv", row.names = FALSE)


library(tidyverse)
map_cod <- read_csv("map_cod.csv")
map_nonfatal <- read_csv("map_nonfatal.csv")

both <- full_join(map_cod, map_nonfatal, by="Cause")

both$ICD10_parsed_union<- NA
for (i in 1:nrow(both)) {
  dg1 <- unlist(strsplit(both$ICD10_parsed.x[i], ","))
  dg2 <- unlist(strsplit(both$ICD10_parsed.y[i], ","))
  if (is.na(dg1[1]) & is.na(dg2[1])) { # both diagnoses missing
    both$ICD10_parsed_union[i] <- NA
  } else if (is.na(dg1[1]) & !is.na(dg2[1])) { # diagnosis 1 missing
    both$ICD10_parsed_union[i] <- paste(dg2, collapse=",")
  } else if (!is.na(dg1[1]) & is.na(dg2[1])) { # diagnosis 2 missing
    both$ICD10_parsed_union[i] <- paste(dg1, collapse=",")
  } else { # both present
    both$ICD10_parsed_union[i] <- paste(union(dg1, dg2), collapse=",")
  }
}

names(both) <- c("Cause", "ICD10_cod", "ICD10_cod_parsed", "ICD10_nonfatal", "ICD10_nonfatal_parsed", "ICD10_union_parsed")

both <- as.data.frame(both)

both_withdot <- both

for (i in 2:ncol(both)) {
  both[,i] <- gsub("\\.", "", both[,i])
}

  

write_tsv(both, "gbd_to_icd10_ontology 130220.tsv")




