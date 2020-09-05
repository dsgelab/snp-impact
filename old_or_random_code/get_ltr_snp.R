library(tidyverse)

outcome_name <- "Breast cancer"
sex <- "Female"

get_ltr_overall <- function(outcome_name, sex) {

  # Get incidences in cases and controls in FinnGen -------------------------
  ## TODO
  

  
  # Load and filter GBD life tables -----------------------------------------
  lt <- read.csv("/Users/sakarijukarainen/Library/Mobile Documents/com~apple~CloudDocs/Documents/FIMM/GWAS Impact/Auxiliary data/IHME_GBD_2017_ABRIDGED_LIFE_TABLES_2017_Y2018M11D08.CSV")
  lt <- lt %>%
    filter(location_name == "Finland" & measure_name == "Probability of death")
  
  names(lt)[c(6,12)] <- c("age_name", "prob_death")
  
  # Load population data ----------------------------------------------------
  pop <- read.csv("/Users/sakarijukarainen/Library/Mobile Documents/com~apple~CloudDocs/Documents/FIMM/GWAS Impact/Auxiliary data/IHME_GBD_2017_POP_2015_2017_Y2018M11D08.CSV")
  
  pop <- pop %>%
    filter(location_name == "Finland" &
             year_id == 2017)
  
  names(pop)[c(6,11)] <- c("age_name", "population")
  
  # Get incidences of outcome & all cause mortality from GBD -----------------------------------------------
  gbdfin <- read.csv("/Users/sakarijukarainen/Library/Mobile Documents/com~apple~CloudDocs/Documents/FIMM/GWAS Impact/Auxiliary data/IHME-GBD_2017_DATA_FINLAND/IHME-GBD_2017_DATA_FINLAND.csv")
  
  #df for all cause mortality
  acm <- gbdfin %>% filter(year == 2017 &
                             cause_name == "All causes" &
                             metric_name == "Number" &
                             #sex_name == "Both" &
                             measure_name %in% c("Deaths"))
  
  acm <- left_join(acm, pop[,c("age_name", "sex_name", "population")], by=c("age_name", "sex_name"))
  
  acm$acm <- acm$val / acm$population # calculate all-cause mortality rate
  
  # df for cause specific incidence, deaths and prevalence
  gbdfin <- gbdfin %>% filter(year == 2017 &
                                cause_name == outcome_name &
                                metric_name == "Number" &
                                measure_name %in% c("Incidence", "Deaths", "Prevalence") &
                                !(age_id %in% c(2,3,4,5,148,45,33,28,199,44, 21, 22, 32)))
  
  
  # Combine estimates to one table ------------------------------------------
  temp <- gbdfin %>% filter(measure_name == "Incidence")
  names(temp)[14] <- "incidence"
  
  # join no. outcome specific deaths
  temp <- left_join(temp, gbdfin[gbdfin$measure_name == "Deaths", c("sex_name", "age_name", "val")],
                    by=c("sex_name", "age_name"))
  
  names(temp)[17] <- "outcome_deaths"
  temp$outcome_deaths[is.na(temp$outcome_deaths)] <- 0
  
  # join prevalence of outcome (no. prevalent cases)
  temp <- left_join(temp, gbdfin[gbdfin$measure_name == "Prevalence", c("sex_name", "age_name", "val")],
                    by=c("sex_name", "age_name"))
  names(temp)[18] <- "prevalence"
  
  # join all cause mortality
  temp <- left_join(temp, acm[,c("age_name", "sex_name", "acm")], by=c("age_name", "sex_name"))
  
  # join probability of death in age group
  temp <- left_join(temp, lt[,c("age_name", "sex_name", "prob_death")], by=c("age_name", "sex_name"))
  
  # join population numbers
  temp <- left_join(temp, pop[,c("age_name", "sex_name", "population")], by=c("age_name", "sex_name"))
  
  # calculate prevalence proportion
  temp$prevalence <- temp$prevalence/temp$population
  
  
  # Estimate lifetime risk of outcome ---------------------------------------
  temp <- temp %>% filter(sex_name == sex)
  
  tab <- data.frame(age_name = temp$age_name, # age group name
                    l = NA, # population of hypothetical cohort
                    m = temp$acm, # death rate (per resident)
                    mo = temp$outcome_deaths/temp$population, #death rate from outcome (per resident)
                    rp = (temp$incidence / temp$population), # incidence rate of outcome in population
                    a = NA, # no. new cases within interval
                    d = NA) # no. non-outcome related deaths among outcome free individuals
  
  tab$gp = 1 - exp(-5*tab$rp) #prob of outcome in total population
  tab$g = tab$gp * (1 / (1 - temp$prevalence)) #prob of outcome in outcome-free population
  tab$r = (-1/5) * log(1 - tab$g) # incidence of outcome in outcome free population
  
  tab$l[1] <- 10*10^6 #set 10 million as population of hypothetical cohort
  
  for (i in 1:nrow(tab)) {
    if (i >= 2) {
      tab$l[i] <- tab$l[i-1] * exp(-5*(tab$m[i-1] + tab$r[i-1] - tab$mo[i-1]))
    }
  }
  tab$a <- tab$l * (1 - exp(-5*(tab$m + tab$r - tab$mo))) * (tab$r/(tab$m + tab$r - tab$mo))
  
  tab$d <- tab$l * (1 - exp(-5*(tab$m + tab$r - tab$mo))) * (tab$m/(tab$m + tab$r - tab$mo))
  
  # Calculate a and d for last open ended age interval
  ind <- nrow(tab)
  tab$a[ind] <- tab$l[ind] * (tab$r[ind]/(tab$r[ind] + tab$m[ind] - tab$mo[ind]))
  tab$d[ind] <- tab$l[ind] * (tab$m[ind]/(tab$r[ind] + tab$m[ind] - tab$mo[ind]))
  
  # Calculate lifetime risk
  ltr <- sum(tab$a)/tab$l[1]
  
  return(ltr)
}

get_ltr_overall("Tracheal, bronchus, and lung cancer", "Male")


