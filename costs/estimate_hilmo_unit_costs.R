setwd("/home/jsjukara/costs/")
library(data.table)
library(tidyverse)
library(lubridate)

# load HILMO unit costs table scraped from the THL report and do some wrangling
unit_costs_hilmo <- read_delim("unit_costs/unit_costs_hilmo_21_4_2020.csv", delim=";")[2:7]
names(unit_costs_hilmo) <- tolower(names(unit_costs_hilmo))
names(unit_costs_hilmo)[c(4,5)] <- c("hospital_type", "pala_class")
unit_costs_hilmo <- unit_costs_hilmo[!is.na(unit_costs_hilmo$ea),]
unit_costs_hilmo$unit_cost <- as.numeric(unit_costs_hilmo$unit_cost)
unit_costs_hilmo$ea <- as.character(unit_costs_hilmo$ea)

startyear <- 1972
endyear <- 2019

# load custom more detailed longitudinal phenotype file for FinnGen and do some wrangling
d <- fread("zcat finngen_R5_v3_custom_detailed_longitudinal_allcodes_REFINERY_ONLY.gz")
dnames <- names(d)
d$year <- decimal_date(as.Date(d$APPROX_EVENT_DAY))
d <- d[year >= startyear]
d <- d[year < endyear]


d$rowindex <- 1:nrow(d)
dhilmo <- d[SOURCE %in% c("INPAT", "OUTPAT")]
dprim <- d[SOURCE == "PRIM_OUT"]

# modify specialty
dhilmo$CODE6 <- as.integer(gsub("[^0-9.-]", "", dhilmo$CODE6))
dhilmo$year <- decimal_date(as.Date(dhilmo$APPROX_EVENT_DAY))
dhilmo <- dhilmo[year >= startyear]
dhilmo <- dhilmo[CATEGORY=="M"]
names(dhilmo) <- tolower(names(dhilmo))
names(dhilmo)[c(5:11)] <- c("pdgo", "pdge", "atc", "hospdays", "pala", "ea", "hospital_type")
# make specialties that are missing from unit costs table "-1"
missing_eas <- unique(dhilmo$ea)[!(unique(dhilmo$ea)%in%unique(unit_costs_hilmo$ea))]
dhilmo$ea[dhilmo$ea %in% missing_eas] <- "-1"
dhilmo$hospdays <- dhilmo$hospdays+1


# impute some missing values for the unit costs 
unit_costs_hilmo$thl_days[is.na(unit_costs_hilmo$thl_days)] <- 4.9
day_ward_costs <- unit_costs_hilmo %>% filter(pala_class=="ward")
day_ward_costs$thl_days[day_ward_costs$ea==15] <- mean(dhilmo$hospdays[dhilmo$ea==15 & dhilmo$source=="INPAT"])

day_ward_costs$hilmo_days_mean <- NA
day_ward_costs$hilmo_days_median <- NA
for (i in 1:nrow(day_ward_costs)) {
    ea <- day_ward_costs$ea[i]
    day_ward_costs$hilmo_days_mean[i] <- mean(dhilmo$hospdays[dhilmo$source=="INPAT" & dhilmo$ea==ea & dhilmo$year>=1998])
    day_ward_costs$hilmo_days_median[i] <- median(dhilmo$hospdays[dhilmo$source=="INPAT" & dhilmo$ea==ea & dhilmo$year>=1998])
}

# assuming 408€ per incremental ward episode day, estimate costs per episode based on a model of:
# costs = fixed costs + 408 * ward_days
day_ward_costs$cost_fixed <- day_ward_costs$unit_cost - 408 * day_ward_costs$thl_days
day_ward_costs$cost_fixed <- ifelse(day_ward_costs$cost_fixed<0, 0, day_ward_costs$cost_fixed)
day_ward_costs$cost_per_thlday <- day_ward_costs$unit_cost/day_ward_costs$thl_days
day_ward_costs$cost_per_day_wofixed <- (day_ward_costs$unit_cost-day_ward_costs$cost_fixed)/day_ward_costs$thl_days


pala_outpatient <- c("83", "92", "93", "94") #päiväsairaanhoito, ajanvaraus, ajanvaraus(uusintak.), konsultaatio
pala_ward <- c("1", "2", "5", "6", "7", "8")

dhilmo <- dhilmo %>%
    mutate(pala_class = case_when(pala == "91" & source == "OUTPAT" ~ "emergency",
                                 source == "OUTPAT" ~ "outpatient",
                                 source == "INPAT" ~ "ward"))

# join unit costs
dhilmo <- left_join(dhilmo, unit_costs_hilmo[unit_costs_hilmo$hospital_type=="all",c("ea", "pala_class", "unit_cost")], by=c("ea", "pala_class"))
# join fixed costs and thl_days
dhilmo <- left_join(dhilmo, day_ward_costs[,c("ea", "pala_class", "cost_fixed", "thl_days")],by=c("ea", "pala_class"))

# adjust ward episode unit cost based on the model of fixed costs + costs/day
dhilmo$unit_cost_w_fixed <- ifelse(dhilmo$pala_class == "ward",
                                   dhilmo$cost_fixed + (dhilmo$unit_cost - dhilmo$cost_fixed) * (dhilmo$hospdays / dhilmo$thl_days),
                                   dhilmo$unit_cost)

# for some specialties, instead use specific costs/day estimates from the report
dhilmo <- dhilmo %>%
    mutate(unit_cost_w_fixed = case_when(ea == "98" & pala_class == "ward" & hospdays<=90 ~hospdays*234,
                                ea == "98" & pala_class == "ward" & hospdays>90 ~ 90*234 + (hospdays-90)*191,
                                ea == "70" & pala_class == "ward" & hospdays<=90 ~hospdays*408,
                                ea == "70" & pala_class == "ward" & hospdays>90 ~ 90*408 + (hospdays-90)*408*0.816,
                                ea == "74" & pala_class == "ward" & hospdays<=90 ~hospdays*568,
                                ea == "74" & pala_class == "ward" & hospdays>90 ~ 90*568 + (hospdays-90)*568*0.816,
                                ea == "75" & pala_class == "ward" & hospdays<=90 ~hospdays*638,
                                ea == "75" & pala_class == "ward" & hospdays>90 ~ 90*638 + (hospdays-90)*638*0.816,
                                TRUE ~ unit_cost_w_fixed),
          unit_cost = case_when(ea == "98" & pala_class == "ward" & hospdays<=90 ~hospdays*234,
                                ea == "98" & pala_class == "ward" & hospdays>90 ~ 90*234 + (hospdays-90)*191,
                                ea == "70" & pala_class == "ward" & hospdays<=90 ~hospdays*408,
                                ea == "70" & pala_class == "ward" & hospdays>90 ~ 90*408 + (hospdays-90)*408*0.816,
                                ea == "74" & pala_class == "ward" & hospdays<=90 ~hospdays*568,
                                ea == "74" & pala_class == "ward" & hospdays>90 ~ 90*568 + (hospdays-90)*568*0.816,
                                ea == "75" & pala_class == "ward" & hospdays<=90 ~hospdays*638,
                                ea == "75" & pala_class == "ward" & hospdays>90 ~ 90*638 + (hospdays-90)*638*0.816,
                                TRUE ~ unit_cost))


#convert 2011 euro estimates to 2018 euros using price index of public expenditure
dhilmo$unit_cost <- dhilmo$unit_cost/0.9341038

names(dhilmo)[1] <- "FINNGENID"

# load additional phenotype file to get year of birth and death
 
endbig <- fread("/home/jsjukara/finngen_R5_V2_endpoint.txt", data.table = FALSE)

keepcols_endbig <- c("FINNGENID", "BL_AGE", "BL_YEAR", "SEX",
                     'DEATH','DEATH_AGE','DEATH_YEAR')

fg_df <- endbig[,keepcols_endbig]


fg_df$YEAR_OF_BIRTH <- fg_df$BL_YEAR - fg_df$BL_AGE

# make end of follow-up be death or 2019 if did not die during follow-up
fg_df$FU_END <- ifelse(is.na(fg_df$DEATH_YEAR), 2019, fg_df$DEATH_YEAR)

## Calculate individual sums of costs for overall data
# calculate sum of costs by id for inpatient and outpatient data from HILMO after 1998
indiv_costs_post1998 <- dhilmo %>%
    filter(year >= 1998) %>%
    group_by(FINNGENID) %>%
    summarize(HILMO_COSTS_POST1998=sum(unit_cost))

# calculate sum of costs by id for inpatient data only from HILMO for whole follow-up
indiv_costs_inpatient <- dhilmo %>%
    filter(source == "INPAT") %>%
    group_by(FINNGENID) %>%
    summarize(HILMO_COSTS_INPATIENT=sum(unit_cost))

# create data frame containing all ids
indiv_costs <- data.frame(FINNGENID = unique(fg_df$FINNGENID))

# join costs
indiv_costs <- left_join(indiv_costs, indiv_costs_post1998)

# make post1998 costs for IDs that died before 1998 NA, make NA costs for others 0.
ids_pre1998 <- fg_df$FINNGENID[fg_df$DEATH_YEAR<=1998]
indiv_costs$HILMO_COSTS_POST1998 <- ifelse(indiv_costs$FINNGENID %in% ids_pre1998, NA,
                                          ifelse(is.na(indiv_costs$HILMO_COSTS_POST1998), 0, indiv_costs$HILMO_COSTS_POST1998))


indiv_costs <- left_join(indiv_costs, indiv_costs_inpatient)
                  
# make inpatient costs for IDs that died before 1972 NA, make NA costs for others 0.
ids_pre1972 <- fg_df$FINNGENID[fg_df$DEATH_YEAR<=1972]
indiv_costs$HILMO_COSTS_INPATIENT <- ifelse(indiv_costs$FINNGENID %in% ids_pre1972, NA,
                                          ifelse(is.na(indiv_costs$HILMO_COSTS_INPATIENT), 0, indiv_costs$HILMO_COSTS_INPATIENT))

# save costs
fwrite(indiv_costs, file="hilmo_costs/finngen_R5_overall_hilmo_costs.txt.gz")


## Calculate sums of costs by age group
# join year of birth to hilmo and calculate age at event
dhilmo <- left_join(dhilmo, fg_df[,c("FINNGENID", "YEAR_OF_BIRTH")], by="FINNGENID")

dhilmo$age_at_event <- dhilmo$year - dhilmo$YEAR_OF_BIRTH



dhilmo_age <- dhilmo %>% mutate(age_group = case_when(age_at_event<10 ~ "0_10",
                                                 age_at_event<20 ~ "10_20",
                                                 age_at_event<30 ~ "20_30",
                                                 age_at_event<40 ~ "30_40",
                                                 age_at_event<50 ~ "40_50",
                                                 age_at_event<60 ~ "50_60",
                                                 age_at_event<70 ~ "60_70",
                                                 age_at_event<80 ~ "70_80",
                                                 age_at_event<90 ~ "80_90"),
                               post1998 = case_when(year>1998 ~ 1,
                                                   year<=1998 ~ 0)) %>%
    filter(!is.na(age_group),
          source == "INPAT") %>% 
    group_by(FINNGENID, age_group) %>%
    summarize(HILMO_COSTS=sum(unit_cost))

dhilmo_age <- dhilmo_age %>%
    pivot_wider(names_from="age_group", values_from="HILMO_COSTS")

dhilmo_age <- dhilmo_age[,c("FINNGENID", sort(names(dhilmo_age[2:10])))]

# make data frame with all FINNGENID's, even those without any costs in the data
dhilmo_age_full <- data.frame(FINNGENID = unique(fg_df$FINNGENID))
dhilmo_age_full <- left_join(dhilmo_age_full, dhilmo_age, by="FINNGENID")
dhilmo_age_full <- left_join(dhilmo_age_full, fg_df[,c("FINNGENID", "YEAR_OF_BIRTH", "DEATH_YEAR")])

dhilmo_age <- dhilmo_age_full

age_low <- seq(0, 80, 10)

# convert temporarily to matrix for fast assignment (at least 1000x faster)
dhilmo_agem <- as.matrix(dhilmo_age)

# loop to make costs that are NA to 0 if patient was followed-up during that age interval
for (i in 1:nrow(dhilmo_age)) {  
    age_start <- 1972 - dhilmo_age$YEAR_OF_BIRTH[i]  # can be negative
    
    # calculate age at end of follow-up/death
    yod <- ifelse(is.na(dhilmo_age$DEATH_YEAR[i]), 2019, dhilmo_age$DEATH_YEAR[i]) # if no death recorded, assume complete follow-up until 2019
    age_end <- yod - dhilmo_age$YEAR_OF_BIRTH[i]

    # set NA costs to 0 if follow-up started before upper bound of age interval
    # and follow-up ended after the lower bound of age interval
    indices <- which(age_start < age_low + 10 & age_end > age_low & is.na(dhilmo_agem[i,2:10])) + 1
    dhilmo_agem[i,indices] <- 0
    
    # set costs to NA for age interval if follow-up started after upper bound of age interval
    indices <- which(age_start > age_low + 10) + 1
    dhilmo_agem[i,indices] <- NA
}

# convert back to data frame
dhilmo_age <- as.data.frame(dhilmo_agem, stringsAsFactors=FALSE)
fwrite(dhilmo_age[,1:10], file="hilmo_costs/finngen_R5_hilmo_inpatient_costs_by_age.txt.gz")



dhilmo_age <- dhilmo %>% mutate(age_group = case_when(age_at_event<10 ~ "0_10",
                                                 age_at_event<20 ~ "10_20",
                                                 age_at_event<30 ~ "20_30",
                                                 age_at_event<40 ~ "30_40",
                                                 age_at_event<50 ~ "40_50",
                                                 age_at_event<60 ~ "50_60",
                                                 age_at_event<70 ~ "60_70",
                                                 age_at_event<80 ~ "70_80",
                                                 age_at_event<90 ~ "80_90"),
                               post1998 = case_when(year>1998 ~ 1,
                                                   year<=1998 ~ 0)) %>%
    filter(!is.na(age_group),
          year >= 1998) %>% 
    group_by(FINNGENID, age_group) %>%
    summarize(HILMO_COSTS=sum(unit_cost))

dhilmo_age <- dhilmo_age %>%
    pivot_wider(names_from="age_group", values_from="HILMO_COSTS")

dhilmo_age <- dhilmo_age[,c("FINNGENID", sort(names(dhilmo_age[2:10])))]

# make data frame with all FINNGENID's, even those without any costs in the data
dhilmo_age_full <- data.frame(FINNGENID = unique(fg_df$FINNGENID))
dhilmo_age_full <- left_join(dhilmo_age_full, dhilmo_age, by="FINNGENID")
dhilmo_age_full <- left_join(dhilmo_age_full, fg_df[,c("FINNGENID", "YEAR_OF_BIRTH", "DEATH_YEAR")])

dhilmo_age <- dhilmo_age_full

# convert temporarily to matrix for fast assignment (at least 1000x faster)
dhilmo_agem <- as.matrix(dhilmo_age)

# loop to make costs that are NA to 0 if patient was followed-up during that age interval
for (i in 1:nrow(dhilmo_age)) {  
    age_start <- 1998 - dhilmo_age$YEAR_OF_BIRTH[i]  # can be negative
    
    # calculate age at end of follow-up/death
    yod <- ifelse(is.na(dhilmo_age$DEATH_YEAR[i]), 2019, dhilmo_age$DEATH_YEAR[i]) # if no death recorded, assume complete follow-up until 2019
    age_end <- yod - dhilmo_age$YEAR_OF_BIRTH[i]

    # set NA costs to 0 if follow-up started before upper bound of age interval
    # and follow-up ended after the lower bound of age interval
    indices <- which(age_start < age_low + 10 & age_end > age_low & is.na(dhilmo_agem[i,2:10])) + 1
    dhilmo_agem[i,indices] <- 0
    
    # set costs to NA for age interval if follow-up started after upper bound of age interval
    indices <- which(age_start > age_low + 10) + 1
    dhilmo_agem[i,indices] <- NA
}
# convert back to data frame
dhilmo_age <- as.data.frame(dhilmo_agem, stringsAsFactors=FALSE)
head(dhilmo_age)
fwrite(dhilmo_age[,1:10], file="hilmo_costs/finngen_R5_hilmo_costs_by_age_post1998.txt.gz")