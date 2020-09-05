#save.image(file = "ukbb_workspace.RData", compress=FALSE)
## Compare incidences in FinnGen vs UKBB
fg_df <- readRDS("/home/jsjukara/ukbb/fg_df.rds")
endpoint_map <- readRDS("/home/jsjukara/ukbb/endpoint_map.rds")
# have to read ukbb data as well

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



get_incidences <- function(endpoint, method, cohort, fg_df) {
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

    age_df$raw_prob_death <- age_df$n_death / age_df$n
    age_df$raw_prob_outcome <- age_df$n_outcome / age_df$n
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