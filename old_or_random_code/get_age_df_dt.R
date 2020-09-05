get_age_df_dt <- function(endpoint, rsid, method, cohort, fg_df) {
  # define data frame by age groups that is to contain probability estimates in age groups for outcome and death
  age_df <- data.frame(age_name = c("Under 5", "5 to 9", "10 to 14", "15 to 19", "20 to 24", "25 to 29", "30 to 34",
                                    "35 to 39", "40 to 44", "45 to 49", "50 to 54", "55 to 59", "60 to 64", "65 to 69",
                                    "70 to 74", "75 to 79", "80 to 84"),
                       age_low = seq(0,80,5),
                       age_high = seq(5,85,5),
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

    if (method == "coxph" | method == "crude_incidence") {
        
        # loop over age groups
        for (i in 1:nrow(age_df)) {
            if (cohort == "post1998") {
              temp <- fg_df[1998 - YEAR_OF_BIRTH <= age_df$age_low[i] & # excl. if age at 1998 higher than start of interval
                       FU_END_AGE >= age_df$age_low[i] & #excl if FU end before start of interval
                       get(endpoint_age) >= age_df$age_low[i]] # excl. if outcome before interval
            } else {
              temp <- fg_df[FU_START_AGE <= age_df$age_low[i] & # excl. if age at 1998 higher than start of interval
                       FU_END_AGE >= age_df$age_low[i] & #excl if FU end before start of interval
                       get(endpoint_age) >= age_df$age_low[i]] # excl. if outcome before interval

            }
            #temp[,outcome1 := fcase(get(endpoint) == 1 &
            #                        get(endpoint_age) > age_df$age_low[i] &
            #                         get(endpoint_age) <= age_df$age_high[i], 1,
            #                         TRUE, 0)]
            # set outcome positive if outcome in interval
            temp[, outcome:= ifelse(get(endpoint) == 1 &
                                     get(endpoint_age) > age_df$age_low[i] &
                                     get(endpoint_age) <= age_df$age_high[i], 1, 0)]

            #temp <- temp %>% mutate(outcome = case_when(get(endpoint) == 1 &
            #                                           get(endpoint_age) > age_df$age_low[i] &
            #                                           get(endpoint_age) <= age_df$age_high[i] ~ 1,
            #                                           TRUE ~ 0)
            #                                           ) # hack to get else condition working

            # set futime

            #temp <- temp %>% mutate(futime = case_when(outcome == 1 ~ get(endpoint_age),
            #                                          TRUE ~ age_df$age_high[i])
            #)

            temp[, futime:= ifelse(outcome == 1, get(endpoint_age), age_df$age_high[i])]

            # set death as positive in interval
            #temp <- temp %>% mutate(death = case_when(DEATH == 1 &
            #                                          DEATH_AGE >= age_df$age_low[i] &
            #                                           DEATH_AGE < age_df$age_high[i] ~ 1,
            #                                           TRUE ~ 0)
            #                                           ) # hack to get else condition working

            temp[, person_time := futime - age_df$age_low[i]]
            
            #temp$person_time <- temp$futime - age_df$age_low[i]
            
            #if (sum(temp$person_time<0) != 0) {
            #  print("person time below 0, check")
            #  print(paste(endpoint, rsid))
            #  print(age_df$age_low[i])
            #  return(temp)
            #}

            #temp %>% mutate(person_time = case_when(person_time > 5 ~ 5,
            #                                        !is.na(FINNGENID) ~ person_time))

            age_df$person_time[i] <- temp[,sum(person_time)]
            #age_df$n_death[i] = sum(temp$death)
            age_df$n[i] = nrow(temp)
            age_df$n_outcome[i] = temp[,sum(outcome)]

            #age_df$person_time[i] <- sum(temp$person_time)
            #age_df$n_death[i] = sum(temp$death)
            #age_df$n[i] = nrow(temp)
            #age_df$n_outcome[i] = sum(temp$outcome)


            if (method == "coxph") {
              # predict outcome
              form <- formula(paste("Surv(futime, outcome) ~ ", rsid, " + ", covariates_formula, sep=""))
              fit_outcome <- suppressWarnings(coxph(formula = form,
                                 data=temp))
              # save HR
              age_df$hr01[i] <- exp(fit_outcome$coefficients[1])
              # estimate 95% CI
              if (!is.na(age_df$hr01[i])) {
                  age_df$hr01loci[i] <- summary(fit_outcome)$conf.int[1,3]
                  age_df$hr01hici[i] <- summary(fit_outcome)$conf.int[1,4]
              }

              age_df$hr02[i] <- exp(fit_outcome$coefficients[1]*2)

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
              if (!is.na(age_df$hr01[i])) {
                se <- sqrt(1/n0+1/n1)

                age_df$hr01loci[i] <- exp(log(age_df$hr01[i]) - 1.96*se)
                age_df$hr01hici[i] <- exp(log(age_df$hr01[i]) + 1.96*se)
              }

            }

            # if missing HRs, impute 1
            if (sum(is.na(age_df[,c("hr01", "hr02")]))!=0) {
                # print("Had to impute HRs to 1")
                age_df$hr01 <- ifelse(is.na(age_df$hr01), 1, age_df$hr01)
                age_df$hr02 <- ifelse(is.na(age_df$hr02), 1, age_df$hr02)
            }

        }
    } else if (method == "coxph_single" & cohort == "all") {
        temp <- fg_df
        
        # set outcome positive if outcome in interval
        temp <- temp %>% mutate(outcome = get(endpoint))
        # set death as positive in interval
        temp <- temp %>% mutate(death = DEATH)

        temp[,"FU_START_AGE"] <- ifelse(temp[,"FU_START_AGE"]<0, 0, temp[,"FU_START_AGE"])

        temp$person_time <- temp[,endpoint_age]
        if (sum(temp$person_time<0) != 0) {
          print("person time below 0, check")
          print(paste(endpoint, rsid))
          print(age_df$age_low[i])
          print(summary(temp$FU_START_AGE))
          print(summary(temp$person_time))
          return(temp)
        }

        age_df$person_time <- sum(temp$person_time)
        age_df$n_death = sum(temp$death)
        age_df$n = nrow(temp)
        age_df$n_outcome = sum(temp$outcome)


        # predict outcome
        form <- formula(paste("Surv(", endpoint_age, ", ", endpoint,") ~ ", rsid, " + ", covariates_formula, sep=""))
        fit_outcome <- suppressWarnings(coxph(formula = form,
                           data=temp))
        # save HRs
        age_df$hr01 <- exp(fit_outcome$coefficients[1])
        age_df$hr01loci <- summary(fit_outcome)$conf.int[1,3]
        age_df$hr01hici <- summary(fit_outcome)$conf.int[1,4]

        age_df$hr02 <- exp(fit_outcome$coefficients[1]*2)
    }
    age_df$raw_prob_death <- age_df$n_death / age_df$n
    age_df$raw_prob_outcome <- age_df$n_outcome / age_df$n
    age_df$outcome_incidence_fg <- age_df$n_outcome/age_df$person_time

    return(age_df)
}

selected_colnames <- c("FINNGENID", "YEAR_OF_BIRTH", "FU_START_AGE",
                       "FU_END_AGE", "AGE_AT_DEATH_OR_NOW", "DEATH",
                       "DEATH_AGE", "FEMALE", "PC1", "PC2", "PC3", "PC4")
fg_df_dt <- data.table(fg_df)
fitlist <- list()
tic()
for (i in start:stop) {
    endpoint <- ltr_df$endpoint_name[i]
    gbd_endpoint <- ltr_df$cause_name[i]
    endpoint_age <- paste(endpoint, "_AGE", sep="")
    rsid <- ltr_df$rsid[i]
    temp_selected_colnames <- c(selected_colnames, endpoint, endpoint_age, rsid)
    temp_fg_df <- fg_df_dt[,..temp_selected_colnames]
    age_df <- get_age_df_dt(endpoint=endpoint,
                         rsid=rsid,
                         method=ltr_df$method[i],
                         cohort=ltr_df$cohort[i],
                         fg_df=temp_fg_df)
    if (ltr_df$cohort[i] == "all" & ltr_df$method[i] == "coxph_single") {
        fitlist[[endpoint]] <- get_fit(endpoint=endpoint,
                             rsid=rsid,
                             method=ltr_df$method[i],
                             cohort=ltr_df$cohort[i],
                             fg_df=temp_fg_df)
        }
    ltr <- get_ltr_incidence(outcome_name = gbd_endpoint,
                             age_df = age_df,
                             gbdfin=gbd_subset,
                             rsid=rsid,
                             genotype_freqs=genotype_freqs)
    
    ltr_df[i, 6:8] <- ltr
    
    age_df$rsid <- rsid
    age_df$endpoint_name <- endpoint
    age_df$cause_name <- gbd_endpoint
    age_df$cohort <- ltr_df$cohort[i]
    age_df$method <- ltr_df$method[i]
    age_df_list[[paste(ltr_df$cohort[i], ltr_df$method[i], rsid, endpoint)]] <- age_df
}
toc()

#age_df_all <- bind_rows(age_df_list)