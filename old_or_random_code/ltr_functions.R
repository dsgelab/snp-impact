## Define functions to estimate hazard ratios and lifetime risk

# function to go from cumulative incidence (probability) to incidence via exponential formula
ci_to_i <- function(x) {
    return((1/5)*log(1/(1-x)))
}

# function that accepts endpoint (FinnGen) name and snp, returns age_df
# function that accepts endpoint (FinnGen) name and snp, returns age_df
get_age_df <- function(endpoint, rsid, method, cohort, fg_df) {
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

        #temp[,"FU_START_AGE"] <- ifelse(temp[,"FU_START_AGE"]<0, 0, temp[,"FU_START_AGE"])

        temp$person_time <- temp[,endpoint_age]

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

get_fit <- function(endpoint, rsid, method, cohort, fg_df) {
  # define data frame by age groups that is to contain probability estimates in age groups for outcome and death
  endpoint_age <- paste(endpoint, "_AGE", sep="")
  fg_df <- fg_df %>% select(FINNGENID, YEAR_OF_BIRTH, FU_START_AGE, FU_END_AGE, AGE_AT_DEATH_OR_NOW, DEATH, DEATH_AGE, endpoint, endpoint_age, rsid, FEMALE, PC1, PC2, PC3, PC4)
  temp <- fg_df
  
  # set outcome positive if outcome in interval
  temp <- temp %>% mutate(outcome = get(endpoint))

  # predict outcome
  form <- formula(paste("Surv(", endpoint_age, ", ", endpoint,") ~ ", rsid, " + ", covariates_formula, sep=""))
  fit_outcome <- suppressWarnings(coxph(formula = form,
                     data=temp))
  return(fit_outcome)
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
    gbd_outcome <- gbdfin %>% filter(cause_name == gbd_endpoint)

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