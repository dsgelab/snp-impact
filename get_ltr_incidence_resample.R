#define function to output lifetime risk given age_df, version for incidence method
get_ltr_incidence_resample <- function(outcome_name,
                              age_df,
                              gbdfin,
                              rsid,
                              genotype_freqs,
                              B=1000) {

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
  
  n_agegroups <- nrow(age_df)
  
  ltr_resample_df <- data.frame(b=1:B, p0=NA, p1=NA, p2=NA)
  
  for (b in 1:B) {
    loghr_resample <- rnorm(n = nrow(age_df), mean = age_df$loghr, sd=age_df$loghr_se)

    hr01_resample <- exp(loghr_resample)
    hr02_resample <- exp(loghr_resample*2)


    # estimate incidence attributable to 0, 1, 2 alleles
    age_df$i0 <- (age_df$outcome_incidence*age_df$population) / (props[1] * age_df$population + hr01_resample * props[2] * age_df$population + hr02_resample * props[3] * age_df$population)
    age_df$i1 <- age_df$i0 * hr01_resample
    age_df$i2 <- age_df$i0 * hr02_resample
    
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
      tab$a[n_agegroups] <- tab$l[n_agegroups] * (tab$r[n_agegroups]/(tab$r[n_agegroups] + tab$m[n_agegroups] - tab$mo[n_agegroups]))
      tab$d[n_agegroups] <- tab$l[n_agegroups] * ((tab$m[n_agegroups] - tab$mo[n_agegroups])/(tab$r[n_agegroups] + tab$m[n_agegroups] - tab$mo[n_agegroups]))

      # Calculate lifetime risk
      ltrvec[dosage+1] <- sum(tab$a)/tab$l[1]
    }

  ltr_resample_df[b,2:4] <- ltrvec
  }
  
  return(ltr_resample_df)
}