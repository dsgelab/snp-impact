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