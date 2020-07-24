get_ltr_boot <- function(data,
                         #endpoint,
                         #gbd_endpoint,
                         #rsid,
                         #method,
                         #cohort,
                         indices,
                         hr_all) {
    d <- data[indices,]
    age_df <- get_age_df_se(endpoint=endpoint,
                     rsid=rsid,
                     method=method,
                     cohort=cohort,
                     fg_df=d)
    
    age_df$hr01 <- ifelse(is.na(age_df$hr01), hr_all, age_df$hr01)
    age_df$hr02 <- ifelse(is.na(age_df$hr02), hr_all^2, age_df$hr02)
    
    age_df$hr01 <- ifelse(age_df$hr01 > 10, hr_all, age_df$hr01)
    age_df$hr01 <- ifelse(age_df$hr01 < 0.1, hr_all, age_df$hr01)
    
    age_df$hr02 <- ifelse(age_df$hr02 > 10^2, hr_all^2, age_df$hr02)
    age_df$hr02 <- ifelse(age_df$hr02 < 0.1^2, hr_all^2, age_df$hr02)
    
    print(gbd_endpoint)
    print(age_df)
    print(gbd_subset)
    print(rsid)
    print(genotype_freqs)
    ltr <- get_ltr_incidence(outcome_name = gbd_endpoint,
                             age_df = age_df,
                             gbdfin=gbd_subset,
                             rsid=rsid,
                             genotype_freqs=genotype_freqs)
    
    return(ltr)
}