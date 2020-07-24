m_vec <- seq(1,12,1)
age_df_all <- bind_rows(age_df_list)
age_df_all$mult <- NA
age_df_all$hr01loci <- exp(age_df_all$loghr - 1.96 * age_df_all$loghr_se)
age_df_all$hr01hici <- exp(age_df_all$loghr + 1.96 * age_df_all$loghr_se)

adjust_hr_bayes <- function(age_df_all, mult) {
    # Bayesian adjustment of age group HR's
    bayes <- age_df_all %>% filter(method=="coxph")

    #get normally distributed prior distribution from coxph_single model
    temp <- age_df_all %>% filter(method=="coxph_single")

    #assign priors
    m.pr <- temp$loghr
    s.pr <- temp$loghr_se * mult
    #assign observed values
    m.obs <- bayes$loghr
    s.obs <- bayes$loghr_se
    
    #calculate posteriors
    s.po <- 1/sqrt(1/s.pr^2 + 1/s.obs^2)
    m.po <- s.po^2 * (m.pr/s.pr^2 + m.obs/s.obs^2)

    #if posteriors are NA (from missing observed values), use priors
    m.po <- ifelse(is.na(m.po), m.pr, m.po)
    s.po <- ifelse(is.na(s.po), s.pr, s.po)
    
    print(m.po)
    print(s.po)
    bayes$loghr <- m.po
    bayes$loghr_se <- s.po

    #calculate HR and 95% CI
    bayes$hr01 <- exp(m.po)
    bayes$hr01loci <- exp(log(bayes$hr01) - 1.96*s.po)
    bayes$hr01hici <- exp(log(bayes$hr01) + 1.96*s.po)
    
    #name the method
    bayes$method <- paste("coxph_bayes", mult, sep="")
    
    return(bayes)
}

