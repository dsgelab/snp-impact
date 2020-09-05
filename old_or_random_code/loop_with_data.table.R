library(tictoc)
selected_colnames <- c("FINNGENID", "YEAR_OF_BIRTH", "FU_START_AGE",
                       "FU_END_AGE", "AGE_AT_DEATH_OR_NOW", "DEATH",
                       "DEATH_AGE", "FEMALE", "PC1", "PC2", "PC3", "PC4")
fg_df_dt <- data.table(fg_df)
fitlist <- list()
tic()
for (i in 1:nrow(ltr_df)) {
    endpoint <- ltr_df$endpoint_name[i]
    gbd_endpoint <- ltr_df$cause_name[i]
    endpoint_age <- paste(endpoint, "_AGE", sep="")
    rsid <- ltr_df$rsid[i]
    temp_selected_colnames <- c(selected_colnames, endpoint, endpoint_age, rsid)
    temp_fg_df <- fg_df_dt[,..temp_selected_colnames]
    age_df <- get_age_df(endpoint=endpoint,
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