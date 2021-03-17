get_attr_individual_dalys_boot_incremental_dementia <- function(data,
                         endpoints,
                         rsid,
                         indices,
                         mult=6) {
    
    bootstrap_ltr_list <- list()
    
    d <- data[indices,]
    
    ltr_df <- expand.grid(rsid=rsid,
                endpoint_name=endpoints,
                stringsAsFactors=FALSE)
    
    temp <- endpoint_map
    temp$endpoint_name <- temp$endpoint
    temp$cause_name <- temp$endpoint_gbd
    ltr_df <- left_join(ltr_df, temp[,c("endpoint_name", "cause_name")], by="endpoint_name")
    
    ltr_df$ltr0 <- NA
    ltr_df$ltr1 <- NA
    ltr_df$ltr2 <- NA
    ltr_df <- left_join(ltr_df, dalys_both[,c("cause_name", "val")], by="cause_name")
    names(ltr_df)[grep("val", names(ltr_df))] <- "population_dalys"
    
    ltr_df <- left_join(ltr_df, genotype_freqs, by="rsid")
    
    for (i in 1:nrow(ltr_df)) {
        endpoint <- ltr_df$endpoint_name[i]
        endpoint_age <- paste(endpoint, "_AGE", sep="")
        gbd_endpoint <- endpoint_map$endpoint_gbd[endpoint_map$endpoint == endpoint]

        selected_colnames <- c("FINNGENID", "YEAR_OF_BIRTH", "FU_START_AGE",
                           "FU_END_AGE", "AGE_AT_DEATH_OR_NOW", "DEATH",
                           "DEATH_AGE", "FEMALE", "PC1", "PC2", "PC3", "PC4")
        
        temp_selected_colnames <- c(selected_colnames, endpoint, endpoint_age, rsid)
        
        dtemp <- d[,temp_selected_colnames]
        
        age_df_coxph <- get_age_df_se(endpoint=ltr_df$endpoint_name[i],
                         rsid=rsid,
                         method="coxph",
                         cohort="all",
                         fg_df=dtemp)
        
        age_df_coxph$method <- "coxph"

        age_df_coxph_single <- get_age_df_se(endpoint=ltr_df$endpoint_name[i],
                         rsid=rsid,
                         method="coxph_single",
                         cohort="all",
                         fg_df=d)

        age_df_coxph_single$method <- "coxph_single"

        age_df_all <- rbind(age_df_coxph, age_df_coxph_single)

        temp <- adjust_hr_bayes(age_df_all, mult)

        age_df_all <- rbind(age_df_all, temp)

        age_df_all$hr02 <- exp(age_df_all$loghr*2)

        age_df <- age_df_all %>% filter(method==paste("coxph_bayes", mult, sep=""))

        ltr <- get_ltr_incidence(outcome_name = gbd_endpoint,
                                 age_df = age_df,
                                 gbdfin=gbd_subset,
                                 rsid=rsid,
                                 genotype_freqs=genotype_freqs)
        ltr_df[i,4:6] <- ltr
        
    }
    
    n_population <- 5511371
    ltr_df$n0 <- n_population*ltr_df$p0
    ltr_df$n1 <- n_population*ltr_df$p1
    ltr_df$n2 <- n_population*ltr_df$p2
    ltr_df$denom <- (ltr_df$ltr0 * ltr_df$p0 + ltr_df$ltr1 * ltr_df$p1 + ltr_df$ltr2 * ltr_df$p2)
    ltr_df$pop_dalys0 <- ltr_df$ltr0*ltr_df$p0 / ltr_df$denom * ltr_df$population_dalys
    ltr_df$pop_dalys1 <- ltr_df$ltr1*ltr_df$p1 / ltr_df$denom * ltr_df$population_dalys
    ltr_df$pop_dalys2 <- ltr_df$ltr2*ltr_df$p2 / ltr_df$denom * ltr_df$population_dalys
    ltr_df$person_dalys0 <- ltr_df$pop_dalys0 / ltr_df$n0
    ltr_df$person_dalys1 <- ltr_df$pop_dalys1 / ltr_df$n1
    ltr_df$person_dalys2 <- ltr_df$pop_dalys2 / ltr_df$n2

    ltr_df$attr_person_dalys1 <- ltr_df$person_dalys1 - ltr_df$person_dalys0
    ltr_df$attr_person_dalys2 <- ltr_df$person_dalys2 - ltr_df$person_dalys0

    ltr_df$attr_population_dalys1 <- ltr_df$attr_person_dalys1 * ltr_df$n1
    ltr_df$attr_population_dalys2 <- ltr_df$attr_person_dalys2 * ltr_df$n2

    ltr_df$person_attr1_lifetime <- ltr_df$attr_person_dalys1*82
    ltr_df$person_attr2_lifetime <- ltr_df$attr_person_dalys2*82
    ltr_df$attr_population_dalys12 <- ltr_df$attr_population_dalys1 + ltr_df$attr_population_dalys2
    
    # permute rows to randomize order of endpoints
    ltr_df <- ltr_df[sample(1:nrow(ltr_df)),]
    # swap dementia to top
    dementiarow <- which(ltr_df$endpoint_name == "KRA_PSY_DEMENTIA")
    save <- ltr_df[1,]
    ltr_df[1,]<- ltr_df[dementiarow,]
    ltr_df[dementiarow,] <- save
    
    
    
    indiv <- cumsum(ltr_df$person_attr1_lifetime)
    pop <- cumsum(ltr_df$attr_population_dalys12)
    
    results <- c(indiv, pop)
    
    return(results)
}

mult <- 6
B <- 1000
endpoints <- unique(ltr_df$endpoint_name)
library(tictoc)
tic()
results <- boot(data=fg_df,
                endpoints=endpoints,
                rsid="chr19_44908684_T_C",
                mult=6,
                statistic=get_attr_individual_dalys_boot_incremental_dementia,
                R=B,
                parallel = "multicore", ncpus=20)
toc()

estimates_df <- data.frame(measure=c(rep(paste("indiv", 1:35, sep="_")),
                                    rep(paste("pop", 1:35, sep="_"))),
                          estimate=NA,
                          loci=NA,
                          hici=NA)

temp <- results$t
for (j in 1:ncol(temp)) {
    temp[,j] <- as.numeric(temp[,j])
}

temp <- temp[!is.na(temp[,1]),]

for (j in 1:ncol(temp)) {
    estimates_df[j, 3:4] <- quantile(as.numeric(temp[,j]), probs=c(0.025, 0.975))
    estimates_df[j, 2] <- mean(as.numeric(temp[,j]))
}

estimates_df$measure <- as.character(estimates_df$measure)
estimates_df$n <- rep(1:35,2)
indiv <- estimates_df[1:35,]
pop <- estimates_df[36:70,]

estimates_df$measure <- as.character(estimates_df$measure)
estimates_df$n <- rep(1:35,2)
indiv <- estimates_df[1:35,]
pop <- estimates_df[36:70,]

options(repr.plot.width=12, repr.plot.height=8)

indiv %>%
    ggplot(aes(x=n, y=estimate)) +
    geom_errorbar(aes(x=n, ymin=loci, ymax=hici)) +
    geom_point() +
    geom_hline(yintercept=0) +
    labs(y="Total attributable individual lifetime DALYS") +
    theme_bw()


options(repr.plot.width=12, repr.plot.height=8)

pop %>%
    ggplot(aes(x=n, y=estimate)) +
    geom_errorbar(aes(x=n, ymin=loci, ymax=hici)) +
    geom_point() +
    geom_hline(yintercept=0) +
    labs(y="Total attributable population DALYs") +
    theme_bw()
    
#######


library(boot)
# individual bootstrap to get attributable individual dalys 

get_attr_individual_dalys_boot <- function(data,
                         endpoints,
                         rsid,
                         indices,
                         mult=6) {
    
    bootstrap_ltr_list <- list()
    
    d <- data[indices,]
    
    ltr_df <- expand.grid(rsid=rsid,
                endpoint_name=endpoints,
                stringsAsFactors=FALSE)
    
    temp <- endpoint_map
    temp$endpoint_name <- temp$endpoint
    temp$cause_name <- temp$endpoint_gbd
    ltr_df <- left_join(ltr_df, temp[,c("endpoint_name", "cause_name")], by="endpoint_name")
    
    ltr_df$ltr0 <- NA
    ltr_df$ltr1 <- NA
    ltr_df$ltr2 <- NA
    ltr_df <- left_join(ltr_df, dalys_both[,c("cause_name", "val")], by="cause_name")
    names(ltr_df)[grep("val", names(ltr_df))] <- "population_dalys"
    
    ltr_df <- left_join(ltr_df, genotype_freqs, by="rsid")
    
    for (i in 1:nrow(ltr_df)) {
        endpoint <- ltr_df$endpoint_name[i]
        endpoint_age <- paste(endpoint, "_AGE", sep="")
        gbd_endpoint <- endpoint_map$endpoint_gbd[endpoint_map$endpoint == endpoint]

        selected_colnames <- c("FINNGENID", "YEAR_OF_BIRTH", "FU_START_AGE",
                           "FU_END_AGE", "AGE_AT_DEATH_OR_NOW", "DEATH",
                           "DEATH_AGE", "FEMALE", "PC1", "PC2", "PC3", "PC4")
        
        temp_selected_colnames <- c(selected_colnames, endpoint, endpoint_age, rsid)
        
        dtemp <- d[,temp_selected_colnames]
        
        age_df_coxph <- get_age_df_se(endpoint=ltr_df$endpoint_name[i],
                         rsid=rsid,
                         method="coxph",
                         cohort="all",
                         fg_df=dtemp)
        
        age_df_coxph$method <- "coxph"

        age_df_coxph_single <- get_age_df_se(endpoint=ltr_df$endpoint_name[i],
                         rsid=rsid,
                         method="coxph_single",
                         cohort="all",
                         fg_df=d)

        age_df_coxph_single$method <- "coxph_single"

        age_df_all <- rbind(age_df_coxph, age_df_coxph_single)

        temp <- adjust_hr_bayes(age_df_all, mult)

        age_df_all <- rbind(age_df_all, temp)

        age_df_all$hr02 <- exp(age_df_all$loghr*2)

        age_df <- age_df_all %>% filter(method==paste("coxph_bayes", mult, sep=""))

        ltr <- get_ltr_incidence(outcome_name = gbd_endpoint,
                                 age_df = age_df,
                                 gbdfin=gbd_subset,
                                 rsid=rsid,
                                 genotype_freqs=genotype_freqs)
        ltr_df[i,4:6] <- ltr
        
    }
    
    n_population <- 5511371
    ltr_df$n0 <- n_population*ltr_df$p0
    ltr_df$n1 <- n_population*ltr_df$p1
    ltr_df$n2 <- n_population*ltr_df$p2
    ltr_df$denom <- (ltr_df$ltr0 * ltr_df$p0 + ltr_df$ltr1 * ltr_df$p1 + ltr_df$ltr2 * ltr_df$p2)
    ltr_df$pop_dalys0 <- ltr_df$ltr0*ltr_df$p0 / ltr_df$denom * ltr_df$population_dalys
    ltr_df$pop_dalys1 <- ltr_df$ltr1*ltr_df$p1 / ltr_df$denom * ltr_df$population_dalys
    ltr_df$pop_dalys2 <- ltr_df$ltr2*ltr_df$p2 / ltr_df$denom * ltr_df$population_dalys
    ltr_df$person_dalys0 <- ltr_df$pop_dalys0 / ltr_df$n0
    ltr_df$person_dalys1 <- ltr_df$pop_dalys1 / ltr_df$n1
    ltr_df$person_dalys2 <- ltr_df$pop_dalys2 / ltr_df$n2

    ltr_df$attr_person_dalys1 <- ltr_df$person_dalys1 - ltr_df$person_dalys0
    ltr_df$attr_person_dalys2 <- ltr_df$person_dalys2 - ltr_df$person_dalys0

    ltr_df$attr_population_dalys1 <- ltr_df$attr_person_dalys1 * ltr_df$n1
    ltr_df$attr_population_dalys2 <- ltr_df$attr_person_dalys2 * ltr_df$n2

    ltr_df$person_attr1_lifetime <- ltr_df$attr_person_dalys1*82
    ltr_df$person_attr2_lifetime <- ltr_df$attr_person_dalys2*82
    ltr_df$attr_population_dalys12 <- ltr_df$attr_population_dalys1 + ltr_df$attr_population_dalys2
    
    indiv <- sum(ltr_df$person_attr1_lifetime)#, na.rm=T)
    pop <- sum(ltr_df$attr_population_dalys12)#, na.rm=T)
    
    results <- c(indiv, pop, ltr_df$person_attr1_lifetime, ltr_df$attr_population_dalys12)
    
    return(results)
}

test_rsids <- c("chr19_44908684_T_C",
                "chr8_142684641_C_G",
                "chr2_55812900_A_T",
                "chr15_78618839_T_C",
                "chr15_78576294_G_A",
                "chr12_20424871_T_C",
                "chr3_146516482_A_G")

mult <- 6
B <- 500
endpoints <- endpoint_map$endpoint

tic()
results <- list(NULL)
for (current_rsid in test_rsids) {
    results[[current_rsid]] <- boot(data=fg_df,
                endpoints=endpoints,
                rsid=current_rsid,
                mult=6,
                statistic=get_attr_individual_dalys_boot,
                R=B,
                parallel = "multicore", ncpus=63)
}

toc()

estimates_df <- data.frame(measure=c("sum_indiv_dalys",
                                    "sum_pop_dalys",
                                    paste(endpoint_map$endpoint, "_", "indiv", sep=""),
                                    paste(endpoint_map$endpoint, "_", "pop", sep="")),
                          estimate=NA,
                          loci=NA,
                          hici=NA)
for (i in 1:length(test_rsids)) {
    temp <- results[[test_rsids[i]]]$t
    for (j in 1:ncol(temp)) {
        temp[,j] <- as.numeric(temp[,j])
    }
    
    temp <- temp[!is.na(temp[,1]),]
    
    for (j in 1:ncol(temp)) {
        estimates_df[j, 3:4] <- quantile(as.numeric(temp[,j]), probs=c(0.025, 0.975))
        estimates_df[j, 2] <- mean(as.numeric(temp[,j]))
    }
    estimates_df$rsid <- test_rsids[i]
    if (i == 1) {
        estimates_df_all <- estimates_df
    } else {
        estimates_df_all <- rbind(estimates_df_all, estimates_df)
    }
}

estimates_df_all %>%
    filter(measure %in% paste(endpoint_map$endpoint, "_", "indiv", sep=""),
          rsid=="chr15_78576294_G_A") %>%
    arrange(desc(estimate))

estimates_df_all %>%
    filter(measure %in% c(paste(endpoint_map$endpoint, "_", "pop", sep=""), "sum_pop_dalys" ),
          rsid=="chr15_78576294_G_A") %>%
    arrange(desc(estimate))


options(repr.plot.width=12, repr.plot.height=8)

estimates_df_all %>%
    filter(measure=="sum_indiv_dalys") %>%
    ggplot(aes(x=rsid, y=estimate)) +
    geom_errorbar(aes(x=rsid, ymin=loci, ymax=hici)) +
    geom_point() +
    geom_hline(yintercept=0) +
    labs(y="Total attributable individual lifetime DALYS") +
    theme_bw()
    

estimates_df_all %>%
    filter(measure=="sum_pop_dalys") %>%
    ggplot(aes(x=rsid, y=estimate)) +
    geom_errorbar(aes(x=rsid, ymin=loci, ymax=hici)) +
    geom_point() +
    geom_hline(yintercept=0) +
    labs(y="Total attributable population DALYS") +
    theme_bw()

options(repr.plot.width=9, repr.plot.height=6)

estimates_df_all %>%
    filter(measure %in% paste(endpoint_map$endpoint, "_", "indiv", sep=""),
          rsid=="chr19_44908684_T_C") %>%
    ggplot(aes(x=reorder(measure, estimate), y=estimate)) +
    geom_errorbar(aes(x=reorder(measure, estimate), ymin=loci, ymax=hici)) +
    geom_point() +
    geom_hline(yintercept=0) +
    labs(y="Attributable individual lifetime DALYs") +
    theme_bw() +
    coord_flip() +
    labs(x="")

options(repr.plot.width=9, repr.plot.height=6)

estimates_df_all %>%
    filter(measure %in% paste(endpoint_map$endpoint, "_", "indiv", sep=""),
          rsid=="chr15_78618839_T_C") %>%
    ggplot(aes(x=reorder(measure, estimate), y=estimate)) +
    geom_errorbar(aes(x=reorder(measure, estimate), ymin=loci, ymax=hici)) +
    geom_point() +
    geom_hline(yintercept=0) +
    labs(y="Attributable individual lifetime DALYs") +
    theme_bw() +
    coord_flip() +
    labs(x="", title=="")

for (rs in test_rsids) {
    plot(
    estimates_df_all %>%
        filter(measure %in% paste(endpoint_map$endpoint, "_", "indiv", sep=""),
              rsid==rs) %>%
        ggplot(aes(x=measure, y=estimate)) +
        geom_errorbar(aes(x=measure, ymin=loci, ymax=hici)) +
        geom_point() +
        geom_hline(yintercept=0) +
        labs(y="Total attributable population DALYS", title=rs) +
        theme_bw()
    )
}

indiv_dalys <- as.numeric(temp$t[,1])
indiv_dalys <- indiv_dalys[!is.na(indiv_dalys)]
indiv_dalys <- indiv_dalys[!is.na(indiv_dalys)]

pop_dalys <- as.numeric(temp$t[,2])
pop_dalys <- pop_dalys[!is.na(pop_dalys)]
pop_dalys <- pop_dalys[!is.na(pop_dalys)]

endpoint_indiv_dalys <- temp$t[,3:(nrow(endpoint_map)+2)]

endpoint_pop_dalys <- temp$t[,(nrow(endpoint_map)+3):(ncol(temp$t))]

indiv_dalys <- as.numeric(results$t[,1])
indiv_dalys <- indiv_dalys[!is.na(indiv_dalys)]

mult <- 6
B <- 300
endpoints <- endpoint_map$endpoint

tic()
results1 <- boot(data=fg_df,
                endpoints="KRA_PSY_DEMENTIA",
                rsid="chr19_44908684_T_C",
                mult=6,
                statistic=get_attr_individual_dalys_boot,
                R=B,
                parallel = "multicore", ncpus=63)
toc()


######

library(boot)

# individual bootstrap  to get lifetime risks
get_ltr_boot <- function(data,
                         #endpoint,
                         #gbd_endpoint,
                         #rsid,
                         #method,
                         #cohort,
                         indices,
                         mult=8) {
    d <- data[indices,]
    age_df_coxph <- get_age_df_se(endpoint=endpoint,
                     rsid=rsid,
                     method="coxph",
                     cohort=cohort,
                     fg_df=d)
    age_df_coxph$method <- "coxph"

    age_df_coxph_single <- get_age_df_se(endpoint=endpoint,
                     rsid=rsid,
                     method="coxph_single",
                     cohort=cohort,
                     fg_df=d)
    
    age_df_coxph_single$method <- "coxph_single"

    age_df_all <- rbind(age_df_coxph, age_df_coxph_single)

    #print(head(age_df_all))
    #print(tail(age_df_all))
    temp <- adjust_hr_bayes(age_df_all, mult)
    #temp$mult <- mult

    age_df_all <- rbind(age_df_all, temp)

    age_df_all$hr02 <- exp(age_df_all$loghr*2)

    age_df <- age_df_all %>% filter(method==paste("coxph_bayes", mult, sep=""))
    
    
    ltr <- get_ltr_incidence(outcome_name = gbd_endpoint,
                             age_df = age_df,
                             gbdfin=gbd_subset,
                             rsid=rsid,
                             genotype_freqs=genotype_freqs)
    
    return(ltr)
}


## Bootstrap

options(repr.plot.width=12, repr.plot.height=7)
mult <- 6
B <- 300
endpoints <- c("KRA_PSY_DEMENTIA", "G6_MS", "C3_PROSTATE")
snps <- c("chr19_44908684_T_C", "chr13_19556224_G_A")

loop_df <- expand.grid(endpoints,snps, stringsAsFactors=FALSE)
names(loop_df) <- c("endpoint", "rsid")

bootstrap_ltr_list <- list()
bootstrap_boot_list <- list()
tic()
for (i in 1:nrow(loop_df)) {
    endpoint <- loop_df$endpoint[i]
    gbd_endpoint <- endpoint_map$endpoint_gbd[endpoint_map$endpoint == endpoint]
    endpoint_age <- paste(endpoint, "_AGE", sep="")
    
    cohort <- "all"
    method <- "coxph"
    rsid <- loop_df$rsid[i]
    rs <- rsid
    
    selected_colnames <- c("FINNGENID", "YEAR_OF_BIRTH", "FU_START_AGE",
                       "FU_END_AGE", "AGE_AT_DEATH_OR_NOW", "DEATH",
                       "DEATH_AGE", "FEMALE", "PC1", "PC2", "PC3", "PC4")
    temp_selected_colnames <- c(selected_colnames, endpoint, endpoint_age, rsid)
    temp_fg_df <- fg_df[,temp_selected_colnames]
    
    results <- boot(data=temp_fg_df,
                #hr_all=hr_all,
                mult=8,
                statistic=get_ltr_boot,
                R=B,
                parallel = "multicore", ncpus=30)
    
    bootstrap_boot_list[[rsid]][[endpoint]] <- results 

    temp <- data.frame(b=1:B,
                       ltr0=results$t[,1],
                       ltr1=results$t[,2],
                       ltr2=results$t[,3])

    
    bootstrap_ltr_list[[rsid]][[endpoint]] <- temp
    
}
toc()

# initialize ltr df version where there are different rows for the estimate, and low/high confidence intervals
# from bootstrapping
temp <- ltr_df %>% filter(method=="coxph_bayes6")
temp$statistic <- "estimate"
ltr_df_ci <- temp

temp$statistic <- "loci"
ltr_df_ci <- rbind(ltr_df_ci, temp)

temp$statistic <- "hici"
ltr_df_ci <- rbind(ltr_df_ci, temp)
head(ltr_df_ci)

for (i in 1:nrow(ltr_df_ci)) {
    endpoint <- ltr_df_ci$endpoint[i]
    rsid <- ltr_df_ci$rsid[i]
    
    boot <- bootstrap_ltr_list[[rsid]][[endpoint]]
    if (ltr_df_ci$statistic[i] == "loci") {
        ltr_df_ci$ltr0[i] <- quantile(boot$ltr0, probs=c(0.025))
        ltr_df_ci$ltr1[i] <- quantile(boot$ltr1, probs=c(0.025))
        ltr_df_ci$ltr2[i] <- quantile(boot$ltr2, probs=c(0.025))
    } else if (ltr_df_ci$statistic[i] == "hici") {
        ltr_df_ci$ltr0[i] <- quantile(boot$ltr0, probs=c(0.975))
        ltr_df_ci$ltr1[i] <- quantile(boot$ltr1, probs=c(0.975))
        ltr_df_ci$ltr2[i] <- quantile(boot$ltr2, probs=c(0.975))
    }
}

# wrangle data frame into long format with estimate and confidence intervals as columns
ltr_df_ci <- ltr_df_ci %>%
    pivot_longer(names_to="ltr_level", values_to="value", c(ltr0, ltr1, ltr2)) %>%
    pivot_wider(names_from="statistic", values_from=c("value"))
head(ltr_df_ci)






