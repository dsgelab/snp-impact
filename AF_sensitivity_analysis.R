library(AF)
compare_af <- function(data=fg_df, variant, cause, cause_age) {
    data$snp_12 <- ifelse(data[,variant]==0, 0, 1)
    data$snp_1 <- ifelse(data[,variant]==0, 0,
                        ifelse(data[,variant]==1, 1, NA))
    
    covariates_formula <- paste("Surv(",  cause_age, ", ", cause, ") ~ ", "snp_12", "+ FEMALE + ", paste(pcs, collapse=" + "), sep="")
    form <- formula(covariates_formula)
    AFest <- AF.ch(form, data = data, exposure = "snp_12", times = seq(10,90,10))
    print(AFest)
    plot(AFest, CI = TRUE, xlab = "time (years)")
     
    fit <- coxph(form, data=data)
    HR <- summary(fit)$coefficients[1,2]
    P <- sum(data[,"snp_12"]==1)/nrow(data)
    print(paste("COXPH AF:", (P*(HR-1))/(1+(P*(HR-1)))))
    
    
    covariates_formula <- paste("Surv(",  cause_age, ", ", cause, ") ~ ", "snp_1", "+ FEMALE + ", paste(pcs, collapse=" + "), sep="")
    form <- formula(covariates_formula)
    AFest <- AF.ch(form, data = data[!is.na(data$snp_1),], exposure = "snp_1", times = seq(10,90,10))
    print(AFest)
    plot(AFest, CI = TRUE, xlab = "time (years)")
}

b2v <- function(x) {
    return(paste("chr", gsub(":", "_", x),sep=""))
}

compare_af(data=fg_df, variant="chr19_44908684_T_C",
           cause="Alzheimers_disease_and_other_dementias",
           cause_age="F5_DEMENTIA_AGE")

compare_af(data=fg_df, variant=b2v("2:43839108:G:C"),
           cause="Gallbladder_and_biliary_diseases",
           cause_age="Gallbladder_and_biliary_diseases_AGE")

compare_af(data=fg_df, variant=b2v("20:44189982:G:T"),
           cause="Diabetes_mellitus_type_2",
           cause_age="T2D_WIDE_AGE")