library(survminer)
test.ph <- cox.zph(fit_outcome)
test.ph
ggcoxzph(test.ph)
fit <- fitlist$KRA_PSY_DEMENTIA
sfit <- survfit(fit,
             newdata=data.frame(chr19_44908684_T_C=c(1),
                                FEMALE=0.5,
                                PC1=0,
                                PC2=0,
                                PC3=0,
                                PC4=0))
plot(sfit, xscale=1, xlab = "Years", ylab="Survival")
       
tab <- data.frame(sfit$time, sfit$cumhaz, sfit$surv, sfit$lower, sfit$upper)
fit <- coxph(Surv(futime, fustat) ~ age, data = ovarian) 
plot(survfit(fit, newdata=data.frame(age=60)),
     xscale=365.25, xlab = "Years", ylab="Survival") 

gehansurv=Surv(temp$KRA_PSY_DEMENTIA_AGE, temp$KRA_PSY_DEMENTIA)
plot(survfit(gehansurv ~ temp$chr19_44908684_T_C), col=c("black", "red", "blue"), fun="cloglog")



##NINAN SETIT
# Malli:
covs <- paste(c(paste0("PC", 1:10), names(data)[names(data) %like% "BATCH"]), collapse = " + ")

formula_prsvariable_25975 <- formula(paste0("Surv(", age, ", event = ", outcome, ") ~ ",
                         "prsvariable_25975 + strata(SEX_IMPUTED) + ", covs))

fit_prsvariable_25975 <- coxph(formula_prsvariable_25975, data = data)

# Adjustoidun selviytymiskäyräobjektin tallennus, tuon kunkin käyrän koordinaatit sisältävän plotdata-objektin tallennan sitten
# write.table():lla:
assign(paste0(outcome, "_plot_25975"), ggadjustedcurves(fit_prsvariable_25975, variable = "prsvariable_25975", fun = outcome,
                                              data = as.data.frame(data))
plotdata <- get(paste0(outcome, "_plot_25975"))
plotdata <- plotdata$data