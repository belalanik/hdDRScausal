library(survival)
library(Publish)
library(jtools)

setwd("~/GitHub/hdDRScausal")

################################# Data ################################
rm(list = ls())
load("Data/simdata.RData")
load("Data/simdata_with_proxy.RData")

# Keep only relevant dataset
rm(list = ls()[!ls() %in% c("dat.with.ec", "top200.lasso")])

#### Full cohort ####
dat.full <- dat.with.ec

# Offset of follow-up time for rate-based models 
dat.full$log.offset <- log(dat.full$follow_up)
dat.full$log.offset8 <- 0

# Full cohort with setting tb infection == "No"
dat.full.unexposed <- dat.full
dat.full.unexposed$tb.infection <- "No"
dat.full.unexposed$log.offset <- 0

#### Unexposed cohort ####
dat.unexposed <- subset(dat.full, tb.infection == "No")

#### Covariates ####
# Investigator-specified confounders 
vars.investigator <- c("age", "sex", "comorbidity")

# Empirical covariates
vars.empirical <- top200.lasso

# Investigator-specified + Empirical
vars <- c(vars.investigator, vars.empirical)


################################# Traditional ################################
# Formula
Formula.out <- as.formula(paste0("Surv(follow_up, cvd) ~ tb.infection + ", 
                                 paste(vars.investigator, collapse = "+")))

# Adjusted for investigator-specified confounders
fit.traditional <- coxph(Formula.out, data = dat.full)
publish(fit.traditional, print = F)$regressionTable[1:2,]


################################# hdPS ################################
# PS model specification
Formula.exp <- as.formula(paste0("I(tb.infection == 'Yes') ~ ", paste(vars, collapse = "+")))

# Fitting the PS model
fit.ps <- glm(Formula.exp, data = dat.full, family = binomial)
fit.ps$coefficients[is.na(fit.ps$coefficients)] <- 0
dat.full$ps <- predict(fit.ps, type = "response")
summary(dat.full$ps)

# Outcome analysis
dat.full$ps.decile <- as.factor(dplyr::ntile(dat.full$ps, 10))
fit.hdps <- coxph(Surv(follow_up, cvd) ~ tb.infection + ps.decile + age + sex + comorbidity, data = dat.full)
publish(fit.hdps, pvalue.method = "robust", confint.method = "robust", print = F)$regressionTable[1:2,]


################################# hdDRS-Full-Logistic ################################
# Outcome model specification
Formula.drs.full.logistic <- as.formula(paste0("cvd ~ tb.infection + ", paste(vars, collapse = "+")))

# Fitting the drs model 
fit.drs.full.logistic <- glm(Formula.drs.full.logistic, data = dat.full, family = binomial)
fit.drs.full.logistic$coefficients[is.na(fit.drs.full.logistic$coefficients)] <- 0
dat.full$drs.full.logistic <- predict(fit.drs.full.logistic, type = "response", newdata = dat.full.unexposed)
summary(dat.full$drs.full.logistic)

# Outcome analysis
dat.full$drs.full.logistic.decile <- as.factor(dplyr::ntile(dat.full$drs.full.logistic, 10))
fit.hddrs.full.logistic <- coxph(Surv(follow_up, cvd) ~ tb.infection + drs.full.logistic.decile + age + sex + 
                                  comorbidity, data = dat.full)
publish(fit.hddrs.full.logistic, pvalue.method = "robust", confint.method = "robust", 
        print = F)$regressionTable[1:2,]


################################# hdDRS-Full-Survival ################################
# Outcome model specification
Formula.drs.full.survival <- as.formula(paste0("Surv(follow_up, cvd) ~ tb.infection + ", 
                                               paste(vars, collapse = "+")))

# Fitting the drs model 
fit.drs.full.survival <- coxph(Formula.drs.full.survival, data = dat.full)
fit.drs.full.survival$coefficients[is.na(fit.drs.full.survival$coefficients)] <- 0

# Predicting drs 
dat.full$drs.full.survival <- predict(fit.drs.full.survival, type = "survival", newdata = dat.full.unexposed)
summary(dat.full$drs.full.survival)

# Outcome analysis 
dat.full$drs.full.survival.decile <- as.factor(dplyr::ntile(dat.full$drs.full.survival, 10))
fit.hddrs.full.survival <- coxph(Surv(follow_up, cvd) ~ tb.infection + drs.full.survival.decile + 
                                   age + sex + comorbidity, data = dat.full)
publish(fit.hddrs.full.survival, pvalue.method = "robust", confint.method = "robust", 
        print = F)$regressionTable[1:2,]


################################# hdDRS-Full-Hazard ################################
# Outcome model specification
Formula.drs.full.hazard <- as.formula(paste0("Surv(follow_up, cvd) ~ tb.infection + ", 
                                               paste(vars, collapse = "+")))

# Fitting the drs model 
fit.drs.full.hazard <- coxph(Formula.drs.full.hazard, data = dat.full)
fit.drs.full.hazard$coefficients[is.na(fit.drs.full.hazard$coefficients)] <- 0

# Predicting drs 
dat.full$drs.full.hazard <- predict(fit.drs.full.hazard, type = "risk", newdata = dat.full.unexposed)
summary(dat.full$drs.full.hazard)

# Outcome analysis with drs.full.hazard
dat.full$drs.full.hazard.decile <- as.factor(dplyr::ntile(dat.full$drs.full.hazard, 10))
fit.hddrs.full.hazard <- coxph(Surv(follow_up, cvd) ~ tb.infection + drs.full.hazard.decile + age + sex + 
                                 comorbidity, data = dat.full)
publish(fit.hddrs.full.hazard, pvalue.method = "robust", confint.method = "robust", 
        print = F)$regressionTable[1:2,]


################################# hdDRS-Full-Rate ################################
# Outcome model specification
Formula.drs.full.rate <- as.formula(paste0("cvd ~ tb.infection + offset(log.offset) + ", 
                                           paste(vars, collapse = "+")))

# Fitting the drs model 
fit.drs.full.rate <- glm(Formula.drs.full.rate, data = dat.full, family = poisson)
fit.drs.full.rate$coefficients[is.na(fit.drs.full.rate$coefficients)] <- 0

# Predicting drs 
dat.full$drs.full.rate <- predict(fit.drs.full.rate, type = "response", newdata = dat.full.unexposed)
summary(dat.full$drs.full.rate)

# Outcome analysis
dat.full$drs.full.rate.decile <- as.factor(dplyr::ntile(dat.full$drs.full.rate, 10))
fit.hddrs.full.rate <- coxph(Surv(follow_up, cvd) ~ tb.infection + drs.full.rate.decile + age + 
                               sex + comorbidity, data = dat.full)
publish(fit.hddrs.full.rate, pvalue.method = "robust", confint.method = "robust", 
        print = F)$regressionTable[1:2,]


################################# hdDRS-Unexposed-Logistic ################################
# Outcome model specification
Formula.drs.unexposed.logistic <- as.formula(paste0("cvd ~ ", paste(vars, collapse = "+")))

# Fitting the drs model 
fit.drs.unexposed.logistic <- glm(Formula.drs.unexposed.logistic, data = dat.unexposed, family = binomial)
fit.drs.unexposed.logistic$coefficients[is.na(fit.drs.unexposed.logistic$coefficients)] <- 0

# Predicting drs 
dat.full$drs.unexposed.logistic <- predict(fit.drs.unexposed.logistic, type = "response", newdata = dat.full)
summary(dat.full$drs.unexposed.logistic)

# Outcome analysis
dat.full$drs.unexposed.logistic.decile <- as.factor(dplyr::ntile(dat.full$drs.unexposed.logistic, 10))
fit.hddrs.unexposed.logistic <- coxph(Surv(follow_up, cvd) ~ tb.infection + drs.unexposed.logistic.decile + 
                                        age + sex + comorbidity, data = dat.full)
publish(fit.hddrs.unexposed.logistic, pvalue.method = "robust", confint.method = "robust", 
        print = F)$regressionTable[1:2,]


################################# hdDRS-Unexposed-Survival ################################
# Outcome model specification
Formula.drs.unexposed.survival <- as.formula(paste0("Surv(follow_up, cvd) ~ ", paste(vars, collapse = "+")))

# Fitting the drs model 
fit.drs.unexposed.survival <- coxph(Formula.drs.unexposed.survival, data = dat.unexposed)
fit.drs.unexposed.survival$coefficients[is.na(fit.drs.unexposed.survival$coefficients)] <- 0

# Predicting drs 
dat.full$drs.unexposed.survival <- predict(fit.drs.full.survival, type = "survival", newdata = dat.full)
summary(dat.full$drs.unexposed.survival)

# Outcome analysis with drs.unexposed.survival
dat.full$drs.unexposed.survival.decile <- as.factor(dplyr::ntile(dat.full$drs.unexposed.survival, 10))
fit.hddrs.unexposed.survival <- coxph(Surv(follow_up, cvd) ~ tb.infection + drs.unexposed.survival.decile +
                                        age + sex + comorbidity, data = dat.full)
publish(fit.hddrs.unexposed.survival, pvalue.method = "robust", confint.method = "robust", 
        print = F)$regressionTable[1:2,]


################################# hdDRS-Unexposed-Hazard ################################
# Outcome model specification
Formula.drs.unexposed.hazard <- as.formula(paste0("Surv(follow_up, cvd) ~ ", paste(vars, collapse = "+")))

# Fitting the drs model 
fit.drs.unexposed.hazard <- coxph(Formula.drs.unexposed.hazard, data = dat.unexposed)
fit.drs.unexposed.hazard$coefficients[is.na(fit.drs.unexposed.hazard$coefficients)] <- 0

# Predicting drs 
dat.full$drs.unexposed.hazard <- predict(fit.drs.unexposed.hazard, type = "risk", newdata = dat.full)
summary(dat.full$drs.unexposed.hazard)

# Outcome analysis with drs.unexposed.hazard
dat.full$drs.unexposed.hazard.decile <- as.factor(dplyr::ntile(dat.full$drs.unexposed.hazard, 10))
fit.hddrs.unexposed.hazard <- coxph(Surv(follow_up, cvd) ~ tb.infection + drs.unexposed.hazard.decile + 
                                      age + sex + comorbidity, data = dat.full)
publish(fit.hddrs.unexposed.hazard, pvalue.method = "robust", confint.method = "robust", 
        print = F)$regressionTable[1:2,]


################################# hdDRS-Unexposed-Rate ################################
# Outcome model specification
Formula.drs.unexposed.rate <- as.formula(paste0("cvd ~ offset(log.offset8) + ", paste(vars, collapse = "+")))

# Fitting the drs model 
fit.drs.unexposed.rate <- glm(Formula.drs.unexposed.rate, data = dat.unexposed, family = poisson)
fit.drs.unexposed.rate$coefficients[is.na(fit.drs.unexposed.rate$coefficients)] <- 0

# Predicting drs 
dat.full$drs.unexposed.rate <- predict(fit.drs.unexposed.rate, type = "response", newdata = dat.full)
summary(dat.full$drs.unexposed.rate)

# Outcome analysis
dat.full$drs.unexposed.rate.decile <- as.factor(dplyr::ntile(dat.full$drs.unexposed.rate, 10))
fit.hddrs.unexposed.rate <- coxph(Surv(follow_up, cvd) ~ tb.infection + drs.unexposed.rate.decile + age + 
                                    sex + comorbidity, data = dat.full)
publish(fit.hddrs.unexposed.rate, pvalue.method = "robust", confint.method = "robust", 
        print = F)$regressionTable[1:2,]

################################# Plot all results ################################
plot_summs(fit.traditional, fit.hdps, fit.hddrs.full.logistic, fit.hddrs.full.survival, 
           fit.drs.full.hazard, fit.drs.full.rate, fit.hddrs.unexposed.logistic, 
           fit.hddrs.unexposed.survival, fit.hddrs.unexposed.hazard,fit.hddrs.unexposed.rate,
           scale = F, robust = TRUE,
           coefs = c("tb.infection" = "TB infection"),
           exp = TRUE,
           legend.title = "Method",
           model.names = c("Traditional", "hdPS", "hdDRS-Full-Logistic", "hdDRS-Full-Survival", 
                           "hdDRS-Full-Hazard", "hdDRS-Full-Rate", "hdDRS-Unexposed-Logistic", 
                           "hdDRS-Unexposed-Survival", "hdDRS-Unexposed-Hazard","hdDRS-Unexposed-Rate"))

