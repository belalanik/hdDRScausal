library(simsurv)
library(tableone)
library(autoCovariateSelection)
library(survival)
library(glmnet)
library(parallel)
library(doParallel)

setwd("~/GitHub/hdDRScausal")

# Analytic data
rm(list = ls())
load("Data/simdata.RData")

#################### Random proxies from different website ######################
# Proxies from different website such as from cms.gov, gov.bc.ca
proxies <- read.csv("Data/Proxy.csv", header = T)
head(proxies)

## 3-digit diagnostic codes from hospital database
dat.diag <- subset(proxies, dim == "diag")
dat.diag$code <- substr(dat.diag$code, start = 1, stop = 3)
dat.diag$code[dat.diag$code==""] <- NA
dat.diag <- na.omit(dat.diag)

## 3-digit procedure codes from hospital database
dat.proc <- subset(proxies, dim == "proc")
dat.proc$code <- substr(dat.proc$code, start = 1, stop = 3)
dat.proc$code[dat.proc$code==""] <- NA
dat.proc <- na.omit(dat.proc)

## 3-digit icd codes from physician claim database
dat.msp <- subset(proxies, dim == "msp")
dat.msp$code <- substr(dat.msp$code, start = 1, stop = 3)
dat.msp$code[dat.msp$code==""] <- NA
dat.msp <- na.omit(dat.msp)

## DINPIN from drug dispensation database
dat.din <- subset(proxies, dim == "din")
dat.din$code[dat.din$code==""] <- NA
dat.din <- na.omit(dat.din)

# Merge all dimensions
dat.proxy <- rbind(dat.diag, dat.proc, dat.msp, dat.din)

# Drop missing proxies 
dat.proxy <- na.omit(dat.proxy)
table(dat.proxy$dim, useNA = "always")

################### Converting proxies into empirically identified variables ################
id <- simdat$studyid

#### Generate candidate empirical covariates ####
step1 <- get_candidate_covariates(df = dat.proxy, domainVarname = "dim", 
                                  eventCodeVarname = "code", 
                                  patientIdVarname = "studyid", 
                                  patientIdVector = id, 
                                  n = 1000, 
                                  min_num_patients = 20)
out1 <- step1$covars_data

#### Assessing recurrence of codes ####
all.equal(id, step1$patientIds)

step2 <- get_recurrence_covariates(df = out1, 
                                   eventCodeVarname = "code", 
                                   patientIdVarname = "studyid",
                                   patientIdVector = id)
out2 <- step2$recurrence_data

#### Prioritizing covariates using Cox-LASSO ####
# Proxy variables
vars.proxy <- names(out2)[-1]

# Merge outcome data with proxies
dat.all <- merge(simdat[,c("studyid", "follow_up", "cvd")], out2, by = "studyid", all.x = T)

# Outcome model
formula.out <- as.formula(paste("Surv(follow_up, cvd) ~ ", paste(vars.proxy, collapse = " + ")))

# Model matrix for fitting LASSO
X <- model.matrix(formula.out, data = dat.all)[,-1]
Y <- as.matrix(data.frame(time = dat.all$follow_up, status = dat.all$cvd))

# Detect the number of cores
n_cores <- parallel::detectCores()

# Create a cluster of cores
cl <- makeCluster(n_cores - 1)

# Register the cluster for parallel processing
registerDoParallel(cl)

# 5-fold cross-validation for selecting hyperparameters
set.seed(123)
fit.lasso <- cv.glmnet(x = X, y = Y, nfolds = 5, parallel = T, alpha = 1, family = "cox")
stopCluster(cl)

plot(fit.lasso)

# Variable ranking based on Cox-LASSO
#lasso.coef <- coef(fit.lasso, s = fit.lasso$lambda.min)
lasso.coef <- coef(fit.lasso, s = exp(-6))
head(lasso.coef)
lasso.coef <- data.frame(as.matrix(lasso.coef))
lasso.coef <- data.frame(vars = rownames(lasso.coef), coef = lasso.coef)
colnames(lasso.coef) <- c("vars", "coef")
rownames(lasso.coef) <- NULL
lasso.coef$coef.abs <- abs(lasso.coef$coef)
lasso.coef <- lasso.coef[order(lasso.coef$coef.abs, decreasing = T),]
head(lasso.coef)

## Top 200 empirical covariates section based on Cox-LASSO
top200.lasso <- lasso.coef$vars[1:200]

# Dataset with all empirical covariates ranked by Cox-LASSO
dat.ec.lasso <- out2[, c("studyid", lasso.coef$vars)]

# Original simulated dataset with empirical covariates
dat.with.ec <- merge(simdat, dat.ec.lasso, by = "studyid")

#### Save ####
save(dat.proxy, dat.ec.lasso, lasso.coef, top200.lasso, dat.with.ec, file = "Data/simdata_with_proxy.RData")
