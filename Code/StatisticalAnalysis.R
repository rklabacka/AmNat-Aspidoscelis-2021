# Load libraries
library(tidyverse)
library(rstatix)
library(broom)
library(factoextra)
library(dplyr)
library(nlme) # for building of linear model using Frequentist approach
library(piecewiseSEM) # for calculating r^2
library(boot) # for bootstrapping
# - Use Bayesian approach to mean-correct variance values and compare
library(MCMCglmm) # For linear models using Bayesian approach
library(wolakR) #<-- handy functions to work with MCMCglmm models
### --- Can use `remotes` package to install from GitHub:
### --- remotes::install_github("matthewwolak/wolakR")


## All phylogenetic network analyses performed in the Julia package PhyloNetworks (see AspidoscelisReticAnalysis.jl)

# Read in physiology data
# This is assuming the code was ran from the directory where it is stored
dat_phys_indiv <- read.csv("../SampleInformation/PhysiologyData_2019_Individuals.csv")

# - Prepare data
## -- Scale and center SVL
dat_phys_indiv$scSVL <- scale(dat_phys_indiv$SVL)  #<-- work with this from now on
## -- Convert Sexual Mode variable to factor
dat_phys_indiv$SexualMode <- as.factor(dat_phys_indiv$SexualMode)
# Get summary statistics:
dat_phys_indiv %>% group_by(Sex) %>% group_by(Species) %>% summarize(m=mean(CII_State3), sd=sd(CII_State3), n=n(), ci=sd / sqrt(n))

# -----------------------------------------------------------------------------------------------------
#    -----------------------------------------------------------------------------------------------------
# ----------------------------------     Section I: LME APPROACH     ----------------------------------
#    ----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

# Examine effect of sexual mode on response variables and perform model selection for one or two residual variance parameters

## **********************
## ****  Endurance  *****
## **********************
# NOTE: Comparisons between models a and b (within both Model 1 and Model 2) is to check whether log-transformed data better fits the model.
## -- Model 1: No differences in variance between sexual modes
endur_lme_1a <- lme(Endurance ~ SexualMode + scSVL, data = dat_phys_indiv, random = ~ 1 | Species)
endur_lme_1b <- lme(Log.Endurance ~ SexualMode + scSVL, data = dat_phys_indiv, random = ~ 1 | Species)
## -- Model 2: Differences in variance between sexual modes
endur_lme_2a <- lme(Endurance ~ SexualMode + scSVL, data = dat_phys_indiv, random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
endur_lme_2b <- lme(Log.Endurance ~ SexualMode + scSVL, data = dat_phys_indiv, random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
# Compare the two models: Hypothesis test for different (residual) variances
anova(endur_lme_1a, endur_lme_2a) # model 2a preferred, AIC = 134.9557, p = <.0001
anova(endur_lme_1b, endur_lme_2b) # model 2b preferred, AIC = -7.647765, p = 0.0227
# Plot residuals of models (notice model 2b shows less heteroscedasticity in plotted residuals)
plot(endur_lme_2a)
plot(endur_lme_2b)
## -- See details of selected model
summary(endur_lme_2b)

bootFunc=function(bootData,repeats){
    tryCatch({
    # Fit model
    boot_lme <- lme(Log.Endurance ~ SexualMode + scSVL, data = bootData[repeats,], random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
    # Extract standard deviation ratios
    boot_sdrs <- coef(boot_lme$modelStruct$varStruct, unconstrained = FALSE)
    # Multiply standard deviation ratios by overall residual standard deviation
    boot_sds <- c(1, boot_sdrs) * boot_lme$sigma
    return(boot_sds)
    },
    error = function(err) {return(NA)}
    )
}

sds <- boot(dat_phys_indiv, bootFunc, R=1000)
#       original      bias    std. error
# t1* 0.18312367 -0.01953424  0.03015752
# t2* 0.09026297 -0.01074882  0.01597630


## ***********************
## ****  CI State 3  *****
## ***********************
## -- Model 1: No differences in variance between sexual modes
CIS3_lme_1 <- lme(CI_State3 ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species)
## -- Model 2: Differences in variance between sexual modes
CIS3_lme_2 <- lme(CI_State3 ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
# Compare the two models: Hypothesis test for different (residual) variances
anova(CIS3_lme_1, CIS3_lme_2) # model 1 preferred, p = 0.8405
## -- See details of selected model
summary(CIS3_lme_1)

bootFunc=function(bootData,repeats){
    tryCatch({
    # Fit model
    boot_lme <- lme(CI_State3 ~ SexualMode, data = bootData[repeats,], random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
    # Extract standard deviation ratios
    boot_sdrs <- coef(boot_lme$modelStruct$varStruct, unconstrained = FALSE)
    # Multiply standard deviation ratios by overall residual standard deviation
    boot_sds <- c(1, boot_sdrs) * boot_lme$sigma
    return(boot_sds)
    },
    error = function(err) {return(NA)}
    )
}

sds <- boot(dat_phys_indiv, bootFunc, R=1000)
#     original     bias    std. error
# t1* 4.772031 -0.2755333   0.5613797
# t2* 4.510147 -0.3632354   0.8413166


## ***********************
## ****  CI State 4  *****
## ***********************
## -- Model 1: No differences in variance between sexual modes
CIS4_lme_1 <- lme(CI_State4 ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species)
## -- Model 2: Differences in variance between sexual modes
CIS4_lme_2 <- gls(CI_State4 ~ SexualMode, data = dat_phys_indiv, weights = varIdent(form = ~ 1 | SexualMode))
# Compare the two models: Hypothesis test for different (residual) variances
anova(CIS4_lme_1, CIS4_lme_2) # model 2 preferred, p = 7e-04
## -- See details of selected model
summary(CIS4_lme_2)

bootFunc=function(bootData,repeats){
    tryCatch({
    # Fit model
    boot_lme <- lme(CI_State4 ~ SexualMode, data = bootData[repeats,], random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
    # Extract standard deviation ratios
    boot_sdrs <- coef(boot_lme$modelStruct$varStruct, unconstrained = FALSE)
    # Multiply standard deviation ratios by overall residual standard deviation
    boot_sds <- c(1, boot_sdrs) * boot_lme$sigma
    return(boot_sds)
    },
    error = function(err) {return(NA)}
    )
}

sds <- boot(dat_phys_indiv, bootFunc, R=1000)
#      original     bias    std. error
# t1* 2.2594527 -0.1482513   0.3608133
# t2* 0.7460734 -0.0938581   0.1872636

## ***********************
## ****    CI RCR    *****
## ***********************
## -- Model 1: No differences in variance between sexual modes
CIRCR_lme_1 <- lme(CI_RCR ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species)
## -- Model 2: Differences in variance between sexual modes
CIRCR_lme_2 <- lme(CI_RCR ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
# Compare the two models: Hypothesis test for different (residual) variances
anova(CIRCR_lme_1, CIRCR_lme_2) # model 1 preferred, p = 0.37
## -- See details of selected model
summary(CIRCR_lme_1)

bootFunc=function(bootData,repeats){
    tryCatch({
    # Fit model
    boot_lme <- lme(CI_RCR ~ SexualMode, data = bootData[repeats,], random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
    # Extract standard deviation ratios
    boot_sdrs <- coef(boot_lme$modelStruct$varStruct, unconstrained = FALSE)
    # Multiply standard deviation ratios by overall residual standard deviation
    boot_sds <- c(1, boot_sdrs) * boot_lme$sigma
    return(boot_sds)
    },
    error = function(err) {return(NA)}
    )
}

sds <- boot(dat_phys_indiv, bootFunc, R=1000)
#      original      bias    std. error
# t1* 0.6864893 -0.07373868   0.1709778
# t2* 0.8782012 -0.06456399   0.209786

## ***********************
## ****  CII State 3 *****
## ***********************
## -- Model 1: No differences in variance between sexual modes
CIIS3_lme_1 <- lme(CII_State3 ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species)
## -- Model 2: Differences in variance between sexual modes
CIIS3_lme_2 <- lme(CII_State3 ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
# Compare the two models: Hypothesis test for different (residual) variances
anova(CIIS3_lme_1, CIIS3_lme_2) # model 2 preferred, p = 0.047
## -- See details of selected model
summary(CIIS3_lme_2)

bootFunc=function(bootData,repeats){
    tryCatch({
    # Fit model
    boot_lme <- lme(CII_State3 ~ SexualMode, data = bootData[repeats,], random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
    # Extract standard deviation ratios
    boot_sdrs <- coef(boot_lme$modelStruct$varStruct, unconstrained = FALSE)
    # Multiply standard deviation ratios by overall residual standard deviation
    boot_sds <- c(1, boot_sdrs) * boot_lme$sigma
    return(boot_sds)
    },
    error = function(err) {return(NA)}
    )
}

sds <- boot(dat_phys_indiv, bootFunc, R=1000)
#     original     bias    std. error
# t1* 4.531504 -0.2930536   0.6174195
# t2* 2.490282 -0.3456246   0.4801312

## ***********************
## ****  CII State 4 *****
## ***********************
## -- Model 1: No differences in variance between sexual modes
CIIS4_lme_1 <- lme(CII_State4 ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species)
## -- Model 2: Differences in variance between sexual modes
CIIS4_lme_2 <- lme(CII_State4 ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
# Compare the two models: Hypothesis test for different (residual) variances
anova(CIIS4_lme_1, CIIS4_lme_2) # model 2 preferred, p = 0.004
## -- See details of selected model
summary(CIIS4_lme_2)

bootFunc=function(bootData,repeats){
    tryCatch({
    # Fit model
    boot_lme <- lme(CII_State4 ~ SexualMode, data = bootData[repeats,], random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
    # Extract standard deviation ratios
    boot_sdrs <- coef(boot_lme$modelStruct$varStruct, unconstrained = FALSE)
    # Multiply standard deviation ratios by overall residual standard deviation
    boot_sds <- c(1, boot_sdrs) * boot_lme$sigma
    return(boot_sds)
    },
    error = function(err) {return(NA)}
    )
}

sds <- boot(dat_phys_indiv, bootFunc, R=1000)
#     original     bias    std. error
# t1* 3.290843 -0.2512533   0.5134656
# t2* 1.321751 -0.1666855   0.3584912

## ***********************
## ****    CII RCR   *****
## ***********************
## -- Model 1: No differences in variance between sexual modes
CIIRCR_lme_1 <- lme(CII_RCR ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species)
## -- Model 2: Differences in variance between sexual modes
CIIRCR_lme_2 <- lme(CII_RCR ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
# Compare the two models: Hypothesis test for different (residual) variances
anova(CIIRCR_lme_1, CIIRCR_lme_2) # model 1 preferred, p = 0.088
## -- See details of selected model
summary(CIIRCR_lme_1)

bootFunc=function(bootData,repeats){
    tryCatch({
    # Fit model
    boot_lme <- lme(CII_RCR ~ SexualMode, data = bootData[repeats,], random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
    # Extract standard deviation ratios
    boot_sdrs <- coef(boot_lme$modelStruct$varStruct, unconstrained = FALSE)
    # Multiply standard deviation ratios by overall residual standard deviation
    boot_sds <- c(1, boot_sdrs) * boot_lme$sigma
    return(boot_sds)
    },
    error = function(err) {return(NA)}
    )
}

sds <- boot(dat_phys_indiv, bootFunc, R=1000)
#      original      bias    std. error
# t1* 0.3590600 -0.05336864  0.08380936
# t2* 0.2141631 -0.01807205  0.03937686

## ***********************
## ****  Endur~Mito  *****
## ***********************
# Examine the relationship between endurance and respiration
all_EndurCIS3_lme<-  lme(Log.Endurance ~ CI_State3, random = ~ 1 | Species, data=dat_phys_indiv)
rsquared(all_EndurCIS3_lme) # Marginal

all_EndurCIS4_lme<-  lme(Log.Endurance ~ CI_State4, random = ~ 1 | Species, data=dat_phys_indiv)
rsquared(all_EndurCIS4_lme) # Marginal

all_EndurCIRCR_lme<-  lme(Log.Endurance ~ CI_RCR, random = ~ 1 | Species, data=dat_phys_indiv)
rsquared(all_EndurCIRCR_lme) # Marginal

all_EndurCIIS3_lme<-  lme(Log.Endurance ~ CII_State3, random = ~ 1 | Species, data=dat_phys_indiv)
rsquared(all_EndurCIIS3_lme) # Marginal

all_EndurCIIS4_lme<-  lme(Log.Endurance ~ CII_State4, random = ~ 1 | Species, data=dat_phys_indiv)
rsquared(all_EndurCIIS4_lme) # Marginal

all_EndurCIIRCR_lme<-  lme(Log.Endurance ~ CII_RCR, random = ~ 1 | Species, data=dat_phys_indiv)
rsquared(all_EndurCIIRCR_lme) # Marginal



# -----------------------------------------------------------------------------------------------------
#    -----------------------------------------------------------------------------------------------------
# ----------------------------------  Section II: MCMCglmm APPROACH   ---------------------------------
#    ----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
# Use Bayesian approach to examine differences in variance between mean-corrected variance (a.k.a., coefficient of variation)
## -- We did this for selected models which two residual variance parameters (see above)

# First set our priors
## Use scaled central F-distribution for Species random effects prior distribution
## Use non-informative improper prior for residual variances (marginal posterior distribution of variance component is equivalent to the REML estimate)
pr1 <- list(R = list(V = 1e-12, nu = -2),
  G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = c(0), alpha.V = diag(1)*1000)))

# multivariate version of the non-informative improper prior for 2 residual variances
pr2 <- list(R = list(V = diag(2)*1e-12, nu = -1),#<-- only thing differing from pr1
  G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = c(0), alpha.V = diag(1)*1000)))

# Setup MCMC chain parameters
nsample <- 1000
BURN <- 3000; THIN <- 50; NITT <- BURN + nsample * THIN

# Set up plotting environment
par(mfrow=c(2,1))

## ***********************
## ****   Endurance  *****
## ***********************
# Create model
endur_b <- MCMCglmm(Log.Endurance ~ SexualMode + scSVL,
  random = ~ Species,
  rcov = ~ idh(SexualMode):units,
  data = dat_phys_indiv,
  prior = pr2,
  nitt = NITT, burnin = BURN, thin = THIN)

# Plot standard deviations for each repro mode
endur_b_sdPost <- sqrt(endur_b$VCV)
# Get mean values
mean(endur_b_sdPost[, 2])
mean(endur_b_sdPost[, 3])
# Get confidence intervals
endur_b_sd_hpd <- HPDinterval(endur_b_sdPost)
# Plot posterior distributions
par(mfrow=c(2,1))
endur_asex_sd <- postPlot(endur_b_sdPost[, 2], xlim = c(0, 0.40))
# add sd determined from frequentist approach
abline(v=0.09, col="red", lwd=3)
endur_sex_sd <- postPlot(endur_b_sdPost[, 3], xlim = c(0, 0.40))
abline(v=0.18, col="red", lwd=3)

# Plot differences between coefficient of variation
endur_cvdiffPost <- ((sqrt(endur_b$VCV[, 3]) / (endur_b$Sol[, 1] + endur_b$Sol[, 2]) - sqrt(endur_b$VCV[, 2]) / endur_b$Sol[, 1]))
endur_hpd <- HPDinterval(endur_cvdiffPost)
endur_b_postPlot <- postPlot(endur_cvdiffPost, xlim=c(-0.5,0.5),plotHist = FALSE)
endur_b_post <- mean(endur_cvdiffPost > 0) # posterior probability that asexual coefficient of variance < sexual
abline(v=0, col="red", lwd=3)

## ***********************
## ****  CI State 3  *****
## ***********************
# Create model
CIS3_b <- MCMCglmm(CI_State3 ~ SexualMode ,
  random = ~ Species,
  rcov = ~ idh(SexualMode):units,
  data = dat_phys_indiv,
  prior = pr2,
  nitt = NITT, burnin = BURN, thin = THIN)

# Plot standard deviations for each repro mode
CIS3_b_sdPost <- sqrt(CIS3_b$VCV)
# Get mean values
mean(CIS3_b_sdPost[, 2]) #asex
mean(CIS3_b_sdPost[, 3]) #sex
# Get confidence intervals
CIS3_b_sd_hpd <- HPDinterval(CIS3_b_sdPost)
# Plot posterior distributions
par(mfrow=c(2,1))
CIS3_asex_sd <- postPlot(CIS3_b_sdPost[, 2], xlim=c(2,12))
# add sd determined from frequentist approach
abline(v=4.51, col="red", lwd=3)
CIS3_sex_sd <- postPlot(CIS3_b_sdPost[, 3], xlim=c(2,12))
# add sd determined from frequentist approach
abline(v=4.77, col="red", lwd=3)

CIS3_cvdiffPost <- ((sqrt(CIS3_b$VCV[, 3]) / (CIS3_b$Sol[, 1] + CIS3_b$Sol[, 2]) - sqrt(CIS3_b$VCV[, 2]) / CIS3_b$Sol[, 1]))
CIS3_hpd <- HPDinterval(CIS3_cvdiffPost)
postPlot(CIS3_cvdiffPost, xlim = c(-0.5, 0.5), plotHist = FALSE)
CIS3_b_post <- mean(CIS3_cvdiffPost > 0) # posterior probability that asexual coefficient of variance < sexual
abline(v=0, col="red", lwd=3)

## ***********************
## ****  CI State 4  *****
## ***********************
# Create model
CIS4_b <- MCMCglmm(CI_State4 ~ SexualMode ,
  random = ~ Species,
  rcov = ~ idh(SexualMode):units,
  data = dat_phys_indiv,
  prior = pr2,
  nitt = NITT, burnin = BURN, thin = THIN)

# Plot standard deviations for each repro mode
CIS4_b_sdPost <- sqrt(CIS4_b$VCV)
# Get mean values
mean(CIS4_b_sdPost[, 2]) #asex
mean(CIS4_b_sdPost[, 3]) #sex
# Get confidence intervals
CIS4_b_sd_hpd <- HPDinterval(CIS4_b_sdPost)
# Plot posterior distributions
CIS4_b_sdPost <- sqrt(CIS4_b$VCV)
par(mfrow=c(2,1))
CIS4_asex_sd <- postPlot(CIS4_b_sdPost[, 2], xlim=c(0,5))
# add sd determined from frequentist approach
abline(v=0.75, col="red", lwd=3)
CIS4_sex_sd <- postPlot(CIS4_b_sdPost[, 3], xlim=c(0,5))
# add sd determined from frequentist approach
abline(v=2.26, col="red", lwd=3)

CIS4_cvdiffPost <- ((sqrt(CIS4_b$VCV[, 3]) / (CIS4_b$Sol[, 1] + CIS4_b$Sol[, 2]) - sqrt(CIS4_b$VCV[, 2]) / CIS4_b$Sol[, 1]))
CIS4_hpd <- HPDinterval(CIS4_cvdiffPost)
postPlot(CIS4_cvdiffPost, xlim = c(-0.5, 0.5), plotHist = FALSE)
CIS4_b_post <- mean(CIS4_cvdiffPost > 0) # posterior probability that asexual coefficient of variance < sexual
abline(v=0, col="red", lwd=3)

## ***********************
## ****    CI RCR    *****
## ***********************
# Create model
CIRCR_b <- MCMCglmm(CI_RCR ~ SexualMode ,
  random = ~ Species,
  rcov = ~ idh(SexualMode):units,
  data = dat_phys_indiv,
  prior = pr2,
  nitt = NITT, burnin = BURN, thin = THIN)

# Plot standard deviations for each repro mode
CIRCR_b_sdPost <- sqrt(CIRCR_b$VCV)
# Get mean values
mean(CIRCR_b_sdPost[, 2]) #asex
mean(CIRCR_b_sdPost[, 3]) #sex
# Get confidence intervals
CIRCR_b_sd_hpd <- HPDinterval(CIRCR_b_sdPost)
# Plot posterior distributions
CIRCR_b_sdPost <- sqrt(CIRCR_b$VCV)
par(mfrow=c(2,1))
CIRCR_asex_sd <- postPlot(CIRCR_b_sdPost[, 2], xlim=c(0,2.5))
# add sd determined from frequentist approach
abline(v=0.88, col="red", lwd=3)
CIRCR_sex_sd <- postPlot(CIRCR_b_sdPost[, 3], xlim=c(0,2.5))
# add sd determined from frequentist approach
abline(v=0.69, col="red", lwd=3)

CIRCR_cvdiffPost <- ((sqrt(CIRCR_b$VCV[, 3]) / (CIRCR_b$Sol[, 1] + CIRCR_b$Sol[, 2]) - sqrt(CIRCR_b$VCV[, 2]) / CIRCR_b$Sol[, 1]))
CIRCR_hpd <- HPDinterval(CIRCR_cvdiffPost)
postPlot(CIRCR_cvdiffPost, xlim = c(-0.5, 0.5), plotHist = FALSE)
CIRCR_b_post <- mean(CIRCR_cvdiffPost > 0) # posterior probability that asexual coefficient of variance < sexual
abline(v=0, col="red", lwd=3)

## ***********************
## ****  CII State 3 *****
## ***********************
# Create model
CIIS3_b <- MCMCglmm(CII_State3 ~ SexualMode ,
  random = ~ Species,
  rcov = ~ idh(SexualMode):units,
  data = dat_phys_indiv,
  prior = pr2,
  nitt = NITT, burnin = BURN, thin = THIN)

# Plot standard deviations for each repro mode
CIIS3_b_sdPost <- sqrt(CIIS3_b$VCV)
# Get mean values
mean(CIIS3_b_sdPost[, 2]) #asex
mean(CIIS3_b_sdPost[, 3]) #sex
# Get confidence intervals
CIIS3_b_sd_hpd <- HPDinterval(CIIS3_b_sdPost)
# Plot posterior distributions
CIIS3_b_sdPost <- sqrt(CIIS3_b$VCV)
par(mfrow=c(2,1))
CIIS3_asex_sd <- postPlot(CIIS3_b_sdPost[, 2], xlim=c(0,10))
# add sd determined from frequentist approach
abline(v=2.49, col="red", lwd=3)
CIIS3_sex_sd <- postPlot(CIIS3_b_sdPost[, 3], xlim=c(0,10))
# add sd determined from frequentist approach
abline(v=4.53, col="red", lwd=3)

CIIS3_cvdiffPost <- ((sqrt(CIIS3_b$VCV[, 3]) / (CIIS3_b$Sol[, 1] + CIIS3_b$Sol[, 2]) - sqrt(CIIS3_b$VCV[, 2]) / CIIS3_b$Sol[, 1]))
CIIS3_hpd <- HPDinterval(CIIS3_cvdiffPost)
postPlot(CIIS3_cvdiffPost, xlim = c(-0.5, 0.5), plotHist = FALSE)
CIIS3_b_post <- mean(CIIS3_cvdiffPost > 0) # posterior probability that asexual coefficient of variance < sexual
abline(v=0, col="red", lwd=3)

## ***********************
## ****  CII State 4 *****
## ***********************
# Create model
CIIS4_b <- MCMCglmm(CII_State4 ~ SexualMode ,
  random = ~ Species,
  rcov = ~ idh(SexualMode):units,
  data = dat_phys_indiv,
  prior = pr2,
  nitt = NITT, burnin = BURN, thin = THIN)


# Plot standard deviations for each repro mode
CIIS4_b_sdPost <- sqrt(CIIS4_b$VCV)
# Get mean values
mean(CIIS4_b_sdPost[, 2]) #asex
mean(CIIS4_b_sdPost[, 3]) #sex
# Get confidence intervals
CIIS4_b_sd_hpd <- HPDinterval(CIIS4_b_sdPost)
# Plot posterior distributions
CIIS4_b_sdPost <- sqrt(CIIS4_b$VCV)
par(mfrow=c(2,1))
CIIS4_asex_sd <- postPlot(CIIS4_b_sdPost[, 2], xlim=c(0,7))
# add sd determined from frequentist approach
abline(v=1.32, col="red", lwd=3)
CIIS4_sex_sd <- postPlot(CIIS4_b_sdPost[, 3], xlim=c(0,7))
# add sd determined from frequentist approach
abline(v=3.29, col="red", lwd=3)

CIIS4_cvdiffPost <- ((sqrt(CIIS4_b$VCV[, 3]) / (CIIS4_b$Sol[, 1] + CIIS4_b$Sol[, 2]) - sqrt(CIIS4_b$VCV[, 2]) / CIIS4_b$Sol[, 1]))
CIIS4_hpd <- HPDinterval(CIIS4_cvdiffPost)
postPlot(CIIS4_cvdiffPost, xlim = c(-0.5, 0.5), plotHist = FALSE)
CIIS4_b_post <- mean(CIIS4_cvdiffPost > 0) # posterior probability that asexual coefficient of variance < sexual
abline(v=0, col="red", lwd=3)

## ***********************
## ****    CII RCR   *****
## ***********************
# Create model
CIIRCR_b <- MCMCglmm(CII_RCR ~ SexualMode ,
  random = ~ Species,
  rcov = ~ idh(SexualMode):units,
  data = dat_phys_indiv,
  prior = pr2,
  nitt = NITT, burnin = BURN, thin = THIN)

# Plot standard deviations for each repro mode
CIIRCR_b_sdPost <- sqrt(CIIRCR_b$VCV)
# Get mean values
mean(CIIRCR_b_sdPost[, 2]) #asex
mean(CIIRCR_b_sdPost[, 3]) #sex
# Get confidence intervals
CIIRCR_b_sd_hpd <- HPDinterval(CIIRCR_b_sdPost)
# Plot posterior distributions
CIIRCR_b_sdPost <- sqrt(CIIRCR_b$VCV)
# Get mean values
mean(CIIRCR_b_sdPost[, 2])
mean(CIIRCR_b_sdPost[, 3])
# Get confidence intervals
CIIRCR_sdhpd <- HPDinterval(CIIRCR_b_sdPost)
# Plot posterior distributions
par(mfrow=c(2,1))
CIIRCR_asex_sd <- postPlot(CIIRCR_b_sdPost[, 2], xlim=c(0,0.7))
# add sd determined from frequentist approach
abline(v=0.21, col="red", lwd=3)
CIIRCR_sex_sd <- postPlot(CIIRCR_b_sdPost[, 3], xlim=c(0,0.7))
# add sd determined from frequentist approach
abline(v=0.36, col="red", lwd=3)

CIIRCR_cvdiffPost <- ((sqrt(CIIRCR_b$VCV[, 3]) / (CIIRCR_b$Sol[, 1] + CIIRCR_b$Sol[, 2]) - sqrt(CIIRCR_b$VCV[, 2]) / CIIRCR_b$Sol[, 1]))
CIIRCR_hpd <- HPDinterval(CIIRCR_cvdiffPost)
postPlot(CIIRCR_cvdiffPost, xlim = c(-0.5, 0.5), plotHist = FALSE)
CIIRCR_b_post <- mean(CIIRCR_cvdiffPost > 0) # posterior probability that asexual coefficient of variance < sexual




# -----------------------------------------------------------------------------------------------------
#    -----------------------------------------------------------------------------------------------------
# -----------------------------------  Section III: Within-group LM  ----------------------------------
#    ----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

# Subset data
dat_phys_tess <- subset(dat_phys_indiv, Species == 'tesselatus' | Species == 'marmoratus' | Species == 'septemvittatus')
dat_phys_neom <- subset(dat_phys_indiv, Species == 'neomexicanus' | Species == 'marmoratus' | Species == 'inornatus')
dat_phys_marm <- subset(dat_phys_indiv, Species == 'neomexicanus' | Species == 'marmoratus' | Species == 'tesselatus')

dat_phys_tess$Species <- factor(dat_phys_tess$Species, ordered=FALSE)
dat_phys_neom$Species <- factor(dat_phys_neom$Species, ordered=FALSE)
dat_phys_marm$Species <- factor(dat_phys_marm$Species, ordered=FALSE)

tess_endur_lm <- lm(Log.Endurance~relevel(Species,ref="tesselatus")+SVL, dat_phys_tess)
tess_CIS3_lm <- lm(CI_State3~relevel(Species,ref="tesselatus"), dat_phys_tess)
tess_CIS4_lm <- lm(CI_State4~relevel(Species,ref="tesselatus"), dat_phys_tess)
tess_CIRCR_lm <- lm(CI_RCR~relevel(Species,ref="tesselatus"), dat_phys_tess)
tess_CIIS3_lm <- lm(CII_State3~relevel(Species,ref="tesselatus"), dat_phys_tess)
tess_CIIS4_lm <- lm(CII_State4~relevel(Species,ref="tesselatus"), dat_phys_tess)
tess_CIIRCR_lm <- lm(CII_RCR~relevel(Species,ref="tesselatus"), dat_phys_tess)

neom_endur_lm <- lm(Log.Endurance~relevel(Species,ref="neomexicanus")+SVL, dat_phys_neom)
neom_CIS3_lm <- lm(CI_State3~relevel(Species,ref="neomexicanus"), dat_phys_neom)
neom_CIS4_lm <- lm(CI_State4~relevel(Species,ref="neomexicanus"), dat_phys_neom)
neom_CIRCR_lm <- lm(CI_RCR~relevel(Species,ref="neomexicanus"), dat_phys_neom)
neom_CIIS3_lm <- lm(CII_State3~relevel(Species,ref="neomexicanus"), dat_phys_neom)
neom_CIIS4_lm <- lm(CII_State4~relevel(Species,ref="neomexicanus"), dat_phys_neom)
neom_CIIRCR_lm <- lm(CII_RCR~relevel(Species,ref="neomexicanus"), dat_phys_neom)

marm_endur_lm <- lm(Log.Endurance~relevel(Species,ref="marmoratus")+SVL, dat_phys_marm)
marm_CIS3_lm <- lm(CI_State3~relevel(Species,ref="marmoratus"), dat_phys_marm)
marm_CIS4_lm <- lm(CI_State4~relevel(Species,ref="marmoratus"), dat_phys_marm)
marm_CIRCR_lm <- lm(CI_RCR~relevel(Species,ref="marmoratus"), dat_phys_marm)
marm_CIIS3_lm <- lm(CII_State3~relevel(Species,ref="marmoratus"), dat_phys_marm)
marm_CIIS4_lm <- lm(CII_State4~relevel(Species,ref="marmoratus"), dat_phys_marm)
marm_CIIRCR_lm <- lm(CII_RCR~relevel(Species,ref="marmoratus"), dat_phys_marm)
