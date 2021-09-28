# Load libraries
library(tidyverse)
library(rstatix)
library(broom)
library(factoextra)
library(dplyr)
library(nlme) # for building of linear model using Frequentist approach
library(piecewiseSEM) # for calculating r^2
# - Use Bayesian approach to mean-correct variance values and compare
library(MCMCglmm) # For linear models using Bayesian approach
library(wolakR) #<-- handy functions to work with MCMCglmm models
### --- Can use `remotes` package to install from GitHub:
### --- remotes::install_github("matthewwolak/wolakR")


## All phylogenetic network analyses performed in the Julia package PhyloNetworks (see AspidoscelisReticAnalysis.jl)

setwd("/path/to/gitHubRepo/Aspidoscelis-AmNat-2021/IndividualData")
# Read in physiology data
dat_phys_indiv <- read.csv("PhysiologyData_2019_Individuals.csv")

# - Prepare data
## -- Scale and center SVL
dat_phys_indiv$scSVL <- scale(dat_phys_indiv$SVL)  #<-- work with this from now on
## -- Convert Sexual Mode variable to factor
dat_phys_indiv$SexualMode <- as.factor(dat_phys_indiv$SexualMode)
# Get summary statistics:
dat_phys_indiv %>% group_by(Sex) %>% group_by(Species) %>% summarize(m=mean(Endurance), sd=sd(Endurance), n=n(), ci=sd / sqrt(n))

# --------------------------------------------------------------------------------------------------------
# Examine effect of sexual mode on response variables and perform model selection for one or two residual variance parameters

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


# *************
#  CI State 3
# *************
## -- Model 1: No differences in variance between sexual modes
CIS3_lme_1 <- lme(CI_State3 ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species)
## -- Model 2: Differences in variance between sexual modes
CIS3_lme_2 <- lme(CI_State3 ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
# Compare the two models: Hypothesis test for different (residual) variances
anova(CIS3_lme_1, CIS3_lme_2) # model 1 preferred, p = 0.8405
## -- See details of selected model
summary(CIS3_lme_1)

# *************
#  CI State 4
# *************
## -- Model 1: No differences in variance between sexual modes
CIS4_lme_1 <- lme(CI_State4 ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species)
## -- Model 2: Differences in variance between sexual modes
CIS4_lme_2 <- lme(CI_State4 ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
# Compare the two models: Hypothesis test for different (residual) variances
anova(CIS4_lme_1, CIS4_lme_2) # model 2 preferred, p = 7e-04
## -- See details of selected model
summary(CIS4_lme_2)

# *************
#  CI RCR
# *************
## -- Model 1: No differences in variance between sexual modes
CIRCR_lme_1 <- lme(CI_RCR ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species)
## -- Model 2: Differences in variance between sexual modes
CIRCR_lme_2 <- lme(CI_RCR ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
# Compare the two models: Hypothesis test for different (residual) variances
anova(CIRCR_lme_1, CIRCR_lme_2) # model 1 preferred, p = 0.37
## -- See details of selected model
summary(CIRCR_lme_1)


# *************
#  CII State 3
# *************
## -- Model 1: No differences in variance between sexual modes
CIIS3_lme_1 <- lme(CII_State3 ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species)
## -- Model 2: Differences in variance between sexual modes
CIIS3_lme_2 <- lme(CII_State3 ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
# Compare the two models: Hypothesis test for different (residual) variances
anova(CIIS3_lme_1, CIIS3_lme_2) # model 2 preferred, p = 0.047
## -- See details of selected model
summary(CIIS3_lme_2)

# *************
#  CII State 4
# *************
## -- Model 1: No differences in variance between sexual modes
CIIS4_lme_1 <- lme(CII_State4 ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species)
## -- Model 2: Differences in variance between sexual modes
CIIS4_lme_2 <- lme(CII_State4 ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
# Compare the two models: Hypothesis test for different (residual) variances
anova(CIIS4_lme_1, CIIS4_lme_2) # model 2 preferred, p = 0.004
## -- See details of selected model
summary(CIIS4_lme_2)

# *************
#  CII RCR
# *************
## -- Model 1: No differences in variance between sexual modes
CIIRCR_lme_1 <- lme(CII_RCR ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species)
## -- Model 2: Differences in variance between sexual modes
CIIRCR_lme_2 <- lme(CII_RCR ~ SexualMode, data = dat_phys_indiv, random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
# Compare the two models: Hypothesis test for different (residual) variances
anova(CIIRCR_lme_1, CIIRCR_lme_2) # model 1 preferred, p = 0.088
## -- See details of selected model
summary(CIIRCR_lme_1)

# --------------------------------------------------------------------------------------------------------
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

# Create model
endur_b <- MCMCglmm(Log.Endurance ~ SexualMode,
  random = ~ Species,
  rcov = ~ idh(SexualMode):units,
  data = dat_phys_indiv,
  prior = pr2,
  nitt = NITT, burnin = BURN, thin = THIN)
endur_cvdiffPost <- (sqrt(endur_b$VCV[, 2]) / endur_b$Sol[, 1]) - (sqrt(endur_b$VCV[, 3]) / (endur_b$Sol[, 1] + endur_b$Sol[, 2]))
endur_hpd <- HPDinterval(endur_cvdiffPost)
postPlot(endur_cvdiffPost, xlim = c(-0.5, 0.5))
endur_b_post <- mean(endur_cvdiffPost < 0) # posterior probability that asexual coefficient of variance < sexual

CIS4_b <- MCMCglmm(CI_State4 ~ SexualMode ,
  random = ~ Species,
  rcov = ~ idh(SexualMode):units,
  data = dat_phys_indiv,
  prior = pr2,
  nitt = NITT, burnin = BURN, thin = THIN)
CIS4_cvdiffPost <- (sqrt(CIS4_b$VCV[, 2]) / CIS4_b$Sol[, 1]) - (sqrt(CIS4_b$VCV[, 3]) / (CIS4_b$Sol[, 1] + CIS4_b$Sol[, 2]))
CIS4_hpd <- HPDinterval(CIS4_cvdiffPost)
postPlot(CIS4_cvdiffPost, xlim = c(-0.5, 0.5))
CIS4_b_post <- mean(CIS4_cvdiffPost < 0) # posterior probability that asexual coefficient of variance < sexual

CIIS3_b <- MCMCglmm(CII_State3 ~ SexualMode ,
  random = ~ Species,
  rcov = ~ idh(SexualMode):units,
  data = dat_phys_indiv,
  prior = pr2,
  nitt = NITT, burnin = BURN, thin = THIN)
CIIS3_cvdiffPost <- (sqrt(CIIS3_b$VCV[, 2]) / CIIS3_b$Sol[, 1]) - (sqrt(CIIS3_b$VCV[, 3]) / (CIIS3_b$Sol[, 1] + CIIS3_b$Sol[, 2]))
CIIS3_hpd <- HPDinterval(CIIS3_cvdiffPost)
postPlot(CIIS3_cvdiffPost, xlim = c(-0.5, 0.5))
CIIS3_b_post <- mean(CIIS3_cvdiffPost < 0) # posterior probability that asexual coefficient of variance < sexual

CIIS4_b <- MCMCglmm(CII_State4 ~ SexualMode ,
  random = ~ Species,
  rcov = ~ idh(SexualMode):units,
  data = dat_phys_indiv,
  prior = pr2,
  nitt = NITT, burnin = BURN, thin = THIN)
CIIS4_cvdiffPost <- (sqrt(CIIS4_b$VCV[, 2]) / CIIS4_b$Sol[, 1]) - (sqrt(CIIS4_b$VCV[, 3]) / (CIIS4_b$Sol[, 1] + CIIS4_b$Sol[, 2]))
CIIS4_hpd <- HPDinterval(CIIS4_cvdiffPost)
postPlot(CIIS4_cvdiffPost, xlim = c(-0.5, 0.5))
CIIS4_b_post <- mean(CIIS4_cvdiffPost < 0) # posterior probability that asexual coefficient of variance < sexual

# --------------------------------------------------------------------------------------------------------
# Examine the relationship between endurance and respiration
all_EndurCIS3_lme<-  lme(Log.Endurance ~ CI_State3  + SVL, random = ~ 1 | Species, data=dat_phys_indiv)
rsquared(all_EndurCIS3_lme) # Marginal: 0.34

all_EndurCIS4_lme<-  lme(Log.Endurance ~ CI_State4  + SVL, random = ~ 1 | Species, data=dat_phys_indiv)
rsquared(all_EndurCIS4_lme) # Marginal: 0.29

all_EndurCIRCR_lme<-  lme(Log.Endurance ~ CI_RCR  + SVL, random = ~ 1 | Species, data=dat_phys_indiv)
rsquared(all_EndurCIRCR_lme) # Marginal: 0.15

all_EndurCIIS3_lme<-  lme(Log.Endurance ~ CII_State3  + SVL, random = ~ 1 | Species, data=dat_phys_indiv)
rsquared(all_EndurCIIS3_lme) # Marginal: 0.32

all_EndurCIIS4_lme<-  lme(Log.Endurance ~ CII_State4  + SVL, random = ~ 1 | Species, data=dat_phys_indiv)
rsquared(all_EndurCIIS4_lme) # Marginal: 0.19

all_EndurCIIRCR_lme<-  lme(Log.Endurance ~ CII_RCR  + SVL, random = ~ 1 | Species, data=dat_phys_indiv)
rsquared(all_EndurCIIRCR_lme) # Marginal: 0.13
# -----------------------------------------------------------------------------------------------------


# ------------
# Within-group ANCOVAs
# ------------

# Subset data
dat_phys_tess <- subset(dat_phys_indiv, Species == 'tesselata' | Species == 'marmorata' | Species == 'septemvittata')
dat_phys_neom <- subset(dat_phys_indiv, Species == 'neomexicana' | Species == 'marmorata' | Species == 'inornata')
dat_phys_marm <- subset(dat_phys_indiv, Species == 'neomexicana' | Species == 'marmorata' | Species == 'tesselata')

dat_phys_tess$Species <- factor(dat_phys_tess$Species, ordered=FALSE)
dat_phys_neom$Species <- factor(dat_phys_neom$Species, ordered=FALSE)
dat_phys_marm$Species <- factor(dat_phys_marm$Species, ordered=FALSE)

tess_endur_lm <- lm(log(Endurance)~relevel(Species,ref="tesselata")+SVL, dat_phys_tess)
tess_CIS3_lm <- lm(CI_State3~relevel(Species,ref="tesselata"), dat_phys_tess)
tess_CIS4_lm <- lm(CI_State4~relevel(Species,ref="tesselata"), dat_phys_tess)
tess_CIRCR_lm <- lm(CI_RCR~relevel(Species,ref="tesselata"), dat_phys_tess)
tess_CIIS3_lm <- lm(CII_State3~relevel(Species,ref="tesselata"), dat_phys_tess)
tess_CIIS4_lm <- lm(CII_State4~relevel(Species,ref="tesselata"), dat_phys_tess)
tess_CIIRCR_lm <- lm(CII_RCR~relevel(Species,ref="tesselata"), dat_phys_tess)

neom_endur_lm <- lm(log(Endurance)~relevel(Species,ref="neomexicana")+SVL, dat_phys_neom)
neom_CIS3_lm <- lm(CI_State3~relevel(Species,ref="neomexicana"), dat_phys_neom)
neom_CIS4_lm <- lm(CI_State4~relevel(Species,ref="neomexicana"), dat_phys_neom)
neom_CIRCR_lm <- lm(CI_RCR~relevel(Species,ref="neomexicana"), dat_phys_neom)
neom_CIIS3_lm <- lm(CII_State3~relevel(Species,ref="neomexicana"), dat_phys_neom)
neom_CIIS4_lm <- lm(CII_State4~relevel(Species,ref="neomexicana"), dat_phys_neom)
neom_CIIRCR_lm <- lm(CII_RCR~relevel(Species,ref="neomexicana"), dat_phys_neom)

marm_endur_lm <- lm(log(Endurance)~relevel(Species,ref="marmorata")+SVL, dat_phys_marm)
marm_CIS3_lm <- lm(CI_State3~relevel(Species,ref="marmorata"), dat_phys_marm)
marm_CIS4_lm <- lm(CI_State4~relevel(Species,ref="marmorata"), dat_phys_marm)
marm_CIRCR_lm <- lm(CI_RCR~relevel(Species,ref="marmorata"), dat_phys_marm)
marm_CIIS3_lm <- lm(CII_State3~relevel(Species,ref="marmorata"), dat_phys_marm)
marm_CIIS4_lm <- lm(CII_State4~relevel(Species,ref="marmorata"), dat_phys_marm)
marm_CIIRCR_lm <- lm(CII_RCR~relevel(Species,ref="marmorata"), dat_phys_marm)
