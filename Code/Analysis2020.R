# Load libraries
library(tidyverse)
library(rstatix)
library(broom)
library(MuMIn)
library(factoextra)
library(lme4)
library(lmerTest)
library(dplyr)

## All phylogenetic network analyses performed in the Julia package PhyloNetworks (see AspidoscelisReticAnalysis.jl)

setwd("/path/to/gitHubRepo/Aspidoscelis-AmNat-2021/IndividualData")
# Read in physiology data
dat_phys_indiv <- read.csv("PhysiologyData_2019_Individuals.csv")
# To get summary statistics:
dat_phys_indiv %>% group_by(Sex) %>% group_by(Species) %>% summarize(m=mean(Endurance), sd=sd(Endurance), n=n(), ci=sd / sqrt(n))
# Add means to dataframe:
dat_phys_indiv <- dat_phys_indiv %>% 
  group_by(Species) %>% 
  mutate(AvgEndur = mean(Endurance))

# ------------
# Mixed-Effects Models
# ------------

# Effect of SexMode on response variables
## Note: we use the 
all_endur_lmer<-lmer(Log.Endurance~SexualMode+SVL +(1|Species), data=dat_phys_indiv, REML=T, control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
# check redidual distribution
plot(all_endur_lmer)
anova(all_endur_lmer)
var.test(residuals(all_endur_lmer) ~ as.factor(dat_phys_indiv$SexualMode))
# Add residuals to dataframe
dat_phys_indiv$endur_resid<-residuals(all_endur_lmer)
# Obtain mean-scaled deviation measure
dat_phys_indiv$endur_resid_scaled<-dat_phys_indiv$endur_resid / dat_phys_indiv$AvgEndur


all_CIS3_lmer<-lmer(CI_State3~SexualMode +(1|Species), data=dat_phys_indiv, REML=T, control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anova(all_CIS3_lmer)
leveneTest(residuals(all_CIS3_lmer) ~ as.factor(dat_phys_indiv$SexualMode))
var.test(residuals(all_CIS3_lmer) ~ as.factor(dat_phys_indiv$SexualMode))

all_CIS4_lmer<-lmer(CI_State4~SexualMode +(1|Species), data=dat_phys_indiv, REML=T, control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anova(all_CIS4_lmer)
var.test(residuals(all_CIS4_lmer) ~ as.factor(dat_phys_indiv$SexualMode))

all_CIRCR_lmer<-lmer(CI_RCR~SexualMode +(1|Species), data=dat_phys_indiv, REML=T, control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anova(all_CIRCR_lmer)
var.test(residuals(all_CIRCR_lmer) ~ as.factor(dat_phys_indiv$SexualMode))

all_CIIS3_lmer<-lmer(CII_State3~SexualMode +(1|Species), data=dat_phys_indiv, REML=T, control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anova(all_CIIS3_lmer)
var.test(residuals(all_CIIS3_lmer) ~ as.factor(dat_phys_indiv$SexualMode))

all_CIIS4_lmer<-lmer(CII_State4~SexualMode +(1|Species), data=dat_phys_indiv, REML=T, control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anova(all_CIIS4_lmer)
var.test(residuals(all_CIIS4_lmer) ~ as.factor(dat_phys_indiv$SexualMode))

all_CIIRCR_lmer<-lmer(CII_RCR~SexualMode +(1|Species), data=dat_phys_indiv, REML=T, control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anova(all_CIIRCR_lmer)
var.test(residuals(all_CIIRCR_lmer) ~ as.factor(dat_phys_indiv$SexualMode))


# Endurance-Respiration Relationship
all_EndurCIS3_lmer<-  lmer(Log.Endurance ~ CI_State3  + SVL + (1|Species), data=dat_phys_indiv, REML=T, control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
all_EndurCIS4_lmer<-  lmer(Log.Endurance ~ CI_State4  + SVL + (1|Species), data=dat_phys_indiv, REML=T, control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
all_EndurCIRCR_lmer<- lmer(Log.Endurance ~ CI_RCR     + SVL + (1|Species), data=dat_phys_indiv, REML=T, control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
all_EndurCIIS3_lmer<- lmer(Log.Endurance ~ CII_State3 + SVL + (1|Species), data=dat_phys_indiv, REML=T, control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
all_EndurCIIS4_lmer<- lmer(Log.Endurance ~ CII_State4 + SVL + (1|Species), data=dat_phys_indiv, REML=T, control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
all_EndurCIIRCR_lmer<-lmer(Log.Endurance ~ CII_RCR    + SVL + (1|Species), data=dat_phys_indiv, REML=T, control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
r.squaredGLMM(<insertModelHere>) # uses MuMIn package to get marginal (takes fixed effects into account) and conditional (which takes random and fixed effects into account) R^2. I'm using conditional.



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
