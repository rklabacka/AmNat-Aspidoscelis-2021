# Load libraries
library(tidyverse)
library(rstatix)
library(broom)
library(MuMIn)
library(factoextra)
library(lme4)
library(lmerTest)

## All phylogenetic network analyses performed in the Julia package PhyloNetworks (see AspidoscelisReticAnalysis.jl)

setwd("/path/to/gitHubRepo/Aspidoscelis-AmNat-2021/IndividualData")
# Read in physiology data
dat_phys_indiv <- read.csv("PhysiologyData_2019_Individuals.csv")

# ------------
# Mixed-Effects Models
# ------------

# Effect of SexMode on response variables
all_endur_lmer<-lmer(Log.Endurance~SexualMode+SVL +(1|Species), data=dat_phys_indiv, REML=T)
anova(all_endur_lmer)
var.test(residuals(all_endur_lmer) ~ as.factor(dat_phys_indiv$SexualMode))


all_CIS3_lmer<-lmer(CI_State3~SexualMode +(1|Species), data=dat_phys_indiv, REML=T)
anova(all_CIS3_lmer)
leveneTest(residuals(all_CIS3_lmer) ~ as.factor(dat_phys_indiv$SexualMode))
var.test(residuals(all_CIS3_lmer) ~ as.factor(dat_phys_indiv$SexualMode))

all_CIS4_lmer<-lmer(CI_State4~SexualMode +(1|Species), data=dat_phys_indiv, REML=T)
anova(all_CIS4_lmer)
var.test(residuals(all_CIS4_lmer) ~ as.factor(dat_phys_indiv$SexualMode))

all_CIRCR_lmer<-lmer(CI_RCR~SexualMode +(1|Species), data=dat_phys_indiv, REML=T)
anova(all_CIRCR_lmer)
var.test(residuals(all_CIRCR_lmer) ~ as.factor(dat_phys_indiv$SexualMode))

all_CIIS3_lmer<-lmer(CII_State3~SexualMode +(1|Species), data=dat_phys_indiv, REML=T)
anova(all_CIIS3_lmer)
var.test(residuals(all_CIIS3_lmer) ~ as.factor(dat_phys_indiv$SexualMode))

all_CIIS4_lmer<-lmer(CII_State4~SexualMode +(1|Species), data=dat_phys_indiv, REML=T)
anova(all_CIIS4_lmer)
var.test(residuals(all_CIIS4_lmer) ~ as.factor(dat_phys_indiv$SexualMode))

all_CIIRCR_lmer<-lmer(CII_RCR~SexualMode +(1|Species), data=dat_phys_indiv, REML=T)
anova(all_CIIRCR_lmer)
var.test(residuals(all_CIIRCR_lmer) ~ as.factor(dat_phys_indiv$SexualMode))


# Endurancce-Respiration Relationship
all_EndurCIS3_lmer<-  lmer(CI_State3  ~ Log.Endurance +(1|Species), data=dat_phys_indiv, REML=T)
all_EndurCIS4_lmer<-  lmer(CI_State4  ~ Log.Endurance +(1|Species), data=dat_phys_indiv, REML=T)
all_EndurCIRCR_lmer<- lmer(CI_RCR     ~ Log.Endurance +(1|Species), data=dat_phys_indiv, REML=T)
all_EndurCIIS3_lmer<- lmer(CII_State3 ~ Log.Endurance +(1|Species), data=dat_phys_indiv, REML=T)
all_EndurCIIS4_lmer<- lmer(CII_State4 ~ Log.Endurance +(1|Species), data=dat_phys_indiv, REML=T)
all_EndurCIIRCR_lmer<-lmer(CII_RCR      ~ Log.Endurance +(1|Species), data=dat_phys_indiv, REML=T)
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
