# AmNat-Aspidoscelis-2021

Klabacka et al. (2021) measured the endurance capacity and mitochondrial respiration of five Aspidoscelis lizards, comparing sexual and asexual species and examining the relationship between endurance and mitochondrial respiration. This repository houses the data and coding files for the analyses.

# Contents
-   [Code Orientation](#code-orientation)
    -   [System Requirements](#system-requirements)
    -   [Data Analysis Walk-through](#walk-through)
-   [Sample Information](#sample-information)   
-   [Mitochondrial Physiology](#mito-physiology)
-   [Phylogenetics](#phylogenetics)

# Code Orientation

Three code files contain the commands used for data analysis

1. StatisticalAnalysis.R
    - Linear mixed-effects models for examining: + Effect of hybrid asexuality on each response variable + Differences in variation between hybrid asexual and parental sexual species + The effect of mitochondrial respiration on endurance capacity - Bayesian models to examine differences in mean-corrected variance - Linear models for subgroups based on mitochondrial history and parentage

1. StatisticalAnalysis.jl
    - Phylogenetic network linear models for examining:
        + Effect of hybrid asexuality on each response variable
        + The effect of mitochondrial respiration on endurance capacity

# System Requirements

In order to run the analyses, [julia](https://julialang.org/downloads/) and [R](https://cran.r-project.org/doc/manuals/r-release/R-admin.html) are required languages.
Within this readme, you will see command prompts '$' (for bash), '>' (for R), and 'julia>' (for julia).

Once you install julia, installing software packages is simply: 

    julia> using Pkg
    julia> Pkg.add("<package-name>")

The julia packages [PhyloNetworks](http://crsl4.github.io/PhyloNetworks.jl/latest/man/installation/), [PhyloPlots](https://github.com/cecileane/PhyloPlots.jl), [CSV](https://juliapackages.com/p/csv), [StatsModels](https://juliastats.org/StatsModels.jl/stable/), [GLM](https://juliapackages.com/p/glm), and [DataFrames](https://dataframes.juliadata.org/stable/) are required. 

Once you install R, installing software packages is typically:

    > install.packages("<package-name>")

However, some packages require a different installation approach, such as DevTools or install_github. Where this is needed, it is specified within the StatisticalAnalysis.R file.

The R packages [tidyverse](https://www.tidyverse.org/), [rstatix](https://www.rdocumentation.org/packages/rstatix/versions/0.7.0), [broom](https://www.rdocumentation.org/packages/broom/versions/0.7.10), [factoextra](https://www.rdocumentation.org/packages/factoextra/versions/1.0.7), [dplyr](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8), [nlme](https://cran.r-project.org/web/packages/nlme/index.html), [piecewiseSEM](https://cran.r-project.org/web/packages/piecewiseSEM/index.html), [boot](https://www.rdocumentation.org/packages/boot/versions/1.3-28), [MCMCglmm](https://www.rdocumentation.org/packages/MCMCglmm/versions/2.32), and [wolakR](https://github.com/matthewwolak/wolakR/) are required.

[PhyloNet](http://old-bioinfo.cs.rice.edu/phylonet/#Downloads) (a standalone software package) is also required.


# Data Analysis Walk-through

1- Within terminal navigate to the Code directory of the repository:

    $ cd AmNat-Aspidoscelis-2021/Code

Within this directory you can see the code files. I recommend opening them and running them within an interactive R environment (as shown in this walk-through)

2- Begin R session:

    $ R

3- Read in and prepare the data: 

    > dat_phys_indiv <- read.csv("../SampleInformation/PhysiologyData_2019_Individuals.csv")
    > dat_phys_indiv$scSVL <- scale(dat_phys_indiv$SVL)
    > dat_phys_indiv$SexualMode <- as.factor(dat_phys_indiv$SexualMode)
    > head(dat_phys_indiv)

4- Get summary statistics from data:
    
    > library(tidyverse)
    > dat_phys_indiv %>% group_by(Sex) %>% group_by(Species) %>% summarize(m=mean(CII_State3), sd=sd(CII_State3), n=n(), ci=sd / sqrt(n))

5- Create linear mixed-effects model and fit data to model:
    
    > library(nlme)
    > endur_lme_1a <- lme(Endurance ~ SexualMode + scSVL, data = dat_phys_indiv, random = ~ 1 | Species) 
    > summary(endur_lme_1a)

This is just one of the models from the code used as an example. 
For all linear mixed-effects models, see StatisticalAnalyses.R

6- See if data better fits model with differing residual variation:
    
    > endur_lme_2a <- lme(Endurance ~ SexualMode + scSVL, data = dat_phys_indiv, random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
    > # compare models
    > anova(endur_lme_1a, endur_lme_2a)
    > # low p-value supports model 2a

7- See if log-transformed data better fits the model:
    
    > endur_lme_2b <- lme(Log.Endurance ~ SexualMode + scSVL, data = dat_phys_indiv, random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
    > # plot residuals from models
    > plot(endur_lme_2a)
    > plot(endur_lme_2b)
    > # clustering noticably decreased in model 2b

8- Perform bootstrapping to get confidence intervals of residual standard deviations:

    > library(boot)
    > # Create boot function containing linear mixed-effects model
    > bootFunc=function(bootData,repeats){
    >   tryCatch({
    >   # Fit model
    >   boot_lme <- lme(Log.Endurance ~ SexualMode + scSVL, data = bootData[repeats,], random = ~ 1 | Species, weights = varIdent(form = ~ 1 | SexualMode))
    >   # Extract standard deviation ratios
    >   boot_sdrs <- coef(boot_lme$modelStruct$varStruct, unconstrained = FALSE)
    >   # Multiply standard deviation ratios by overall residual standard deviation
    >   boot_sds <- c(1, boot_sdrs) * boot_lme$sigma
    >   return(boot_sds)
    >   },
    >   error = function(err) {return(NA)}
    >   )
    > }
    > # call function
    > sds <- boot(dat_phys_indiv, bootFunc, R=1000)
    > sds
  
Steps 6 & 8 should be performed for each variable

9- Examine the relationship between endurance and mitochondrial respiration: 

    > all_EndurCIS3_lme<-  lme(Log.Endurance ~ CI_State3, random = ~ 1 | Species, data=dat_phys_indiv)
    > summary(all_EndurCIS3_lme)
    > rsquared(all_EndurCIS3_lme)
    > # use marginal 

This step should be performed for each variable of mito respiration

10- Now load librares to use a Bayesian approach to examine differences in variability between mean-corrected values of variation:

    > library(MCMCglmm)
    > library(wolakR)

11 - Set non-informative priors for residual variances:

    > # single variance value
    > pr1 <- list(R = list(V = 1e-12, nu = -2), G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = c(0), alpha.V = diag(1)*1000)))
    > # two variance values
    > pr2 <- list(R = list(V = diag(2)*1e-12, nu = -1), G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = c(0), alpha.V = diag(1)*1000)))

12- Set MCMC chain parameters:

    > nsample <- 1000
    > BURN <- 3000; THIN <- 50; NITT <- BURN + nsample * THIN

13- Create model:

    > endur_b <- MCMCglmm(Log.Endurance ~ SexualMode + scSVL,
    >   random = ~ Species,
    >   rcov = ~ idh(SexualMode):units,
    >   data = dat_phys_indiv,
    >   prior = pr2,
    >   nitt = NITT, burnin = BURN, thin = THIN)

14- Plot standard deviations for each repro mode:

    > endur_b_sdPost <- sqrt(endur_b$VCV)
    > # Get mean values
    > mean(endur_b_sdPost[, 2])
    > mean(endur_b_sdPost[, 3])
    > # Get confidence intervals
    > endur_b_sd_hpd <- HPDinterval(endur_b_sdPost)
    > # Plot posterior distributions
    > par(mfrow=c(2,1))
    > endur_asex_sd <- postPlot(endur_b_sdPost[, 2], xlim = c(0, 0.40))
    > # add sd determined from frequentist approach
    > abline(v=0.09, col="red", lwd=3)
    > endur_sex_sd <- postPlot(endur_b_sdPost[, 3], xlim = c(0, 0.40))
    > abline(v=0.18, col="red", lwd=3)

Look at plot!
The two plots show the posterior distribution of standard deviations for asexual species (top) and sexual species (bottom)


15- Plot differences between coefficient of variation:

    > endur_cvdiffPost <- ((sqrt(endur_b$VCV[, 3]) / (endur_b$Sol[, 1] + endur_b$Sol[, 2]) - sqrt(endur_b$VCV[, 2]) / endur_b$Sol[, 1]))
    > endur_hpd <- HPDinterval(endur_cvdiffPost)
    > endur_b_postPlot <- postPlot(endur_cvdiffPost, xlim=c(-0.5,0.5),plotHist = FALSE)
    > endur_b_post <- mean(endur_cvdiffPost > 0) # posterior probability that asexual coefficient of variance < sexual
    > abline(v=0, col="red", lwd=3)

Look at plot!
This is the difference in posterior distributions between the coefficient of variation.
A distribution whose center is greater than 0 shows greater mean-corrected variation in the sexual species

Steps 13-15 should be performed for each variable of mito respiration

16- Subset and prepare data for within-group linear modeling:

    > dat_phys_tess <- subset(dat_phys_indiv, Species == 'tesselatus' | Species == 'marmoratus' | Species == 'septemvittatus')
    > # convert species to factor
    > dat_phys_tess$Species <- factor(dat_phys_tess$Species, ordered=FALSE)

17- Examine effect of species on log.endurance:

    > tess_endur_lm <- lm(Log.Endurance~relevel(Species,ref="tesselatus")+SVL, dat_phys_tess) 
    > summary(tess_endur_lm)

Steps 17 and 18 should be done for each subgroup and each variable

18- Close R interactive session:

    > quit()

If given option to save workspace, you can say N

19- Begin interactive julia session:
 
    $ julia

20- Import species network:

    julia> cd("Phylogenetics/PhyloNet")
    julia> using PhyloNetworks
    julia> species_network = "network.tre";
    julia> spe_net = readTopology(species_network)
    julia> # root tree
    julia> spe_net_rooted = rootatnode!(spe_net, "septemvittatus")

21- Import physiology data:

    julia> using CSV
    julia> using DataFrames
    julia> cd("../../SampleInformation")
    julia> dat_means = CSV.read("PhysiologyData_2019_Means.csv", DataFrame)

22- Examine effect of hybrid asexuality on endurance:
    
    julia> using StatsModels
    julia> PhyNetLM_endur = phylolm(@formula(LogEndurance ~ SexMode + SVL), dat_means, spe_net_rooted, y_mean_std=true, reml=false)

Step 22 should be done for each variable

23- Examine effect of mitochondrial respiration on endurance:

    julia> using GLM
    julia> PhyNetLM_CIS3endur  = phylolm(@formula(LogEndurance ~ CI_State3), dat_means, spe_net_rooted, y_mean_std=true, reml=false)
    julia> r2(PhyNetLM_CIS3endur)

Step 23 should be done for each variable of mito respiration

# Sample Information

1.  PhysiologyData_2019_Individuals.csv
    - This file contains all individuals (rows) and their values for each respective variable (columns).
    - These data were used for the linear mixed-effects models.
    - Variables:
        + **FIELD** : Field ID
        + **Species** : Specific epithet (all individuals belong to genus Aspidoscelis) 
        + **Sex**
        + **SexualMode** : Mode of reproduction (sexual / asexual)
        + **STATE**
        + **COUNTY**
        + **CITY**
        + **LOCAL** : Specific location where individual was found
        + **LAT** : Latitude (decimal coordinate)
        + **LONG** : Longitude (decimal coordinate)
        + **Order** : Order lizards were ran within a day
        + **DateRan** : Assigned day lizards were ran
        + **Endurance** : Length of time (minutes) lizard ran
        + **Log.Endrance** : Log(10)-transformed Endurance value
        + **MitoDate** : Assigned day lizards were euthanized and mito respiration measured 
        + **EuthanasiaTime** : Time (24 hr) lizard was euthanized
        + **DisLength** : Length of time (minutes) from euthanasia to removal of all skeletal muscle from fore- and hindlimbs
        + **MitoMuscleMass** : Mass (g) of muscle used for mitochondrial respiration
        + **TotalMuscleMass** : Total mass (g) of skeletal muscle from fore- and hindlimbs
        + **TotalMass** : Total mass (g) of lizard
        + **SVL** : Snout-vent length
        + **CI_State3** : State 3 respiration (nmoles oxygen consumed per minute per milligram of protein) initiated through mitochondrial Complex I
        + **CI_State4** : State 4 respiration (nmoles oxygen consumed per minute per milligram of protein) initiated through mitochondrial Complex I
        + **CI_RCR** : State 3 / State 4 for respiration initiated through mitochondrial Complex I
        + **CII_State3** : State 3 respiration (nmoles oxygen consumed per minute per milligram of protein) initiated through mitochondrial Complex II
        + **CII_State4** : State 4 respiration (nmoles oxygen consumed per minute per milligram of protein) initiated through mitochondrial Complex II
        + **CII_RCR** : State 3 / State 4 for respiration initiated through mitochondrial Complex II

1.  PhysiologyData_2019_Means.csv 
    - This file contains all species (rows) and their variable values (columns).
    - These data were used for the phylogenetic network linear models.
    - Where values are means of variables, see the variable descriptions above
    - Variables:
        + **tipNames** : Specific epithet
        + **Endurance** : Mean Endurance
        + **Log.Endurance** : Log(10)-transformed mean endurance capacity
        + **DisLength** : Mean DisLength 
        + **MitoMuscleMass** : Mean MitoMuscleMass
        + **TotalMuscleMass** : Mean TotalMuscleMass
        + **TotalMass** : Mean TotalMass
        + **SVL** : Mean SVL
        + **CI_State3** : Mean CI_State3
        + **CI_State3_n** : Number of individuals measured for CI_State3
        + **CI_State3_sd** : Standard deviation for CI_State3
        + **CI_State4** : Mean CI_State4
        + **CI_State4_n** : Number of individuals measured for CI_State3
        + **CI_State4_sd** : Standard deviation for CI_State3
        + **CI_RCR** : Mean CI_RCR
        + **CI_RCR_n** : Number of individuals measured for CI_State3
        + **CI_RCR_sd** : Standard deviation for CI_State3
        + **CII_State3** : Mean CII_State3
        + **CII_State3_n** : Number of individuals measured for CI_State3
        + **CII_State3_sd** : Standard deviation for CI_State3
        + **CII_State4** : Mean CII_State4
        + **CII_State4_n** : Number of individuals measured for CII_State3
        + **CII_State4_sd** : Standard deviation for CII_State3
        + **CII_RCR** : Mean CII_RCR
        + **CII_RCR_n** : Number of individuals measured for CII_State3
        + **CII_RCR_sd** : Standard deviation for CII_State3

# Mitochondrial Physiology

1.  RawMitoData/\*.csv 
    - Plotted oxygen consumption for isolated mitochondria of each individual

1.  Mitochondria Bradford Data.xlsx
    - Protein concentrations from Bradford assay of isolated mitochondria (for respiration normalization)

1.  RCR calculations_Lizard Endruance Study 2019.xlsx
    - Normalization of respiration data and calculation of RCR

# Phylogenetics

## IQ-Tree

1. mito.nex
    - Sequence alignment of mitochondrial data used for phylogeny estimation

1.  mito.tre 
    - Mitochondrial consensus tree estimated from IQ-Tree

## PhyloNet

1. PhyloNet_ML.nex
    - Input for Phylonet with gene trees and parameters

1. network.tre
    - Phylogenetic network estimated from PhyloNet
