# AmNat-Aspidoscelis-2021

Klabacka et al. (2021) measured the endurance capacity and mitochondrial respiration of five Aspidoscelis lizards, comparing sexual and asexual species and examining the relationship between endurance and mitochondrial respiration. This repository houses the data and coding files for the analyses.

## Code Orientation

Three code files contain the commands used for data analysis

1. StatisticalAnalysis.R
    - Linear mixed-effects models for examining: + Effect of hybrid asexuality on each response variable + Differences in variation between hybrid asexual and parental sexual species + The effect of mitochondrial respiration on endurance capacity - Bayesian models to examine differences in mean-corrected variance - Linear models for subgroups based on mitochondrial history and parentage

1. StatisticalAnalysis.jl
    - Phylogenetic network linear models for examining:
        + Effect of hybrid asexuality on each response variable
        + The effect of mitochondrial respiration on endurance capacity

### Analyses System Requirements

In order to run the analyses, [julia](https://julialang.org/downloads/) and [R](https://cran.r-project.org/doc/manuals/r-release/R-admin.html) are required languages.

Once you install julia, installing software packages is simply: 

    $ using Pkg
    $ Pkg.add("<package-name>")

The julia packages [PhyloNetworks](http://crsl4.github.io/PhyloNetworks.jl/latest/man/installation/), [PhyloPlots](https://github.com/cecileane/PhyloPlots.jl), [CSV](https://juliapackages.com/p/csv), [StatsModels](https://juliastats.org/StatsModels.jl/stable/), [GLM](https://juliapackages.com/p/glm), and [DataFrames](https://dataframes.juliadata.org/stable/) are required. 

Once you install R, installing software packages is typically:

    $ install.packages("<package-name>")

However, some packages require a different installation approach, such as DevTools or install_github. Where this is needed, it is specified within the StatisticalAnalysis.R file.

The R packages [tidyverse](https://www.tidyverse.org/), [rstatix](https://www.rdocumentation.org/packages/rstatix/versions/0.7.0), [broom](https://www.rdocumentation.org/packages/broom/versions/0.7.10), [factoextra](https://www.rdocumentation.org/packages/factoextra/versions/1.0.7), [dplyr](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8), [nlme](https://cran.r-project.org/web/packages/nlme/index.html), [piecewiseSEM](https://cran.r-project.org/web/packages/piecewiseSEM/index.html), [boot](https://www.rdocumentation.org/packages/boot/versions/1.3-28), [MCMCglmm](https://www.rdocumentation.org/packages/MCMCglmm/versions/2.32), and [wolakR](https://github.com/matthewwolak/wolakR/) are required.

[PhyloNet](http://old-bioinfo.cs.rice.edu/phylonet/#Downloads) (a standalone software package) is also required.


### Data Analysis Walk-through

1.  Within terminal navigate to the Code directory of the repository

    $ cd AmNat-Aspidoscelis-2021/Code

    Within this directory you can see the code files. I recommend opening them and running them within an interactive R environment (as shown in this walk-through)

1.  Begin R session

    $ R

1.  Read in and prepare the data 

    > dat_phys_indiv <- read.csv("../SampleInformation/PhysiologyData_2019_Individuals.csv")
    > dat_phys_indiv$scSVL <- scale(dat_phys_indiv$SVL)
    > dat_phys_indiv$SexualMode <- as.factor(dat_phys_indiv$SexualMode)
    > head(dat_phys_indiv)

1.  Get summary statistics from data
    
    > library(tidyverse)
    > dat_phys_indiv %>% group_by(Sex) %>% group_by(Species) %>% summarize(m=mean(CII_State3), sd=sd(CII_State3), n=n(), ci=sd / sqrt(n))

1.  Create linear mixed-effects model and fit data to model
    note: This is just one of the models from the code used as an example. 
    For all linear mixed-effects models, see StatisticalAnalyses.R
    
    > library(nlme)
    > endur_lme_1a <- lme(Endurance ~ SexualMode + scSVL, data = dat_phys_indiv, random = ~ 1 | Species) 
    > summary(endur_lme_1a)

1.  See if log-transformed data better fits the model
    
    > endur_lme_1b <- lme(Log.Endurance ~ SexualMode + scSVL, data = dat_phys_indiv, random = ~ 1 | Species)

## SampleInformation

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

## Mitochondrial Physiology

1.  RawMitoData/\*.csv 
    - Plotted oxygen consumption for isolated mitochondria of each individual

1.  Mitochondria Bradford Data.xlsx
    - Protein concentrations from Bradford assay of isolated mitochondria (for respiration normalization)

1.  RCR calculations_Lizard Endruance Study 2019.xlsx
    - Normalization of respiration data and calculation of RCR

## Phylogenetics

### IQ-Tree

1. mito.nex
    - Sequence alignment of mitochondrial data used for phylogeny estimation

1.  mito.tre 
    - Mitochondrial consensus tree estimated from IQ-Tree

### PhyloNet

1. PhyloNet_ML.nex
    - Input for Phylonet with gene trees and parameters

1. network.tre
    - Phylogenetic network estimated from PhyloNet
