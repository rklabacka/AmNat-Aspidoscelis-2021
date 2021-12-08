# AmNat-Aspidoscelis-2021

Klabacka et al. (2021) measured the endurance capacity and mitochondrial respiration of five Aspidoscelis lizards, comparing sexual and asexual species and examining the relationship between endurance and mitochondrial respiration. This repository houses the data and coding files for the analyses.

# Data Orientation

The data are contained within several files/directories:

1.  PhysiologyData_2019_Individuals.csv
    This file contains all individuals (rows) and their values for each respective variable (columns).
1.  PhysiologyData_2019_Means.csv 
1.  RawMitoData/\*.csv 
1.  Mitochondria Bradford Data.xlsx
1.  RCR calculations_Lizard Endruance Study 2019.xlsx

# Analyses System Requirements

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



