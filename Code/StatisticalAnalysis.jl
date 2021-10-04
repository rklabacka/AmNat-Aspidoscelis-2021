# These analyses were completed using Julia 1.1.0 and PhyloNetworks v0.14.0

# After opening Julia, prepare the environment
## -- If needed, install the packages that will be used
## using Pkg
## Pkg.add("PhyloNetworks")
## Pkg.add("PhyloPlots")
## Pkg.add("CSV")
## Pkg.add("StatsModels")
using PhyloNetworks
using PhyloPlots
using CSV
using StatsModels
using GLM
using DataFrames

# Set working directory
cd("/path/to/gitHubRepo/Aspidoscelis-AmNat-2021/Trees")

# Import species network
species_network = "AspNet_ML_Species.tre";
spe_net = readTopology(species_network)
spe_net_rooted = rootatnode!(spe_net, "septemvittatus")
plot(spe_net_rooted, :R, showGamma=false, showEdgeLength=true);


# Read in physiology data
cd("/path/to/gitHubRepo/Aspidoscelis-AmNat-2021/IndividualDat")
dat_means = CSV.read("PhysiologyData_2019_Means.csv", DataFrame)


# -----------------------------------------------------------------
# Analyze Between Sexual Modes
## Pagel's Lambda in theory captures within-species variation in the (1 - lambda) portion of the variance (see https://github.com/crsl4/PhyloNetworks.jl/issues/101)

# Physiology phylolm
PhyNetLM_endur = phylolm(@formula(LogEndurance ~ SexMode + SVL), dat_means, spe_net_rooted, y_mean_std=true, reml=false)
PhyNetLM_CIS3 = phylolm(@formula(CI_State3 ~ SexMode), dat_means, spe_net_rooted, y_mean_std=true, reml=false)
PhyNetLM_CIS4 = phylolm(@formula(CI_State4 ~ SexMode), dat_means, spe_net_rooted, y_mean_std=true, reml=false)
PhyNetLM_CIRCR = phylolm(@formula(CI_RCR ~ SexMode), dat_means, spe_net_rooted, y_mean_std=true, reml=false)
PhyNetLM_CIIS3 = phylolm(@formula(CII_State3 ~ SexMode), dat_means, spe_net_rooted, y_mean_std=true, reml=false)
PhyNetLM_CIIS4 = phylolm(@formula(CII_State4 ~ SexMode), dat_means, spe_net_rooted, y_mean_std=true, reml=false)
PhyNetLM_CIIRCR = phylolm(@formula(CII_RCR ~ SexMode), dat_means, spe_net_rooted, y_mean_std=true, reml=false)

# Endurance ~ Mito
PhyNetLM_CIS3endur  = phylolm(@formula(LogEndurance ~ CI_State3), dat_means, spe_net_rooted, y_mean_std=true, reml=false)
PhyNetLM_CIS4endur  = phylolm(@formula(LogEndurance ~ CI_State4), dat_means, spe_net_rooted, y_mean_std=true, reml=false)
PhyNetLM_CIRCRendur = phylolm(@formula(LogEndurance ~ CI_RCR), dat_means, spe_net_rooted, y_mean_std=true, reml=false)
PhyNetLM_CIIS3endur = phylolm(@formula(LogEndurance ~ CII_State3), dat_means, spe_net_rooted, y_mean_std=true, reml=false)
PhyNetLM_CIIS4endur = phylolm(@formula(LogEndurance ~ CII_State4), dat_means, spe_net_rooted, y_mean_std=true, reml=false)
PhyNetLM_CIIRCR     = phylolm(@formula(LogEndurance ~ CII_RCR), dat_means, spe_net_rooted, y_mean_std=true, reml=false)
# get r^2
r2(<insertModelHere>)
