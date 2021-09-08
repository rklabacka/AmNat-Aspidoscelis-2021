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

# Set working directory
cd("/path/to/gitHubRepo/Aspidoscelis-AmNat-2021/Trees")

# Import species network
species_network = "AspNet_ML_Species.tre";
spe_net = readTopology(species_network)
spe_net_rooted = rootatnode!(spe_net, "septemvittatus")
plot(spe_net_rooted, :R, showGamma=false, showEdgeLength=true);


# Read in physiology data
cd("/path/to/gitHubRepo/Aspidoscelis-AmNat-2021/IndividualDat")
dat_means = CSV.read("PhysiologyData_2019_Means.csv")


# -----------------------------------------------------------------
# Analyze Between Sexual Modes
## Pagel's Lambda in theory captures within-species variation in the (1 - lambda) portion of the variance (see https://github.com/crsl4/PhyloNetworks.jl/issues/101)

# Physiology Phylonetworklm
PhyNetLM_endur = phyloNetworklm(@formula(LogEndurance ~ SexMode + SVL), dat_means, spe_net_rooted, model="lambda")
PhyNetLM_CIS3 = phyloNetworklm(@formula(CI_State3 ~ SexMode), dat_means, spe_net_rooted, model="lambda")
PhyNetLM_CIS4 = phyloNetworklm(@formula(CI_State4 ~ SexMode), dat_means, spe_net_rooted, model="lambda")
PhyNetLM_CIRCR = phyloNetworklm(@formula(CI_RCR ~ SexMode), dat_means, spe_net_rooted, model="lambda")
PhyNetLM_CIIS3 = phyloNetworklm(@formula(CII_State3 ~ SexMode), dat_means, spe_net_rooted, model="lambda")
PhyNetLM_CIIS4 = phyloNetworklm(@formula(CII_State4 ~ SexMode), dat_means, spe_net_rooted, model="lambda")
PhyNetLM_CIIRCR = phyloNetworklm(@formula(CII_RCR ~ SexMode), dat_means, spe_net_rooted, model="lambda")

# Endurance ~ Mito
PhyNetLM_CIS3endur  = phyloNetworklm(@formula(LogEndurance ~ CI_State3  + SVL), dat_means, spe_net_rooted, model="lambda")
PhyNetLM_CIS4endur  = phyloNetworklm(@formula(LogEndurance ~ CI_State4  + SVL), dat_means, spe_net_rooted, model="lambda")
PhyNetLM_CIRCRendur = phyloNetworklm(@formula(LogEndurance ~ CI_RCR     + SVL), dat_means, spe_net_rooted, model="lambda")
PhyNetLM_CIIS3endur = phyloNetworklm(@formula(LogEndurance ~ CII_State3 + SVL), dat_means, spe_net_rooted, model="lambda")
PhyNetLM_CIIS4endur = phyloNetworklm(@formula(LogEndurance ~ CII_State4 + SVL), dat_means, spe_net_rooted, model="lambda")
PhyNetLM_CIIRCR     = phyloNetworklm(@formula(LogEndurance ~ CII_RCR    + SVL), dat_means, spe_net_rooted, model="lambda")
# get r^2
r2(<insertModelHere>)
