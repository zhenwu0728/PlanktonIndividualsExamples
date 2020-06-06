# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Julia 1.4.0
#     language: julia
#     name: julia-1.4
# ---

# + [markdown] slideshow={"slide_type": "slide"}
# # 0-Dimensional 1-species example
#
# - Introduce equations used in the model
#
# - Go through the whole processes of runing the model and visualization
#

# + [markdown] slideshow={"slide_type": "fragment"}
# ## For more documentation
# - <https://zhenwu0728.github.io/PlanktonIndividuals.jl/dev/> for PlanktonIndividuals.jl
# - <https://github.com/JuliaPy/PyPlot.jl> for PyPlot.jl
# - <https://docs.julialang.org/en/v1/> for Julia

# + [markdown] slideshow={"slide_type": "slide"}
# ![skematic](PI_Quota.jpeg)

# + [markdown] slideshow={"slide_type": "subslide"}
# ## Photosynthesis
#
# Basically, we follow Geider et al. 1998 for the parameterization of photosynthesis, but without nutrient limitation.
#
#
# $PP=PP_{max}\cdot (1-e^{\frac{-\alpha \cdot I\cdot Chl}{PP_{max}\cdot C}})$
#
# **Unit: mmolC/cell/s**
#
# $PP_{max}$ is scaled by a power-law relationship of cell size
#
# ## Nutrient Uptake
#
# Include intracellular nutrient limitation:
#
# $RegQ_i=\bigg[\frac{R_{iC}^{max}-Q_i}{R_{iC}^{max}-R_{iC}^{min}}\bigg]_0^1$
#
# $V_i=V_i^{max}\cdot regQ_i\cdot\frac{[i]}{[i]+K_i^{sat}}$
#
# $i$ denotes $NH_4$, $NO_3$, $PO_4$.
#
# **Unit: mmolN/cell/s**
#
#

# + [markdown] slideshow={"slide_type": "subslide"}
# ### Metabolic Partitioning
# $\beta=\frac{a\cdot Vol_{cell}^b}{1+a\cdot Vol_{cell}^b}$
#
# $BioSynC = \beta\cdot k_{mtb}\cdot Q_C^R$
#
# $MaintenC=(1-\beta)\cdot k_{mtb}\cdot Q_C^R$
#
# $BioSynN = k_{mtb}\cdot Q_N^R/R_{NC}$
#
# $BioSynP = k_{mtb}\cdot Q_P^R/R_{PC}$
#
# ### Compute biosynthesis rate and excretion rate
# $BioSyn=min(BioSynC,BioSynN,BioSynP)$
#
# $ExcretC=BioSynC - BioSyn$

# + [markdown] slideshow={"slide_type": "slide"}
# ## 1. Import Software, i.e. `Julia` Packages

# + slideshow={"slide_type": "fragment"}
using PlanktonIndividuals, Oceananigans, YAML, KernelDensity, DelimitedFiles, PyPlot

# + [markdown] slideshow={"slide_type": "slide"}
# ## 2. Initialize the model
# 1. Initialize model grid
# 2. read parameters (Stored in `params.yml`)
# 3. Initialize nutrient field
# 4. Initialize individuals
# 5. User-defined diagnostics

# + slideshow={"slide_type": "fragment"}
Nx = 25       # Number of grid points in x and y
Nz = 50       # Number of grid points in z
Δz = 4.0      # Grid spacing in x, y, z (meters)
grid = RegularCartesianGrid(size=(1,1,1), extent=(Δz*Nx, Δz*Nx, Δz*Nz))

# + slideshow={"slide_type": "subslide"}
resultpath = PrepRunDir()
phy_grid = read_Ogrids(grid);
nday = 2
RunParam.nTime = 1440*nday # simulated days in minute
RunParam.ΔT = 60 # model time step: 60 seconds

# update paramters
phy_params = YAML.load(open("params.yml"))
update_params!(RunParam.params,phy_params)

Nsp = 1
RunParam.params["P_Nsp"] = Nsp

# Initial conditions
#           DIC  NH4   NO3   PO4   DOC  DON  DOP  POC  PON  POP  ZOO
nut_init = [2.0, 0.05, 0.10, 0.02, 1.0, 0.1, 0.1, 1.0, 0.1, 0.1, 2.0];
phy_model = PI_Model(phy_grid, RunParam;
                     nutrients = setup_nutrients(phy_grid,nut_init));

# User-defined size distribution diagnostics
size_dens = zeros(45,RunParam.nTime,phy_model.params["P_Nsp"])
xInd = collect(0.8:0.05:3.0);

# + [markdown] slideshow={"slide_type": "slide"}
# ## 3. Run the model
# - `PI_TimeStep!` updates the model at each time step
# - with user-defined diagnostics

# + slideshow={"slide_type": "fragment"}
for i in 1:RunParam.nTime
    PI_TimeStep!(phy_model, RunParam.ΔT, resultpath)
    phyts_sp = sort_species(phy_model.individuals.phytos, phy_model.params["P_Nsp"])
    write_species_dynamics(phy_model.t, phyts_sp, resultpath)

    for j in 1:size(phyts_sp,1)
        phyts = phyts_sp[j]
        ksd =kde(phyts[4,:])
        iksd = InterpKDE(ksd)
        size_dens[:,i,j] = pdf(iksd,xInd)
    end
    # write_output(phyts_sp, resultpath, phy_model.t)
end

# + slideshow={"slide_type": "skip"}
# post-processing model results
include("post_process.jl")
c = ["r","b","g","c","m"]
titles = ["Population (cell)", "Division (per day)", "Average Cell Size", "Grazing (per day)", "Generation",
          "Biomass (mmol/cell)", "C Reserve (mmol/cell)", "N Reserve (mmol/cell)", "P Reserve (mmol/cell)", "Average Age (hour)"]
tcks = collect(0:86400:86400*nday);
lbs = ["0","1","2","3","4","5","6","7","8","9","10"]
labels = ["Species 1","Species 2"];

# + [markdown] slideshow={"slide_type": "slide"}
# ## 4. Visualization using `PyPlot.jl`
# Similar to `matplotlib.pyplot` in Python

# + slideshow={"slide_type": "fragment"}
fig, axs = PyPlot.subplots(5,2,sharex=true,figsize=(12,10))
for i in 1:Nsp
    axs[1].plot(collect(1:60:86400*nday),rawdata[:,2,i] .* 1e11,label =labels[i],color=c[i])
    axs[2].plot(collect(1:3600:86400*nday),pops[:,i,1] ./ pops[:,i,2] .* 24,label =labels[i],color=c[i])
    axs[3].plot(collect(1:60:86400*nday),rawdata[:,5,i],label =labels[i],color=c[i])
    axs[4].plot(collect(1:3600:86400*nday),pops[:,1,3] ./ pops[:,1,2] .* 24,label =labels[i],color=c[i])
    axs[5].plot(collect(1:60:86400*nday),rawdata[:,3,i],label =labels[i],color=c[i])
    axs[6].plot(collect(1:60:86400*nday),rawdata[:,6,i] ./ 1e11,label =labels[i],color=c[i])
    axs[7].plot(collect(1:60:86400*nday),rawdata[:,7,i] ./ 1e11,label =labels[i],color=c[i])
    axs[8].plot(collect(1:60:86400*nday),rawdata[:,8,i] ./ 1e11,label =labels[i],color=c[i])
    axs[9].plot(collect(1:60:86400*nday),rawdata[:,9,i] ./ 1e11,label =labels[i],color=c[i])
    axs[10].plot(collect(1:60:86400*nday),rawdata[:,4,i],label =labels[i],color=c[i])
end
for ax in axs
    ax.tick_params("both",direction="in")
    ax.set_xticks(tcks[1:nday+1])
    ax.set_xticklabels(lbs[1:nday+1])
end
for i in 1:10
    axs[i].set_title(titles[i])
end
axs[1].legend(loc=2, fontsize=9)
axs[5].set_xlabel("Time(day)")
axs[10].set_xlabel("Time(day)")
fig.subplots_adjust(hspace=0.3)
fig.set_visible(false)
nothing

# + slideshow={"slide_type": "subslide"}
fig.set_visible(true)
fig

# + slideshow={"slide_type": "slide"}
fig,ax = plt.subplots(1,1,figsize=(8,3),sharex=true)
p = ax.pcolormesh(collect(0:60:86400*nday),collect(0.8:0.05:3),size_dens[:,:,1] ./ 20)
ax.set_ylim(0.9,3)
ax.set_ylabel("Relative Cell Size")
ax.set_xticks(tcks[1:nday+1])
ax.set_xticklabels(lbs[1:nday+1])
p1b = fig.colorbar(p,ax=ax, pad = 0.01)
p1b.set_label("Size Class Proportion")
ax.set_xlabel("Time(day)");
#plt.savefig("size_distribution.png")

# + [markdown] slideshow={"slide_type": "fragment"}
# ### for >=2 species, use code below to plot

# + slideshow={"slide_type": "fragment"}
fig,axs = plt.subplots(Nsp,1,figsize=(8,3*Nsp),sharex=true)
for i in 1:Nsp
    p = axs[i].pcolormesh(collect(0:60:86400*nday),collect(0.8:0.05:3),size_dens[:,:,i] ./ 20)
    axs[i].set_ylim(0.9,3)
    axs[i].set_ylabel("Relative Cell Size")
    axs[i].set_xticks(tcks[1:nday+1])
    axs[i].set_xticklabels(lbs[1:nday+1])
    p1b = fig.colorbar(p,ax=axs[i], pad = 0.01)
    p1b.set_label("Size Class Proportion")
end
axs[Nsp].set_xlabel("Time(day)");
#plt.savefig("size_distribution.png")

# + [markdown] slideshow={"slide_type": "slide"}
# ## Thing to do next...
# 1. Multiple species with different growth rate etc..
#     - Change parameters in `params.yml`
#     - `PCmax`: maximum photosynthesis rate
#     - `VNO3`, `VNH4`, `VP`: maximum nutrient uptake rates
#     - $\alpha$, $\Phi$: light absorption coefficient and maximum quantum yield of photosynthesis
# 2. 2-Dimensional setup with submesoscale motions
#     - Single cell example
#     - Heat induced vertical convection
#     - Random location of the individual
# -


