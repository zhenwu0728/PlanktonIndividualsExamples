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
# # 2-Dimensional single cell example
#
# - How to couple the physical model with IBM
#
# - Go through the whole processes of runing the model and visualization

# + [markdown] slideshow={"slide_type": "fragment"}
# ## For more documentation
# - <https://zhenwu0728.github.io/PlanktonIndividuals.jl/dev/> for PlanktonIndividuals.jl
# - <https://clima.github.io/Oceananigans.jl/stable/> for Oceananigans.jl
# - <http://docs.juliaplots.org/latest/> for Plots.jl
# - <https://docs.julialang.org/en/v1/> for Julia

# + [markdown] slideshow={"slide_type": "slide"}
# ## 1. Load `Julia` packages

# + slideshow={"slide_type": "-"}
using PlanktonIndividuals, Oceananigans, YAML, Plots

# + [markdown] slideshow={"slide_type": "slide"}
# ## 2. Model initialization
# ### 2.1 Physical model initialization

# + slideshow={"slide_type": "-"}
Nx = 25       # Number of grid points in x and y
Nz = 50       # Number of grid points in z
Δz = 4.0      # Grid spacing in x, y, z (meters)
Qᵀ = 1e-5     # Temperature flux at surface
∂T∂z = 0.005  # Initial vertical temperature gradient
f = 1e-4      # Coriolis parameter
α = 2e-4      # Thermal expansion coefficient
β = 8e-4      # Haline contraction coefficient

grids = RegularCartesianGrid(size=(Nx, 1, Nz), extent=(Δz*Nx, Δz*Nx, Δz*Nz))
T_bcs = TracerBoundaryConditions(grids, top = BoundaryCondition(Flux, Qᵀ), bottom = BoundaryCondition(Gradient, ∂T∂z))

model = IncompressibleModel(
                 architecture = CPU(),
                         grid = grids,
                     coriolis = FPlane(f=f),
                     buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(α=α, β=β)),
                      closure = AnisotropicMinimumDissipation(),
          boundary_conditions = (T=T_bcs,)
)

## Random noise damped at top and bottom
Ξ(z) = randn() * z / model.grid.Lz * (1 + z / model.grid.Lz) # noise

## Temperature initial condition: a stable density tradient with random noise superposed.
T₀(x, y, z) = 20 + ∂T∂z * z + ∂T∂z * model.grid.Lz * 1e-6 * Ξ(z)

set!(model, T=T₀)
wizard = TimeStepWizard(cfl=0.2, Δt=1.0, max_change=1.1, max_Δt=60.0);

# + [markdown] slideshow={"slide_type": "subslide"}
# ### 2.1 Physical model warm up

# + slideshow={"slide_type": "-"}
# warm up of physical model
simulation = Simulation(model, Δt=wizard, stop_iteration=400*40, progress_frequency=200)
run!(simulation)

Nsimulation = Simulation(model, Δt=10.0, stop_iteration=16002, progress_frequency=2)

# + [markdown] slideshow={"slide_type": "subslide"}
# ### 2.2 IBM initialization

# + slideshow={"slide_type": "-"}
resultpath = PrepRunDir()
phy_grid = read_Ogrids(grids);
nday = 2
RunParam.nTime = 1440*nday # simulated days in minute
RunParam.ΔT = 60 # model time step: 60 seconds

# update paramters
phy_params = YAML.load(open("params.yml"))
update_params!(RunParam.params,phy_params)
RunParam.params["P_Nind"] = 1
RunParam.params["Grz_P"] = 0

# Initial conditions
#           DIC  NH4   NO3   PO4   DOC  DON  DOP  POC  PON  POP  ZOO
nut_init = [2.0, 0.05, 0.10, 0.02, 1.0, 0.1, 0.1, 1.0, 0.1, 0.1, 2.0];
phy_model = PI_Model(phy_grid, RunParam;
                     nutrients = setup_nutrients(phy_grid,nut_init));

# the depth of the individual
phy_model.individuals.phytos[3,:] = rand(collect(-30.0:0.1:-2.0),1);

# + [markdown] slideshow={"slide_type": "slide"}
# ### 3. Run the model
# - **Velocities will be passed to IBM to advect individuals and nutrient tracers**
# - RK4 is used for the advection of individuals

# + jupyter={"outputs_hidden": true} slideshow={"slide_type": "-"}
anim = @animate for i in 1:1440*nday
    vel_field =[]
    for j in 1:3
        run!(Nsimulation)
        Nsimulation.stop_iteration += 2
        u = model.velocities.u.data.parent
        v = model.velocities.v.data.parent
        w = model.velocities.w.data.parent
        vel = PlanktonIndividuals.velocity(u, v, w)
        push!(vel_field,vel)
    end
    vel_itps = (generate_vel_itp(phy_model.grid, vel_field[1]),
                generate_vel_itp(phy_model.grid, vel_field[2]),
                generate_vel_itp(phy_model.grid, vel_field[3]))
    PI_advectRK4!(phy_model, RunParam.ΔT, vel_itps)
    PI_TimeStep!(phy_model, RunParam.ΔT, vel_field[end], resultpath)

    phyts_sp = sort_species(phy_model.individuals.phytos, phy_model.params["P_Nsp"])
    write_species_dynamics(phy_model.t, phyts_sp, resultpath)
    
    # ploting
    w = Array(interior(model.velocities.w))[:, 1, :]
    p1 = Plots.heatmap(phy_model.grid.xC[2:end-1], phy_model.grid.zF[2:end-1], w',xlabel="x (m)", ylabel="z (m)",clims=(-4e-2, 4e-2),color=:balance, fill=true)
    Plots.scatter!(p1, phyts_sp[1][1,:],phyts_sp[1][3,:], markersize=7,markercolor=:lime, xlims=(0,100), ylims=(-200,0), xlabel="x (m)", ylabel="z (m)", cbar=false, legend=:none)
    Plots.plot(p1, size=(300, 600),title=lpad(i÷1440,2,"0")*"day "*lpad(i÷60-24*(i÷1440),2,"0")*"hour")
end

# + slideshow={"slide_type": "slide"}
mp4(anim, "motions.mp4", fps = 60)

# + [markdown] slideshow={"slide_type": "fragment"}
# ## 4. Things to do next...
# 1. Different heat flux $Q^T$
# 2. Different location of the individual
