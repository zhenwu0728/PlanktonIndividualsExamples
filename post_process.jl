# Read in model results from text file
Ndata = readdlm("results/cons_N.txt")
Pdata = readdlm("results/cons_P.txt")
rawdata = zeros(RunParam.nTime,10,5)
for i in 1:Nsp
    rawdata[:,:,i] = readdlm("results/dynamic_species"*lpad(i,3,"0")*".txt")
end

# Reorganize model diagnostics in a 2D array (depth,time)
nT_diags = Int(RunParam.nTime*RunParam.ΔT/3600)
CNP = zeros(nT_diags,5,4)     # cellular CNP concentrations
counts = zeros(nT_diags,5)    # cell counts during each diagnostic period
rates = zeros(nT_diags,5,2)   # photosynthesis rate, biosynthesis rates, nutrient uptake rates etc.
pops = zeros(nT_diags,5,3)    # population dynamics: division, grazing etc.
diags_spcs = phy_model.diags.spcs
diags_pop  = phy_model.diags.pop
diags_tr   = phy_model.diags.tr
for i in 1:Nsp
    CNP[:,i,1] = diags_spcs[1,1,1,:,i,10] .+ diags_spcs[1,1,1,:,i,11]
    CNP[:,i,2] = diags_spcs[1,1,1,:,i,10] ./ 106 .* 16 .+ diags_spcs[1,1,1,:,i,12]
    CNP[:,i,3] = diags_spcs[1,1,1,:,i,10] ./ 106 .+ diags_spcs[1,1,1,:,i,13]
    CNP[:,i,4] = diags_spcs[1,1,1,:,i,14]
    counts[:,i]= diags_spcs[1,1,1,:,i,1]
    rates[:,i,1]=diags_spcs[1,1,1,:,i,2]
    rates[:,i,2]=diags_spcs[1,1,1,:,i,7]
    pops[:,i,1]= diags_pop[1,1,1,:,i,1]
    pops[:,i,2]= diags_tr[1,1,1,:,i+1] ./ 60
    pops[:,i,3]= diags_pop[1,1,1,:,i,3]
end