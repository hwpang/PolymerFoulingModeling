# save_directory = "1,3-BUTADIENE_1.0"
# import Pkg
# Pkg.activate("/home/gridsan/hwpang/Jobs/Dow")
# Pkg.instantiate()

using CSV
using DataFrames
using NonlinearSolve

using ReactionMechanismSimulator
using ReactionMechanismSimulator.SciMLBase
using ReactionMechanismSimulator.Sundials
using ReactionMechanismSimulator.YAML
using ReactionMechanismSimulator.SparseArrays
using ReactionMechanismSimulator.ForwardDiff

# read in command line arguments
film_mech_path = ARGS[1]
model_name = ARGS[2]
liquid_simulation_results_path = ARGS[3]

if model_name in ["basecase_debutanizer_model", "trace_oxygen_perturbed_debutanizer_model"]
    aspen_condition_path = ARGS[4]
    tray = parse(Int, ARGS[5])
elseif model_name == "QCMD_cell_model"
    tray = 1
    nothing
end

println("film_mech_path: $(film_mech_path)")
println("model_name: $(model_name)")
println("liquid_simulation_results_path: $(liquid_simulation_results_path)")
if model_name in ["basecase_debutanizer_model", "trace_oxygen_perturbed_debutanizer_model"]
    println("aspen_condition_path: $(aspen_condition_path)")
end

# set up paths
film_mech_directory = dirname(film_mech_path)
println("film_mech_directory: $(film_mech_directory)")
save_directory = dirname(liquid_simulation_results_path)
println("save_directory: $(save_directory)")
liquid_species_mapping_path = joinpath(film_mech_directory, "liquid_species_mapping.yml")
println("liquid_species_mapping_path: $(liquid_species_mapping_path)")
ASFWnparams_path = joinpath(save_directory, "ASFWnparams.yml")
println("ASFWnparams_path: $(ASFWnparams_path)")

# set up parameters
if model_name in ["basecase_debutanizer_model", "trace_oxygen_perturbed_debutanizer_model"]
    if model_name == "basecase_debutanizer_model"
        include_oxygen = false
        abstol = 1e-18
        reltol = 1e-6
    elseif model_name == "trace_oxygen_perturbed_debutanizer_model"
        include_oxygen = true
        abstol = 1e-20
        reltol = 1e-6
    end

    d = 2.5
    h = 0.3
    A = (d / 2)^2 * pi
    Vliq = A * h
    spacing = 0.6
    Vgas = A * (spacing - h)
    hfilm0 = 1e-5
    epsilon = 0.2 #vol% of liquid in swollen film
    rho = 900.0 #density of solid
    Vfilm0 = A * hfilm0
    Vsolidinfilm0 = Vfilm0 * (1 - epsilon)
    Vliqinfilm0 = Vfilm0 * epsilon

    tf0 = 3600 * 24 * 365 * 500
    tf = 3600 * 24 * 365

    trays = 1:40
    aspen_condition = YAML.load_file(aspen_condition_path)
    Ts = aspen_condition["T"]

elseif model_name == "QCMD_cell_model"
    include_oxygen = true

    Vreactor = 40 * 1e-9
    d = 0.55 / 2.54 / 100 #in to meter
    A = (d / 2)^2 * pi
    h = Vreactor / A
    Vliq = Vreactor
    epsilon = 0.2 #vol% of liquid in swollen film
    rho = 900.0 #density of solid
    hsolid0 = 100 * 1e-9
    Vsolidinfilm0 = A * hsolid0
    Vliqinfilm0 = Vsolidinfilm0 / (1 - epsilon) * epsilon

    tf0 = 10000.0

    trays = 1:1
    Ts = [90.0 + 273.15]

    abstol = 1e-20 # smaller tol is needed for oxygen chemistry
    reltol = 1e-8
end
mass0 = Vsolidinfilm0 * rho

println("Loading film phase submodel...")

input_file = "chem_liquid_film_phase.rms"
phaseDict = readinput(input_file)

liqspcs = phaseDict["liquid"]["Species"]
liqrxns = phaseDict["liquid"]["Reactions"]

filmspcs = phaseDict["film"]["Species"]
filmrxns = phaseDict["film"]["Reactions"]

interfacerxns = phaseDict[Set(["liquid", "film"])]["Reactions"]
solvent = phaseDict["Solvents"][1]

liq = IdealDiluteSolution(liqspcs, liqrxns, solvent; name="liquid", diffusionlimited=true);
liqspcnames = getfield.(liq.species, :name)
liqspcmws = getfield.(liq.species, :molecularweight)
liq_name_mw = Dict(zip(liqspcnames, liqspcmws))
liq_name_smiles = Dict(zip(liqspcnames, getfield.(liq.species, :smiles)))
liqrxnstrs = getrxnstr.(liq.reactions)
liqrxncomments = getfield.(liq.reactions, :comment)

film = FragmentBasedIdealFilm(filmspcs, filmrxns; name="film", diffusionlimited=false);
println("Number of film species: $(length(filmspcs))")
println("Number of film reactions: $(length(filmrxns))")
filmspcnames = getfield.(film.fragments, :name)
filmrxnstrs = getrxnstr.(film.reactions)
filmrxncomments = getfield.(film.reactions, :comment)
fragmentnames = [spc.name for spc in filmspcs if spc.isfragment]

println("Loading liquid simulation results...")
if model_name in ["basecase_debutanizer_model", "trace_oxygen_perturbed_debutanizer_model"]
    liquid_steady_state_mols = DataFrame(CSV.File(liquid_simulation_results_path))
elseif model_name == "QCMD_cell_model"
    liquid_steady_state_mols = DataFrame(CSV.File(liquid_simulation_results_path))
    liquid_steady_state_mols = liquid_steady_state_mols[[end], 2:end]
end

println("Loading ASF parameters...")
liquid_species_mapping = YAML.load_file(liquid_species_mapping_path)
ASFWnparams = YAML.load_file(ASFWnparams_path)

T = Ts[tray]
println("Estimating oligomer concentrations...")
dimer_weight = sum(Float64[liquid_steady_state_mols[tray, name] * liq_name_mw[name*"(L)"] for name in liquid_species_mapping["dimer(L)"]])
dimer_allylic_C_weight = sum(Float64[liquid_steady_state_mols[tray, name] * liq_name_mw[name*"(L)"] for name in liquid_species_mapping["dimer_allylic_C.(L)"]])
dimer_alkyl_C_weight = sum(Float64[liquid_steady_state_mols[tray, name] * liq_name_mw[name*"(L)"] for name in liquid_species_mapping["dimer_alkyl_C.(L)"]])
if include_oxygen
    dimer_COO_weight = sum(Float64[liquid_steady_state_mols[tray, name] * liq_name_mw[name*"(L)"] for name in liquid_species_mapping["dimer_COO.(L)"]])
    dimer_COOH_weight = sum(Float64[liquid_steady_state_mols[tray, name] * liq_name_mw[name*"(L)"] for name in liquid_species_mapping["dimer_COOH(L)"]])
end
monomer_Mw = sum(Float64[liquid_steady_state_mols[tray, name] * liq_name_mw[name*"(L)"] for name in liquid_species_mapping["monomer(L)"]]) / sum(Float64[liquid_steady_state_mols[tray, name] for name in liquid_species_mapping["monomer(L)"]])
oligomer_Mw = ASFWnparams["n oligomer"][tray] * monomer_Mw
oligomer_mol = dimer_weight / ASFWnparams["Wn dimer"][tray] * ASFWnparams["Wn oligomer"][tray] / (oligomer_Mw)
oligomer_allylic_C_mol = dimer_allylic_C_weight / ASFWnparams["Wn dimer"][tray] * ASFWnparams["Wn oligomer"][tray] / (oligomer_Mw - 1.008 / 1000)
oligomer_alkyl_C_mol = dimer_alkyl_C_weight / ASFWnparams["Wn dimer"][tray] * ASFWnparams["Wn oligomer"][tray] / (oligomer_Mw - 1.008 / 1000)
if include_oxygen
    oligomer_COO_mol = dimer_COO_weight / ASFWnparams["Wn dimer"][tray] * ASFWnparams["Wn oligomer"][tray] / (oligomer_Mw - 1.008 / 1000)
    oligomer_COOH_mol = dimer_COOH_weight / ASFWnparams["Wn dimer"][tray] * ASFWnparams["Wn oligomer"][tray] / (oligomer_Mw)
end

oligomer_init_mol = Dict(
    "C=C(L,oligomer)" => oligomer_mol * ASFWnparams["#db oligomer"][tray],
    "allylic_CH(L,oligomer)" => oligomer_mol * ASFWnparams["#db oligomer"][tray] * 4,
    "allylic_C.(L,oligomer)" => oligomer_allylic_C_mol,
    "alkyl_C.(L,oligomer)" => oligomer_alkyl_C_mol,
)
if include_oxygen
    if !occursin("OXYGEN_0.0", liquid_simulation_results_path)
        oligomer_init_mol["COO.(L,oligomer)"] = oligomer_COO_mol >= 0.0 ? oligomer_COO_mol : 0.0 # can be a very small negative number due to floating point error for bottom trays
        oligomer_init_mol["COOH(L,oligomer)"] = oligomer_COOH_mol >= 0.0 ? oligomer_COO_mol : 0.0 # can be a very small negative number due to floating point error for bottom trays
    else
        oligomer_init_mol["COO.(L,oligomer)"] = 0.0
        oligomer_init_mol["COOH(L,oligomer)"] = 0.0
    end
end

for (oligomer, mol) in oligomer_init_mol
    println(oligomer, " (mol/m^3) = ", mol / Vliq)
end

liqtotalmass = sum(liquid_steady_state_mols[tray, name] * liq_name_mw[name*"(L)"] for name in names(liquid_steady_state_mols))

liqinitialconds = Dict{String,Float64}()
liqinitialconds["T"] = T
liqinitialconds["V"] = Vliqinfilm0
for name in names(liquid_steady_state_mols)
    if name * "(L)" in keys(liq_name_smiles)
        smiles = liq_name_smiles[name*"(L)"]
        if occursin("O", smiles) && occursin("OXYGEN_0.0", liquid_simulation_results_path)
            conc = 0.0
        else
            conc = liquid_steady_state_mols[tray, name] / Vliq
            conc = conc >= 0.0 ? conc : 0.0 # can be a very small negative number due to floating point error for oxygen species in bottom trays
        end
        liqinitialconds[name*"(L)"] = conc * Vliqinfilm0
    end
end

for oligomer_label in keys(oligomer_init_mol)
    liqinitialconds[oligomer_label] = oligomer_init_mol[oligomer_label] / Vliq * Vliqinfilm0
end

filminitialconds = Dict{String,Float64}()
filminitialconds["T"] = T
filminitialconds["rho"] = rho
filminitialconds["mass"] = mass0
filminitialconds["A"] = A
if model_name == "QCMD_cell_model"
    filminitialconds["AH"] = 4.0 / (56.0 / 1000.0) * mass0
    filminitialconds["CDB"] = 1.0 / (56.0 / 1000.0) * mass0
    for fragment in fragmentnames
        if fragment != "AH" && fragment != "CDB"
            filminitialconds[fragment] = 0.0
        end
    end
else
    filminitialconds["AH"] = sum(Float64[liquid_steady_state_mols[tray, name] for name in liquid_species_mapping["allylic_CH(L)"]]) / liqtotalmass * mass0
    filminitialconds["AR"] = sum(Float64[liquid_steady_state_mols[tray, name] for name in liquid_species_mapping["allylic_C.(L)"]]) / liqtotalmass * mass0
    filminitialconds["KR"] = sum(Float64[liquid_steady_state_mols[tray, name] for name in liquid_species_mapping["alkyl_C.(L)"]]) / liqtotalmass * mass0
    filminitialconds["CDB"] = sum(Float64[liquid_steady_state_mols[tray, name] for name in liquid_species_mapping["C=C(L)"]]) / liqtotalmass * mass0

    if include_oxygen
        if !occursin("OXYGEN_0.0", liquid_simulation_results_path)
            filminitialconds["PR"] = sum(Float64[liquid_steady_state_mols[tray, name] > 0.0 ? liquid_steady_state_mols[tray, name] : 0.0 for name in liquid_species_mapping["COO.(L)"]]) / liqtotalmass * mass0
            filminitialconds["CP"] = sum(Float64[liquid_steady_state_mols[tray, name] > 0.0 ? liquid_steady_state_mols[tray, name] : 0.0 for name in liquid_species_mapping["COOC(L)"]]) / liqtotalmass * mass0
            filminitialconds["HP"] = sum(Float64[liquid_steady_state_mols[tray, name] > 0.0 ? liquid_steady_state_mols[tray, name] : 0.0 for name in liquid_species_mapping["COOH(L)"]]) / liqtotalmass * mass0
            filminitialconds["OR"] = sum(Float64[liquid_steady_state_mols[tray, name] > 0.0 ? liquid_steady_state_mols[tray, name] : 0.0 for name in liquid_species_mapping["CO.(L)"]]) / liqtotalmass * mass0
            filminitialconds["OH"] = 0.0
        else
            filminitialconds["PR"] = 0.0
            filminitialconds["CP"] = 0.0
            filminitialconds["HP"] = 0.0
            filminitialconds["OR"] = 0.0
            filminitialconds["OH"] = 0.0
        end
    end
end

for fragment in fragmentnames
    if fragment in keys(filminitialconds)
        println(fragment, " (mol/kg) = ", filminitialconds[fragment] / mass0)
    end
end

println("Solving film phase submodel...")

domainfilm, y0film, pfilm = FragmentBasedConstantTrhoDomain(phase=film, initialconds=filminitialconds);

domainliq, y0liq, pliq = ConstantTVDomain(phase=liq, initialconds=liqinitialconds, constantspecies=liqspcnames);

inter, pinter = FragmentBasedReactiveFilmGrowthInterfaceConstantT(domainfilm, domainliq, interfacerxns);

react, y0, p = Reactor((domainfilm, domainliq), (y0film, y0liq), (0.0, tf0), (inter,), (pfilm, pliq, pinter));

interrxnstrs = getrxnstr.(inter.reactions)
interrxncomments = getfield.(inter.reactions, :comment)

allrxnstrs = [filmrxnstrs; liqrxnstrs; interrxnstrs]
allrxncomments = [filmrxncomments; liqrxncomments; interrxncomments]

@time sol = solve(react.ode, react.recommendedsolver, abstol=abstol, reltol=reltol)

df = DataFrame(sol)
rename!(df, names(df)[domainliq.indexes[1]+1:domainliq.indexes[2]+1] .=> liqspcnames)
rename!(df, names(df)[domainfilm.indexes[1]+1:domainfilm.indexes[2]+1] .=> filmspcnames)
rename!(df, names(df)[domainfilm.indexes[3]+1] .=> "mass")

println("Initializing fragment concentrations with steady state fragment concentrations...")

for fragment in fragmentnames
    if fragment != "inert(S)"
        new_value = df[end, fragment] / df[end, "mass"] * mass0
        println(fragment, " (mol/kg): ", filminitialconds[fragment] / mass0, " ", new_value / mass0)
        if new_value > 0.0
            filminitialconds[fragment] = new_value
        else
            filminitialconds[fragment] = 0.0
        end
    end
end

# println("Initializing liquid concentrations in film with steady state liquid concentrations in film...")

# for spc in keys(liqinitialconds)
#     if spc != "T" && spc != "V"
#         new_value = df[end, spc] / df[end, "mass"] * mass0
#         println(spc, " (mol/kg): ", liqinitialconds[spc] / mass0, " ", new_value / mass0)
#         if new_value > 0.0
#             liqinitialconds[spc] = new_value
#         else
#             liqinitialconds[spc] = 0.0
#         end
#     end
# end

println("Solving film phase submodel with steady state fragment concentrations...")

num_cells = 1
dtheta = 1.0 / num_cells
Cbulk = y0 ./ Vliqinfilm0

for (key, value) in filminitialconds
    if key != "A" && key != "T" && key != "rho"
        filminitialconds[key] = value / num_cells
    end
end

for (key, value) in liqinitialconds
    if key != "T"
        liqinitialconds[key] = value / num_cells
    end
end

liqinitialconds["epsilon"] = epsilon

domainfilm, y0film, pfilm = FragmentBasedConstantTrhoDomain(phase=film, initialconds=filminitialconds)
domainliq, y0liq, pliq = ConstantTLiqFilmDomain(phase=liq, initialconds=liqinitialconds)
# domainliq, y0liq, pliq = ConstantTVDomain(phase=liq, initialconds=liqinitialconds)
inter, pinter = FragmentBasedReactiveFilmGrowthInterfaceConstantT(domainfilm, domainliq, interfacerxns)
react, y0, p = Reactor((domainfilm, domainliq), (y0film, y0liq), (0.0, tf), (inter,), (pfilm, pliq, pinter))

mu = liq.solvent.mu(T)
diffs = [x(T=T, mu=mu, P=1e8) for x in getfield.(liq.species,:diffusion)]

function f_film_growth!(dy, y, p, t, react, num_cells, diffs, dtheta, Cbulk, epsilon, A)
    dy .= 0.0
    domainfilm = react.domain[1]
    domainliq = react.domain[2]
    liq_inds = domainliq.indexes[1]:domainliq.indexes[2]
    @views hhats = y[domainfilm.indexes[3], :] / domainfilm.rho / domainfilm.A / (1 - epsilon)

    for j in 1:num_cells

        # reaction terms
        @views react.ode.f(dy[:, j], y[:, j], p, t)

        # diffusion terms and boundary conditions

        if j == num_cells
            # no flux condition at tray surface
            Jjp1half = 0.0
        else
            @views Jjp1half = -diffs .* (y[liq_inds, j+1] / y[domainliq.indexes[3], j+1] .- y[liq_inds, j] / y[domainliq.indexes[3], j]) ./ (hhats[j+1] / 2 + hhats[j] / 2)
        end
        if j == 1
            # C = C bulk in liquid
            @views Jjm1half = -diffs .* (y[liq_inds, j] / y[domainliq.indexes[3], j] .- Cbulk[liq_inds]) ./ (hhats[j] / 2 + hhats[j] / 2)
        else
            @views Jjm1half = -diffs .* (y[liq_inds, j] / y[domainliq.indexes[3], j] .- y[liq_inds, j-1] / y[domainliq.indexes[3], j-1]) ./ (hhats[j] / 2 + hhats[j-1] / 2)
        end
        @views dy[liq_inds, j] .+= -(Jjp1half .- Jjm1half) .* A # normalized x by film thickness

        # change due to moving control surface
        if j == 1
            @views dy[liq_inds, j] .+= dy[domainliq.indexes[3], j] .* Cbulk[liq_inds]
        else
            rates = dy[domainliq.indexes[3], j] .* (y[liq_inds, j-1] / y[domainliq.indexes[3], j-1])
            @views dy[liq_inds, j] .+= rates
            @views dy[liq_inds, j-1] .-= rates
        end
    end
end

function f!(dy, y, p, t)
    f_film_growth!(unflatten(dy), unflatten(y), p, t, react, num_cells, diffs, dtheta, Cbulk, epsilon, A)
end

function jacobianyforwarddiff!(J, y, p, t)
    ForwardDiff.jacobian!(J,(dy, y) -> f!(dy, y, p, t), zeros(length(y)), y)
end

function unflatten(y::AbstractVector)
    return reshape(y, num_variables, :)
end

function rops(ssys, t, cell_ind)
    domains = getfield.(ssys.sims, :domain)
    Nrxns = sum([length(sim.domain.phase.reactions) for sim in ssys.sims]) + sum([length(inter.reactions) for inter in ssys.interfaces if hasproperty(inter, :reactions)])
    Nspcs = sum([length(getphasespecies(sim.domain.phase)) for sim in ssys.sims])
    cstot = zeros(Nspcs)
    vns = Array{Any,1}(undef, length(domains))
    vcs = Array{Any,1}(undef, length(domains))
    vT = Array{Any,1}(undef, length(domains))
    vP = Array{Any,1}(undef, length(domains))
    vV = Array{Any,1}(undef, length(domains))
    vC = Array{Any,1}(undef, length(domains))
    vN = Array{Any,1}(undef, length(domains))
    vmu = Array{Any,1}(undef, length(domains))
    vkfs = Array{Any,1}(undef, length(domains))
    vkrevs = Array{Any,1}(undef, length(domains))
    vHs = Array{Any,1}(undef, length(domains))
    vUs = Array{Any,1}(undef, length(domains))
    vGs = Array{Any,1}(undef, length(domains))
    vdiffs = Array{Any,1}(undef, length(domains))
    vCvave = Array{Any,1}(undef, length(domains))
    vphi = Array{Any,1}(undef, length(domains))
    ropmat = spzeros(Nrxns, Nspcs)
    start = 1
    for (k, sim) in enumerate(ssys.sims)
        vns[k], vcs[k], vT[k], vP[k], vV[k], vC[k], vN[k], vmu[k], vkfs[k], vkrevs[k], vHs[k], vUs[k], vGs[k], vdiffs[k], vCvave[k], vphi[k] = calcthermo(sim.domain, unflatten(ssys.sol(t))[:, cell_ind], t)
        cstot[sim.domain.indexes[1]:sim.domain.indexes[2]] = vcs[k]
        rops!(ropmat, sim.domain.rxnarray, cstot, vkfs[k], vkrevs[k], vV[k], start)
        start += length(vkfs[k])
    end
    for inter in ssys.interfaces
        if inter isa FragmentBasedReactiveFilmGrowthInterfaceConstantT
            kfs, krevs = getkfskrevs(inter)
            rops!(ropmat, inter.rxnarray, inter.fragmentbasedrxnarray, cstot, kfs, krevs, vV[inter.domaininds[1]], start)
            start += length(kfs)
        elseif hasproperty(inter, :reactions)
            kfs, krevs = getkfskrevs(inter, vT[inter.domaininds[1]], vT[inter.domaininds[2]], vphi[inter.domaininds[1]], vphi[inter.domaininds[2]], vGs[inter.domaininds[1]], vGs[inter.domaininds[2]], cstot)
            rops!(ropmat, inter.rxnarray, cstot, kfs, krevs, inter.A, start)
            start += length(kfs)
        end
    end
    return ropmat
end


num_variables = length(y0)
y0z = zeros(num_variables, num_cells)

for j in 1:num_cells
    y0z[:, j] .= y0
end

y0z = reshape(y0z, :)

dy0z = zeros(num_variables * num_cells)
f!(dy0z, y0z, p, 0.0)
println("y0z")
display(unflatten(y0z))
println("dy0z")
display(unflatten(dy0z))

odefcn = ODEFunction(f!; jac=jacobianyforwarddiff!)
odeprob = ODEProblem(odefcn, y0z, (0.0, 3600 * 24 * 365), p)

@time sol = solve(odeprob, react.recommendedsolver, abstol=1e-18, reltol=1e-6)

println("Saving results...")
df = DataFrame(sol)
cell_names = Array{String,2}(undef, num_variables, num_cells)
for j in 1:num_cells
    cell_names[domainliq.indexes[1]:domainliq.indexes[2], j] .= liqspcnames .* "_cell_" .* string(j)
    cell_names[domainliq.indexes[3], j] = "Vliqinfilm_cell_" .* string(j)
    cell_names[domainfilm.indexes[1]:domainfilm.indexes[2], j] .= filmspcnames .* "_cell_" .* string(j)
    cell_names[domainfilm.indexes[3], j] = "mass_cell_" .* string(j)
end

cell_names = reshape(cell_names, :)
rename!(df, names(df)[2:end] .=> cell_names)

path = "$(save_directory)/simulation_film_1D_$(tray)_$(num_cells)cells.csv"
println("Saving results to $(path)")
CSV.write(path, df)

# println("Saving rops...")

# function save_rop(sol)

#     ssys = SystemSimulation(sol, (domainfilm, domainliq,), (inter,), p)
#     spcnames = ssys.names
#     spcmassnames = spcnames[:]
#     push!(spcmassnames, "Vliqinfilm")
#     push!(spcmassnames, "mass")

#     for cell_ind in 1:num_cells

#         count = 0

#         while true

#             ropmat = rops(ssys, sol.t[end]/2, cell_ind)

#             rxnind, spcind, rop = findnz(ropmat)
#             no_nan = any((isnan).(rop)) == false
#             println("No NaN ", no_nan)
#             if model_name == "basecase_debutanizer_model"
#                 no_large = all((abs).(rop) .< 1) == true
#             elseif model_name == "trace_oxygen_perturbed_debutanizer_model"
#                 no_large = true
#             end
#             println("No large ", no_large)

#             if no_nan && no_large

#                 rxnind = trunc.(Int, rxnind)
#                 spcind = trunc.(Int, spcind)
#                 rop_rxnstrs = allrxnstrs[rxnind]
#                 rop_rxncomments = allrxncomments[rxnind]
#                 rop_spcnames = spcmassnames[spcind]

#                 df_rop = DataFrame()
#                 df_rop.rop_rxnstr = rop_rxnstrs
#                 df_rop.rop_rxncomment = rop_rxncomments
#                 df_rop.rop_spcname = rop_spcnames
#                 df_rop.rop = rop

#                 ####
#                 CSV.write("$(save_directory)/simulation_film_rop_1D_$(tray).csv", df_rop, quotestrings=true)

#                 rops, rop_rxnstrs, rop_rxncomments = get_rops(df_rop, "mass")
#                 no_neg_rop = all((rops .>= 0) .|| occursin.("<=>", rop_rxncomments)) # positive rop except for reversible reactions
#                 println("No negative rop ", no_neg_rop)

#                 if no_neg_rop

#                     CSV.write("$(save_directory)/simulation_film_rop_$(tray)_cell$(cell)_$(num_cells)cells.csv", df_rop, quotestrings=true)

#                     break

#                 end
#             end

#             count += 1

#             if count > 10
#                 break
#             end
#         end
#     end
# end

# save_rop(sol)

println("Done!")