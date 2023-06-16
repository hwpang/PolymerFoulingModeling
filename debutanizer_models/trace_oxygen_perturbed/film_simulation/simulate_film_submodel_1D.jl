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
rms_mech_path = ARGS[1]
model_name = ARGS[2]
liquid_simulation_results_path = ARGS[3]

if model_name in ["basecase_debutanizer_model", "trace_oxygen_perturbed_debutanizer_model"]
    aspen_condition_path = ARGS[4]
    tray = parse(Int, ARGS[5])
elseif model_name == "QCMD_cell_model"
    tray = 1
    nothing
end

println("rms_mech_path: $(rms_mech_path)")
println("model_name: $(model_name)")
println("liquid_simulation_results_path: $(liquid_simulation_results_path)")
if model_name in ["basecase_debutanizer_model", "trace_oxygen_perturbed_debutanizer_model"]
    println("aspen_condition_path: $(aspen_condition_path)")
end

# set up paths
rms_mech_directory = dirname(rms_mech_path)
println("rms_mech_directory: $(rms_mech_directory)")
save_directory = dirname(liquid_simulation_results_path)
println("save_directory: $(save_directory)")
liquid_species_mapping_path = joinpath(rms_mech_directory, "liquid_species_mapping.yml")
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

input_file = joinpath(rms_mech_directory, "chem_film_phase.rms")
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

domainliq, y0liq, pliq = ConstantTVDomain(phase=liq, initialconds=liqinitialconds);

inter, pinter = FragmentBasedReactiveFilmGrowthInterfaceConstantT(domainfilm, domainliq, interfacerxns);

react, y0, p = Reactor((domainfilm, domainliq), (y0film, y0liq), (0.0, tf0), (inter,), (pfilm, pliq, pinter));

interrxnstrs = getrxnstr.(inter.reactions)
interrxncomments = getfield.(inter.reactions, :comment)

allrxnstrs = [filmrxnstrs; liqrxnstrs; interrxnstrs]
allrxncomments = [filmrxncomments; liqrxncomments; interrxncomments]

# @time sol = solve(react.ode, react.recommendedsolver, abstol=abstol, reltol=reltol);

num_cells = 5
dtheta = 1.0 / num_cells
mu = liq.solvent.mu(T)
diffs = [x(T=T, mu=mu, P=1e8) for x in getfield.(liq.species,:diffusion)]

function f_film_growth!(dy, y, p, t, react, num_cells, diffs, dtheta, rho, A)
    dy .= 0.0
    liq_inds = react.domain[2].indexes[1]:react.domain[2].indexes[2]
    @views h = (sum(y[end, :]) / rho / A)

    for j in 1:num_cells

        # reaction terms
        @views react.ode.f(dy[:, j], y[:, j], p, t)

        # diffusion terms and boundary conditions
        if j == 1
            # C = Cbulk at boundary
            nothing
        else
            if j == num_cells
                # no flux condition at tray surface
                Jjp1half = 0.0
            else
                Jjp1half = - diffs .* (y[liq_inds, j+1] .- y[liq_inds, j]) ./ dtheta
            end
            Jjm1half = - diffs .* (y[liq_inds, j] .- y[liq_inds, j-1]) ./ dtheta
            @views dy[liq_inds, j] .+= (Jjp1half .- Jjm1half) ./ dtheta .* h^2 # normalized x by film thickness
        end
    end
end

function f!(dy, y, p, t)
    f_film_growth!(unflatten(dy), unflatten(y), p, t, react, num_cells, diffs, dtheta, rho, A)
end

function jacobianyforwarddiff!(J, y, p, t)
    ForwardDiff.jacobian!(J,(dy, y) -> f!(dy, y, p, t), zeros(length(y)), y)
end

function unflatten(y::AbstractVector)
    return reshape(y, num_variables, :)
end

num_variables = length(y0)
y0z = zeros(num_variables, num_cells)

for j in 1:num_cells
    y0z[:, j] .= y0 ./ num_cells
end

y0z = reshape(y0z, :)

dy0z = zeros(num_variables * num_cells)
f!(dy0z, y0z, p, 0.0)
println("y0z")
display(unflatten(y0z))
println("dy0z")
display(unflatten(dy0z))

odefcn = ODEFunction(f!; jac=jacobianyforwarddiff!)
odeprob = ODEProblem(odefcn, y0z, (0.0, 1.0), p)

@time sol = solve(odeprob, react.recommendedsolver, abstol=1e-18, reltol=1e-6)

# save results

df = DataFrame(sol)
cell_names = Array{String,2}(undef, num_variables, num_cells)
for j in 1:num_cells
    cell_names[domainliq.indexes[1]+1:domainliq.indexes[2]+1, j] .= liqspcnames .* "_cell_" .* string(j)
    cell_names[domainfilm.indexes[1]+1:domainfilm.indexes[2]+1, j] .= filmspcnames .* "_cell_" .* string(j)
    cell_names[domainfilm.indexes[3]+1, j] .= "mass_cell_" .* string(j)
end
rename!(df, names(df) .=> cell_names)
CSV.write("$(save_directory)/simulation_film_1D_$(tray).csv", df)
