using CSV
using DataFrames
using NonlinearSolve

using ReactionMechanismSimulator
using ReactionMechanismSimulator.SciMLBase
using ReactionMechanismSimulator.Sundials
using ReactionMechanismSimulator.YAML
using ReactionMechanismSimulator.SparseArrays
using ReactionMechanismSimulator.ForwardDiff
import ReactionMechanismSimulator: getphasespecies, rops, rops!, getkfskrevs

# read in command line arguments
film_mech_directory = ARGS[1]
model_name = ARGS[2]
liquid_simulation_results_path = ARGS[3]

if model_name in ["basecase_debutanizer_model", "trace_oxygen_perturbed_debutanizer_model"]
    aspen_condition_path = ARGS[4]
    tray = parse(Int, ARGS[5])
    num_cells = parse(Int, ARGS[6])
    perturb_factor_string = ARGS[7]
    method = ARGS[8]
elseif model_name == "QCMD_cell_model"
    tray = 1
    nothing
end



println("film_mech_directory: $(film_mech_directory)")
println("model_name: $(model_name)")
println("liquid_simulation_results_path: $(liquid_simulation_results_path)")
if model_name in ["basecase_debutanizer_model", "trace_oxygen_perturbed_debutanizer_model"]
    println("aspen_condition_path: $(aspen_condition_path)")
end
println("tray: $(tray)")
println("perturb_factor_string: $(perturb_factor_string)")
println("method: $(method)")

# set up paths
# save_directory = dirname(liquid_simulation_results_path)
save_directory = "discretization_studies"
println("save_directory: $(save_directory)")
liquid_species_mapping_path = joinpath(film_mech_directory, "liquid_species_mapping.yml")
println("liquid_species_mapping_path: $(liquid_species_mapping_path)")
liquid_results_directory = dirname(liquid_simulation_results_path)
ASFWnparams_path = joinpath(liquid_results_directory, "ASFWnparams.yml")
println("ASFWnparams_path: $(ASFWnparams_path)")

# set up parameters
if model_name in ["basecase_debutanizer_model", "trace_oxygen_perturbed_debutanizer_model"]
    if model_name == "basecase_debutanizer_model"
        include_oxygen = false
        tf0 = 3600 * 24 * 365 * 50
        abstol = 1e-18
        reltol = 1e-6
    elseif model_name == "trace_oxygen_perturbed_debutanizer_model"
        include_oxygen = true
        tf0 = 3600 * 24 * 365 * 0.5
        abstol = 1e-18
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

input_file = joinpath(film_mech_directory, "chem_liquid_film_phase.rms")
phaseDict = readinput(input_file)

liqspcs = phaseDict["liquid"]["Species"]
liqrxns = phaseDict["liquid"]["Reactions"]

filmspcs = phaseDict["film"]["Species"]
filmrxns = phaseDict["film"]["Reactions"]

interfacerxns = phaseDict[Set(["liquid", "film"])]["Reactions"]
solvent = phaseDict["Solvents"][1]

liq_no_rxns = IdealDiluteSolution(liqspcs, [], solvent; name="liquid", diffusionlimited=true);
liq_with_rxns = IdealDiluteSolution(liqspcs, liqrxns, solvent; name="liquid", diffusionlimited=true);
liqspcnames = getfield.(liq_with_rxns.species, :name)
liqspcmws = getfield.(liq_with_rxns.species, :molecularweight)
liq_name_mw = Dict(zip(liqspcnames, liqspcmws))
liq_name_smiles = Dict(zip(liqspcnames, getfield.(liq_with_rxns.species, :smiles)))
liqrxnstrs = getrxnstr.(liq_with_rxns.reactions)
liqrxncomments = getfield.(liq_with_rxns.reactions, :comment)

film = FragmentBasedIdealFilm(filmspcs, filmrxns; name="film", diffusionlimited=false);
println("Number of film species: $(length(filmspcs))")
println("Number of film reactions: $(length(filmrxns))")
filmspcnames = getfield.(film.fragments, :name)
filmrxnstrs = getrxnstr.(film.reactions)
filmrxncomments = getfield.(film.reactions, :comment)
fragmentnames = [spc.name for spc in filmspcs if spc.isfragment && spc.name != "CD"]

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

domainliq, y0liq, pliq = ConstantTVDomain(phase=liq_no_rxns, initialconds=liqinitialconds, constantspecies=liqspcnames);

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

println("Solving film phase submodel with steady state fragment concentrations...")

function sparsity_pattern(dy, y, p, t, react, num_cells, diffs, Cbulk, epsilon, A)
    liq_inds = domainliq.indexes[1]:domainliq.indexes[2]
    liq_volume_ind = domainliq.indexes[3]
    film_inds = domainfilm.indexes[1]:domainfilm.indexes[2]
    mass_ind = domainfilm.indexes[3]

    J = spzeros(length(y), length(y)length(y))

    @views hhats = y[mass_ind, :] / domainfilm.rho / domainfilm.A / (1 - epsilon)

    for j in 1:num_cells

        # reaction terms
        @views react.ode.f(dy[:, j], y[:, j], p, t)
        J[:, j, :] .= 1

        # diffusion terms and boundary conditions

        if j == num_cells
            # no flux condition at tray surface
            Jjp1half = 0.0
        else
            @views Jjp1half = -diffs .* (y[liq_inds, j+1] / y[liq_volume_ind, j+1] .- y[liq_inds, j] / y[liq_volume_ind, j]) ./ (hhats[j+1] / 2 + hhats[j] / 2)
            J[liq_inds, j, :] .= 1
            J[liq_inds, j+1, :] .= 1
            J[liq_volume_ind, j, :] .= 1
            J[liq_volume_ind, j+1, :] .= 1
            J[mass_ind, j, :] .= 1
            J[mass_ind, j+1, :] .= 1
        end
        if j == 1
            # C = C bulk in liquid
            @views Jjm1half = -diffs .* (y[liq_inds, j] / y[liq_volume_ind, j] .- Cbulk[liq_inds]) ./ (hhats[j] / 2 + hhats[j] / 2)
            J[liq_inds, j, :] .= 1
        else
            @views Jjm1half = -diffs .* (y[liq_inds, j] / y[liq_volume_ind, j] .- y[liq_inds, j-1] / y[liq_volume_ind, j-1]) ./ (hhats[j] / 2 + hhats[j-1] / 2)
            J[liq_inds, j, :] .= 1
            J[liq_inds, j-1, :] .= 1
        end
        @views dy[liq_inds, j] .+= -(Jjp1half .- Jjm1half) .* A


        upper_surface_volume_change = sum(dy[mass_ind, j:num_cells]) / domainfilm.rho / domainfilm.A / (1 - domainliq.epsilon)
        lower_surface_volume_change = sum(dy[mass_ind, j+1:num_cells]) / domainfilm.rho / domainfilm.A / (1 - domainliq.epsilon)
        if upper_surface_volume_change >= 0.0
            if j == 1
                @views dy[liq_inds, j] .+= upper_surface_volume_change * domainliq.epsilon * A .* Cbulk[liq_inds]
                J[liq_inds, j:num_cells, mass_ind] .= 1
            else
                @views dy[liq_inds, j] .+= upper_surface_volume_change * domainliq.epsilon * A .* (y[liq_inds, j-1] / y[liq_volume_ind, j-1])
                J[liq_inds, j-1, :] .= 1
            end
        else
            @views dy[liq_inds, j] .-= abs(upper_surface_volume_change) * domainliq.epsilon * A .* (y[liq_inds, j] / y[liq_volume_ind, j])
            J[liq_inds, j, :] .= 1
        end
        if lower_surface_volume_change >= 0.0
            @views dy[liq_inds, j] .-= lower_surface_volume_change * domainliq.epsilon * A .* (y[liq_inds, j] / y[liq_volume_ind, j])
            J[liq_inds, j, :] .= 1
        else
            if j != num_cells
                @views dy[liq_inds, j] .+= abs(lower_surface_volume_change) * domainliq.epsilon * A .* (y[liq_inds, j+1] / y[liq_volume_ind, j+1])
                J[liq_inds, j, :] .= 1
            end
        end
    end
end
function sparsity_2!(dy, y, p, t, react, num_cells, diffs, Cbulk, epsilon, A)
    dy .= 0.0
    liq_inds = domainliq.indexes[1]:domainliq.indexes[2]
    liq_volume_ind = domainliq.indexes[3]
    film_inds = domainfilm.indexes[1]:domainfilm.indexes[2]
    mass_ind = domainfilm.indexes[3]

    @views hhats = y[mass_ind, :]

    for j in 1:num_cells

        # reaction terms
        dy[:, j] .+= y[:, j]

        # diffusion terms and boundary conditions

        if j == num_cells
            # no flux condition at tray surface
            Jjp1half = 0.0
        else
            @views Jjp1half = y[liq_inds, j+1] .+ y[domainliq.indexes[3], j+1] .+ y[liq_inds, j] + y[domainliq.indexes[3], j] .+ hhats[j+1] + hhats[j]
        end
        if j == 1
            # C = C bulk in liquid
            @views Jjm1half = y[liq_inds, j] .+ y[domainliq.indexes[3], j] .+ Cbulk[liq_inds] .+ hhats[j] .+ hhats[j]
        else
            @views Jjm1half = y[liq_inds, j] .+ y[domainliq.indexes[3], j] .+ y[liq_inds, j-1] / y[domainliq.indexes[3], j-1] .+ (hhats[j] + hhats[j-1])
        end
        @views dy[liq_inds, j] .+= Jjp1half .+ Jjm1half

        upper_surface_volume_change = sum(dy[mass_ind, j:num_cells]) / domainfilm.rho / domainfilm.A / (1 - domainliq.epsilon)
        lower_surface_volume_change = sum(dy[mass_ind, j+1:num_cells]) / domainfilm.rho / domainfilm.A / (1 - domainliq.epsilon)
        if upper_surface_volume_change >= 0.0
            if j == 1
                @views dy[liq_inds, j] .+= upper_surface_volume_change * domainliq.epsilon * A .* Cbulk[liq_inds]
            else
                @views dy[liq_inds, j] .+= upper_surface_volume_change * domainliq.epsilon * A .* (y[liq_inds, j-1] / y[domainliq.indexes[3], j-1])
            end
        else
            @views dy[liq_inds, j] .-= abs(upper_surface_volume_change) * domainliq.epsilon * A .* (y[liq_inds, j] / y[domainliq.indexes[3], j])
        end
        if lower_surface_volume_change >= 0.0
            @views dy[liq_inds, j] .-= lower_surface_volume_change * domainliq.epsilon * A .* (y[liq_inds, j] / y[domainliq.indexes[3], j])
        else
            if j != num_cells
                @views dy[liq_inds, j] .+= abs(lower_surface_volume_change) * domainliq.epsilon * A .* (y[liq_inds, j+1] / y[domainliq.indexes[3], j+1])
            end
        end
    end
end
function f_film_growth!(dy, y, p, t, react, num_cells, diffs, Cbulk, epsilon, A)
    dy .= 0.0
    liq_inds = domainliq.indexes[1]:domainliq.indexes[2]
    liq_volume_ind = domainliq.indexes[3]
    film_inds = domainfilm.indexes[1]:domainfilm.indexes[2]
    mass_ind = domainfilm.indexes[3]

    @views hhats = y[mass_ind, :] / domainfilm.rho / domainfilm.A / (1 - epsilon)

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
        @views dy[liq_inds, j] .+= -(Jjp1half .- Jjm1half) .* A

        upper_surface_volume_change = sum(dy[mass_ind, j:num_cells]) / domainfilm.rho / domainfilm.A / (1 - domainliq.epsilon)
        lower_surface_volume_change = sum(dy[mass_ind, j+1:num_cells]) / domainfilm.rho / domainfilm.A / (1 - domainliq.epsilon)
        if upper_surface_volume_change >= 0.0
            if j == 1
                @views dy[liq_inds, j] .+= upper_surface_volume_change * domainliq.epsilon * A .* Cbulk[liq_inds]
            else
                @views dy[liq_inds, j] .+= upper_surface_volume_change * domainliq.epsilon * A .* (y[liq_inds, j-1] / y[domainliq.indexes[3], j-1])
            end
        else
            @views dy[liq_inds, j] .-= abs(upper_surface_volume_change) * domainliq.epsilon * A .* (y[liq_inds, j] / y[domainliq.indexes[3], j])
        end
        if lower_surface_volume_change >= 0.0
            @views dy[liq_inds, j] .-= lower_surface_volume_change * domainliq.epsilon * A .* (y[liq_inds, j] / y[domainliq.indexes[3], j])
        elseif j != num_cells
            @views dy[liq_inds, j] .+= abs(lower_surface_volume_change) * domainliq.epsilon * A .* (y[liq_inds, j+1] / y[domainliq.indexes[3], j+1])
        end
    end
end

dy0 = zeros(num_variables * num_cells)
function jacobianyforwarddiff!(J, y, p, t)
    ForwardDiff.jacobian!(J, (dy, y) -> f!(dy, y, p, t), dy0, y)
end

import Base: size, getindex, setindex!

struct MyMatrix{T} <: AbstractMatrix{T}
    vector::Array{T}
    num_variables::Int
    num_cells::Int
end

size(A::MyMatrix{T}) where {T} = (A.num_variables, A.num_cells)

function getindex(A::MyMatrix{T}, i::Int, j::Int) where {T}
    return A.vector[i+(j-1)*A.num_variables]
end

function setindex!(A::MyMatrix{T}, value::T, i::Int, j::Int) where {T}
    A.vector[i+(j-1)*A.num_variables] = value
end

struct Unflatten
    num_variables::Int64
    num_cells::Int64
end

(unflatten::Unflatten)(y::T) where {T<:AbstractArray} = MyMatrix(y, unflatten.num_variables, unflatten.num_cells)

unflatten = Unflatten(num_variables, num_cells)

function rops_1D(ssys, t, cell_ind)
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
    Nthermos = sum([length(sim.domain.thermovariabledict) for sim in ssys.sims])
    if any([domain isa FragmentBasedConstantTrhoDomain for domain in domains])
        ropmat = spzeros(Nrxns, Nspcs + Nthermos)
    else
        ropmat = spzeros(Nrxns, Nspcs)
    end
    start = 0
    for (k, sim) in enumerate(ssys.sims)
        vns[k], vcs[k], vT[k], vP[k], vV[k], vC[k], vN[k], vmu[k], vkfs[k], vkrevs[k], vHs[k], vUs[k], vGs[k], vdiffs[k], vCvave[k], vphi[k] = calcthermo(sim.domain, unflatten(ssys.sol(t))[:, cell_ind], t)
        cstot[sim.domain.indexes[1]:sim.domain.indexes[2]] = vcs[k]
        rops!(ropmat, sim.domain.rxnarray, cstot, vkfs[k], vkrevs[k], vV[k], start)
        start += length(vkfs[k])
    end
    for inter in ssys.interfaces
        if inter isa FragmentBasedReactiveFilmGrowthInterfaceConstantT
            kfs, krevs = getkfskrevs(inter)
            rops!(ropmat, inter.rxnarray, inter.fragmentbasedrxnarray, cstot, kfs, krevs, vV[inter.domaininds[1]], inter.Mws, inter.domainfilm.indexes[1]:inter.domainfilm.indexes[2], inter.domainfilm.indexes[3], start)
            start += length(kfs)
        elseif hasproperty(inter, :reactions)
            kfs, krevs = getkfskrevs(inter, vT[inter.domaininds[1]], vT[inter.domaininds[2]], vphi[inter.domaininds[1]], vphi[inter.domaininds[2]], vGs[inter.domaininds[1]], vGs[inter.domaininds[2]], cstot)
            rops!(ropmat, inter.rxnarray, cstot, kfs, krevs, inter.A, start)
            start += length(kfs)
        end
    end
    return ropmat
end

function get_rops(df, rop_name; loss_only=false, production_only=false, N=5)
    name_inds = (df[!, "rop_spcname"] .== rop_name)
    rop_rxnstrs = df[name_inds, "rop_rxnstr"]
    rop_rxncomments = df[name_inds, "rop_rxncomment"]
    rops_selected = df[name_inds, "rop"]
    if loss_only
        loss_inds = (rops_selected .< 0)
        rops_selected = abs.(rops_selected[loss_inds])
        rop_rxnstrs = rop_rxnstrs[loss_inds]
        rop_rxncomments = rop_rxncomments[loss_inds]
    elseif production_only
        prod_inds = (rops_selected .> 0)
        rops_selected = rops_selected[prod_inds]
        rop_rxnstrs = rop_rxnstrs[prod_inds]
        rop_rxncomments = rop_rxncomments[prod_inds]
    end
    sorted_inds = sortperm(abs.(rops_selected), rev=true)[1:N]
    return rops_selected[sorted_inds], rop_rxnstrs[sorted_inds], rop_rxncomments[sorted_inds]
end

function save_rop(sol; suffix="")

    ssys = SystemSimulation(sol, (domainfilm, domainliq,), (inter,), p)
    spcnames = ssys.names
    spcmassnames = spcnames[:]
    push!(spcmassnames, "Vliqinfilm")
    push!(spcmassnames, "mass")

    for cell_ind in 1:num_cells

        count = 0

        while true

            ropmat = rops_1D(ssys, sol.t[end-1], cell_ind)

            rxnind, spcind, rop = findnz(ropmat)
            no_nan = any((isnan).(rop)) == false
            println("No NaN ", no_nan)
            if model_name == "basecase_debutanizer_model"
                no_large = all((abs).(rop) .< 1) == true
            elseif model_name == "trace_oxygen_perturbed_debutanizer_model"
                no_large = true
            end
            println("No large ", no_large)

            if no_nan && no_large

                rxnind = trunc.(Int, rxnind)
                spcind = trunc.(Int, spcind)
                rop_rxnstrs = allrxnstrs[rxnind]
                rop_rxncomments = allrxncomments[rxnind]
                rop_spcnames = spcmassnames[spcind]

                df_rop = DataFrame()
                df_rop.rop_rxnstr = rop_rxnstrs
                df_rop.rop_rxncomment = rop_rxncomments
                df_rop.rop_spcname = rop_spcnames
                df_rop.rop = rop

                ####

                rops, rop_rxnstrs, rop_rxncomments = get_rops(df_rop, "mass")
                no_neg_rop = all((rops .>= 0) .|| occursin.("<=>", rop_rxncomments)) # positive rop except for reversible reactions
                println("No negative rop ", no_neg_rop)

                if no_neg_rop

                    CSV.write("$(save_directory)/simulation_film_rop_$(tray)_$(num_cells)cells_cell$(cell_ind)$(suffix).csv", df_rop, quotestrings=true)

                    break

                end
            end

            count += 1

            if count > 10
                break
            end
        end
    end
end

if method == "same_initial_size"

    print("Using same initial size...")

    cell_size_fractions = ones(num_cells)
    cell_size_fractions ./= sum(cell_size_fractions)
    Cbulk = y0 ./ Vliqinfilm0

    liqinitialconds["epsilon"] = epsilon

    domainfilm, y0film, pfilm = FragmentBasedConstantTrhoDomain(phase=film, initialconds=filminitialconds)
    domainliq, y0liq, pliq = ConstantTLiqFilmDomain(phase=liq_with_rxns, initialconds=liqinitialconds)
    inter, pinter = FragmentBasedReactiveFilmGrowthInterfaceConstantT(domainfilm, domainliq, interfacerxns)
    react, y0, p = Reactor((domainfilm, domainliq), (y0film, y0liq), (0.0, tf), (inter,), (pfilm, pliq, pinter))

    mu = liq_with_rxns.solvent.mu(T)
    diffs = [x(T=T, mu=mu, P=1e8) for x in getfield.(liq_with_rxns.species, :diffusion)]

    num_variables = length(y0)
    y0z = zeros(num_variables, num_cells)

    for j in 1:num_cells
        y0z[:, j] .= y0 * cell_size_fractions[j]
    end

    y0z = reshape(y0z, :)
    struct model_info{T1,T2,T3,T4,T5,T6,T7}
        react::T1
        num_cells::T2
        diffs::T3
        Cbulk::T4
        epsilon::T5
        A::T6
        unflatten::T7
    end
    function (model::model_info)(dy, y, p, t)
        (; react, num_cells, diffs, Cbulk, epsilon, A, unflatten) = model
        f_film_growth!(unflatten(dy), unflatten(y), p, t, react, num_cells, diffs, Cbulk, epsilon, A)
    end

    f! = model_info(react, num_cells, diffs, Cbulk, epsilon, A, unflatten)

    odefcn = ODEFunction(f!; jac=jacobianyforwarddiff!)
    odeprob = ODEProblem(odefcn, y0z, (0.0, tf), p)

    @time sol = solve(odeprob, react.recommendedsolver, abstol=1e-18, reltol=1e-6)

elseif method == "different_initial_size"

    print("Using different initial size...")

    cell_size_fractions = 10 .^ range(0 - (num_cells - 1), 0, length=num_cells)
    cell_size_fractions ./= sum(cell_size_fractions)
    Cbulk = y0 ./ Vliqinfilm0

    liqinitialconds["epsilon"] = epsilon

    domainfilm, y0film, pfilm = FragmentBasedConstantTrhoDomain(phase=film, initialconds=filminitialconds)
    domainliq, y0liq, pliq = ConstantTLiqFilmDomain(phase=liq_with_rxns, initialconds=liqinitialconds)
    inter, pinter = FragmentBasedReactiveFilmGrowthInterfaceConstantT(domainfilm, domainliq, interfacerxns)
    react, y0, p = Reactor((domainfilm, domainliq), (y0film, y0liq), (0.0, tf), (inter,), (pfilm, pliq, pinter))

    mu = liq_with_rxns.solvent.mu(T)
    diffs = [x(T=T, mu=mu, P=1e8) for x in getfield.(liq_with_rxns.species, :diffusion)]

    num_variables = length(y0)
    y0z = zeros(num_variables, num_cells)

    for j in 1:num_cells
        y0z[:, j] .= y0 * cell_size_fractions[j]
    end

    y0z = reshape(y0z, :)

    function f!(dy, y, p, t)
        f_film_growth!(unflatten(dy), unflatten(y), p, t, react, num_cells, diffs, Cbulk, epsilon, A)
    end

    odefcn = ODEFunction(f!; jac=jacobianyforwarddiff!)
    odeprob = ODEProblem(odefcn, y0z, (0.0, tf), p)

    @time sol = solve(odeprob, react.recommendedsolver, abstol=1e-18, reltol=1e-6)

elseif method == "callback"

    print("Using callback...")

    cell_size_fractions = ones(num_cells)
    cell_size_fractions ./= sum(cell_size_fractions)
    Cbulk = y0 ./ Vliqinfilm0

    liqinitialconds["epsilon"] = epsilon

    domainfilm, y0film, pfilm = FragmentBasedConstantTrhoDomain(phase=film, initialconds=filminitialconds)
    domainliq, y0liq, pliq = ConstantTLiqFilmDomain(phase=liq_with_rxns, initialconds=liqinitialconds)
    inter, pinter = FragmentBasedReactiveFilmGrowthInterfaceConstantT(domainfilm, domainliq, interfacerxns)
    react, y0, p = Reactor((domainfilm, domainliq), (y0film, y0liq), (0.0, tf), (inter,), (pfilm, pliq, pinter))

    mu = liq_with_rxns.solvent.mu(T)
    diffs = [x(T=T, mu=mu, P=1e8) for x in getfield.(liq_with_rxns.species, :diffusion)]

    num_variables = length(y0)
    y0z = zeros(num_variables, num_cells)

    for j in 1:num_cells
        y0z[:, j] .= y0 * cell_size_fractions[j]
    end

    y0z = reshape(y0z, :)

    function f!(dy, y, p, t)
        f_film_growth!(unflatten(dy), unflatten(y), p, t, react, num_cells, diffs, Cbulk, epsilon, A)
    end

    struct my_condition{F<:Function,D<:FragmentBasedReactiveFilmGrowthInterfaceConstantT} <: Function
        condition::F
        num_variables::Int64
        num_cells::Int64
        domainfilm::D
    end
    function (c::my_condition{F,D})(u, t, integrator) where {F,D}
        (; condition, num_variables, num_cells) = c
        u = unflatten(integrator.u, num_variables, num_cells)
        mass_ind = domainfilm.indexes[3]

        min, max = minmax(u[mass_ind, i] for i in axes(u, 2))

        return max / min > 2.0
    end

    temp_u = zeros(num_variables, num_cells)

    function affect!(integrator)
        u = unflatten(integrator.u)
        temp_u .= u
        liq_inds = domainliq.indexes[1]:domainliq.indexes[2]
        liq_volume_ind = domainliq.indexes[3]
        film_inds = domainfilm.indexes[1]:domainfilm.indexes[2]
        mass_ind = domainfilm.indexes[3]

        @views total_mass = sum(u[mass_ind, :])
        new_mass = total_mass / num_cells
        new_liquid_volume = new_mass / domainfilm.rho / (1 - domainliq.epsilon) * domainliq.epsilon
        @views temp_u[mass_ind, :] .= new_mass
        @views temp_u[liq_volume_ind, :] .= new_liquid_volume

        for j in 1:num_cells
            @views upper_surface_volume_change = (sum(temp_u[mass_ind, j:num_cells]) - sum(u[mass_ind, j:num_cells])) / domainfilm.rho / (1 - domainliq.epsilon)
            @views lower_surface_volume_change = (sum(temp_u[mass_ind, j+1:num_cells]) - sum(u[mass_ind, j+1:num_cells])) / domainfilm.rho / (1 - domainliq.epsilon)

            if upper_surface_volume_change >= 0.0
                if j == 1
                    @views temp_u[liq_inds, j] .+= upper_surface_volume_change * domainliq.epsilon * Cbulk[liq_inds]
                else
                    @views temp_u[liq_inds, j] .+= upper_surface_volume_change * domainliq.epsilon * (u[liq_inds, j-1] / u[liq_volume_ind, j-1])
                    @views temp_u[film_inds, j] .+= upper_surface_volume_change * (1 - domainliq.epsilon) * (u[film_inds, j-1] / (u[mass_ind, j-1] / domainfilm.rho))
                end
            else
                @views temp_u[liq_inds, j] .+= -abs(upper_surface_volume_change) * domainliq.epsilon * (u[liq_inds, j] / u[liq_volume_ind, j])
                @views temp_u[film_inds, j] .+= -abs(upper_surface_volume_change) * (1 - domainliq.epsilon) * (u[film_inds, j] / (u[mass_ind, j] / domainfilm.rho))
            end
            if lower_surface_volume_change >= 0.0
                @views temp_u[liq_inds, j] .+= -lower_surface_volume_change * domainliq.epsilon * (u[liq_inds, j] / u[liq_volume_ind, j])
                @views temp_u[film_inds, j] .+= -lower_surface_volume_change * (1 - domainliq.epsilon) * (u[film_inds, j] / (u[mass_ind, j] / domainfilm.rho))
            else
                if j != num_cells
                    @views temp_u[liq_inds, j] .+= abs(lower_surface_volume_change) * domainliq.epsilon * (u[liq_inds, j+1] / u[liq_volume_ind, j+1])
                    @views temp_u[film_inds, j] .+= abs(lower_surface_volume_change) * (1 - domainliq.epsilon) * (u[film_inds, j+1] / (u[mass_ind, j+1] / domainfilm.rho))
                end
            end
        end

        u .= temp_u
    end
    cb = DiscreteCallback(condition, affect!)
    odefcn = ODEFunction(f!; jac=jacobianyforwarddiff!)
    odeprob = ODEProblem(odefcn, y0z, (0.0, tf), p)

    @time sol = solve(odeprob, react.recommendedsolver, abstol=1e-18, reltol=1e-6, callback=cb)

elseif method == "different_initial_size_and_callback"

    print("Using different initial size and callback...")

    cell_size_fractions = 2 .^ range(0 - (num_cells - 1), 0, length=num_cells)
    cell_size_fractions ./= sum(cell_size_fractions)
    Cbulk = y0 ./ Vliqinfilm0

    liqinitialconds["epsilon"] = epsilon

    domainfilm, y0film, pfilm = FragmentBasedConstantTrhoDomain(phase=film, initialconds=filminitialconds)
    domainliq, y0liq, pliq = ConstantTLiqFilmDomain(phase=liq_with_rxns, initialconds=liqinitialconds)
    inter, pinter = FragmentBasedReactiveFilmGrowthInterfaceConstantT(domainfilm, domainliq, interfacerxns)
    react, y0, p = Reactor((domainfilm, domainliq), (y0film, y0liq), (0.0, tf), (inter,), (pfilm, pliq, pinter))

    mu = liq_with_rxns.solvent.mu(T)
    diffs = [x(T=T, mu=mu, P=1e8) for x in getfield.(liq_with_rxns.species, :diffusion)]

    num_variables = length(y0)
    y0z = zeros(num_variables, num_cells)

    for j in 1:num_cells
        y0z[:, j] .= y0 * cell_size_fractions[j]
    end

    y0z = reshape(y0z, :)

    function f!(dy, y, p, t)
        f_film_growth!(unflatten(dy), unflatten(y), p, t, react, num_cells, diffs, Cbulk, epsilon, A)
    end

    function condition(u, t, integrator)
        u = unflatten(integrator.u)
        mass_ind = domainfilm.indexes[3]

        return @views argmin(u[mass_ind, :]) != 1 # want the stencil near the bulk liquid to be the smallest
    end

    temp_u = zeros(num_variables, num_cells)

    function affect!(integrator)
        u = unflatten(integrator.u)
        temp_u .= u
        liq_inds = domainliq.indexes[1]:domainliq.indexes[2]
        liq_volume_ind = domainliq.indexes[3]
        film_inds = domainfilm.indexes[1]:domainfilm.indexes[2]
        mass_ind = domainfilm.indexes[3]

        @views total_mass = sum(u[mass_ind, :])
        new_masses = total_mass .* cell_size_fractions
        new_liquid_volumes = new_masses / domainfilm.rho / (1 - domainliq.epsilon) * domainliq.epsilon
        @views temp_u[mass_ind, :] .= new_masses
        @views temp_u[liq_volume_ind, :] .= new_liquid_volumes

        for j in 1:num_cells
            @views upper_surface_volume_change = (sum(temp_u[mass_ind, j:num_cells]) - sum(u[mass_ind, j:num_cells])) / domainfilm.rho / (1 - domainliq.epsilon)
            @views lower_surface_volume_change = (sum(temp_u[mass_ind, j+1:num_cells]) - sum(u[mass_ind, j+1:num_cells])) / domainfilm.rho / (1 - domainliq.epsilon)

            if upper_surface_volume_change >= 0.0
                if j == 1
                    @views temp_u[liq_inds, j] .+= upper_surface_volume_change * domainliq.epsilon * Cbulk[liq_inds]
                else
                    @views temp_u[liq_inds, j] .+= upper_surface_volume_change * domainliq.epsilon * (u[liq_inds, j-1] / u[liq_volume_ind, j-1])
                    @views temp_u[film_inds, j] .+= upper_surface_volume_change * (1 - domainliq.epsilon) * (u[film_inds, j-1] / (u[mass_ind, j-1] / domainfilm.rho))
                end
            else
                @views temp_u[liq_inds, j] .+= -abs(upper_surface_volume_change) * domainliq.epsilon * (u[liq_inds, j] / u[liq_volume_ind, j])
                @views temp_u[film_inds, j] .+= -abs(upper_surface_volume_change) * (1 - domainliq.epsilon) * (u[film_inds, j] / (u[mass_ind, j] / domainfilm.rho))
            end
            if lower_surface_volume_change >= 0.0
                @views temp_u[liq_inds, j] .+= -lower_surface_volume_change * domainliq.epsilon * (u[liq_inds, j] / u[liq_volume_ind, j])
                @views temp_u[film_inds, j] .+= -lower_surface_volume_change * (1 - domainliq.epsilon) * (u[film_inds, j] / (u[mass_ind, j] / domainfilm.rho))
            else
                if j != num_cells
                    @views temp_u[liq_inds, j] .+= abs(lower_surface_volume_change) * domainliq.epsilon * (u[liq_inds, j+1] / u[liq_volume_ind, j+1])
                    @views temp_u[film_inds, j] .+= abs(lower_surface_volume_change) * (1 - domainliq.epsilon) * (u[film_inds, j+1] / (u[mass_ind, j+1] / domainfilm.rho))
                end
            end
        end

        u .= temp_u
    end
    cb = DiscreteCallback(condition, affect!)
    odefcn = ODEFunction(f!; jac=jacobianyforwarddiff!)
    odeprob = ODEProblem(odefcn, y0z, (0.0, tf), p)

    @time sol = solve(odeprob, react.recommendedsolver, abstol=1e-18, reltol=1e-6, callback=cb)

end

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

if !isdir(save_directory)
    mkdir(save_directory)
end
path = "$(save_directory)/simulation_film_1D_$(perturb_factor_string)_$(tray)_$(num_cells)cells_$(method).csv"
println("Saving results to $(path)")
CSV.write(path, df)

# println("Saving rops...")

# save_rop(sol; suffix="_$(method)")

println("Done!")