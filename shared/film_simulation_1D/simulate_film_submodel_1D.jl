using CSV
using DataFrames

using ReactionMechanismSimulator
using ReactionMechanismSimulator.SciMLBase
using ReactionMechanismSimulator.Sundials
using ReactionMechanismSimulator.YAML
using ReactionMechanismSimulator.SparseArrays
using ReactionMechanismSimulator.ForwardDiff
using ReactionMechanismSimulator.PreallocationTools
using ReactionMechanismSimulator.PyPlot
import ReactionMechanismSimulator: getphasespecies, rops, rops!, getkfskrevs

# read in command line arguments
film_mech_directory = ARGS[1]
model_name = ARGS[2]
liquid_simulation_results_path = ARGS[3]

if model_name in ["basecase_debutanizer_model", "trace_oxygen_perturbed_debutanizer_model"]
    aspen_condition_path = ARGS[4]
    tray = parse(Int64, ARGS[5])
    num_cells = parse(Int64, ARGS[6])
    perturb_factor_string = ARGS[7]
    method = ARGS[8]
    is_discretization_study = ARGS[9] == "true"
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
if is_discretization_study
    save_directory = "discretization_studies"
else
    save_directory = dirname(liquid_simulation_results_path)
end
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
        tf0 = 3600 * 24 * 365
        abstol = 1e-20
        reltol = 1e-7
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
println("Number of liquid species: $(length(liqspcs))")
println("Number of liquid reactions: $(length(liqrxns))")

filmspcs = phaseDict["film"]["Species"]
filmrxns = phaseDict["film"]["Reactions"]
println("Number of film species: $(length(filmspcs))")
println("Number of film reactions: $(length(filmrxns))")

interfacerxns = phaseDict[Set(["liquid", "film"])]["Reactions"]
println("Number of interface reactions: $(length(interfacerxns))")
solvent = phaseDict["Solvents"][1]

liq_without_rxns = IdealDiluteSolution(liqspcs, [], solvent; name="liquid", diffusionlimited=true);
liq_with_rxns = IdealDiluteSolution(liqspcs, liqrxns, solvent; name="liquid", diffusionlimited=true);
liqspcnames = getfield.(liq_with_rxns.species, :name)
liqspcmws = getfield.(liq_with_rxns.species, :molecularweight)
liq_name_mw = Dict(zip(liqspcnames, liqspcmws))
liq_name_smiles = Dict(zip(liqspcnames, getfield.(liq_with_rxns.species, :smiles)))
liqrxnstrs = getrxnstr.(liq_with_rxns.reactions)
liqrxncomments = getfield.(liq_with_rxns.reactions, :comment)

film = FragmentBasedIdealFilm(filmspcs, filmrxns; name="film", diffusionlimited=false);
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
    filminitialconds["CD"] = sum(Float64[liquid_steady_state_mols[tray, name] for name in liquid_species_mapping["conjugated_diene(L)"]]) / liqtotalmass * mass0

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

domainliq, y0liq, pliq = ConstantTVDomain(phase=liq_without_rxns, initialconds=liqinitialconds, constantspecies=liqspcnames);

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

if !is_discretization_study
    println("Saving asymptotic steady state film phase concentrations...")
    path = "$(save_directory)/simulation_film_$(tray)_asymptotic.csv"
    CSV.write(path, df)
end

println("Initializing fragment concentrations with steady state fragment concentrations...")

for fragment in fragmentnames
    if fragment != "inert(S)"
        new_value = df[end, fragment] / df[end, "mass"] * mass0
        println(fragment, " (mol/kg): ", filminitialconds[fragment] / mass0, " ", new_value / mass0)
        if new_value > 0.0
            filminitialconds[fragment] = new_value
        else
            filminitialconds[fragment] = 0.0 # can be small negative value for oxygen species in certain trays 
        end
    end
end

for (key, value) in liqinitialconds
    println(key, ": ", value)
end

for (key, value) in filminitialconds
    println(key, ": ", value)
end

println("Solving film phase submodel with steady state fragment concentrations...")

import Base: size, getindex, setindex!

struct MyMatrix{T} <: AbstractMatrix{T}
    vector::Array{T}
    num_variables::Int64
    num_cells::Int64
end

size(A::MyMatrix{T}) where {T} = (A.num_variables, A.num_cells)

function getindex(A::MyMatrix{T}, i::Int64, j::Int64) where {T}
    return A.vector[i+(j-1)*A.num_variables]
end

function setindex!(A::MyMatrix{T}, value::T, i::Int64, j::Int64) where {T}
    A.vector[i+(j-1)*A.num_variables] = value
end

struct Unflatten
    num_variables::Int64
    num_cells::Int64
end

(unflatten::Unflatten)(y::T) where {T<:AbstractArray} = MyMatrix(y, unflatten.num_variables, unflatten.num_cells)

struct MyGetIndex
    num_variables::Int64
    num_cells::Int64
end

function (mygetindex::MyGetIndex)(y::T, i::Int64, j::Int64) where {T<:Vector}
    return y[i+(j-1)*mygetindex.num_variables]
end

struct OneDimensionalFilmGrowth{F1<:Function}
    reaction!::F1
    liq_inds::UnitRange{Int64}
    liq_volume_ind::Int64
    film_inds::UnitRange{Int64}
    mass_ind::Int64
    rho::Float64
    A::Float64
    epsilon::Float64
    num_cells::Int64
    num_variables::Int64
    diffs::Vector{Float64}
    csbulk::Vector{Float64}
end

function (one_d_film_growth::OneDimensionalFilmGrowth)(dy::T1, y::T2, p::T3, t::T4) where {T1<:AbstractArray,T2<:AbstractArray,T3<:AbstractArray,T4<:Number}

    dy = unflatten(dy)
    y = unflatten(y)
    dy .= 0.0
    (; reaction!, liq_inds, liq_volume_ind, film_inds, mass_ind, rho, A, epsilon, num_cells, num_variables, diffs, csbulk) = one_d_film_growth

    for j in 1:num_cells

        # reaction terms
        @views reaction!(dy[:, j], y[:, j], p, t)

        # diffusion terms and boundary conditions

        for ind in liq_inds
            if j == num_cells
                # no flux condition at tray surface
                Jjp1half = 0.0
            else
                delta_z_j = y[mass_ind, j] / rho / A / (1 - epsilon)
                delta_z_jp1 = y[mass_ind, j+1] / rho / A / (1 - epsilon)
                @views Jjp1half = -diffs[ind-liq_inds[1]+1] * (y[ind, j+1] / y[liq_volume_ind, j+1] - y[ind, j] / y[liq_volume_ind, j]) / (delta_z_jp1 / 2 + delta_z_j / 2)
            end
            if j == 1
                # C = C bulk in liquid
                delta_z_j = y[mass_ind, j] / rho / A / (1 - epsilon)
                @views Jjm1half = -diffs[ind-liq_inds[1]+1] * (y[ind, j] / y[liq_volume_ind, j] - csbulk[ind]) / (delta_z_j / 2 + delta_z_j / 2)
            else
                delta_z_j = y[mass_ind, j] / rho / A / (1 - epsilon)
                delta_z_jm1 = y[mass_ind, j-1] / rho / A / (1 - epsilon)
                @views Jjm1half = -diffs[ind-liq_inds[1]+1] * (y[ind, j] / y[liq_volume_ind, j] - y[ind, j-1] / y[liq_volume_ind, j-1]) / (delta_z_j / 2 + delta_z_jm1 / 2)
            end
            @views dy[ind, j] += -(Jjp1half - Jjm1half) * A
        end

        @views upper_surface_volume_change = sum(dy[mass_ind, j:num_cells]) / rho / (1 - epsilon)
        @views lower_surface_volume_change = sum(dy[mass_ind, j+1:num_cells]) / rho / (1 - epsilon)

        for ind in liq_inds
            if upper_surface_volume_change >= 0.0
                if j == 1
                    dy[ind, j] += upper_surface_volume_change * epsilon * csbulk[ind]
                else
                    dy[ind, j] += upper_surface_volume_change * epsilon * (y[ind, j-1] / y[liq_volume_ind, j-1])
                end
            else
                dy[ind, j] -= abs(upper_surface_volume_change) * epsilon * (y[ind, j] / y[liq_volume_ind, j])
            end
            if lower_surface_volume_change >= 0.0
                dy[ind, j] -= lower_surface_volume_change * epsilon * (y[ind, j] / y[liq_volume_ind, j])
            else
                if j != num_cells
                    dy[ind, j] += abs(lower_surface_volume_change) * epsilon * (y[ind, j+1] / y[liq_volume_ind, j+1])
                end
            end
        end
    end
end

function get_jacobian_prototype(one_d_film_growth::OneDimensionalFilmGrowth)

    (; reaction!, liq_inds, liq_volume_ind, film_inds, mass_ind, rho, A, epsilon, num_cells, num_variables, diffs, csbulk) = one_d_film_growth
    J = zeros(num_variables, num_cells, num_variables, num_cells)

    for j in 1:num_cells

        # reaction terms
        # @views reaction!(dy[:, j], y[:, j], p, t)
        J[:, j, :, j] .= 1.0

        # diffusion terms and boundary conditions

        for ind in liq_inds
            if j == num_cells
                # no flux condition at tray surface
                Jjp1half = 0.0
            else
                # delta_z_j = y[mass_ind, j] / rho / A / (1 - epsilon)
                J[ind, j, mass_ind, j] = 1.0
                # delta_z_jp1 = y[mass_ind, j+1] / rho / A / (1 - epsilon)
                J[ind, j, mass_ind, j+1] = 1.0
                # @views Jjp1half = -diffs[ind] * (y[ind, j+1] / y[liq_volume_ind, j+1] - y[ind, j] / y[liq_volume_ind, j]) ./ (delta_z_jp1 / 2 + delta_z_j / 2)
                J[ind, j, ind, j+1] = 1.0
                J[ind, j, liq_volume_ind, j] = 1.0
            end
            if j == 1
                # C = C bulk in liquid
                # delta_z_j = y[mass_ind, j] / rho / A / (1 - epsilon)
                J[ind, j, mass_ind, j] = 1.0
                # @views Jjm1half = -diffs[ind] * (y[ind, j] / y[liq_volume_ind, j] .- csbulk[ind]) ./ (delta_z_j / 2 + delta_z_j / 2)
                J[ind, j, ind, j] = 1.0
                J[ind, j, liq_volume_ind, j] = 1.0
            else
                # delta_z_j = y[mass_ind, j] / rho / A / (1 - epsilon)
                J[ind, j, mass_ind, j] = 1.0
                # delta_z_jm1 = y[mass_ind, j-1] / rho / A / (1 - epsilon)
                J[ind, j, mass_ind, j] = 1.0
                # @views Jjm1half = -diffs * (y[ind, j] / y[liq_volume_ind, j] .- y[ind, j-1] / y[liq_volume_ind, j-1]) ./ (delta_z_j / 2 + delta_z_jm1 / 2)
                J[ind, j, ind, j] = 1.0
                J[ind, j, liq_volume_ind, j] = 1.0
                J[ind, j, ind, j-1] = 1.0
                J[ind, j, liq_volume_ind, j-1] = 1.0
            end
            # @views dy[ind, j] += -(Jjp1half - Jjm1half) * A
        end

        # @views upper_surface_volume_change = sum(dy[mass_ind, j:num_cells]) / rho / A / (1 - epsilon)
        # @views lower_surface_volume_change = sum(dy[mass_ind, j+1:num_cells]) / rho / A / (1 - epsilon)

        for ind in liq_inds
            # if upper_surface_volume_change >= 0.0
            if j == 1
                # dy[ind, j] += upper_surface_volume_change * epsilon * A * csbulk[ind]
                J[ind, j, mass_ind, j:num_cells] .= 1.0
            else
                # dy[ind, j] += upper_surface_volume_change * epsilon * A * (y[ind, j-1] / y[liq_volume_ind, j-1])
                J[ind, j, mass_ind, j:num_cells] .= 1.0
                J[ind, j, ind, j-1] = 1.0
                J[ind, j, liq_volume_ind, j-1] = 1.0
            end

            # else
            # dy[ind, j] -= abs(upper_surface_volume_change) * epsilon * A * (y[ind, j] / y[liq_volume_ind, j])
            J[ind, j, mass_ind, j:num_cells] .= 1.0
            J[ind, j, ind, j] = 1.0
            J[ind, j, liq_volume_ind, j] = 1.0
            # end

            # if lower_surface_volume_change >= 0.0
            # dy[ind, j] -= lower_surface_volume_change * epsilon * A * (y[ind, j] / y[liq_volume_ind, j])
            J[ind, j, mass_ind, j+1:num_cells] .= 1.0
            J[ind, j, ind, j] = 1.0
            J[ind, j, liq_volume_ind, j] = 1.0
            # else
            if j != num_cells
                # dy[ind, j] += abs(lower_surface_volume_change) * epsilon * A * (y[ind, j+1] / y[liq_volume_ind, j+1])
                J[ind, j, mass_ind, j+1:num_cells] .= 1.0
                J[ind, j, ind, j+1] = 1.0
                J[ind, j, liq_volume_ind, j+1] = 1.0
            end
            # end
        end
    end
    return J
end

struct MyJacobianYForwardDiff
    dy0::Vector{Float64}
end

function (myjacobianyforwarddiff::MyJacobianYForwardDiff)(J, y, p, t)
    ForwardDiff.jacobian!(J, (dy, y) -> f!(dy, y, p, t), myjacobianyforwarddiff.dy0, y)
end

struct MyAffect
    liq_inds::UnitRange{Int64}
    liq_volume_ind::Int64
    film_inds::UnitRange{Int64}
    mass_ind::Int64
    cell_size_fractions::Vector{Float64}
    epsilon::Float64
    rho::Float64
    csbulk::Vector{Float64}
    num_cells::Int64
    temp_u::Matrix{Float64}
    lower_cell_surface_locations::Vector{Float64}
end

function (myaffect::MyAffect)(integrator)
    t = integrator.t
    # println("t = $t")

    (; liq_inds, liq_volume_ind, film_inds, mass_ind, cell_size_fractions, epsilon, rho, csbulk, num_cells, temp_u, lower_cell_surface_locations) = myaffect

    u = unflatten(integrator.u)
    temp_u .= 0.0

    @views total_mass = sum(u[mass_ind, :])
    new_masses = total_mass .* cell_size_fractions
    new_liquid_volumes = new_masses / rho / (1 - epsilon) * epsilon
    for j in 1:num_cells
        temp_u[mass_ind, j] = new_masses[j]
        temp_u[liq_volume_ind, j] = new_liquid_volumes[j]
    end

    # println("u[mass_ind, :]: ", u[mass_ind, :])
    # println("u[liq_volume_ind, :]: ", u[liq_volume_ind, :])
    # println("sum(u[mass_ind, :]): ", sum(u[mass_ind, :]))
    # println("sum(u[liq_volume_ind, :]): ", sum(u[liq_volume_ind, :]))
    # println("u[liq_inds[1], :]: ", u[liq_inds[1], :])
    # println("sum(u[liq_inds[1], :]): ", sum(u[liq_inds[1], :]))

    @views lower_cell_surface_locations .= cumsum(u[mass_ind, :]) / rho / (1 - epsilon) / A

    # println("lower_cell_surface_locations: ", lower_cell_surface_locations)

    for j in 1:num_cells

        @views new_upper_surface_location = sum(temp_u[mass_ind, 1:j-1]) / rho / (1 - epsilon) / A
        @views new_lower_surface_location = sum(temp_u[mass_ind, 1:j]) / rho / (1 - epsilon) / A

        # println("new_upper_surface_location: ", new_upper_surface_location)
        # println("new_lower_surface_location: ", new_lower_surface_location)

        new_upper_surface_cell_ind = findfirst(x -> x >= new_upper_surface_location, lower_cell_surface_locations)
        new_lower_surface_cell_ind = findfirst(x -> x >= new_lower_surface_location, lower_cell_surface_locations)
        if new_lower_surface_cell_ind === nothing
            new_lower_surface_cell_ind = num_cells
        end

        if new_upper_surface_cell_ind == new_lower_surface_cell_ind
            for ind in liq_inds
                temp_u[ind, j] += (new_lower_surface_location - new_upper_surface_location) * A * epsilon * (u[ind, new_upper_surface_cell_ind] / u[liq_volume_ind, new_upper_surface_cell_ind])
                # if ind == liq_inds[1]
                #     println("new_lower_surface_cell_ind: ", new_lower_surface_cell_ind)
                #     println("new_upper_surface_cell_ind: ", new_upper_surface_cell_ind)
                #     println("temp_u[ind, j]: ", temp_u[ind, j])
                # end
            end
            for ind in film_inds
                temp_u[ind, j] += (new_lower_surface_location - new_upper_surface_location) * A * (1 - epsilon) * (u[ind, new_upper_surface_cell_ind] / (u[mass_ind, new_upper_surface_cell_ind] / rho))
            end
        else
            for ind in liq_inds
                temp_u[ind, j] += (lower_cell_surface_locations[new_upper_surface_cell_ind] - new_upper_surface_location) * A * epsilon * (u[ind, new_upper_surface_cell_ind] / u[liq_volume_ind, new_upper_surface_cell_ind])
                # if ind == liq_inds[1]
                #     println("new_upper_surface_cell_ind: ", new_upper_surface_cell_ind)
                #     println("temp_u[ind, j]: ", temp_u[ind, j])
                # end
            end
            for ind in film_inds
                temp_u[ind, j] += (lower_cell_surface_locations[new_upper_surface_cell_ind] - new_upper_surface_location) * A * (1 - epsilon) * (u[ind, new_upper_surface_cell_ind] / (u[mass_ind, new_upper_surface_cell_ind] / rho))
            end
            for middle_cell_ind in new_upper_surface_cell_ind+1:new_lower_surface_cell_ind-1
                for ind in liq_inds
                    temp_u[ind, j] += u[ind, middle_cell_ind]
                    # if ind == liq_inds[1]
                    #     println("middle_cell_ind: ", middle_cell_ind)
                    #     println("temp_u[ind, j]: ", temp_u[ind, j])
                    # end
                end
                for ind in film_inds
                    temp_u[ind, j] += u[ind, middle_cell_ind]
                end
            end
            for ind in liq_inds
                temp_u[ind, j] += (new_lower_surface_location - lower_cell_surface_locations[new_lower_surface_cell_ind-1]) * A * epsilon * (u[ind, new_lower_surface_cell_ind] / u[liq_volume_ind, new_lower_surface_cell_ind])
                # if ind == liq_inds[1]
                #     println("new_lower_surface_cell_ind: ", new_lower_surface_cell_ind)
                #     println("temp_u[ind, j]: ", temp_u[ind, j])
                # end
            end
            for ind in film_inds
                temp_u[ind, j] += (new_lower_surface_location - lower_cell_surface_locations[new_lower_surface_cell_ind-1]) * A * (1 - epsilon) * (u[ind, new_lower_surface_cell_ind] / (u[mass_ind, new_lower_surface_cell_ind] / rho))
            end
        end
    end

    u .= temp_u

    # println("u[mass_ind, :]: ", u[mass_ind, :])
    # println("u[liq_volume_ind, :]: ", u[liq_volume_ind, :])
    # println("sum(u[mass_ind, :]): ", sum(u[mass_ind, :]))
    # println("sum(u[liq_volume_ind, :]): ", sum(u[liq_volume_ind, :]))
    # println("u[liq_inds[1], :]: ", u[liq_inds[1], :])
    # println("sum(u[liq_inds[1], :]): ", sum(u[liq_inds[1], :]))

end

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

            ropmat = rops_1D(ssys, sol.t[end], cell_ind)

            rxnind, spcind, rop = findnz(ropmat)
            no_nan = any((isnan).(rop)) == false
            println("No NaN ", no_nan)

            if no_nan

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
                CSV.write("$(save_directory)/simulation_film_rop_$(tray).csv", df_rop, quotestrings=true)

                rops, rop_rxnstrs, rop_rxncomments = get_rops(df_rop, "mass")
                no_neg_rop = all((rops .>= 0) .|| occursin.("<=>", rop_rxncomments)) # positive rop except for reversible reactions
                println("No negative rop ", no_neg_rop)

                ind_max = argmax(rops)
                println("Max rop ", rops[ind_max])
                println("Max rop rxn ", rop_rxnstrs[ind_max])
                println("Max rop rxn comment ", rop_rxncomments[ind_max])
                if model_name == "basecase_debutanizer_model"
                    check_max_rxn = occursin("CYCLOPENTADIENE(L) + CDB Diels-Alder addition", rop_rxncomments[ind_max]) || occursin("1,3-BUTADIENE(L) + AR radical addition", rop_rxncomments[ind_max]) || occursin("AR + CYCLOPENTADIENE(L) radical addition", rop_rxncomments[ind_max])
                elseif model_name == "trace_oxygen_perturbed_debutanizer_model"
                    check_max_rxn = occursin("CYCLOPENTADIENE(L) + CDB Diels-Alder addition", rop_rxncomments[ind_max]) || occursin("1,3-BUTADIENE(L) + AR radical addition", rop_rxncomments[ind_max]) || occursin("AR + CYCLOPENTADIENE(L) radical addition", rop_rxncomments[ind_max])
                end

                if no_neg_rop && check_max_rxn

                    CSV.write("$(save_directory)/simulation_film_1D_rop_$(tray)_$(num_cells)cells_cell$(cell_ind)$(suffix).csv", df_rop, quotestrings=true)

                    break

                end
            end

            count += 1

            if count > 100
                break
            end
        end

        count = 0

        while true

            ropmat = rops_1D(ssys, 0.0, cell_ind)

            rxnind, spcind, rop = findnz(ropmat)
            no_nan = any((isnan).(rop)) == false
            println("No NaN ", no_nan)

            if no_nan

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

                rops, rop_rxnstrs, rop_rxncomments = get_rops(df_rop, "mass")
                no_neg_rop = all((rops .>= 0) .|| occursin.("<=>", rop_rxncomments)) # positive rop except for reversible reactions
                println("No negative rop ", no_neg_rop)

                ind_max = argmax(rops)
                println("Max rop ", rops[ind_max])
                println("Max rop rxn ", rop_rxnstrs[ind_max])
                println("Max rop rxn comment ", rop_rxncomments[ind_max])
                if model_name == "basecase_debutanizer_model"
                    check_max_rxn = occursin("CYCLOPENTADIENE(L) + CDB Diels-Alder addition", rop_rxncomments[ind_max]) || occursin("1,3-BUTADIENE(L) + AR radical addition", rop_rxncomments[ind_max]) || occursin("AR + CYCLOPENTADIENE(L) radical addition", rop_rxncomments[ind_max])
                elseif model_name == "trace_oxygen_perturbed_debutanizer_model"
                    check_max_rxn = occursin("CYCLOPENTADIENE(L) + CDB Diels-Alder addition", rop_rxncomments[ind_max]) || occursin("1,3-BUTADIENE(L) + AR radical addition", rop_rxncomments[ind_max]) || occursin("AR + CYCLOPENTADIENE(L) radical addition", rop_rxncomments[ind_max]) || occursin("1,3-BUTADIENE(L) + PR radical addition", rop_rxncomments[ind_max])
                elseif model_name == "QCMD_cell_model"
                    check_max_rxn = occursin("AR + 2-methylcyclohexadiene(L) radical addition", rop_rxncomments[ind_max]) || occursin("[O][O](L) + AR <=> PR", rop_rxncomments[ind_max])
                end

                if no_neg_rop && check_max_rxn

                    CSV.write("$(save_directory)/simulation_film_1D_rop_$(tray)_$(num_cells)cells_cell$(cell_ind)$(suffix)_t0.csv", df_rop, quotestrings=true)

                    break

                end
            end

            count += 1

            if count > 100
                break
            end
        end
    end
end

if method == "same_initial_size"

    println("Using same initial size...")

    csbulk = y0 ./ Vliqinfilm0

    cell_size_fractions = ones(num_cells)
    cell_size_fractions ./= sum(cell_size_fractions)

    liqinitialconds["epsilon"] = epsilon

    domainfilm, y0film, pfilm = FragmentBasedConstantTrhoDomain(phase=film, initialconds=filminitialconds)
    domainliq, y0liq, pliq = ConstantTLiqFilmDomain(phase=liq_with_rxns, initialconds=liqinitialconds)
    inter, pinter = FragmentBasedReactiveFilmGrowthInterfaceConstantT(domainfilm, domainliq, interfacerxns)
    react, y0, p = Reactor((domainfilm, domainliq), (y0film, y0liq), (0.0, tf), (inter,), (pfilm, pliq, pinter))

    num_variables = length(y0)
    unflatten = Unflatten(num_variables, num_cells)
    liq_inds = domainliq.indexes[1]:domainliq.indexes[2]
    liq_volume_ind = domainliq.indexes[3]
    film_inds = domainfilm.indexes[1]:domainfilm.indexes[2]
    mass_ind = domainfilm.indexes[3]

    mu = liq_with_rxns.solvent.mu(T)
    diffs = [x(T=T, mu=mu, P=1e8) for x in getfield.(liq_with_rxns.species, :diffusion)]

    y0z = zeros(num_variables, num_cells)

    for j in 1:num_cells
        @views y0z[:, j] .= y0 * cell_size_fractions[j]
    end

    y0z = reshape(y0z, :)

    f! = OneDimensionalFilmGrowth(react.ode.f, liq_inds, liq_volume_ind, film_inds, mass_ind, rho, A, epsilon, num_cells, num_variables, diffs, csbulk)
    jacobianyforwarddiff! = MyJacobianYForwardDiff(zeros(num_variables * num_cells))

    odefcn = ODEFunction(f!; jac=jacobianyforwarddiff!)
    odeprob = ODEProblem(odefcn, y0z, (0.0, tf), p)

    println("Using default :Dense Solver...")
    solver = Sundials.CVODE_BDF()
    @time sol = solve(odeprob, solver, abstol=abstol, reltol=reltol)

    println(sol.stats)

    # println("Using :GMRES Solver without jac_prototype...")
    # solver = Sundials.CVODE_BDF(linear_solver=:GMRES)
    # @time sol = solve(odeprob, solver, abstol=abstol, reltol=reltol)

    # jac_prototype = get_jacobian_prototype(f!)
    # jac_prototype = reshape(jac_prototype, (num_variables * num_cells, num_variables * num_cells))
    # jac_prototype = sparse(jac_prototype)
    # println("Nonzero entries in jacobian prototype: ", length(findnz(jac_prototype))/length(jac_prototype))

    # odefcn = ODEFunction(f!; jac=jacobianyforwarddiff!, jac_prototype=float.(jac_prototype))
    # odeprob = ODEProblem(odefcn, y0z, (0.0, tf), p)

    # println("Using :GMRES Solver with jac_prototype...")
    # solver = Sundials.CVODE_BDF(linear_solver=:GMRES)
    # @time sol = solve(odeprob, solver, abstol=abstol, reltol=reltol)

    # println("Using :KLU Solver with jac_prototype...")
    # solver = Sundials.CVODE_BDF(linear_solver=:KLU)
    # @time sol = solve(odeprob, solver, abstol=abstol, reltol=reltol)

elseif method == "callback"

    println("Using callback...")

    csbulk = y0 ./ Vliqinfilm0

    cell_size_fractions = ones(num_cells)
    cell_size_fractions ./= sum(cell_size_fractions)

    liqinitialconds["epsilon"] = epsilon

    domainfilm, y0film, pfilm = FragmentBasedConstantTrhoDomain(phase=film, initialconds=filminitialconds)
    domainliq, y0liq, pliq = ConstantTLiqFilmDomain(phase=liq_with_rxns, initialconds=liqinitialconds)
    inter, pinter = FragmentBasedReactiveFilmGrowthInterfaceConstantT(domainfilm, domainliq, interfacerxns)
    react, y0, p = Reactor((domainfilm, domainliq), (y0film, y0liq), (0.0, tf), (inter,), (pfilm, pliq, pinter))

    num_variables = length(y0)
    unflatten = Unflatten(num_variables, num_cells)
    liq_inds = domainliq.indexes[1]:domainliq.indexes[2]
    liq_volume_ind = domainliq.indexes[3]
    film_inds = domainfilm.indexes[1]:domainfilm.indexes[2]
    mass_ind = domainfilm.indexes[3]

    mu = liq_with_rxns.solvent.mu(T)
    diffs = [x(T=T, mu=mu, P=1e8) for x in getfield.(liq_with_rxns.species, :diffusion)]

    y0z = zeros(num_variables, num_cells)

    for j in 1:num_cells
        y0z[:, j] .= y0 * cell_size_fractions[j]
    end

    y0z = reshape(y0z, :)

    struct MyCondition
        mass_ind::Int64
        num_variables::Int64
        cell_inds::UnitRange{Int64}
    end

    function (mycondition::MyCondition)(u, t, integrator)
        (; mass_ind, num_variables, cell_inds) = mycondition
        min_mass = 1e100
        for cell_ind in cell_inds
            @views min_mass = min(min_mass, minimum(u[mass_ind+(cell_ind-1)*num_variables]))
        end
        max_mass = 0.0
        for cell_ind in cell_inds
            @views max_mass = max(max_mass, maximum(u[mass_ind+(cell_ind-1)*num_variables]))
        end
        return max_mass / min_mass > 2.0
    end

    condition = MyCondition(mass_ind, num_variables, 1:num_cells)
    temp_u = zeros(num_variables, num_cells)
    lower_cell_surface_locations = zeros(num_cells)
    affect! = MyAffect(liq_inds, liq_volume_ind, film_inds, mass_ind, cell_size_fractions, epsilon, rho, csbulk, num_cells, temp_u, lower_cell_surface_locations)

    cb = DiscreteCallback(condition, affect!)
    f! = OneDimensionalFilmGrowth(react.ode.f, liq_inds, liq_volume_ind, film_inds, mass_ind, rho, A, epsilon, num_cells, num_variables, diffs, csbulk)
    jacobianyforwarddiff! = MyJacobianYForwardDiff(zeros(num_variables * num_cells))

    odefcn = ODEFunction(f!; jac=jacobianyforwarddiff!)
    odeprob = ODEProblem(odefcn, y0z, (0.0, tf), p)

    println("Using default :Dense Solver...")
    solver = Sundials.CVODE_BDF()
    @time sol = solve(odeprob, solver, abstol=abstol, reltol=reltol, callback=cb)

    println(sol.stats)

    # println("Using :GMRES Solver without jac_prototype...")
    # solver = Sundials.CVODE_BDF(linear_solver=:GMRES)
    # @time sol = solve(odeprob, solver, abstol=abstol, reltol=reltol, callback=cb)

    # jac_prototype = get_jacobian_prototype(f!)
    # jac_prototype = reshape(jac_prototype, (num_variables * num_cells, num_variables * num_cells))
    # jac_prototype = sparse(jac_prototype)
    # println("Nonzero entries in jacobian prototype: ", length(findnz(jac_prototype))/length(jac_prototype))

    # odefcn = ODEFunction(f!; jac=jacobianyforwarddiff!, jac_prototype=float.(jac_prototype))
    # odeprob = ODEProblem(odefcn, y0z, (0.0, tf), p)

    # println("Using :GMRES Solver with jac_prototype...")
    # solver = Sundials.CVODE_BDF(linear_solver=:GMRES)
    # @time sol = solve(odeprob, solver, abstol=abstol, reltol=reltol, callback=cb)

    # println("Using :KLU Solver with jac_prototype...")
    # solver = Sundials.CVODE_BDF(linear_solver=:KLU)
    # @time sol = solve(odeprob, solver, abstol=abstol, reltol=reltol, callback=cb)

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

if is_discretization_study
    path = "$(save_directory)/simulation_film_1D_$(perturb_factor_string)_$(tray)_$(num_cells)cells_$(method).csv"
else
    path = "$(save_directory)/simulation_film_1D_$(tray)_$(num_cells)cells_$(method).csv"
end

println("Saving results to $(path)")
CSV.write(path, df)

println("Saving rops...")

save_rop(sol; suffix="_$(method)")

println("Done!")