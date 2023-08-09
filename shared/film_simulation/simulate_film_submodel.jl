# save_directory = "1,3-BUTADIENE_1.0"
# import Pkg
# Pkg.activate("/home/gridsan/hwpang/Jobs/Dow")
# Pkg.instantiate()

using CSV
using DataFrames

using ReactionMechanismSimulator
using ReactionMechanismSimulator.SciMLBase
using ReactionMechanismSimulator.Sundials
using ReactionMechanismSimulator.YAML
using ReactionMechanismSimulator.SparseArrays

# read in command line arguments
rms_mech_directory = ARGS[1]
model_name = ARGS[2]
liquid_simulation_results_path = ARGS[3]

if model_name in ["basecase_debutanizer_model", "trace_oxygen_perturbed_debutanizer_model"]
    aspen_condition_path = ARGS[4]
    tray = parse(Int, ARGS[5])
elseif model_name == "QCMD_cell_model"
    tray = 1
    nothing
end

println("rms_mech_directory: $(rms_mech_directory)")
println("model_name: $(model_name)")
println("liquid_simulation_results_path: $(liquid_simulation_results_path)")
if model_name == "basecase_debutanizer_model"
    println("aspen_condition_path: $(aspen_condition_path)")
end

# set up paths
save_directory = dirname(liquid_simulation_results_path)
println("save_directory: $(save_directory)")
liquid_species_mapping_path = joinpath(rms_mech_directory, "liquid_species_mapping.yml")
println("liquid_species_mapping_path: $(liquid_species_mapping_path)")
ASFWnparams_path = joinpath(save_directory, "ASFWnparams.yml")
println("ASFWnparams_path: $(ASFWnparams_path)")

# set up parameters
if model_name in ["basecase_debutanizer_model", "trace_oxygen_perturbed_debutanizer_model"]
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

    if model_name == "basecase_debutanizer_model"
        include_oxygen = false
        tf0 = 3600 * 24 * 365 * 50
        abstol = 1e-18
        reltol = 1e-6
    else
        include_oxygen = true
        tf0 = 3600 * 24 * 365
        abstol = 1e-18
        reltol = 1e-6
    end

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

    tf0 = 3600.0

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
    oligomer_init_mol["COO.(L,oligomer)"] = oligomer_COO_mol
    oligomer_init_mol["COOH(L,oligomer)"] = oligomer_COOH_mol
end

for (oligomer, mol) in oligomer_init_mol
    println(oligomer, " (mol/m^3) = ", mol / Vliq)
end

liqtotalmass = sum(liquid_steady_state_mols[tray, name] * liq_name_mw[name*"(L)"] for name in names(liquid_steady_state_mols))

liqinitialconds = Dict{String,Float64}()
liqinitialconds["T"] = T
liqinitialconds["V"] = Vliqinfilm0
for name in names(liquid_steady_state_mols)
    if name * "(L)" in liqspcnames
        liqinitialconds[name*"(L)"] = liquid_steady_state_mols[tray, name] / Vliq * Vliqinfilm0
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
        filminitialconds["PR"] = sum(Float64[liquid_steady_state_mols[tray, name] for name in liquid_species_mapping["COO.(L)"]]) / liqtotalmass * mass0
        filminitialconds["CP"] = sum(Float64[liquid_steady_state_mols[tray, name] for name in liquid_species_mapping["COOC(L)"]]) / liqtotalmass * mass0
        filminitialconds["HP"] = sum(Float64[liquid_steady_state_mols[tray, name] for name in liquid_species_mapping["COOH(L)"]]) / liqtotalmass * mass0
        filminitialconds["OR"] = sum(Float64[liquid_steady_state_mols[tray, name] for name in liquid_species_mapping["CO.(L)"]]) / liqtotalmass * mass0
        filminitialconds["OH"] = 0.0
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

@time sol = solve(react.ode, react.recommendedsolver, abstol=abstol, reltol=reltol);

df = DataFrame(sol)
rename!(df, names(df)[domainliq.indexes[1]+1:domainliq.indexes[2]+1] .=> liqspcnames)
rename!(df, names(df)[domainfilm.indexes[1]+1:domainfilm.indexes[2]+1] .=> filmspcnames)
rename!(df, names(df)[domainfilm.indexes[3]+1] .=> "mass")

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

function save_rop(sol)

    count = 0

    while true
        ssys = SystemSimulation(sol, (domainfilm, domainliq,), (inter,), p)
        spcnames = ssys.names

        ropmat = ReactionMechanismSimulator.rops(ssys, sol.t[end])

        rxnind, spcind, rop = findnz(ropmat)
        no_nan = any((isnan).(rop)) == false
        println("No NaN ", no_nan)

        if no_nan

            rxnind = trunc.(Int, rxnind)
            spcind = trunc.(Int, spcind)
            rop_rxnstrs = allrxnstrs[rxnind]
            rop_rxncomments = allrxncomments[rxnind]
            spcmassnames = spcnames[:]
            push!(spcmassnames, "mass")
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
                check_max_rxn = occursin("CYCLOPENTADIENE(L) + CDB Diels-Alder addition", rop_rxncomments[ind_max]) || occursin("1,3-BUTADIENE(L) + AR radical addition", rop_rxncomments[ind_max]) || occursin("AR + CYCLOPENTADIENE(L) radical addition", rop_rxncomments[ind_max]) || occursin("1,3-BUTADIENE(L) + PR radical addition", rop_rxncomments[ind_max])
            elseif model_name == "QCMD_cell_model"
                check_max_rxn = occursin("AR + 2-methylcyclohexadiene(L) radical addition", rop_rxncomments[ind_max]) || occursin("[O][O](L) + AR <=> PR", rop_rxncomments[ind_max])
            end

            if no_neg_rop && check_max_rxn

                CSV.write("$(save_directory)/simulation_film_rop_$(tray).csv", df_rop, quotestrings=true)

                break

            end
        end

        count += 1

        if count > 100
            break
        end
    end
end

if model_name in ["basecase_debutanizer_model", "trace_oxygen_perturbed_debutanizer_model"]
    CSV.write("$(save_directory)/simulation_film_$(tray)_asymptotic.csv", df)
end

if model_name == "QCMD_cell_model"
    CSV.write("$(save_directory)/simulation_film_$(tray).csv", df)
    save_rop(sol)
else

    println("Initializing fragment concentrations with steady state fragment concentrations...")

    for fragment in fragmentnames
        if fragment != "inert(S)"
            println(fragment, " (mol/mass): ", filminitialconds[fragment] / mass0, " ", df[end, fragment] / df[end, "mass"])
            filminitialconds[fragment] = df[end, fragment] / df[end, "mass"] * mass0
        end
    end

    println("Solving film phase submodel with steady state fragment concentrations...")

    domainfilm, y0film, pfilm = FragmentBasedConstantTrhoDomain(phase=film, initialconds=filminitialconds)

    domainliq, y0liq, pliq = ConstantTVDomain(phase=liq, initialconds=liqinitialconds, constantspecies=liqspcnames)

    inter, pinter = FragmentBasedReactiveFilmGrowthInterfaceConstantT(domainfilm, domainliq, interfacerxns)

    react, y0, p = Reactor((domainfilm, domainliq), (y0film, y0liq), (0.0, tf), (inter,), (pfilm, pliq, pinter))

    @time sol = solve(react.ode, react.recommendedsolver, abstol=abstol, reltol=reltol)

    println("Saving film phase submodel results...")

    df = DataFrame(sol)
    rename!(df, names(df)[domainliq.indexes[1]+1:domainliq.indexes[2]+1] .=> liqspcnames)
    rename!(df, names(df)[domainfilm.indexes[1]+1:domainfilm.indexes[2]+1] .=> filmspcnames)
    rename!(df, names(df)[domainfilm.indexes[3]+1] .=> "mass")

    CSV.write("$(save_directory)/simulation_film_$(tray).csv", df)
    save_rop(sol)
end

println("Done!")


