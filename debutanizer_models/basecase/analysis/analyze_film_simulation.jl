# -*- coding: utf-8 -*-
# +
using CSV
using DataFrames
using ReactionMechanismSimulator.PyPlot
using ReactionMechanismSimulator.YAML
using ReactionMechanismSimulator.Statistics

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
font0 = Dict(
    "font.size" => 12,
)
merge!(rcParams, font0)

# model = "QCMD_cell_model"
# model = "basecase_debutanizer_model"
model = ARGS[1]

function load_film_simulations(perturb_species_list, perturb_factor_list, trays; asymptotic=false)
    film_simulations = Dict()
    for perturb_species in perturb_species_list
        for perturb_factor in perturb_factor_list
            film_simulations[perturb_species, perturb_factor] = Dict()
            if model == "basecase_debutanizer_model"
                simulation_result_folder = "../simulation_results/$(perturb_species)_$(perturb_factor)_3600.0_64.0"
            elseif model == "QCMD_cell_model"
                simulation_result_folder = "../simulation_results/$(perturb_species)_$(perturb_factor)"
            end
            for tray in trays
                if asymptotic
                    file = joinpath(simulation_result_folder, "simulation_film_$(tray)_asymptotic.csv")
                else
                    file = joinpath(simulation_result_folder, "simulation_film_$(tray).csv")
                end
                df = DataFrame(CSV.File(file))
                film_simulations[perturb_species, perturb_factor][tray] = df
            end
        end
    end
    return film_simulations
end

function load_film_rops(perturb_species_list, perturb_factor_list, trays)
    film_rops = Dict()
    for perturb_species in perturb_species_list
        for perturb_factor in perturb_factor_list
            film_rops[perturb_species, perturb_factor] = Dict()
            if model == "basecase_debutanizer_model"
                simulation_result_folder = "../simulation_results/$(perturb_species)_$(perturb_factor)_3600.0_64.0"
            elseif model == "QCMD_cell_model"
                simulation_result_folder = "../simulation_results/$(perturb_species)_$(perturb_factor)"
            end
            for tray in trays
                file = joinpath(simulation_result_folder, "simulation_film_rop_$(tray).csv")
                df = DataFrame(CSV.File(file))
                film_rops[perturb_species, perturb_factor][tray] = df
            end
        end
    end
    return film_rops
end

function calculate_film_growth_time_constant(df; use_average=false)
    if use_average
        dmdt = (df[end, "mass"] - df[1, "mass"]) / (df[end, "timestamp"] - df[1, "timestamp"])
        m = df[end, "mass"]
    else
        dmdt = (df[end, "mass"] - df[end-2, "mass"]) / (df[end, "timestamp"] - df[end-2, "timestamp"])
        m = df[end-1, "mass"]
    end
    return m / dmdt / 3600 / 24 / 365
end

function calculate_film_growth_rate(df)
    inds = 1:nrow(df)-1
    film_growth_rates = zeros(length(inds))
    for (i, ind) in enumerate(inds)
        dmdt = (df[end, "mass"] - df[ind, "mass"]) / (df[end, "timestamp"] - df[ind, "timestamp"])
        film_growth_rates[i] = dmdt
    end
    return mean(film_growth_rates), std(film_growth_rates)
end

function get_rops(df, rop_name; loss_only=false, production_only=false, N=5)
    name_inds = (df[!, "rop_spcname"] .== rop_name)
    rop_rxncomments = df[name_inds, "rop_rxncomment"]
    rop_rxncomments = [replace(rxnstr, "\n" => " ") for rxnstr in rop_rxncomments]
    rop_rxncomments = [replace(rxnstr, " H abstraction" => "") for rxnstr in rop_rxncomments]
    rop_rxnstrs = df[name_inds, "rop_rxnstr"]
    rops = df[name_inds, "rop"]
    if loss_only
        loss_inds = (rops .< 0)
        rops = abs.(rops[loss_inds])
        rop_rxncomments = rop_rxncomments[loss_inds]
        rop_rxnstrs = rop_rxnstrs[loss_inds]
    elseif production_only
        prod_inds = (rops .> 0)
        rops = rops[prod_inds]
        rop_rxncomments = rop_rxncomments[prod_inds]
        rop_rxnstrs = rop_rxnstrs[prod_inds]
    end
    sorted_inds = sortperm(abs.(rops), rev=true)
    if length(sorted_inds) > N
        sorted_inds = sorted_inds[1:N]
    end
    return rops[sorted_inds], rop_rxncomments[sorted_inds], rop_rxnstrs[sorted_inds]
end

function calculate_fragment_per_mass(name, df)
    mass = df[end, "mass"]
    return df[end, name] / mass
end

if model == "basecase_debutanizer_model"

    d = 2.5
    h = 0.3
    A = (d / 2)^2 * pi
    Vliq = A * h
    spacing = 0.6
    Vgas = A * (spacing - h)
    hfilm0 = 1e-5
    Vfilm0 = A * hfilm0
    epsilon = 0.2 #vol% of liquid in swollen film
    Vsolidinfilm0 = Vfilm0 * (1 - epsilon)
    rho = 900.0 #density of solid in swollen film
    mass = Vsolidinfilm0 * rho
    tf0 = 3600 * 24 * 365 * 200
    tf = 3600 * 24 * 365
    trays = 1:40

    perturb_species_list = ["1,3-BUTADIENE", "CYCLOPENTADIENE"]
    perturb_factor_list = ["0.5", "0.7", "0.9", "1.0", "1.1", "1.3", "1.5", "1.7", "1.9"]
    factor_num_list = [parse(Float64, perturb_factor) for perturb_factor in perturb_factor_list]

    perturb_species = "1,3-BUTADIENE"
    perturb_factor = "1.0"

    factorcmap = get_cmap(:PuRd)

    ########## plot asymptotic
    # plot fragments vs. time

    film_simulations = load_film_simulations([perturb_species], [perturb_factor], trays; asymptotic=true)

    selected_trays = trays
    selected_fragments = ["AR", "KR", "CDB", "AH"]
    nrows = length(selected_fragments)
    ncols = 1
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6, 6), sharex=true)

    cmap = get_cmap(:RdPu)
    for (fragment_ind, fragment) in enumerate(selected_fragments)
        for (tray_ind, tray) in enumerate(selected_trays)
            df = film_simulations[perturb_species, perturb_factor][tray]
            axs[fragment_ind].plot(df[!, "timestamp"] / 3600 / 24 / 365, df[!, fragment] ./ df[!, "mass"], label=tray, color=cmap(tray / length(selected_trays)))
            axs[fragment_ind].set_yscale("log")
            axs[fragment_ind].set_ylabel("$(fragment)/mass\n(mol/kg)")
            axs[fragment_ind].set_xlim([0, 100])
        end
    end

    axs[end].set_xlabel("Time (year)")

    fig.tight_layout()
    fig.savefig("basecase_debutanizer_model_film_asymptotic.pdf", bbox_inches="tight")
    plt.close()

    ########## plot real simulation
    # plot film growth rop

    film_rops = load_film_rops([perturb_species], [perturb_factor], trays)

    # plot mass rop

    selected_trays = [1, 5, 10, 15, 20, 25, 30, 35, 40]
    nrows = length(selected_trays)
    ncols = 1

    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9, 12), sharex=true)

    all_rop_rxncomments = Array{String,1}()

    for (ind, tray) in enumerate(selected_trays)
        rops, rop_rxncomments, rop_rxnstrs = get_rops(film_rops[perturb_species, perturb_factor][tray], "mass")
        append!(all_rop_rxncomments, rop_rxnstrs)
        local xs = 1:length(rop_rxncomments)
        axs[ind, 1].barh(xs, rops, align="center")
        axs[ind, 1].set_yticks(xs)
        axs[ind, 1].set_yticklabels(rop_rxncomments)
        axs[ind, 1].set_ylabel("Tray $(tray)")
        axs[ind, 1].set_xscale("log")
        axs[ind, 1].invert_yaxis()
        axs[ind, 1].set_xlim([1e-16, 1e-9])
    end

    axs[end, 1].set_xlabel("Rate of film growth (kg/s)")
    fig.tight_layout()
    fig.savefig("basecase_debutanizer_model_film_rop_mass.pdf", bbox_inches="tight")
    plt.close()

    # plot AR rop
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9, 12), sharex=true)

    for (ind, tray) in enumerate(selected_trays)
        rops, rop_rxncomments = get_rops(film_rops[perturb_species, perturb_factor][tray], "AR")
        local xs = 1:length(rop_rxncomments)
        axs[ind, 1].barh(xs, rops, align="center")
        axs[ind, 1].set_yticks(xs)
        axs[ind, 1].set_yticklabels(rop_rxncomments)
        axs[ind, 1].set_ylabel("Tray $(tray)")
        axs[ind, 1].set_xscale("symlog", linthresh=1e-16)
        axs[ind, 1].invert_yaxis()
        axs[ind, 1].set_xlim([-1e-9, 1e-9])
        axs[ind, 1].set_xticks([-1e-9, -1e-15, 1e-15, 1e-9])
        axs[ind, 1].set_xticklabels([-1e-9, -1e-15, 1e-15, 1e-9])
    end

    axs[end, 1].set_xlabel("ROP of AR (mol/s)")
    fig.tight_layout()
    fig.savefig("basecase_debutanizer_model_film_rop_AR.pdf", bbox_inches="tight")
    plt.close()

    # plot AR rop loss
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9, 12), sharex=true)

    for (ind, tray) in enumerate(selected_trays)
        rops, rop_rxncomments = get_rops(film_rops[perturb_species, perturb_factor][tray], "AR", loss_only=true)
        local xs = 1:length(rop_rxncomments)
        axs[ind, 1].barh(xs, rops, align="center")
        axs[ind, 1].set_yticks(xs)
        axs[ind, 1].set_yticklabels(rop_rxncomments)
        axs[ind, 1].set_ylabel("Tray $(tray)")
        axs[ind, 1].set_xscale("log")
        axs[ind, 1].set_xlim([1e-16, 1e-9])
        axs[ind, 1].invert_yaxis()
        axs[ind, 1].invert_xaxis()
        axs[ind, 1].set_xticks([1e-16, 1e-12, 1e-9])
        axs[ind, 1].set_xticklabels([-1e-16, -1e-12, -1e-9])
    end

    axs[end, 1].set_xlabel("Rate of AR loss (mol/s)")
    fig.tight_layout()
    fig.savefig("basecase_debutanizer_model_film_rop_AR_loss.pdf", bbox_inches="tight")
    plt.close()

    # plot KR rop
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9, 12), sharex=true)

    for (ind, tray) in enumerate(selected_trays)
        rops, rop_rxncomments = get_rops(film_rops[perturb_species, perturb_factor][tray], "KR")
        local xs = 1:length(rop_rxncomments)
        axs[ind, 1].barh(xs, rops, align="center")
        axs[ind, 1].set_yticks(xs)
        axs[ind, 1].set_yticklabels(rop_rxncomments)
        axs[ind, 1].set_ylabel("Tray $(tray)")
        axs[ind, 1].set_xscale("symlog", linthresh=1e-16)
        axs[ind, 1].invert_yaxis()
        axs[ind, 1].set_xlim([-1e-9, 1e-9])
        axs[ind, 1].set_xticks([-1e-9, -1e-15, 1e-15, 1e-9])
        axs[ind, 1].set_xticklabels([-1e-9, -1e-15, 1e-15, 1e-9])
    end

    axs[end, 1].set_xlabel("ROP of KR (mol/s)")
    fig.tight_layout()
    fig.savefig("basecase_debutanizer_model_film_rop_KR.pdf", bbox_inches="tight")
    plt.close()

    # plot KR rop loss
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9, 12), sharex=true)

    for (ind, tray) in enumerate(selected_trays)
        rops, rop_rxncomments = get_rops(film_rops[perturb_species, perturb_factor][tray], "KR", loss_only=true)
        local xs = 1:length(rop_rxncomments)
        axs[ind, 1].barh(xs, rops, align="center")
        axs[ind, 1].set_yticks(xs)
        axs[ind, 1].set_yticklabels(rop_rxncomments)
        axs[ind, 1].set_ylabel("Tray $(tray)")
        axs[ind, 1].set_xscale("log")
        axs[ind, 1].set_xlim([1e-16, 1e-9])
        axs[ind, 1].invert_yaxis()
        axs[ind, 1].invert_xaxis()
        axs[ind, 1].set_xticks([1e-16, 1e-12, 1e-9])
        axs[ind, 1].set_xticklabels([-1e-16, -1e-12, -1e-9])
    end

    axs[end, 1].set_xlabel("Rate of KR loss (mol/s)")
    fig.tight_layout()
    fig.savefig("basecase_debutanizer_model_film_rop_KR_loss.pdf", bbox_inches="tight")
    plt.close()

    # plot simulation results

    function calculate_film_chemistry_contribution(df_dict)
        chem_contributions_dict = Dict("Diels-Alder addition" => zeros(length(trays)), "radical addition" => zeros(length(trays)))
        for tray in trays
            df = df_dict[tray]
            name_inds = (df[!, "rop_spcname"] .== "mass")
            rop_rxncomments = df[name_inds, "rop_rxncomment"]
            rops = df[name_inds, "rop"]
            for (rop_rxncomment, rop) in zip(rop_rxncomments, rops)
                for family in keys(chem_contributions_dict)
                    if occursin(family, rop_rxncomment)
                        chem_contributions_dict[family][tray] += rop
                    end
                end
            end
        end
        all_rop = sum(values(chem_contributions_dict))
        normalized_chem_contributions_dict = Dict(key => value ./ all_rop .* 100.0 for (key, value) in chem_contributions_dict)
        return chem_contributions_dict, normalized_chem_contributions_dict
    end

    film_simulations = load_film_simulations([perturb_species], [perturb_factor], trays)

    s = 15
    nrows = 1
    ncols = 3
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9, 3.5), sharex=true)
    for label in ["AR", "KR", "CDB", "AH"]
        axs[1].scatter(trays, [calculate_fragment_per_mass(label, film_simulations[perturb_species, perturb_factor][tray]) for tray in trays], label=label, s=s)
    end
    axs[1].set_xlabel("Tray")
    axs[1].set_ylabel("Fragments/mass (mol/kg)")
    axs[1].set_yscale("log")
    axs[1].set_title("(a)", loc="left")
    axs[1].legend(ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.3),)

    axs[2].scatter(trays, [calculate_film_growth_time_constant(film_simulations[perturb_species, perturb_factor][tray]) for tray in trays], s=s)
    axs[2].set_xlabel("Tray")
    axs[2].set_ylabel("Film growth " * L"\tau" * " (yr)")
    axs[2].set_yscale("log")
    axs[2].set_ylim([1e0, 1e2])
    axs[2].set_title("(b)", loc="left")

    local chem_contributions_dict, normalized_chem_contributions_dict = calculate_film_chemistry_contribution(film_rops[perturb_species, perturb_factor])
    for (chem, ratios) in normalized_chem_contributions_dict
        axs[3].scatter(trays, ratios, label=chem, s=s)
    end
    axs[3].set_xlabel("Tray")
    axs[3].set_ylabel("Film growth\ncontribution (%)")
    axs[3].set_title("(c)", loc="left")
    axs[3].legend(loc="upper center", bbox_to_anchor=(0.5, -0.3))

    fig.tight_layout()
    fig.savefig("basecase_debutanizer_model_film_all.pdf", bbox_inches="tight")
    plt.close()

    # plot sensitivity to monomer perturbation

    liquid_mech = YAML.load_file("/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/liquid_mechanism/chem.rms")
    liquid_radical_names = [spc["name"] for spc in liquid_mech["Phases"][1]["Species"] if spc["radicalelectrons"] > 0]

    film_simulations = load_film_simulations(perturb_species_list, perturb_factor_list, trays)

    vapor_liquid_simulations = Dict()

    for perturb_species in perturb_species_list
        for perturb_factor in perturb_factor_list
            if model == "basecase_debutanizer_model"
                simulation_result_folder = "../simulation_results/$(perturb_species)_$(perturb_factor)_3600.0_64.0"
            elseif model == "QCMD_cell_model"
                simulation_result_folder = "../simulation_results/$(perturb_species)_$(perturb_factor)"
            end
            file = joinpath(simulation_result_folder, "simulation_vapor_liquid_yliqn_3648.0.csv")
            vapor_liquid_simulations[perturb_species, perturb_factor] = DataFrame(CSV.File(file))
        end
    end

    function calculate_liquid_radical_concentration(df, tray, liquid_radical_names)
        return sum([df[tray, radical] for radical in liquid_radical_names]) / Vliq
    end

    nrows = 4
    ncols = length(perturb_species_list)
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6, 8), sharex=true, sharey="row")
    for (species_ind, perturb_species) in enumerate(perturb_species_list)
        for (factor_ind, perturb_factor) in enumerate(perturb_factor_list)
            liquid_radical_concs = [calculate_liquid_radical_concentration(vapor_liquid_simulations[perturb_species, perturb_factor], tray, liquid_radical_names) for tray in trays]
            axs[1, species_ind].scatter(trays, liquid_radical_concs, color=factorcmap(factor_ind / length(perturb_factor_list)))
            axs[1, species_ind].set_yscale("log")
            # axs[1,species_ind].set_ylim([1e-11,1e-5])

            local film_growth_time_constants = [calculate_film_growth_time_constant(film_simulations[perturb_species, perturb_factor][tray]) for tray in trays]
            axs[2, species_ind].scatter(trays, film_growth_time_constants, color=factorcmap(factor_ind / length(perturb_factor_list)))
            axs[2, species_ind].set_yscale("log")
            # axs[2,species_ind].set_ylim([1e0,2e3])

            AR_concs = [calculate_fragment_per_mass("AR", film_simulations[perturb_species, perturb_factor][tray]) for tray in trays]
            axs[3, species_ind].scatter(trays, AR_concs, color=factorcmap(factor_ind / length(perturb_factor_list)))
            axs[3, species_ind].set_yscale("log")

            KR_concs = [calculate_fragment_per_mass("KR", film_simulations[perturb_species, perturb_factor][tray]) for tray in trays]
            axs[4, species_ind].scatter(trays, KR_concs, color=factorcmap(factor_ind / length(perturb_factor_list)))
            axs[4, species_ind].set_yscale("log")

        end
    end

    axs[1, 1].set_ylabel("R.(L) (mol/m³)")
    axs[2, 1].set_ylabel("Film growth " * L"\tau" * " (yr)")
    axs[3, 1].set_ylabel("AR/mass (mol/kg)")
    axs[4, 1].set_ylabel("KR/mass (mol/kg)")

    for (species_ind, perturb_species) in enumerate(perturb_species_list)
        axs[1, species_ind].set_title(perturb_species)
        axs[4, species_ind].set_xlabel("Trays")
    end

    cbar_ax = fig.add_axes([1.0, 0.15, 0.02, 0.7])
    sm = plt.cm.ScalarMappable(cmap=:RdPu, norm=plt.Normalize(vmin=0.5, vmax=1.9))
    cbar = fig.colorbar(sm, ticks=factor_num_list, orientation="vertical", label="Perturbation", cax=cbar_ax)

    # fig.add_subplot(111, frameon=false)
    # plt.tick_params(labelcolor="none", which="both", top=false, bottom=false, left=false, right=false)
    # plt.xlabel("Tray", fontsize=14)

    fig.tight_layout()
    fig.savefig("basecase_debutanizer_model_sens_all.pdf", bbox_inches="tight")
    plt.close()

    # plot chemistry contribution

    film_rops = Dict()
    for perturb_species in perturb_species_list
        for perturb_factor in perturb_factor_list
            film_rops[perturb_species, perturb_factor] = Dict()
            if model == "basecase_debutanizer_model"
                simulation_result_folder = "../simulation_results/$(perturb_species)_$(perturb_factor)_3600.0_64.0"
            elseif model == "QCMD_cell_model"
                simulation_result_folder = "../simulation_results/$(perturb_species)_$(perturb_factor)"
            end
            for tray in trays
                file = joinpath(simulation_result_folder, "simulation_film_rop_$(tray).csv")
                df = DataFrame(CSV.File(file))
                film_rops[perturb_species, perturb_factor][tray] = df
            end
        end
    end

    nrows = 1
    ncols = length(perturb_species_list)

    chem_colormap = Dict("Diels-Alder addition" => factorcmap, "radical addition" => get_cmap(:Greys))

    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6, 3), sharex=true, sharey=true)
    for (species_ind, perturb_species) in enumerate(perturb_species_list)
        for (factor_ind, perturb_factor) in enumerate(perturb_factor_list)
            local chem_contributions_dict, normalized_chem_contributions_dict = calculate_film_chemistry_contribution(film_rops[perturb_species, perturb_factor])
            for (chem, chem_contributions) in chem_contributions_dict
                inds = findall(x -> x > 1, chem_contributions)
                if length(inds) > 0
                    println("Chemistry contribution is too large for $(perturb_species) $(perturb_factor) $(inds) $(chem)")
                end
                inds = findall(x -> x < 0, chem_contributions)
                if length(inds) > 0
                    println("Chemistry contribution is too small for $(perturb_species) $(perturb_factor) $(inds) $(chem)")
                end
            end
            for (chem, chem_contributions) in normalized_chem_contributions_dict
                axs[species_ind].scatter(trays, chem_contributions, color=chem_colormap[chem](factor_ind / length(perturb_factor_list)))
            end
        end
        axs[species_ind].set_xlabel("Tray")
        axs[species_ind].set_title(perturb_species)
    end

    axs[1].set_ylabel("Film growth\nchemistry (%)")

    cbar_ax = fig.add_axes([1.0, 0.1, 0.02, 0.7])
    sm = plt.cm.ScalarMappable(cmap=:RdPu, norm=plt.Normalize(vmin=0.5, vmax=1.9))
    cbar = fig.colorbar(sm, ticks=factor_num_list, orientation="vertical", label="Diels-Alder addition", cax=cbar_ax)

    cbar_ax = fig.add_axes([1.13, 0.1, 0.02, 0.7])
    sm = plt.cm.ScalarMappable(cmap=:Greys, norm=plt.Normalize(vmin=0.5, vmax=1.9))
    cbar = fig.colorbar(sm, ticks=[], orientation="vertical", label="radical addition", cax=cbar_ax)

    fig.tight_layout()
    fig.savefig("basecase_debutanizer_model_sens_chem_contribution.pdf", bbox_inches="tight")
    plt.close()

elseif model == "QCMD_cell_model"

    Vreactor = 40 * 1e-9
    d = 0.55 #in
    d = d / 2.54 / 100 #m
    A = (d / 2)^2 * pi
    h = Vreactor / A
    Vliq = Vreactor
    hsolid0 = 100 * 1e-9
    Vsolidinfilm0 = A * hsolid0
    epsilon = 0.2 #vol% of liquid in swollen film
    rho = 900.0 # density of solid
    rho_liq = 970.2 # density of liquid 30vol% MCHD in Dowtherm A
    tf0 = 10000.0
    # tf=10000.0
    T = 90.0 + 273.15
    trays = 1:1

    # plot QCM-D simulation

    expt_rate_dict = Dict()
    expt_rate_dict["Unsparged"] = 9.4 #ng/hr
    expt_rate_dict["N2 sparged"] = 6.1
    expt_rate_dict["O2 sparged"] = 927.3
    cases = ["Unsparged", "N2 sparged"]

    perturb_species = "O2"
    perturb_factor_list = ["0.0", "1e-3", "1e-2", "1e-1", "1e0"]
    factor_num_list = [parse(Float64, perturb_factor) for perturb_factor in perturb_factor_list]
    tray = 1
    film_simulations = load_film_simulations([perturb_species], perturb_factor_list, trays)

    # pure_O2_sat_conc = 0.008093871600706785 #M estimated by RMG
    pure_O2_sat_molfrac = 8.1 * 1e-4 #mol fraction measured experimentally from Battino, R., Rettich, T. R., & Tominaga, T. (1983). The Solubility of Oxygen and Ozone in Liquids. Journal of Physical and Chemical Reference Data, Vol. 12, pp. 163–178. https://doi.org/10.1063/1.555680
    benzene_density = 876 #kg/m^3 from Lide, D. R., Data, S. R., Board, E. A., Baysinger, G., Chemistry, S., Library, C. E., … Zwillinger, D. (2004). CRC Handbook of Chemistry and Physics. 2660.
    benzene_mw = 78.11 / 1000 #kg/mol from Lide, D. R., Data, S. R., Board, E. A., Baysinger, G., Chemistry, S., Library, C. E., … Zwillinger, D. (2004). CRC Handbook of Chemistry and Physics. 2660.
    benzene_sat_conc = benzene_density / benzene_mw #mol/m^3
    pure_O2_sat_conc = pure_O2_sat_molfrac * benzene_sat_conc #mol/m^3
    air_O2_sat_conc = 0.21 * pure_O2_sat_conc #mol/m^3
    air_O2_sat_conc /= 1000 #mol/L

    plt.figure(figsize=(4, 3))
    xs = factor_num_list .* air_O2_sat_conc

    # simulated film growth rate
    film_growth_rates = [calculate_film_growth_rate(film_simulations[perturb_species, perturb_factor][tray]) for perturb_factor in perturb_factor_list]
    ys = [film_growth_rate[1] for film_growth_rate in film_growth_rates] * 1e9 * 1e3 * 3600 # kg/s to ng/hr
    ys .+= ys / rho * epsilon * rho_liq # converting from mass of solid to mass of film by adding mass of liquid in film
    yerr = 10.0 # a perturb_factor of x error
    label = "Prediction"
    plt.errorbar(xs, ys, yerr=[ys - ys * 1 / yerr, -ys + ys * yerr], color="C0", marker="o", capsize=5, label=label)
    plt.fill_between(xs, ys * 1 / yerr, ys * yerr, color="C0", alpha=0.5, linewidth=0)

    expt_rate_std = 1.635621377
    # unsparged
    label = "Unsparged"
    expt_rate_mean = expt_rate_dict[label] #ng/hr
    xerr = 3
    expt_o2 = air_O2_sat_conc / xerr #M
    yerr = expt_rate_std
    color = "C1"
    plt.errorbar([expt_o2], [expt_rate_mean], yerr=[yerr], xerr=[[expt_o2 - expt_o2 / xerr], [-expt_o2 + expt_o2 * xerr]], color=color, marker="o", capsize=5, label=label)
    plt.fill_between([expt_o2 / xerr, expt_o2 * xerr], [expt_rate_mean - yerr, expt_rate_mean - yerr], [expt_rate_mean + yerr, expt_rate_mean + yerr], color=color, alpha=0.5, linewidth=0)

    # N2 sparged
    label = "N2 sparged"
    expt_rate_mean = expt_rate_dict[label] #ng/hr
    xerr = 10
    reduction_factor = 10
    expt_o2 = expt_o2 / reduction_factor #M
    yerr = expt_rate_std
    color = "C2"
    plt.errorbar([expt_o2], [expt_rate_mean], yerr=[yerr], xerr=[[expt_o2 - expt_o2 / xerr], [-expt_o2 + expt_o2 * xerr]], color=color, marker="o", capsize=5, label=label)
    plt.fill_between([expt_o2 / xerr, expt_o2 * xerr], [expt_rate_mean - yerr, expt_rate_mean - yerr], [expt_rate_mean + yerr, expt_rate_mean + yerr], color=color, alpha=0.5, linewidth=0)

    plt.xlabel("O2 (M)")
    plt.ylabel("Film growth rate (ng/hr)")
    plt.xscale("symlog", linthresh=1e-6)
    plt.yscale(:log)
    # plt.ylim([1e-3, 1e2])
    plt.tight_layout()
    plt.legend(fontsize=9)
    plt.savefig("QCMD_cell_model_film_growth_rates.pdf", bbox_inches="tight")
    plt.close()


    # load film rops
    film_rops = load_film_rops([perturb_species], perturb_factor_list, trays)


    # plot mass rop
    nrows = length(perturb_factor_list)
    ncols = 1
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9, 12), sharex=true)
    for (ind, perturb_factor) in enumerate(perturb_factor_list)
        rops, rop_rxncomments, rop_rxnstrs = get_rops(film_rops[perturb_species, perturb_factor][tray], "mass")
        local xs = 1:length(rop_rxncomments)
        axs[ind, 1].barh(xs, rops, align="center")
        axs[ind, 1].set_yticks(xs)
        axs[ind, 1].set_yticklabels(rop_rxncomments)
        axs[ind, 1].set_ylabel("$(perturb_factor) M O2")
        axs[ind, 1].set_xscale("symlog", linthresh=1e-18)
        axs[ind, 1].invert_yaxis()
    end
    axs[end, 1].set_xlabel("Rate of film growth (kg/s)")
    fig.tight_layout()
    fig.savefig("QCMD_cell_model_film_rop_mass.pdf", bbox_inches="tight")
    plt.close(fig)


    # plot fragment loss
    selected_fragments = ["AR", "KR", "PR", "OR"]
    nrows = length(perturb_factor_list)
    ncols = 1
    for fragment in selected_fragments
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9, 12), sharex=true)

        for (ind, perturb_factor) in enumerate(perturb_factor_list)
            rops, rop_rxncomments = get_rops(film_rops[perturb_species, perturb_factor][tray], fragment; loss_only=true)
            local xs = 1:length(rop_rxncomments)
            axs[ind, 1].barh(xs, rops, align="center")
            axs[ind, 1].set_yticks(xs)
            axs[ind, 1].set_yticklabels(rop_rxncomments)
            axs[ind, 1].set_ylabel("$(perturb_factor) M O2")
            axs[ind, 1].set_xscale("symlog", linthresh=1e-18)
            axs[ind, 1].invert_yaxis()
        end

        axs[end, 1].set_xlabel("Rate of loss of $(fragment) (mol/s)")
        fig.tight_layout()
        fig.savefig("QCMD_cell_model_film_loss_$(fragment).pdf", bbox_inches="tight")
        plt.close(fig)
    end


    # plot fragment/mass vs. oxygen
    xs = factor_num_list .* air_O2_sat_conc
    plt.figure(figsize=(5, 4))
    for fragment in selected_fragments
        fragment_per_mass_list = [calculate_fragment_per_mass(fragment, film_simulations[perturb_species, perturb_factor][tray]) for perturb_factor in perturb_factor_list]
        plt.plot(xs, fragment_per_mass_list, marker="o", label=fragment)
    end
    plt.yscale("symlog", linthresh=1e-18)
    plt.xscale("symlog", linthresh=1e-6)
    plt.xlabel("O2 (M)")
    plt.ylabel("Fragment/mass (mol/kg)")
    plt.tight_layout()
    plt.legend()
    plt.savefig("QCMD_cell_model_film_fragment_per_mass.pdf", bbox_inches="tight")
    plt.close()

end

println("Done!")
# -


