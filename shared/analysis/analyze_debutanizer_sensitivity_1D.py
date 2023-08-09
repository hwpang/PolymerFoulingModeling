import os
import yaml
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils import extrapolate_conc

# change default font size to 12
plt.rcParams.update({"font.size": 12})


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model_name",
        type=str,
        required=True,
        help="The name of the model.",
    )
    parser.add_argument(
        "--all_simulation_directory",
        type=str,
        required=True,
        help="The path to all simulation results.",
    )

    args = parser.parse_args()
    model_name = args.model_name
    all_simulation_directory = args.all_simulation_directory

    return (
        model_name,
        all_simulation_directory,
    )


(
    model_name,
    all_simulation_directory,
) = parse_arguments()

print("model_name: ", model_name)
print("all_simulation_directory: ", all_simulation_directory)

d = 2.5
hliq = 0.3
A = (d / 2) ** 2 * np.pi
Vliq = A * hliq
rho = 900  # solid density kg/m3
epsilon = 0.2

trays = range(1, 41)
selected_trays = [1, 10, 20, 30, 40]

if model_name == "basecase_debutanizer_model":
    perturbed_species_list = ["1,3-BUTADIENE", "CYCLOPENTADIENE"]
    perturbed_factor_string_list = [
        "0.5",
        "0.7",
        "0.9",
        "1.1",
        "1.3",
        "1.5",
        "1.7",
        "1.9",
    ]
    delta_t = 64.0
    reaction_types = ["Diels-Alder addition", "radical addition"]
else:
    perturbed_species_list = ["OXYGEN"]
    perturbed_factor_string_list = [
        "0.0",
        "1e-6",
        "1e-5",
        "1e-4",
        "1e-3",
        "1e-2",
        "1e-1",
        "1e0",
    ]
    delta_t = 32.0
    # reaction_types = [
    #     "carbon-center radical addition",
    #     "oxygen-center radical addition",
    # ]
    reaction_types = ["Diels-Alder addition", "radical addition"]
    num_cells = 3
    # num_cells = 1
    cell_inds = range(1, num_cells + 1)
    # method = "same_initial_size"
    method = "callback"

perturbed_factor_num_list = [
    float(perturbed_factor) for perturbed_factor in perturbed_factor_string_list
]
tf_liq = np.arange(0, 3600.0 + delta_t, delta_t)[-1]

print("Loading all liquid simulation results...")
liquid_simulations = dict()
for perturbed_species in perturbed_species_list:
    for perturbed_factor in perturbed_factor_string_list:
        simulation_path = os.path.join(
            all_simulation_directory,
            f"{perturbed_species}_{perturbed_factor}_3600.0_{delta_t}",
            f"simulation_vapor_liquid_yliqn_{tf_liq}.csv",
        )
        liquid_simulations[perturbed_species, perturbed_factor] = pd.read_csv(
            simulation_path
        )

print("Loading all film simulation results...")

film_simulations = dict()
for perturbed_species in perturbed_species_list:
    for perturbed_factor in perturbed_factor_string_list:
        for tray in trays:
            simulation_path = os.path.join(
                all_simulation_directory,
                f"{perturbed_species}_{perturbed_factor}_3600.0_{delta_t}",
                f"simulation_film_1D_{tray}_{num_cells}cells_{method}.csv",
            )
            film_simulations[perturbed_species, perturbed_factor, tray] = pd.read_csv(
                simulation_path
            )

print("Loading radical labels...")

alpha_rates_path = os.path.join(
    all_simulation_directory,
    f"{perturbed_species_list[0]}_{perturbed_factor_string_list[0]}_3600.0_{delta_t}",
    "alpha_rates.yml",
)

with open(alpha_rates_path, "r") as f:
    results = yaml.load(f, Loader=yaml.FullLoader)
    all_alphas, consumption_rates, production_rates, rxn_rates, radical_labels = results
    alpha1, alpha2, alphas, alphas_DA = all_alphas
    R_labels, ROO_labels, RO_labels = radical_labels

carbon_center_radical_names = [label + "(L)" for label in R_labels]
peroxyl_radical_names = [label + "(L)" for label in ROO_labels]
# define helper functions


def calculate_film_growth_time_constant(df):
    n_rows = len(df.index)
    m = np.sum(
        [df.loc[n_rows - 2, f"mass_cell_{cell_ind}"] for cell_ind in cell_inds], axis=0
    )
    dmdt = (
        np.sum(
            [df.loc[n_rows - 1, f"mass_cell_{cell_ind}"] for cell_ind in cell_inds],
            axis=0,
        )
        - np.sum(
            [df.loc[n_rows - 3, f"mass_cell_{cell_ind}"] for cell_ind in cell_inds],
            axis=0,
        )
    ) / (df.loc[n_rows - 1, "timestamp"] - df.loc[n_rows - 3, "timestamp"])
    return m / dmdt / 3600.0 / 24.0 / 365.0  # year


def calculate_fragment_per_mass(df, label):
    n_rows = len(df.index)
    if df.loc[n_rows - 1, "timestamp"] < 3e7:
        return None
    return np.sum(
        [df.loc[n_rows - 1, f"{label}_cell_{cell_ind}"] for cell_ind in cell_inds]
    ) / np.sum([df.loc[n_rows - 1, f"mass_cell_{cell_ind}"] for cell_ind in cell_inds])


def calculate_fouling_chemistry_contribution(film_rate_of_productions, normalize=True):
    fouling_chemistry_contribution = dict()
    for perturbed_species in perturbed_species_list:
        for perturbed_factor in perturbed_factor_string_list:
            fouling_chemistry_contribution[(perturbed_species, perturbed_factor)] = {
                reaction_type: np.zeros(len(trays)) for reaction_type in reaction_types
            }
            for tray in trays:
                df = film_rate_of_productions[perturbed_species, perturbed_factor, tray]
                name_inds = df["rop_spcname"] == "mass"
                rop_rxncomments = df.loc[name_inds, "rop_rxncomment"]
                rops = df.loc[name_inds, "rop"]

                for rop_rxncomment, rop in zip(rop_rxncomments, rops):
                    for reaction_type in reaction_types:
                        # if model_name == "basecase_debutanizer_model":
                        if reaction_type in rop_rxncomment:
                            fouling_chemistry_contribution[
                                perturbed_species, perturbed_factor
                            ][reaction_type][tray - 1] += rop
                    # elif model_name == "trace_oxygen_perturbed_debutanizer_model":
                    #     if "PR" in rop_rxncomment or "OR" in rop_rxncomment:
                    #         fouling_chemistry_contribution[
                    #             perturbed_species, perturbed_factor
                    #         ]["oxygen-center radical addition"][tray - 1] += rop
                    #     elif "AR" in rop_rxncomment or "KR" in rop_rxncomment:
                    #         fouling_chemistry_contribution[
                    #             perturbed_species, perturbed_factor
                    #         ]["carbon-center radical addition"][tray - 1] += rop

    if normalize:
        for perturbed_species in perturbed_species_list:
            for perturbed_factor in perturbed_factor_string_list:
                rop_sum = np.sum(
                    [
                        fouling_chemistry_contribution[
                            perturbed_species, perturbed_factor
                        ][reaction_type]
                        for reaction_type in reaction_types
                    ],
                    axis=0,
                )
                for reaction_type in reaction_types:
                    fouling_chemistry_contribution[perturbed_species, perturbed_factor][
                        reaction_type
                    ] /= rop_sum
                    fouling_chemistry_contribution[perturbed_species, perturbed_factor][
                        reaction_type
                    ] *= 100.0

    return fouling_chemistry_contribution


factorcmap = plt.get_cmap("PuRd")


def get_from_1D_simulation(df, name):
    n_rows = len(df.index)
    if df.loc[n_rows - 1, "timestamp"] < 3e7:
        return None
    return np.sum(
        [df.loc[n_rows - 1, f"{name}_cell_{cell_ind}"] for cell_ind in cell_inds],
        axis=0,
    )


print("Plotting h vs. t for selected trays and all sensitivity simulations...")

nrows = len(selected_trays)
ncols = 1

min_thickness = 1e10
max_thickness = 0.0

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6, 10), sharex=True)
for i, tray in enumerate(selected_trays):
    ax = axs[i]
    for perturbed_species in perturbed_species_list:
        for perturbed_factor_ind, perturbed_factor in enumerate(
            perturbed_factor_string_list
        ):
            film_simulation = film_simulations[
                perturbed_species, perturbed_factor, tray
            ]
            ts = np.array(film_simulation["timestamp"]) / 3600 / 24 / 365
            hs = (
                np.sum(
                    [
                        film_simulation[f"mass_cell_{cell_ind}"]
                        for cell_ind in cell_inds
                    ],
                    axis=0,
                )
                / rho
                / (1 - epsilon)
                / A
            )
            ax.plot(
                ts,
                hs,
                color=factorcmap(
                    perturbed_factor_ind / len(perturbed_factor_string_list)
                ),
            )
            min_thickness = min(min_thickness, np.min(hs))
            max_thickness = max(max_thickness, np.max(hs))
            ax.set_yscale("log")
            ax.set_ylabel(f"Tray {tray}")

axs[-1].set_xlabel("Time (yr)")
for ax in axs.flat:
    ax.set_ylim([min_thickness, max_thickness])

fig.align_labels()

fig.add_subplot(111, frameon=False)
plt.tick_params(
    labelcolor="none", which="both", top=False, bottom=False, left=False, right=False
)
plt.ylabel("Film thickness (m)", x=-1.0)

cbar_ax = fig.add_axes([1.0, 0.15, 0.02, 0.7])
if model_name == "basecase_debutanizer_model":
    sm = plt.cm.ScalarMappable(cmap="RdPu", norm=plt.Normalize(vmin=0.5, vmax=1.9))
else:
    sm = plt.cm.ScalarMappable(
        cmap="RdPu", norm=plt.matplotlib.colors.LogNorm(vmin=1e-7, vmax=1e0)
    )
cbar = fig.colorbar(
    sm,
    ticks=perturbed_factor_num_list,
    orientation="vertical",
    label="Perturbation",
    cax=cbar_ax,
)

fig.tight_layout()
fig.savefig(
    f"Figures/{model_name}_1D_sens_film_thickness_vs_t.pdf",
    bbox_inches="tight",
)

print("Plotting all sensitivity simulation results...")

# normalized_fouling_chemistry_contribution = calculate_fouling_chemistry_contribution(
#     film_rate_of_productions, normalize=True
# )
# print(normalized_fouling_chemistry_contribution)

# markers = ["x", "v", "D", "P", "X", "8", "p", "*", "h", "H", "d"]

if model_name == "basecase_debutanizer_model":
    nrows = 5
    ncols = 2
    fig, axs = plt.subplots(
        nrows=nrows, ncols=ncols, figsize=(6, 10), sharex=True, sharey="row"
    )
elif model_name == "trace_oxygen_perturbed_debutanizer_model":
    nrows = 4
    ncols = 2
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6, 8), sharex=True)

if model_name == "basecase_debutanizer_model":
    for species_ind, perturbed_species in enumerate(perturbed_species_list):
        for factor_ind, perturbed_factor in enumerate(perturbed_factor_string_list):
            ax = axs[0, species_ind]
            liquid_carbon_center_radical_concs = (
                liquid_simulations[perturbed_species, perturbed_factor]
                .loc[:, R_labels]
                .sum(axis=1)
                / Vliq
            )
            ax.plot(
                trays,
                liquid_carbon_center_radical_concs,
                "-o",
                color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
            )
            if species_ind == 0:
                ax.set_ylabel("R.(L) (mol/m³)")
            ax.set_yscale("log")
            ax.set_title(perturbed_species)

            ax = axs[1, species_ind]
            film_growth_time_constants = [
                calculate_film_growth_time_constant(
                    film_simulations[perturbed_species, perturbed_factor, tray]
                )
                for tray in trays
            ]
            ax.plot(
                trays,
                film_growth_time_constants,
                "-o",
                color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
            )
            if species_ind == 0:
                ax.set_ylabel("Film growth τ (year)")
            ax.set_yscale("log")

            ax = axs[2, species_ind]
            film_fragment_per_mass = [
                calculate_fragment_per_mass(
                    film_simulations[perturbed_species, perturbed_factor, tray],
                    "AR",
                )
                for tray in trays
            ]
            ax.plot(
                trays,
                film_fragment_per_mass,
                "-o",
                color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
            )
            if species_ind == 0:
                ax.set_ylabel("AR (mol/kg)")
            ax.set_yscale("log")

            ax = axs[3, species_ind]
            film_fragment_per_mass = [
                calculate_fragment_per_mass(
                    film_simulations[perturbed_species, perturbed_factor, tray],
                    "KR",
                )
                for tray in trays
            ]
            ax.plot(
                trays,
                film_fragment_per_mass,
                "-o",
                color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
            )
            if species_ind == 0:
                ax.set_ylabel("KR (mol/kg)")
            ax.set_yscale("log")

            ax = axs[4, species_ind]
            contribution = normalized_fouling_chemistry_contribution[
                perturbed_species, perturbed_factor
            ]
            label = None
            for ind, reaction_type in enumerate(reaction_types):
                if (
                    perturbed_factor == perturbed_factor_string_list[-1]
                    and perturbed_species == perturbed_species_list[-1]
                ):
                    label = reaction_type
                ax.scatter(
                    trays,
                    contribution[reaction_type],
                    # color=plt.get_cmap(cmaps[ind])(factor_ind / len(perturbed_factor_string_list)),
                    color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
                    marker=markers[ind],
                    edgecolors="black",
                    label=label,
                    zorder=ind,
                )
            if species_ind == 0:
                ax.set_ylabel("Importance of\nfilm growth pathways (%)")

            if (
                perturbed_factor == perturbed_factor_string_list[-1]
                and perturbed_species == perturbed_species_list[-1]
            ):
                ax.legend(bbox_to_anchor=(1.05, -0.3), loc="upper right")

else:
    for species_ind, perturbed_species in enumerate(perturbed_species_list):
        for factor_ind, perturbed_factor in enumerate(perturbed_factor_string_list):
            ax = axs[0, 0]
            liquid_carbon_center_radical_concs = (
                liquid_simulations[perturbed_species, perturbed_factor]
                .loc[:, R_labels]
                .sum(axis=1)
                / Vliq
            )
            ax.plot(
                trays,
                liquid_carbon_center_radical_concs,
                "-o",
                color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
            )
            ax.set_ylabel("RC.(L,bulk) (mol/m³)")
            ax.set_yscale("log")
            ax.set_title("(a)", loc="left")

            if perturbed_factor != "0.0":
                ax = axs[0, 1]
                liquid_peroxyl_radical_concs = (
                    liquid_simulations[perturbed_species, perturbed_factor]
                    .loc[:, ROO_labels]
                    .sum(axis=1)
                    / Vliq
                )
                ax.plot(
                    trays,
                    liquid_peroxyl_radical_concs,
                    "-o",
                    color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
                )
                ax.set_ylabel("ROO.(L,bulk) (mol/m³)")
                ax.set_yscale("log")
                ax.set_title("(b)", loc="left")

                ax = axs[1, 0]
                liquid_alkoxyl_radical_concs = (
                    liquid_simulations[perturbed_species, perturbed_factor]
                    .loc[:, RO_labels]
                    .sum(axis=1)
                    / Vliq
                )
                ax.plot(
                    trays,
                    liquid_alkoxyl_radical_concs,
                    "-o",
                    color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
                )
                ax.set_ylabel("RO.(L,bulk) (mol/m³)")
                ax.set_yscale("log")
                ax.set_title("(c)", loc="left")

            ax = axs[1, 1]
            ms = [
                get_from_1D_simulation(
                    film_simulations[perturbed_species, perturbed_factor, tray], "mass"
                )
                for tray in trays
            ]
            trays_plot = [trays[ind] for ind, m in enumerate(ms) if m is not None]
            hs_plot = (
                np.array([ms[ind] for ind, m in enumerate(ms) if m is not None])
                / rho
                / (1 - epsilon)
                / A
            )
            ax.plot(
                trays_plot,
                hs_plot,
                "-o",
                color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
            )
            ax.set_ylabel("Film thickness (m)")
            ax.set_yscale("log")
            ax.set_title("(d)", loc="left")

            ax = axs[2, 0]
            concs = [
                calculate_fragment_per_mass(
                    film_simulations[perturbed_species, perturbed_factor, tray], "AR"
                )
                for tray in trays
            ]
            trays_plot = [
                trays[ind] for ind, conc in enumerate(concs) if conc is not None
            ]
            concs_plot = [
                concs[ind] for ind, conc in enumerate(concs) if conc is not None
            ]
            ax.plot(
                trays_plot,
                concs_plot,
                "-o",
                color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
            )
            ax.set_ylabel("AR/mass (mol/kg)")
            ax.set_yscale("log")
            ax.set_title("(e)", loc="left")

            ax = axs[2, 1]
            concs = [
                calculate_fragment_per_mass(
                    film_simulations[perturbed_species, perturbed_factor, tray], "KR"
                )
                for tray in trays
            ]
            trays_plot = [
                trays[ind] for ind, conc in enumerate(concs) if conc is not None
            ]
            concs_plot = [
                concs[ind] for ind, conc in enumerate(concs) if conc is not None
            ]
            ax.plot(
                trays_plot,
                concs_plot,
                "-o",
                color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
            )
            ax.set_ylabel("KR/mass (mol/kg)")
            ax.set_yscale("log")
            ax.set_title("(f)", loc="left")

            if perturbed_factor != "0.0":
                ax = axs[3, 0]
                concs = [
                    calculate_fragment_per_mass(
                        film_simulations[perturbed_species, perturbed_factor, tray],
                        "PR",
                    )
                    for tray in trays
                ]
                trays_plot = [
                    trays[ind] for ind, conc in enumerate(concs) if conc is not None
                ]
                concs_plot = [
                    concs[ind] for ind, conc in enumerate(concs) if conc is not None
                ]
                ax.plot(
                    trays_plot,
                    concs_plot,
                    "-o",
                    color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
                )
                ax.set_ylabel("PR/mass (mol/kg)")
                ax.set_yscale("log")
                ax.set_title("(g)", loc="left")

                ax = axs[3, 1]
                concs = [
                    calculate_fragment_per_mass(
                        film_simulations[perturbed_species, perturbed_factor, tray],
                        "OR",
                    )
                    for tray in trays
                ]
                trays_plot = [
                    trays[ind] for ind, conc in enumerate(concs) if conc is not None
                ]
                concs_plot = [
                    concs[ind] for ind, conc in enumerate(concs) if conc is not None
                ]
                ax.plot(
                    trays_plot,
                    concs_plot,
                    "-o",
                    color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
                )
                ax.set_ylabel("OR/mass (mol/kg)")
                ax.set_yscale("log")
                ax.set_title("(h)", loc="left")


axs[-1, 0].set_xlabel("Trays")
axs[-1, 1].set_xlabel("Trays")

cbar_ax = fig.add_axes([1.0, 0.15, 0.02, 0.7])
if model_name == "basecase_debutanizer_model":
    sm = plt.cm.ScalarMappable(cmap="RdPu", norm=plt.Normalize(vmin=0.5, vmax=1.9))
else:
    sm = plt.cm.ScalarMappable(
        cmap="RdPu", norm=plt.matplotlib.colors.LogNorm(vmin=1e-7, vmax=1e0)
    )
cbar = fig.colorbar(
    sm,
    ticks=perturbed_factor_num_list,
    orientation="vertical",
    label="Perturbation",
    cax=cbar_ax,
)
fig.align_ylabels(axs[:, 0])
fig.align_ylabels(axs[:, 1])
fig.tight_layout()
fig.savefig(f"Figures/{model_name}_1D_sens_all.pdf", bbox_inches="tight")
plt.close()

print("Plot h_i, [O2_i], ... vs. z")

for tray in selected_trays:
    min_thickness = 1e10
    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(8, 9), sharex=True)
    for perturbed_species in perturbed_species_list:
        for factor_ind, perturbed_factor in enumerate(perturbed_factor_string_list):
            simulation = film_simulations[perturbed_species, perturbed_factor, tray]
            row_ind = len(simulation.index) - 1

            masshats = np.array(
                [
                    simulation.loc[row_ind, f"mass_cell_{cell_ind}"]
                    for cell_ind in cell_inds
                ]
            )
            masshats_init = np.array(
                [simulation.loc[0, f"mass_cell_{cell_ind}"] for cell_ind in cell_inds]
            )
            Vhatsolidinfilms = masshats / rho
            Vhatliqinfilms = np.array(
                [
                    simulation.loc[row_ind, f"Vliqinfilm_cell_{cell_ind}"]
                    for cell_ind in cell_inds
                ]
            )
            Vhatliqinfilms_init = np.array(
                [
                    simulation.loc[0, f"Vliqinfilm_cell_{cell_ind}"]
                    for cell_ind in cell_inds
                ]
            )
            Vhatfilms = Vhatsolidinfilms / (1 - epsilon)
            t = simulation.loc[row_ind, "timestamp"]
            zs = list(np.cumsum(Vhatfilms / A) - Vhatfilms / A / 2)
            zs.insert(0, 0.0)
            zs.append(zs[-1] + Vhatfilms[-1] / A / 2)
            zs = np.array(zs)
            min_thickness = min(min_thickness, zs[-1])

            # zs /= zs[-1]  # normalized with thickness of film

            if perturbed_factor != "0.0":
                ax = axs[0, 0]
                concs = list(
                    np.array(
                        [
                            simulation.loc[row_ind, f"OXYGEN(L)_cell_{cell_ind}"]
                            for cell_ind in cell_inds
                        ]
                    )
                    / Vhatliqinfilms
                )
                extrapolated_conc = extrapolate_conc(concs, zs)
                concs.insert(0, extrapolated_conc)
                concs.append(concs[-1])
                ax.plot(
                    zs,
                    concs,
                    "-o",
                    color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
                )
                ax.set_ylabel("[O2(L,film)] (mol/m^3)")
                if max(concs) / min(concs) > 10:
                    ax.set_yscale("log")
                ax.set_title("(a)", loc="left", y=1.08)

            ax = axs[1, 0]
            concs = list(
                np.array(
                    [
                        np.sum(
                            [
                                simulation.loc[row_ind, f"{name}_cell_{cell_ind}"]
                                for name in carbon_center_radical_names
                            ]
                        )
                        for cell_ind in cell_inds
                    ]
                )
                / Vhatliqinfilms
            )

            extrapolated_conc = extrapolate_conc(concs, zs)
            concs.insert(0, extrapolated_conc)
            concs.append(concs[-1])
            ax.plot(
                zs,
                concs,
                "-o",
                color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
            )
            ax.set_ylabel("[RC.(L,film)] (mol/m^3)")
            if max(concs) / min(concs) > 10:
                ax.set_yscale("log")
            ax.set_title("(b)", loc="left", y=1.08)

            if perturbed_factor != "0.0":
                ax = axs[2, 0]
                concs = list(
                    np.array(
                        [
                            np.sum(
                                [
                                    simulation.loc[row_ind, f"{name}_cell_{cell_ind}"]
                                    for name in peroxyl_radical_names
                                ]
                            )
                            for cell_ind in cell_inds
                        ]
                    )
                    / Vhatliqinfilms
                )

                extrapolated_conc = extrapolate_conc(concs, zs)
                concs.insert(0, extrapolated_conc)
                concs.append(concs[-1])
                ax.plot(
                    zs,
                    concs,
                    "-o",
                    color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
                )
                ax.set_ylabel("[ROO.(L,film)] (mol/m^3)")
                if max(concs) / min(concs) > 10:
                    ax.set_yscale("log")
                ax.set_title("(c)", loc="left", y=1.08)
                # ax.set_xlabel("z/h (-)")
                ax.set_xlabel("z (m)")

            ax = axs[0, 1]
            concs = list(
                np.array(
                    [
                        simulation.loc[row_ind, f"AR_cell_{cell_ind}"]
                        for cell_ind in cell_inds
                    ]
                )
                / masshats
            )
            extrapolated_conc = extrapolate_conc(concs, zs)
            concs.insert(0, extrapolated_conc)
            concs.append(concs[-1])
            ax.plot(
                zs,
                concs,
                "-o",
                color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
            )
            ax.set_ylabel("[AR] (mol/kg)")
            if max(concs) / min(concs) > 10:
                ax.set_yscale("log")
            ax.set_title("(d)", loc="left", y=1.08)

            ax = axs[1, 1]
            concs = list(
                np.array(
                    [
                        simulation.loc[row_ind, f"KR_cell_{cell_ind}"]
                        for cell_ind in cell_inds
                    ]
                )
                / masshats
            )
            extrapolated_conc = extrapolate_conc(concs, zs)
            concs.insert(0, extrapolated_conc)
            concs.append(concs[-1])
            ax.plot(
                zs,
                concs,
                "-o",
                color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
            )
            ax.set_ylabel("[KR] (mol/kg)")
            if max(concs) / min(concs) > 10:
                ax.set_yscale("log")
            ax.set_title("(e)", loc="left", y=1.08)

            if perturbed_factor != "0.0":
                ax = axs[2, 1]
                concs = list(
                    np.array(
                        [
                            simulation.loc[row_ind, f"PR_cell_{cell_ind}"]
                            for cell_ind in cell_inds
                        ]
                    )
                    / masshats
                )
                extrapolated_conc = extrapolate_conc(concs, zs)
                concs.insert(0, extrapolated_conc)
                concs.append(concs[-1])
                ax.plot(
                    zs,
                    concs,
                    "-o",
                    color=factorcmap(factor_ind / len(perturbed_factor_string_list)),
                )
                ax.set_ylabel("[PR] (mol/kg)")
                if max(concs) / min(concs) > 10:
                    ax.set_yscale("log")
                ax.set_title("(f)", loc="left", y=1.08)
                # ax.set_xlabel("z/h (-)")
                ax.set_xlabel("z (m)")

    # for ax in axs.flatten():
    #     ax.set_xlim([0, min_thickness])
    # axs[2, 1].legend(bbox_to_anchor=(0.0, -0.25), loc="upper left")

    cbar_ax = fig.add_axes([1.0, 0.15, 0.02, 0.7])
    if model_name == "basecase_debutanizer_model":
        sm = plt.cm.ScalarMappable(cmap="RdPu", norm=plt.Normalize(vmin=0.5, vmax=1.9))
    else:
        sm = plt.cm.ScalarMappable(
            cmap="RdPu", norm=plt.matplotlib.colors.LogNorm(vmin=1e-7, vmax=1e0)
        )
    cbar = fig.colorbar(
        sm,
        ticks=perturbed_factor_num_list,
        orientation="vertical",
        label="Perturbation",
        cax=cbar_ax,
    )
    fig.align_ylabels(axs[:, 0])
    fig.align_ylabels(axs[:, 1])
    fig.tight_layout()
    fig.savefig(
        f"Figures/{model_name}_{tray}_sens_all_vs_z.pdf",
        bbox_inches="tight",
    )
