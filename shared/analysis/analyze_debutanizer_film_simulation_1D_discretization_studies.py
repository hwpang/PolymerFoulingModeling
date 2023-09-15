import os
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
        "--simulation_directory",
        type=str,
        required=True,
        help="The path to the film simulation results.",
    )

    args = parser.parse_args()
    model_name = args.model_name
    simulation_directory = args.simulation_directory

    return (
        model_name,
        simulation_directory,
    )


(
    model_name,
    simulation_directory,
) = parse_arguments()

print("model_name: ", model_name)
print("simulation_directory: ", simulation_directory)

rho = 900.0
d = 2.5
h = 0.3
A = (d / 2) ** 2 * np.pi
Vliq = A * h
spacing = 0.6
Vgas = A * (spacing - h)
hfilm0 = 1e-5
epsilon = 0.2  # vol% of liquid in swollen film
rho = 900.0  # density of solid
Vfilm0 = A * hfilm0
Vsolidinfilm0 = Vfilm0 * (1 - epsilon)
Vliqinfilm0 = Vfilm0 * epsilon

print("Load 1-D film simulation results")
num_cells_list = [
    1,
    2,
    3,
    4,
    # 5,
]
tray = 17
perturbed_factor_string_list = [
    "1e0",
    "1e-1",
    "1e-2",
    "1e-3",
    "1e-4",
    "1e-5",
    "1e-6",
]
all_discretization_methods = [
    "same_initial_size",
    "different_initial_size",
    "callback",
    "different_initial_size_and_callback",
]
discretization_methods = [
    # "same_initial_size",
    # "different_initial_size",
    "callback",
    # "different_initial_size_and_callback",
]
markers = ["o", "s", "^", "v"]
discretization_method_markers = dict(zip(all_discretization_methods, markers))
colors = [
    "tab:blue",
    "tab:orange",
    "tab:green",
    "tab:red",
    "tab:purple",
    "tab:brown",
    "tab:pink",
]
num_cells_colors = dict(zip(num_cells_list, colors))
linestyles = ["-", "--", "-.", ":"]
discretization_method_linestyles = dict(zip(all_discretization_methods, linestyles))

one_d_simulations = dict()
for perturbed_factor_string in perturbed_factor_string_list:
    for num_cells in num_cells_list:
        for discretization_method in discretization_methods:
            one_d_simulation_path = os.path.join(
                simulation_directory,
                f"simulation_film_1D_{perturbed_factor_string}_{tray}_{num_cells}cells_{discretization_method}.csv",
            )
            one_d_simulations[
                perturbed_factor_string, num_cells, discretization_method
            ] = pd.read_csv(one_d_simulation_path)

save_name = "discretization_studies"


def get_film_growth_time_constants(hfilms, ts):
    dhdt = (hfilms[2:] - hfilms[:-2]) / (ts[2:] - ts[:-2])
    h = hfilms[1:-1]
    return h / dhdt / 3600 / 24 / 365


print("Plot h, O2, ... vs. t")

for perturbed_factor_string in perturbed_factor_string_list:
    fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(8, 9), sharex=True)

    num_cells_fine = num_cells_list[-1]
    discretization_method_fine = discretization_methods[-1]
    cell_inds_fine = range(1, num_cells_fine + 1)
    simulation = one_d_simulations[
        perturbed_factor_string, num_cells_fine, discretization_method_fine
    ]
    hfilms_fine = (
        np.sum(
            [simulation[f"mass_cell_{cell_ind}"] for cell_ind in cell_inds_fine], axis=0
        )
        / rho
        / A
        / (1 - epsilon)
    )
    film_growth_time_constants_fine = get_film_growth_time_constants(
        hfilms_fine, np.array(simulation["timestamp"])
    )

    for num_cells in num_cells_list:
        cell_inds = range(1, num_cells + 1)

        for discretization_method in discretization_methods:
            simulation = one_d_simulations[
                perturbed_factor_string, num_cells, discretization_method
            ]

            # ax = axs[0, 0]
            # hfilms = (
            #     np.sum(
            #         [simulation[f"mass_cell_{cell_ind}"] for cell_ind in cell_inds],
            #         axis=0,
            #     )
            #     / rho
            #     / A
            #     / (1 - epsilon)
            # )
            # label = f"{num_cells} cells ({discretization_method}, diff={(hfilms[-1] - hfilms_fine[-1])/hfilms_fine[-1]*100:.0f} %)"
            # ax.plot(
            #     simulation["timestamp"],
            #     hfilms,
            #     label=label,
            #     linestyle=discretization_method_linestyles[discretization_method],
            #     color=num_cells_colors[num_cells],
            # )
            # if max(hfilms) / min(hfilms) > 10:
            #     ax.set_yscale("log")
            # ax.set_ylabel("Film thickness (m)")
            # ax.set_title("(a)", loc="left", y=1.08)

            ax = axs[0, 0]
            ts = np.array(simulation["timestamp"])
            hfilms = (
                np.sum(
                    [simulation[f"mass_cell_{cell_ind}"] for cell_ind in cell_inds],
                    axis=0,
                )
                / rho
                / A
                / (1 - epsilon)
            )
            film_growth_time_constants = get_film_growth_time_constants(hfilms, ts)
            # label = f"{num_cells} cells ({discretization_method}, diff={(film_growth_time_constants[-1] - film_growth_time_constants_fine[-1])/film_growth_time_constants_fine[-1]*100:.0f} %)"
            label = f"{num_cells} cells ({discretization_method})"
            ax.plot(
                ts[1:-1],
                film_growth_time_constants,
                label=label,
                linestyle=discretization_method_linestyles[discretization_method],
                color=num_cells_colors[num_cells],
            )
            if max(film_growth_time_constants) / min(film_growth_time_constants) > 10:
                ax.set_yscale("log")
            ax.set_ylabel(r"$\tau$ (year)")
            ax.set_title("(a)", loc="left", y=1.08)

            ax = axs[1, 0]
            concs = np.sum(
                [simulation[f"OXYGEN(L)_cell_{cell_ind}"] for cell_ind in cell_inds],
                axis=0,
            ) / np.sum(
                [simulation[f"Vliqinfilm_cell_{cell_ind}"] for cell_ind in cell_inds],
                axis=0,
            )
            ax.plot(
                simulation["timestamp"],
                concs,
                linestyle=discretization_method_linestyles[discretization_method],
                color=num_cells_colors[num_cells],
            )
            if max(concs) / min(concs) > 10:
                ax.set_yscale("log")
            ax.set_ylabel("[O2(L,film)] (mol/m^3)")
            # ax.set_ylabel("O2(L) (mol)")
            ax.set_title("(b)", loc="left", y=1.08)

            species_names = [
                name.replace("_cell_1", "")
                for name in simulation.columns.values
                if "_cell_1" in name
            ]
            carbon_center_radical_names = [
                name for name in species_names if "[C" in name or "[c" in name
            ]
            oxygen_center_radical_names = [
                name for name in species_names if "[O" in name or "[o" in name
            ]

            ax = axs[2, 0]
            concs = np.sum(
                [
                    simulation[f"{name}_cell_{cell_ind}"]
                    for name in carbon_center_radical_names
                    for cell_ind in cell_inds
                ],
                axis=0,
            ) / np.sum(
                [simulation[f"Vliqinfilm_cell_{cell_ind}"] for cell_ind in cell_inds],
                axis=0,
            )
            ax.plot(
                simulation["timestamp"],
                concs,
                linestyle=discretization_method_linestyles[discretization_method],
                color=num_cells_colors[num_cells],
            )
            if max(concs) / min(concs) > 10:
                ax.set_yscale("log")
            ax.set_ylabel("[RC.(L,film)] (mol/m^3)")
            ax.set_title("(c)", loc="left", y=1.08)

            ax = axs[3, 0]
            concs = np.sum(
                [
                    simulation[f"{name}_cell_{cell_ind}"]
                    for name in oxygen_center_radical_names
                    for cell_ind in cell_inds
                ],
                axis=0,
            ) / np.sum(
                [simulation[f"Vliqinfilm_cell_{cell_ind}"] for cell_ind in cell_inds],
                axis=0,
            )
            ax.plot(
                simulation["timestamp"],
                concs,
                linestyle=discretization_method_linestyles[discretization_method],
                color=num_cells_colors[num_cells],
            )
            if max(concs) / min(concs) > 10:
                ax.set_yscale("log")
            ax.set_ylabel("[ROO.(L,film)] (mol/m^3)")
            ax.set_title("(d)", loc="left", y=1.08)
            ax.set_xlabel("Time (s)")

            ax = axs[1, 1]
            concs = np.sum(
                [simulation[f"AR_cell_{cell_ind}"] for cell_ind in cell_inds], axis=0
            ) / np.sum(
                [simulation[f"mass_cell_{cell_ind}"] for cell_ind in cell_inds], axis=0
            )
            ax.plot(
                simulation["timestamp"],
                concs,
                linestyle=discretization_method_linestyles[discretization_method],
                color=num_cells_colors[num_cells],
            )
            if max(concs) / min(concs) > 10:
                ax.set_yscale("log")
            ax.set_ylabel("[AR] (mol/kg)")
            ax.set_title("(e)", loc="left", y=1.08)

            ax = axs[2, 1]
            concs = np.sum(
                [simulation[f"KR_cell_{cell_ind}"] for cell_ind in cell_inds], axis=0
            ) / np.sum(
                [simulation[f"mass_cell_{cell_ind}"] for cell_ind in cell_inds], axis=0
            )
            ax.plot(
                simulation["timestamp"],
                concs,
                linestyle=discretization_method_linestyles[discretization_method],
                color=num_cells_colors[num_cells],
            )
            if max(concs) / min(concs) > 10:
                ax.set_yscale("log")
            ax.set_ylabel("[KR] (mol/kg)")
            ax.set_title("(f)", loc="left", y=1.08)

            ax = axs[3, 1]
            concs = np.sum(
                [simulation[f"PR_cell_{cell_ind}"] for cell_ind in cell_inds], axis=0
            ) / np.sum(
                [simulation[f"mass_cell_{cell_ind}"] for cell_ind in cell_inds], axis=0
            )
            ax.plot(
                simulation["timestamp"],
                concs,
                linestyle=discretization_method_linestyles[discretization_method],
                color=num_cells_colors[num_cells],
            )
            if max(concs) / min(concs) > 10:
                ax.set_yscale("log")
            ax.set_ylabel("[PR] (mol/kg)")
            ax.set_title("(g)", loc="left", y=1.08)
            ax.set_xlabel("Time (s)")

    # get legend handles and labels
    handles, labels = axs[0, 0].get_legend_handles_labels()
    axs[0, 1].legend(
        handles=handles, labels=labels, loc="upper left", bbox_to_anchor=(0.0, 1.0)
    )
    axs[0, 1].axis("off")
    plt.tight_layout()
    plt.savefig(
        f"Figures/{save_name}_{perturbed_factor_string}_{tray}_film_simulation_1D_h_O2_RC_ROO_vs_t.pdf",
        bbox_inches="tight",
    )

print("Plot h_i, [O2_i], ... vs. z")

for perturbed_factor_string in perturbed_factor_string_list:
    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(8, 9), sharex=True)

    for num_cells in num_cells_list:
        for discretization_method in discretization_methods:
            simulation = one_d_simulations[
                perturbed_factor_string, num_cells, discretization_method
            ]
            cell_inds = range(1, num_cells + 1)
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
            print("zs: ", zs)
            zs = np.array(zs)
            zs /= zs[-1]

            label = f"{num_cells} cells ({discretization_method})"

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
            # concs_init = list(np.array([simulation.loc[0, f"OXYGEN(L)_cell_{cell_ind}"] for cell_ind in cell_inds]) / Vhatliqinfilms_init)
            # concs.insert(0, (concs[0] + concs_init[0]) / 2)
            extrapolated_conc = extrapolate_conc(concs, zs)
            concs.insert(0, extrapolated_conc)
            concs.append(concs[-1])
            print("concs: ", concs)
            ax.plot(
                zs,
                concs,
                "-o",
                label=label,
                marker=discretization_method_markers[discretization_method],
                color=num_cells_colors[num_cells],
            )
            ax.set_ylabel("[O2(L)] (mol/m^3)")
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
            # concs_init = list(np.array([np.sum([simulation.loc[0, f"{name}_cell_{cell_ind}"] for name in carbon_center_radical_names]) for cell_ind in cell_inds]) / Vhatliqinfilms_init)
            # concs.insert(0, (concs[0] + concs_init[0]) / 2)
            extrapolated_conc = extrapolate_conc(concs, zs)
            concs.insert(0, extrapolated_conc)
            concs.append(concs[-1])
            ax.plot(
                zs,
                concs,
                "-o",
                label=label,
                marker=discretization_method_markers[discretization_method],
                color=num_cells_colors[num_cells],
            )
            ax.set_ylabel("[RC.(L)] (mol/m^3)")
            if max(concs) / min(concs) > 10:
                ax.set_yscale("log")
            ax.set_title("(b)", loc="left", y=1.08)

            ax = axs[2, 0]
            concs = list(
                np.array(
                    [
                        np.sum(
                            [
                                simulation.loc[row_ind, f"{name}_cell_{cell_ind}"]
                                for name in oxygen_center_radical_names
                            ]
                        )
                        for cell_ind in cell_inds
                    ]
                )
                / Vhatliqinfilms
            )
            # concs_init = list(np.array([np.sum([simulation.loc[0, f"{name}_cell_{cell_ind}"] for name in oxygen_center_radical_names]) for cell_ind in cell_inds]) / Vhatliqinfilms_init)
            # concs.insert(0, (concs[0] + concs_init[0]) / 2)
            extrapolated_conc = extrapolate_conc(concs, zs)
            concs.insert(0, extrapolated_conc)
            concs.append(concs[-1])
            ax.plot(
                zs,
                concs,
                "-o",
                label=label,
                marker=discretization_method_markers[discretization_method],
                color=num_cells_colors[num_cells],
            )
            ax.set_ylabel("[ROO.(L)] (mol/m^3)")
            if max(concs) / min(concs) > 10:
                ax.set_yscale("log")
            ax.set_title("(c)", loc="left", y=1.08)
            ax.set_xlabel("z/h (-)")

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
                label=label,
                marker=discretization_method_markers[discretization_method],
                color=num_cells_colors[num_cells],
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
                label=label,
                marker=discretization_method_markers[discretization_method],
                color=num_cells_colors[num_cells],
            )
            ax.set_ylabel("[KR] (mol/kg)")
            if max(concs) / min(concs) > 10:
                ax.set_yscale("log")
            ax.set_title("(e)", loc="left", y=1.08)

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
                label=label,
                marker=discretization_method_markers[discretization_method],
                color=num_cells_colors[num_cells],
            )
            ax.set_ylabel("[PR] (mol/kg)")
            if max(concs) / min(concs) > 10:
                ax.set_yscale("log")
            ax.set_title("(f)", loc="left", y=1.08)
            ax.set_xlabel("z/h (-)")

    axs[2, 1].legend(bbox_to_anchor=(0.0, -0.25), loc="upper left")
    plt.tight_layout()
    plt.savefig(
        f"Figures/{save_name}_{perturbed_factor_string}_{tray}_film_simulation_1D_conc_vs_z.pdf",
        bbox_inches="tight",
    )

print("Plotting max concentrations of liquid-phase radical in film vs. z/h...")

for perturbed_factor_string in perturbed_factor_string_list:
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(8, 9), sharex=True)

    for num_cells in num_cells_list:
        for discretization_method in discretization_methods:
            simulation = one_d_simulations[
                perturbed_factor_string, num_cells, discretization_method
            ]
            cell_inds = range(1, num_cells + 1)
            row_ind = len(simulation.index) - 1

            masshats = np.array(
                [
                    simulation.loc[row_ind, f"mass_cell_{cell_ind}"]
                    for cell_ind in cell_inds
                ]
            )
            Vhatsolidinfilms = masshats / rho
            Vhatliqinfilms = np.array(
                [
                    simulation.loc[row_ind, f"Vliqinfilm_cell_{cell_ind}"]
                    for cell_ind in cell_inds
                ]
            )
            Vhatfilms = Vhatsolidinfilms / (1 - epsilon)
            t = simulation.loc[row_ind, "timestamp"]
            zs = list(np.cumsum(Vhatfilms / A) - Vhatfilms / A / 2)
            zs.insert(0, 0.0)
            zs.append(zs[-1] + Vhatfilms[-1] / A / 2)
            print("zs: ", zs)
            zs = np.array(zs)
            zs /= zs[-1]

            label = f"{num_cells} cells ({discretization_method})"

            max_radical_names = set()
            for cell_ind in cell_inds:
                concs = np.array(
                    [
                        simulation.loc[row_ind, f"{name}_cell_{cell_ind}"]
                        for name in carbon_center_radical_names
                    ]
                )
                print("concs: ", concs)
                print("max: ", np.max(concs))
                print("argmax: ", np.argmax(concs))
                max_radical_name = carbon_center_radical_names[np.argmax(concs)]
                max_radical_names.add(max_radical_name)

            ax = axs[0]
            for name in max_radical_names:
                concs = list(
                    np.array(
                        [
                            simulation.loc[row_ind, f"{name}_cell_{cell_ind}"]
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
                    label=label + f" ({name})",
                    marker=discretization_method_markers[discretization_method],
                    color=num_cells_colors[num_cells],
                )
            ax.legend(bbox_to_anchor=(1.05, 1.05), loc="upper left")
            ax.set_ylabel("[max(RC.(L))] (mol/m^3)")
            if max(concs) / min(concs) > 10:
                ax.set_yscale("log")

            # ax = axs[1]
            # max_radical_names = set()
            # for cell_ind in cell_inds:
            #     concs = np.array(
            #         [
            #             simulation.loc[row_ind, f"{name}_cell_{cell_ind}"]
            #             for name in oxygen_center_radical_names
            #         ]
            #     )
            #     max_radical_name = oxygen_center_radical_names[np.argmax(concs)]
            #     max_radical_names.add(max_radical_name)

            # for name in max_radical_names:
            #     concs = list(
            #         np.array(
            #             [
            #                 simulation.loc[row_ind, f"{name}_cell_{cell_ind}"]
            #                 for cell_ind in cell_inds
            #             ]
            #         )
            #         / Vhatliqinfilms
            #     )
            #     extrapolated_conc = extrapolate_conc(concs, zs)
            #     concs.insert(0, extrapolated_conc)
            #     concs.append(concs[-1])
            #     ax.plot(
            #         zs,
            #         concs,
            #         "-o",
            #         label=label + f" ({name})",
            #         marker=discretization_method_markers[discretization_method],
            #         color=num_cells_colors[num_cells],
            #     )
            # ax.legend(bbox_to_anchor=(1.05, 1.05), loc="upper left")
            # ax.set_ylabel("[max(ROO.(L))] (mol/m^3)")
            # if max(concs) / min(concs) > 10:
            #     ax.set_yscale("log")

            ax = axs[1]
            concs = list(
                np.array(
                    [
                        simulation.loc[row_ind, f"CDB_cell_{cell_ind}"]
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
                label=label + f" ({name})",
                marker=discretization_method_markers[discretization_method],
                color=num_cells_colors[num_cells],
            )
            if max(concs) / min(concs) > 10:
                ax.set_yscale("log")
            ax.set_ylabel("[CDB] (mol/kg)")

            ax = axs[2]
            concs = list(
                np.array(
                    [
                        simulation.loc[row_ind, f"CD_cell_{cell_ind}"]
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
                label=label + f" ({name})",
                marker=discretization_method_markers[discretization_method],
                color=num_cells_colors[num_cells],
            )
            if max(concs) / min(concs) > 10:
                ax.set_yscale("log")
            ax.set_ylabel("[CD] (mol/kg)")

    axs[-1].set_xlabel("z/h (-)")

    plt.tight_layout()
    plt.savefig(
        f"Figures/{save_name}_{perturbed_factor_string}_{tray}_film_simulation_1D_max_conc_vs_z.pdf",
        bbox_inches="tight",
    )
