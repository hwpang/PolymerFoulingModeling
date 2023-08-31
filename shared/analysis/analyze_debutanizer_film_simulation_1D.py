import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from utils import get_film_rops_1D

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

trays = range(1, 41)
# num_cells = 1
num_cells = 3
cell_inds = range(1, num_cells + 1)
method = "callback"
# method = "same_initial_size"
d = 2.5
hliq = 0.3
A = (d / 2) ** 2 * np.pi
Vliq = A * hliq
rho = 900  # solid density kg/m3
epsilon = 0.2

print("Loading asymptotic film simulation results...")

asymptotic_simulations = dict()
for tray in trays:
    asymptotic_simulation_path = os.path.join(
        simulation_directory, f"simulation_film_{tray}_asymptotic.csv"
    )
    asymptotic_simulations[tray] = pd.read_csv(asymptotic_simulation_path)

os.makedirs("Figures", exist_ok=True)

print("Plot asymptotic film simulation results")

selected_fragments = ["AR", "KR", "CDB", "AH", "CD"]

if model_name != "basecase_debutanizer_model":
    selected_fragments += ["CP", "HP", "PR", "OR"]

nrows = len(selected_fragments)
ncols = 1
figsize = (6, 12)

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, sharex=True)

cmap = plt.get_cmap("RdPu")
for fragment_ind, fragment in enumerate(selected_fragments):
    for tray_ind, tray in enumerate(trays):
        df = asymptotic_simulations[tray]
        axs[fragment_ind].plot(
            df.loc[:, "timestamp"] / 3600 / 24 / 365,
            df.loc[:, fragment] / df.loc[:, "mass"],
            label=tray,
            color=cmap(tray / len(trays)),
        )
        axs[fragment_ind].set_yscale("log")
        axs[fragment_ind].set_ylabel(f"{fragment}" + "/mass\n(mol/kg)")

axs[-1].set_xlabel("Time (year)")

fig.tight_layout()
fig.savefig(f"Figures/{model_name}_film_asymptotic.pdf", bbox_inches="tight")
plt.close()

print("Loading film simulation results...")

film_simulations = dict()
for tray in trays:
    simulation_path = os.path.join(
        simulation_directory,
        f"simulation_film_1D_{tray}_{num_cells}cells_{method}.csv",
        # f"simulation_film_1D_{tray}_callback.csv",
    )
    film_simulations[tray] = pd.read_csv(simulation_path)

traycmap = plt.get_cmap("plasma")
selected_trays = [1, 10, 20, 30, 40]


def get_film_growth_time_constants(ms, ts):
    dmdt = (ms[2:] - ms[:-2]) / (ts[2:] - ts[:-2])
    m = ms[1:-1]
    return m / dmdt / 3600 / 24 / 365


fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 4), sharey=True)

for tray in selected_trays:
    film_simulation = film_simulations[tray]

    ts = np.array(film_simulation["timestamp"]) / 3600.0 / 24.0 / 365.0
    ms = np.sum(
        [film_simulation[f"mass_cell_{cell_ind}"] for cell_ind in cell_inds], axis=0
    )
    hs = ms / rho / (1 - epsilon) / A

    ax = axs[0]
    ax.plot(
        ts,
        hs,
        color=traycmap(tray / len(trays)),
        label=f"Tray {tray}",
    )
    ax.set_yscale("log")
    ax.set_title("(a) Film thickness", loc="left")
    ax.set_ylabel("(m)")
    ax.set_xlabel("Time (yr)")

    ax.plot([0, 0], [1e-5, 1e0], "k--")
    ax.plot([0, 1 / 12 / 2], [1e-5, 1e-5], "k--")
    ax.plot([1 / 12 / 2, 1 / 12 / 2], [1e-5, 1e0], "k--")
    ax.plot([0, 1 / 12 / 2], [1e0, 1e0], "k--")

    ax = axs[1]
    ax.plot(
        ts * 12 * 4,
        hs,
        color=traycmap(tray / len(trays)),
        label=f"Tray {tray}",
    )
    ax.set_yscale("log")
    ax.set_xlim([0, 2])  # show 2 weeks of growth
    ax.set_title("Zoomed in on (a)", loc="left")
    ax.set_xlabel("Time (week)")

sm = plt.cm.ScalarMappable(cmap="plasma", norm=plt.Normalize(vmin=1, vmax=40))
cbar_ax = fig.add_axes([1.0, 0.15, 0.02, 0.7])
cbar = fig.colorbar(
    sm,
    ticks=selected_trays,
    orientation="vertical",
    label="Trays",
    cax=cbar_ax,
)

fig.tight_layout()
fig.savefig(f"Figures/{model_name}_1D_film_thickness_vs_t.pdf", bbox_inches="tight")

print("Loading rate of loss results...")

film_rate_of_productions = dict()
for tray in trays:
    for i, cell_ind in enumerate(cell_inds):
        rate_of_production_path = os.path.join(
            simulation_directory,
            f"simulation_film_1D_rop_{tray}_{num_cells}cells_cell{cell_ind}_{method}.csv",
        )
        film_rate_of_productions[tray, cell_ind] = pd.read_csv(rate_of_production_path)

print("Loading t0 rate of loss results...")

film_rate_of_productions_t0 = dict()
for tray in trays:
    for i, cell_ind in enumerate(cell_inds):
        rate_of_production_path = os.path.join(
            simulation_directory,
            f"simulation_film_1D_rop_{tray}_{num_cells}cells_cell{cell_ind}_{method}_t0.csv",
        )
        film_rate_of_productions_t0[tray, cell_ind] = pd.read_csv(
            rate_of_production_path
        )

print("Plotting rate of film growth at t0 and tf...")


def select_bar_color_pattern(comment):
    if "Diels-Alder" in comment:
        return "tab:green", "//"
    elif (
        "radical addition" in comment
        and ("PR" in comment or "OR" in comment or "[O]" in comment or "O." in comment)
    ) or "[O][O](L)" in comment:
        return "tab:red", "\\\\"
    elif "radical addition" in comment:
        return "tab:gray", "||"
    else:
        return "tab:blue", "--"


colors = {
    "Diels-Alder addition": "tab:green",
    "Oxygen-center radical pathways": "tab:red",
    "Carbon-center radical pathways": "tab:gray",
    # "Other": "tab:blue",
}
patterns = {
    "Diels-Alder addition": "//",
    "Oxygen-center radical pathways": "\\\\",
    "Carbon-center radical pathways": "||",
    # "Other": "--",
}
labels = list(colors.keys())
handles = [
    mpatches.Patch(facecolor=colors[label], hatch=patterns[label]) for label in labels
]

selected_trays = [1, 10, 20, 30, 40]
nrows = len(selected_trays) * 2
ncols = 1

fig = plt.figure(figsize=(6, 14))
gs = fig.add_gridspec(nrows, ncols)
axs = []

min_rop = 1e10
max_rop = 0
for ind, tray in enumerate(selected_trays):
    ax = fig.add_subplot(gs[ind * 2, 0])
    axs.append(ax)

    rops, rop_rxncomments, rop_rxnstrs = get_film_rops_1D(
        film_rate_of_productions_t0, tray, cell_inds, "mass"
    )
    colors_patterns = [select_bar_color_pattern(comment) for comment in rop_rxncomments]
    colors = [color for color, pattern in colors_patterns]
    patterns = [pattern for color, pattern in colors_patterns]
    df = film_simulations[tray]
    mass = np.sum(
        [df.loc[0, f"mass_cell_{cell_ind}"] for cell_ind in cell_inds], axis=0
    )
    normalized_rops = np.array(rops / mass)
    min_rop = min(min_rop, min(normalized_rops))
    max_rop = max(max_rop, max(normalized_rops))
    xs = np.arange(len(rop_rxncomments))
    for bar_ind, x in enumerate(xs):
        ax.barh(
            [x],
            [normalized_rops[bar_ind]],
            align="center",
            color=colors[bar_ind],
            hatch=patterns[bar_ind],
        )
    ax.set_yticks(xs)
    ax.set_yticklabels(rop_rxncomments)
    ax.set_xscale("log")
    ax.invert_yaxis()
    ax.set_ylabel("($t_0$)")

    ax = fig.add_subplot(gs[ind * 2 + 1, 0])
    axs.append(ax)

    rops, rop_rxncomments, rop_rxnstrs = get_film_rops_1D(
        film_rate_of_productions, tray, cell_inds, "mass"
    )
    colors_patterns = [select_bar_color_pattern(comment) for comment in rop_rxncomments]
    colors = [color for color, pattern in colors_patterns]
    patterns = [pattern for color, pattern in colors_patterns]
    df = film_simulations[tray]
    mass = np.sum(
        [df.loc[len(df.index) - 1, f"mass_cell_{cell_ind}"] for cell_ind in cell_inds],
        axis=0,
    )
    normalized_rops = np.array(rops / mass)
    min_rop = min(min_rop, min(normalized_rops))
    max_rop = max(max_rop, max(normalized_rops))
    xs = np.arange(len(rop_rxncomments))
    for bar_ind, x in enumerate(xs):
        ax.barh(
            [x],
            [normalized_rops[bar_ind]],
            align="center",
            color=colors[bar_ind],
            hatch=patterns[bar_ind],
        )
    ax.set_yticks(xs)
    ax.set_yticklabels(rop_rxncomments)
    ax.set_xscale("log")
    ax.invert_yaxis()
    ax.set_ylabel("($t_f$)")

    ax = fig.add_subplot(gs[(ind * 2) : (ind * 2 + 2), :], frameon=False)
    ax.set_ylabel(f"Tray {tray}", labelpad=20)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_yticklabels([])

print("min_rop: ", min_rop)
print("max_rop: ", max_rop)

for ax in axs:
    ax.set_xlim(min_rop, max_rop)
    if ax != axs[-1]:
        ax.set_xticks([])
        ax.set_xticklabels([])

axs[-1].set_xlabel("Rate of film growth (kg/(kg*s))")
axs[-1].legend(
    handles=handles, labels=labels, bbox_to_anchor=(1.00, -0.5), loc="upper right"
)
fig.align_labels()
fig.tight_layout()
fig.savefig(f"Figures/{model_name}_1D_film_rop_mass_t0_tf.pdf", bbox_inches="tight")

# print("Plotting rate of film growth at tf near bulk liquid and tray surface...")

# fig = plt.figure(figsize=(6, 14))
# gs = fig.add_gridspec(nrows, ncols)
# axs = []

# min_rop = 1e10
# max_rop = 0
# for ind, tray in enumerate(selected_trays):
#     ax = fig.add_subplot(gs[ind * 2, 0])
#     axs.append(ax)

#     rops, rop_rxncomments, rop_rxnstrs = get_film_rops_1D(
#         film_rate_of_productions, tray, [cell_inds[0]], "mass"
#     )
#     colors_patterns = [select_bar_color_pattern(comment) for comment in rop_rxncomments]
#     colors = [color for color, pattern in colors_patterns]
#     patterns = [pattern for color, pattern in colors_patterns]
#     df = film_simulations[tray]
#     mass = np.sum(
#         [
#             df.loc[len(df.index) - 1, f"mass_cell_{cell_ind}"]
#             for cell_ind in [cell_inds[0]]
#         ],
#         axis=0,
#     )
#     normalized_rops = np.array(rops / mass)
#     min_rop = min(min_rop, min(normalized_rops))
#     max_rop = max(max_rop, max(normalized_rops))
#     xs = np.arange(len(rop_rxncomments))
#     for bar_ind, x in enumerate(xs):
#         ax.barh(
#             [x],
#             [normalized_rops[bar_ind]],
#             align="center",
#             color=colors[bar_ind],
#             hatch=patterns[bar_ind],
#         )
#     ax.set_yticks(xs)
#     ax.set_yticklabels(rop_rxncomments)
#     ax.set_xscale("log")
#     ax.invert_yaxis()

#     ax = fig.add_subplot(gs[ind * 2 + 1, 0])
#     axs.append(ax)

#     cell_ind = cell_inds[-1]
#     rops, rop_rxncomments, rop_rxnstrs = get_film_rops_1D(
#         film_rate_of_productions, tray, [cell_inds[0]], "mass"
#     )
#     colors_patterns = [select_bar_color_pattern(comment) for comment in rop_rxncomments]
#     colors = [color for color, pattern in colors_patterns]
#     patterns = [pattern for color, pattern in colors_patterns]
#     df = film_simulations[tray]
#     mass = np.sum(
#         [
#             df.loc[len(df.index) - 1, f"mass_cell_{cell_ind}"]
#             for cell_ind in [cell_inds[0]]
#         ],
#         axis=0,
#     )
#     normalized_rops = np.array(rops / mass)
#     min_rop = min(min_rop, min(normalized_rops))
#     max_rop = max(max_rop, max(normalized_rops))
#     xs = np.arange(len(rop_rxncomments))
#     for bar_ind, x in enumerate(xs):
#         ax.barh(
#             [x],
#             [normalized_rops[bar_ind]],
#             align="center",
#             color=colors[bar_ind],
#             hatch=patterns[bar_ind],
#         )
#     ax.set_yticks(xs)
#     ax.set_yticklabels(rop_rxncomments)
#     ax.set_xscale("log")
#     ax.invert_yaxis()

#     ax = fig.add_subplot(gs[(ind * 2) : (ind * 2 + 2), :], frameon=False)
#     ax.set_ylabel(f"Tray {tray}\n(tray surface)      (bulk liquid)")
#     ax.set_xticks([])
#     ax.set_xticklabels([])
#     ax.set_yticks([])
#     ax.set_yticklabels([])

# print("min_rop: ", min_rop)
# print("max_rop: ", max_rop)

# for ax in axs:
#     ax.set_xlim(min_rop, max_rop)
#     if ax != axs[-1]:
#         ax.set_xticks([])
#         ax.set_xticklabels([])

# axs[-1].set_xlabel("Local rate of film growth (kg/(kg*s))")
# axs[-1].legend(
#     handles=handles, labels=labels, bbox_to_anchor=(1.00, -0.5), loc="upper right"
# )
# fig.align_labels()
# fig.tight_layout()
# fig.savefig(
#     f"Figures/{model_name}_1D_film_rop_mass_bulk_liquid_tray_surface.pdf",
#     bbox_inches="tight",
# )


def plot_fragment_rops_1D(species_label):
    fig = plt.figure(figsize=(6, 14))
    gs = fig.add_gridspec(nrows, ncols)
    axs = []

    min_rop = 1e10
    max_rop = 0
    for ind, tray in enumerate(selected_trays):
        ax = fig.add_subplot(gs[ind * 2, 0])
        axs.append(ax)

        rops, rop_rxncomments, rop_rxnstrs = get_film_rops_1D(
            film_rate_of_productions_t0, tray, cell_inds, species_label, loss_only=True
        )
        colors_patterns = [
            select_bar_color_pattern(comment) for comment in rop_rxncomments
        ]
        colors = [color for color, pattern in colors_patterns]
        patterns = [pattern for color, pattern in colors_patterns]
        df = film_simulations[tray]
        mass = np.sum(
            [df.loc[0, f"mass_cell_{cell_ind}"] for cell_ind in cell_inds], axis=0
        )
        normalized_rops = np.array(np.abs(rops) / mass)
        min_rop = min(min_rop, min(normalized_rops))
        max_rop = max(max_rop, max(normalized_rops))
        xs = np.arange(len(rop_rxncomments))
        for bar_ind, x in enumerate(xs):
            ax.barh(
                [x],
                [normalized_rops[bar_ind]],
                align="center",
                color=colors[bar_ind],
                hatch=patterns[bar_ind],
            )
        ax.set_yticks(xs)
        ax.set_yticklabels(rop_rxncomments)
        ax.set_xscale("log")
        ax.invert_yaxis()
        ax.set_ylabel("($t_0$)")

        ax = fig.add_subplot(gs[ind * 2 + 1, 0])
        axs.append(ax)

        rops, rop_rxncomments, rop_rxnstrs = get_film_rops_1D(
            film_rate_of_productions, tray, cell_inds, species_label, loss_only=True
        )
        colors_patterns = [
            select_bar_color_pattern(comment) for comment in rop_rxncomments
        ]
        colors = [color for color, pattern in colors_patterns]
        patterns = [pattern for color, pattern in colors_patterns]
        df = film_simulations[tray]
        mass = np.sum(
            [
                df.loc[len(df.index) - 1, f"mass_cell_{cell_ind}"]
                for cell_ind in cell_inds
            ],
            axis=0,
        )
        normalized_rops = np.array(np.abs(rops) / mass)
        min_rop = min(min_rop, min(normalized_rops))
        max_rop = max(max_rop, max(normalized_rops))
        xs = np.arange(len(rop_rxncomments))
        for bar_ind, x in enumerate(xs):
            ax.barh(
                [x],
                [normalized_rops[bar_ind]],
                align="center",
                color=colors[bar_ind],
                hatch=patterns[bar_ind],
            )
        ax.set_yticks(xs)
        ax.set_yticklabels(rop_rxncomments)
        ax.set_xscale("log")
        ax.invert_yaxis()
        ax.set_ylabel("($t_f$)")

        ax = fig.add_subplot(gs[(ind * 2) : (ind * 2 + 2), :], frameon=False)
        ax.set_ylabel(f"Tray {tray}", labelpad=20)
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_yticks([])
        ax.set_yticklabels([])

    print("min_rop: ", min_rop)
    print("max_rop: ", max_rop)

    for ax in axs:
        ax.set_xlim(min_rop, max_rop)
        if ax != axs[-1]:
            ax.set_xticks([])
            ax.set_xticklabels([])

    axs[-1].set_xlabel(f"Rate of {species_label} loss (mol/(kg*s))")
    axs[-1].legend(
        handles=handles, labels=labels, bbox_to_anchor=(1.00, -0.5), loc="upper right"
    )
    fig.align_labels()
    fig.tight_layout()
    fig.savefig(
        f"Figures/{model_name}_1D_film_rop_{species_label}_t0_tf.pdf",
        bbox_inches="tight",
    )


print("Plotting rate of loss for AR...")

species_label = "AR"
plot_fragment_rops_1D(species_label)

print("Plotting rate of loss for KR...")

species_label = "KR"
plot_fragment_rops_1D(species_label)

simulation = film_simulations[trays[0]]
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

# print("Plotting rate of loss for carbon center radicals...")

# fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9, 12), sharex=True)

# min_rop = 1e10
# max_rop = 0
# for ind, tray in enumerate(selected_trays):
#     rops, rop_rxncomments, rop_rxnstrs = get_film_rops_1D(
#         film_rate_of_productions,
#         tray,
#         cell_inds,
#         carbon_center_radical_names,
#         loss_only=True,
#     )
#     rop_rxnlabels = [
#         comment if isinstance(comment, str) else rxnstr
#         for comment, rxnstr in zip(rop_rxncomments, rop_rxnstrs)
#     ]
#     df = film_simulations[tray]
#     Vliqinfilm = df.loc[len(df.index) - 1, f"Vliqinfilm_cell_{cell_ind}"]
#     normalized_rops = np.abs(rops) / Vliqinfilm
#     min_rop = min(min_rop, min(normalized_rops))
#     max_rop = max(max_rop, max(normalized_rops))
#     xs = np.arange(len(rop_rxncomments))
#     axs[ind].barh(xs, normalized_rops, align="center")
#     axs[ind].set_yticks(xs)
#     axs[ind].set_yticklabels(rop_rxnlabels)
#     axs[ind].set_ylabel(f"Tray {tray}")
#     axs[ind].set_xscale("log")
#     axs[ind].invert_yaxis()

# print("min_rop: ", min_rop)
# print("max_rop: ", max_rop)

# for ind, tray in enumerate(selected_trays):
#     axs[ind].set_xlim(min_rop, max_rop)

# axs[-1].set_xlabel("Rate of loss for RC.(L in film) (mol/(m^3*s))")

# fig.tight_layout()
# fig.savefig(f"Figures/{model_name}_1D_film_rop_RC..pdf", bbox_inches="tight")

# print("Plotting rate of loss for oxygen center radicals...")

# fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9, 12), sharex=True)

# min_rop = 1e10
# max_rop = 0
# for ind, tray in enumerate(selected_trays):
#     rops, rop_rxncomments, rop_rxnstrs = get_film_rops_1D(
#         film_rate_of_productions,
#         tray,
#         cell_inds,
#         oxygen_center_radical_names,
#         loss_only=True,
#     )
#     rop_rxnlabels = [
#         comment if isinstance(comment, str) else rxnstr
#         for comment, rxnstr in zip(rop_rxncomments, rop_rxnstrs)
#     ]
#     df = film_simulations[tray]
#     Vliqinfilm = df.loc[len(df.index) - 1, f"Vliqinfilm_cell_{cell_ind}"]
#     normalized_rops = np.abs(rops) / Vliqinfilm
#     min_rop = min(min_rop, min(normalized_rops))
#     max_rop = max(max_rop, max(normalized_rops))
#     xs = np.arange(len(rop_rxncomments))
#     axs[ind].barh(xs, normalized_rops, align="center")
#     axs[ind].set_yticks(xs)
#     axs[ind].set_yticklabels(rop_rxnlabels)
#     axs[ind].set_ylabel(f"Tray {tray}")
#     axs[ind].set_xscale("log")
#     axs[ind].invert_yaxis()

# print("min_rop: ", min_rop)
# print("max_rop: ", max_rop)

# for ind, tray in enumerate(selected_trays):
#     axs[ind].set_xlim(min_rop, max_rop)

# axs[-1].set_xlabel("Rate of loss for ROO.(L in film) (mol/(m^3*s))")

# fig.tight_layout()
# fig.savefig(f"Figures/{model_name}_1D_film_rop_ROO..pdf", bbox_inches="tight")

if model_name != "basecase_debutanizer_model":
    print("Plotting rate of loss for PR...")

    species_label = "PR"
    plot_fragment_rops_1D(species_label)

    print("Plotting rate of loss for OR...")

    species_label = "OR"
    plot_fragment_rops_1D(species_label)
