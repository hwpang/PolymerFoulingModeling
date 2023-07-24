import os
import yaml
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils import get_film_rops

#change default font size to 12
plt.rcParams.update({'font.size': 12})

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model_name", type=str, required=True, help="The name of the model.",
    )
    parser.add_argument(
        "--simulation_directory", type=str, required=True, help="The path to the film simulation results.",
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

print("Loading asymptotic film simulation results...")

asymptotic_simulations = dict()
for tray in trays:
    asymptotic_simulation_path = os.path.join(simulation_directory, f"simulation_film_{tray}_asymptotic.csv")
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
        axs[fragment_ind].plot(df.loc[:, "timestamp"] / 3600 / 24 / 365, df.loc[:, fragment] / df.loc[:, "mass"], label=tray, color=cmap(tray / len(trays)))
        axs[fragment_ind].set_yscale("log")
        axs[fragment_ind].set_ylabel(f"{fragment}"+"/mass\n(mol/kg)")

axs[-1].set_xlabel("Time (year)")

fig.tight_layout()
fig.savefig(f"Figures/{model_name}_film_asymptotic.pdf", bbox_inches="tight")
plt.close()

print("Loading rate of production restuls...")

rate_of_productions = dict()
for tray in trays:
    rate_of_production_path = os.path.join(simulation_directory, f"simulation_film_rop_{tray}.csv")
    rate_of_productions[tray] = pd.read_csv(rate_of_production_path)

print("Loading film simulation results...")

simulations = dict()
for tray in trays:
    simulation_path = os.path.join(simulation_directory, f"simulation_film_{tray}.csv")
    simulations[tray] = pd.read_csv(simulation_path)

print("Plotting rate of film growth...")

selected_trays = [1, 5, 10, 15, 20, 25, 30, 35, 40]
nrows = len(selected_trays)
ncols = 1

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9, 12), sharex=True)

min_rop = 1e10
max_rop = 0
for ind, tray in enumerate(selected_trays):
    rops, rop_rxncomments, rop_rxnstrs = get_film_rops(rate_of_productions[tray], "mass")
    df = simulations[tray]
    mass = df.loc[len(df.index)-1, "mass"]
    normalized_rops = rops / mass
    min_rop = min(min_rop, min(normalized_rops))
    max_rop = max(max_rop, max(normalized_rops))
    xs = np.arange(len(rop_rxncomments))
    axs[ind].barh(xs, normalized_rops, align="center")
    axs[ind].set_yticks(xs)
    axs[ind].set_yticklabels(rop_rxncomments)
    axs[ind].set_ylabel(f"Tray {tray}")
    axs[ind].set_xscale("log")
    axs[ind].invert_yaxis()

print("min_rop: ", min_rop)
print("max_rop: ", max_rop)

for ax in axs:
    ax.set_xlim(min_rop, max_rop)

axs[-1].set_xlabel("Rate of film growth (kg/(kg*s))")
fig.tight_layout()
fig.savefig(f"Figures/{model_name}_film_rop_mass.pdf", bbox_inches="tight")

print("Plot rate of production for AR")

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9, 12), sharex=True)

min_rop = 1e10
max_rop = 0
for ind, tray in enumerate(selected_trays):
    rops, rop_rxncomments, rop_rxnstrs = get_film_rops(rate_of_productions[tray], "AR", loss_only=True)
    df = simulations[tray]
    mass = df.loc[len(df.index)-1, "mass"]
    normalized_rops = np.abs(rops) / mass
    min_rop = min(min_rop, min(normalized_rops))
    max_rop = max(max_rop, max(normalized_rops))
    xs = np.arange(len(rop_rxncomments))
    axs[ind].barh(xs, normalized_rops, align="center")
    axs[ind].set_yticks(xs)
    axs[ind].set_yticklabels(rop_rxncomments)
    axs[ind].set_ylabel(f"Tray {tray}")
    axs[ind].set_xscale("log")
    axs[ind].invert_yaxis()

for ax in axs:
    ax.set_xlim(min_rop, max_rop)

axs[-1].set_xlabel("Rate of loss of AR (mol/(kg*s))")
fig.tight_layout()
fig.savefig(f"Figures/{model_name}_film_rop_AR.pdf", bbox_inches="tight")
plt.close()

print("Plot rate of production for KR")

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9, 12), sharex=True)

min_rop = 1e10
max_rop = 0
for ind, tray in enumerate(selected_trays):
    rops, rop_rxncomments, rop_rxnstrs = get_film_rops(rate_of_productions[tray], "KR", loss_only=True)
    df = simulations[tray]
    mass = df.loc[len(df.index)-1, "mass"]
    normalized_rops = np.abs(rops) / mass
    min_rop = min(min_rop, min(normalized_rops))
    max_rop = max(max_rop, max(normalized_rops))
    xs = np.arange(len(rop_rxncomments))
    axs[ind].barh(xs, normalized_rops, align="center")
    axs[ind].set_yticks(xs)
    axs[ind].set_yticklabels(rop_rxncomments)
    axs[ind].set_ylabel(f"Tray {tray}")
    axs[ind].set_xscale("log")
    axs[ind].invert_yaxis()

for ax in axs:
    ax.set_xlim(min_rop, max_rop)

axs[-1].set_xlabel("Rate of loss of KR (mol/(kg*s))")
fig.tight_layout()
fig.savefig(f"Figures/{model_name}_film_rop_KR.pdf", bbox_inches="tight")
plt.close()

if model_name != "basecase_debutanizer_model":

    print("Plot rate of production for PR")

    min_rop = 1e10
    max_rop = 0
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9, 12), sharex=True)

    for ind, tray in enumerate(selected_trays):
        rops, rop_rxncomments, rop_rxnstrs = get_film_rops(rate_of_productions[tray], "PR", loss_only=True)
        df = simulations[tray]
        mass = df.loc[len(df.index)-1, "mass"]
        normalized_rops = np.abs(rops) / mass
        min_rop = min(min_rop, min(normalized_rops))
        max_rop = max(max_rop, max(normalized_rops))
        xs = np.arange(len(rop_rxncomments))
        axs[ind].barh(xs, normalized_rops, align="center")
        axs[ind].set_yticks(xs)
        axs[ind].set_yticklabels(rop_rxncomments)
        axs[ind].set_ylabel(f"Tray {tray}")
        axs[ind].set_xscale("log")
        axs[ind].invert_yaxis()

    for ax in axs:
        ax.set_xlim(min_rop, max_rop)

    axs[-1].set_xlabel("Rate of loss of PR (mol/(kg*s))")
    fig.tight_layout()
    fig.savefig(f"Figures/{model_name}_film_rop_PR.pdf", bbox_inches="tight")
    plt.close()

    print("Plot rate of production for OR")

    min_rop = 1e10
    max_rop = 0
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9, 12), sharex=True)

    for ind, tray in enumerate(selected_trays):
        rops, rop_rxncomments, rop_rxnstrs = get_film_rops(rate_of_productions[tray], "OR", loss_only=True)
        df = simulations[tray]
        mass = df.loc[len(df.index)-1, "mass"]
        normalized_rops = np.abs(rops) / mass
        min_rop = min(min_rop, min(normalized_rops))
        max_rop = max(max_rop, max(normalized_rops))
        xs = np.arange(len(rop_rxncomments))
        axs[ind].barh(xs, normalized_rops, align="center")
        axs[ind].set_yticks(xs)
        axs[ind].set_yticklabels(rop_rxncomments)
        axs[ind].set_ylabel(f"Tray {tray}")
        axs[ind].set_xscale("log")
        axs[ind].invert_yaxis()

    for ax in axs:
        ax.set_xlim(min_rop, max_rop)

    axs[-1].set_xlabel("Rate of loss of OR (mol/(kg*s))")
    fig.tight_layout()
    fig.savefig(f"Figures/{model_name}_film_rop_OR.pdf", bbox_inches="tight")
    plt.close()
