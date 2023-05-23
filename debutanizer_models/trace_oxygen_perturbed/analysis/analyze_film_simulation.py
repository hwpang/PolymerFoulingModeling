import os
import yaml
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
selected_trays = [1, 10, 20, 30, 40]

print("Load asymptotic film simulation results")

asymptotic_simulations = dict()
for tray in trays:
    asymptotic_simulation_path = os.path.join(simulation_directory, f"simulation_film_{tray}_asymptotic.csv")
    asymptotic_simulations[tray] = pd.read_csv(asymptotic_simulation_path)

print("Plot asymptotic film simulation results")

selected_fragments = ["AR", "KR", "CDB", "AH", "CP", "HP", "PR", "OR"]

nrows = len(selected_fragments)
ncols = 1
figsize = (6, 12)

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, sharex=True)

cmap = plt.get_cmap("RdPu")
for fragment_ind, fragment in enumerate(selected_fragments):
    for tray_ind, tray in enumerate(selected_trays):
        df = asymptotic_simulations[tray]
        axs[fragment_ind].plot(df.loc[:, "timestamp"] / 3600 / 24 / 365, df.loc[:, fragment] / df.loc[:, "mass"], label=tray, color=cmap(tray / len(selected_trays)))
        axs[fragment_ind].set_yscale("log")
        axs[fragment_ind].set_ylabel(f"{fragment}"+"/mass\n(mol/kg)")
        axs[fragment_ind].set_xlim([0, 100])

axs[-1].set_xlabel("Time (year)")

fig.tight_layout()
fig.savefig(f"{model_name}_film_asymptotic.pdf", bbox_inches="tight")
plt.close()

print("Load rate of production restuls")

rate_of_productions = dict()
for tray in trays:
    rate_of_production_path = os.path.join(simulation_directory, f"simulation_film_rop_{tray}.csv")
    rate_of_productions[tray] = pd.read_csv(rate_of_production_path)

print("Plot rate of film growth")

def get_rops(df, rop_name, loss_only=False, production_only=False, N=5):
    name_inds = df["rop_spcname"] == rop_name
    rop_rxncomments = df.loc[name_inds, "rop_rxncomment"]
    # rop_rxncomments = [rxnstr.replace("\n", " ") for rxnstr in rop_rxncomments]
    # rop_rxncomments = [rxnstr.replace(" H abstraction", "") for rxnstr in rop_rxncomments]
    rop_rxnstrs = df.loc[name_inds, "rop_rxnstr"]
    rops = df.loc[name_inds, "rop"]
    
    if loss_only:
        loss_inds = rops < 0
        rops = rops[loss_inds]
        rop_rxncomments = rop_rxncomments[loss_inds]
        rop_rxnstrs = rop_rxnstrs[loss_inds]
    elif production_only:
        prod_inds = rops > 0
        rops = rops[prod_inds]
        rop_rxncomments = rop_rxncomments[prod_inds]
        rop_rxnstrs = rop_rxnstrs[prod_inds]
    
    sorted_inds = sorted(range(len(rops)), key=lambda i: abs(rops[i]), reverse=True)
    if len(sorted_inds) > N:
        sorted_inds = sorted_inds[:N]
    
    return rops[sorted_inds], rop_rxncomments[sorted_inds], rop_rxnstrs[sorted_inds]

nrows = len(selected_trays)
ncols = 1

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9, 12), sharex=True)

for ind, tray in enumerate(selected_trays):
    rops, rop_rxncomments, rop_rxnstrs = get_rops(rate_of_productions[tray], "mass")
    xs = np.arange(len(rop_rxncomments))
    mass = film_simulations[perturb_species, perturb_factor][tray][-1, "mass"]
    axs[ind].barh(xs, rops / mass, align="center")
    axs[ind].set_yticks(xs)
    axs[ind].set_yticklabels(rop_rxncomments)
    axs[ind].set_ylabel(f"Tray {tray}")
    axs[ind].set_xscale("log")
    axs[ind].invert_yaxis()
    axs[ind].set_xlim([1e-16, 1e-9])

axs[-1, 0].set_xlabel("Rate of film growth (kg/(kg*s))")
fig.tight_layout()
fig.savefig(f"{model_name}_film_rop_mass.pdf", bbox_inches="tight")
