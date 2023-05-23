import os
import yaml
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils import get_rops

#change default font size to 12
plt.rcParams.update({'font.size': 12})

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model_name", type=str, required=True, help="The name of the model.",
    )
    parser.add_argument(
        "--all_simulation_directory", type=str, required=True, help="The path to all simulation results.",
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
h = 0.3
A = (d / 2) ** 2 * np.pi
Vliq = A * h

trays = range(1, 41)
perturbed_species_list = ["OXYGEN"]
perturbed_factor_list = ["0.0", "1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1", "1e0"]
perturbed_factor_num_list = [float(perturbed_factor) for perturbed_factor in perturbed_factor_list]
delta_t = 32.0
tf_liq = np.arange(0, 3600.0 + delta_t, delta_t)[-1]

print("Load all liquid simulation results")
liquid_simulations = dict()
for perturbed_species in perturbed_species_list:
    for perturbed_factor in perturbed_factor_list:
        simulation_path = os.path.join(all_simulation_directory, f"{perturbed_species}_{perturbed_factor}_3600.0_{delta_t}", f"simulation_vapor_liquid_yliqn_{tf_liq}.csv")
        liquid_simulations[(perturbed_species, perturbed_factor)] = pd.read_csv(simulation_path)

print("Load all film simulation results")

film_simulations = dict()
for perturbed_species in perturbed_species_list:
    for perturbed_factor in perturbed_factor_list:
        for tray in trays:
            simulation_path = os.path.join(all_simulation_directory, f"{perturbed_species}_{perturbed_factor}_3600.0_{delta_t}", f"simulation_film_{tray}.csv")
            film_simulations[(perturbed_species, perturbed_factor, tray)] = pd.read_csv(simulation_path)

print("Load radical labels")

alpha_rates_path = os.path.join(all_simulation_directory, f"{perturbed_species_list[0]}_{perturbed_factor_list[0]}_3600.0_{delta_t}", "alpha_rates.yml")
with open(alpha_rates_path, "r") as f:
    results = yaml.load(f, Loader=yaml.FullLoader)
    all_alphas, consumption_rates, production_rates, rxn_rates, radical_labels = results
    alpha1, alpha2, alphas, alphas_DA = all_alphas
    R_labels, ROO_labels, RO_labels = radical_labels

print("Plot all sensitivity simulation results")

def calculate_film_growth_time_constant(df):
    m = df.loc[-2, "mass"]
    dmdt = (df.loc[-1, "mass"] - df.loc[-3, "mass"])/(df.loc[-1, "timestamp"] - df.loc[-3, "timestamp"])
    return m / dmdt / 3600.0 / 24.0 / 365.0 # year

def calculate_fragment_per_mass(df, label):
    return df.loc[-1, label] / df.loc[-1, "mass"] 

factorcmap = plt.get_cmap("PuRd")

nrows = 4
ncols = 2

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6, 8), sharex=True, sharey="row")
for species_ind, perturbed_species in enumerate(perturbed_species_list):
    for factor_ind, perturbed_factor in enumerate(perturbed_factor_list):
        ax = axs[1, 1]
        liquid_carbon_center_radical_concs = liquid_simulations[perturbed_species, perturbed_factor].loc[:, R_labels].sum(axis=1) / Vliq
        ax.scatter(trays, liquid_carbon_center_radical_concs, color=factorcmap(factor_ind / len(perturbed_factor_list)))
        ax.set_ylabel("RC.(L) (mol/m³)")
        ax.set_yscale("log")
        ax.set_title("(a)")

        ax = axs[1, 2]
        liquid_peroxyl_radical_concs = liquid_simulations[perturbed_species, perturbed_factor].loc[:, ROO_labels].sum(axis=1) / Vliq
        ax.scatter(trays, liquid_peroxyl_radical_concs, color=factorcmap(factor_ind / len(perturbed_factor_list)))
        ax.set_ylabel("ROO.(L) (mol/m³)")
        ax.set_yscale("log")
        ax.set_title("(b)")

        ax = axs[2, 1]
        liquid_alkoxyl_radical_concs = liquid_simulations[perturbed_species, perturbed_factor].loc[:, RO_labels].sum(axis=1) / Vliq
        ax.scatter(trays, liquid_alkoxyl_radical_concs, color=factorcmap(factor_ind / len(perturbed_factor_list)))
        ax.set_ylabel("RO.(L) (mol/m³)")
        ax.set_yscale("log")
        ax.set_title("(c)")

        ax = axs[2, 2]
        film_growth_time_constants = [calculate_film_growth_time_constant(film_simulations[perturbed_species, perturbed_factor, tray]) for tray in trays]
        ax.scatter(trays, film_growth_time_constants, color=factorcmap(factor_ind / len(perturbed_factor_list)))
        ax.set_ylabel("Film growth τ (yr)")
        ax.set_yscale("log")
        ax.set_title("(d)")

        ax = axs[3, 1]
        AR_concs = [calculate_fragment_per_mass(film_simulations[perturbed_species, perturbed_factor][tray], "AR") for tray in trays]
        ax.scatter(trays, AR_concs, color=factorcmap(factor_ind / len(perturbed_factor_list)))
        ax.set_ylabel("AR/mass (mol/kg)")
        ax.set_yscale("log")
        ax.set_title("(e)")

        ax = axs[3, 2]
        KR_concs = [calculate_fragment_per_mass(film_simulations[perturbed_species, perturbed_factor][tray], "KR") for tray in trays]
        ax.scatter(trays, KR_concs, color=factorcmap(factor_ind / len(perturbed_factor_list)))
        ax.set_ylabel("KR/mass (mol/kg)")
        ax.set_yscale("log")
        ax.set_title("(f)")

        ax = axs[4, 1]
        PR_concs = [calculate_fragment_per_mass(film_simulations[perturbed_species, perturbed_factor][tray], "PR") for tray in trays]
        ax.scatter(trays, PR_concs, color=factorcmap(factor_ind / len(perturbed_factor_list)))
        ax.set_ylabel("PR/mass (mol/kg)")
        ax.set_yscale("log")
        ax.set_title("(g)")

        ax = axs[4, 2]
        OR_concs = [calculate_fragment_per_mass(film_simulations[perturbed_species, perturbed_factor][tray], "OR") for tray in trays]
        ax.scatter(trays, OR_concs, color=factorcmap(factor_ind / len(perturbed_factor_list)))
        ax.set_ylabel("OR/mass (mol/kg)")
        ax.set_yscale("log")
        ax.set_title("(h)")

axs[4, 1].set_xlabel("Trays")
axs[4, 2].set_xlabel("Trays")

cbar_ax = fig.add_axes([1.0, 0.15, 0.02, 0.7])
sm = plt.cm.ScalarMappable(cmap="RdPu", norm=plt.Normalize(vmin=0.5, vmax=1.9))
cbar = fig.colorbar(sm, ticks=perturbed_factor_num_list, orientation="vertical", label="Perturbation", cax=cbar_ax)

fig.tight_layout()
fig.savefig(f"{model_name}_sens_all.pdf", bbox_inches="tight")
plt.close()