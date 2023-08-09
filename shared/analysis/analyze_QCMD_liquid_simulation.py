import os
import yaml
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils import get_last_row, get_liquid_rops

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
print("alls_imulation_directory: ", all_simulation_directory)

Vreactor = 40 * 1e-9
d = 0.55  # in
d = d / 2.54 / 100  # m
A = (d / 2) ** 2 * np.pi
h = Vreactor / A
Vliq = Vreactor
hsolid0 = 100 * 1e-9
Vsolidinfilm0 = A * hsolid0
epsilon = 0.2  # vol% of liquid in swollen film
rho = 900.0  # density of solid
rho_liq = 970.2  # density of liquid 30vol% MCHD in Dowtherm A
tf0 = 10000.0
T = 90.0 + 273.15
tray = 1

perturb_species = "O2"
perturb_factor_list = ["0.0", "1e-3", "1e-2", "1e-1", "1e0"]
factor_num_list = [float(perturb_factor) for perturb_factor in perturb_factor_list]

print("Loading liquid simulation results...")
liquid_simulations = dict()
for perturb_factor in perturb_factor_list:
    liquid_simulation_path = os.path.join(
        all_simulation_directory,
        f"{perturb_species}_{perturb_factor}",
        f"simulation_liquid_{tray}.csv",
    )
    liquid_simulations[perturb_factor] = pd.read_csv(liquid_simulation_path)

print("Loading liquid rops...")
liquid_rop_results = dict()
for perturb_factor in perturb_factor_list:
    liquid_rop_path = os.path.join(
        all_simulation_directory,
        f"{perturb_species}_{perturb_factor}",
        f"simulation_liquid_liqrop_{tray}.csv",
    )
    liquid_rop_results[perturb_factor] = pd.read_csv(liquid_rop_path)

print("Loading alpha rates...")
with open(
    os.path.join(
        all_simulation_directory,
        f"{perturb_species}_{perturb_factor_list[0]}",
        f"alpha_rates.yml",
    ),
    "r",
) as f:
    results = yaml.load(f, Loader=yaml.FullLoader)
    all_alphas, consumption_rates, production_rates, rxn_rates, radical_labels = results
    alpha1, alpha2, alphas, alphas_DA = all_alphas
    R_labels, ROO_labels, RO_labels = radical_labels
    alpha1 = np.array(alpha1)
    alpha2 = np.array(alpha2)
    alphas = np.array(alphas)
    alphas_DA = np.array(alphas_DA)
    for key in consumption_rates.keys():
        consumption_rates[key] = np.array(consumption_rates[key])

# pure_O2_sat_conc = 0.008093871600706785 #M estimated by RMG
pure_O2_sat_molfrac = (
    8.1 * 1e-4
)  # mol fraction measured experimentally from Battino, R., Rettich, T. R., & Tominaga, T. (1983). The Solubility of Oxygen and Ozone in Liquids. Journal of Physical and Chemical Reference Data, Vol. 12, pp. 163–178. https://doi.org/10.1063/1.555680
benzene_density = 876  # kg/m^3 from Lide, D. R., Data, S. R., Board, E. A., Baysinger, G., Chemistry, S., Library, C. E., … Zwillinger, D. (2004). CRC Handbook of Chemistry and Physics. 2660.
benzene_mw = (
    78.11 / 1000
)  # kg/mol from Lide, D. R., Data, S. R., Board, E. A., Baysinger, G., Chemistry, S., Library, C. E., … Zwillinger, D. (2004). CRC Handbook of Chemistry and Physics. 2660.
benzene_sat_conc = benzene_density / benzene_mw  # mol/m^3
pure_O2_sat_conc = pure_O2_sat_molfrac * benzene_sat_conc  # mol/m^3
air_O2_sat_conc = 0.21 * pure_O2_sat_conc  # mol/m^3
air_O2_sat_conc /= 1000  # mol/L

xs = np.array(factor_num_list) * air_O2_sat_conc

os.makedirs("Figures", exist_ok=True)

print("Plotting max radical concentrations vs. [O2]...")


def get_max_conc_liq_radical_label(df, labels):
    radical_mols = get_last_row(df, labels)
    return labels[np.argmax(radical_mols)]


max_conc_liq_radical_labels = set(
    [
        get_max_conc_liq_radical_label(
            liquid_simulations[perturb_factor], R_labels + ROO_labels + RO_labels
        )
        for perturb_factor in perturb_factor_list
    ]
)
print("max_conc_liq_radical_labels: ", max_conc_liq_radical_labels)

fig, axs = plt.subplots(1, 1, figsize=(4, 3))
ax = axs
for label in max_conc_liq_radical_labels:
    ax.plot(
        xs,
        [
            get_last_row(liquid_simulations[perturb_factor], label) / Vliq
            for perturb_factor in perturb_factor_list
        ],
        marker="o",
        label=label,
    )
ax.set_yscale("log")
ax.set_ylabel("Concentration (mol/m^3)")
ax.set_xscale("symlog", linthresh=1e-6)
ax.set_xlabel("[O2] (M)")
ax.legend()

plt.tight_layout()
plt.savefig("Figures/QCMD_cell_model_max_liquid_radical_concs.pdf", bbox_inches="tight")


def plot_liquid_rop(name, labels, loss_only=False, production_only=False):
    print(f"Plotting {name} ROPs vs. [O2]...")

    fig, axs = plt.subplots(
        nrows=len(perturb_factor_list), ncols=1, figsize=(9, 12), sharex=True
    )

    min_rop = 1e10
    max_rop = 0
    for ind, perturb_factor in enumerate(perturb_factor_list):
        rops, rop_sourcestrings = get_liquid_rops(
            liquid_rop_results[perturb_factor],
            labels,
            radicals_only=True,
            loss_only=loss_only,
            production_only=production_only,
        )
        normalized_rops = rops.abs() / Vliq
        if normalized_rops.empty:
            continue
        min_rop = min(min_rop, min(normalized_rops))
        max_rop = max(max_rop, max(normalized_rops))
        xs = np.arange(len(rop_sourcestrings))
        axs[ind].barh(xs, normalized_rops, align="center")
        axs[ind].set_yticks(xs)
        axs[ind].set_yticklabels(rop_sourcestrings)
        x = xs[ind]
        axs[ind].set_ylabel(f"{x:.1e} M [O2]")
        axs[ind].set_xscale("log")
        axs[ind].invert_yaxis()

    print("min_rop: ", min_rop)
    print("max_rop: ", max_rop)

    for ax in axs:
        ax.set_xlim(min_rop, max_rop)

    if loss_only:
        label = "loss"
    elif production_only:
        label = "production"

    axs[-1].set_xlabel(f"Rate of {name} {label} (mol/(m^3*s))")
    fig.tight_layout()
    fig.savefig(
        f"Figures/QCMD_cell_model_liquid_rop_{label}_{name}.pdf", bbox_inches="tight"
    )


plot_liquid_rop("RC.", R_labels, loss_only=True)
plot_liquid_rop("ROO.", ROO_labels, loss_only=True)
plot_liquid_rop("RO.", RO_labels, loss_only=True)
plot_liquid_rop("RC.", R_labels, production_only=True)
plot_liquid_rop("ROO.", ROO_labels, production_only=True)
plot_liquid_rop("RO.", RO_labels, production_only=True)
