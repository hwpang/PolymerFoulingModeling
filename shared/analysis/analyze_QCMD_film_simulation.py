import os
import yaml
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils import get_liq_radicals_conc, get_film_radical_reactive_site_conc, get_film_rops

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
print("alls_imulation_directory: ",all_simulation_directory)

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

# experimental film growth rates
expt_rate_dict = {}
expt_rate_dict["Unsparged"] = 9.4  # ng/hr
expt_rate_dict["N2 sparged"] = 6.1
expt_rate_dict["O2 sparged"] = 927.3

print("Loading liquid simulation results...")
liquid_simulations = dict()
for perturb_factor in perturb_factor_list:
    liquid_simulation_path = os.path.join(all_simulation_directory, f"{perturb_species}_{perturb_factor}", f"simulation_liquid_{tray}.csv")
    liquid_simulations[perturb_factor] = pd.read_csv(liquid_simulation_path)

print("Loading liquid rops...")
liquid_rops = dict()
for perturb_factor in perturb_factor_list:
    liquid_rop_path = os.path.join(all_simulation_directory, f"{perturb_species}_{perturb_factor}", f"simulation_liquid_liqrop_{tray}.csv")
    liquid_rops[perturb_factor] = pd.read_csv(liquid_rop_path)

print("Loading alpha rates...")
with open(os.path.join(all_simulation_directory, f"{perturb_species}_{perturb_factor_list[0]}", f"alpha_rates.yml"), "r") as f:
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

print("Loading film simulation results...")
film_simulations = dict()
for perturb_factor in perturb_factor_list:
    film_simulation_path = os.path.join(all_simulation_directory, f"{perturb_species}_{perturb_factor}", f"simulation_film_{tray}.csv")
    film_simulations[perturb_factor] = pd.read_csv(film_simulation_path)

print("Loading film rops...")
film_rops = dict()
for perturb_factor in perturb_factor_list:
    film_rop_path = os.path.join(all_simulation_directory, f"{perturb_species}_{perturb_factor}", f"simulation_film_rop_{tray}.csv")
    film_rops[perturb_factor] = pd.read_csv(film_rop_path)

print("Plotting film growth rates vs. [O2]...")

# pure_O2_sat_conc = 0.008093871600706785 #M estimated by RMG
pure_O2_sat_molfrac = 8.1 * 1e-4  # mol fraction measured experimentally from Battino, R., Rettich, T. R., & Tominaga, T. (1983). The Solubility of Oxygen and Ozone in Liquids. Journal of Physical and Chemical Reference Data, Vol. 12, pp. 163–178. https://doi.org/10.1063/1.555680
benzene_density = 876  # kg/m^3 from Lide, D. R., Data, S. R., Board, E. A., Baysinger, G., Chemistry, S., Library, C. E., … Zwillinger, D. (2004). CRC Handbook of Chemistry and Physics. 2660.
benzene_mw = 78.11 / 1000  # kg/mol from Lide, D. R., Data, S. R., Board, E. A., Baysinger, G., Chemistry, S., Library, C. E., … Zwillinger, D. (2004). CRC Handbook of Chemistry and Physics. 2660.
benzene_sat_conc = benzene_density / benzene_mw  # mol/m^3
pure_O2_sat_conc = pure_O2_sat_molfrac * benzene_sat_conc  # mol/m^3
air_O2_sat_conc = 0.21 * pure_O2_sat_conc  # mol/m^3
air_O2_sat_conc /= 1000  # mol/L

xs = np.array(factor_num_list) * air_O2_sat_conc

def calculate_film_growth_rate(df):
    inds = range(1, len(df.index))
    dmdts = np.array([(df.loc[ind, "mass"] - df.loc[0, "mass"]) / (df.loc[ind, "timestamp"] - df.loc[0, "timestamp"])  for ind in inds])
    return np.mean(dmdts), np.std(dmdts)

plt.figure(figsize=(4, 3))

# simulated film growth rate
film_growth_rates = [calculate_film_growth_rate(film_simulations[perturb_factor]) for perturb_factor in perturb_factor_list]
print("film_growth_rates: ", film_growth_rates)
ys = np.array([film_growth_rate[0] for film_growth_rate in film_growth_rates])
ys += ys / rho * epsilon * rho_liq  # converting from mass of solid to mass of film by adding mass of liquid in film
ys *= 1e9 * 1e3 * 3600  # kg/s to ng/hr
yerr = np.exp(
    2 * 1000 * 4.184 / 8.314 / (273.15 + 90)
)  # 2 kcal/mol uncertainty in activation energy
label = "Prediction"
plt.errorbar(xs, ys, yerr=[ys - ys * 1 / yerr, -ys + ys * yerr], color="C0", marker="o", capsize=5, label=label)
plt.fill_between(xs, ys * 1 / yerr, ys * yerr, color="C0", alpha=0.5, linewidth=0)

expt_rate_std = 1.635621377
# unsparged
label = "Unsparged"
expt_rate_mean = expt_rate_dict[label]  # ng/hr
xerr = 3
expt_o2 = air_O2_sat_conc / xerr  # M
yerr = expt_rate_std
color = "C1"
plt.errorbar([expt_o2], [expt_rate_mean], yerr=[yerr], xerr=[[expt_o2 - expt_o2 / xerr], [-expt_o2 + expt_o2 * xerr]], color=color, marker="o", capsize=5, label=label)
plt.fill_between([expt_o2 / xerr, expt_o2 * xerr], [expt_rate_mean - yerr, expt_rate_mean - yerr], [expt_rate_mean + yerr, expt_rate_mean + yerr], color=color, alpha=0.5, linewidth=0)

# N2 sparged
label = "N2 sparged"
expt_rate_mean = expt_rate_dict[label]  # ng/hr
xerr = 10
reduction_factor = 10
expt_o2 = expt_o2 / reduction_factor  # M
yerr = expt_rate_std
color = "C2"
plt.errorbar(
    [expt_o2],
    [expt_rate_mean],
    yerr=[yerr],
    xerr=[[expt_o2 - expt_o2 / xerr], [-expt_o2 + expt_o2 * xerr]],
    color=color,
    marker="o",
    capsize=5,
    label="$N_2$ sparged",
)
plt.fill_between(
    [expt_o2 / xerr, expt_o2 * xerr],
    [expt_rate_mean - yerr, expt_rate_mean - yerr],
    [expt_rate_mean + yerr, expt_rate_mean + yerr],
    color=color,
    alpha=0.5,
    linewidth=0,
)
plt.ylim([1e-2, 1e4])
plt.xlabel("[$O_2$] (M)")
plt.ylabel("Film growth rate (ng/hr)")
plt.xscale("symlog", linthresh=1e-6)
plt.yscale("log")
# plt.ylim([1e-3, 1e2])
plt.tight_layout()
plt.legend(fontsize=9)

os.makedirs("Figures", exist_ok=True)

plt.savefig("Figures/QCMD_cell_model_film_growth_rates.pdf", bbox_inches="tight")
plt.close()


print("Plotting radical concentrations vs. [O2]...")

fig, axs = plt.subplots(1, 2, figsize=(8, 3))

ax = axs[0]
ax.plot(xs, [get_liq_radicals_conc(liquid_simulations[perturb_factor], R_labels, Vliq) for perturb_factor in perturb_factor_list], marker="o", label="R.")
ax.plot(xs, [get_liq_radicals_conc(liquid_simulations[perturb_factor], ROO_labels, Vliq) for perturb_factor in perturb_factor_list], marker="o", label="ROO.")
ax.plot(xs, [get_liq_radicals_conc(liquid_simulations[perturb_factor], RO_labels, Vliq) for perturb_factor in perturb_factor_list], marker="o", label="RO.")
ax.set_yscale("log")
ax.set_ylabel("Concentration (mol/m^3)")
ax.set_xscale("symlog", linthresh=1e-6)
ax.set_xlabel("[O2] (M)")
ax.legend()

ax = axs[1]
ax.plot(xs, [get_film_radical_reactive_site_conc(film_simulations[perturb_factor], "AR") for perturb_factor in perturb_factor_list], marker="o", label="AR")
ax.plot(xs, [get_film_radical_reactive_site_conc(film_simulations[perturb_factor], "KR") for perturb_factor in perturb_factor_list], marker="o", label="KR")
ax.plot(xs, [get_film_radical_reactive_site_conc(film_simulations[perturb_factor], "PR") for perturb_factor in perturb_factor_list], marker="o", label="PR")
ax.plot(xs, [get_film_radical_reactive_site_conc(film_simulations[perturb_factor], "OR") for perturb_factor in perturb_factor_list], marker="o", label="OR")
ax.set_yscale("log")
ax.set_ylabel("Concentration (mol/kg)")
ax.set_xscale("symlog", linthresh=1e-6)
ax.set_xlabel("[O2] (M)")
ax.legend()

plt.tight_layout()
plt.savefig("Figures/QCMD_cell_model_radical_concs.pdf", bbox_inches="tight")

print("Plotting mass ROPs vs. [O2]...")

fig, axs = plt.subplots(nrows=len(perturb_factor_list), ncols=1, figsize=(9, 12), sharex=True)

min_rop = 1e10
max_rop = 0
for ind, perturb_factor in enumerate(perturb_factor_list):
    rops, rop_rxncomments, rop_rxnstrs = get_film_rops(film_rops[perturb_factor], "mass")
    df = film_simulations[perturb_factor]
    mass = df.loc[len(df.index)-1, "mass"]
    normalized_rops = rops / mass
    min_rop = min(min_rop, min(normalized_rops))
    max_rop = max(max_rop, max(normalized_rops))
    ys = np.arange(len(rop_rxncomments))
    axs[ind].barh(xs, normalized_rops, align="center")
    axs[ind].set_yticks(ys)
    axs[ind].set_yticklabels(rop_rxncomments)
    x = xs[ind]
    axs[ind].set_ylabel(f"[$O_2$] = {x:.1e} M")
    axs[ind].set_xscale("log")
    axs[ind].invert_yaxis()

for ax in axs:
    ax.set_xlim(min_rop, max_rop)

axs[-1].set_xlabel("Rate of film growth (kg/(kg*s))")
fig.tight_layout()
fig.savefig(f"Figures/QCMD_cell_model_film_rop_mass.pdf", bbox_inches="tight")


def plot_film_rop(name, loss_only=False, production_only=False):
    print(f"Plotting {name} ROPs vs. [O2]...")

    fig, axs = plt.subplots(nrows=len(perturb_factor_list), ncols=1, figsize=(9, 12), sharex=True)

    min_rop = 1e10
    max_rop = 0
    for ind, perturb_factor in enumerate(perturb_factor_list):
        if (name == "PR" or name == "OR") and perturb_factor == "0.0":
            continue
        rops, rop_rxncomments, rop_rxnstrs = get_film_rops(film_rops[perturb_factor], name, loss_only=loss_only, production_only=production_only)
        df = film_simulations[perturb_factor]
        mass = df.loc[len(df.index)-1, "mass"]
        normalized_rops = rops.abs() / mass
        min_rop = min(min_rop, min(normalized_rops))
        max_rop = max(max_rop, max(normalized_rops))
        ys = np.arange(len(rop_rxncomments))
        axs[ind].barh(ys, normalized_rops, align="center")
        axs[ind].set_yticks(ys)
        axs[ind].set_yticklabels(rop_rxncomments)
        x = xs[ind]
        axs[ind].set_ylabel(f"{x:.1e} M [$O_2$]")
        axs[ind].set_xscale("log")
        axs[ind].invert_yaxis()

    for ax in axs:
        ax.set_xlim(min_rop, max_rop)
        
    if loss_only:
        label = "loss"
    elif production_only:
        label = "production"

    axs[-1].set_xlabel(f"Rate of {name} {label} (kg/(kg*s))")
    fig.tight_layout()
    fig.savefig(f"Figures/QCMD_cell_model_film_rop_{label}_{name}.pdf", bbox_inches="tight")

plot_film_rop("AR", loss_only=True)
plot_film_rop("KR", loss_only=True)
plot_film_rop("PR", loss_only=True)
plot_film_rop("OR", loss_only=True)

plot_film_rop("AR", production_only=True)
plot_film_rop("KR", production_only=True)
plot_film_rop("PR", production_only=True)
plot_film_rop("OR", production_only=True)