# %%
import sys
sys.path.insert(0, "/home/gridsan/hwpang/Software/RMG-Py/")

import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"

import yaml
import argparse
import numpy as np
import matplotlib.pyplot as plt
#change default font size to 12
plt.rcParams.update({'font.size': 12})

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--alpha_rates_path", type=str, required=True, help="The path to the yml file containing alpha and rates.",
    )
    parser.add_argument(
        "--model_name", type=str, required=True, help="The name of the model.",
    )

    args = parser.parse_args()
    alpha_rates_path = args.alpha_rates_path
    model_name = args.model_name

    return (
        alpha_rates_path,
        model_name,
    )

(
    alpha_rates_path,
    model_name,
) = parse_arguments()

with open(alpha_rates_path, 'r') as f:
    results = yaml.load(f, Loader=yaml.FullLoader)
    all_alphas, consumption_rates, production_rates, rxn_rates, radical_concs = results
    alpha1, alpha2, alphas, alphas_DA = all_alphas

d = 2.5
h = 0.3
A = (d / 2) ** 2 * np.pi
Vliq = A * h

trays = np.arange(0, 40, 1)
carbon_center_radical_mols = np.sum(list(radical_concs[0].values()), axis=0)
peroxyl_radical_mols = np.sum(list(radical_concs[1].values()), axis=0)
alkoxyl_radical_mols = np.sum(list(radical_concs[2].values()), axis=0)

carbon_center_radical_transport_rates = np.array(consumption_rates["R._outlet"]) + np.array(consumption_rates["R._evap"])
carbon_center_radical_reaction_rates = np.array(consumption_rates["R._Add"]) + np.array(consumption_rates["R._Habs"]) + np.array(consumption_rates["R._Recomb"]) + np.array(consumption_rates["R._Disprop"]) + np.array(consumption_rates["R._CycEther"])
carbon_center_radical_chemical_lifetime = carbon_center_radical_mols / carbon_center_radical_reaction_rates
carbon_center_radical_residence_time = carbon_center_radical_mols / carbon_center_radical_transport_rates

peroxyl_radical_transport_rates = np.array(consumption_rates["ROO._outlet"]) + np.array(consumption_rates["ROO._evap"])
peroxyl_radical_reaction_rates = np.array(consumption_rates["ROO._Add"]) + np.array(consumption_rates["ROO._Habs"]) + np.array(consumption_rates["ROO._Recomb"]) + np.array(consumption_rates["ROO._Disprop"]) + np.array(consumption_rates["ROO._eli"])
peroxyl_radical_chemical_lifetime = peroxyl_radical_mols / peroxyl_radical_reaction_rates
peroxyl_radical_residence_time = peroxyl_radical_mols / peroxyl_radical_transport_rates

alkoxyl_radical_transport_rates = np.array(consumption_rates["RO._outlet"]) + np.array(consumption_rates["RO._evap"])
alkoxyl_radical_reaction_rates = np.array(consumption_rates["RO._Add"]) + np.array(consumption_rates["RO._Habs"]) + np.array(consumption_rates["RO._Recomb"]) + np.array(consumption_rates["RO._Disprop"]) + np.array(consumption_rates["RO._CycEther"])
alkoxyl_radical_chemical_lifetime = alkoxyl_radical_mols / alkoxyl_radical_reaction_rates
alkoxyl_radical_residence_time = alkoxyl_radical_mols / alkoxyl_radical_transport_rates

#%% Plot the radical concentration, radical lifetime to radical radical residence time ratio, consumption rate distribution, and production rate distribution
fig = plt.figure(figsize=(9, 9))
gs = fig.add_gridspec(3, 2)

ax = fig.add_subplot(gs[0, 0])
ax.plot(trays, carbon_center_radical_mols / Vliq, "-o", label="RC.")
concs = peroxyl_radical_mols / Vliq
ax.plot(trays, concs, "-x", label="ROO.")
concs = alkoxyl_radical_mols / Vliq
ax.plot(trays, concs, "-s", label="RO.")
ax.set_ylabel("R.(liq) (mol/m^3)", fontsize=12)
ax.set_xlabel("Tray", fontsize=12)
ax.set_yscale("log")
ax.set_title("(a)", loc="left")

ax = fig.add_subplot(gs[0, 1])
ax.plot(trays, carbon_center_radical_chemical_lifetime / carbon_center_radical_residence_time, "-o", label="RC.")
ax.plot(trays, peroxyl_radical_chemical_lifetime / peroxyl_radical_residence_time, "-x", label="ROO.")
ax.plot(trays, alkoxyl_radical_chemical_lifetime / alkoxyl_radical_residence_time, "-s", label="RO.")
ax.set_ylabel("R.(liq) $\tau_\mathrm{chem}$/$\tau_\mathrm{res}$", fontsize=12)
ax.set_xlabel("Tray", fontsize=12)
ax.set_yscale("log")
ax.set_title("(b)", loc="left")

ax = fig.add_subplot(gs[1, :])
bottom = np.zeros(len(trays))
consumption_paths = consumption_rates.keys()
total_consumption_rates = np.sum(list(consumption_rates.values()), axis=0)
for consumption_path in consumption_paths:
    if any(consumption_rates[consumption_path]/total_consumption_rates > 0.02):
        if "R." in consumption_path:
            label = consumption_path.replace("R.", "RC.")
        else:
            label = consumption_path
        ax.bar(trays, consumption_rates[consumption_path], bottom=bottom, label=label)
        bottom += consumption_rates[consumption_path]
ax.bar(trays, 100.0 - bottom, bottom=bottom, label="other")
ax.set_ylabel("R.(liq) consumption (%)")
ax.set_ylim([0, 100])
ax.set_title("(c)", loc="left")
ax.set_xticks([])
ax.set_xticklabels([])
ax.legend(bbox_to_anchor=(1, 1))

ax = fig.add_subplot(gs[2, :])
bottom = np.zeros(len(trays))
production_paths = production_rates.keys()
total_production_rates = np.sum(list(production_rates.values()), axis=0)
for production_path in production_paths:
    if any(production_rates[production_path]/total_production_rates > 0.02):
        if "R." in production_path:
            label = production_path.replace("R.", "RC.")
        else:
            label = production_path
        ax.bar(trays, production_rates[production_path], bottom=bottom, label=label)
        bottom += production_rates[production_path]
ax.bar(trays, 100.0 - bottom, bottom=bottom, label="other")
ax.set_ylabel("R.(liq) production (%)")
ax.set_ylim([0, 100])
ax.set_title("(d)", loc="left")
ax.legend(bbox_to_anchor=(1, 1))

fig.tight_layout()
fig.savefig(f"{model_name}_liquid_radical.pdf", bbox_inches="tight")

#%% Plot alpha_DA, alpha_R, and Wn
cmap = plt.get_cmap("plasma")
cs = cmap(trays / len(trays))

fig = plt.figure(figsize=(10, 5))
gs = fig.add_gridspec(3, 3)

# ax = fig.add_subplot(gs[0, 0])
# ax.scatter(trays, alphas_DA, c=cs, zorder=2)
# ax.plot(trays, alphas_DA, zorder=1, color="grey")
# ax.set_ylabel("$\alpha_\mathrm{DA}$", fontsize=13)
# ax.set_xlabel("Tray")
# ax.set_title("(a)", loc="left")
# ax.set_yscale("log")

ax = fig.add_subplot(gs[0, 0])
ax.scatter(trays, alpha1, c=cs, zorder=2)
ax.plot(trays, alphas_DA, zorder=1, color="grey")
ax.set_ylabel("$\alpha_\mathrm{RC.}$", fontsize=13)
ax.set_xlabel("Tray")
ax.set_title("(a)", loc="left")
ax.set_yscale("log")

ax = fig.add_subplot(gs[1, 0])
ax.scatter(trays, alpha2, c=cs, zorder=2)
ax.plot(trays, alphas, zorder=1, color="grey")
ax.set_ylim([0, 1])
ax.set_ylabel("$\alpha_\mathrm{ROO.}$", fontsize=13)
ax.set_xlabel("Tray")
ax.set_title("(b)", loc="left")

ax = fig.add_subplot(gs[2, 0])
ax.scatter(trays, alphas, c=cs, zorder=2)
ax.plot(trays, alphas, zorder=1, color="grey")
ax.set_ylim([0, 1])
ax.set_ylabel("$\alpha_\mathrm{R.}$", fontsize=13)
ax.set_xlabel("Tray")
ax.set_title("(c)", loc="left")

ax = fig.add_subplot(gs[:, 1:], projection="3d")
ns = np.arange(1, 40)
def Wn(n, alpha):
    return n * (1 - alpha) ** 2 * alpha ** (n - 1)

for n in ns:
    cs = cmap(trays / len(trays))
    ax.bar(trays, Wn(n, alphas), zs=n, zdir="x", color=cs)
ax.set_xlabel("$k$")
ax.set_ylabel("Tray", labelpad=10)
ax.set_zlabel("W($k$)", labelpad=10)
ax.set_title("(d)", loc="left")

fig.tight_layout()
fig.savefig(f"{model_name}_ASF_distribution.pdf", bbox_inches="tight")







