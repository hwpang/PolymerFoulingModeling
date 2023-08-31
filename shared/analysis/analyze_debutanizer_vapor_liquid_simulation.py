# %%
import os
import yaml
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils import get_liquid_rops

# change default font size to 12
plt.rcParams.update({"font.size": 12})


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--simulation_results_directory",
        type=str,
        required=True,
        help="The path to the csv file containing the simulation results.",
    )
    parser.add_argument(
        "--aspen_condition_path",
        type=str,
        required=True,
        help="The path to the yml file containing the Aspen conditions.",
    )
    parser.add_argument(
        "--model_name",
        type=str,
        required=True,
        help="The name of the model.",
    )

    args = parser.parse_args()
    model_name = args.model_name
    aspen_condition_path = args.aspen_condition_path
    simulation_results_directory = args.simulation_results_directory

    return (
        model_name,
        aspen_condition_path,
        simulation_results_directory,
    )


(
    model_name,
    aspen_condition_path,
    simulation_results_directory,
) = parse_arguments()

d = 2.5
h = 0.3
A = (d / 2) ** 2 * np.pi
Vliq = A * h

trays = np.arange(1, 41, 1)

if model_name == "basecase_debutanizer_model":
    delta_t = 64.0
else:
    delta_t = 32.0

tf_liq = np.arange(0, 3600.0 + delta_t, delta_t)[-1]
alpha_rates_path = os.path.join(simulation_results_directory, "alpha_rates.yml")
liquid_simulation_results_path = os.path.join(
    simulation_results_directory, f"simulation_vapor_liquid_yliqn_{tf_liq}.csv"
)
selected_trays = [1, 5, 10, 15, 20, 25, 30, 35, 40]

print("Loading Aspen conditions...")

with open(aspen_condition_path, "r") as f:
    results = yaml.load(f, Loader=yaml.FullLoader)
    Ts = results["T"]
    Ps = results["P"]

print("Loading liquid simulation results...")

liquid_simulation_df = pd.read_csv(liquid_simulation_results_path)

print("Loading liquid rop results...")
liquid_rop_results = {}
for tray in trays:
    path = os.path.join(
        simulation_results_directory, f"simulation_vapor_liquid_liqrop_{tray}.csv"
    )
    liquid_rop_results[tray] = pd.read_csv(path)

print("Loading alpha and rates...")

with open(alpha_rates_path, "r") as f:
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

print(
    "Plotting radical concentrations, chemical lifetime to residence time ratio, consumption rates, and production rates..."
)

carbon_center_radical_mols = liquid_simulation_df.loc[:, R_labels].sum(axis=1)
peroxyl_radical_mols = liquid_simulation_df.loc[:, ROO_labels].sum(axis=1)
alkoxyl_radical_mols = liquid_simulation_df.loc[:, RO_labels].sum(axis=1)

carbon_center_radical_transport_rates = np.array(
    consumption_rates["R._outlet"]
) + np.array(consumption_rates["R._evap"])
carbon_center_radical_reaction_rates = (
    np.array(consumption_rates["R._Add"])
    + np.array(consumption_rates["R._Habs"])
    + np.array(consumption_rates["R._Recomb"])
    + np.array(consumption_rates["R._Disprop"])
    + np.array(consumption_rates["R._CycEther"])
)
carbon_center_radical_chemical_lifetime = (
    carbon_center_radical_mols / carbon_center_radical_reaction_rates
)
carbon_center_radical_residence_time = (
    carbon_center_radical_mols / carbon_center_radical_transport_rates
)

peroxyl_radical_transport_rates = np.array(consumption_rates["ROO._outlet"]) + np.array(
    consumption_rates["ROO._evap"]
)
peroxyl_radical_reaction_rates = (
    np.array(consumption_rates["ROO._Add"])
    + np.array(consumption_rates["ROO._Habs"])
    + np.array(consumption_rates["ROO._Recomb"])
    + np.array(consumption_rates["ROO._Disprop"])
    + np.array(consumption_rates["ROO._eli"])
)
peroxyl_radical_chemical_lifetime = (
    peroxyl_radical_mols / peroxyl_radical_reaction_rates
)
peroxyl_radical_residence_time = peroxyl_radical_mols / peroxyl_radical_transport_rates

alkoxyl_radical_transport_rates = np.array(consumption_rates["RO._outlet"]) + np.array(
    consumption_rates["RO._evap"]
)
alkoxyl_radical_reaction_rates = (
    np.array(consumption_rates["RO._Add"])
    + np.array(consumption_rates["RO._Habs"])
    + np.array(consumption_rates["RO._Recomb"])
    + np.array(consumption_rates["RO._Disprop"])
    + np.array(consumption_rates["RO._CycEther"])
)
alkoxyl_radical_chemical_lifetime = (
    alkoxyl_radical_mols / alkoxyl_radical_reaction_rates
)
alkoxyl_radical_residence_time = alkoxyl_radical_mols / alkoxyl_radical_transport_rates

os.makedirs("Figures", exist_ok=True)

fig = plt.figure(figsize=(10, 9))
gs = fig.add_gridspec(3, 2)

ax = fig.add_subplot(gs[0, 0])
concs = carbon_center_radical_mols / Vliq
ax.plot(trays, concs, "-o", label="RC.")
concs = peroxyl_radical_mols / Vliq
ax.plot(trays, concs, "-x", label="ROO.")
concs = alkoxyl_radical_mols / Vliq
ax.plot(trays, concs, "-s", label="RO.")
ax.set_ylabel(r"R.(liq) (mol/m$^3$)", fontsize=12)
ax.set_xlabel("Tray", fontsize=12)
ax.set_yscale("log")
ax.set_title("(a)", loc="left")

ax2 = ax.twinx()
ax2.plot(trays, Ts, "m-v")
ax2.set_ylabel("Temperature (K)", fontsize=12, color="m")

ax = fig.add_subplot(gs[0, 1])
ax.plot(
    trays,
    carbon_center_radical_chemical_lifetime / carbon_center_radical_residence_time,
    "-o",
    label="RC.",
)
ax.plot(
    trays,
    peroxyl_radical_chemical_lifetime / peroxyl_radical_residence_time,
    "-x",
    label="ROO.",
)
ax.plot(
    trays,
    alkoxyl_radical_chemical_lifetime / alkoxyl_radical_residence_time,
    "-s",
    label="RO.",
)
ax.set_ylabel(
    "R.(liq) " + r"$\tau_\mathrm{chem}$" + "/" + r"$\tau_\mathrm{res}$", fontsize=12
)
ax.set_xlabel("Tray", fontsize=12)
ax.set_yscale("log")
ax.set_title("(b)", loc="left")

patterns = ["//", "\\\\", "||", "--", "++", "xx", "oo", "OO", "..", "**"]

ax = fig.add_subplot(gs[1, :])
bottom = np.zeros(len(trays))
consumption_paths = [
    "R._outlet",
    "R._evap",
    "R._Add",
    "R._Habs",
    "R._Recomb",
    "R._Disprop",
    "R.+O2",
    "R._CycEther",
    "ROO._outlet",
    "ROO._evap",
    "ROO._Add",
    "ROO._Habs",
    "ROO._Recomb",
    "ROO._Disprop",
    "ROO._eli",
    "RO._outlet",
    "RO._evap",
    "RO._Add",
    "RO._Habs",
    "RO._Recomb",
    "RO._Disprop",
    "RO._CycEther",
]
total_consumption_rates = np.sum(
    [consumption_rates[path] for path in consumption_paths], axis=0
)
count = 0
for consumption_path in consumption_paths:
    percentages = consumption_rates[consumption_path] / total_consumption_rates * 100
    if any(percentages > 5):
        if "R." in consumption_path and model_name != "basecase_debutanizer_model":
            label = consumption_path.replace("R.", "RC.")
        else:
            label = consumption_path
        ax.bar(
            trays,
            percentages,
            bottom=bottom,
            label=label,
            hatch=patterns[count],
            width=1.0,
        )
        count += 1
        bottom += percentages
ax.bar(trays, 100.0 - bottom, bottom=bottom, label="other")
ax.set_ylabel("R.(liq) consumption (%)")
ax.set_ylim([0, 100])
ax.set_title("(c)", loc="left")
ax.set_xticks([])
ax.set_xticklabels([])
ax.legend(bbox_to_anchor=(1, 1))

ax = fig.add_subplot(gs[2, :])
bottom = np.zeros(len(trays))
production_paths = [
    "R._inlet",
    "R._cond",
    "R._RevDisprop",
    "ROO._inlet",
    "ROO._cond",
    "RO._inlet",
    "RO._cond",
    "RO._BondDiss",
]
total_production_rates = np.sum(
    [production_rates[path] for path in production_paths], axis=0
)
count = 0
for production_path in production_paths:
    percentages = production_rates[production_path] / total_production_rates * 100
    if any(percentages > 5):
        if "R." in production_path and model_name != "basecase_debutanizer_model":
            label = production_path.replace("R.", "RC.")
        else:
            label = production_path
        ax.bar(
            trays,
            percentages,
            bottom=bottom,
            label=label,
            hatch=patterns[count],
            width=1.0,
        )
        count += 1
        bottom += percentages
if model_name != "basecase_debutanizer_model":
    ax.bar(trays, 100.0 - bottom, bottom=bottom, label="other")
ax.set_ylabel("R.(liq) production (%)")
ax.set_ylim([0, 100])
ax.set_title("(d)", loc="left")
ax.legend(bbox_to_anchor=(1, 1))
ax.set_xlabel("Tray", fontsize=12)

plt.subplots_adjust(wspace=0, hspace=0)
fig.tight_layout()
fig.savefig(f"Figures/{model_name}_liquid_radical.pdf", bbox_inches="tight")

print("Plotting alphas and ASF distributions...")

cmap = plt.get_cmap("plasma")
cs = cmap(trays / len(trays))

fig = plt.figure(figsize=(10, 5))

if model_name == "basecase_debutanizer_model":
    gs = fig.add_gridspec(2, 3)

    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(trays, alphas_DA, c=cs, zorder=2)
    ax.plot(trays, alphas_DA, zorder=1, color="grey")
    ax.set_ylabel(r"$\alpha_\mathrm{DA}$", fontsize=13)
    ax.set_xlabel("Tray")
    ax.set_title("(a)", loc="left")
    ax.set_yscale("log")

    ax = fig.add_subplot(gs[1, 0])
    ax.scatter(trays, alphas, c=cs, zorder=2)
    ax.plot(trays, alphas, zorder=1, color="grey")
    ax.set_ylim([0, 1])
    ax.set_ylabel(r"$\alpha_\mathrm{R.}$", fontsize=13)
    ax.set_xlabel("Tray")
    ax.set_title("(b)", loc="left")

else:
    gs = fig.add_gridspec(3, 3)

    ax = fig.add_subplot(gs[0, 0])
    alpha1_RC_add = alpha1 * (
        consumption_rates["R._Add"]
        / (consumption_rates["R._Add"] + consumption_rates["R.+O2"])
    )
    alpha1_RC_O2 = alpha1 * (
        consumption_rates["R.+O2"]
        / (consumption_rates["R._Add"] + consumption_rates["R.+O2"])
    )
    ax.scatter(
        trays,
        alpha1_RC_add,
        c=cs,
        zorder=2,
        label=r"$\alpha_\mathrm{RC. add}$",
        marker="s",
    )
    ax.plot(trays, alpha1_RC_add, zorder=1, color="grey")
    ax.scatter(
        trays,
        alpha1_RC_O2,
        c=cs,
        zorder=2,
        label=r"$\alpha_\mathrm{RC.+O2}$",
        marker="x",
    )
    ax.plot(trays, alpha1_RC_O2, zorder=1, color="grey")
    ax.set_ylim([0, 1])
    ax.set_ylabel(r"$\alpha_\mathrm{RC.}$", fontsize=13)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_title("(a)", loc="left")
    ax.legend(bbox_to_anchor=(1, 1))

    ax = fig.add_subplot(gs[1, 0])
    ax.scatter(trays, alpha2, c=cs, zorder=2)
    ax.plot(trays, alpha2, zorder=1, color="grey")
    ax.set_ylim([0, 1])
    ax.set_ylabel(r"$\alpha_\mathrm{ROO.}$", fontsize=13)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_title("(b)", loc="left")

    ax = fig.add_subplot(gs[2, 0])
    ax.scatter(trays, alphas, c=cs, zorder=2)
    ax.plot(trays, alphas, zorder=1, color="grey")
    ax.set_ylim([0, 1])
    ax.set_ylabel(r"$\alpha_\mathrm{R.}$", fontsize=13)
    ax.set_xlabel("Tray")
    ax.set_title("(c)", loc="left")

ax = fig.add_subplot(gs[:, 1:], projection="3d")

if model_name == "basecase_debutanizer_model":
    ax.set_title("(c)", loc="left")
    ns = np.arange(1, 16)
else:
    ax.set_title("(d)", loc="left")
    ns = np.arange(1, 40)


def Wn(n, alpha):
    return n * (1 - alpha) ** 2 * alpha ** (n - 1)


for n in ns:
    cs = cmap(trays / len(trays))
    ax.bar(trays, Wn(n, alphas), zs=n, zdir="x", color=cs)
ax.set_xlabel(r"$k$")
ax.set_ylabel("Tray", labelpad=10)
ax.set_zlabel(r"W($k$)", labelpad=10)

fig.tight_layout()
fig.savefig(f"Figures/{model_name}_ASF_distribution.pdf", bbox_inches="tight")


def get_max_conc_liq_radical_label(df, tray, labels, N=1):
    radical_mols = df.loc[tray - 1, labels]
    inds = radical_mols.argsort()[::-1][:N]
    return [labels[ind] for ind in inds]


print("Plotting max liquid radical concentration...")
max_conc_liq_radical_labels = set(
    [
        label
        for tray in trays
        for label in get_max_conc_liq_radical_label(
            liquid_simulation_df, tray, R_labels + ROO_labels + RO_labels, N=1
        )
    ]
)
print(max_conc_liq_radical_labels)

fig, ax = plt.subplots()
for label in max_conc_liq_radical_labels:
    ax.plot(trays, liquid_simulation_df.loc[:, label], label=label)
ax.set_yscale("log")
ax.set_ylabel("Concentration (mol/m$^3$)")
ax.set_xlabel("Tray")
ax.legend(bbox_to_anchor=(1, 1))
fig.tight_layout()
fig.savefig(f"Figures/{model_name}_max_conc_liq_radical.pdf", bbox_inches="tight")

print("Plotting liquid radical rop...")


def plot_liquid_rop(name, labels, loss_only=False, production_only=False):
    print(f"Plotting {name} ROPs vs. tray...")

    fig, axs = plt.subplots(
        nrows=len(selected_trays), ncols=1, figsize=(9, 12), sharex=True
    )

    min_rop = 1e10
    max_rop = 0
    for ind, tray in enumerate(selected_trays):
        rops, rop_sourcestrings = get_liquid_rops(
            liquid_rop_results[tray],
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
        axs[ind].set_ylabel(f"Tray {tray}")
        axs[ind].set_xscale("log")
        axs[ind].invert_yaxis()

    print("min_rop: ", min_rop)
    print("max_rop: ", max_rop)

    for ax in axs:
        ax.set_xlim(min_rop, max_rop)

    axs[-1].set_xlabel(f"Rate of {name} loss (mol/(m^3*s))")
    fig.tight_layout()

    if loss_only:
        label = "loss"
    elif production_only:
        label = "production"

    fig.savefig(
        f"Figures/{model_name}_liquid_rop_{label}_{name}.pdf", bbox_inches="tight"
    )


plot_liquid_rop("RC.", R_labels, loss_only=True)
plot_liquid_rop("ROO.", ROO_labels, loss_only=True)
plot_liquid_rop("RO.", RO_labels, loss_only=True)

plot_liquid_rop("RC.", R_labels, production_only=True)
plot_liquid_rop("ROO.", ROO_labels, production_only=True)
plot_liquid_rop("RO.", RO_labels, production_only=True)
