import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# change default font size to 12
plt.rcParams.update({"font.size": 20})


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model_name",
        type=str,
        required=True,
        help="The name of the model.",
    )
    parser.add_argument(
        "--aspen_results_path",
        type=str,
        required=True,
        help="The path to the Aspen results file.",
    )

    args = parser.parse_args()
    model_name = args.model_name
    aspen_results_path = args.aspen_results_path

    return (
        model_name,
        aspen_results_path,
    )


(
    model_name,
    aspen_results_path,
) = parse_arguments()

# Read TPFQ data
TPFQ = pd.read_excel(aspen_results_path, engine='openpyxl', sheet_name="TPFQ")

# Read liquid molar fraction data
liquid_molfraction = pd.read_excel(aspen_results_path, engine='openpyxl', sheet_name="LiquidCompositionMoleBasis")

# Read vapor molar fraction data
vapor_molfraction = pd.read_excel(aspen_results_path, engine='openpyxl', sheet_name="VaporCompositionMoleBasis")

# Read stream data
stream = pd.read_excel(aspen_results_path, engine='openpyxl', sheet_name="StreamResults")

# Define species name dictionary
spcnamedict = {
    "N-BUT-01": "N-BUTANE",
    "2-BUT-01": "2-BUTENE",
    "1:3-B-01": "1,3-BUTADIENE",
    "CYCLO-01": "CYCLOPENTADIENE",
    "BENZE-01": "BENZENE",
    "1:3-C-01": "1,3-CYCLOHEXADIENE",
    "TOLUE-01": "TOLUENE",
    "STYRE-01": "STYRENE"
}
spcnames = [
    "N-BUTANE",
    "2-BUTENE",
    "1,3-BUTADIENE",
    "CYCLOPENTADIENE",
    "BENZENE",
    "1,3-CYCLOHEXADIENE",
    "TOLUENE",
    "STYRENE"
]

# Add "OXYGEN" to the species name dictionary if model_name is "trace_oxygen_perturbed_debutanizer_model"
if model_name == "trace_oxygen_perturbed_debutanizer_model":
    spcnamedict["OXYGEN"] = "OXYGEN"
    spcnames.append("OXYGEN")

# Initialize initial conditions dictionary
initial_conditions = {}

# Calculate liquid outlet volumetric flowrate
liquid_molar_density = stream.loc[13, "DC4-L"] * 1000  # mol/m^3
liquid_concentration = liquid_molfraction.iloc[:, 1:] * liquid_molar_density
liquid_outlet_molar_flowrate = TPFQ.loc[1:, "Liquid from (Mole)"] * 1000  # mol/s
liquid_outlet_volumetric_flowrate = liquid_outlet_molar_flowrate / liquid_molar_density
initial_conditions["liquid_outlet_volumetric_flowrate"] = liquid_outlet_volumetric_flowrate.tolist()

# Calculate vapor outlet volumetric flowrate
vapor_molar_density = stream.loc[13, "DC4-V"] * 1000  # mol/m^3
vapor_concentration = vapor_molfraction.iloc[:, 1:] * vapor_molar_density
vapor_outlet_molar_flowrate = TPFQ.loc[1:, "Vapor from (Mole)"] * 1000  # mol/s
vapor_outlet_volumetric_flowrate = vapor_outlet_molar_flowrate / vapor_molar_density
initial_conditions["vapor_outlet_volumetric_flowrate"] = vapor_outlet_volumetric_flowrate.tolist()

# Extract temperature and pressure data
T = TPFQ.loc[1:, "Temperature"].tolist()
initial_conditions["T"] = T
P = TPFQ.loc[1:, "Pressure"].tolist()
initial_conditions["P"] = P

# Initialize liquid and vapor concentration dictionaries
initial_conditions["liquid_concentration"] = {}
initial_conditions["vapor_concentration"] = {}

# Extract liquid and vapor concentrations for each species
for aspenname, spcname in spcnamedict.items():
    conc = liquid_concentration[aspenname]
    conc[conc < 1e-18] = 0.0
    initial_conditions["liquid_concentration"][spcname] = conc.tolist()
    conc = vapor_concentration[aspenname]
    conc[conc < 1e-18] = 0.0
    initial_conditions["vapor_concentration"][spcname] = conc.tolist()

# Set save_name based on model_name
if model_name == "basecase_debutanizer_model":
    save_name = "aspen_conditions"
elif model_name == "trace_oxygen_perturbed_debutanizer_model":
    save_name = "aspen_conditions_oxygen"

# Write initial conditions to YAML file
import yaml
with open(f"{save_name}.yml", "w") as f:
    yaml.dump(initial_conditions, f)

# Define plot parameters
d = 2.5
h = 0.3
A = (d / 2) ** 2 * np.pi
Vliq = A * h
spacing = 0.6
Vvap = A * (spacing - h)

lines = ["-v", "-^", "-<", "->", "-1", "-2", "-3", "-4", "-8", "-s", "-p", "-P", "-*", "-h", "-H", "-+", "-x", "-X", "-D", "-d"]

# Create subplots
if model_name == "basecase_debutanizer_model":
    fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(10, 14))
elif model_name == "trace_oxygen_perturbed_debutanizer_model":
    fig, axs = plt.subplots(nrows=5, ncols=2, figsize=(10, 17))

# Adjust subplot spacing
fig.subplots_adjust(wspace=0, hspace=0)

# Plot temperature data
axs[0, 0].plot(range(1, 41), T, color="k", marker="o")
axs[0, 0].set_ylabel("Temperature (K)", size=20)
axs[0, 0].set_title("(a)", loc="left", size=20)
axs[0, 0].set_xticklabels([])

# Plot liquid concentration data
for i, spc in enumerate(spcnames):
    label = spc
    axs[1, 0].plot(range(1, 41), initial_conditions["liquid_concentration"][spc], lines[i], label=label)

axs[1, 0].set_ylabel("Conc. (mol/m^3)", size=20)
axs[1, 0].set_title("(b) Liquid phase", loc="left", size=20)
axs[1, 0].set_xticklabels([])
axs[1, 0].plot([0, 0], [-400, 800], "k--")
axs[1, 0].plot([0, 41], [-400, -400], "k--")
axs[1, 0].plot([41, 41], [-400, 800], "k--")
axs[1, 0].plot([0, 41], [800, 800], "k--")

# Plot zoomed-in liquid concentration data
for i, spc in enumerate(spcnames):
    label = spc
    axs[2, 0].plot(range(1, 41), initial_conditions["liquid_concentration"][spc], lines[i], label=label)

axs[2, 0].set_ylabel("Conc. (mol/m^3)", size=20)
axs[2, 0].set_ylim([-20, 600])
axs[2, 0].set_title("Zoomed in on (b)", loc="left", size=20)
axs[2, 0].set_xticklabels([])

# Plot vapor concentration data
for i, spc in enumerate(spcnames):
    label = spc
    axs[1, 1].plot(range(1, 41), initial_conditions["vapor_concentration"][spc], lines[i], label=label)

axs[1, 1].set_title("(c) Vapor phase", loc="left", size=20)
axs[1, 1].set_xticklabels([])
axs[1, 1].plot([0, 0], [-10, 20], "k--")
axs[1, 1].plot([0, 41], [-10, -10], "k--")
axs[1, 1].plot([41, 41], [-10, 20], "k--")
axs[1, 1].plot([0, 41], [20, 20], "k--")

handles, labels = axs[1, 1].get_legend_handles_labels()

# Hide the second subplot in the first row
axs[0, 1].legend(handles=handles, loc="center", fontsize=14)
axs[0, 1].axis("off")

# Plot zoomed-in vapor concentration data
for i, spc in enumerate(spcnames):
    label = spc
    axs[2, 1].plot(range(1, 41), initial_conditions["vapor_concentration"][spc], lines[i], label=label)

axs[2, 1].set_ylim([-1, 20])
axs[2, 1].set_title("Zoomed in on (c)", loc="left", size=20)
axs[2, 1].set_xticklabels([])

# Plot liquid phase residence time
flowrate = np.array(initial_conditions["liquid_outlet_volumetric_flowrate"])
flowrate[-1] = flowrate[-2]
residence_times = Vliq / flowrate
axs[3, 0].plot(range(1, 41), residence_times, color="k", marker="o")
axs[3, 0].set_ylabel("Residence time (s)", size=20)
axs[3, 0].set_title("(d) Liquid phase", loc="left", size=20)
axs[3, 0].set_ylim([0, np.max(residence_times) * 1.1])
if model_name == "basecase_debutanizer_model":
    axs[3, 0].set_xlabel("Tray", size=20)
elif model_name == "trace_oxygen_perturbed_debutanizer_model":
    axs[3, 0].set_xticklabels([])

# Plot vapor phase residence time
flowrate = np.array(initial_conditions["vapor_outlet_volumetric_flowrate"])
flowrate[0] = flowrate[1]
residence_times = Vvap / flowrate
axs[3, 1].plot(range(1, 41), residence_times, color="k", marker="o")
axs[3, 1].set_title("(e) Vapor phase", loc="left", size=20)
axs[3, 1].set_ylim([0, np.max(residence_times) * 1.1])
if model_name == "basecase_debutanizer_model":
    axs[3, 1].set_xlabel("Tray", size=20)
elif model_name == "trace_oxygen_perturbed_debutanizer_model":
    axs[3, 1].set_xticklabels([])

# Plot additional plots for "trace_oxygen_perturbed_debutanizer_model"
if model_name == "trace_oxygen_perturbed_debutanizer_model":
    spc = "OXYGEN"
    i = spcnames.index(spc)
    label = spc
    axs[4, 0].plot(range(1, 41), initial_conditions["liquid_concentration"][spc], lines[i], label=spc, color="C8")
    axs[4, 0].set_ylabel("Conc. (mol/m^3)", size=20)
    axs[4, 0].set_title("(f) Liquid phase", loc="left", size=20)
    axs[4, 0].set_xlabel("Tray", size=20)
    axs[4, 0].set_yscale("log")
    axs[4, 0].set_ylim([1e-18, 1e0])

    axs[4, 1].plot(range(1, 41), initial_conditions["vapor_concentration"][spc], lines[i], label=spc, color="C8")
    axs[4, 1].set_title("(g) Vapor phase", loc="left", size=20)
    axs[4, 1].set_xlabel("Tray", size=20)
    axs[4, 1].set_yscale("log")
    axs[4, 0].set_ylim([1e-18, 1e0])

fig.tight_layout()

# Save the figure
if model_name == "basecase_debutanizer_model":
    save_name = "monomer_conc_trays"
elif model_name == "trace_oxygen_perturbed_debutanizer_model":
    save_name = "monomer_conc_trays_oxygen"
plt.savefig(f"Figures/{save_name}.pdf", bbox_inches="tight")
