import os
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

# plot QCM-D simulation

expt_rate_dict = {}
expt_rate_dict["Unsparged"] = 9.4  # ng/hr
expt_rate_dict["N2 sparged"] = 6.1
expt_rate_dict["O2 sparged"] = 927.3

perturb_species = "O2"
perturb_factor_list = ["0.0", "1e-3", "1e-2", "1e-1", "1e0"]
factor_num_list = [float(perturb_factor) for perturb_factor in perturb_factor_list]
tray = 1

print("Loading film simulation results...")
film_simulations = dict()
for perturb_factor in perturb_factor_list:
    film_simulation_path = os.path.join(all_simulation_directory, f"{perturb_species}_{perturb_factor}", f"simulation_film_{tray}.csv")
    film_simulations[perturb_factor] = pd.read_csv(film_simulation_path)

# pure_O2_sat_conc = 0.008093871600706785 #M estimated by RMG
pure_O2_sat_molfrac = 8.1 * 1e-4  # mol fraction measured experimentally from Battino, R., Rettich, T. R., & Tominaga, T. (1983). The Solubility of Oxygen and Ozone in Liquids. Journal of Physical and Chemical Reference Data, Vol. 12, pp. 163–178. https://doi.org/10.1063/1.555680
benzene_density = 876  # kg/m^3 from Lide, D. R., Data, S. R., Board, E. A., Baysinger, G., Chemistry, S., Library, C. E., … Zwillinger, D. (2004). CRC Handbook of Chemistry and Physics. 2660.
benzene_mw = 78.11 / 1000  # kg/mol from Lide, D. R., Data, S. R., Board, E. A., Baysinger, G., Chemistry, S., Library, C. E., … Zwillinger, D. (2004). CRC Handbook of Chemistry and Physics. 2660.
benzene_sat_conc = benzene_density / benzene_mw  # mol/m^3
pure_O2_sat_conc = pure_O2_sat_molfrac * benzene_sat_conc  # mol/m^3
air_O2_sat_conc = 0.21 * pure_O2_sat_conc  # mol/m^3
air_O2_sat_conc /= 1000  # mol/L

def calculate_film_growth_rate(df):
    inds = range(1, len(df.index))
    dmdts = np.array([(df.loc[ind, "mass"] - df.loc[0, "mass"]) / (df.loc[ind, "timestamp"] - df.loc[0, "timestamp"])  for ind in inds])
    return np.mean(dmdts), np.std(dmdts)


plt.figure(figsize=(4, 3))
xs = np.array(factor_num_list) * air_O2_sat_conc

# simulated film growth rate
film_growth_rates = [calculate_film_growth_rate(film_simulations[perturb_factor]) for perturb_factor in perturb_factor_list]
print("film_growth_rates: ", film_growth_rates)
ys = np.array([film_growth_rate[0] for film_growth_rate in film_growth_rates])
ys += ys / rho * epsilon * rho_liq  # converting from mass of solid to mass of film by adding mass of liquid in film
ys *= 1e9 * 1e3 * 3600  # kg/s to ng/hr
yerr = 10.0  # a perturb_factor of x error
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
plt.errorbar([expt_o2], [expt_rate_mean], yerr=[yerr], xerr=[[expt_o2 - expt_o2 / xerr], [-expt_o2 + expt_o2 * xerr]], color=color, marker="o", capsize=5, label=label)
plt.fill_between([expt_o2 / xerr, expt_o2 * xerr], [expt_rate_mean - yerr, expt_rate_mean - yerr], [expt_rate_mean + yerr, expt_rate_mean + yerr], color=color, alpha=0.5, linewidth=0)

# # O2 sparged
# label = "O2 sparged"
# expt_rate_mean = expt_rate_dict[label]  # ng/hr
# xerr = 1
# reduction_factor = 10
# expt_o2 = air_O2_sat_conc  # M
# yerr = expt_rate_std
# color = "C3"
# plt.errorbar([expt_o2], [expt_rate_mean], yerr=[yerr], xerr=[[expt_o2 - expt_o2 / xerr], [-expt_o2 + expt_o2 * xerr]], color=color, marker="o", capsize=5, label=label)
# plt.fill_between([expt_o2 / xerr, expt_o2 * xerr], [expt_rate_mean - yerr, expt_rate_mean - yerr], [expt_rate_mean + yerr, expt_rate_mean + yerr], color=color, alpha=0.5, linewidth=0)

plt.xlabel("O2 (M)")
plt.ylabel("Film growth rate (ng/hr)")
plt.xscale("symlog", linthresh=1e-6)
plt.yscale("log")
# plt.ylim([1e-3, 1e2])
plt.tight_layout()
plt.legend(fontsize=9)

os.makedirs("Figures", exist_ok=True)

plt.savefig("Figures/QCMD_cell_model_film_growth_rates.pdf", bbox_inches="tight")
plt.close()
