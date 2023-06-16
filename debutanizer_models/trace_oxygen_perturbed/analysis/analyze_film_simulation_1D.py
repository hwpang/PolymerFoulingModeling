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

print("Load 1-D film simulation results")

one_d_simulations = dict()
for tray in trays:
    one_d_simulation_path = os.path.join(simulation_directory, f"simulation_film_1D_{tray}.csv")
    one_d_simulations[tray] = pd.read_csv(one_d_simulation_path)

print("Plot 1-D film simulation results")
nrows = len(selected_trays)
ncols = 1

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9, 12), sharex=True)

for ind, tray in enumerate(selected_trays):
    axs[ind].plot()