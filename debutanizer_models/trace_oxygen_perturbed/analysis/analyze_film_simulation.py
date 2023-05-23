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
        "--simulation_results_path", type=str, required=True, help="The path to the csv file containing the film simulation results.",
    )
    parser.add_argument(
        "--rop_results_path", type=str, required=True, help="The path to the csv file containing the film rop results.",
    )

    args = parser.parse_args()
    model_name = args.model_name
    simulation_results_path = args.simulation_results_path
    rop_results_path = args.rop_results_path

    return (
        model_name,
        simulation_results_path,
        rop_results_path,
    )

(
    model_name,
    simulation_results_path,
    rop_results_path,
) = parse_arguments()

