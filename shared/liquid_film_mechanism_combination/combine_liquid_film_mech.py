import os
import yaml
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--liquid_rms_path",
        type=str,
        required=True,
        help="The path to liquid phase mechanism rms file",
    )
    parser.add_argument(
        "--film_rms_path",
        type=str,
        required=True,
        help="The path to film phase mechanism rms file",
    )
    parser.add_argument(
        "--save_directory",
        type=str,
        default="film_mechanism",
        help="The directory to save the generated film phase mechanism.",
    )

    args = parser.parse_args()

    liquid_rms_path = args.liquid_rms_path
    film_rms_path = args.film_rms_path
    save_directory = args.save_directory

    return liquid_rms_path, film_rms_path, save_directory


liquid_rms_path, film_rms_path, save_directory = parse_arguments()

os.makedirs(save_directory, exist_ok=True)

with open(liquid_rms_path, "r") as f:
    liquid_mech_dict = yaml.load(stream=f, Loader=yaml.Loader)

with open(film_rms_path, "r") as f:
    film_mech_dict = yaml.load(stream=f, Loader=yaml.Loader)

for spc_dict in liquid_mech_dict["Phases"][0]["Species"]:
    spc_dict["name"] = spc_dict["name"] + "(L)"

liquid_liq_spc_names = set(
    [spc_dict["name"] for spc_dict in liquid_mech_dict["Phases"][0]["Species"]]
)
film_liq_spc_names = set(
    [spc_dict["name"] for spc_dict in film_mech_dict["Phases"][1]["Species"]]
)
extra_liq_spc_names = liquid_liq_spc_names - film_liq_spc_names

for spc_dict in liquid_mech_dict["Phases"][0]["Species"]:
    if spc_dict["name"] in extra_liq_spc_names:
        film_mech_dict["Phases"][1]["Species"].append(spc_dict)

for rxn_dict in liquid_mech_dict["Reactions"]:
    for i, spc_name in enumerate(rxn_dict["reactants"]):
        rxn_dict["reactants"][i] = spc_name + "(L)"
    for i, spc_name in enumerate(rxn_dict["products"]):
        rxn_dict["products"][i] = spc_name + "(L)"
    if "efficiencies" in rxn_dict["kinetics"]:
        for spc_name in list(rxn_dict["kinetics"]["efficiencies"]):
            if "(L)" not in spc_name:
                rxn_dict["kinetics"]["efficiencies"][spc_name + "(L)"] = rxn_dict[
                    "kinetics"
                ]["efficiencies"][spc_name]
                del rxn_dict["kinetics"]["efficiencies"][spc_name]

film_mech_dict["Reactions"] += liquid_mech_dict["Reactions"]

path = os.path.join(save_directory, "chem_liquid_film_phase.rms")
with open(path, "w+") as f:
    yaml.dump(film_mech_dict, stream=f)

print("Done!")
