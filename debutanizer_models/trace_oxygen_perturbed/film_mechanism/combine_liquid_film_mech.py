import yaml

liquid_mech_path = "debutanizer_models/trace_oxygen_perturbed/liquid_mechanism/chem.rms"
film_mech_path = "debutanizer_models/trace_oxygen_perturbed/film_mechanism/chem_film_phase.rms"

with open(liquid_mech_path, "r") as f:
    liquid_mech_dict = yaml.load(stream=f, Loader=yaml.Loader)

with open(film_mech_path, "r") as f:
    film_mech_dict = yaml.load(stream=f, Loader=yaml.Loader)

liquid_liq_spc_names = set([spc["name"] for spc in liquid_mech_dict["Phases"][0]["Species"]])
film_liq_spc_names = set([spc["name"] for spc in film_mech_dict["Phases"][1]["Species"]])
extra_liq_spc_names = liquid_liq_spc_names - film_liq_spc_names

for spc_dict in liquid_mech_dict["Phases"][0]["Species"]:
    if spc_dict["name"] in extra_liq_spc_names:
        film_mech_dict["Phases"][1]["Species"].append(spc_dict)

film_mech_dict["Reactions"] += liquid_mech_dict["Reactions"]

with open("chem_liquid_film_phase.rms","w+") as f:
    yaml.dump(film_mech_dict,stream=f)

print("Done!")