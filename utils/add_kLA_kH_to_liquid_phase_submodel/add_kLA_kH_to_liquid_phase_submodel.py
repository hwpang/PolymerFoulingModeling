"""
This script can be used to add/update the kLA and kH parameters to a given liquid phase mechanism in .rms format.
"""
import os
import yaml
import argparse

from rmgpy import settings
from rmgpy.chemkin import load_chemkin_file
from rmgpy.quantity import Quantity
from rmgpy.rmg.main import RMG
from rmgpy.kinetics.diffusionLimited import diffusion_limiter
from rmgpy.data.vaporLiquidMassTransfer import vapor_liquid_mass_transfer, liquidVolumetricMassTransferCoefficientPowerLaw

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--chemkin_path",
        type=str,
        required=True,
        help="The path to liquid phase mechanism chemkin file",
    )
    parser.add_argument(
        "--species_dict_path",
        type=str,
        required=True,
        help="The path to species dict file",
    )
    parser.add_argument(
        "--rms_path", type=str, required=True, help="The path to rms file",
    )
    parser.add_argument(
        "--save_directory",
        type=str,
        default=".",
        help="The directory to save the generated film phase mechanism.",
    )

    args = parser.parse_args()

    chemkin_path = args.chemkin_path
    species_dict_path = args.species_dict_path
    rms_path = args.rms_path
    save_directory = args.save_directory

    return (
        chemkin_path,
        species_dict_path,
        rms_path,
        save_directory,
    )


# +

(
    chemkin_path,
    species_dict_path,
    rms_path,
    save_directory,
) = parse_arguments()

liqspcs, liqrxns = load_chemkin_file(chemkin_path,
                        dictionary_path=species_dict_path)

rmg = RMG()

solvent = "benzene"
rmg.database_directory = settings["database.directory"]
rmg.thermo_libraries = ['primaryThermoLibrary']
rmg.kinetics_families = 'default'
rmg.kinetics_depositories = ['training']
rmg.kinetics_estimator = 'rate rules'
rmg.solvent = solvent
rmg.reaction_libraries = []

rmg.load_database()

solvent_data = rmg.database.solvation.get_solvent_data(rmg.solvent)
diffusion_limiter.enable(solvent_data, rmg.database.solvation)

with open(rms_path,"r") as f:
    mech_dict = yaml.load(stream=f,Loader=yaml.Loader)

prefactor = 14.5
diffusion_coefficient_power = 1/2
solvent_viscosity_power = -1/6
solvent_density_power = -1/6

liquid_volumetric_mass_transfer_coefficient_power_law = liquidVolumetricMassTransferCoefficientPowerLaw(
    prefactor=Quantity(prefactor).value_si,
    diffusion_coefficient_power=diffusion_coefficient_power,
    solvent_viscosity_power=solvent_viscosity_power,
    solvent_density_power=solvent_density_power,
)

vapor_liquid_mass_transfer.enable(solvent_data, rmg.database.solvation, liquid_volumetric_mass_transfer_coefficient_power_law)

label_to_species = {}
for spc in liqspcs:
    label_to_species[spc.label] = spc

for spc_dict in mech_dict["Phases"][0]["Species"]:
    try:
        spc = label_to_species[spc_dict["name"]]
    except:
        print(spc["name"])
    else:
        spc.get_liquid_volumetric_mass_transfer_coefficient_data()
        Ts = spc.liquid_volumetric_mass_transfer_coefficient_data.Ts
        kLAs = spc.liquid_volumetric_mass_transfer_coefficient_data.kLAs
        spc_dict["liquidvolumetricmasstransfercoefficient"] = dict()
        spc_dict["liquidvolumetricmasstransfercoefficient"]["type"] = "TemperatureDependentLiquidVolumetricMassTransferCoefficient"
        spc_dict["liquidvolumetricmasstransfercoefficient"]["Ts"] = Ts
        spc_dict["liquidvolumetricmasstransfercoefficient"]["kLAs"] = kLAs
        spc.get_henry_law_constant_data()
        Ts = spc.henry_law_constant_data.Ts
        kHs = spc.henry_law_constant_data.kHs
        spc_dict["henrylawconstant"] = dict()
        spc_dict["henrylawconstant"]["type"] = "TemperatureDependentHenryLawConstant"
        spc_dict["henrylawconstant"]["Ts"] = Ts
        spc_dict["henrylawconstant"]["kHs"] = kHs

with open(os.path.join(save_directory, "chem.rms"),"w+") as f:
    yaml.dump(mech_dict,stream=f)