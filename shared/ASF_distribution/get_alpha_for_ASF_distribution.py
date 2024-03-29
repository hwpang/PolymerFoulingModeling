# %%
import sys

sys.path.insert(0, "/home/gridsan/hwpang/Software/RMG-Py/")
import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
import re
import yaml
import argparse
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from joblib import Parallel, delayed

from rmgpy import settings
from rmgpy.rmg.main import RMG
from rmgpy.molecule.group import Group
from rmgpy.chemkin import load_chemkin_file
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.kinetics.kineticsdata import KineticsData
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.depository import DepositoryReaction
from rmgpy.kinetics.diffusionLimited import diffusion_limiter


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--chemkin_path", type=str, required=True, help="The path including chemkin file",
    )
    parser.add_argument(
        "--species_dict_path",
        type=str,
        required=True,
        help="The path including species dict file",
    )
    parser.add_argument(
        "--rms_path", type=str, required=True, help="The path including rms file."
    )
    parser.add_argument(
        "--aspen_condition_path",
        type=str,
        default=None,
        help="The path including aspen condition file.",
    )
    parser.add_argument(
        "--results_directory",
        type=str,
        required=True,
        help="The path including rop csvs from simulations.",
    )
    parser.add_argument(
        "--n_jobs", type=int, default=1, help="The number of jobs to run in parallel.",
    )
    parser.add_argument(
        "--model_name",
        type=str,
        required=True,
        help="The name of the model, i.e., basecase_debutanizer_model or QCMD_cell_model.",
    )

    args = parser.parse_args()

    chemkin_path = args.chemkin_path
    species_dict_path = args.species_dict_path
    rms_path = args.rms_path
    aspen_condition_path = args.aspen_condition_path
    results_directory = args.results_directory
    n_jobs = args.n_jobs
    model_name = args.model_name

    return (
        chemkin_path,
        species_dict_path,
        rms_path,
        aspen_condition_path,
        results_directory,
        n_jobs,
        model_name,
    )


# Create a helper function to parallelize the calculation for trimerization rates
def generate_trimerization(spc_i, spc_j):
    # Generate trimerization DA reactions
    reaction_list = rmg.database.kinetics.generate_reactions_from_families(
        reactants=[spc_i, spc_j], only_families=["Diels_alder_addition"]
    )
    # Assign kinetics and thermochemistries
    cerm1 = CoreEdgeReactionModel()
    cerm1.solvent_name = "benzene"
    for rxn0 in reaction_list:
        for spc in rxn0.reactants + rxn0.products:
            cerm1.generate_thermo(spc)
        cerm1.apply_kinetics_to_reaction(rxn0)
        if isinstance(rxn0.kinetics, KineticsData):
            rxn0.kinetics = rxn0.kinetics.to_arrhenius()
        if isinstance(rxn0, (TemplateReaction, DepositoryReaction)):
            rxn0.fix_barrier_height()
    return reaction_list


(
    chemkin_path,
    species_dict_path,
    rms_path,
    aspen_condition_path,
    results_directory,
    n_jobs,
    model_name,
) = parse_arguments()

# %%
### You can manually set the path

# chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/liquid_mechanism/chem_annotated.inp"
# species_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/liquid_mechanism/species_dictionary.txt"
# rms_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/liquid_mechanism/chem.rms"
# aspen_condition_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/aspen_simulation/aspen_conditions.yml"
# results_directory = "../../simulation_results/1,3-BUTADIENE_1.0_3600.0_64.0/"
# n_jobs = 4
# model_name = "basecase_debutanizer_model"

if model_name in ["basecase_debutanizer_model", "trace_oxygen_perturbed_debutanizer_model"]:
    initial_monomer_labels = [
        "N-BUTANE",
        "2-BUTENE",
        "1,3-BUTADIENE",
        "CYCLOPENTADIENE",
        "BENZENE",
        "1,3-CYCLOHEXADIENE",
        "TOLUENE",
        "STYRENE",
    ]
    d = 2.5
    h = 0.3
    A = (d / 2) ** 2 * np.pi
    Vliq = A * h
    with open(aspen_condition_path, "r") as f:
        aspen_condition = yaml.load(f, Loader=yaml.FullLoader)
    temp = aspen_condition["T"]

    trays = np.arange(0, 40, 1)
    ROP_CSV = r"^simulation_vapor_liquid_liqrop_(\d+).csv"
    
    if model_name == "basecase_debutanizer_model":
        time_step = 64.0
    elif model_name == "trace_oxygen_perturbed_debutanizer_model":
        time_step = 32.0

    tf = np.arange(0.0, 3600.0+time_step, time_step)[-1]
    SS_MOL_CSV = f"simulation_vapor_liquid_yvapn_{tf}.csv"
elif model_name == "QCMD_cell_model":
    initial_monomer_labels = [
        "5-methylcyclohexadiene",
        "1-methylcyclohexadiene",
        "2-methylcyclohexadiene",
        "methylenecyclohexene",
    ]
    Vliq = 40 * 1e-9
    temp = [90.0 + 273.15]

    trays = np.arange(0, 1)
    ROP_CSV = r"^simulation_liquid_liqrop_(\d+).csv"
    SS_MOL_CSV = "simulation_liquid_1.csv"

RXN_TYPES = [
    "R._Add",
    "R.+O2",
    "R._Habs",
    "R._Recomb",
    "R._Disprop",
    "R._CycEther",
    "ROO._Add",
    "ROO._Habs_ROO.",
    "ROO._Habs",
    "ROO._Recomb",
    "ROO._Disprop",
    "ROO._R._Habs",
    "ROO._eli",
    "RO._Add",
    "RO._Habs",
    "RO._Recomb",
    "RO._Disprop",
    "RO._CycEther",
    "diels_alder",
    "retroene",
]

ROO = Group().from_adjacency_list(
    """
1  O u1 p2 c0 {2,S}
2  O u0 p2 c0 {1,S} {3,S}
3  C u0 p0 c0 {2,S}
"""
)

RO = Group().from_adjacency_list(
    """
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S}
"""
)


def Wn(n, alpha):
    return n * ((1 - alpha) ** 2) * (alpha ** (n - 1))


# %% [markdown]
# ## Load RMS model file

# %%
print("loading RMS model...")
with open(rms_path) as f:
    phase_dict = yaml.load(f, Loader=yaml.Loader)
rms_spcs = phase_dict["Phases"][0]["Species"]
rms_rxns = phase_dict["Reactions"]

# %% [markdown]
# ## Load RMG database and model file

# %%
print("loading RMG model...")
rmg = RMG()
rmg.database_directory = settings["database.directory"]
rmg.thermo_libraries = ["primaryThermoLibrary"]
rmg.kinetics_families = "default"
rmg.kinetics_depositories = ["training"]
rmg.kinetics_estimator = "rate rules"
rmg.solvent = "benzene"
rmg.reaction_libraries = []

rmg.load_database()
# solvent_data = rmg.database.solvation.get_solvent_data(rmg.solvent)
# diffusion_limiter.enable(solvent_data, rmg.database.solvation)

liqspcs, liqrxns = load_chemkin_file(chemkin_path, species_dict_path)

# %% [markdown]
# ## Mark reaction families

def get_reaction_family(rxn):
    if rxn.family in rmg.database.kinetics.families:
        return rxn.family
    else:
        rxns = rmg.database.kinetics.generate_reactions_from_families(
            reactants=rxn.reactants, products=rxn.products
        )
        if rxns:
            print(rxns[0].family)
            return rxns[0].family
        else:
            return rxn.family

families = Parallel(n_jobs=n_jobs, verbose=5, backend="multiprocessing")(
    delayed(get_reaction_family)(rxn) for rxn in liqrxns
)

for rxn, family in zip(liqrxns, families):
    rxn.family = family

# %%
# Check if the RMG and RMS model match each other
if len(rms_spcs) != len(liqspcs):
    raise ValueError("The RMG and RMS model have different number of species.")
if len(rms_rxns) != len(liqrxns):
    raise ValueError("The RMG and RMS model have different number of reactions.")

# %% [markdown]
# ## Sort Species according to RMS indexes

# %%
print("Updating RMG species and reaction list...")
rms_idx_to_spcs = {}
label_spc_dict = {}
for spc in liqspcs:
    label_spc_dict[spc.label] = spc

for idx, rms_spc in enumerate(rms_spcs):
    rms_label = rms_spc["name"]
    spc = label_spc_dict.get(rms_label)
    if spc:
        rms_idx_to_spcs[idx] = spc
    else:
        raise RuntimeError(
            f"Unexpected Error: RMS model has a species {rms_label} not inside the RMG model. Contact Xiaorui!"
        )

liqspcs = [item[1] for item in sorted(rms_idx_to_spcs.items(), key=lambda x: x[0])]
spc_label_index_dict = {}
for ind, spc in enumerate(liqspcs):
    spc_label_index_dict[spc.label] = ind
source_label_index_dict = {}
for ind, rxn in enumerate(liqrxns):
    label = ""
    for spc in rxn.reactants:
        label += spc.label
        label += "+"
    label = label[:-1]
    label += "<=>"
    for spc in rxn.products:
        label += spc.label
        label += "+"
    label = label[:-1]
    source_label_index_dict[label] = ind
source_label_index_dict["cond"] = len(liqrxns)
source_label_index_dict["evap"] = len(liqrxns) + 1
source_label_index_dict["inlet"] = len(liqrxns) + 2
source_label_index_dict["outlet"] = len(liqrxns) + 3
# %% [markdown]
# ## Sort Reactions according to RMS indexes

# %%
rms_idx_to_rxns = {}
rxns_not_match = []
for idx, rms_rxn in enumerate(rms_rxns):
    rmg_rxn = liqrxns[idx]
    rms_kin = rms_rxn["kinetics"]
    rmg_kin = rmg_rxn.kinetics
    match = True
    # Match reactants and products
    for r in rmg_rxn.reactants:
        if r.label not in rms_rxn["reactants"]:
            match = False
            break
    for p in rmg_rxn.products:
        if p.label not in rms_rxn["products"]:
            match = False
            break
    if rms_kin["type"] != rmg_kin.__repr__().split("(")[0]:
        match = False
    if not match:
        rxns_not_match.append((idx, rms_rxn))
        continue
    rms_idx_to_rxns[idx] = rmg_rxn

if rxns_not_match:
    raise RuntimeError(
        f"Unexpected Error: RMS model reaction list is different than RMG model. Contact Xiaorui!"
    )

liqrxns = [item[1] for item in sorted(rms_idx_to_rxns.items(), key=lambda x: x[0])]

# %% [markdown]
# ## Identify R and ROO

# %%
print("Categorizing reactions...")
Rs_id = []
ROOs_id = []
ROs_id = []

# Find Rs and ROOs
for i, spc in enumerate(liqspcs):
    if not spc.molecule[0].is_radical():
        continue

    if spc.molecule[0].is_subgraph_isomorphic(ROO):
        ROOs_id.append(i)
    elif spc.molecule[0].is_subgraph_isomorphic(RO):
        ROs_id.append(i)
    elif "C" in spc.smiles or "c" in spc.smiles:
        Rs_id.append(i)

# %%
def find_R_ROO_and_RO(spc_list):
    R, ROO, RO = False, False, False
    for spc in spc_list:
        spc_id = spc_label_index_dict.get(spc.label)
        if spc_id in Rs_id:
            R = True
        elif spc_id in ROOs_id:
            ROO = True
        elif spc_id in ROs_id:
            RO = True
    return R, ROO, RO


# %%
################# Classify reactions #########################
family_categories = {rxn_type: [] for rxn_type in RXN_TYPES}

for i, rxn in enumerate(liqrxns):
    rxn_family = rxn.family

    if rxn_family == "Retroene":
        family_categories["retroene"].append(i)
        continue

    elif rxn_family == "Diels_alder_addition":
        family_categories["diels_alder"].append(i)
        continue

    # Find if R or ROO in the reaction
    r_in, roo_in, ro_in = find_R_ROO_and_RO(rxn.reactants + rxn.products)
    if not r_in and not roo_in and not ro_in:
        continue  # irrelavant reaction

    if rxn_family == "R_Recombination":
        # if O2 and R + O2
        has_O2 = any([spc.smiles == "[O][O]" for spc in rxn.reactants + rxn.products])
        if has_O2:
            if r_in and roo_in:
                family_categories["R.+O2"].append(i)
            else:
                pass
        else:
            if len(rxn.reactants) == 2:
                r_in, roo_in, ro_in = find_R_ROO_and_RO(rxn.reactants)
            else:
                r_in, roo_in, ro_in = find_R_ROO_and_RO(rxn.products)
            if r_in:
                family_categories["R._Recomb"].append(i)
            if roo_in:
                family_categories["ROO._Recomb"].append(i)
            if ro_in:
                family_categories["RO._Recomb"].append(i)

    elif rxn_family == "Disproportionation":
        has_O2 = any([spc.smiles == "[O][O]" for spc in rxn.reactants + rxn.products])
        if not has_O2:
            if r_in:
                family_categories["R._Disprop"].append(i)
            if roo_in:
                family_categories["ROO._Disprop"].append(i)
            if ro_in:
                family_categories["RO._Disprop"].append(i)

    elif rxn_family == "H_Abstraction":
        if r_in:
            family_categories["R._Habs"].append(i)
        elif roo_in:
            family_categories["ROO._Habs_ROO."].append(i)
        if roo_in:
            family_categories["ROO._Habs"].append(i)
        if ro_in:
            family_categories["RO._Habs"].append(i)

    elif rxn_family == "R_Addition_MultipleBond":
        if len(rxn.reactants) == 2:
            r_in, roo_in, ro_in = find_R_ROO_and_RO(rxn.reactants)
        else:
            r_in, roo_in, ro_in = find_R_ROO_and_RO(rxn.products)
        if r_in:
            family_categories["R._Add"].append(i)
        if roo_in:
            family_categories["ROO._Add"].append(i)
        if ro_in:
            family_categories["RO._Add"].append(i)

    elif rxn_family == "HO2_Elimination_from_PeroxyRadical":
        family_categories["ROO._eli"].append(i)

    elif rxn_family == "Cyclic_Ether_Formation":
        family_categories["R._CycEther"].append(i)
        family_categories["RO._CycEther"].append(i)

# %%
rop_files = {}
for file in os.listdir(results_directory):
    if re.match(ROP_CSV, file):
        rop_files[int(re.split(ROP_CSV, file)[1]) - 1] = file

# %%
print("loading all ROP results...")
rop_matrix = {}
conc = {}

if model_name in ["basecase_debutanizer_model", "trace_oxygen_perturbed_debutanizer_model"]:
    ss_mol_df = pd.read_csv(os.path.join(results_directory, SS_MOL_CSV))
elif model_name == "QCMD_cell_model":
    ss_mol_df = pd.read_csv(os.path.join(results_directory, SS_MOL_CSV))
    ss_mol_df = ss_mol_df.iloc[[-1], :]
    ss_mol_df.reset_index(inplace=True)

for tray in trays:
    try:
        rop_df = pd.read_csv(os.path.join(results_directory, rop_files[tray]))
    except KeyError:
        print(results_directory)
        raise RuntimeError(f"Missing rop csv or outlet csv at tray {tray}.")
    spcinds = np.array([spc_label_index_dict[spcname] for spcname in rop_df["rop_spcname"]])
    sourceinds = np.array(
        [source_label_index_dict[spcname] for spcname in rop_df["rop_sourcestring"]]
    )
    if np.max(sourceinds) > len(liqrxns) - 1 + 4:
        raise ValueError(
            f"The ROP result {rop_files[tray]} contains reaction index exceeding the number of reactions in the model."
        )
    if np.max(spcinds) > len(liqspcs):
        raise ValueError(
            f"The ROP result {rop_files[tray]} contains species index exceeding the number of species in the model."
        )
    sparse_mat = coo_matrix(
        (rop_df["rop"].to_numpy(), (spcinds, sourceinds)),
        shape=(len(liqspcs), len(liqrxns) + 4),
    )
    try:
        rop_matrix[tray] = sparse_mat.toarray()
    except AttributeError:
        rop_matrix[tray] = sparse_mat.to_array()

    conc[tray] = np.zeros(len(liqspcs))
    for label, idx in spc_label_index_dict.items():
        conc[tray][idx] = ss_mol_df.loc[tray, label] / Vliq

Rs_labels = [liqspcs[i].label for i in Rs_id]
ROOs_labels = [liqspcs[i].label for i in ROOs_id]
ROs_labels = [liqspcs[i].label for i in ROs_id]

# %%
print("Getting detailed rates...")
rxn_rates = {}
for family in family_categories:
    rxn_rates[family] = {}
    for i in family_categories[family]:
        rxn = liqrxns[i]
        rxn_label = "+".join([spc.label for spc in rxn.reactants]) + "<=>" + "+".join(spc.label for spc in rxn.products)
        rxn_rates[family][rxn_label] = np.zeros_like(trays, dtype=float)
        for tray in trays:
            rops = rop_matrix[tray][:, [i]][Rs_id, :]
            rops = rops[rops < 0]
            rxn_rates[family][rxn_label][tray] = np.sum(rops)
        rxn_rates[family][rxn_label] = rxn_rates[family][rxn_label].tolist()

# %%
print("Calculating alpha...")
rates = {}
production_rates = {}
Habs = np.zeros_like(trays, dtype=float)

for family in family_categories.keys():
    rates[family] = np.zeros_like(trays, dtype=float)

rates["R._outlet"] = np.zeros_like(trays, dtype=float)
rates["ROO._outlet"] = np.zeros_like(trays, dtype=float)
rates["RO._outlet"] = np.zeros_like(trays, dtype=float)
rates["R._evap"] = np.zeros_like(trays, dtype=float)
rates["ROO._evap"] = np.zeros_like(trays, dtype=float)
rates["RO._evap"] = np.zeros_like(trays, dtype=float)

production_rates["R._inlet"] = np.zeros_like(trays, dtype=float)
production_rates["ROO._inlet"] = np.zeros_like(trays, dtype=float)
production_rates["RO._inlet"] = np.zeros_like(trays, dtype=float)
production_rates["R._cond"] = np.zeros_like(trays, dtype=float)
production_rates["ROO._cond"] = np.zeros_like(trays, dtype=float)
production_rates["RO._cond"] = np.zeros_like(trays, dtype=float)
production_rates["R._RevDisprop"] = np.zeros_like(trays, dtype=float)
production_rates["RO._BondDiss"] = np.zeros_like(trays, dtype=float)

for tray in trays:
    for family, rxn_indices in family_categories.items():
        if "R." in family:
            rops = rop_matrix[tray][:, rxn_indices][Rs_id, :]
            rops = rops[rops < 0]
            rates[family][tray] = np.sum(rops)
        elif "ROO." in family:
            rops = rop_matrix[tray][:, rxn_indices][ROOs_id, :]
            rops = rops[rops < 0]
            rates[family][tray] = np.sum(rops)
        elif "RO." in family:
            rops = rop_matrix[tray][:, rxn_indices][ROs_id, :]
            rops = rops[rops < 0]
            rates[family][tray] = np.sum(rops)

        if family == "R._Disprop":
            rops = rop_matrix[tray][:, rxn_indices][Rs_id, :]
            rops = rops[rops > 0]
            production_rates["R._RevDisprop"][tray] = np.sum(rops)
        elif family == "RO._Recomb":
            rops = rop_matrix[tray][:, rxn_indices][ROs_id, :]
            rops = rops[rops > 0]
            production_rates["RO._BondDiss"][tray] = np.sum(rops)

    outlet_ind = source_label_index_dict["outlet"]
    rops = rop_matrix[tray][Rs_id, outlet_ind]
    rops = rops[rops < 0]
    rates["R._outlet"][tray] = np.sum(rops)
    rops = rop_matrix[tray][ROOs_id, outlet_ind]
    rops = rops[rops < 0]
    rates["ROO._outlet"][tray] = np.sum(rops)
    rops = rop_matrix[tray][ROs_id, outlet_ind]
    rops = rops[rops < 0]
    rates["RO._outlet"][tray] = np.sum(rops)

    evap_ind = source_label_index_dict["evap"]
    rops = rop_matrix[tray][Rs_id, evap_ind]
    rops = rops[rops < 0]
    rates["R._evap"][tray] = np.sum(rops)
    rops = rop_matrix[tray][ROOs_id, evap_ind]
    rops = rops[rops < 0]
    rates["ROO._evap"][tray] = np.sum(rops)
    rops = rop_matrix[tray][ROs_id, evap_ind]
    rops = rops[rops < 0]
    rates["RO._evap"][tray] = np.sum(rops)

    inlet_ind = source_label_index_dict["inlet"]
    rops = rop_matrix[tray][Rs_id, inlet_ind]
    rops = rops[rops > 0]
    production_rates["R._inlet"][tray] = np.sum(rops)
    rops = rop_matrix[tray][ROOs_id, inlet_ind]
    rops = rops[rops > 0]
    production_rates["ROO._inlet"][tray] = np.sum(rops)
    rops = rop_matrix[tray][ROs_id, inlet_ind]
    rops = rops[rops > 0]
    production_rates["RO._inlet"][tray] = np.sum(rops)

    cond_ind = source_label_index_dict["cond"]
    rops = rop_matrix[tray][Rs_id, cond_ind]
    rops = rops[rops > 0]
    production_rates["R._cond"][tray] = np.sum(rops)
    rops = rop_matrix[tray][ROOs_id, cond_ind]
    rops = rops[rops > 0]
    production_rates["ROO._cond"][tray] = np.sum(rops)
    rops = rop_matrix[tray][ROs_id, cond_ind]
    rops = rops[rops > 0]
    production_rates["RO._cond"][tray] = np.sum(rops)

# %%
####### Calculate alphas ########
alpha1 = (rates["R._Add"] + rates["R.+O2"]) / (
    rates["R._Add"]
    + rates["R.+O2"]
    + rates["R._outlet"]
    + rates["R._evap"]
    + rates["R._Habs"]
    + rates["R._Recomb"]
    + rates["R._Disprop"]
    + rates["R._CycEther"]
)
alpha2 = (rates["ROO._Add"]) / (
    rates["ROO._Add"]
    + rates["ROO._outlet"]
    + rates["ROO._evap"]
    + rates["ROO._Habs_ROO."]
    + rates["ROO._Recomb"]
    + rates["ROO._Disprop"]
    + rates["ROO._eli"]
)

if np.all(np.isnan(alpha2)):
    alphas = alpha1
else:
    alphas = alpha1 * (
        1 - (1 - alpha2) * (rates["R.+O2"]) / (rates["R._Add"] + rates["R.+O2"])
    )

# ###### Save results #########
# with open(os.path.join(results_directory, f'alpha_rates.yml'), 'w+') as f:
#     yaml.dump([[alpha1.tolist(), alpha2.tolist(), alphas.tolist()], {key: value.tolist() for key, value in rates.items()}, {key: value.tolist() for key, value in production_rates.items()}], f)


def Wn(n, alpha):
    return n * (1 - alpha) ** 2 * alpha ** (n - 1)


print(alpha1)
print(alpha2)
print(alphas)

ASFparams = dict()
ASFparams["Wn oligomer"] = [0.0 for i in trays]
ASFparams["Wn dimer"] = [0.0 for i in trays]
ASFparams["n oligomer"] = [0.0 for i in trays]
ASFparams["#db oligomer"] = [0.0 for i in trays]

for tray in trays:
    ASFparams["Wn oligomer"][tray] = float(
        sum(Wn(n, alphas[tray]) for n in range(3, 100))
    )
    ASFparams["Wn dimer"][tray] = float(Wn(2, alphas[tray]))
for tray in trays:
    ASFparams["n oligomer"][tray] = float(
        sum(n * Wn(n, alphas[tray]) for n in range(3, 100))
        / ASFparams["Wn oligomer"][tray]
    )
    ASFparams["#db oligomer"][tray] = float(
        sum(n * Wn(n, alphas[tray]) for n in range(3, 100))
        / ASFparams["Wn oligomer"][tray]
    )

with open(os.path.join(results_directory, "ASFWnparams.yml"), "w+") as f:
    yaml.dump(ASFparams, stream=f)


############### Diels-Alder alpha calculation ##############
print("Calculating the alpha of Diels-Alder dimers ...")
rates["dimer_DA"] = np.zeros(len(trays), dtype=np.float64)
rates["dimer_side_rxn"] = np.zeros(len(trays), dtype=np.float64)
rates["dimer_outlet"] = np.zeros(len(trays), dtype=np.float64)
rates["dimer_evap"] = np.zeros(len(trays), dtype=np.float64)

conjugated_diene = Group().from_adjacency_list(
    """
1  C u0 p0 c0 {2,D}
2  C u0 p0 c0 {1,D} {3,S}
3  C u0 p0 c0 {2,S} {4,D}
4  C u0 p0 c0 {3,D}
"""
)
print("Identifying monomers, dimers and trimers...")
# Find all monomers (C=C-C=C)
target_species = conjugated_diene
target_species_index = []
for i, spc in enumerate(liqspcs):
    if spc.label in initial_monomer_labels:
        target_species_index.append(i)

# Find all Diels Alder reactions in the model
certain_family = "Diels_alder_addition"
da_rxn_idxs = []
for i, rxn in enumerate(liqrxns):
    if rxn.family == certain_family:
        da_rxn_idxs.append(i)

dimer_idxs = []
# From all DA reactions find all dimers
# Note, there is a mix of dimerization and trimerization.
for da_idx in da_rxn_idxs:
    rxn = liqrxns[da_idx]
    if all(reactant.label in initial_monomer_labels for reactant in rxn.reactants):
        dimer_idxs.append(spc_label_index_dict[rxn.products[0].label])

non_da_rxn_idxs = list(set(range(len(liqrxns))) - set(da_rxn_idxs))
for tray in trays:
    dimer_rops = rop_matrix[tray][dimer_idxs, :]
    rops = dimer_rops[:, da_rxn_idxs]
    rates["dimer_DA"][tray] = rops[rops < 0].sum()
    rops = dimer_rops[:, non_da_rxn_idxs]
    rates["dimer_side_rxn"][tray] = rops[rops < 0].sum()
    # Evaporation and outlet
    outlet_ind = source_label_index_dict["outlet"]
    rates["dimer_outlet"][tray] = dimer_rops[:, outlet_ind].sum()
    evap_ind = source_label_index_dict["evap"]
    rates["dimer_evap"][tray] = dimer_rops[:, evap_ind].sum()

# Generate Trimerization reactions
print("Generating Trimerization reactions...")
all_pair_indexes = set(
    (spc_i, spc_j) for spc_i in dimer_idxs for spc_j in target_species_index
)
print("Number of trimerization pairs: {}".format(len(all_pair_indexes)))

raw_results = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=50)(
    delayed(generate_trimerization)(liqspcs[spc_i], liqspcs[spc_j])
    for spc_i, spc_j in all_pair_indexes
)

print("Summing up all ROPs for trimerization...")

for (spc_i, spc_j), reaction_list in zip(all_pair_indexes, raw_results):
    #     pair = tuple(sorted([spc_i, spc_j]))
    # if not trimerization_pairs.get(pair):  # To avoid double counting
    for tray in trays:
        conc_prod = conc[tray][spc_i] * conc[tray][spc_j]
        for rxn0 in reaction_list:
            rates["dimer_DA"][tray] += (
                -(rxn0.get_rate_coefficient(temp[tray]) * conc_prod) * Vliq
            )

alphas_DA = (rates["dimer_DA"]) / (
    rates["dimer_DA"] + rates["dimer_side_rxn"] + rates["dimer_outlet"]
)

print(alphas_DA)

###### Save results #########
with open(os.path.join(results_directory, "alpha_rates.yml"), "w+") as f:
    yaml.dump(
        [
            [alpha1.tolist(), alpha2.tolist(), alphas.tolist(), alphas_DA.tolist()],
            {key: value.tolist() for key, value in rates.items()},
            {key: value.tolist() for key, value in production_rates.items()},
            rxn_rates,
            [Rs_labels, ROOs_labels, ROs_labels],
        ],
        f,
    )
