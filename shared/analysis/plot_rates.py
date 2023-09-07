import numpy as np
import matplotlib.pyplot as plt

from rmgpy.chemkin import load_chemkin_file

chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/trace_oxygen_perturbed/film_mechanism/chem_annotated_film.inp"
species_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/trace_oxygen_perturbed/film_mechanism/species_dictionary_film.txt"
filmspcs, filmrxns = load_chemkin_file(chemkin_path, species_dict_path)

Ts = np.linspace(300, 400, 50)

CPD_radical_CDB_RAdds = []
CPDOO_radical_CDB_RAdds = []

for rxn in filmrxns:
    if len(rxn.reactants) == 2 and len(rxn.products) == 1:
        spcs = rxn.reactants
        if any(spc.label == "CDB" for spc in spcs):
            if any(spc.label == "[CH]1C=CC=C1(L)" for spc in spcs):
                CPD_radical_CDB_RAdds.append(rxn)
            elif any(spc.label == "[O]OC1C=CC=C1(L)" for spc in spcs):
                CPDOO_radical_CDB_RAdds.append(rxn)

plt.figure()

for rxn in CPD_radical_CDB_RAdds:
    ks = [rxn.get_rate_coefficient(T) for T in Ts]
    plt.plot(1000 / Ts, ks, label="cyclopentadienyl + CDB radical addition")
    print(rxn)
    for spc in rxn.reactants + rxn.products:
        print(spc.smiles)
    print(rxn.kinetics)

for rxn in CPDOO_radical_CDB_RAdds:
    ks = [rxn.get_rate_coefficient(T) for T in Ts]
    plt.plot(1000 / Ts, ks, label="cyclopentadienyl peroxyl + CDB radical addition")
    print(rxn)
    for spc in rxn.reactants + rxn.products:
        print(spc.smiles)
    print(rxn.kinetics)

plt.yscale("log")
plt.xlabel("1000 / T (1000/K)")
plt.ylabel("k ($\mathrm{m}^3$/(mol*s))")
plt.legend()
plt.tight_layout()
plt.savefig("Figures/reaction_rates_CPD.+CDB_CPDOO.+CDB.pdf", bbox_inches="tight")

AR_BD_RAdds = []
AR_CPD_HAbs = []
AR_O2_RRecomb = []
PR_BD_RAdds = []
PR_CPD_HAbs = []

for rxn in filmrxns:
    if len(rxn.reactants) == 2 and len(rxn.products) == 1:
        spcs = rxn.reactants
        if any(spc.label == "AR" for spc in spcs):
            if any(spc.label == "1,3-BUTADIENE(L)" for spc in spcs):
                if any(spc.label == "C=C[CH]CCC(=C)C" for spc in rxn.products):
                    AR_BD_RAdds.append(rxn)
            elif any(spc.smiles == "[O][O]" for spc in spcs):
                AR_O2_RRecomb.append(rxn)
        elif any(spc.label == "PR" for spc in spcs):
            if any(spc.label == "1,3-BUTADIENE(L)" for spc in spcs):
                if any(spc.label == "C=C[CH]COOCC(C)C" for spc in rxn.products):
                    PR_BD_RAdds.append(rxn)
    if len(rxn.reactants) == 2 and len(rxn.products) == 2:
        spcs = rxn.reactants + rxn.products
        if any(spc.label == "AR" for spc in spcs):
            if any(spc.label == "CYCLOPENTADIENE(L)" for spc in spcs):
                AR_CPD_HAbs.append(rxn)
        elif any(spc.label == "PR" for spc in spcs):
            if any(spc.label == "CYCLOPENTADIENE(L)" for spc in spcs):
                PR_CPD_HAbs.append(rxn)

plt.figure()

for rxn in AR_BD_RAdds:
    ks = [rxn.get_rate_coefficient(T) for T in Ts]
    plt.plot(1000 / Ts, ks, label="AR + BD radical addition")
    print(rxn)
    for spc in rxn.reactants + rxn.products:
        print(spc.smiles)
    print(rxn.kinetics)

for rxn in AR_CPD_HAbs:
    ks = [rxn.get_rate_coefficient(T) for T in Ts]
    plt.plot(1000 / Ts, ks, label="AR + CPD hydrogen abstraction")
    print(rxn)
    for spc in rxn.reactants + rxn.products:
        print(spc.smiles)
    print(rxn.kinetics)

for rxn in AR_O2_RRecomb:
    ks = [rxn.get_rate_coefficient(T) for T in Ts]
    plt.plot(1000 / Ts, ks, label="AR + $\mathrm{O}_2$ => PR")
    print(rxn)
    for spc in rxn.reactants + rxn.products:
        print(spc.smiles)
    print(rxn.kinetics)

for rxn in PR_BD_RAdds:
    ks = [rxn.get_rate_coefficient(T) for T in Ts]
    plt.plot(1000 / Ts, ks, label="PR + BD radical addition")
    print(rxn)
    for spc in rxn.reactants + rxn.products:
        print(spc.smiles)
    print(rxn.kinetics)

for rxn in PR_CPD_HAbs:
    ks = [rxn.get_rate_coefficient(T) for T in Ts]
    plt.plot(1000 / Ts, ks, label="PR + CPD hydrogen abstraction")
    print(rxn)
    for spc in rxn.reactants + rxn.products:
        print(spc.smiles)
    print(rxn.kinetics)

plt.yscale("log")
plt.xlabel("1000 / T (1000/K)")
plt.ylabel("k ($\mathrm{m}^3$/(mol*s))")
plt.legend()
plt.tight_layout()
plt.savefig("Figures/reaction_rates_AR+BD_PR+BD.pdf", bbox_inches="tight")
