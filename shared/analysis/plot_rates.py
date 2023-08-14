import numpy as np
import matplotlib.pyplot as plt

from rmgpy.chemkin import load_chemkin_file

chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/trace_oxygen_perturbed/film_mechanism/chem_annotated_film.inp"
species_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/trace_oxygen_perturbed/film_mechanism/species_dictionary_film.txt"
filmspcs, filmrxns = load_chemkin_file(chemkin_path, species_dict_path)

AR_BD_RAdds = []
AR_O2_RRecomb = []
PR_BD_RAdds = []

for rxn in filmrxns:
    if len(rxn.reactants) == 2 and len(rxn.products) == 1:
        spcs = rxn.reactants
        if any(spc.label == "AR" for spc in spcs):
            if any(spc.label == "1,3-BUTADIENE(L)" for spc in spcs):
                AR_BD_RAdds.append(rxn)
            elif any(spc.smiles == "[O][O]" for spc in spcs):
                AR_O2_RRecomb.append(rxn)
        elif any(spc.label == "PR" for spc in spcs):
            if any(spc.label == "1,3-BUTADIENE(L)" for spc in spcs):
                PR_BD_RAdds.append(rxn)

Ts = np.linspace(300, 400, 50)
plt.figure()

for rxn in AR_BD_RAdds:
    ks = [rxn.get_rate_coefficient(T) for T in Ts]
    plt.plot(1000 / Ts, ks, label="AR + BD radical addition")
    print(rxn.kinetics)

for rxn in AR_O2_RRecomb:
    ks = [rxn.get_rate_coefficient(T) for T in Ts]
    plt.plot(1000 / Ts, ks, label="AR + $\mathrm{O}_2$ => PR")
    print(rxn.kinetics)

for rxn in PR_BD_RAdds:
    ks = [rxn.get_rate_coefficient(T) for T in Ts]
    plt.plot(1000 / Ts, ks, label="PR + BD radical addition")
    print(rxn.kinetics)

plt.yscale("log")
plt.xlabel("1000 / T (1000/K)")
plt.ylabel("k ($\mathrm{m}^3$/(mol*s))")
plt.legend()
plt.tight_layout()
plt.savefig("Figures/reaction_rates.pdf", bbox_inches="tight")
