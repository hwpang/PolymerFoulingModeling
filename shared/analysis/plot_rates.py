import numpy as np
import matplotlib.pyplot as plt

from rmgpy.chemkin import load_chemkin_file

chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/trace_oxygen_perturbed/film_mechanism/chem_annotated_film.inp"
species_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/trace_oxygen_perturbed/film_mechanism/species_dictionary_film.txt"
filmspcs, filmrxns = load_chemkin_file(chemkin_path, species_dict_path)

def get_rxn_smiles(rxn):
    return "+".join([spc.smiles for spc in rxn.reactants]) + "->" + "+".join([spc.smiles for spc in rxn.products])

def plot_rxn_rates(rxns, Ts, ax, label, linestyle="-"):
    for rxn in rxns:
        ks = [rxn.get_rate_coefficient(T) for T in Ts]
        rxn_smiles = get_rxn_smiles(rxn)
        ax.plot(1000 / Ts, ks, label=f"{label}", linestyle=linestyle)
        print(rxn.kinetics)
Ts = np.linspace(300, 400, 50)

fig, axs = plt.subplots(2, 3, figsize=(12, 7), sharex=True, sharey=True)

print("Plotting AR/PR + BD radical addition rates...")
AR_BD_RAdds = []
PR_BD_RAdds = []

for rxn in filmrxns:
    if len(rxn.reactants) == 2 and len(rxn.products) == 1:
        spcs = rxn.reactants
        if any(spc.label == "AR" for spc in spcs):
            if any(spc.label == "1,3-BUTADIENE(L)" for spc in spcs):
                if any(spc.label == "C=C[CH]CCC(=C)C" for spc in rxn.products):
                    AR_BD_RAdds.append(rxn)
        elif any(spc.label == "PR" for spc in spcs):
            if any(spc.label == "1,3-BUTADIENE(L)" for spc in spcs):
                if any(spc.label == "C=C[CH]COOCC(C)C" for spc in rxn.products):
                    PR_BD_RAdds.append(rxn)

ax = axs[0, 0]
ax.set_title("(a)", loc="left")
plot_rxn_rates(AR_BD_RAdds, Ts, ax, "AR + BD")
plot_rxn_rates(PR_BD_RAdds, Ts, ax, "PR + BD", linestyle="--")

ax.set_yscale("log")
ax.legend()

print("Plotting AR/PR + CPD radical addition rates...")
AR_CPD_RAdds = []
PR_CPD_RAdds = []

for rxn in filmrxns:
    if len(rxn.reactants) == 2 and len(rxn.products) == 1:
        spcs = rxn.reactants
        if any(spc.label == "AR" for spc in spcs):
            if any(spc.label == "CYCLOPENTADIENE(L)" for spc in spcs):
                if any(spc.smiles == "C=C(C)CC1C=C[CH]C1" for spc in rxn.products):
                    AR_CPD_RAdds.append(rxn)
        elif any(spc.label == "PR" for spc in spcs):
            if any(spc.label == "CYCLOPENTADIENE(L)" for spc in spcs):
                if any(spc.smiles == "CC(C)COOC1[CH]C=CC1" for spc in rxn.products):
                    PR_CPD_RAdds.append(rxn)

ax = axs[0, 1]
ax.set_title("(b)", loc="left")
plot_rxn_rates(AR_CPD_RAdds, Ts, ax, "AR + CPD")
plot_rxn_rates(PR_CPD_RAdds, Ts, ax, "PR + CPD", linestyle="--")

ax.set_yscale("log")
ax.legend()

print("Plotting cyclopentenyl RC./RCOO. + CDB radical addition rates...")
cyclopentenyl_radical_CDB_RAdds = []
cyclopentenylperoxyl_radical_CDB_RAdds = []

for rxn in filmrxns:
    if len(rxn.reactants) == 2 and len(rxn.products) == 1:
        spcs = rxn.reactants
        if any(spc.label == "CDB" for spc in spcs):
            if any(spc.label == "[CH]1C=CCC1(L)" for spc in spcs):
                cyclopentenyl_radical_CDB_RAdds.append(rxn)
            elif any(spc.label == "[O]OC1C=CCC1(L)" for spc in spcs):
                cyclopentenylperoxyl_radical_CDB_RAdds.append(rxn)

ax = axs[0, 2]
ax.axis("off")

ax = axs[1, 0]
ax.set_title("(c)", loc="left")
plot_rxn_rates(cyclopentenyl_radical_CDB_RAdds, Ts, ax, "1-cyclopenten-3-yl + CDB")
plot_rxn_rates(cyclopentenylperoxyl_radical_CDB_RAdds, Ts, ax, "1-cyclopenten-3-ylperoxyl + CDB", linestyle="--")

ax.set_yscale("log")
ax.legend(bbox_to_anchor=(0.0, -0.2), loc="upper left")

print("Plotting butenyl1 RC./RCOO. + CDB radical addition rates...")
butenyl1_radical_CDB_RAdds = []
butenyl1peroxyl_radical_CDB_RAdds = []

for rxn in filmrxns:
    if len(rxn.reactants) == 2 and len(rxn.products) == 1:
        spcs = rxn.reactants
        if any(spc.label == "CDB" for spc in spcs):
            if any(spc.label == "C=C[CH]C(L)" for spc in spcs):
                if any(spc.label == "C=CC(C)C([CH]CC)CC" for spc in rxn.products):
                    butenyl1_radical_CDB_RAdds.append(rxn)
            elif any(spc.label == "C=CC(C)O[O](L)" for spc in spcs):
                butenyl1peroxyl_radical_CDB_RAdds.append(rxn)

ax = axs[1, 1]
ax.set_title("(d)", loc="left")
plot_rxn_rates(butenyl1_radical_CDB_RAdds, Ts, ax, "1-buten-3-yl + CDB")
plot_rxn_rates(butenyl1peroxyl_radical_CDB_RAdds, Ts, ax, "1-buten-3-ylperoxyl + CDB", linestyle="--")

ax.set_yscale("log")
ax.legend(bbox_to_anchor=(0.0, -0.2), loc="upper left")

print("Plotting butenyl2 RC./RCOO. + CDB radical addition rates...")
butenyl2_radical_CDB_RAdds = []
butenyl2peroxyl_radical_CDB_RAdds = []

for rxn in filmrxns:
    if len(rxn.reactants) == 2 and len(rxn.products) == 1:
        spcs = rxn.reactants
        if any(spc.label == "CDB" for spc in spcs):
            if any(spc.label == "C=C[CH]C(L)" for spc in spcs):
                if any(spc.label == "CC=CCC([CH]CC)CC" for spc in rxn.products):
                    butenyl2_radical_CDB_RAdds.append(rxn)
            elif any(spc.label == "CC=CCO[O](L)" for spc in spcs):
                butenyl2peroxyl_radical_CDB_RAdds.append(rxn)

ax = axs[1, 2]
ax.set_title("(e)", loc="left")
plot_rxn_rates(butenyl2_radical_CDB_RAdds, Ts, ax, "2-buten-3-yl + CDB")
plot_rxn_rates(butenyl2peroxyl_radical_CDB_RAdds, Ts, ax, "2-buten-3-ylperoxyl + CDB", linestyle="--")

ax.set_yscale("log")
ax.legend(bbox_to_anchor=(0.0, -0.2), loc="upper left")

for ax in axs[:, 0]:
    ax.set_ylabel("k ($\mathrm{m}^3$/(mol*s))")
for ax in axs[-1, :]:
    ax.set_xlabel("1000 / T (1000/K)")

fig.subplots_adjust(wspace=0, hspace=0)
fig.align_labels()
fig.tight_layout()
fig.savefig("Figures/reaction_rates_AR+BD_PR+BD_RC.+CDB_RCOO.+CDB.pdf", bbox_inches="tight")
