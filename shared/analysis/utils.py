import numpy as np


def get_last_row(df, name):
    return df.loc[len(df.index) - 1, name]


def get_liq_radicals_conc(df, labels, Vliq):
    return sum([get_last_row(df, name) for name in labels]) / Vliq


def get_film_radical_reactive_site_conc(df, name):
    return get_last_row(df, name) / get_last_row(df, "mass")


def get_film_rops(df, rop_names, loss_only=False, production_only=False, N=5):
    if isinstance(rop_names, str):
        rop_names = [rop_names]
    name_inds = df["rop_spcname"].isin(rop_names)
    rop_rxncomments = df.loc[name_inds, "rop_rxncomment"]
    rop_rxnstrs = df.loc[name_inds, "rop_rxnstr"]
    rops = df.loc[name_inds, "rop"]

    if loss_only:
        inds = rops < 0
    elif production_only:
        inds = rops > 0
    else:
        inds = rops != 0

    rops = rops[inds]
    rop_rxncomments = rop_rxncomments[inds]
    rop_rxnstrs = rop_rxnstrs[inds]

    sorted_inds = rops.abs().argsort()[::-1]
    if len(sorted_inds) > N:
        sorted_inds = sorted_inds[:N]

    return (
        rops.iloc[sorted_inds],
        rop_rxncomments.iloc[sorted_inds],
        rop_rxnstrs.iloc[sorted_inds],
    )


def get_film_rops_1D(
    film_rate_of_productions,
    tray,
    cell_inds,
    rop_names,
    loss_only=False,
    production_only=False,
    N=3,
):
    if isinstance(rop_names, str):
        rop_names = [rop_names]

    cell_ind = cell_inds[0]
    df = film_rate_of_productions[tray, cell_ind]
    df.set_index(["rop_rxnstr", "rop_spcname"], inplace=True, drop=False)

    name_inds = df["rop_spcname"].isin(rop_names)
    rop_rxncomments = df.loc[name_inds, "rop_rxncomment"]
    rop_rxnstrs = df.loc[name_inds, "rop_rxnstr"]
    rops = df.loc[name_inds, "rop"]

    for cell_ind in cell_inds[1:]:
        df = film_rate_of_productions[tray, cell_ind]
        df.set_index(["rop_rxnstr", "rop_spcname"], inplace=True, drop=False)
        name_inds = df["rop_spcname"].isin(rop_names)
        rops += df.loc[name_inds, "rop"]

    if loss_only:
        inds = rops < 0
    elif production_only:
        inds = rops > 0
    else:
        inds = rops != 0

    rops = rops[inds]
    rop_rxncomments = rop_rxncomments[inds]
    rop_rxnstrs = rop_rxnstrs[inds]

    sorted_inds = rops.abs().argsort()[::-1]
    if len(sorted_inds) > N:
        sorted_inds = sorted_inds[:N]

    return (
        rops.iloc[sorted_inds],
        rop_rxncomments.iloc[sorted_inds],
        rop_rxnstrs.iloc[sorted_inds],
    )


def get_liquid_rops(
    df, rop_names, loss_only=False, production_only=False, radicals_only=False, N=5
):
    name_inds = df["rop_spcname"].isin(rop_names)
    rop_sourcestrings = df.loc[name_inds, "rop_sourcestring"]
    rops = df.loc[name_inds, "rop"]

    if loss_only:
        inds = rops < 0
    elif production_only:
        inds = rops > 0
    else:
        inds = rops != 0

    rops = rops[inds]
    rop_sourcestrings = rop_sourcestrings[inds]

    if radicals_only:
        inds = []
        for ind, sourcestring in enumerate(rop_sourcestrings):
            if "<=>" in sourcestring:
                r_smis, p_smis = sourcestring.split("<=>")

                p_rad = 0
                for p_smi in p_smis.split("+"):
                    if p_smi.count("[") == 1:
                        p_rad += 1

                r_rad = 0
                for r_smi in r_smis.split("+"):
                    if r_smi.count("[") == 1:
                        r_rad += 1

                if p_rad - r_rad != 0:
                    inds.append(ind)
            else:
                inds.append(ind)

        rops = rops.iloc[inds]
        rop_sourcestrings = rop_sourcestrings.iloc[inds]

    sorted_inds = rops.abs().argsort()[::-1]
    if len(sorted_inds) > N:
        sorted_inds = sorted_inds[:N]

    return rops.iloc[sorted_inds], rop_sourcestrings.iloc[sorted_inds]


def get_reaction_family(rmg, rxn):
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


def extrapolate_conc(concs, zs):
    if len(concs) == 1:
        return concs[0]
    else:
        return 10 ** (
            np.log10(concs[0])
            - (zs[1] - zs[0])
            * (np.log10(concs[1]) - np.log10(concs[0]))
            / (zs[2] - zs[1])
        )
