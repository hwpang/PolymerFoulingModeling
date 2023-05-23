def get_rops(df, rop_name, loss_only=False, production_only=False, N=5):
    name_inds = df["rop_spcname"] == rop_name
    rop_rxncomments = df.loc[name_inds, "rop_rxncomment"]
    # rop_rxncomments = [rxnstr.replace("\n", " ") for rxnstr in rop_rxncomments]
    # rop_rxncomments = [rxnstr.replace(" H abstraction", "") for rxnstr in rop_rxncomments]
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
    
    sorted_inds = rops.argsort()[::-1]
    print(rops)
    print(sorted_inds)
    if len(sorted_inds) > N:
        sorted_inds = sorted_inds[:N]
    
    return rops[sorted_inds], rop_rxncomments[sorted_inds], rop_rxnstrs[sorted_inds]
