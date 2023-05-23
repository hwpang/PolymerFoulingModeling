def get_rops(df, rop_name, loss_only=False, production_only=False, N=5):
    name_inds = df["rop_spcname"] == rop_name
    rop_rxncomments = df.loc[name_inds, "rop_rxncomment"]
    # rop_rxncomments = [rxnstr.replace("\n", " ") for rxnstr in rop_rxncomments]
    # rop_rxncomments = [rxnstr.replace(" H abstraction", "") for rxnstr in rop_rxncomments]
    rop_rxnstrs = df.loc[name_inds, "rop_rxnstr"]
    rops = df.loc[name_inds, "rop"]
    
    if loss_only:
        loss_inds = rops < 0
        rops = rops[loss_inds]
        rop_rxncomments = rop_rxncomments[loss_inds]
        rop_rxnstrs = rop_rxnstrs[loss_inds]
    elif production_only:
        prod_inds = rops > 0
        rops = rops[prod_inds]
        rop_rxncomments = rop_rxncomments[prod_inds]
        rop_rxnstrs = rop_rxnstrs[prod_inds]
    
    sorted_inds = rops.argsort()
    if len(sorted_inds) > N:
        sorted_inds = sorted_inds[:N]
    
    return rops[sorted_inds], rop_rxncomments[sorted_inds], rop_rxnstrs[sorted_inds]
