import numpy as np

def get_rops(df, rop_name, loss_only=False, production_only=False, N=5):
    name_inds = df["rop_spcname"] == rop_name
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
    
    sorted_inds = np.argsort(np.abs(rops))[::-1]
    if len(sorted_inds) > N:
        sorted_inds = sorted_inds[:N]
    
    return rops.iloc[sorted_inds], rop_rxncomments.iloc[sorted_inds], rop_rxnstrs.iloc[sorted_inds]
