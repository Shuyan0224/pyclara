import numpy as _np

def covarianceFromBeamDict(beam_dict) :
    x     = beam_dict["x"]
    xp    = beam_dict["xp"]
    y     = beam_dict["y"]
    yp    = beam_dict["yp"]
    zeta  = beam_dict["zeta"]  # zeta = (t-t0)*c
    delta = beam_dict["delta"] # (p-p0)/p0

    cov_matrix = _np.cov([x, xp, y, yp, zeta, delta])
    return cov_matrix