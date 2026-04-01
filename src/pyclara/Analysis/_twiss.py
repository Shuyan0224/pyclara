import numpy as _np

from ._covariance import covarianceFromBeamDict as _covarianceFromBeamDict
from ._emittance import emittance as _emittanceFromBeamDict

def twissFromBeamDict(beam_dict):
    covm = _covarianceFromBeamDict(beam_dict)
    emit = _emittanceFromBeamDict(beam_dict)

    beta_x  =  covm[0][0]/emit["emit_x"]
    alpha_x = -covm[0][1]/emit["emit_x"]
    gamma_x =  covm[1][1]/emit["emit_x"]

    beta_y  =  covm[2][2]/emit["emit_y"]
    alpha_y = -covm[2][3]/emit["emit_y"]
    gamma_y =  covm[3][3]/emit["emit_y"]

    beta_z  =  covm[4][4]/emit["emit_z"]
    alpha_z = -covm[4][5]/emit["emit_z"]
    gamma_z =  covm[5][5]/emit["emit_z"]

    return {"beta_x":beta_x,
            "alpha_x":alpha_x,
            "gamma_x":gamma_x,
            "beta_y":beta_y,
            "alpha_y":alpha_y,
            "gamma_y":gamma_y,
            "beta_z":beta_z,
            "alpha_z":alpha_z,
            "gamma_z":gamma_z} | emit

def twissFromCovariance(covariance_matrix) :
    pass

