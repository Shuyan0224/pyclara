import numpy as _np
from ._covariance import covarianceFromBeamDict as _covarianceFromBeamDict

def emittance(beam_dict) :
    cm = _covarianceFromBeamDict(beam_dict)

    emit_x = _np.sqrt(_np.linalg.det(cm[0:2,0:2]))
    emit_y = _np.sqrt(_np.linalg.det(cm[2:4,2:4]))
    emit_z = _np.sqrt(_np.linalg.det(cm[4:6,4:6]))

    return {"emit_x":emit_x,
            "emit_y":emit_y,
            "emit_z":emit_z}


