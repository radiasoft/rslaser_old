u"""Default values for highdimensional inputs
Copyright (c) 2021 RadiaSoft LLC. All rights reserved
"""
import math
from pykern.pkcollections import PKDict
import scipy.constants as const

PHE_DEFAULT = const.h * const.c / 1e-6
Z_WAIST_DEFAULT = 0
Z_CENTER_DEFAULT = 0
LASER_PULSE_SLICE_DEFAULTS = PKDict(
    sigrW=0.000186,
    propLen=15,
    pulseE=0.001,
    poltype=1,
    sampFact=5,
    numsig=3.,
    mx=0,
    my=0
)
LASER_PULSE_DEFAULTS = PKDict(
        phE=PHE_DEFAULT,
        nslice=3,
        chirp=0,
        w0=.1,
        a0=.01,
        dw0x=0.0,
        dw0y=0.0,
        z_waist=Z_WAIST_DEFAULT,
        dzwx=0.0,
        dzwy=0.0,
        tau_fwhm=0.1 / const.c / math.sqrt(2.),
        z_center=Z_CENTER_DEFAULT,
        x_shift = 0.,
        y_shift=0.,
        d_to_w=Z_WAIST_DEFAULT - Z_CENTER_DEFAULT,
        slice_params=LASER_PULSE_SLICE_DEFAULTS,
)
LASER_CAVITY_DEFAULTS = PKDict(
    drift_right_length=0.5,
    drift_left_length=0.5,
    lens_left_focal_length=0.2,
    lens_right_focal_length=0.2,
    n0 = 1.75,
    n2 = 0.001,
    L_half_cryst=0.2,
)