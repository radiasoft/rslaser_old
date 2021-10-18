import math
import numpy as np
import srwlib
from array import array
from pykern.pkcollections import PKDict

class LaserCavity:
    """
    create laser cavity
    """
    def __init__(self,kwargs):
        k=PKDict(kwargs).pksetdefault(
            n0 = 1.75,
            n2 = 0.001,
            L_half_cryst=0.2,
            laser_pulse_length=0.1,
            nslice=11,
            drift_right_length=0.5,
            drift_left_length=0.5,
            lens_left_focal_length=0.2,
            lens_right_focal_length=0.2,
            sigrW=0.000437,
            propLen=15,
            sig_s=0.1,
            pulseE=0.001,
            poltype=1,
            phE=1.55,
            energyChirp=.01,
            sampFact=5,
            mx=0,
            my=0
        )
        self.laser_pulse = LaserPulse(k)
        self.crystal_right = Crystal(k.n0,k.n2,k.L_half_cryst)

        self.crystal_left = Crystal(k.n0,k.n2,k.L_half_cryst)
        self.drift_right = Drift(k.drift_right_length)
        self.drift_left = Drift(k.drift_left_length)
        self.lens_right = Lens(k.lens_right_focal_length)
        self.lens_left  = Lens(k.lens_left_focal_length)
