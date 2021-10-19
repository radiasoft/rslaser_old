import math
import numpy as np

import rslaser.rsoptics
from rslaser.rsoptics.element import *
from rslaser.rsoptics.wavefront import *
import rslaser.rspulse
from rslaser.rspulse.pulse import *

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

    def propagate(self, num_cycles, callback=None):
        l = self.laser_pulse
        l._sxvals = []
        l._syvals = []
        current_position = 0

        # initial wavefront rms values
        vals = l.compute_middle_slice_intensity()
        positions = [current_position]
        if callback:
            callback(current_position, vals)

        for n in range(num_cycles):
            for element in (
                self.crystal_right,
                self.drift_right,
                self.lens_right,
                self.drift_right,
                self.crystal_right,
                self.crystal_left,
                self.drift_left,
                self.lens_left,
                self.drift_left,
                self.crystal_left,
            ):
                element.propagate(l)
                current_position += element.length
                vals = l.compute_middle_slice_intensity()
                positions.append(current_position)
                if callback:
                    callback(current_position, vals)

        return(positions, l._sxvals, l._syvals)
