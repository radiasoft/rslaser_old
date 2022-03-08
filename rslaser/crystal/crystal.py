import math
import numpy as np
import srwlib
from array import array
from pykern.pkcollections import PKDict

import rslaser.optics
from rslaser.optics.element import *

class CrystalSlice:
    """
    This class represents a slice of a crystal in a laser cavity.
    Args:
        length
        n0: on-axis index of refraction
        n2: transverse variation of index of refraction
            n(r) = n0 - 0.5 n2 r^2

    To be added: alpha0, alpha2 laser gain parameters

    Note: Initially, these parameters are fixed. Later we will update
    these parameters as the laser passes through.
    """
    def __init__(self,n0,n2):
        self.n0=n0
        self.n2=n2

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
