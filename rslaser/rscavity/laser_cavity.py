import math
import numpy as np

import rslaser.rsoptics
from rslaser.rsoptics.element import *
from rslaser.rsoptics.wavefront import *
import rslaser.rspulse
from rslaser.rspulse.pulse import *

from array import array
from pykern.pkcollections import PKDict


_REQUIRED_CAVITY_PARAMS = ['drift_right_length', 'drift_left_length',
    'lens_left_focal_length', 'lens_right_focal_length', 'n0', 'n2',
    'L_half_cryst', 'pulse_params']


class InvalidLaserCavityInputError(Exception):
    pass


def _validate_params(input_params):
    check_fields(input_params, _REQUIRED_CAVITY_PARAMS, InvalidLaserCavityInputError, 'LaserCavity')


class LaserCavity:
    """
    create laser cavity

    Args:
        params (PKDict):
            required fields:
                drift_right_length
                drift_left_length
                lens_left_focal_length
                lens_right_focal_length
                n0
                n2
                L_half_cryst
                pulse_params (PKDict): see LaserPulse docs
    """
    __LASER_CAVITY_DEFAULTS = PKDict(
        drift_right_length=0.5,
        drift_left_length=0.5,
        lens_left_focal_length=0.2,
        lens_right_focal_length=0.2,
        n0 = 1.75,
        n2 = 0.001,
        L_half_cryst=0.2,
        pulse_params=PKDict()
    )

    def __init__(self, params=None):
        params = self.get_params(params)
        self.laser_pulse = LaserPulse(params.pulse_params)
        self.crystal_right = Crystal(params.n0,params.n2,params.L_half_cryst)
        self.crystal_left = Crystal(params.n0,params.n2,params.L_half_cryst)
        self.drift_right = Drift(params.drift_right_length)
        self.drift_left = Drift(params.drift_left_length)
        self.lens_right = Lens(params.lens_right_focal_length)
        self.lens_left  = Lens(params.lens_left_focal_length)

    def get_params(self, params):
        if params == None:
            return self.__LASER_CAVITY_DEFAULTS.copy()
        validate_type(params, PKDict, LaserCavity, 'params', InvalidLaserCavityInputError)
        for k in self.__LASER_CAVITY_DEFAULTS:
            if k not in params:
                params[k] = self.__LASER_CAVITY_DEFAULTS[k]
        _validate_params(params)
        return params

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
