# -*- coding: utf-8 -*-
u"""Test that module imports.

You should delete the test once you have real tests.
Only necessary if you have no other tests so that
tox will work.
"""
from __future__ import absolute_import, division, print_function
from curses import echo
import math
from pykern.pkdebug import pkdp
from pykern.pkcollections import PKDict
import pytest
from rslaser.rspulse.pulse import LaserPulse, LaserPulseSlice, InvalidLaserPulseInputError
import scipy.constants as const


_PHE_DEFAULT = const.h * const.c / 1e-6
_Z_WAIST_DEFAULT = 0
_Z_CENTER_DEFAULT = 0
_LASER_PULSE_SLICE_DEFAULTS = PKDict(
    sigrW=0.000186,
    propLen=15,
    sig_s=0.1,
    pulseE=0.001,
    poltype=1,
    sampFact=5,
    mx=0,
    my=0
)
_LASER_PULSE_DEFAULTS = PKDict(
        phE=_PHE_DEFAULT,
        nslice=3,
        chirp=0,
        # TODO (gurhar1133): format k {kv pairs ..., hermite kv pairs: {}, slice kv pairs: {}} ?
        w0=.1,
        a0=.01,
        dw0x=0.0,
        dw0y=0.0,
        z_waist=_Z_WAIST_DEFAULT,
        dzwx=0.0,
        dzwy=0.0,
        tau_fwhm=0.1 / const.c / math.sqrt(2.),
        z_center=_Z_CENTER_DEFAULT,
        x_shift = 0.,
        y_shift=0.,
        d_to_w=_Z_WAIST_DEFAULT - _Z_CENTER_DEFAULT,
        slice_params=_LASER_PULSE_SLICE_DEFAULTS,
)


def test_basic_instantiation1():
    l = LaserPulse(_LASER_PULSE_DEFAULTS)
    for s in l.slice:
        if s.phE != l.phE:
            raise AssertionError('LaserPulseSlice has different wavelength than pulse as a whole')


def test_basic_instantiation2():
    k = _LASER_PULSE_DEFAULTS.copy()
    k.chirp = 0.01*_PHE_DEFAULT
    l = LaserPulse(k)
    a = [s.phE for s in l.slice]
    assert len(set(a)) == len(a)


def test_basic_pulse_slice_instantiation():
    s = [LaserPulseSlice(i, _LASER_PULSE_DEFAULTS) for i in range(10)]


def test_slice_input_validators():
    k = _LASER_PULSE_DEFAULTS.copy()
    k.pkdel('slice_params')
    try:
        l = LaserPulse(k)
    except InvalidLaserPulseInputError as e:
        return e
    assert False


def test_pulse_input_validators_type():
    try:
        l = LaserPulse([])
    except InvalidLaserPulseInputError as e:
        return e
    assert False

def test_pulse_input_validators_fields():
    k=PKDict(
        phE=_PHE_DEFAULT,
        nslice=3,
        chirp=0,
        # TODO (gurhar1133): format k {kv pairs ..., hermite kv pairs: {}, slice kv pairs: {}} ?
        w0=.1,
        a0=.01,
        dw0x=0.0,
        dw0y=0.0,
        z_waist=_Z_WAIST_DEFAULT,
        dzwx=0.0,
        dzwy=0.0,
        tau_fwhm=0.1 / const.c / math.sqrt(2.),
        z_center=_Z_CENTER_DEFAULT,
        x_shift = 0.,
        y_shift=0.,
        d_to_w=_Z_WAIST_DEFAULT - _Z_CENTER_DEFAULT,
    )
    try:
        l = LaserPulse(k)
    except InvalidLaserPulseInputError as e:
        return e
    assert False
