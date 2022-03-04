# -*- coding: utf-8 -*-
u"""Test that module imports.

You should delete the test once you have real tests.
Only necessary if you have no other tests so that
tox will work.
"""
from __future__ import absolute_import, division, print_function
import math
from pykern.pkdebug import pkdp, pkdlog
from pykern.pkcollections import PKDict
import pytest
from rslaser.rspulse.pulse import LaserPulse, LaserPulseSlice, InvalidLaserPulseInputError
from rslaser.rscavity.laser_cavity import LaserCavity
import scipy.constants as const



_PHE_DEFAULT = const.h * const.c / 1e-6
_Z_WAIST_DEFAULT = 0
_Z_CENTER_DEFAULT = 0
_LASER_PULSE_SLICE_DEFAULTS = PKDict(
    sigrW=0.000186,
    propLen=15,
    pulseE=0.001,
    poltype=1,
    sampFact=5,
    numsig=3.,
    mx=0,
    my=0
)
_LASER_PULSE_DEFAULTS = PKDict(
        phE=_PHE_DEFAULT,
        nslice=3,
        chirp=0,
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
_LASER_CAVITY_DEFAULTS = PKDict(
    drift_right_length=0.5,
    drift_left_length=0.5,
    lens_left_focal_length=0.2,
    lens_right_focal_length=0.2,
    n0 = 1.75,
    n2 = 0.001,
    L_half_cryst=0.2,
)


def pulse_instantiation_test(pulse, field):
    for s in pulse.slice:
        if getattr(s, field) != getattr(pulse, field):
            raise AssertionError(f'LaserPulseSlice has different {field} than pulse as a whole')


def slice_instantiation_test(pulse, field):
    a = [getattr(s, field) for s in pulse.slice]
    assert len(set(a)) == len(a)


def test_basic_instantiation1():
    l = LaserPulse(_LASER_PULSE_DEFAULTS)
    pulse_instantiation_test(l, 'phE')


def test_basic_instantiation2():
    k = _LASER_PULSE_DEFAULTS.copy()
    k.chirp = 0.01*_PHE_DEFAULT
    l = LaserPulse(k)
    slice_instantiation_test(l, 'phE')


def test_basic_instantiation3():
    l = LaserPulse(_LASER_PULSE_DEFAULTS)
    pulse_instantiation_test(l, '_lambda0')


def test_basic_instantiation4():
    k = _LASER_PULSE_DEFAULTS.copy()
    k.chirp = 0.01*_PHE_DEFAULT
    l = LaserPulse(k)
    slice_instantiation_test(l, '_lambda0')


def test_basic_pulse_slice_instantiation():
    s = [LaserPulseSlice(i, _LASER_PULSE_DEFAULTS) for i in range(10)]


def test_slice_input_validators():
    k = _LASER_PULSE_DEFAULTS.copy()
    k.pkdel('slice_params')
    c = [(i, k) for i in range(10)]
    for a in c:
        trigger_exception_test(LaserPulseSlice, *a)


def test_pulse_slice_input_validators_type():
    c = [(i, 0) for i in range(10)]
    for a in c:
        trigger_exception_test(LaserPulseSlice, *a)


def test_pulse_input_validators_type():
    trigger_exception_test(LaserPulse, [])


def test_pulse_input_validators_fields():
    k = _LASER_PULSE_DEFAULTS.copy()
    k.pkdel('slice_params')
    trigger_exception_test(LaserPulse, k)


def test_correct_slice_params_type():
    k = _LASER_PULSE_DEFAULTS.copy()
    k.slice_params = []
    trigger_exception_test(LaserPulse, k)


def test_correct_slice_params():
    k = _LASER_PULSE_DEFAULTS.copy()
    k.slice_params = PKDict(foo='bar', hello='world')
    trigger_exception_test(LaserPulse, k)


def test_laser_pulse_slice_index():
    a = ('10', _LASER_PULSE_DEFAULTS)
    trigger_exception_test(LaserPulseSlice, *a)


def test_laser_cavity():
    c = LaserCavity(
        PKDict(
            **_LASER_CAVITY_DEFAULTS,
            pulse_params=_LASER_PULSE_DEFAULTS,
        )
    )

def test_laser_cavity_fail():
    k = _LASER_CAVITY_DEFAULTS.copy()
    k.pkdel('n0')
    trigger_exception_test(LaserCavity, k)


def trigger_exception_test(call, *args):
    try:
        l = call(*args)
    except Exception as e:
        pkdlog('EXCEPTION:{}, with message "{}" triggered by call: {} and args {}', type(e), e, call, args)
        return e
    assert False
