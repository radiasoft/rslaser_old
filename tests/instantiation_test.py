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
from rslaser.utils import defaults
import scipy.constants as const


# TODO (gurhar1133): rewrite tests to reflect __DEFAULTS set on class now
def pulse_instantiation_test(pulse, field):
    for s in pulse.slice:
        if getattr(s, field) != getattr(pulse, field):
            raise AssertionError(f'LaserPulseSlice has different {field} than pulse as a whole')


def test_no_params_instantiation():
    l = LaserPulse()


def slice_instantiation_test(pulse, field):
    a = [getattr(s, field) for s in pulse.slice]
    assert len(set(a)) == len(a)


def test_basic_instantiation1():
    l = LaserPulse(defaults.LASER_PULSE_DEFAULTS)
    pulse_instantiation_test(l, 'phE')


def test_basic_instantiation2():
    k = defaults.LASER_PULSE_DEFAULTS.copy()
    k.chirp = 0.01*defaults.PHE_DEFAULT
    l = LaserPulse(k)
    slice_instantiation_test(l, 'phE')


def test_basic_instantiation3():
    l = LaserPulse(defaults.LASER_PULSE_DEFAULTS)
    pulse_instantiation_test(l, '_lambda0')


def test_basic_instantiation4():
    k = defaults.LASER_PULSE_DEFAULTS.copy()
    k.chirp = 0.01*defaults.PHE_DEFAULT
    l = LaserPulse(k)
    slice_instantiation_test(l, '_lambda0')


def test_basic_pulse_slice_instantiation():
    s = [LaserPulseSlice(i, defaults.LASER_PULSE_DEFAULTS) for i in range(10)]


def test_slice_input_validators():
    k = defaults.LASER_PULSE_DEFAULTS.copy()
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
    k = defaults.LASER_PULSE_DEFAULTS.copy()
    k.pkdel('slice_params')
    trigger_exception_test(LaserPulse, k)


def test_correct_slice_params_type():
    k = defaults.LASER_PULSE_DEFAULTS.copy()
    k.slice_params = []
    trigger_exception_test(LaserPulse, k)


def test_correct_slice_params():
    k = defaults.LASER_PULSE_DEFAULTS.copy()
    k.slice_params = PKDict(foo='bar', hello='world')
    trigger_exception_test(LaserPulse, k)


def test_laser_pulse_slice_index():
    a = ('10', defaults.LASER_PULSE_DEFAULTS)
    trigger_exception_test(LaserPulseSlice, *a)


def test_laser_cavity():
    c = LaserCavity(
        PKDict(
            **defaults.LASER_CAVITY_DEFAULTS,
            pulse_params=defaults.LASER_PULSE_DEFAULTS,
        )
    )

def test_laser_cavity_fail():
    k = defaults.LASER_CAVITY_DEFAULTS.copy()
    k.pkdel('n0')
    trigger_exception_test(LaserCavity, k)


def trigger_exception_test(call, *args):
    try:
        l = call(*args)
    except Exception as e:
        pkdlog('EXCEPTION:{}, with message "{}" triggered by call: {} and args {}', type(e), e, call, args)
        return e
    assert False
