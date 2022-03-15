# -*- coding: utf-8 -*-
u"""Tests for instantiation of LaserCavity, LaserPulse
and LaserPulseSlice
"""
from __future__ import absolute_import, division, print_function
import math
import numpy as np
from pykern.pkdebug import pkdp, pkdlog
from pykern.pkcollections import PKDict
import pytest
from rslaser.pulse import pulse
from rslaser.cavity import laser_cavity
import scipy.constants as const


def pulse_instantiation_test(pulse, field):
    for s in pulse.slice:
        if getattr(s, field) != getattr(pulse, field):
            raise AssertionError(f'LaserPulseSlice has different {field} than pulse as a whole')


def test_no_params_instantiation():
    l = pulse.LaserPulse()


def slice_instantiation_test(pulse, field):
    a = [getattr(s, field) for s in pulse.slice]
    assert len(set(a)) == len(a)


def test_basic_instantiation1():
    l = pulse.LaserPulse()
    pulse_instantiation_test(l, 'phE')


def test_basic_instantiation2():
    k = PKDict(chirp=0.01*(const.h * const.c / const.e / 1e-6))
    l = pulse.LaserPulse(k)
    slice_instantiation_test(l, 'phE')


def test_basic_instantiation3():
    l = pulse.LaserPulse()
    pulse_instantiation_test(l, '_lambda0')


def test_basic_instantiation4():
    k = k = PKDict(chirp=0.01*(const.h * const.c / const.e / 1e-6))
    l = pulse.LaserPulse(k)
    slice_instantiation_test(l, '_lambda0')


def test_basic_pulse_slice_instantiation():
    s = [pulse.LaserPulseSlice(i) for i in range(10)]


def test_pulse_slice_input_validators_type():
    c = [[i, 0] for i in range(10)]
    for a in c:
        trigger_exception_test(pulse.LaserPulseSlice, a)


def test_wrong_input():
    p = PKDict(PHHe=0.1)
    trigger_exception_test(pulse.LaserPulse, p)


def test_wrong_input2():
    p = PKDict(slice_params=PKDict(blonk=9))
    trigger_exception_test(pulse.LaserPulse, p)


def test_wrong_filling():
    p = PKDict(slice_params=PKDict(pulseE=9))
    l = pulse.LaserPulse(p)


def test_wrong_filling2():
    p = PKDict(slice_params=PKDict())
    l = pulse.LaserPulse(p)


def test_pulse_input_validators_type():
    trigger_exception_test(pulse.LaserPulse, [])


def test_correct_slice_params_type():
    k = PKDict(slice_params=[])
    trigger_exception_test(pulse.LaserPulse, k)


def test_correct_slice_params():
    k = PKDict(
        slice_params=PKDict(
            foo='bar', hello='world'
        )
    )
    trigger_exception_test(pulse.LaserPulse, k)


def test_laser_pulse_slice_index():
    a = '10'
    trigger_exception_test(pulse.LaserPulseSlice, a)


def test_laser_cavity():
    c = laser_cavity.LaserCavity()


def test_laser_cavity_fail():
    k = PKDict(n3='fail')
    trigger_exception_test(laser_cavity.LaserCavity, k)


def test_cavity_partial_pulse_params():
    k = PKDict(pulse_params=PKDict(chirp=0.01*(const.h * const.c / 1e-6)))
    c = laser_cavity.LaserCavity(k)


def test_cavity_partial_pulse_params2():
    k = PKDict(
        pulse_params=PKDict(
            chirp=0.01*(const.h * const.c / 1e-6),
            slice_params=PKDict(
                propLen=2,
            )
        )
    )
    c = laser_cavity.LaserCavity(k)


def test_cavity_partial_pulse_params3():
    k = PKDict(
        pulse_params=PKDict(
            slice_params=PKDict(
                propLen=5,
            )
        )
    )
    c = laser_cavity.LaserCavity(k)

def test_cavity_partial_pulse_params4():
    k = PKDict()
    c = laser_cavity.LaserCavity(k)


def test_cavity_wrong_in_type():
    trigger_exception_test(laser_cavity.LaserCavity, 'wrong')


def trigger_exception_test(call, args):
    try:
        if len(args) > 0 and type(args) != PKDict and type(args) != str:
            l = call(*args)
        else:
            l = call(args)
    except Exception as e:
        pkdlog('EXCEPTION:{}, with message "{}" triggered by call: {} and args {}', type(e), e, call, args)
        return e
    assert False
