# -*- coding: utf-8 -*-
u"""Test that module imports.

You should delete the test once you have real tests.
Only necessary if you have no other tests so that
tox will work.
"""
from __future__ import absolute_import, division, print_function
import math
from pykern.pkdebug import pkdp
from pykern.pkcollections import PKDict
import pytest
from rslaser.rspulse.pulse import LaserPulse
import scipy.constants as const


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

def test_basic_instantiation1():
    lambda0 = 1e-6
    phE = const.h * const.c / lambda0
    chirp = 0.0
    z_waist = 0
    z_center = 0.0
    k=PKDict(
        phE=phE,
        nslice=3,
        chirp=chirp,
        # TODO (gurhar1133): format k {kv pairs ..., hermite kv pairs: {}, slice kv pairs: {}} ?
        w0=.1,
        a0=.01,
        dw0x=0.0,
        dw0y=0.0,
        z_waist=z_waist,
        dzwx=0.0,
        dzwy=0.0,
        tau_fwhm=0.1 / const.c / math.sqrt(2.),
        z_center=z_center,
        x_shift = 0.,
        y_shift=0.,
        d_to_w=z_waist - z_center,
    ).pkupdate(_LASER_PULSE_SLICE_DEFAULTS)
    l = LaserPulse(k)
    for s in l.slice:
        if s.phE != l.phE:
            raise AssertionError('LaserPulseSlice has different wavelength than pulse as a whole')


def test_basic_instantiation2():
    lambda0 = 1e-6
    phE = const.h * const.c / lambda0
    chirp = 0.01*phE
    z_waist = 0
    z_center = 0.0
    k=PKDict(
        phE=phE,
        nslice=3,
        chirp=chirp,
        # TODO (gurhar1133): format k {kv pairs ..., hermite kv pairs: {}, slice kv pairs: {}} ?
        w0=.1,
        a0=.01,
        dw0x=0.0,
        dw0y=0.0,
        z_waist=z_waist,
        dzwx=0.0,
        dzwy=0.0,
        tau_fwhm=0.1 / const.c / math.sqrt(2.),
        z_center=z_center,
        x_shift = 0.,
        y_shift=0.,
        d_to_w=z_waist - z_center,
    ).pkupdate(_LASER_PULSE_SLICE_DEFAULTS)
    l = LaserPulse(k)
    a = [s.phE for s in l.slice]
    assert len(set(a)) == len(a)


def test_slice_input_validators():
    assert False


def test_pulse_input_validators():

    assert False


def test_pulse_slice_instantiation():
    assert False
