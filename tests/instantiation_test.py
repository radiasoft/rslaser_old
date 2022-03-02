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


def test_basic_instantiation1():
    # instantiate LaserPulse instance with zero chirp (i.e. d_lambda = 0),
    # central wavelength lambda0=1e-6 meters and with 3 slices;
    # query the wavelength (i.e. _lambda0) associated with each slice;
    # it should be the same as for the pulse as a whole
    k=PKDict(
        d_lambda=0,
        lambda0=1e-6,
        nslice=3,
        # TODO (gurhar1133): instantiation of LaserPulse fails without
        w0=.1, # fields below this comment. Needed for calling GaussHermite(). What are reasonable values?
        a0=.01, # current values pulled from TestCavityFull.ipynb
        dw0x=0.0,
        dw0y=0.0,
        z_waist=0,
        dzwx=0.0,
        dzwy=0.0,
        tau_fwhm=0.1 / const.c / math.sqrt(2.),
        z_center=0.0,
        x_shift = 0.,
        y_shift=0.,
    )
    k.d_to_w = k.z_waist - k.z_center

    l = LaserPulse(k)
    for s in l.slice:
        if s._lambda0 != l._lambda0:
            raise AssertionError('LaserPulseSlice has different wavelength than pulse as a whole')


def test_basic_instantiation2():
    # do the same with d_lambda = lambda0/10; in this case,
    # each slice should have a different wavelength,
    # but probably won't because of a bug
    k=PKDict(
        lambda0=1e-6,
        nslice=3,
        # TODO (gurhar1133): instantiation of LaserPulse fails without
        w0=.1, # fields below this comment. Needed for calling GaussHermite(). What are reasonable values?
        a0=.01, # current values pulled from TestCavityFull.ipynb
        dw0x=0.0,
        dw0y=0.0,
        z_waist=0,
        dzwx=0.0,
        dzwy=0.0,
        tau_fwhm=0.1 / const.c / math.sqrt(2.),
        z_center=0.0,
        x_shift = 0.,
        y_shift=0.,
    )
    k.d_to_w = k.z_waist - k.z_center
    k.d_lambda = k.lambda0/10
    l = LaserPulse(k)
    a = [s._lambda0 for s in l.slice]
    assert len(set(a)) == len(a)


def test_slice_input_validators():
    assert False


def test_pulse_input_validators():
    assert False


def test_pulse_slice_instantiation():
    assert False
