# -*- coding: utf-8 -*-
u"""Test that module imports.

You should delete the test once you have real tests.
Only necessary if you have no other tests so that
tox will work.
"""
from __future__ import absolute_import, division, print_function
import pytest
from pykern.pkdebug import pkdp
from pykern.pkcollections import PKDict
import rslaser.rspulse.pulse as rsp

def test_hello_world():
    return 'hello world'

def test_basic_instantiation():
    import rslaser.rspulse.pulse
    k=PKDict(
        d_lambda=0,
        lambda0=1e-6,
        nslice=3,
    )
    lp = rsp.LaserPulse(k)
    # for i, s in enumerate(lp.slice):
    #     print(f'slice {i} energy {s._phE}')
    # print(f'PulseEnergy was {lp.phE}')

    # instantiate LaserPulse instance with zero chirp (i.e. d_lambda = 0),
    # central wavelength lambda0=1e-6 meters and with 3 slices;
    # query the wavelength (i.e. _lambda0) associated with each slice;
    # it should be the same as for the pulse as a whole
