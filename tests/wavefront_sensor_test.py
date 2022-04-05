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
from rslaser.optics.wavefront import WavefrontSensor
from rslaser.cavity import laser_cavity
import scipy.constants as const
import test_utils


# TODO (gurhar1133): check logs, the exceptions being triggered arent the right ones


def test_wfs_instantiation():
    WavefrontSensor('w1', 2.0)


def test_propagate():
    p = pulse.LaserPulse()
    wfs = WavefrontSensor('w1', 2.0)
    # TODO (gurhar1133): check attrs other than avgPhotEn
    assert wfs.propagate(p).avgPhotEn == 1.55


def test_fail_instantiate():
    test_utils.trigger_exception_test(WavefrontSensor, ('w1', '2'))


def test_propagate_fail():
    p = pulse.LaserPulseSlice(0)
    wfs = WavefrontSensor('w1', 2.0)
    test_utils.trigger_exception_test(wfs.propagate, p)


def test_propagation_vals():
    '''IN: 1 slice and distance_from_pulse_center = 0, OUT: laserPulseSlice.wfr should remain the same'''
    # TODO (gurhar1133): get clarity on this one, which attr of srw object are we comparing with res?
    # rn checking avgPhotEn
    p = pulse.LaserPulse(PKDict(nslice=1))
    wfs = WavefrontSensor('w1', 0.0)
    res = wfs.propagate(p)
    assert res.avgPhotEn == p.slice[0].wfr.avgPhotEn

