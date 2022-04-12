# -*- coding: utf-8 -*-
u"""Tests for instantiation of LaserCavity, LaserPulse
and LaserPulseSlice
"""
from __future__ import absolute_import, division, print_function
from pykern.pkdebug import pkdp, pkdlog
from pykern.pkcollections import PKDict
import pykern.pkunit
import array
import pytest
import copy
import srwlib
from rslaser.pulse import pulse
from rslaser.optics.wavefront import WavefrontSensor, InvalidWaveFrontSensorInputError


EPSILON = 1e-2
WFR_ATTRS_LIST = [
    'arElecPropMatr',
    'arEx',
    'arEy',
    'arMomX',
    'arMomY',
    'arWfrAuxData',
    'avgPhotEn',
    'dRx',
    'dRy',
    'presCA',
    'presFT',
    'unitElFld',
    'unitElFldAng',
    'xc',
    'yc',
    ]


def _check_epsilon_diff(val1, val2, epsilon, message):
    if val1 != 0:
        if not (val1 - val2)/val1 < epsilon:
            pykern.pkunit.pkfail(message)
    else:
        if not abs(val2) < epsilon:
            pykern.pkunit.pkfail(message)


def test_instantiation01():
    WavefrontSensor('w1', 2.0)
    with pykern.pkunit.pkexcept(InvalidWaveFrontSensorInputError):
        WavefrontSensor('w1', '2.0')


def test_propagation01():
    from pykern import pkunit, pkjson
    p = pulse.LaserPulse()
    w = WavefrontSensor('w1', 2.0)
    r = w.propagate(p)
    actual = PKDict()
    for a in WFR_ATTRS_LIST:
        v = getattr(r, a)
        actual.update({a: v})
    pkunit.empty_work_dir()
    pkunit.file_eq('res.json', actual=actual)


def test_propagation02():
    p = pulse.LaserPulse()
    wfs = WavefrontSensor('w1', 2.0)
    if not type(wfs.propagate(p)) == srwlib.SRWLWfr:
        pykern.pkunit.pkfail(f'WavefrontSensor.propagate failed to return type {srwlib.SRWLWfr}')
    p = pulse.LaserPulseSlice(0)
    with pykern.pkunit.pkexcept(InvalidWaveFrontSensorInputError):
        wfs.propagate(p)


def test_propagation043():
    p = pulse.LaserPulse(PKDict(nslice=1))
    pc = copy.deepcopy(p)
    wfs = WavefrontSensor('w1', 0.0)
    res = wfs.propagate(p)
    for a in WFR_ATTRS_LIST:
        o = getattr(pc.slice[0].wfr, a)
        n = getattr(res, a)
        if type(o) == array.array:
            for i, v in enumerate(o):
                _check_epsilon_diff(v, n[i], EPSILON,
                    f'epsilon check failed with vals: {v}, {n[i]} on attr: {a}, index: {i}')
        else:
            _check_epsilon_diff(o, n, EPSILON,
                f'epsilon check failed with vals: {v}, {n[i]} on attr: {a}')
