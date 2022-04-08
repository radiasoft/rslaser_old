# -*- coding: utf-8 -*-
u"""Tests for instantiation of LaserCavity, LaserPulse
and LaserPulseSlice
"""
from __future__ import absolute_import, division, print_function
from tabnanny import check
from tkinter import E
from pykern.pkdebug import pkdp, pkdlog
from pykern.pkcollections import PKDict
from pykern.pkunit import pkfail
import array
import pytest
import copy
import srwlib
from rslaser.pulse import pulse
from rslaser.optics.wavefront import WavefrontSensor
import test_utils


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


def test_wfs_instantiation():
    WavefrontSensor('w1', 2.0)


def test_propagate():
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


def check_epsilon_diff(val1, val2, epsilon, message):
    if val1 != 0:
        if not (val1 - val2)/val1 < epsilon:
            pkfail(message)
    else:
        if not abs(val2) < epsilon:
            pkfail(message)


def test_propagate_ret_type():
    p = pulse.LaserPulse()
    wfs = WavefrontSensor('w1', 2.0)
    assert type(wfs.propagate(p)) == srwlib.SRWLWfr


def test_fail_instantiate():
    test_utils.trigger_exception_test(WavefrontSensor, ['w1', '2.0'])


def test_propagate_fail():
    p = pulse.LaserPulseSlice(0)
    wfs = WavefrontSensor('w1', 2.0)
    test_utils.trigger_exception_test(wfs.propagate, p)


def test_propagation_vals():
    p = pulse.LaserPulse(PKDict(nslice=1))
    pc = copy.deepcopy(p)
    wfs = WavefrontSensor('w1', 0.0)
    res = wfs.propagate(p)
    for a in WFR_ATTRS_LIST:
        o = getattr(pc.slice[0].wfr, a)
        n = getattr(res, a)
        if type(o) == array.array:
            for i, v in enumerate(o):
                check_epsilon_diff(v, n[i], EPSILON,
                    f'epsilon check failed with vals: {v}, {n[i]} on attr: {a}, index: {i}')
        else:
            check_epsilon_diff(o, n, EPSILON,
                f'epsilon check failed with vals: {v}, {n[i]} on attr: {a}')

