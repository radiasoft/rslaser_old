# -*- coding: utf-8 -*-
u"""Tests for instantiation of LaserCavity, LaserPulse
and LaserPulseSlice
"""
from __future__ import absolute_import, division, print_function
from tabnanny import check
from tkinter import E
from pykern.pkdebug import pkdp, pkdlog
from pykern.pkcollections import PKDict
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
    wfs = WavefrontSensor('w1', 2.0)
    res = wfs.propagate(p)
    res_attrs = PKDict()
    for attr in WFR_ATTRS_LIST:
        res_attrs.update({attr: getattr(res, attr)})
    # TODO (gurhar1133): for some reason diff is empty, fix test
    # so that it is not (and that there is no diff)
    pkunit.file_eq('res.json', actual=res_attrs)


def check_epsilon_diff(val1, val2, epsilon, attr):
    # TODO (gurhar1133): make sure this is correct
    def _check(x, y, epsilon):
        assert (x - y)/x < epsilon, f'epsilon check failed with vals: {x}, {y} on attr: {attr}'

    if val1 != 0:
        _check(val1, val2, epsilon)
    else:
        assert abs(val2) < epsilon, f'epsilon check failed with vals: {val1}, {val2} on attr: {attr}'


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
                check_epsilon_diff(v, n[i], EPSILON, a)
        else:
            check_epsilon_diff(o, n, EPSILON, a)

