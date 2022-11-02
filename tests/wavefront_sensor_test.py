# -*- coding: utf-8 -*-
u"""Tests for instantiation of LaserCavity, LaserPulse
and LaserPulseSlice
"""
from __future__ import absolute_import, division, print_function
from pykern.pkdebug import pkdp, pkdlog
from pykern.pkcollections import PKDict
import pykern.pksubprocess
import pykern.pkunit
import pykern.pkio
import re
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


def test_instantiation01():
    WavefrontSensor('w1', 2.0)
    with pykern.pkunit.pkexcept(InvalidWaveFrontSensorInputError):
        WavefrontSensor('w1', '2.0')


def test_propagation01():
    from pykern import pkunit
    from pykern import pkio

    data_dir = pkunit.data_dir()
    work_dir = pkunit.empty_work_dir()
    p = pulse.LaserPulse()
    w = WavefrontSensor('w1', 2.0)
    r = w.propagate(p)
    actual = PKDict()
    for a in WFR_ATTRS_LIST:
        v = getattr(r, a)
        actual.update({a: v})
    ndiff_files(
        data_dir.join("res.txt"),
        pkio.write_text(
            work_dir.join("res_actual.txt"),
            str(actual),
        ),
        work_dir.join("ndiff.out"),
        data_dir
    )


def test_propagation02():
    p = pulse.LaserPulse()
    wfs = WavefrontSensor('w1', 2.0)
    if not type(wfs.propagate(p)) == srwlib.SRWLWfr:
        pykern.pkunit.pkfail(f'WavefrontSensor.propagate failed to return type {srwlib.SRWLWfr}')
    p = pulse.LaserPulseSlice(0)
    with pykern.pkunit.pkexcept(InvalidWaveFrontSensorInputError):
        wfs.propagate(p)


# TODO (gurhar1133): this will be scrapped for file_eq when file_q
# supports ndiff
def ndiff_files(expect_path, actual_path, diff_file, data_dir):
    pykern.pksubprocess.check_call_with_signals(
        [
            "ndiff",
            actual_path,
            expect_path,
            data_dir.join("ndiff_conf.txt"),
        ],
        output=str(diff_file),
    )

    d = pykern.pkio.read_text(diff_file)
    if re.search("diffs have been detected", d):
        raise AssertionError(f"{d}")
