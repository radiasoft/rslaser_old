# -*- coding: utf-8 -*-
u"""Tests for LaserPulseEnvelope
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
import test_utils

def test_envelope():
    e = pulse.LaserPulseEnvelope()


def test_envelope2():
    test_utils.trigger_exception_test(pulse.LaserPulseEnvelope, 1)


def test_envelope3():
    test_utils.trigger_exception_test(
        pulse.LaserPulseEnvelope,
        PKDict(
            test='test'
        )
    )


def test_envelope4():
    e = pulse.LaserPulseEnvelope()
    e.evaluate_envelope_ex(
        np.random.rand(12),
        np.random.rand(12),
        0.1)


def test_envelope5():
    e = pulse.LaserPulseEnvelope()
    e.evaluate_envelope_ex(1, 2, 3)


def test_envelope6():
    e = pulse.LaserPulseEnvelope()
    test_utils.trigger_exception_test(e.evaluate_envelope_ex, 'should fail')




