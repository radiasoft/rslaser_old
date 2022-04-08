# -*- coding: utf-8 -*-
u"""Tests for LaserPulseEnvelope
"""
from __future__ import absolute_import, division, print_function
import math
import numpy as np
from pykern.pkdebug import pkdp, pkdlog
from pykern.pkcollections import PKDict
from pykern.pkunit import pkexcept
import pytest
from rslaser.pulse import pulse
from rslaser.cavity import laser_cavity
import scipy.constants as const


def test_envelope():
    e = pulse.LaserPulseEnvelope()


def test_envelope2():
    with pkexcept(pulse.InvalidLaserPulseInputError):
        pulse.LaserPulseEnvelope(1)


def test_envelope3():
    with pkexcept(pulse.InvalidLaserPulseInputError):
        pulse.LaserPulseEnvelope(PKDict(test='test'))


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
    with pkexcept(TypeError):
        e.evaluate_envelope_ex('should fail')
