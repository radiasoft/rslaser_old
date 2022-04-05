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


def test_wfs_instantiation():
    wfs = WavefrontSensor('w1', 2.0)


def test_propagate():
    p = pulse.LaserPulse()
    wfs = WavefrontSensor('w1', 2.0)
    wfs.propagate(p)