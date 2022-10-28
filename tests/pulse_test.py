# -*- coding: utf-8 -*-
u"""Tests for instantiation of LaserCavity, LaserPulse
and LaserPulseSlice
"""
from __future__ import absolute_import, division, print_function
import math
import numpy as np
from pykern.pkdebug import pkdp, pkdlog
from pykern.pkcollections import PKDict
import pykern.pkunit
import pytest
from rslaser.pulse import pulse
from rslaser.cavity import laser_cavity
import scipy


def pulse_instantiation_test(pulse, field):
    for s in pulse.slice:
        if getattr(s, field) != getattr(pulse, field):
            pykern.pkunit.pkfail(f'LaserPulseSlice has different {field} than pulse as a whole')


def test_no_params_instantiation():
    l = pulse.LaserPulse()


def slice_instantiation_test(pulse, field):
    a = [getattr(s, field) for s in pulse.slice]

    # TODO (gurhar1133): do we still want to do this?
    if not len(set(a)) == len(a):
        pykern.pkunit.pkfail(f'Expected pulse.slice field {field} to all be unique vals')


def test_instantiation01():
    l = pulse.LaserPulse()
    pulse_instantiation_test(l, 'phE')
    pulse_instantiation_test(l, '_lambda0')
    k = PKDict(chirp=0.01*(scipy.constants.h * scipy.constants.c / scipy.constants.e / 1e-6))
    l = pulse.LaserPulse(k)
    slice_instantiation_test(l, '_lambda0')

def test_instantiation02():
    k = PKDict(chirp=0.01*(scipy.constants.h * scipy.constants.c / scipy.constants.e / 1e-6))
    l = pulse.LaserPulse(k)
    slice_instantiation_test(l, 'phE')
    s = [pulse.LaserPulseSlice(i) for i in range(10)]
    c = [[i, 0] for i in range(10)]
    for a in c:
        with pykern.pkunit.pkexcept(pulse.InvalidLaserPulseInputError):
            pulse.LaserPulseSlice(a)


def test_instantiation03():
    p = PKDict(PHHe=0.1)
    with pykern.pkunit.pkexcept(pulse.InvalidLaserPulseInputError):
        pulse.LaserPulse(p)
    p = PKDict(slice_params=PKDict(**pulse._LASER_PULSE_SLICE_DEFAULTS))
    p.slice_params.update(PKDict(blonk=9))
    with pykern.pkunit.pkexcept(pulse.InvalidLaserPulseInputError):
        pulse.LaserPulse(p)


def test_instantiation04():
    e = "'PKDict' object has no attribute 'sigx_waist'"
    with pykern.pkunit.pkexcept(pulse.InvalidLaserPulseInputError):
        pulse.LaserPulse([])
    p = PKDict(slice_params=PKDict(**pulse._LASER_PULSE_SLICE_DEFAULTS))
    p.slice_params = PKDict(
        slice_params=PKDict(
            foo='bar', hello='world'
        )
    )
    with pykern.pkunit.pkexcept(e):
        pulse.LaserPulse(p)


def test_instantiation05():
    laser_cavity.LaserCavity()

    k = PKDict(n3='fail')
    with pykern.pkunit.pkexcept(laser_cavity.InvalidLaserCavityInputError):
        laser_cavity.LaserCavity(k)

    k = PKDict(pulse_params=PKDict(chirp=0.01*(scipy.constants.h * scipy.constants.c / 1e-6)))
    laser_cavity.LaserCavity(k)

    k = PKDict(
        pulse_params=PKDict(
            chirp=0.01*(scipy.constants.h * scipy.constants.c / 1e-6),
            slice_params=PKDict(
                **pulse._LASER_PULSE_SLICE_DEFAULTS
            )
        )
    )
    k.pulse_params.slice_params.update(PKDict(poltype=5))
    laser_cavity.LaserCavity(k)
    k = PKDict()
    laser_cavity.LaserCavity(k)
    with pykern.pkunit.pkexcept(
        laser_cavity.InvalidLaserCavityInputError,
        'Invalid instantiation of LaserCavity with in="should raise"'
        ):
        laser_cavity.LaserCavity('should raise')


def test_cavity_propagation():
    L_cav = 8
    dfL = 1
    dfR = 1

    L_cryst = 2 * 1e-2
    n0 = 2
    n2 = 0.02

    wavefrontEnergy = 1.55
    lam = scipy.constants.c * scipy.constants.value('Planck constant in eV/Hz') / wavefrontEnergy

    L_eff = L_cav + (1/n0 - 1) * L_cryst
    beta0 = math.sqrt(L_eff * (L_cav / 4 + dfL) - L_eff ** 2 / 4)
    sigx0 = math.sqrt(lam * beta0 / 4 / math.pi)

    lc = laser_cavity.LaserCavity(PKDict(
        drift_right_length=L_cav / 2 - L_cryst / 2,
        drift_left_length=L_cav / 2 - L_cryst / 2,
        lens_left_focal_length=L_cav / 4 + dfR,
        lens_right_focal_length=L_cav / 4 + dfL,
        n0=n0,
        n2=n2,
        L_half_cryst=L_cryst / 2,
        pulse_params=PKDict(
            phE=wavefrontEnergy,
            nslice=11,
            slice_params=PKDict(**pulse._LASER_PULSE_SLICE_DEFAULTS),
        ),
    ))

    results = []

    def intensity_callback(position, vals):
        (x, y) = lc.laser_pulse.rmsvals()
        results.append([
            lc.laser_pulse.pulsePos(),
            x,
            y,
            lc.laser_pulse.intensity_vals(),
            lc.laser_pulse.energyvals(),
        ])

    lc.propagate(num_cycles=4, callback=intensity_callback)
    print(results[-1])

    # TODO (gurhar1133): make this an ndiff
    pykern.pkunit.pkeq(
        str([[-0.09007729696643918, -0.07369960660890479, -0.05732191625137039, -0.04094422589383599, -0.0245665355363016, -0.00818884517876721, 0.008188845178767196, 0.024566535536301587, 0.04094422589383598, 0.05732191625137037, 0.07369960660890476], [0.001267401706551219, 0.0010762227299791487, 0.0008145963753205132, 0.0007290811958715039, 0.00046594193873525304, 0.00018620557753150257, 0.00018483984936330976, 0.0005109487588028025, 0.000650165346892296, 0.0009435856089303615, 0.0010309987887040782], [0.0012674016346121467, 0.0010762227439608087, 0.0008145964242898868, 0.0007290812090433151, 0.0004659419224893989, 0.00018620555092352633, 0.00018483983760892106, 0.0005109488069060423, 0.0006501653585135822, 0.0009435856849413944, 0.001030998795197547], [3360553700000.0, 6520655300000.0, 22691694000000.0, 100779076000000.0, 271603130000000.0, 848986700000000.0, 1652040400000000.0, 163319100000000.0, 76051780000000.0, 16465634000000.0, 49943532000000.0], [1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55]]),
        str(results[-1]),
    )
