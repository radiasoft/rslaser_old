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
    # TODO (gurhar1133): looks like the validator checks that
    # things that arent supposed to be there are not there.
    # But does not trigger InvalidLaserPulseInputError when
    # there are missing params that are essential.
    # check if this is happening with both LaserPulse and
    # LaserPulseSlice validators. May want to reconsider
    # validation behavior

    p = PKDict(PHHe=0.1)
    with pykern.pkunit.pkexcept(pulse.InvalidLaserPulseInputError):
        pulse.LaserPulse(p)
    p = PKDict(slice_params=PKDict(blonk=9))
    with pykern.pkunit.pkexcept(pulse.InvalidLaserPulseInputError):
        pulse.LaserPulse(p)
    p = PKDict(slice_params=PKDict(pulseE=9))
    pulse.LaserPulse(p)
    p = PKDict(slice_params=PKDict())
    pulse.LaserPulse(p)


def test_instantiation04():
    with pykern.pkunit.pkexcept(pulse.InvalidLaserPulseInputError):
        pulse.LaserPulse([])
    k = PKDict(slice_params=[])
    with pykern.pkunit.pkexcept(pulse.InvalidLaserPulseInputError):
        pulse.LaserPulse(k)
    k = PKDict(
        slice_params=PKDict(
            foo='bar', hello='world'
        )
    )
    with pykern.pkunit.pkexcept(pulse.InvalidLaserPulseInputError):
        pulse.LaserPulse(k)
    with pykern.pkunit.pkexcept(pulse.InvalidLaserPulseInputError):
        pulse.LaserPulseSlice('10')


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
                propLen=2,
            )
        )
    )
    laser_cavity.LaserCavity(k)
    k = PKDict(
        pulse_params=PKDict(
            slice_params=PKDict(
                propLen=5,
            )
        )
    )
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
            slice_params=PKDict(
                sigrW=sigx0,
                propLen=15,
                #sig_s=0.1,
                pulseE=0.001,
                poltype=1,
                sampFact=1,
            ),
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
    pykern.pkunit.pkeq(
        str([
            [-0.09007729696643918, -0.07369960660890479, -0.05732191625137039, -0.04094422589383599, -0.0245665355363016, -0.00818884517876721, 0.008188845178767196, 0.024566535536301587, 0.04094422589383598, 0.05732191625137037, 0.07369960660890476],
            [0.00042060832381247547, 0.00042148869279589853, 0.00042237938068257706, 0.00042328053915324966, 0.0004241917709675302, 0.000425113259460114, 0.0004260448258076334, 0.0004269864714645692, 0.00042793823411483965, 0.0004288998097451079, 0.00042987124363644844],
            [0.000420608251804938, 0.00042148867911241705, 0.00042237936420157245, 0.00042328052228606146, 0.0004241917082775849, 0.0004251133507806255, 0.00042604486544507405, 0.00042698643452929934, 0.00042793819893895553, 0.00042889977737009993, 0.00042987124878598516],
            [768158600000.0, 3386736400000.0, 11088713000000.0, 26961888000000.0, 48684220000000.0, 65282380000000.0, 65008927000000.0, 48075112000000.0, 26401999000000.0, 10767730000000.0, 3261234200000.0],
            [1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55]]),
        str(results[-1]),
    )
