# -*- coding: utf-8 -*-
u"""Test LaserCavity propagation
"""
from pykern.pkcollections import PKDict
from pykern.pkdebug import pkdp, pkdlog
from rslaser.cavity.laser_cavity import LaserCavity
import math
import scipy.constants
from pykern.pkunit import pkeq, pkexcept, pkok, pkre


def test_propagation():
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

    lc = LaserCavity(PKDict(
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
    pkeq(
        str([
            [-0.09007729696643918, -0.07369960660890479, -0.05732191625137039, -0.04094422589383599, -0.0245665355363016, -0.00818884517876721, 0.008188845178767196, 0.024566535536301587, 0.04094422589383598, 0.05732191625137037, 0.07369960660890476],
            [0.00042060832381247547, 0.00042148869279589853, 0.00042237938068257706, 0.00042328053915324966, 0.0004241917709675302, 0.000425113259460114, 0.0004260448258076334, 0.0004269864714645692, 0.00042793823411483965, 0.0004288998097451079, 0.00042987124363644844],
            [0.000420608251804938, 0.00042148867911241705, 0.00042237936420157245, 0.00042328052228606146, 0.0004241917082775849, 0.0004251133507806255, 0.00042604486544507405, 0.00042698643452929934, 0.00042793819893895553, 0.00042889977737009993, 0.00042987124878598516],
            [768158600000.0, 3386736400000.0, 11088713000000.0, 26961888000000.0, 48684220000000.0, 65282380000000.0, 65008927000000.0, 48075112000000.0, 26401999000000.0, 10767730000000.0, 3261234200000.0],
            [1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55]]),
        str(results[-1]),
    )
