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
            [-0.09007729696643918, -0.07369960660890479, -0.05732191625137039, -0.04094422589383599,
             -0.0245665355363016, -0.00818884517876721, 0.008188845178767196, 0.024566535536301587,
              0.04094422589383598, 0.05732191625137037, 0.07369960660890476], [0.0004206083571113905,
              0.00042148866696306385, 0.0004223794433054218, 0.0004232805397432438, 0.0004241917784428871,
              0.000425113299639614, 0.00042604485049852994, 0.0004269865220384565, 0.0004279381999672457,
              0.00042889980587495087, 0.0004298713099470023], [0.0004206082780286289, 0.00042148868408716775,
              0.00042237944785462596, 0.00042328054426942787, 0.00042419175847993286, 0.0004251132953753882,
               0.00042604485714848017, 0.0004269865118091315, 0.0004279381772178525, 0.00042889976698366305,
               0.00042987128270201806], [768161600000.0, 3386752600000.0, 11088761000000.0, 26961980000000.0,
               48684415000000.0, 65282614000000.0, 65009240000000.0, 48075314000000.0, 26402127000000.0,
               10767783000000.0, 3261246500000.0], [1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55]
               ]),
        str(results[-1]),
    )
