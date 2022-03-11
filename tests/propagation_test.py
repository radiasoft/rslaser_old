# -*- coding: utf-8 -*-
u"""Test LaserCavity propagation
"""
from pykern.pkcollections import PKDict
from pykern.pkdebug import pkdp, pkdlog
from rslaser.rscavity.laser_cavity import LaserCavity
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
            [-0.09007729696643918, -0.07206183757315135, -0.054046378179863505, -0.03603091878657567, -0.01801545939328783, 0.0, 0.018015459393287844, 0.03603091878657569, 0.05404637817986352, 0.07206183757315135, 0.09007729696643918],
            [0.0004206083571113905, 0.0004215773580955795, 0.000422558801457325, 0.0004235528433050065, 0.00042455921404788854, 0.00042557776027345, 0.00042660872344884643, 0.0004276516539122094, 0.00042870669687600845, 0.0004297736312473442, 0.0004308525472033005],
            [0.0004206082780286289, 0.00042157735198943637, 0.0004225588410242694, 0.0004235528647439696, 0.0004245591094213043, 0.0004255777869700445, 0.00042660868997989037, 0.0004276516954501241, 0.0004287067054042391, 0.0004297736687977002, 0.0004308525513948011],
            [768161600000.0, 3864636400000.0, 13564228000000.0, 33213339000000.0, 56736250000000.0, 67614554000000.0, 56214814000000.0, 32605664000000.0, 13193718000000.0, 3724545000000.0, 733520100000.0],
            [1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55],
        ]),
        str(results[-1]),
    )
