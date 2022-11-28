# -*- coding: utf-8 -*-
"""Tests for instantiation of LaserCavity, LaserPulse
and LaserPulseSlice
"""
from __future__ import absolute_import, division, print_function
import math
import numpy as np
from pykern.pkdebug import pkdp, pkdlog
from pykern.pkcollections import PKDict
from wavefront_sensor_test import ndiff_files
import pykern.pkunit
import pytest
from rslaser.pulse import pulse
from rslaser.cavity import laser_cavity
import scipy
import rslaser
import srwlib


_PACKAGE_DATA_DIR = rslaser.pkg_resources.resource_filename("rslaser", "package_data")


def pulse_instantiation_test(pulse, field):
    for s in pulse.slice:
        if getattr(s, field) != getattr(pulse, field):
            pykern.pkunit.pkfail(
                f"LaserPulseSlice has different {field} than pulse as a whole"
            )


def test_instantiation01():
    pulse.LaserPulse()


def test_instantiation02():
    p = PKDict(PHHe=0.1)
    with pykern.pkunit.pkexcept(pulse.InvalidLaserPulseInputError):
        pulse.LaserPulse(p)
    p = PKDict(slice_params=PKDict(**pulse._LASER_PULSE_SLICE_DEFAULTS))
    p.slice_params.update(PKDict(blonk=9))
    with pykern.pkunit.pkexcept(pulse.InvalidLaserPulseInputError):
        pulse.LaserPulse(p)


def test_instantiation03():
    e = "'PKDict' object has no attribute 'sigx_waist'"
    with pykern.pkunit.pkexcept(pulse.InvalidLaserPulseInputError):
        pulse.LaserPulse([])
    p = PKDict(slice_params=PKDict(**pulse._LASER_PULSE_SLICE_DEFAULTS))
    p.slice_params = PKDict(slice_params=PKDict(foo="bar", hello="world"))
    with pykern.pkunit.pkexcept(e):
        pulse.LaserPulse(p)


def test_instantiation04():
    laser_cavity.LaserCavity()

    k = PKDict(n3="fail")
    with pykern.pkunit.pkexcept(laser_cavity.InvalidLaserCavityInputError):
        laser_cavity.LaserCavity(k)

    k = PKDict(
        pulse_params=PKDict(chirp=0.01 * (scipy.constants.h * scipy.constants.c / 1e-6))
    )
    laser_cavity.LaserCavity(k)

    k = PKDict(
        pulse_params=PKDict(
            chirp=0.01 * (scipy.constants.h * scipy.constants.c / 1e-6),
            slice_params=PKDict(**pulse._LASER_PULSE_SLICE_DEFAULTS),
        )
    )
    k.pulse_params.slice_params.update(PKDict(poltype=5))
    laser_cavity.LaserCavity(k)
    k = PKDict()
    laser_cavity.LaserCavity(k)
    with pykern.pkunit.pkexcept(
        laser_cavity.InvalidLaserCavityInputError,
        'Invalid instantiation of LaserCavity with in="should raise"',
    ):
        laser_cavity.LaserCavity("should raise")

# TODO (gurhar1133): propagation is a work in progress.
# def test_cavity_propagation():
#     from pykern import pkunit
#     from pykern import pkio

#     data_dir = pkunit.data_dir()
#     work_dir = pkunit.empty_work_dir()
#     L_cav = 8
#     dfL = 1
#     dfR = 1

#     L_cryst = 2 * 1e-2
#     n0 = 2
#     n2 = 0.02

#     wavefrontEnergy = 1.55
#     lam = (
#         scipy.constants.c
#         * scipy.constants.value("Planck constant in eV/Hz")
#         / wavefrontEnergy
#     )

#     L_eff = L_cav + (1 / n0 - 1) * L_cryst
#     beta0 = math.sqrt(L_eff * (L_cav / 4 + dfL) - L_eff**2 / 4)
#     sigx0 = math.sqrt(lam * beta0 / 4 / math.pi)

#     lc = laser_cavity.LaserCavity(
#         PKDict(
#             drift_right_length=L_cav / 2 - L_cryst / 2,
#             drift_left_length=L_cav / 2 - L_cryst / 2,
#             lens_left_focal_length=L_cav / 4 + dfR,
#             lens_right_focal_length=L_cav / 4 + dfL,
#             n0=n0,
#             n2=n2,
#             L_half_cryst=L_cryst / 2,
#             pulse_params=PKDict(
#                 phE=wavefrontEnergy,
#                 nslice=11,
#                 slice_params=PKDict(**pulse._LASER_PULSE_SLICE_DEFAULTS),
#             ),
#         )
#     )

#     results = []

#     def intensity_callback(position, vals):
#         (x, y) = lc.laser_pulse.rmsvals()
#         results.append(
#             [
#                 lc.laser_pulse.pulsePos(),
#                 x,
#                 y,
#                 lc.laser_pulse.intensity_vals(),
#                 lc.laser_pulse.energyvals(),
#             ]
#         )

#     lc.propagate(num_cycles=4, callback=intensity_callback)

#     ndiff_files(
#         data_dir.join("res.txt"),
#         pkio.write_text(
#             work_dir.join("res_actual.txt"),
#             str(results[-1]),
#         ),
#         work_dir.join("ndiff.out"),
#         data_dir,
#     )


def test_from_file():
    from pykern import pkunit
    from pykern import pkio

    data_dir = pkunit.data_dir()
    work_dir = pkunit.empty_work_dir()
    pulse_inputs = pulse._LASER_PULSE_DEFAULTS.copy()
    pulse_inputs.nslice = 1
    f = PKDict(
        ccd=pkio.py_path(_PACKAGE_DATA_DIR).join("ccd_pump_off.txt"),
        wfs=pkio.py_path(_PACKAGE_DATA_DIR).join("wfs_pump_off.txt"),
        meta=pkio.py_path(_PACKAGE_DATA_DIR).join("wfs_meta.dat"),
    )
    wavefront = pulse.LaserPulse(
        pulse_inputs,
        files=f,
    ).slice_wfr(0)
    intensity = srwlib.array('f', [0]*wavefront.mesh.nx*wavefront.mesh.ny)
    srwlib.srwl.CalcIntFromElecField(intensity, wavefront, 6, 0, 3, wavefront.mesh.eStart, 0, 0)
    ndiff_files(
        data_dir.join("2d_wf_intensity.txt"),
        pkio.write_text(
            work_dir.join("2d_wf_intensity_actual.txt"),
            str(intensity),
        ),
        work_dir.join("ndiff.out"),
        data_dir,
    )
    pulse_inputs.nslice = 2
    with pykern.pkunit.pkexcept(
        pulse.InvalidLaserPulseInputError,
        "cannot use file inputs with more than one slice",
    ):
        pulse.LaserPulse(
            pulse_inputs,
            f,
        )
