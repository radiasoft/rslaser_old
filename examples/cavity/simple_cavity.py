# -*- coding: utf-8 -*-
u"""Propagation of a wavefront in an optical cavity

The initial wavefront is a Hermite-Gaussian mode located at the center of a cavity
of length L_cav. There is a crystal located at the center of the cavity of length L_cryst.
The wavefront is propagated through the crystal and then to the end of the cavity.

:copyright: Copyright (c) 2016 RadiaSoft LLC.  All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""

from __future__ import absolute_import, division, print_function

import math
import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np

from pykern.pkcollections import PKDict
from pykern import pkdebug

import rslaser.cavity.laser_cavity as rslc
import rslaser.optics.element as rse
import rsmath.const as rsc

import scipy.constants as const
import time

# This is necessary on the JupyterHub commandline
# Delete elsewhere to automatically use the local defaul GUI backend
mpl.use("agg")

# This is helpful for dealing with agg vs GUI backends
def agg_plot(filename):
    if mpl.get_backend() == 'agg':
        plt.savefig(filename)
        plt.close()
    else:
        plt.show()

print(" ")
print("Define laser pulse, crystal and cavity parameters...")

# specify all default values of system parameters in a Pykern dictionary object

# photon energy corresponding to lambda = 1 micron
_PHE_DEFAULT = const.h * const.c / 1e-6
_Z_WAIST_DEFAULT = 0
_Z_CENTER_DEFAULT = 0

_LASER_PULSE_SLICE_DEFAULTS = PKDict(
    sigrW=0.000186,
    propLen=15,
    pulseE=0.001,
    poltype=1,
    sampFact=5,
    numsig=3.,
    mx=0,
    my=0
)
_LASER_PULSE_DEFAULTS = PKDict(
        phE=_PHE_DEFAULT,
        nslice=3,
        chirp=0,
        w0=.1,
        a0=.01,
        dw0x=0.0,
        dw0y=0.0,
        z_waist=_Z_WAIST_DEFAULT,
        dzwx=0.0,
        dzwy=0.0,
        tau_fwhm=0.1 / const.c / math.sqrt(2.),
        z_center=_Z_CENTER_DEFAULT,
        x_shift = 0.,
        y_shift=0.,
        d_to_w=_Z_WAIST_DEFAULT - _Z_CENTER_DEFAULT,
        slice_params=_LASER_PULSE_SLICE_DEFAULTS,
)
_LASER_CAVITY_DEFAULTS = PKDict(
    drift_right_length=0.5,
    drift_left_length=0.5,
    lens_left_focal_length=0.2,
    lens_right_focal_length=0.2,
    n0 = 1.75,
    n2 = 0.001,
    L_half_cryst=0.2,
)

# the laser pulse parameters are loaded into a Pykern dictionary object
# the slice parameters are included hierarchically
in_pulse = PKDict(
    **_LASER_PULSE_DEFAULTS
    )

# Specify non-default values for laser pulse parameters at center of crystal
in_pulse.nslice = 5

# the laser cavity parameters are loaded into a Pykern dictionary object
# the laser pulse parameters are included hierarchically
in_cavity = PKDict(
    **_LASER_CAVITY_DEFAULTS,
    pulse_params=in_pulse,
    )

# Specify non-default values for laser cavity
# Also, specify/print other values for consideration by the user

# cavity length [m]
L_cav = 1.

# Define right and left mirror focal lengths
# Focal length difference from confocal case [m]
df = 0.3
dfR = df
dfL = df

# focal length
f = L_cav/4. + df

in_cavity.lens_left_focal_length = f
in_cavity.lens_right_focal_length = f

# Length of crystal [m]
L_cryst = 0.1
# index of refraction on axis
in_cavity.n0=1.75
# radial variation of index of refraction: n(r) = n0 - 0.5 n2 r^2
in_cavity.n2=0.0
in_cavity.drift_right_length = L_cav/2. - L_cryst/2.
in_cavity.drift_left_length = in_cavity.drift_right_length

# Define effective length as path length reduced by index of refraction n0
L_eff = L_cav + (1./in_cavity.n0 - 1.) * L_cryst
print("L_eff = %f m" % L_eff)

beta0 = np.sqrt(L_eff*f - L_eff**2/4.)
print("beta0 = %f m" % beta0)

lambda0 = const.h * const.c / in_pulse.phE
sigx0 = np.sqrt(lambda0*beta0/4./np.pi)
print("sigx0 = %f m" % sigx0)

sigrW = sigx0
w0 = sigx0 * math.sqrt(2.)
in_cavity.L_half_cryst = L_cryst/2.

print(" ")
print("Initialize the laser cavity object...")
LC = rslc.LaserCavity(in_cavity)

svals = LC.laser_pulse.pulsePos()
(lpsxvals,lpsyvals) = LC.laser_pulse.rmsvals()
ivals = LC.laser_pulse.intensity_vals()
evals = LC.laser_pulse.energyvals()

print(" ")
print("Plot the RMS values along the initial laser pulse...")

plt.title('RMS Beam size along laser pulse')
plt.plot(svals,lpsyvals)
plt.xlabel('s [m]')
plt.ylabel('rms x [m]')

agg_plot('img01.png')

print(" ")
print("Plot the pulse intensity along the initial laser pulse...")

plt.plot(svals,ivals)
plt.ylabel('pulse intensity []')
plt.xlabel('s [m]')

agg_plot('img02.png')

print(" ")
print("Plot the slice energy along the initial laser pulse...")

plt.title('Slice energy along laser pulse')
plt.plot(svals,evals)
plt.xlabel('s [m]')
plt.ylabel('Energy [ev]')

agg_plot('img03.png')

print(" ")
print("Propagate the laser pulse several times through the cavity....")

start_time = time.time()
(svals, sxvals, syvals) = LC.propagate(num_cycles=4)

print("Simulation time: %s seconds" %(round((time.time() - start_time),5)))

print(" ")
print("Plot the transverse size of the pulse as it oscillates back and forth ...")

fig, ax = plt.subplots()
ax.plot(svals, sxvals)

agg_plot('img04.png')

print(" ")
print("Plot the RMS values along the final laser pulse...")

svals = LC.laser_pulse.pulsePos()
(lpsxvals, lpsyvals) = LC.laser_pulse.rmsvals()
ivals = LC.laser_pulse.intensity_vals()
evals = LC.laser_pulse.energyvals()

plt.title('RMS Beam size along laser pulse')
plt.plot(svals, lpsyvals)
plt.xlabel('s [m]')
plt.ylabel('rms x [m]')

agg_plot('img05.png')

print(" ")
print("Plot the laser intensity along the final laser pulse...")
plt.title('Intensity along laser pulse')
plt.plot(svals,ivals)
plt.ylabel('pulse intensity []')
plt.xlabel('s [m]')

agg_plot('img06.png')

print(" ")
print("Plot the slice energy along the final laser pulse...")
plt.title('Slice energy along laser pulse')
plt.plot(svals, evals)
plt.xlabel('s [m]')
plt.ylabel('Energy [ev]')

agg_plot('img07.png')
