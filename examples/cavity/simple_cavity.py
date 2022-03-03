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

import rslaser.rscavity.laser_cavity as rslc
import rslaser.rsoptics.element as rse
import rslaser.utils.constants as rsc

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

# specify the parameters in a PKDict dictionary object
_PHE_DEFAULT = const.h * const.c / 1e-6   # photon energy corresponding to lambda = 1 micron
_Z_WAIST_DEFAULT = 0
_Z_CENTER_DEFAULT = 0
_LASER_PULSE_SLICE_DEFAULTS = PKDict(
    sigrW=0.000186,
    propLen=15,
    sig_s=0.1,
    pulseE=0.001,
    poltype=1,
    sampFact=5,
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
        z_center=_Z_CENTER_DEFAULT,
        x_shift = 0.,
        y_shift=0.,
        d_to_w=_Z_WAIST_DEFAULT - _Z_CENTER_DEFAULT,
        slice_params=_LASER_PULSE_SLICE_DEFAULTS,
)

# central laser pulse values
a0 = 0.01             # 0.85e-9 * lambda [microns] * Sqrt(Intensity [W/cm^2])
lambda0 = 8.e-7       # wavelength [m]

# numerical values
num_slices = 5        # desired number of slices (i.e. SRW wavefronts) to represent the pulse

# local Pykern dictionary object to hold all physical parameters
k=_LASER_PULSE_DEFAULTS.copy()
k.L_cav = 1 #Length of cavity [m]
k.df = 0.3 #Focal length difference from confocal case [m]
k.phE = const.h * const.c / lambda0
k.nslice = num_slices

# Define right and left mirror focal lengths
k.dfR = k.df
k.dfL = k.df
k.f=k.L_cav/4+k.df #focal length

k.lens_left_focal_length = k.f
k.lens_right_focal_length = k.f

#Define Crystal parameters
k.L_cryst = 0.1   # Length of crystal [m]
k.n0=1.75         # index of refraction on axis
k.n2=0.0          # radial variation of index of refraction: n(r) = n0 - 0.5 n2 r^2
k.drift_right_length=k.L_cav/2-k.L_cryst/2
k.drift_left_length=k.drift_right_length

k.L_eff = k.L_cav+(1/k.n0 - 1)*k.L_cryst #Define effective length as path length reduced by index of refraction n0
print("L_eff = %f m" % k.L_eff)
k.beta0 = np.sqrt(k.L_eff*k.f-k.L_eff**2/4)
print("beta0 = %f m" % k.beta0)
k.sigx0 = np.sqrt(lambda0*k.beta0/4/np.pi)
print("sigx0 = %f m" % k.sigx0)
k.sigrW = k.sigx0
k.w0 = k.sigx0 * math.sqrt(2.)

k.sig_s = 0.1  # rms length of Gaussian laser pulse [m]
tau_fwhm = k.sig_s / const.c / math.sqrt(2.)
k.tau_fwhm = tau_fwhm

k.L_half_cryst=k.L_cryst/2

print(" ")
print("Initialize the laser cavity object...")
LC = rslc.LaserCavity(k)

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
