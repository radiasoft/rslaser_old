# -*- coding: utf-8 -*-
u"""Definition of a laser pulse
Copyright (c) 2021 RadiaSoft LLC. All rights reserved
"""

import math
import numpy as np
import scipy.constants as const

from pykern.pkcollections import PKDict
import rslaser.rsoptics.wavefront as rswf
import rslaser.rspulse.gauss_hermite as rsgh
import rslaser.utils.constants as rsc
import rslaser.utils.unit_conversion as units
# import rslaser.utils.srwl_uti_data as rsdata

import srwlib
from srwlib import srwl


def _validate_input(input_params):
    if type(input_params) != PKDict:
        raise InvalidLaserPulseInputError('invalid LaserPulse inputs: input_parameters to LaserPulse class must be of type PKDict')
    if 'slice_params' not in input_params:
        raise InvalidLaserPulseInputError('invalid LaserPulse inputs: does not include slice_params')


class InvalidLaserPulseInputError(Exception):
    pass


class LaserPulse:
    """
    The LaserPulse contains a GaussHermite object to represent the initial envelope,
    as well as an array of LaserPulseSlice instances, which track details of the evolution in time.

    """
    def __init__(self, params):

        _validate_input(params)
        # instantiate the laser envelope
        self.envelope = rsgh.GaussHermite(params)

        # instantiate the array of slices
        self.slice = []
        self.nslice = params.nslice
        self.phE = params.phE
        self._lambda0 = abs(units.calculate_lambda0_from_phE(params.phE))
        self.phE -= 0.5*params.chirp           # so central slice has the central photon energy
        _de = params.chirp / self.nslice   # photon energy shift from slice to slice
        s = params
        s.phE = self.phE
        for i in range(params.nslice):
            # add the slices; each (slowly) instantiates an SRW wavefront object
            self.slice.append(LaserPulseSlice(i, s))
            s.phE += _de
        self._sxvals = []  # horizontal slice data
        self._syvals = []  # vertical slice data

    def compute_middle_slice_intensity(self):
        wfr = self.slice[len(self.slice) // 2].wfr
        (ar2d, sx, sy, xavg, yavg) = rswf.rmsWavefrontIntensity(wfr)
        self._sxvals.append(sx)
        self._syvals.append(sy)
        return (wfr, sx, sy)

    def rmsvals(self):
        sx = []
        sy = []
        for sl in self.slice:
            (_, sigx,sigy, _, _) = rswf.rmsWavefrontIntensity(sl.wfr)
            sx.append(sigx)
            sy.append(sigy)

        return(sx,sy)

    def intensity_vals(self):
        return [rswf.maxWavefrontIntensity(s.wfr) for s in self.slice]

    def pulsePos(self):
        return [s._pulse_pos for s in self.slice]

    def energyvals(self):
        return [s.phE for s in self.slice]

    def slice_wfr(self,slice_index):
        return self.slice(slice_index).wfr


class LaserPulseSlice:
    """
    This class represents a longitudinal slice in a laser pulse.
    There will be a number of wavefronts each with different wavelengths (energy).
    The slice is composed of an SRW wavefront object, which is defined here:
    https://github.com/ochubar/SRW/blob/master/env/work/srw_python/srwlib.py#L2048
    """
    def __init__(self, slice_index, params):
        """
        #nslice: number of slices of laser pulse
        #slice_index: index of slice
        #d_to_w: distance from the pulse center at t = 0 to the intended waist location [m]
        #sigrW: beam size at waist [m]
        #propLen: propagation length [m] required by SRW to create numerical Gaussian
        #propLen=15,
        #sig_s RMS pulse length [m]
        #pulseE: energy per pulse [J]
        #poltype: polarization type (0=linear horizontal, 1=linear vertical, 2=linear 45 deg, 3=linear 135 deg, 4=circular right, 5=circular left, 6=total)
        #phE: photon energy [eV]
        #sampFact: sampling factor to increase mesh density
        """
        #print([sigrW,propLen,pulseE,poltype])
        _validate_input(params)
        self._lambda0 = units.calculate_lambda0_from_phE(params.phE)
        self.slice_index = slice_index
        self.phE = params.phE
        constConvRad = 1.23984186e-06/(4*3.1415926536)  ##conversion from energy to 1/wavelength
        rmsAngDiv = constConvRad/(self.phE*params.slice_params.sigrW)             ##RMS angular divergence [rad]
        #  if at t=0 distance to the waist location d_to_w < d_to_w_cutoff, initialization in SRW involves/requires propagation
        #  from the distance-to-waist > d_to_w_cutoff to the actual z(t=0) for which d_to_w < d_to_w_cutoff
        d_to_w_cutoff = 0.001  # [m] - verify that this is a reasonable value
        if params.d_to_w > d_to_w_cutoff:  params.slice_params.propLen = params.d_to_w  #  d_to_w = L_d1 +0.5*L_c in the single-pass example
        sigrL=math.sqrt(params.slice_params.sigrW**2+(params.slice_params.propLen*rmsAngDiv)**2)  ##required RMS size to produce requested RMS beam size after propagation by propLen


        #***********Gaussian Beam Source
        GsnBm = srwlib.SRWLGsnBm() #Gaussian Beam structure (just parameters)
        GsnBm.x = 0 #Transverse Positions of Gaussian Beam Center at Waist [m]
        GsnBm.y = 0
        numsig = 3. #Number of sigma values to track. Total range is 2*numsig*sig_s
        ds = 2*numsig*params.slice_params.sig_s/(params.nslice - 1)
        self._pulse_pos = -numsig*params.slice_params.sig_s+slice_index*ds
        GsnBm.z = params.slice_params.propLen + self._pulse_pos #Longitudinal Position of Waist [m]
        GsnBm.xp = 0 #Average Angles of Gaussian Beam at Waist [rad]
        GsnBm.yp = 0
        GsnBm.avgPhotEn = self.phE #Photon Energy [eV]
        GsnBm.pulseEn = params.slice_params.pulseE*np.exp(-self._pulse_pos**2/(2*params.slice_params.sig_s**2)) #Energy per Pulse [J] - to be corrected
        GsnBm.repRate = 1 #Rep. Rate [Hz] - to be corrected
        GsnBm.polar = params.slice_params.poltype #1- linear horizontal?
        GsnBm.sigX = params.slice_params.sigrW #Horiz. RMS size at Waist [m]
        GsnBm.sigY = GsnBm.sigX #Vert. RMS size at Waist [m]

        GsnBm.sigT = 10e-15 #Pulse duration [s] (not used?)
        GsnBm.mx = params.slice_params.mx #Transverse Gauss-Hermite Mode Orders
        GsnBm.my = params.slice_params.my

        #***********Initial Wavefront
        _wfr = srwlib.SRWLWfr() #Initial Electric Field Wavefront
        _wfr.allocate(1, 1000, 1000) #Numbers of points vs Photon Energy (1), Horizontal and Vertical Positions (dummy)
        _wfr.mesh.zStart = 0.0 #Longitudinal Position [m] at which initial Electric Field has to be calculated, i.e. the position of the first optical element
        _wfr.mesh.eStart = GsnBm.avgPhotEn #Initial Photon Energy [eV]
        _wfr.mesh.eFin = GsnBm.avgPhotEn #Final Photon Energy [eV]

        _wfr.unitElFld = 1 #Electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)

        distSrc = _wfr.mesh.zStart - GsnBm.z
        #Horizontal and Vertical Position Range for the Initial Wavefront calculation
        #can be used to simulate the First Aperture (of M1)
        #firstHorAp = 8.*rmsAngDiv*distSrc #[m]
        xAp = 8.*sigrL
        yAp = xAp #[m]

        _wfr.mesh.xStart = -0.5*xAp #Initial Horizontal Position [m]
        _wfr.mesh.xFin = 0.5*xAp #Final Horizontal Position [m]
        _wfr.mesh.yStart = -0.5*yAp #Initial Vertical Position [m]
        _wfr.mesh.yFin = 0.5*yAp #Final Vertical Position [m]

        sampFactNxNyForProp = params.slice_params.sampFact #sampling factor for adjusting nx, ny (effective if > 0)
        arPrecPar = [sampFactNxNyForProp]

        srwlib.srwl.CalcElecFieldGaussian(_wfr, GsnBm, arPrecPar)

        ##Beamline to propagate to waist ( only if d_to_w(t=0) < d_to_w_cutoff )
        if params.d_to_w < d_to_w_cutoff:
          optDriftW=srwlib.SRWLOptD(params.slice_params.propLen)
          propagParDrift = [0, 0, 1., 0, 0, 1.1, 1.2, 1.1, 1.2, 0, 0, 0]
          optBLW = srwlib.SRWLOptC([optDriftW],[propagParDrift])
          srwlib.srwl.PropagElecField(_wfr, optBLW)
        self.wfr = _wfr

    def validate_input(self, input_params):
        pass
