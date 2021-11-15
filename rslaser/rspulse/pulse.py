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
# import rslaser.utils.srwl_uti_data as rsdata

import srwlib
from srwlib import srwl

class LaserPulse:
    """
    The LaserPulse contains a GaussHermite object to represent the initial envelope,
    as well as an array of LaserPulseSlice instances, which track details of the evolution in time. 
    
    """
    def __init__(self,kwargs):
        _k = kwargs.copy()
        
        # instantiate the laser envelope
        self.envelope = rsgh.GaussHermite(_k)
        
        # instantiate the array of slices
        self.slice = []
        self.nslice = _k.nslice
        
        _lambda0 = abs(_k.lambda0)
        _phE = const.h * const.c / _lambda0
        _lambda_p = _lambda0 + 0.5 * _k.d_lambda
        _lambda_m = _lambda0 - 0.5 * _k.d_lambda
        _chirp = const.h * const.c * (1./_lambda_m - 1./_lambda_p)
        _phE -= 0.5*_chirp           # so central slice has the central photon energy
        _de = _chirp / self.nslice   # photon energy shift from slice to slice
        
        for i in range(_k.nslice):
            # add the slices; each (slowly) instantiates an SRW wavefront object
            self.slice.append(LaserPulseSlice(i,**_k))
            _phE += _de
        self.sxvals = []  # horizontal slice data
        self.syvals = []  # vertical slice data

    def compute_middle_slice_intensity(self):
        wfr = self.slice[len(self.slice) // 2]._wfr
        (ar2d, sx, sy, xavg, yavg) = rswf.rmsWavefrontIntensity(wfr)
        self.sxvals.append(sx)
        self.syvals.append(sy)
        return (wfr, sx, sy)

    def rmsvals(self):
        sx = []
        sy = []
        for sl in self._slice:
            (_, sigx,sigy, _, _) = rswf.rmsWavefrontIntensity(sl._wfr)
            sx.append(sigx)
            sy.append(sigy)

        return(sx,sy)

    def intensity_vals(self):
        return [rswf.maxWavefrontIntensity(s.wfr) for s in self.slice]

    def pulsePos(self):
        return [s.pulse_pos for s in self.slice]

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
    def __init__(self,slice_index,nslice,sigrW=0.00043698412731784714,propLen=15,sig_s=0.1,pulseE=0.001,poltype=1,phE=1.55,sampFact=5,mx=0,my=0,**_ignore_kwargs):
        """
        #nslice: number of slices of laser pulse
        #slice_index: index of slice
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
        self.slice_index = slice_index
        self.phE = phE
        constConvRad = 1.23984186e-06/(4*3.1415926536)  ##conversion from energy to 1/wavelength
        rmsAngDiv = constConvRad/(phE*sigrW)             ##RMS angular divergence [rad]
        sigrL=math.sqrt(sigrW**2+(propLen*rmsAngDiv)**2)  ##required RMS size to produce requested RMS beam size after propagation by propLen


        #***********Gaussian Beam Source
        GsnBm = srwlib.SRWLGsnBm() #Gaussian Beam structure (just parameters)
        GsnBm.x = 0 #Transverse Positions of Gaussian Beam Center at Waist [m]
        GsnBm.y = 0
        numsig = 5. #Number of sigma values to track. Total range is 2*numsig*sig_s
        ds = 2*numsig*sig_s/(nslice - 1)
        self._pulse_pos = -numsig*sig_s+slice_index*ds
        GsnBm.z = propLen + self._pulse_pos #Longitudinal Position of Waist [m]
        GsnBm.xp = 0 #Average Angles of Gaussian Beam at Waist [rad]
        GsnBm.yp = 0
        GsnBm.avgPhotEn = phE #Photon Energy [eV]
        GsnBm.pulseEn = pulseE*np.exp(-self._pulse_pos**2/(2*sig_s**2)) #Energy per Pulse [J] - to be corrected
        GsnBm.repRate = 1 #Rep. Rate [Hz] - to be corrected
        GsnBm.polar = poltype #1- linear horizontal?
        GsnBm.sigX = sigrW #Horiz. RMS size at Waist [m]
        GsnBm.sigY = GsnBm.sigX #Vert. RMS size at Waist [m]

        GsnBm.sigT = 10e-15 #Pulse duration [s] (not used?)
        GsnBm.mx = mx #Transverse Gauss-Hermite Mode Orders
        GsnBm.my = my

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

        sampFactNxNyForProp = sampFact #sampling factor for adjusting nx, ny (effective if > 0)
        arPrecPar = [sampFactNxNyForProp]

        srwlib.srwl.CalcElecFieldGaussian(_wfr, GsnBm, arPrecPar)

        ##Beamline to propagate to waist

        optDriftW=srwlib.SRWLOptD(propLen)
        propagParDrift = [0, 0, 1., 0, 0, 1.1, 1.2, 1.1, 1.2, 0, 0, 0]
        optBLW = srwlib.SRWLOptC([optDriftW],[propagParDrift])
        srwlib.srwl.PropagElecField(_wfr, optBLW)
        self.wfr = _wfr
