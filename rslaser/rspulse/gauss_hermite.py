# -*- coding: utf-8 -*-
u"""A Hermite-Gaussian laser field of order (m,n).
Copyright (c) 2021 RadiaSoft LLC. All rights reserved
"""

# SciPy imports
import numpy as np
from numpy.polynomial.hermite import hermval

from pykern.pkcollections import PKDict

# get some physical and mathematical constants ready to go
# this code snippet is adapted from rsbeams.rsphysics.rsconst.py
import math, cmath
import scipy.constants as const

TWO_PI = 2 * math.pi
RT_TWO_PI = math.sqrt(2*math.pi)
RT_2_OVER_PI = math.sqrt(2/math.pi)

C_SQ = const.c**2
C_INV  = 1./const.c
MKS_FACTOR = 1./(4.*math.pi*const.epsilon_0)
M_E_EV = const.m_e * C_SQ / (-const.e)

class GaussHermite:
    """Module defining a Hermite-Gaussian laser field of order (m,n).
    
    For now, we assume linear polarization of E along, x.
    Also, the evaluation is done for z=0 (for simplicity)
    The time variable t is ignored.
    """
    
    def __init__(self, kwargs):
        """Four input parametersrequired; waist is assumed to be round and at the origin.

        Args:
            a0:         normalized amplitude
            w0:         waist size           [m]
            lambda0:    central wavelength   [m]
            tau_fhwm:   FWHM pulse length    [s]
            
        Returns:
            instance of the class
            
        Raises:
            NA
        """
        _k=kwargs.copy()
        self.lambda0 = _k.lambda0     # central wavelength [m]
        
        # useful derived quantities
        self.k0 = 1. / self.lambda0
        self.f0 = self.k0 * const.c
        self.omega0 = TWO_PI * self.f0

        # Peak electric field [V/m]
        self.a0 = _k.a0               # amplitude [dimensionless]
        self.efield0 = self.a0 * const.m_e * self.omega0 * const.c / (const.e)

        self.set_waist_x(_k.w0)       # horizontal waist size [m]
        self.set_waist_y(_k.w0)       # vertical waist size [m]
        self.tau_fwhm = _k.tau_fwhm   # FWHM laser pulse length [s]

        # set defaults
        self.xShift = 0.     # horizontal shift of the spot center
        self.yShift = 0.     # vertical shift of the spot center
        self.zWaistX = 0.    # location (in z) of horiz. waist
        self.zWaistY = 0.    # location (in z) of  vert. waist

        self.setCoeffSingleModeX(0, 1.)       # horiz. coeff of Hermite expansion 
        self.setCoeffSingleModeY(0, 1.)       # vert.  coeff of Hermite expansion 

        self.wRotAngle = 0.

        return

    # set the horiz. waist size [m]
    def set_waist_x(self, _waistX):
        # error handling; very small waist will cause performance issues
        wFac = 20.
        minSize = wFac*self.lambda0
        if _waistX >= minSize:
            self.waistX  = _waistX
        else: 
            message = 'waistX = ' + str(_waistX) + '; must be > ' + str(minSize)
            raise Exception(message) 
            
        self.piWxFac = math.sqrt(RT_2_OVER_PI/self.waistX)
        self.zRx = 0.5*self.k0*self.waistX**2  # horiz. Rayleigh range [m]
        self.qx0 = 0.0 + self.zRx*1j
        return

    # set the vert.  waist size [m]
    def set_waist_y(self, _waistY):
        # error handling; very small waist will cause performance issues
        wFac = 20.
        minSize = wFac*self.lambda0
        if _waistY >= minSize:
            self.waistY  = _waistY
        else: 
            message = 'waistY = ' + str(_waistY) + '; must be > ' + str(minSize)
            raise Exception(message) 

        self.piWyFac = math.sqrt(RT_2_OVER_PI/self.waistY)
        self.zRy = 0.5*self.k0*self.waistY**2  #  vert. Rayleigh range [m]
        self.qy0 = 0.0 + self.zRy*1j
        return

    # set array of horizontal coefficients (complex)
    def setMCoef(self,hCoefs):
        self.mMax = hCoefs.size
        self.hCoefs = hCoefs
        return

    # set array of vertical coefficients (complex)
    def setNCoef(self,vCoefs):
        self.nMax = vCoefs.size
        self.vCoefs = vCoefs
        return

    # set horiz. mode number & coeff for single mode
    def setCoeffSingleModeX(self,mMode,mCoef):
        self.mMax = mMode + 1
        self.hCoefs = np.zeros(self.mMax) + np.zeros(self.mMax)*1j
        self.hCoefs[mMode] = mCoef
        return

    # set vertical mode num. & coeff for single mode
    def setCoeffSingleModeY(self,nMode,nCoef):
        self.nMax = nMode + 1
        self.vCoefs = np.zeros(self.nMax) + np.zeros(self.nMax)*1j
        self.vCoefs[nMode] = nCoef
        return

    # For now, we assume this is the polarization direction
    # Handling of arguments is flexible:
    #   x,y,z,t can all be scalar, to evaluate at a single point.
    #   x,y can both be arrays (same length) to evaluate on a mesh.
    #   t can be an array, to evaluate at a sequence of times.
    #   x,y,t can all be arrays, for a particle distribution with fixed z
    #   z, the longitudinal coordinate, must always be a scalar.
    def evaluate_ex(self,xArray,yArray,z,tArray):

        # get the complex-valued envelope function
        result = self.evaluate_envelope_ex(xArray,yArray,z)

        # multiply by the time-dependent term
        result *= np.exp((self.omega0 * tArray - self.k0 * z)*1j)

        # return only the real part
        return np.real(result)

    # For now, we assume this is the polarization direction
    # Handling of arguments is flexible:
    #   x,y,z can all be scalar, to evaluate at a single point.
    #   x,y can both be arrays (same length) to evaluate on a mesh.
    #   z, the longitudinal coordinate, must always be a scalar.
    def evaluate_envelope_ex(self,xArray,yArray,z):

        # assume array input; try to create temporary array
        try:
            numVals = xArray.size
            result  = np.zeros(numVals, complex)
        except AttributeError:
            # above failed, so input must be a float
            result = 0.


        # rotate and shift the x,y coordinates as necessary
        rotArg = (xArray-self.xShift)*math.cos(self.wRotAngle) \
               + (yArray-self.yShift)*math.sin(self.wRotAngle)

        # z-dependent temporary variables
        qxz   = (z-self.zWaistX) + self.qx0
        xrFac = 0.5 * self.k0 / qxz
        xzFac = math.sqrt(2)/self.waistX/math.sqrt(1+((z-self.zWaistX)/self.zRx)**2)

        # load up array of mode-dependent factors
        xCoefs = np.zeros(self.mMax, complex)
        for mMode in range(self.mMax):
            xCoefs[mMode] = self.hCoefs[mMode] * self.piWxFac * cmath.sqrt(self.qx0/qxz)     \
                          / math.sqrt(math.factorial(mMode)*(2.**mMode))                     \
                          * (self.qx0*qxz.conjugate()/self.qx0.conjugate()/qxz)**(0.5*mMode)

        # evaluate the product of Hermite series
        result = hermval(xzFac*rotArg, xCoefs)
        result *= np.exp(-(xrFac*rotArg**2)*1j)

#
# rinse and repeat:  do the same for the y-dependent Hermite series --
#

        # rotate and shift the x,y coordinates as necessary
        rotArg = (yArray-self.yShift)*math.cos(self.wRotAngle) \
               - (xArray-self.xShift)*math.sin(self.wRotAngle)

        # z-dependent temporary variables
        qyz   = (z-self.zWaistY) + self.qy0
        yrFac = 0.5 * self.k0 / qyz
        yzFac = math.sqrt(2)/self.waistY/math.sqrt(1+((z-self.zWaistY)/self.zRy)**2)

        # load up array of mode-dependent factors
        xCoefs = np.zeros(self.mMax, complex)
        for mMode in range(self.mMax):
            xCoefs[mMode] = self.hCoefs[mMode] * self.piWxFac * cmath.sqrt(self.qx0/qxz)     \
                          / math.sqrt(math.factorial(mMode)*(2.**mMode))                     \
                          * (self.qx0*qxz.conjugate()/self.qx0.conjugate()/qxz)**(0.5*mMode)

        # load up array of mode-dependent factors (y-dependence)
        yCoefs = np.zeros(self.nMax, complex)
        for nMode in range(self.nMax):
            yCoefs[nMode] = self.vCoefs[nMode] * self.piWyFac * cmath.sqrt(self.qy0/qyz)     \
                          / math.sqrt(math.factorial(nMode)*(2.**nMode))                     \
                          * (self.qy0*qyz.conjugate()/self.qy0.conjugate()/qyz)**(0.5*nMode)

        # evaluate product of Hermite series (multiplying into previous value)
        result *= hermval(yzFac*rotArg, yCoefs)
        result *= np.exp(-(yrFac*rotArg**2)*1j)

        # return the complex valued result
        return result

    def evaluate_ey(self,xArray,yArray,z,t):
        return 0. 

    def evaluate_ez(self,x,y,z,t):
        return 0.

    def evaluate_bx(self,x,y,z,t):
        return 0.

    def evaluate_by(self,x,y,z,t):
        return 0.

    def evaluate_bz(self,x,y,z,t):
        return 0.

    def tbd_or_delete():
        self._z_center = k.z_center   # longitudinal location of laser pulse center [m]
        self._z_waist = k.z_waist     # longitidunal location of nearest focus
