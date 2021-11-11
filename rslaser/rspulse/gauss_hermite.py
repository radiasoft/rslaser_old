"""
Module defining a Hermite-Gaussian laser field of order (m,n).
Copyright (c) 2021 RadiaSoft LLC. All rights reserved
"""

# SciPy imports
import numpy as np
from numpy.polynomial.hermite import hermval

from pykern.pkcollections import PKDict

# get some physical and mathematical constants ready to go
# this code snippet is adapted from rsbeams.rsphysics.rsconst.py
import math
import scipy.constants as const

TWO_PI = 2 * math.pi
RT_TWO_PI = math.sqrt(2*math.pi)
RT_2_OVER_PI = math.sqrt(2/math.pi)

c_SQ = const.c**2
c_INV  = 1./const.c
MKS_factor = 1./(4.*math.pi*const.epsilon_0)
m_e_EV = const.m_e * c_SQ / (-const.e)

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
        k=kwargs.copy()
        self._a0 = k.a0               # amplitude [dimensionless]
        self._w0x = k.w0              # horizontal waist size [m]
        self._w0y = k.w0              # vertical waist size [m]
        self._lambda0 = k.lambda0     # central wavelength [m]
        self._tau_fwhm = k.tau_fwhm   # FWHM laser pulse length [s]

        # set defaults
        self.setXShift(0.)   # horiz. shift of the spot center
        self.setYShift(0.)   #  vert. shift of the spot center
        self.setZWaistX(0.)  # location (in z) of horiz. waist
        self.setZWaistY(0.)  # location (in z) of  vert. waist

        self.setCoeffSingleModeX(0, 1.)       # horiz. coeff of Hermite expansion 
        self.setCoeffSingleModeY(0, 1.)       # vert.  coeff of Hermite expansion 

        return

    # set the horiz. waist size [m]
    def setWaistX(self, kwargs):
        # error handling; very small waist will cause performance issues
        wFac = 20.
        minSize = wFac*self.lambda0
        if waistX >= minSize:
            self.waistX  = waistX
        else: 
            message = 'waistX = ' + str(waistX) + '; must be > ' + str(minSize)
            raise Exception(message) 
            
        self.piWxFac = math.sqrt(self.rt2opi/waistX)
        self.zRx = 0.5*self.k0*waistX**2  # horiz. Rayleigh range [m]
        self.qx0 = 0.0 + self.zRx*1j
        return

    # set the vert.  waist size [m]
    def setWaistY(self,waistY):
        # error handling; very small waist will cause performance issues
        wFac = 20.
        minSize = wFac*self.lambda0
        if waistY >= minSize:
            self.waistY  = waistY
        else: 
            message = 'waistY = ' + str(waistY) + '; must be > ' + str(minSize)
            raise Exception(message) 

        self.piWyFac = math.sqrt(self.rt2opi/waistY)
        self.zRy = 0.5*self.k0*waistY**2  #  vert. Rayleigh range [m]
        self.qy0 = 0.0 + self.zRy*1j
        return

    # spot rotation angle [Rad]
    def setWRotAngle(self,wRotAngle):
        self.wRotAngle = wRotAngle
        return

    # horiz. shift of the spot center
    def setXShift(self,xShift):
        self.xShift = xShift
        return

    # vert. shift of the spot center
    def setYShift(self,yShift):
        self.yShift = yShift
        return

    # set location (in z) of horiz. waist
    def setZWaistX(self,zWaistX):
        self.zWaistX = zWaistX
        return

    # set location (in z) of  vert. waist
    def setZWaistY(self,zWaistY):
        self.zWaistY = zWaistY
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
    def evaluateEx(self,xArray,yArray,z,tArray):

        # get the complex-valued envelope function
        result = self.evalEnvelopeEx(xArray,yArray,z)

        # multiply by the time-dependent term
        result *= np.exp((self.omega * tArray - self.k0 * z)*1j)

        # return only the real part
        return np.real(result)

    # For now, we assume this is the polarization direction
    # Handling of arguments is flexible:
    #   x,y,z can all be scalar, to evaluate at a single point.
    #   x,y can both be arrays (same length) to evaluate on a mesh.
    #   z, the longitudinal coordinate, must always be a scalar.
    def evalEnvelopeEx(self,xArray,yArray,z):

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

    def evaluateEy(self,xArray,yArray,z,t):
        return 0. 

    def evaluateEz(self,x,y,z,t):
        return 0.

    def evaluateBx(self,x,y,z,t):
        return 0.

    def evaluateBy(self,x,y,z,t):
        return 0.

    def evaluateBz(self,x,y,z,t):
        return 0.

    def tbd():
        self._z_center = k.z_center   # longitudinal location of laser pulse center [m]
        self._z_waist = k.z_waist     # longitidunal location of nearest focus

        # useful derived quantities
        self._k_0 = 1. / self._lambda_0
        self._f_0 = self._k_0 * const.c
        self._omega_0 = TWO_PI * self._f_0

        # Peak electric field [V/m]
        self._efield_0 = self._a_0 * const.m_e * self._omega_0 * const.c / (const.e)

        # FWHM pulse length
        self._tau_fwhm = k.tau_fwhm   # FWHM laser pulse length [s]
        self._L_fwhm = self._tau_fwhm * const.c

        # set input data (cannot change lambda0, k0, omega)
        self.lambda0 = lambda0         # central wavelength [m]
        self.k0 = 2.*math.pi/lambda0   # central wavenumber [Rad/m]
        self.omega = self.k0 * self.c  # central wavenumber [Rad/m]
        self.setWaistX(waistX)         # horiz. waist size [m]
        self.setWaistY(waistY)         #  vert. waist size [m]
        self.setWRotAngle(wRotAngle)   # spot rotation angle [Rad]
