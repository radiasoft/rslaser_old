# -*- coding: utf-8 -*-
u"""A Hermite-Gaussian laser field of order (m,n).
Copyright (c) 2021 RadiaSoft LLC. All rights reserved
"""
import math, cmath
import numpy as np
from numpy.polynomial.hermite import hermval
from pykern.pkcollections import PKDict
from pykern.pkdebug import pkdc, pkdexc, pkdp
import rslaser.utils.constants as rsc
import rslaser.utils.unit_conversion as units
import scipy.constants as const

class GaussHermite:
    """Module defining a Hermite-Gaussian laser field of order (m,n).

    For now, we assume linear polarization of E along, x.
    Also, the evaluation is done for z=0 (for simplicity)
    The time variable t is ignored.

    Four input parametersrequired; waist is assumed to be round and at the origin.

        Args:
            kwargs:  this is a PKDict object

        Returns:
            instance of the class

        Raises:
            NA
    """

    def __init__(self, input_params):
        params = input_params.copy()
        self.lambda0 = abs(units.calculate_lambda0_from_phE(params.phE))            # central wavelength [m]
        # useful derived quantities
        self.k0 = rsc.TWO_PI / self.lambda0       # central wavenumber [radians/m]
        self.f0 = const.c / self.lambda0          # central frequency  [Hz]
        self.omega0 = rsc.TWO_PI * self.f0        # central angular frequency [radians/s]

        # Peak electric field [V/m]
        self.a0 = abs(params.a0)                      # amplitude [dimensionless]
                                                  # peak electric field [V/m]
        self.efield0 = self.a0 * const.m_e * self.omega0 * const.c / (const.e)

        # waist sizes and locations
        self.w0 = abs(params.w0)                      # the waist size of the pulse
        self.set_waist_x(params.w0 + params.dw0x)         # horizontal waist size [m]
        self.set_waist_y(params.w0 + params.dw0y)         # vertical waist size [m]

        self.z_waist = params.z_waist                 # the longitudinal location of the waist
        self.z_waist_x = params.dzwx                  # location (in z) of horizontal waist
        self.z_waist_y = params.dzwy                  # location (in z) of vertical waist

        # Rayleigh range
        self.zR = 0.5*self.k0*(self.w0)**2        # Rayleigh range, ignoring horizontal/vertical differences
#        print('\n ****** \n zR = ', self.zR)

        # pulse length
        self.tau_fwhm = params.tau_fwhm               # FWHM laser pulse length [s]
        self.L_fwhm = self.tau_fwhm * const.c     # FWHM laser pulse length [m]

        # longitudinal location of pulse center
        self.z_center = params.z_center

        # bulk transverse offsets of the laser pulse
        self.x_shift = params.x_shift                 # horizontal shift of the spot center
        self.y_shift = params.y_shift                 # vertical shift of the spot center

        # for now, we set the higher Hermite modes to zero
        self.setCoeffSingleModeX(0, 1.)           # horizontal coefficients of Hermite expansion
        self.setCoeffSingleModeY(0, 1.)           # vertical coefficients of Hermite expansion

        # for now, we set the rotation angle to zero
        self.wRotAngle = 0.

        return

    # set the horizontal waist size [m]
    def set_waist_x(self, _waistX):
        # error handling; very small waist will cause performance issues
        wFac = 4.0
        minSize = wFac*self.lambda0
        if _waistX >= minSize:
            self.waist_x  = _waistX
        else:
            message = 'waistX = ' + str(_waistX) + '; must be >= ' + str(minSize)
            raise Exception(message)

        # error handling; require that deviations from w0 are small
        self.dw0x = _waistX - self.w0
        if abs(self.dw0x) > 0.1 * self.w0:
            message = 'dw0x/w0 = ' + str(self.dw0x) + '; must be < 0.1 '
            raise Exception(message)

        self.piWxFac = math.sqrt(rsc.RT_2_OVER_PI/self.waist_x)
        self.zRx = 0.5*self.k0*self.waist_x**2     # horizintal Rayleigh range [m]
        self.qx0 = 0.0 + self.zRx*1j
        return

    # set the vertical waist size [m]
    def set_waist_y(self, _waistY):
        # error handling; very small waist will cause performance issues
        wFac = 4.0
        minSize = wFac*self.lambda0
        if _waistY >= minSize:
            self.waist_y  = _waistY
        else:
            message = 'waistY = ' + str(_waistY) + '; must be >= ' + str(minSize)
            raise Exception(message)

        # error handling; require that deviations from w0 are small
        self.dw0y = _waistY - self.w0
        if abs(self.dw0y) > 0.1 * self.w0:
            message = 'dw0y/w0 = ' + str(self.dw0y) + '; must be < 0.1 '
            raise Exception(message)

        self.piWyFac = math.sqrt(rsc.RT_2_OVER_PI/self.waist_y)
        self.zRy = 0.5*self.k0*self.waist_y**2      #  vertical Rayleigh range [m]
        self.qy0 = 0.0 + self.zRy*1j

        return

    # set longitudinal position of the horizontal waist [m]
    def set_z_waist_x(self, _zWaistX):
        self.z_waist_x = _zWaistX
        self.dzwx = _waistX - self.z_waist
        return

    # set longitudinal position of the vertical waist [m]
    def set_z_waist_y(self, _zWaistY):
        self.z_waist_y = _zWaistY
        self.dzwy = _waistY - self.z_waist
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
    #   _z, the longitudinal coordinate, must always be a scalar.
    def evaluate_ex(self,xArray,yArray,_z,tArray):

        # account for location of pulse center
        z_local = _z - self.z_center

        # get the complex-valued envelope function
        result = self.evaluate_envelope_ex(xArray, yArray, _z)

        # multiply by the time-dependent term
        return result * np.exp((self.omega0*tArray - self.k0*z_local)*1j)

    # For now, we assume this is the polarization direction
    # Handling of arguments is flexible:
    #   x,y,z can all be scalar, to evaluate at a single point.
    #   x,y can both be arrays (same length) to evaluate on a mesh.
    #   _z, the longitudinal coordinate, must always be a scalar.
    #
    # We ignore x/y differences in the waist size and location
    # Also, we ignore the higher-order Hermite modes here.

    def evaluate_envelope_ex(self,xArray,yArray,_z):
        # account for location of pulse center
        z_local = _z - self.z_center
        # account for the waist location
        _z -= self.z_waist

        # determine whether xArray is really a Numpy array
        try:
            num_vals_x = xArray.size
            x_is_array = True
        except AttributeError:
            # above failed, so input must be a float
            x_is_array = False

        # determine whether yArray is really a Numpy array
        try:
            num_vals_y = yArray.size
            y_is_array = True
        except AttributeError:
            # above failed, so input must be a float
            y_is_array = False

        if (x_is_array and y_is_array):
            rSq = np.zeros(num_vals_x, complex)
            exp_1 = np.zeros(num_vals_x, complex)
            exp_2 = np.zeros(num_vals_x, complex)
            arg_2 = np.zeros(num_vals_x, complex)

        # radius at which the field amplitudes fall to exp(-1) of their axial values
        #     i.e., where the intensity values fall to exp(-2)
        wZ = self.w0 * math.sqrt(1+(_z/self.zR)**2)
#        pkdc('w(z)/w0 = ' + str(wZ/self.w0))

        # the radius squared
        rSq = np.power(xArray,2) + np.power(yArray,2)

        # the radius of curvature of wavefronts at location z
        invR = _z / (_z**2 + self.zR**2)

        # first exponential
        exp_1 = np.exp(-rSq / wZ**2)

        # Gouy phase at position z
        psi_z = np.arctan(_z / self.zR)

        # 2nd exponential
        arg_2 = 0.5*self.k0*invR*rSq
        exp_2 = np.exp(-1j*(arg_2 - psi_z))

#        pkdc(' k0 = ' + str(self.k0))
#        pkdc(' invR = ' + str(invR))
#        pkdc(' rSq = ' + str(rSq))
#        pkdc(' arg_2 = ' + str(arg_2))
#        pkdc(' psi_z = ' + str(psi_z))
#        pkdc(' Re[exp_2] = ' + str(np.real(exp_2)))

        # return the complex valued result
        # here, we apply a longitudinal Gaussian profile
        return (self.w0 / wZ) * exp_1 * np.exp(-(z_local/self.L_fwhm)**2) * self.efield0 * exp_2


    # Evaluate the radial electric field of a circularly polarized laser pulse in r-z geometry.
    # Handling of arguments is flexible:
    #   rArray can be a scalar, to evaluate at a single point.
    #   rArray can be a Numpy array to evaluate along a line.
    #   _z, the longitudinal coordinate, must always be a scalar.
    #   tArray can be a scalar (works) or a Numpy array (not tested)
    def evaluate_er(self,rArray,_z,tArray):

        # account for location of pulse center
        z_local = _z - self.z_center

        # get the complex-valued envelope function
        result = self.evaluate_envelope_er(rArray, _z)

        # multiply by the time-dependent term
        return result * np.exp((self.omega0*tArray - self.k0*z_local)*1j)


    # Calculate the laser pulse envelope in radial r-z coordinates
    # Handling of arguments is flexible:
    #   rArray can be a scalar, to evaluate at a single point.
    #   rArray can be a Numpy array to evaluate along a line
    #   _z, the longitudinal coordinate, must always be a scalar.
    #
    # We ignore x/y differences in the waist size and location
    # Also, we ignore the higher-order Hermite modes here.

    def evaluate_envelope_er(self,rArray,_z):
        # account for the waist location
        z_local = _z - self.z_center
        _z -= self.z_waist
#        pkdc('z_local, _z = ' + str(z_local) + ', ' + str(_z))

        # determine whether xArray is really a Numpy array
        try:
            num_vals_r = rArray.size
            r_is_array = True
            rSq = np.zeros(num_vals_r, complex)
            exp_1 = np.zeros(num_vals_r, complex)
            exp_2 = np.zeros(num_vals_r, complex)
            arg_2 = np.zeros(num_vals_r, complex)
        except AttributeError:
            # above failed, so input must be a float
            r_is_array = False

        # radius at which the field amplitudes fall to exp(-1) of their axial values
        #     i.e., where the intensity values fall to exp(-2)
        wZ = self.w0 * math.sqrt(1+(_z/self.zR)**2)
#        pkdc('w(z)/w0 = ' + str(wZ/self.w0))

        # the radius squared
        rSq = np.power(rArray,2)

        # the radius of curvature of wavefronts at location z
        invR = _z / (_z**2 + self.zR**2)

        # first exponential
        exp_1 = np.exp(-rSq / wZ**2)

        # Gouy phase at position z
        psi_z = np.arctan(_z / self.zR)

        # 2nd exponential
        arg_2 = 0.5*self.k0*invR*rSq - psi_z
        exp_2 = np.exp(-1j*arg_2)

#        pkdc(' k0 = ' + str(self.k0))
#        pkdc(' invR = ' + str(invR))
#        pkdc(' rSq = ' + str(rSq))
#        pkdc(' arg_2 min/max = ' + str(np.min(arg_2)) + ', ' + str(np.max(arg_2)))
#        pkdc(' psi_z = ' + str(psi_z))
#        pkdc(' Re[exp_2] = ' + str(np.real(exp_2)))

        # return the complex valued result
        # here, we apply a longitudinal Gaussian profile
        return (self.w0 / wZ) * exp_1 * np.exp(-(z_local/self.L_fwhm)**2) * self.efield0 * exp_2

    # This is old, untested code with no concept of a longitudinal envelope
    #   x,y,z can all be scalar, to evaluate at a single point.
    #   x,y can both be arrays (same length) to evaluate on a mesh.
    #   z, the longitudinal coordinate, must always be a scalar.
    #
    # The complicated logic below requires further testing.

    def eval_gh_ex(self,xArray,yArray,z):

        # assume array input; try to create temporary array
        try:
            numVals = xArray.size
            result  = np.zeros(numVals, complex)
        except AttributeError:
            # above failed, so input must be a float
            result = 0.


        # rotate and shift the x,y coordinates as necessary
        rotArg = (xArray-self.x_shift)*math.cos(self.wRotAngle) \
               + (yArray-self.y_shift)*math.sin(self.wRotAngle)

        # z-dependent temporary variables
        qxz   = (z-self.z_waist_x) + self.qx0
        xrFac = 0.5 * self.k0 / qxz
        xzFac = math.sqrt(2)/self.waist_x/math.sqrt(1+((z-self.z_waist_x)/self.zRx)**2)

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
        rotArg = (yArray-self.y_shift)*math.cos(self.wRotAngle) \
               - (xArray-self.x_shift)*math.sin(self.wRotAngle)

        # z-dependent temporary variables
        qyz   = (z-self.z_waist_y) + self.qy0
        yrFac = 0.5 * self.k0 / qyz
        yzFac = math.sqrt(2)/self.waist_y/math.sqrt(1+((z-self.z_waist_y)/self.zRy)**2)

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
