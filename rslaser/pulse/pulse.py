# -*- coding: utf-8 -*-
u"""Definition of a laser pulse
Copyright (c) 2021 RadiaSoft LLC. All rights reserved
"""
import math, cmath
import numpy as np
from pykern.pkdebug import pkdp, pkdlog
from pykern.pkcollections import PKDict
from numpy.polynomial.hermite import hermval
import rslaser.optics.wavefront as rswf
import rslaser.utils.constants as rsc
import rslaser.utils.unit_conversion as units
import scipy.constants as const
import srwlib
from srwlib import srwl
from rslaser.validation.validator import ValidatorBase


_LASER_PULSE_SLICE_DEFAULTS = PKDict(
        sigrW=0.000186,
        propLen=15,
        pulseE=0.001,
        poltype=1,
        sampFact=1,
        numsig=3.,
        mx=0,
        my=0
    )
_ENVELOPE_DEFAULTS = PKDict(
    phE=1.55,
    w0=.1,
    a0=.01,
    dw0x=0.0,
    dw0y=0.0,
    z_waist=0,
    dzwx=0.0,
    dzwy=0.0,
    tau_fwhm=0.1 / const.c / math.sqrt(2.),
    z_center=0,
    x_shift = 0.,
    y_shift=0.,
)
_LASER_PULSE_DEFAULTS = PKDict(
        **_ENVELOPE_DEFAULTS,
        nslice=3,
        chirp=0,
        d_to_w=0,
        slice_params=_LASER_PULSE_SLICE_DEFAULTS,
)





class InvalidLaserPulseInputError(Exception):
    pass


class LaserPulse(ValidatorBase):
    """
    The LaserPulse contains a LaserPulseEnvelope object to represent the initial envelope,
    as well as an array of LaserPulseSlice instances, which track details of the evolution in time.

    Args:
        params (PKDict):
            required fields:
                phE
                nslice
                chirp
                w0
                a0
                dw0x
                dw0y
                z_waist
                dzwx
                dzwy
                tau_fwhm
                z_center
                x_shift
                y_shift
                d_to_w
                slice_params (PKDict):
                    required fields:
                        sigrW
                        propLen
                        sig_s
                        pulseE
                        poltype
                        sampFact
                        mx
                        my

    Returns:
        instance of class
    """
    _INPUT_ERROR = InvalidLaserPulseInputError
    _DEFAULTS = _LASER_PULSE_DEFAULTS

    def __init__(self, params=None):
        params = self._get_params(params)
        self._validate_params(params)
        # instantiate the laser envelope
        e = self._get_envelope_params(params)
        self.envelope = LaserPulseEnvelope(e)
        # instantiate the array of slices
        self.slice = []
        self.nslice = params.nslice
        self.phE = params.phE
        self._lambda0 = abs(units.calculate_lambda0_from_phE(params.phE))
        self.phE -= 0.5*params.chirp           # so central slice has the central photon energy
        _de = params.chirp / self.nslice   # photon energy shift from slice to slice
        s = params.copy()
        for i in range(params.nslice):
            # add the slices; each (slowly) instantiates an SRW wavefront object
            self.slice.append(LaserPulseSlice(i, s))
            s.phE += _de
        self._sxvals = []  # horizontal slice data
        self._syvals = []  # vertical slice data

    def _get_envelope_params(self, params):
        e = PKDict()
        for k in params:
            if k in _ENVELOPE_DEFAULTS:
                e[k] = params[k]
        return e

    def compute_middle_slice_intensity(self):
        wfr = self.slice[len(self.slice) // 2].wfr
        (ar2d, sx, sy, xavg, yavg) = rswf.rmsWavefrontIntensity(wfr)
        self._sxvals.append(sx)
        self._syvals.append(sy)
        return (wfr, ar2d, sx, sy, xavg, yavg)

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


class LaserPulseSlice(ValidatorBase):
    """
    This class represents a longitudinal slice in a laser pulse.
    There will be a number of wavefronts each with different wavelengths (energy).
    The slice is composed of an SRW wavefront object, which is defined here:
    https://github.com/ochubar/SRW/blob/master/env/work/srw_python/srwlib.py#L2048

    Args:
        slice_index (int): index of slice
        params (PKDict): see slice_params field in input params to LaserPulse class __init__

    Returns:
        instance of class
    """
    _INPUT_ERROR = InvalidLaserPulseInputError
    _DEFAULTS = _LASER_PULSE_DEFAULTS

    def __init__(self, slice_index, params=None):
        #print([sigrW,propLen,pulseE,poltype])
        self._validate_type(slice_index, int, 'slice_index')
        params = self._get_params(params)
        self._validate_params(params)
        self._lambda0 = units.calculate_lambda0_from_phE(params.phE)
        self.slice_index = slice_index
        self.phE = params.phE
        constConvRad = 1.23984186e-06/(4*3.1415926536)  ##conversion from energy to 1/wavelength
        rmsAngDiv = constConvRad/(self.phE*params.slice_params.sigrW)             ##RMS angular divergence [rad]
        #  if at t=0 distance to the waist location d_to_w < d_to_w_cutoff, initialization in SRW involves/requires propagation
        #  from the distance-to-waist > d_to_w_cutoff to the actual z(t=0) for which d_to_w < d_to_w_cutoff
        d_to_w_cutoff = 0.001  # [m] - verify that this is a reasonable value
        if params.d_to_w > d_to_w_cutoff:
            params.slice_params.propLen = params.d_to_w  #  d_to_w = L_d1 +0.5*L_c in the single-pass example
        sigrL=math.sqrt(params.slice_params.sigrW**2+(params.slice_params.propLen*rmsAngDiv)**2)  ##required RMS size to produce requested RMS beam size after propagation by propLen


        #***********Gaussian Beam Source
        GsnBm = srwlib.SRWLGsnBm() #Gaussian Beam structure (just parameters)
        GsnBm.x = 0 #Transverse Positions of Gaussian Beam Center at Waist [m]
        GsnBm.y = 0
        sig_s = params.tau_fwhm * const.c / 2.355
        ds = 2*params.slice_params.numsig*sig_s/params.nslice
        self._pulse_pos = -params.slice_params.numsig*sig_s+slice_index*ds
        GsnBm.z = params.slice_params.propLen + self._pulse_pos #Longitudinal Position of Waist [m]
        GsnBm.xp = 0 #Average Angles of Gaussian Beam at Waist [rad]
        GsnBm.yp = 0
        GsnBm.avgPhotEn = self.phE #Photon Energy [eV]
        GsnBm.pulseEn = params.slice_params.pulseE*np.exp(-self._pulse_pos**2/(2*sig_s**2)) #Energy per Pulse [J] - to be corrected
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


    def _get_params(self, params):
        return self.__fixup_slice_params(super()._get_params(params))

    def __fixup_slice_params(self, params):
        self._validate_type(params.slice_params, PKDict, 'params.slice_params')
        for s in self._DEFAULTS.slice_params:
            if s not in params.slice_params:
                params.slice_params[s] = self._DEFAULTS.slice_params[s]
        return params

    def _validate_params(self, input_params):
        super()._validate_params(input_params)
        for p in input_params.slice_params:
            if p not in self._DEFAULTS.slice_params.keys():
                raise self._INPUT_ERROR(f'invalid inputs: {p} is not a parameter to {self.__class__}')


class LaserPulseEnvelope(ValidatorBase):
    """Module defining a Hermite-Gaussian laser field of order (m,n).

    For now, we assume linear polarization of E along, x.
    Also, the evaluation is done for z=0 (for simplicity)
    The time variable t is ignored.

        Args:
            params (PKDict):
                required fields:
                    phE
                    w0
                    a0
                    dw0x
                    dw0y
                    z_waist
                    dzwx
                    dzwy
                    tau_fwhm
                    z_center
                    x_shift
                    y_shift

    Returns:
        instance of class

    Note:
        waist is assumed to be round and at the origin.
    """
    _INPUT_ERROR = InvalidLaserPulseInputError
    _DEFAULTS = _ENVELOPE_DEFAULTS

    def __init__(self, params=None):
        params = self._get_params(params)
        self._validate_params(params)
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

    def set_waist_x(self, _waistX):
        """
            set the horizontal waist size [m]

            Note:
                error handling; very small waist will cause performance issues
        """
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

    def set_waist_y(self, _waistY):
        """
            set the vertical waist size [m]

            Note:
                error handling; very small waist will cause performance issues
        """
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

    def set_z_waist_x(self, _zWaistX):
        """
            set longitudinal position of the horizontal waist [m]
        """
        self.z_waist_x = _zWaistX
        self.dzwx = _waistX - self.z_waist
        return

    def set_z_waist_y(self, _zWaistY):
        """
            set longitudinal position of the vertical waist [m]
        """
        self.z_waist_y = _zWaistY
        self.dzwy = _waistY - self.z_waist
        return

    def setMCoef(self,hCoefs):
        """
            set array of horizontal coefficients (complex)
        """
        self.mMax = hCoefs.size
        self.hCoefs = hCoefs
        return

    def setNCoef(self,vCoefs):
        """
            set array of vertical coefficients (complex)
        """
        self.nMax = vCoefs.size
        self.vCoefs = vCoefs
        return

    def setCoeffSingleModeX(self,mMode,mCoef):
        """
            set horiz. mode number & coeff for single mode
        """
        self.mMax = mMode + 1
        self.hCoefs = np.zeros(self.mMax) + np.zeros(self.mMax)*1j
        self.hCoefs[mMode] = mCoef
        return

    def setCoeffSingleModeY(self,nMode,nCoef):
        """
            set vertical mode num. & coeff for single mode
        """
        self.nMax = nMode + 1
        self.vCoefs = np.zeros(self.nMax) + np.zeros(self.nMax)*1j
        self.vCoefs[nMode] = nCoef
        return

    def evaluate_ex(self,xArray,yArray,_z,tArray):
        """
            For now, we assume this is the polarization direction

            Args:
                x,y,z,t can all be scalar, to evaluate at a single point.
                x,y can both be arrays (same length) to evaluate on a mesh.
                t can be an array, to evaluate at a sequence of times.
                x,y,t can all be arrays, for a particle distribution with fixed z
                _z, the longitudinal coordinate, must always be a scalar.
        """

        # account for location of pulse center
        z_local = _z - self.z_center

        # get the complex-valued envelope function
        result = self.evaluate_envelope_ex(xArray, yArray, _z)

        # multiply by the time-dependent term
        return result * np.exp((self.omega0*tArray - self.k0*z_local)*1j)


    def evaluate_envelope_ex(self,xArray,yArray,_z):
        """
            For now, we assume this is the polarization direction

            Args:
                x,y,z can all be scalar, to evaluate at a single point.
                x,y can both be arrays (same length) to evaluate on a mesh.
                _z, the longitudinal coordinate, must always be a scalar.

            Note:
                We ignore x/y differences in the waist size and location
                Also, we ignore the higher-order Hermite modes here.
        """
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

    def evaluate_er(self,rArray,_z,tArray):
        """
            Evaluate the radial electric field of a circularly polarized laser pulse in r-z geometry.

            Args:
               rArray can be a scalar, to evaluate at a single point.
               rArray can be a Numpy array to evaluate along a line.
               _z, the longitudinal coordinate, must always be a scalar.
               tArray can be a scalar (works) or a Numpy array (not tested)
        """

        # account for location of pulse center
        z_local = _z - self.z_center

        # get the complex-valued envelope function
        result = self.evaluate_envelope_er(rArray, _z)

        # multiply by the time-dependent term
        return result * np.exp((self.omega0*tArray - self.k0*z_local)*1j)

    def evaluate_envelope_er(self,rArray,_z):
        """
            Calculate the laser pulse envelope in radial r-z coordinates

            Args:
                rArray can be a scalar, to evaluate at a single point.
                rArray can be a Numpy array to evaluate along a line
                _z, the longitudinal coordinate, must always be a scalar.

            Note:
                We ignore x/y differences in the waist size and location
                Also, we ignore the higher-order Hermite modes here.
        """
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

    def eval_gh_ex(self,xArray,yArray,z):
        """
        Note: This is old, untested code with no concept of a longitudinal envelope
        The complicated logic requires further testing.

        Args:
            x,y,z can all be scalar, to evaluate at a single point.
            x,y can both be arrays (same length) to evaluate on a mesh.
            z, the longitudinal coordinate, must always be a scalar.
        """

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
