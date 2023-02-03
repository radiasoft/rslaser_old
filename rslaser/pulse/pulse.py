# -*- coding: utf-8 -*-
"""Definition of a laser pulse
Copyright (c) 2021 RadiaSoft LLC. All rights reserved
"""
import array
import math
import cmath
import numpy as np
from pykern.pkdebug import pkdp, pkdlog
from pykern.pkcollections import PKDict
from numpy.polynomial.hermite import hermval
import rslaser.optics.wavefront as rswf
import rsmath.const as rsc
import rslaser.utils.unit_conversion as units
import rslaser.utils.srwl_uti_data as srwutil
import scipy.constants as const
import srwlib
from srwlib import srwl
from rslaser.utils.validator import ValidatorBase


_LASER_PULSE_DEFAULTS = PKDict(
        nslice = 3,
        chirp = 0,
        photon_e_ev = 1.5, #1e3,
        num_sig_long=3.,
        dist_waist = 0,
        tau_fwhm= 0.1 / const.c / math.sqrt(2.),
        pulseE = 0.001,
        sigx_waist = 1.0e-3,
        sigy_waist = 1.0e-3,
        num_sig_trans = 6,
        nx_slice = 500,
        ny_slice = 500,
        poltype = 1,
        mx = 0,
        my = 0,
)
_ENVELOPE_DEFAULTS = PKDict(
    w0=.1,
    a0=.01,
    dw0x=0.0,
    dw0y=0.0,
    dzwx=0.0,
    dzwy=0.0,
    z_center=0,
    x_shift = 0.,
    y_shift=0.,
    photon_e_ev = 1.5, #1e3,
    tau_fwhm= 0.1 / const.c / math.sqrt(2.),
)


class InvalidLaserPulseInputError(Exception):
    pass


class LaserPulse(ValidatorBase):
    """
    The LaserPulse contains an array of LaserPulseSlice instances, which track
    details of the evolution in time.
    Args:
        params (PKDict):
                photon_e_ev (float): Photon energy [eV]
                nslice (int): number of slices
                chirp (float): energy variation from first to last slice in laser pulse [eV]
                dist_waist (float): distance from waist at which initial wavefront is calculated [m]
                w0 (float): beamsize of laser pulse at waist [m]
                a0 (float): laser amplitude, a=0.85e-9 lambda[micron] sqrt(I[W/cm^2])
                dw0x (float): horizontal variation in waist size [m]
                dw0y (float): vertical variation in waist size [m]
                dzwx (float): location (in z) of horizontal waist [m]
                dzwy (float): location (in z) of vertical waist [m]
                tau_fwhm (float): FWHM laser pulse length [s]
                pulseE (float): total laser pulse energy [J]
                z_center (float): # longitudinal location of pulse center [m]
                x_shift (float): horizontal shift of the spot center [m]
                y_shift (float): vertical shift of the spot center [m]
                sigx_waist (float): horizontal RMS waist size [m]
                sigy_waist (float): vertical RMS waist size [m]
                nx_slice (int): no. of horizontal mesh points in slice
                ny_slice (int): no. of vertical mesh points in slice
                num_sig_trans (int): no. of sigmas for transverse Gsn range
                pulseE (float): maximum pulse energy for SRW Gaussian wavefronts [J]
                poltype (int): polarization 1- lin. hor., 2- lin. vert., 3- lin. 45 deg., 4- lin.135 deg., 5- circ. right, 6- circ. left
                mx (int): transverse Gauss-Hermite mode order in horizontal direction
                my (int): transverse Gauss-Hermite mode order in vertical direction

    Returns:
        instance of class with attributes:
            slice: list of LaserPulseSlices each with an SRW wavefront object
            nslice: number of slices
            photon_e_ev: Photon energy [eV]
            sig_s: RMS bunch length [m]
            _lambda0: central wavelength [m]
            _sxvals: RMS horizontal beam size of each slice [m]
            _syvals: RMS vertical beam size of each slice [m]
    """
    _INPUT_ERROR = InvalidLaserPulseInputError
    _DEFAULTS = _LASER_PULSE_DEFAULTS

    def __init__(self, params=None, files=None):
        params = self._get_params(params)

        self._validate_params(params, files)
        # instantiate the array of slices
        self.slice = []
        self.files = files
        self.sigx_waist = params.sigx_waist
        self.sigy_waist = params.sigy_waist
        self.num_sig_trans = params.num_sig_trans
        self.nslice = params.nslice
        self.photon_e_ev = params.photon_e_ev
        self.sig_s = params.tau_fwhm * const.c / 2.355
        self.num_sig_long = params.num_sig_long
        self._lambda0 = abs(units.calculate_lambda0_from_phE(params.photon_e_ev *const.e)) # Function requires energy in J
        self.pulseE = params.pulseE
        # self.photon_e -= 0.5*params.chirp           # so central slice has the central photon energy
        # _de = params.chirp / self.nslice   # photon energy shift from slice to slice
        for i in range(params.nslice):
            # add the slices; each (slowly) instantiates an SRW wavefront object
            self.slice.append(LaserPulseSlice(i, params.copy(), files=self.files))
        self._sxvals = []  # horizontal slice data
        self._syvals = []  # vertical slice data

    def _validate_params(self, input_params, files):
        if files and input_params.nslice > 1:
            raise self._INPUT_ERROR("cannot use file inputs with more than one slice")
        super()._validate_params(input_params)

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
        return [s.photon_e_ev for s in self.slice]

    def slice_wfr(self,slice_index):
        return self.slice[slice_index].wfr


class LaserPulseSlice(ValidatorBase):
    """
    This class represents a longitudinal slice in a laser pulse.
    There will be a number of wavefronts each with different wavelengths (energy).
    The slice is composed of an SRW wavefront object, which is defined here:
    https://github.com/ochubar/SRW/blob/master/env/work/srw_python/srwlib.py#L2048

    Args:
        slice_index (int): index of slice
        params (PKDict): accepts input params from LaserPulse class __init__

    Returns:
        instance of class
    """

    _INPUT_ERROR = InvalidLaserPulseInputError
    _DEFAULTS = _LASER_PULSE_DEFAULTS

    def __init__(self, slice_index, params=None, files=None):
        self._validate_type(slice_index, int, 'slice_index')
        params = self._get_params(params)
        self._validate_params(params)
        self._lambda0 = units.calculate_lambda0_from_phE(params.photon_e_ev *const.e) # Function requires energy in J
        self.slice_index = slice_index
        self.sigx_waist = params.sigx_waist
        self.sigy_waist = params.sigy_waist
        self.num_sig_trans = params.num_sig_trans
        # self.z_waist = params.z_waist
        self.nslice = params.nslice
        self.nx_slice = params.nx_slice
        self.ny_slice = params.ny_slice
        self.dist_waist = params.dist_waist

        #  (Note KW: called this pulseE_slice because right now LPS is also passed pulseE for the whole pulse)
        self.pulseE_slice = (params.pulseE / self.nslice) # currently assumes consistent length and energy across all slices

        # compute slice photon energy from central energy, chirp, and slice index
        self.photon_e_ev = params.photon_e_ev          # check that this is being properly incremented in the correct place (see LaserPulse class)
        _de = params.chirp / self.nslice   # photon energy shift from slice to slice
        self.photon_e_ev -= 0.5*params.chirp + (self.nslice * _de)   # so central slice has the central photon energy
        
        self.sig_s = params.tau_fwhm * const.c / 2.355
        self.num_sig_long = params.num_sig_long
        constConvRad = 1.23984186e-06/(4*3.1415926536)  ##conversion from energy to 1/wavelength
        rmsAngDiv_x = constConvRad/(self.photon_e_ev * self.sigx_waist)             ##RMS angular divergence [rad]
        rmsAngDiv_y = constConvRad/(self.photon_e_ev * self.sigy_waist)

        sigrL_x = math.sqrt(self.sigx_waist**2 + (self.dist_waist * rmsAngDiv_x)**2)
        sigrL_y = math.sqrt(self.sigy_waist**2 + (self.dist_waist * rmsAngDiv_y)**2)

        # *************begin function below**********

        # sig_s = params.tau_fwhm * const.c / 2.355
        self.ds = 2 * params.num_sig_long * self.sig_s / params.nslice    # longitudinal spacing between slices
        # self._pulse_pos = self.dist_waist - params.num_sig_long * self.sig_s + (slice_index + 0.5) * self.ds
        self._pulse_pos = -params.num_sig_long * self.sig_s + (slice_index + 0.5) * self.ds
        self._wavefront(params, files)

        # Calculate the initial number of photons in 2d grid of each slice from pulseE_slice
        self.n_photons_2d = self.calc_init_n_photons() # 2d array


    def _wavefront(self, params, files):
        if files:
            with open(files.meta) as fh:
                for line in fh:
                    if line.startswith("pixel_size_h_microns"):
                        pixel_size_h = float(line.split(":")[-1].split(",")[0])  # microns
                    if line.startswith("pixel_size_v_microns"):
                        pixel_size_v = float(line.split(":")[-1].split(",")[0])  # microns

            # central wavelength of the laser pulse
            lambda0_micron = 0.8
            ccd_data = np.genfromtxt(files.ccd, skip_header=1)
            ccd_data = _reshape_data(ccd_data)

            # specify the mesh size
            nx = ccd_data.shape[1]
            ny = ccd_data.shape[0]

            # create the x,y arrays with physical units based on the diagnostic pixel dimensions
            x_max = 0.002    # [m]
            x_min = -x_max
            y_max = x_max
            y_min = -y_max

            # parse the measured phases of the wavefront
            wfs_data = np.genfromtxt(files.wfs, skip_header=1, skip_footer=0)

            # clean up any NaN's
            indices = np.isnan(wfs_data)
            wfs_data = _reshape_data(_array_cleaner(wfs_data, indices))

            # convert from microns to radians
            rad_per_micron = math.pi / lambda0_micron
            wfs_data *= rad_per_micron
            assert np.shape(wfs_data) == np.shape(ccd_data), 'ERROR -- WFS and CCD data have diferent shapes!!'

            # Calulate the real and imaginary parts of the Ex,Ey electric field components
            e_norm = np.sqrt(ccd_data)
            ex_real = np.multiply(e_norm, np.cos(wfs_data)).flatten(order='C')
            ex_imag = np.multiply(e_norm, np.sin(wfs_data)).flatten(order='C')

            ex_numpy = np.zeros(2*len(ex_real))
            for i in range(len(ex_real)):
                ex_numpy[2*i] = ex_real[i]
                ex_numpy[2*i+1] = ex_imag[i]
            ex = array.array('f', ex_numpy.tolist())
            ey = array.array('f', len(ex)*[0.])
            self.wfr = srwlib.SRWLWfr(_arEx=ex, _arEy=ey, _typeE='f',
                    _eStart=1.55, _eFin=1.55, _ne=1,
                    _xStart=x_min, _xFin=x_max, _nx=nx,
                    _yStart=y_min, _yFin=y_max, _ny=ny,
                    _zStart=0., _partBeam=None)
            return
 
        # Adjust for the length of the pulse + a constant factor to make pulseE = sum(energy_2d)
        constant_factor = 2.94e-2

        # Since pulseE = fwhm_tau * spot_size * intensity, new_pulseE = old_pulseE / fwhm_tau
        length_factor = constant_factor / self.ds

        # calculate field energy in this slice
        sliceEnInt = length_factor * self.pulseE_slice*np.exp(-self._pulse_pos**2/(2*self.sig_s**2))
        
        self.wfr = srwutil.createGsnSrcSRW(self.sigx_waist, self.sigy_waist, self.num_sig_trans, self._pulse_pos, sliceEnInt, params.poltype, \
                                           self.nx_slice, self.ny_slice, self.photon_e_ev, params.mx, params.my)

    def calc_init_n_photons(self):

        intensity = srwlib.array('f', [0]*self.wfr.mesh.nx*self.wfr.mesh.ny) # "flat" array to take 2D intensity data
        srwl.CalcIntFromElecField(intensity, self.wfr, 0, 0, 3, self.wfr.mesh.eStart, 0, 0) #extracts intensity

        # Reshaping intensity data from flat to 2D array
        intens_2d = np.array(intensity).reshape((self.wfr.mesh.nx, self.wfr.mesh.ny), order='C').astype(np.float64)

        efield_abs_sqrd_2d = np.sqrt(const.mu_0 / const.epsilon_0) * 2.0 * intens_2d # [V^2/m^2]

        dx = (self.wfr.mesh.xFin - self.wfr.mesh.xStart)/self.wfr.mesh.nx
        dy = (self.wfr.mesh.yFin - self.wfr.mesh.yStart)/self.wfr.mesh.ny

        # Field energy per grid cell is the area of that cell times the energy density
        cell_volume = self.ds * dx * dy
        energy_2d = cell_volume * (const.epsilon_0 / 2.0) * efield_abs_sqrd_2d

        # Get slice value of photon_e (will be in eV)
        photon_e = self.photon_e_ev * const.e

        # Number of photons in each grid cell can be found by dividing the
        # total energy of the laser in that grid cell by the energy of a photon
        return energy_2d / photon_e


class LaserPulseEnvelope(ValidatorBase):
    """Module defining a Hermite-Gaussian laser field of order (m,n).

    For now, we assume linear polarization of E along, x.
    Also, the evaluation is done for z=0 (for simplicity)
    The time variable t is ignored.

        Args:
            params (PKDict):
                required fields:
                    photon_e_ev (float): Photon energy [eV]
                    w0 (float): beamsize of laser pulse at waist [m]
                    a0 (float): laser amplitude, a=0.85e-9 lambda[micron] sqrt(I[W/cm^2])
                    dw0x (float): horizontal variation in waist size [m]
                    dw0y (float): vertical variation in waist size [m]
                    z_waist (float): longitudinal location of the waist [m]
                    dzwx (float): location (in z) of horizontal waist [m]
                    dzwy (float): location (in z) of vertical waist [m]
                    tau_fwhm (float): FWHM laser pulse length [s]
                    z_center (float): # longitudinal location of pulse center [m]
                    x_shift (float): horizontal shift of the spot center [m]
                    y_shift (float): vertical shift of the spot center [m]

    Returns:
        instance of class

    """
    _INPUT_ERROR = InvalidLaserPulseInputError
    _DEFAULTS = _ENVELOPE_DEFAULTS

    def __init__(self, params=None):
        params = self._get_params(params)
        self._validate_params(params)
        self.lambda0 = abs(units.calculate_lambda0_from_phE(params.photon_e_ev *const.e)) # Function requires energy in J # central wavelength [m]
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

        # self.z_waist = params.z_waist                 # the longitudinal location of the waist
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
        # _z -= self.z_waist

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


def _nan_helper(_arr):
    """
    Clean unwanted NaNs from a numpy array, replacing them via interpolation.

    Args:
        _arr, numpy array with NaNs

    Returns:
        nans, logical indices of NaNs
        index, a function with signature indices = index(logical_indices)
               to convert logical indices of NaNs to 'equivalent' indices

    Example:
        >>> nans, x = nan_helper(my_array)
        >>> my_array[nans] = np.interp(x(nans), x(~nans), my_array[~nans])
    """
    return np.isnan(_arr), lambda z: z.nonzero()[0]

def _array_cleaner(_arr, _ind):
    """
    Clean unwanted values from a numpy array, replacing them via interpolation.

    Args:
        _arr, numpy array with bad values
        _ind, precalculated indices of these bad values

    Returns:
        _arr, cleaned version of the input array

    Example:
        >>> indices = np.isnan(my_array)
        >>> my_array = array_cleaner(my_array, indices)
    """
    _arr[_ind] = np.nan
    nans, x = _nan_helper(_arr)
    _arr[nans] = np.interp(x(nans), x(~nans), _arr[~nans])
    return _arr


def _reshape_data(data):
    data = np.delete(data, 0, axis=1)
    data = np.delete(data, 1, axis=1)
    data = np.delete(data, 2, axis=1)
    data = np.delete(data, 3, axis=1)
    data = np.delete(data, -4, axis=1)
    data = np.delete(data, -3, axis=1)
    data = np.delete(data, -2, axis=1)
    data = np.delete(data, -1, axis=1)
    return data
