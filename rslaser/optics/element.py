import numpy as np
from pykern.pkcollections import PKDict
import srwlib
from rslaser.utils.validator import ValidatorBase
from rslaser.utils import srwl_uti_data as srwutil
from rsmath import lct as rslct


_CRYSTAL_SLICE_DEFAULTS = PKDict(
    n0=1.75,
    n2=0.001,
    length=0.2,
)

_CRYSTAL_DEFAULTS = PKDict(
        n0= _CRYSTAL_SLICE_DEFAULTS.n0,
        n2= _CRYSTAL_SLICE_DEFAULTS.n2,
        length=_CRYSTAL_SLICE_DEFAULTS.length*3,
        nslice=3,
    )

class ElementException(Exception):
    pass


class Element(ValidatorBase):
    def propagate(self, laser_pulse, prop_type='default'):
        if prop_type != 'default':
            raise ElementException(f'Non default prop_type "{prop_type}" passed to propagation')
        if not hasattr(self, '_srwc'):
            raise ElementException(f'_srwc field is expected to be set on {self}')
        for w in laser_pulse.slice:
            srwlib.srwl.PropagElecField(w.wfr,self._srwc)
        return laser_pulse


class Crystal(Element):
    """
    Args:
        params (PKDict) with fields:
            n0 (float): on axis index of refraction in crystal slice
            n2 (float): quadratic variation of index of refraction, with n(r) = n0 - 1/2 n2 r^2  [m^-2]
            length (float): total length of crystal [m]
            nslice (int): number of crystal slices
            l_scale: length scale factor for LCT propagation
    """
    _DEFAULTS = _CRYSTAL_DEFAULTS.copy()
    _INPUT_ERROR = ElementException
    def __init__(self, params=None):
        params = self._get_params(params)
        # self._validate_params(params)
        self.length = params.length
        self.nslice = params.nslice
        self.l_scale = params.l_scale
        self.slice = []
        for _ in range(self.nslice):
            self.slice.append(
                CrystalSlice(
                    PKDict(
                        n0=params.n0,
                        n2=params.n2,
                        length=params.length / params.nslice,
                        l_scale = params.l_scale
                    )
                )
            )

    def propagate(self, laser_pulse, prop_type='default'):
        # TODO (gurhar1133): should this take laser_pulse and prop_type?
        # also, should pass the same pulse through each slice and return
        # the final pulse result?
        for s in self.slice:
            laser_pulse = s.propagate(laser_pulse, prop_type)
        return laser_pulse


class CrystalSlice(Element):
    """
    This class represents a slice of a crystal in a laser cavity.

    Args:
        params (PKDict) with fields:
            length
            n0 (float): on-axis index of refraction
            n2 (float): transverse variation of index of refraction [1/m^2]
            n(r) = n0 - 0.5 n2 r^2
            l_scale: length scale factor for LCT propagation

    To be added: alpha0, alpha2 laser gain parameters

    Note: Initially, these parameters are fixed. Later we will update
    these parameters as the laser passes through.
    """

    _DEFAULTS = _CRYSTAL_SLICE_DEFAULTS.copy()
    _INPUT_ERROR = ElementException
    def __init__(self, params=None):
        params = self._get_params(params)
        # self._validate_params(params)
        self.length = params.length
        self.n0 = params.n0
        self.n2 = params.n2
        self.l_scale = params.l_scale
        # self.pop_inv = params._pop_inv

        #  Assuming wfr0 exsts, created e.g. via
        #  wfr0=createGsnSrcSRW(sigrW,propLen,pulseE,poltype,phE,sampFact,mx,my)
        #n_x = wfr0.mesh.nx  #  nr of grid points in x
        #n_y = wfr0.mesh.ny  #  nr of grid points in y
        #sig_cr_sec = np.ones((n_x, n_y), dtype=np.float32)


    def propagate(self, laser_pulse, prop_type):
        if prop_type == 'attenuate':
            raise NotImplementedError(f'{self}.propagate() with prop_type="attenuate" is not currently supported')
            # n_x = wfront.mesh.nx  #  nr of grid points in x
            # n_y = wfront.mesh.ny  #  nr of grid points in y
            # sig_cr_sec = np.ones((n_x, n_y), dtype=np.float32)
            # pop_inv = self.pop_inv
            # n0_phot = 0.0 *sig_cr_sec # incident photon density (3D), at a given transv. loc-n
            # eta = n0_phot *c_light *tau_pulse
            # gamma_degen = 1.0
            # en_gain = np.log( 1. +np.exp(sig_cr_sec *pop_inv *element.length) *(
            #             np.exp(gamma_degen *sig_cr_sec *eta) -1.0) ) /(gamma_degen *sig_cr_sec *eta)
            # return laser_pulse
        if prop_type == 'placeholder':
            raise NotImplementedError(f'{self}.propagate() with prop_type="placeholder" is not currently supported')
            # nslices = len(laser_pulse.slice)
            # for i in np.arange(nslices):
            #     print ('Pulse slice ', i+1, ' of ', nslices, ' propagated through crystal slice.')
            # return laser_pulse
            
        if prop_type == 'abcd_lct':
            nslices_pulse = len(laser_pulse.slice)
            L_cryst = self.length
            n0 = self.n0
            n2 = self.n2
            l_scale = self.l_scale
            
            phE = laser_pulse.phE
            
            ##Convert energy to wavelength
            hc_ev_um = 1.23984198   # hc [eV*um]
            phLambda = hc_ev_um / phE * 1e-6 # wavelength corresponding to phE in meters
            print("Wavelength corresponding to %g keV: %g microns" %(phE * 1e-3, phLambda / 1e-6))
            
            # calculate components of ABCD matrix corrected with wavelength and scale factor
            gamma = np.sqrt(n2/n0)
            A = np.cos(gamma*L_cryst)
            B = (1/gamma)*np.sin(gamma*L_cryst) * phLambda / (l_scale**2)
            C = -gamma*np.sin(gamma*L_cryst) / phLambda * (l_scale**2)
            D = np.cos(gamma*L_cryst)
            abcd_mat_cryst = np.array([[ A,  B  ],
                     [ C, D ]])
            print('A: %g' %A)
            print('B: %g' %B)
            print('C: %g' %C)
            print('D: %g' %D)
            
            # for i in np.arange(nslices_pulse):
            i = 0
            thisSlice = laser_pulse.slice[i]
            #print(type(thisSlice))

            # construct 2d numpy complex E_field from pulse wfr object
            wfr0 = thisSlice.wfr
            re0, re0_mesh = srwutil.calc_int_from_wfr(wfr0, _pol=0, _int_type=5, _det=None, _fname='', _pr=True)
            im0, im0_mesh = srwutil.calc_int_from_wfr(wfr0, _pol=0, _int_type=6, _det=None, _fname='', _pr=True)
            re0_2d = np.array(re0).reshape((wfr0.mesh.nx, wfr0.mesh.ny), order='C')
            im0_2d = np.array(im0).reshape((wfr0.mesh.nx, wfr0.mesh.ny), order='C')

            Etot0_2d = re0_2d + 1j*im0_2d

            xvals_slice = np.linspace(wfr0.mesh.xStart,wfr0.mesh.xFin,wfr0.mesh.nx)
            yvals_slice = np.linspace(wfr0.mesh.yStart,wfr0.mesh.yFin,wfr0.mesh.ny)

            dX = xvals_slice[1] - xvals_slice[0]                       # horizontal spacing [m]
            dX_scale = dX / l_scale
            dY = yvals_slice[1] - yvals_slice[0]                       # vertical spacing [m]
            dY_scale = dY / l_scale            

            # define input signal
            in_signal_2d = (dX_scale, dY_scale, Etot0_2d)

            # calculate 2D LCT
            dX_out, dY_out, out_signal_2d = rslct.apply_lct_2d_sep(abcd_mat_cryst, abcd_mat_cryst, in_signal_2d)                

            # extract propagated complex field and calculate corresponding x and y mesh arrays
            hx = dX_out * l_scale 
            hy = dY_out * l_scale
            sig_arr = out_signal_2d
            ny, nx = np.shape(sig_arr)
            local_xv = rslct.lct_abscissae(nx, hx)
            local_yv = rslct.lct_abscissae(ny, hy)
                
                
            return local_xv, local_yv, sig_arr
        # return laser_pulse            
            
        if prop_type == 'n0n2':
            nslices = len(laser_pulse.slice)
            L_cryst = self.length
            n0 = self.n0
            n2 = self.n2
            #n2 = 0.001

            for i in np.arange(nslices):
                thisSlice = laser_pulse.slice[i]
                #print(type(thisSlice))

                if n2 == 0:
                    #print('n2 = 0')
                    #A = 1.0
                    #B = L_cryst
                    #C = 0.0
                    #D = 1.0
                    optDrift = srwlib.SRWLOptD(L_cryst/n0)
                    propagParDrift = [0, 0, 1., 0, 0, 1., 1., 1., 1., 0, 0, 0]
                    #propagParDrift = [0, 0, 1., 0, 0, 1.1, 1.2, 1.1, 1.2, 0, 0, 0]
                    optBL = srwlib.SRWLOptC([optDrift],[propagParDrift])
                    #print("L_cryst/n0=",L_cryst/n0)
                else:
                    #print('n2 .ne. 0')
                    gamma = np.sqrt(n2/n0)
                    A = np.cos(gamma*L_cryst)
                    B = (1/gamma)*np.sin(gamma*L_cryst)
                    C = -gamma*np.sin(gamma*L_cryst)
                    D = np.cos(gamma*L_cryst)
                    f1= B/(1-A)
                    L = B
                    f2 = B/(1-D)

                    optLens1 = srwlib.SRWLOptL(f1, f1)
                    optDrift = srwlib.SRWLOptD(L)
                    optLens2 = srwlib.SRWLOptL(f2, f2)

                    propagParLens1 = [0, 0, 1., 0, 0, 1, 1, 1, 1, 0, 0, 0]
                    propagParDrift = [0, 0, 1., 0, 0, 1, 1, 1, 1, 0, 0, 0]
                    propagParLens2 = [0, 0, 1., 0, 0, 1, 1, 1, 1, 0, 0, 0]

                    optBL = srwlib.SRWLOptC([optLens1,optDrift,optLens2],[propagParLens1,propagParDrift,propagParLens2])
                    #optBL = createABCDbeamline(A,B,C,D)

                    srwlib.srwl.PropagElecField(thisSlice.wfr, optBL) # thisSlice s.b. a pointer, not a copy
                    print('Propagated pulse slice ', i+1, ' of ', nslices)
            return laser_pulse
        else:
            return super().propagate(laser_pulse, prop_type)


class Drift(Element):

    def __init__(self,length):
        self.length = length
        self._srwc = srwlib.SRWLOptC(
            [srwlib.SRWLOptD(length)],
            [[0, 0, 1., 0, 0, 1., 1., 1., 1., 0, 0, 0]],
        )


class Lens(Element):
    """
    Create lens element

    Args:
        f (float): focal length [m]
        
    Returns:
        SRW beamline element representing lens
    """

    def __init__(self,f):
        self.length = 0
        self._srwc = srwlib.SRWLOptC(
            [srwlib.SRWLOptL(f, f)],
            [[0, 0, 1., 0, 0, 1., 1., 1., 1., 0, 0, 0]]
        )

