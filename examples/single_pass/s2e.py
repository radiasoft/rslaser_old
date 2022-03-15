from __future__ import division, print_function, absolute_import
from pykern.pkcollections import PKDict

def main():
    #  Instantiate the objects corresponding to the optical line elements and
    #  specify the propagation sequence:

    dr1 = elsdr.Drift(0.2)  # drift of length 1.0 m with a unique label 'dr1'
    dr2 = elsdr.Drift(0.02)
    crystal = elsdr.CrystalSlice('cryst_slice1', 0.1)  # a single-slice crystal
    #lattice = [(dr1,'default'), (crystal,'placeholder'), (dr2,'default')]
    lattice = [(dr1,'default'), (crystal,'abcd'), (dr2,'default')]

    current_position = 0.0
    print('Initial z position of the pulse:', current_position, 'm')

  #  Initialize the laser pulse:
    _PHE_DEFAULT = const.h * const.c / 1e-6
    _Z_WAIST_DEFAULT = 0
    _Z_CENTER_DEFAULT = 0
    _LASER_PULSE_SLICE_DEFAULTS = PKDict(
        sigrW=0.000186,
        # propLen=15,
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
    #  Pulse parameters:
    pupa = _LASER_PULSE_DEFAULTS.copy()
    pupa.d_to_w = 0.10  #  [m] distance from the initial pulse location to the loc-n of the beam waist, > 0 if converging
    pupa.slice_params.poltype = 1  #  0 = linear horizontal, 1 = linear vertical, ...
    pupa.phE = 1.55  # Wavefront energy [eV]. 1.55 eV is 800 nm wavelength (seed laser); 532 nm for the pump laser
    pupa.chirp = 0.0
    pupa.nslice = 2  #  the number of slices the pulse is divided into

    #wfront = rso.wavefront.createGsnSrcSRW(sigrW,propLen,pulseE,poltype,phE,sampFact,mx,my)  # creates Gaussian wavefront in SRW
    #wfront = rso.wavefront.createGsnSrcSRW(sigrW,propLen,pulseE,poltype,phE)  # defualt values for omitted arguments
    thisPulse = plsdv.LaserPulse(pupa)

    #  Propagate the pulse through the optical beamline:

    for i in lattice:
        current_elem, prop_type = i
    #print (current_elem, prop_type)
        thisPulse = current_elem.propagate(thisPulse, prop_type)
        current_position += current_elem.length
        print('Current position in the beamline:', current_position, ' m')

    #  Diagnostics, visualization, saving the output, etc.:

if __name__=="__main__":
    #from __future__ import division, print_function, absolute_import
    import numpy as np
    import matplotlib.pyplot as plt

    import srwlib
    from pykern import pkcli
    from array import array
    from pykern.pkcollections import PKDict

    import rslaser.rscavity as rscav
    import rslaser.crystal as rscr
    import rslaser.optics as rso
    import rslaser.pulse as rsp
    import rslaser.pulse.pulse as plsdv
    import rslaser.optics.element as elsdr

    import scipy
    import scipy.constants as const

    import sys
    import time
    import math

    main()
