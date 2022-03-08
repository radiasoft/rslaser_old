from __future__ import division, print_function, absolute_import

def main():
    #  Instantiate the objects corresponding to the optical line elements and
    #  specify the propagation sequence:

    dr1 = elsdr.drift('dr1', 0.2)  # drift of length 1.0 m with a unique label 'dr1'
    dr2 = elsdr.drift('dr2', 0.02)
    crystal = elscs.crystalSlice('cryst_slice1', 0.1)  # a single-slice crystal
    #lattice = [(dr1,'default'), (crystal,'placeholder'), (dr2,'default')]
    lattice = [(dr1,'default'), (crystal,'abcd'), (dr2,'default')]

    current_position = 0.0
    print('Initial z position of the pulse:', current_position, 'm')

  #  Initialize the laser pulse:

    pupa = PKDict(
        slice_params=PKDict(

        )
    )
    pupa.slice_params.propLen = 15.0  # [m]
    pupa.slice_params.sigrW = 0.000186  # 0.000436984 #  radial rms size at the waist [m]
    pupa.d_to_w = 0.10  #  [m] distance from the initial pulse location to the loc-n of the beam waist, > 0 if converging
    pupa.slice_params.pulseE = 0.001  # pulse energy [J]
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
    import rslaser.rscrystal as rscr
    import rslaser.rsoptics as rso
    import rslaser.rspulse as rsp
    import rslaser.rspulse.pulse as plsdv
    import rslaser.elements.drift as elsdr
    import rslaser.elements.crystalSlice as elscs

    import scipy
    import scipy.constants as const

    import sys
    import time
    import math

    main()
