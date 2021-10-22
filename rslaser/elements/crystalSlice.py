import math
import numpy as np
import srwlib
from array import array
from pykern.pkcollections import PKDict

import rslaser.rsoptics as rso
#from rslaser.rsoptics.element import *

c_light = 299792458.0  # m/s, in vacuum 

class crystalSlice:
  """
  This class represents a slice of a crystal in a laser cavity.
  Parameters:
  #label: a unique tag labeling a physical beamline element 
  #length
  #n0: on-axis index of refraction
  #n2: transverse variation of index of refraction
  n(r) = n0 - 0.5 n2 r^2
  #pop_inv:  population inversion in the pumped crystal  
  #To be added: alpha0, alpha2 laser gain parameters
  Initially, these parameters are fixed. Later we will update
  these parameters as the laser passes through.
  """
  def __init__(self, _label, _length, _n0=1.0, _n2=0.0, _pop_inv=1.0): 
    self.label = _label 
    self.length = _length 
    self.n0 = _n0
    self.n2 = _n2
    self.pop_inv = _pop_inv 
  
  #  Assuming wfr0 exsts, created e.g. via  
  #  wfr0=createGsnSrcSRW(sigrW,propLen,pulseE,poltype,phE,sampFact,mx,my) 
  #n_x = wfr0.mesh.nx  #  nr of grid points in x 
  #n_y = wfr0.mesh.ny  #  nr of grid points in y 
  #sig_cr_sec = np.ones((n_x, n_y), dtype=np.float32) 
  
  def propagate(self, wfront, prop_type): 
    #if prop_type == 'attenuate': 
    n_x = wfront.mesh.nx  #  nr of grid points in x 
    n_y = wfront.mesh.ny  #  nr of grid points in y 
    sig_cr_sec = np.ones((n_x, n_y), dtype=np.float32)
    pop_inv = self.pop_inv
    n0_phot = 0.0 *sig_cr_sec # incident photon density (3D), at a given transv. loc-n 
    eta = n0_phot *c_light *tau_pulse
    gamma_degen = 1.0
    en_gain = np.log( 1. +np.exp(sig_cr_sec *pop_inv *element.length) *(
              np.exp(gamma_degen *sig_cr_sec *eta) -1.0) ) /(gamma_degen *sig_cr_sec *eta)


