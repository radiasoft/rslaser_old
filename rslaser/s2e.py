def main(): 
  #  Instantiate the objects corresponding to the optical line elements and 
  #  specify the propagation sequence: 
  
  dr1 = els.drift('dr1', 1.0)  # drift of length 1.0 m with a unique label 'dr1'  
  dr2 = els.drift('dr2', 1.0) 
  crystal = els.crystalSlice('cryst_slice1', 0.1)  # a single-slice crystal 
  lattice = [(dr1,'default'), (crystal,'attenuate'), (dr2,'default')] 
  
  current_position = 0.0 
  
  #  Initialize the laser pulse: 
  
  sigrW = ...
  propLen = ... 
  pulseE = ... 
  poltype = ... 
  phE = ...  # units - ? 
  #wfront = rso.wavefront.createGsnSrcSRW(sigrW,propLen,pulseE,poltype,phE,sampFact,mx,my)  # creates Gaussian wavefront in SRW 
  wfront = rso.wavefront.createGsnSrcSRW(sigrW,propLen,pulseE,poltype,phE)  # defualt values for omitted arguments 
  
  #  Propagate the pulse through the optical beamline: 
  
  for i in lattice: 
    current_elem, prop_type = i 
    wfront = current_elem.propagate(wfront, prop_type) 
    current_position += current_elem.length 

  #  Diagnostics, visualization, saving the output, etc.: 









if __name__=="__main__":
  from __future__ import division, print_function, absolute_import 
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
  import rslaser.rselements as els 
  
  import sys 
  import time 
  import math 
  
  main() 
