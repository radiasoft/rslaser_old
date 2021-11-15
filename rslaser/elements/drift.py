import math
import numpy as np
from array import array
from pykern.pkcollections import PKDict
import srwlib 

class drift():
  def __init__(self,_label,_length): 
    self.label = _label
    self.length = _length
    self._srwc = srwlib.SRWLOptC(
    [srwlib.SRWLOptD(_length)],
    [[0, 0, 1., 0, 0, 1., 1., 1., 1., 0, 0, 0]],
    )

  def propagate(self, laser_pulse, prop_type): 
    if prop_type == 'default': 
      nslices = len(laser_pulse._slice) 
      for i in np.arange(nslices):
        thisSlice = laser_pulse._slice[i] 
        srwlib.srwl.PropagElecField(thisSlice._wfr, self._srwc) 
    return laser_pulse 
