import math
import numpy as np
from array import array
from pykern.pkcollections import PKDict
import srwlib

class Element:
    def propagate(self,laser_pulse):

class Crystal(Element):
    def __init__(self,n0,n2,L_cryst):
        self.length = L_cryst

    def propagate(self,laser_pulse):
        for w in laser_pulse._slice:
            srwlib.srwl.PropagElecField(w._wfr,self._srwc)

    def _createABCDbeamline(A,B,C,D):
        """
        #Use decomposition of ABCD matrix into kick-drift-kick Pei-Huang 2017 (https://arxiv.org/abs/1709.06222)
        #Construct corresponding SRW beamline container object
        #A,B,C,D are 2x2 matrix components.
        """

        f1= B/(1-A)
        L = B
        f2 = B/(1-D)

        optLens1 = srwlib.SRWLOptL(f1, f1)
        optDrift= srwlib.SRWLOptD(L)
        optLens2 = srwlib.SRWLOptL(f2, f2)

        propagParLens1 = [0, 0, 1., 0, 0, 1, 1, 1, 1, 0, 0, 0]
        propagParDrift = [0, 0, 1., 0, 0, 1, 1, 1, 1, 0, 0, 0]
        propagParLens2 = [0, 0, 1., 0, 0, 1, 1, 1, 1, 0, 0, 0]

        return srwlib.SRWLOptC([optLens1,optDrift,optLens2],[propagParLens1,propagParDrift,propagParLens2])

        def _createDriftBL(Lc):
            optDrift=srwlib.SRWLOptD(Lc)
            propagParDrift = [0, 0, 1., 0, 0, 1., 1., 1., 1., 0, 0, 0]
            #propagParDrift = [0, 0, 1., 0, 0, 1.1, 1.2, 1.1, 1.2, 0, 0, 0]
            return srwlib.SRWLOptC([optDrift],[propagParDrift])

        if n2==0:
            self._srwc=_createDriftBL(L_cryst)
            #print("L_cryst/n0=",L_cryst/n0)
        else:
            gamma = np.sqrt(n2/n0)
            A = np.cos(gamma*L_cryst)
            B = (1/(gamma))*np.sin(gamma*L_cryst)
            C = -gamma*np.sin(gamma*L_cryst)
            D = np.cos(gamma*L_cryst)
            self._srwc=_createABCDbeamline(A,B,C,D)

class Drift(Element):
    def __init__(self,length):
        self.length = length
        self._srwc = srwlib.SRWLOptC(
            [srwlib.SRWLOptD(length)],
            [[0, 0, 1., 0, 0, 1., 1., 1., 1., 0, 0, 0]],
        )

    def propagate(self,laser_pulse):
        for w in laser_pulse._slice:
            srwlib.srwl.PropagElecField(w._wfr,self._srwc)

class Lens(Element):
    """
    #Create lens element
    #f: focal length
    """
    def __init__(self,f):
        self.length = 0
        self._srwc = srwlib.SRWLOptC(
            [srwlib.SRWLOptL(f, f)],
            [[0, 0, 1., 0, 0, 1., 1., 1., 1., 0, 0, 0]]
        )

    def propagate(self,laser_pulse):
        for w in laser_pulse._slice:
            srwlib.srwl.PropagElecField(w._wfr,self._srwc)
