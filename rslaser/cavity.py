import math
import numpy as np
import srwlib
from array import array
from pykern.pkcollections import PKDict

class Element:
    def propagate(self,laser_pulse):
        x = len(laser_pulse._slice)//2
        for w in laser_pulse._slice:
            srwlib.srwl.PropagElecField(w._wfr,self._srwc)
            if w.slice_index != x:
                continue
            (sx,sy) = rmsWavefrontIntensity(w._wfr)
            laser_pulse._sxvals.append(sx)
            laser_pulse._syvals.append(sy)
            #print(f'Element RMS sizes:sx={sx} sy={sy}')

class LaserPulseSlice:
    """
    This class represents a longitudinal slice in a laser pulse.
    There will be a number of wavefronts each with different wavelengths (energy).
    The slice is composed of an SRW wavefront object, which is defined here:
    https://github.com/ochubar/SRW/blob/master/env/work/srw_python/srwlib.py#L2048
    """
    def __init__(self,nslice,slice_index,sigrW=0.00043698412731784714,propLen=15,sig_s=0.1,pulseE=0.001,poltype=1,phE=1.55,sampFact=5,mx=0,my=0,**_ignore_kwargs):
        """
        #nslice: number of slices of laser pulse
        #slice_index: index of slice
        #sigrW: beam size at waist [m]
        #propLen: propagation length [m] required by SRW to create numerical Gaussian
        #propLen=15,
        #sig_s RMS pulse length [m]
        #pulseE: energy per pulse [J]
        #poltype: polarization type (0=linear horizontal, 1=linear vertical, 2=linear 45 deg, 3=linear 135 deg, 4=circular right, 5=circular left, 6=total)
        #phE: photon energy [eV]
        #sampFact: sampling factor to increase mesh density
        """
        #print([sigrW,propLen,pulseE,poltype])
        self.slice_index = slice_index
        constConvRad = 1.23984186e-06/(4*3.1415926536)  ##conversion from energy to 1/wavelength
        rmsAngDiv = constConvRad/(phE*sigrW)             ##RMS angular divergence [rad]
        sigrL=math.sqrt(sigrW**2+(propLen*rmsAngDiv)**2)  ##required RMS size to produce requested RMS beam size after propagation by propLen
    
        
        #***********Gaussian Beam Source
        GsnBm = srwlib.SRWLGsnBm() #Gaussian Beam structure (just parameters)
        GsnBm.x = 0 #Transverse Positions of Gaussian Beam Center at Waist [m]
        GsnBm.y = 0
        numsig = 3. #Number of sigma values to track. Total range is 2*numsig*sig_s
        ds = 2*numsig*sig_s/nslice
        self._pulse_pos = -numsig*sig_s+slice_index*ds
        GsnBm.z = propLen + self._pulse_pos #Longitudinal Position of Waist [m]
        GsnBm.xp = 0 #Average Angles of Gaussian Beam at Waist [rad]
        GsnBm.yp = 0
        GsnBm.avgPhotEn = phE #Photon Energy [eV]
        GsnBm.pulseEn = pulseE*np.exp(-self._pulse_pos**2/(2*sig_s**2)) #Energy per Pulse [J] - to be corrected
        GsnBm.repRate = 1 #Rep. Rate [Hz] - to be corrected
        GsnBm.polar = poltype #1- linear horizontal?
        GsnBm.sigX = sigrW #Horiz. RMS size at Waist [m]
        GsnBm.sigY = GsnBm.sigX #Vert. RMS size at Waist [m]

        GsnBm.sigT = 10e-15 #Pulse duration [s] (not used?)
        GsnBm.mx = mx #Transverse Gauss-Hermite Mode Orders
        GsnBm.my = my

        #***********Initial Wavefront
        wfr = srwlib.SRWLWfr() #Initial Electric Field Wavefront
        wfr.allocate(1, 1000, 1000) #Numbers of points vs Photon Energy (1), Horizontal and Vertical Positions (dummy)
        wfr.mesh.zStart = 0.0 #Longitudinal Position [m] at which initial Electric Field has to be calculated, i.e. the position of the first optical element
        wfr.mesh.eStart = GsnBm.avgPhotEn #Initial Photon Energy [eV]
        wfr.mesh.eFin = GsnBm.avgPhotEn #Final Photon Energy [eV]

        wfr.unitElFld = 1 #Electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)

        distSrc = wfr.mesh.zStart - GsnBm.z
        #Horizontal and Vertical Position Range for the Initial Wavefront calculation
        #can be used to simulate the First Aperture (of M1)
        #firstHorAp = 8.*rmsAngDiv*distSrc #[m]
        xAp = 8.*sigrL
        yAp = xAp #[m]

        wfr.mesh.xStart = -0.5*xAp #Initial Horizontal Position [m]
        wfr.mesh.xFin = 0.5*xAp #Final Horizontal Position [m]
        wfr.mesh.yStart = -0.5*yAp #Initial Vertical Position [m]
        wfr.mesh.yFin = 0.5*yAp #Final Vertical Position [m]

        sampFactNxNyForProp = sampFact #sampling factor for adjusting nx, ny (effective if > 0)
        arPrecPar = [sampFactNxNyForProp]

        srwlib.srwl.CalcElecFieldGaussian(wfr, GsnBm, arPrecPar)

        ##Beamline to propagate to waist

        optDriftW=srwlib.SRWLOptD(propLen)
        propagParDrift = [0, 0, 1., 0, 0, 1.1, 1.2, 1.1, 1.2, 0, 0, 0]
        optBLW = srwlib.SRWLOptC([optDriftW],[propagParDrift])
        #wfrW=deepcopy(wfr)
        srwlib.srwl.PropagElecField(wfr, optBLW)

        self._wfr = wfr
    
            
class LaserPulse:
    """
    A laserPulse is a collection of laserSlices.
    """
    def __init__(self,length,wavelength,nslice,**kwargs):
        self._slice = []
        for i in range(nslice):
            #Creation of laser slices i=0...nslice-1
            self._slice.append(LaserPulseSlice(nslice,i,**kwargs))
            #(sx,sy) = rmsWavefrontIntensity(self._slice[-1]._wfr)
            #print(f'Laser pulse RMS sizes:sx={sx} sy={sy}')
        self._sxvals = []
        self._syvals = []
        
    def rmsvals(self):
        sx = []
        sy = []
        s = []
        for sl in self._slice:
            (sigx,sigy) = rmsWavefrontIntensity(sl._wfr)
            sx.append(sigx)
            sy.append(sigy)
            s.append(sl._pulse_pos)
            
        return(sx,sy,s)
    
    def intensity_vals(self):
        intensity = []
        s = []
        for sl in self._slice:
            slint = maxWavefrontIntensity(sl._wfr)
            intensity.append(slint)
            s.append(sl._pulse_pos)
            
        return(intensity,s)
    
    def slice_wfr(self,slice_index):
        return self(slice_index)._slice._wfr
    
class CrystalSlice:
    """
    This class represents a slice of a crystal in a laser cavity.
    Parameters:
    #length
    #n0: on-axis index of refraction
    #n2: transverse variation of index of refraction 
    n(r) = n0 - 0.5 n2 r^2
    #To be added: alpha0, alpha2 laser gain parameters
    Initially, these parameters are fixed. Later we will update
    these parameters as the laser passes through.
    """
    def __init__(self,n0,n2):
        self.n0=n0
        self.n2=n2
        
class Crystal(Element):
    def __init__(self,n0,n2,L_cryst):
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
            optDrift=srwlib.SRWLOptD(Lc/2)
            propagParDrift = [0, 0, 1., 0, 0, 1., 1., 1., 1., 0, 0, 0]
            #propagParDrift = [0, 0, 1., 0, 0, 1.1, 1.2, 1.1, 1.2, 0, 0, 0]
            return srwlib.SRWLOptC([optDrift],[propagParDrift])
        
        if n2==0:
            self._srwc=_createDriftBL(2*L_cryst) #Note that this drift function divides length by 2
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
        self._srwc = srwlib.SRWLOptC(
            [srwlib.SRWLOptD(length)],
            [[0, 0, 1., 0, 0, 1., 1., 1., 1., 0, 0, 0]],
        )
            
class Lens(Element):
    """
    #Create lens element
    #f: focal length
    """
    def __init__(self,f):
        self._srwc = srwlib.SRWLOptC(
            [srwlib.SRWLOptL(f, f)],
            [[0, 0, 1., 0, 0, 1., 1., 1., 1., 0, 0, 0]]
        )
       
class LaserCavity:
    """
    create laser cavity
    """
    def __init__(self,**kwargs):
        k=PKDict(kwargs).pksetdefault(
            n0 = 1.75,
            n2 = 0.001,
            L_half_cryst=0.2,
            laser_pulse_length=0.1,
            wavelength=800e-9,
            nslice=11,
            drift_right_length=0.5,
            drift_left_length=0.5,
            lens_left_focal_length=0.2,
            lens_right_focal_length=0.2,
            
            sigrW=0.000437,
            propLen=15,
            sig_s=0.1,
            pulseE=0.001,
            poltype=1,
            phE=1.55,
            sampFact=5,
            mx=0,
            my=0
        )
        
        self.laser_pulse = LaserPulse(length = k.laser_pulse_length,**k)
        self.crystal_right = Crystal(k.n0,k.n2,k.L_half_cryst)
       
        self.crystal_left = Crystal(k.n0,k.n2,k.L_half_cryst)
        self.drift_right = Drift(k.drift_right_length)
        self.drift_left = Drift(k.drift_left_length)
        self.lens_right = Lens(k.lens_right_focal_length)
        self.lens_left  = Lens(k.lens_left_focal_length)
    
        
        
    def propagate(self,num_cycles):
        l = self.laser_pulse
        l._sxvals = []
        l._syvals = []
        for n in range(num_cycles):
            self.crystal_right.propagate(l)
            self.drift_right.propagate(l)
            self.lens_right.propagate(l)
            self.drift_right.propagate(l)
            self.crystal_right.propagate(l)
            
            self.crystal_left.propagate(l)
            self.drift_left.propagate(l)
            self.lens_left.propagate(l)
            self.drift_left.propagate(l)
            self.crystal_left.propagate(l)
        
        return(l._sxvals,l._syvals)
            
def rmsWavefrontIntensity(wfr):
    """
    #Compute rms values from a wavefront object
    """
    IntensityArray2D = array('f', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" array to take 2D intensity data
    srwlib.srwl.CalcIntFromElecField(IntensityArray2D, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0) #extracts intensity
    ##Reshaping electric field data from flat to 2D array
    IntensityArray2D = np.array(IntensityArray2D).reshape((wfr.mesh.nx, wfr.mesh.ny), order='C')
    xvals=np.linspace(wfr.mesh.xStart,wfr.mesh.xFin,wfr.mesh.nx)
    yvals=np.linspace(wfr.mesh.yStart,wfr.mesh.yFin,wfr.mesh.ny)
    (sx,sy) = rmsIntensity(IntensityArray2D,xvals,yvals)
    return sx, sy

def rmsIntensity(IntArray,xvals,yvals):
    """
    Compute rms values in x and y from array 
    #IntArray is a 2D array representation of a function
    #xvals represents the horizontal coordinates
    #yvals represents the vertical coordinates
    """
    datax=np.sum(IntArray,axis=1) 
    datay=np.sum(IntArray,axis=0)
    sxsq=sum(datax*xvals*xvals)/sum(datax) 
    xavg=sum(datax*xvals)/sum(datax)
    sx=math.sqrt(sxsq-xavg*xavg)

    sysq=sum(datay*yvals*yvals)/sum(datay) 
    yavg=sum(datay*yvals)/sum(datay)
    sy=math.sqrt(sysq-yavg*yavg)
    return sx, sy

def maxWavefrontIntensity(wfr):
    """
    Compute maximum value of wavefront intensity
    """
    IntensityArray2D = array('f', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" array to take 2D intensity data
    srwlib.srwl.CalcIntFromElecField(IntensityArray2D, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0) #extracts intensity
    return(np.max(IntensityArray2D))
    