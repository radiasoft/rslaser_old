import math
import srwlib
from srwlib import *

def createGsnSrcSRW(sigrW,propLen,pulseE,poltype,phE=10e3,sampFact=15,mx=0,my=0):
    
    """
    #sigrW: beam size at waist [m]
    #propLen: propagation length [m] required by SRW to create numerical Gaussian
    #pulseE: energy per pulse [J]
    #poltype: polarization type (0=linear horizontal, 1=linear vertical, 2=linear 45 deg, 3=linear 135 deg, 4=circular right, 5=circular left, 6=total)
    #phE: photon energy [eV]
    #sampFact: sampling factor to increase mesh density
    """
    
    constConvRad = 1.23984186e-06/(4*3.1415926536)  ##conversion from energy to 1/wavelength
    rmsAngDiv = constConvRad/(phE*sigrW)             ##RMS angular divergence [rad]
    sigrL=math.sqrt(sigrW**2+(propLen*rmsAngDiv)**2)  ##required RMS size to produce requested RMS beam size after propagation by propLen
    
        
    #***********Gaussian Beam Source
    GsnBm = SRWLGsnBm() #Gaussian Beam structure (just parameters)
    GsnBm.x = 0 #Transverse Positions of Gaussian Beam Center at Waist [m]
    GsnBm.y = 0
    GsnBm.z = propLen #Longitudinal Position of Waist [m]
    GsnBm.xp = 0 #Average Angles of Gaussian Beam at Waist [rad]
    GsnBm.yp = 0
    GsnBm.avgPhotEn = phE #Photon Energy [eV]
    GsnBm.pulseEn = pulseE #Energy per Pulse [J] - to be corrected
    GsnBm.repRate = 1 #Rep. Rate [Hz] - to be corrected
    GsnBm.polar = poltype #1- linear horizontal?
    GsnBm.sigX = sigrW #Horiz. RMS size at Waist [m]
    GsnBm.sigY = GsnBm.sigX #Vert. RMS size at Waist [m]

    GsnBm.sigT = 10e-15 #Pulse duration [s] (not used?)
    GsnBm.mx = mx #Transverse Gauss-Hermite Mode Orders
    GsnBm.my = my

    #***********Initial Wavefront
    wfr = SRWLWfr() #Initial Electric Field Wavefront
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
    
    srwl.CalcElecFieldGaussian(wfr, GsnBm, arPrecPar)
    
    ##Beamline to propagate to waist
    
    optDriftW=SRWLOptD(propLen)
    propagParDrift = [0, 0, 1., 0, 0, 1.1, 1.2, 1.1, 1.2, 0, 0, 0]
    optBLW = SRWLOptC([optDriftW],[propagParDrift])
    #wfrW=deepcopy(wfr)
    srwl.PropagElecField(wfr, optBLW)
    
    return wfr

def createBL1to1(L,dfof=0):

    """
    ##Define beamline geometric variables.
    #L: drift length before and after lens
    #dfof: focal length variation factor (=0 for no variation; can be positive or negative)
    """
    
    
    ##Drift lengths between elements beginning with source to 1st crystal and ending with last crystal to start of undulator.
    
    
    ##focal length in meters
    f=(L/2)*(1+dfof)
    
    #Lens    
    optLens = SRWLOptL(f, f) 
    #Drift spaces
    optDrift1=SRWLOptD(L)
    optDrift2=SRWLOptD(L)
       

    #***********Wavefront Propagation Parameters:
    #[0]: Auto-Resize (1) or not (0) Before propagation
    #[1]: Auto-Resize (1) or not (0) After propagation
    #[2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
    #[3] Type of the propagator:
        #0 - Standard - Fresnel (it uses two FFTs);
        #1 - Quadratic Term - with semi-analytical treatment of the quadratic (leading) phase terms (it uses two FFTs);
        #2 - Quadratic Term - Special - special case;
        #3 - From Waist - good for propagation from "waist" over a large distance (it uses one FFT);
        #4 - To Waist - good for propagation to a "waist" (e.g. some 2D focus of an optical system) over some distance (it uses one FFT).
    #[4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
    #[5]: Horizontal Range modification factor at Resizing (1. means no modification)
    #[6]: Horizontal Resolution modification factor at Resizing
    #[7]: Vertical Range modification factor at Resizing
    #[8]: Vertical Resolution modification factor at Resizing
    #[9]: Type of wavefront Shift before Resizing (not yet implemented)
    #[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
    #[11]: New Vertical wavefront Center position after Shift (not yet implemented)
    
    #propagParLens = [0, 0, 1., 0, 0, 1., 1.5, 1., 1.5, 0, 0, 0]
    #propagParDrift = [0, 0, 1., 0, 0, 1., 1., 1., 1., 0, 0, 0]
    
    propagParLens = [0, 0, 1., 0, 0, 1.4, 2., 1.4, 2., 0, 0, 0]
    propagParDrift = [0, 0, 1., 0, 0, 1.1, 1.2, 1.1, 1.2, 0, 0, 0]
    
    ##Beamline consruction
    optBL1to1 = SRWLOptC([optDrift1,optLens,optDrift2],[propagParDrift,propagParLens,propagParDrift])
    
    return optBL1to1