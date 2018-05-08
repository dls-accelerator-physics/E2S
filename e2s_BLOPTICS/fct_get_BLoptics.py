#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun 18 Mar 2018

@author: MA

description: returns a BLopt container to be used in a SRW_*.py file    
"""

import os
import sys
import subprocess
import numpy as np


SRWLIB      = '/dls/physics/xph53246/source_to_beamline/SRWLIB/' # MA 12/03/2018 - repository created for pure SRWlib files 

sys.path.insert(0, SRWLIB)
from srwlib import *
from uti_plot import *

os.system(' pwd ')


#***********Optical Elements and Propagation Parameters
#fx = 1e+23 #Focal Length in Horizontal plane
#fy = 1e+23 # 19.0939 #Focal Length in Vertical plane
#optLens = SRWLOptL(fx, fy) #Ideal Lens

#delta = 4.3712962E-06 #Refractive index decrement of Be at 8830 eV
#attenLen = 6946.13E-06 #[m] Attenuation length of Be at 8830 eV
#geomApertF = 1E-03 #[m] Geometrical aparture of 1D CRL in the Focusing plane
#geomApertNF = 1E-03 #[m] Geometrical aparture of 1D CRL in the plane where there is no focusing
#rMin = 0.5E-03 #[m] radius at tip of parabola of CRL
#nCRL = 3
#wallThick = 50E-06 #[m] wall thickness of CRL

#optCRL = srwl_opt_setup_CRL(2, delta, attenLen, 1, geomApertNF, geomApertF, rMin, nCRL, wallThick, 0, 0) #1D CRL
#print('Saving CRL transmission data to files (for viewing/debugging)...', end='')
#optTrIntCRL = optCRL.get_data(2, 3)
#srwl_uti_save_intens_ascii(optTrIntCRL, optCRL.mesh, os.path.join(os.getcwd(), strExDataFolderName, strOpTrFileName), 0, ['', 'Horizontal Position', 'Vertical Position', 'Intensity Transmission'], _arUnits=['', 'm', 'm', 'r.u.'])

#optPathDifCRL = optCRL.get_data(3, 3)
#srwl_uti_save_intens_ascii(optPathDifCRL, optCRL.mesh, os.path.join(os.getcwd(), strExDataFolderName, strOpPathDifFileName), 0, ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'], _arUnits=['', 'm', 'm', 'm'])
#print('done')

#####optDrift =  SRWLOptD(20.) # SRWLOptD(0.1) #Drift space
#####optApert = SRWLOptA('r', 'a', slitDX*1e-6, slitDY*1e-6) #Aperture

#optDrift = SRWLOptD(12.9) #Drift space

# propagParApert = [0, 0, 1., 0, 0, 1.5, 1.0, 1.1, 8., 0, 0, 0]
#####propagParApert = [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1., 0, 0, 0]
#####propagParLens =  [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1., 0, 0, 0]
#####propagParDrift = [0, 0, 1., 1, 0, 1.0, 1.1, 1.0, 1., 0, 0, 0]

#Wavefront Propagation Parameters:
#[0]: Auto-Resize (1) or not (0) Before propagation
#[1]: Auto-Resize (1) or not (0) After propagation
#[2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
#[3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
#[4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
#[5]: Horizontal Range modification factor at Resizing (1. means no modification)
#[6]: Horizontal Resolution modification factor at Resizing
#[7]: Vertical Range modification factor at Resizing
#[8]: Vertical Resolution modification factor at Resizing
#[9]: Type of wavefront Shift before Resizing (not yet implemented)
#[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
#[11]: New Vertical wavefront Center position after Shift (not yet implemented)
#optBL = SRWLOptC([optApert, optCRL, optDrift], [propagParApert, propagParLens, propagParDrift]) #"Beamline" - Container of Optical Elements (together with the 
#optBL = SRWLOptC([optApert, optLens, optDrift], [propagParApert, propagParLens, propagParDrift]) #"Beamline" - Container of Optical Elements (together with the corresponding wavefront propagation instructions)

###optBL = SRWLOptC([optApert, optDrift], [propagParApert, propagParDrift]) #"Beamline" - Container of Optical Elements (together with the corresponding wavefront propagation instructions)



def DefineBLOptics(BLname,slitDX,slitDY):

#
#   I13-coherence branch
#
    if BLname == 'I13d_ENTRY':
        # ----------------------------------------------------------------------
        # - entry slits
        # ----------------------------------------------------------------------
        
        #
        # A) define the optical elements
        #
        optApert   = SRWLOptA('r', 'a', slitDX*1e-6, slitDY*1e-6) # aperture     
        #optDrift_1 = SRWLOptD(1e-3)                                # drift
        # B) define the propagation parameters
        #
        propagParApert = [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1., 0, 0, 0]
        #propagParDrift = [0, 0, 1., 1, 0, 1.0, 1.1, 1.0, 1., 0, 0, 0]

        optBL = SRWLOptC([optApert], [propagParApert]) #

    if BLname == 'I13d':
        # ----------------------------------------------------------------------
        # - entry slits
        # - drift.1
        # - CRL lens 
        # - drift.2
        # - image plane
        # ----------------------------------------------------------------------
        
        #
        # A) define the optical elements
        #
        optApert   = SRWLOptA('r', 'a', slitDX*1e-6, slitDY*1e-6) # aperture     
        optDrift_1 = SRWLOptD(10.)                                # drift
        optDrift_2 = SRWLOptD(4.7)                                # drift
        optDrift_3 = SRWLOptD(5.3)                                # drift
        optDrift_4 = SRWLOptD(187.1)                              # drift
        #
        # CRL lens parameter + definition
        #
        # CASE E = 8247 (TEST)
#       delta       = 5.011848e-6 #RID @ 8247 eV    ###  4.3712962E-06 #Refractive index decrement of Be at 8830 eV
#       attenLen    = 5749E-06 ### AL@ 8247 eV      ###  6946.13E-06 #[m] Attenuation length of Be at 8830 eV
#        delta       = 2.22e-4       #refractive index decrement of Be @ 1243 eV
#        attenLen    = 17e-6       #[m] Attenuation length of Be at 1243 eV
        delta       = 9.758814e-7 #@18.6 keV  5.011848e-6 #RID @ 8247 eV    ###  4.3712962E-06 #Refractive index decrement of Be at 8830 eV
        attenLen    = 26828e-6    #AL@18.6keV  ## 5749E-06 = AL@ 8247 eV      ###  6946.13E-06 #[m] Attenuation length of Be at 8830 eV
#       
        geomApertH  = 2E-03 #[m] Geometrical aperture of 1D CRL in the H plane
        geomApertV  = 2E-03 #[m] Geometrical aperture of 1D CRL in the V plane
        rMin        = 273.5e-6 ### 1.0525E-03 ###  0.5E-03 #[m] radius at tip of parabola of CRL
        nCRL        = 7 ###  3
        wallThick = 50E-06 #[m] wall thickness of CRL
       
        optCRL = srwl_opt_setup_CRL(3, delta, attenLen, 1, geomApertH, geomApertV, rMin, nCRL, wallThick, 0, 0) #1D CRL
        #
        # M1 - planar mirror deflecting in the Horizontal Plane
        #
        tG, sG  = 0.3, 0.3  
        grazG   = 1. * np.pi / 180; 
        optM1   = SRWLOptMirPl(_size_tang=tG, _size_sag=sG, _ap_shape='r',
                    _nvx=cos(grazG), _nvy=0, _nvz=-sin(grazG), _tvx=sin(grazG), _tvy=0)

        #
        # B) define the propagation parameters
        #
        propagParApert = [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1., 0, 0, 0]
        propagParLens =  [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1., 0, 0, 0]
        propagParDrift = [0, 0, 1., 1, 0, 1.0, 1.1, 1.0, 1., 0, 0, 0]
        propagParMirP  = [0, 0, 1., 1, 0, 1.0, 1.1, 1.0, 1., 0, 0, 0]

#        optBL = SRWLOptC([optApert, optDrift_1, optDrift_2], [propagParApert, propagParDrift, propagParDrift]) #
#        optBL = SRWLOptC([optApert, optDrift_1, optCRL, optDrift_2], [propagParApert, propagParDrift, propagParLens, propagParDrift]) #
        optBL = SRWLOptC([optApert, optDrift_1, optCRL, optDrift_2, optM1, optDrift_3, optDrift_4], 
         [propagParApert, propagParDrift, propagParLens, propagParDrift, propagParMirP, propagParDrift, propagParDrift]) #







#
#   I20-SCAnning branch
#
    if BLname == 'I20_SCA':
        optDrift = SRWLOptD(15.)
        optApert = SRWLOptA('r', 'a', slitDX*1e-6, slitDY*1e-6) #Aperture       
        propagParApert = [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1., 0, 0, 0]
        propagParLens =  [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1., 0, 0, 0]
        propagParDrift = [0, 0, 1., 1, 0, 1.0, 1.1, 1.0, 1., 0, 0, 0]
        
        optBL = SRWLOptC([optApert, optDrift], [propagParApert, propagParDrift]) #

#
#   I20-EDE branch
#
    if BLname == 'I20_EDE':    
        optDrift = SRWLOptD(15.)
        optApert = SRWLOptA('r', 'a', slitDX*1e-6, slitDY*1e-6) #Aperture       
        propagParApert = [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1., 0, 0, 0]
        propagParLens =  [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1., 0, 0, 0]
        propagParDrift = [0, 0, 1., 1, 0, 1.0, 1.1, 1.0, 1., 0, 0, 0]
        
        optBL = SRWLOptC([optApert, optDrift], [propagParApert, propagParDrift]) #



    return optBL
