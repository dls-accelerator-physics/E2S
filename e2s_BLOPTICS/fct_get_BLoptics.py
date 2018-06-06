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

    if BLname == 'I13d_test':
        # ----------------------------------------------------------------------
        # - entry slits
        # - drift.1
        # - CRL
        # - drift.1
        # ----------------------------------------------------------------------
        #
        # A) define the optical elements
        #
        optApertoo    = SRWLOptA('r', 'a', slitDX*1e-3, slitDY*1e-3) # huge aperture, to see the beamspot at entrance 
        optApert      = SRWLOptA('r', 'a', slitDX*1e-6, slitDY*1e-6) # aperture     
        optDrift_mini = SRWLOptD(0.00001)                            # infinitesimal drift
        optDrift_1    = SRWLOptD(9.999)                              # 10m - 1mm drift 
        optDrift_2    = SRWLOptD(21.55)                              # 21.55m drift 
        optDrift_22   = SRWLOptD(4.7)                                # 21.55m drift 
        optDrift_3    = SRWLOptD(100.0-4.7)                          # 100m drift from CRL 
        optDrift_33   = SRWLOptD(187.1-4.7)                          # 182.5m drift from CRL (position of mono) 
        optDrift_4    = SRWLOptD(10)                                 # 10m drift 
        optDrift_44   = SRWLOptD(1)                                  # 1m drift 
        optDrift_Sam  = SRWLOptD(10)                                 # 10m drift to sample 

        # CASE E = 11200 eV
        delta        = 2.712180e-6 #    @11.2keV (SIREPO)
        attenLen     = 12542e-6    #[m] @11.2keV (SIREPO) 
        geomApertH   = 1.1E-03       #[m] Geometrical aperture of 1D CRL in the H plane
        geomApertV   = 1.1E-03       #[m] Geometrical aperture of 1D CRL in the 

        rMin         = 351.0e-6
        nCRL         = 3
        wallThick    = 50E-06 #[m] wall thickness of CRL
        
        ftheo        = rMin / 2 / delta / nCRL 
        print('theoretical CRL focal length = '+str(ftheo))       
        optCRL = srwl_opt_setup_CRL(3, delta, attenLen, 1, geomApertH, geomApertV, rMin, nCRL, wallThick, 0, 0) #1D CRL
        
        # -----------------------------------------------------
        # M1 - planar mirror deflecting in the Horizontal Plane
        # -----------------------------------------------------
        tG, sG  = 0.2, 0.2  
        grazG   = 10e-3 # 10mrad  #1. * np.pi / 180; 
        optM1   = SRWLOptMirPl(_size_tang=tG, _size_sag=sG, _ap_shape='r',
                    _nvx=cos(grazG), _nvy=0, _nvz=-sin(grazG), _tvx=sin(grazG), _tvy=0)
        print('>>>>>>>>>>>>>>>>>>>> _nvx = ', cos(grazG))
        #Si(111) Crystal Constants:
        dSpSi111 = 3.1355713563754857 # Crystal reflecting planes d-spacing for Si(111) crystal
        psi0rSi111 = -7.757245827e-6; psi0iSi111 = 9.506848329e-8 #Real and imaginary parts of 0-th Fourier component of crystal polarizability
        psihrSi111 = -4.095903022e-6; psihiSi111 = 6.637430983e-8 #Real and imaginary parts of h-th Fourier component of crystal polarizability
        psihbrSi111 = psihrSi111; psihbiSi111 = psihiSi111 #Real and imaginary parts of -h-th Fourier component of crystal polarizability
        thickCryst = 10.e-03 #0.5e-03 #Thickness of each crystal [m]
        angAsCryst = 0 #Asymmetry angle of each crystal [rad]


        dE_error =  0.0 # energy error on the mono
        # -----------------------------------------------------
        # monoC1 - 1st crystal of the I13(coh) monochromator
        # -----------------------------------------------------
        optCr1 = SRWLOptCryst(_d_sp=dSpSi111, _psi0r=psi0rSi111,
            _psi0i=psi0iSi111, _psi_hr=psihrSi111, _psi_hi=psihiSi111, _psi_hbr=psihbrSi111, 
            _psi_hbi=psihbiSi111,_tc=thickCryst, _ang_as=angAsCryst)
        #Find appropriate orientation of the 1st crystal and the corresponding output beam frame (in the incident beam frame):
        orientDataCr1 = optCr1.find_orient(11208+dE_error, 3.1415926535897932384626433832795/2 ) # (GsnBm.avgPhotEn) # ,0): deflect in the v-plane (default) / 3.1415926535897932384626433832795/2): deflect in the h-plane 
        orientCr1 = orientDataCr1[0] #1st crystal orientation
        tCr1 = orientCr1[0]; nCr1 = orientCr1[2] # Tangential and Normal vectors to crystal surface
        print('   1st crystal orientation:'); print('   t=', tCr1, 's=', orientCr1[1], 'n=', nCr1)
        #Set crystal orientation:
        optCr1.set_orient(nCr1[0], nCr1[1], nCr1[2], tCr1[0], tCr1[1])
        orientOutFrCr1 = orientDataCr1[1] #Orientation (coordinates of base vectors) of the output beam frame 
        rxCr1 = orientOutFrCr1[0]; ryCr1 = orientOutFrCr1[1]; rzCr1 = orientOutFrCr1[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
        print('   1st crystal output beam frame:'); print('   ex=', rxCr1, 'ey=', ryCr1, 'ez=', rzCr1)
        TrM = [rxCr1, ryCr1, rzCr1] #Input/Output beam transformation matrix (for debugging)
        print('   Beam frame transformation matrix (from the begining of opt. scheme to output of current element):')
        uti_math.matr_print(TrM)

        # -----------------------------------------------------
        # monoC2 - 2nd crystal of the I13(coh) monochromator
        # -----------------------------------------------------
        optCr2 = SRWLOptCryst(_d_sp=dSpSi111, _psi0r=psi0rSi111,
                     _psi0i=psi0iSi111, _psi_hr=psihrSi111, _psi_hi=psihiSi111, _psi_hbr=psihbrSi111, _psi_hbi=psihbiSi111,_tc=thickCryst, _ang_as=angAsCryst)
        #Find appropriate orientation of the 2nd crystal and the corresponding output beam frame (in the incident beam frame):
        orientDataCr2 = optCr2.find_orient(11208+dE_error, _ang_dif_pl= 3.1415926535897932384626433832795 + 3.1415926535897932384626433832795/2) # (GsnBm.avgPhotEn) 
        orientCr2 = orientDataCr2[0] #2nd crystal orientation
        tCr2 = orientCr2[0]; nCr2 = orientCr2[2] # Tangential and Normal vectors to crystal surface
        print('   2nd crystal orientation:'); print('   t=', tCr2, 's=', orientCr2[1], 'n=', nCr2)
        #Set crystal orientation:
        optCr2.set_orient(nCr2[0], nCr2[1], nCr2[2], tCr2[0], tCr2[1])
        orientOutFrCr2 = orientDataCr2[1] #Orientation (coordinates of base vectors) of the output beam frame 
        rxCr2 = orientOutFrCr2[0]; ryCr2 = orientOutFrCr2[1]; rzCr2 = orientOutFrCr2[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
        print('   2nd crystal output beam frame:'); print('   ex=', rxCr2, 'ey=', ryCr2, 'ez=', rzCr2)
        TrM = uti_math.matr_prod(TrM, [rxCr2, ryCr2, rzCr2]) #Input/Output beam transformation matrix (for debugging)
        print('   Beam frame transformation matrix (from the begining of opt. scheme to output of current element):')
        uti_math.matr_print(TrM)
        print('   After the two crystals of DCM, the transformation matrix should be close to the unit matrix.')

        # -----------------------------------------------------
        # monoC3 - 3rd crystal of the I13(coh) monochromator
        # -----------------------------------------------------
        optCr3 = SRWLOptCryst(_d_sp=dSpSi111, _psi0r=psi0rSi111,
                     _psi0i=psi0iSi111, _psi_hr=psihrSi111, _psi_hi=psihiSi111, _psi_hbr=psihbrSi111, _psi_hbi=psihbiSi111,_tc=thickCryst, _ang_as=angAsCryst)
        #Find appropriate orientation of the 2nd crystal and the corresponding output beam frame (in the incident beam frame):
        orientDataCr3 = optCr3.find_orient(11208+dE_error, _ang_dif_pl= + 3.1415926535897932384626433832795/2) # (GsnBm.avgPhotEn) 
        orientCr3 = orientDataCr3[0] #2nd crystal orientation
        tCr3 = orientCr3[0]; nCr3 = orientCr3[2] # Tangential and Normal vectors to crystal surface
        print('   3rd crystal orientation:'); print('   t=', tCr3, 's=', orientCr3[1], 'n=', nCr3)
        #Set crystal orientation:
        optCr3.set_orient(nCr3[0], nCr3[1], nCr3[2], tCr3[0], tCr3[1])
        orientOutFrCr3 = orientDataCr3[1] #Orientation (coordinates of base vectors) of the output beam frame 
        rxCr3 = orientOutFrCr3[0]; ryCr3 = orientOutFrCr3[1]; rzCr3 = orientOutFrCr3[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
        print('   3rd crystal output beam frame:'); print('   ex=', rxCr3, 'ey=', ryCr3, 'ez=', rzCr3)
        TrM = uti_math.matr_prod(TrM, [rxCr3, ryCr3, rzCr3]) #Input/Output beam transformation matrix (for debugging)
        print('   Beam frame transformation matrix (from the begining of opt. scheme to output of current element):')
        uti_math.matr_print(TrM)

        # -----------------------------------------------------
        # monoC4 - 4th crystal of the I13(coh) monochromator
        # -----------------------------------------------------
        optCr4 = SRWLOptCryst(_d_sp=dSpSi111, _psi0r=psi0rSi111,
                     _psi0i=psi0iSi111, _psi_hr=psihrSi111, _psi_hi=psihiSi111, _psi_hbr=psihbrSi111, _psi_hbi=psihbiSi111,_tc=thickCryst, _ang_as=angAsCryst)
        #Find appropriate orientation of the 2nd crystal and the corresponding output beam frame (in the incident beam frame):
        orientDataCr4 = optCr4.find_orient(11208+dE_error, _ang_dif_pl= 3.1415926535897932384626433832795 + 3.1415926535897932384626433832795/2) # (GsnBm.avgPhotEn) 
        orientCr4 = orientDataCr4[0] #2nd crystal orientation
        tCr4 = orientCr4[0]; nCr4 = orientCr4[2] # Tangential and Normal vectors to crystal surface
        print('   4th crystal orientation:'); print('   t=', tCr4, 's=', orientCr4[1], 'n=', nCr4)
        #Set crystal orientation:
        optCr4.set_orient(nCr4[0], nCr4[1], nCr4[2], tCr4[0], tCr4[1])
        orientOutFrCr4 = orientDataCr4[1] #Orientation (coordinates of base vectors) of the output beam frame 
        rxCr4 = orientOutFrCr4[0]; ryCr4 = orientOutFrCr4[1]; rzCr4 = orientOutFrCr4[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
        print('   4th crystal output beam frame:'); print('   ex=', rxCr4, 'ey=', ryCr4, 'ez=', rzCr4)
        TrM = uti_math.matr_prod(TrM, [rxCr4, ryCr4, rzCr4]) #Input/Output beam transformation matrix (for debugging)
        print('   Beam frame transformation matrix (from the begining of opt. scheme to output of current element):')
        uti_math.matr_print(TrM)
        print('   After the two crystals of DCM, the transformation matrix should be close to the unit matrix.')


        # ------------------------------------
        # B) define the propagation parameters
        # ------------------------------------
        propagParApert =  [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1., 0, 0, 0]
        propagParLens  =  [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1., 0, 0, 0]
        propagParDrift =  [0, 0, 1., 1, 0, 1.0, 1.1, 1.0, 1., 0, 0, 0]
        propagParM1    =  [0, 0, 1., 1, 0, 1.0, 1.1, 1.0, 1., 0, 0, 0]
        propagParCryst =  [0, 0, 1., 1, 0, 1.0, 1.1, 1.0, 1., 0, 0, 0]

        #optBL = SRWLOptC([optApertoo, optDrift_mini],   # TEST-0 
        #                 [propagParApert, propagParDrift])               # beam after a very large FE slit 
 

        #optBL = SRWLOptC([optApert, optDrift_mini],   # TEST-1 passed 
        #                 [propagParApert, propagParDrift])             # beam after FE slit 

        optBL = SRWLOptC([optApert, optDrift_1, optCRL, optDrift_3],              # TEST-2 passed   
                         [propagParApert, propagParDrift, propagParLens, propagParDrift])         # FE slit + CRL 

        # optBL = SRWLOptC([optApert, optDrift_1, optCRL, optDrift_22, optM1, optDrift_3],   # TEST-3 passed  
        #                 [propagParApert, propagParDrift, propagParLens, propagParDrift])                   # FE slit + CRL + M1 


        #optBL = SRWLOptC([optApert, optDrift_1, optCRL, optDrift_22, optM1, optDrift_33, optCr1, optDrift_4],    
        #                 [propagParApert, propagParDrift, propagParLens, propagParDrift, propagParM1, propagParDrift, propagParCryst, propagParDrift])  # FE slit + CRL + M1 + monoC1 
        

        #optBL = SRWLOptC([optApert, optDrift_1, optCRL, optDrift_22, optM1, optDrift_33, optCr1, optCr2, optDrift_4],     # TEST-4  passed 
        #                 [propagParApert, propagParDrift, propagParLens, propagParDrift, propagParM1, propagParDrift, propagParCryst, propagParCryst, propagParDrift])  # FE slit + CRL + M1 + monoC1 + monoC2 

#        optBL = SRWLOptC([optApert, optDrift_1, optCRL, optDrift_22, optM1, optDrift_33, 
#                          optCr1, optCr2, optDrift_44,
#                          optCr3, optDrift_Sam],                                                          # TEST-5 passed 
#                         [propagParApert, propagParDrift, propagParLens, propagParDrift, propagParM1, propagParDrift, 
#                          propagParCryst, propagParCryst, propagParDrift, 
#                          propagParCryst, propagParDrift])                                                # FE slit + CRL + M1 + monoC1 + monoC2 + monoC3 

#        optBL = SRWLOptC([optApert, optDrift_1, optCRL, optDrift_22, optM1, optDrift_33, 
#                          optCr1, optCr2, optDrift_44,
#                          optCr3, optCr4, optDrift_Sam],                                                          # TEST-5 passed 
#                         [propagParApert, propagParDrift, propagParLens, propagParDrift, propagParM1, propagParDrift, 
#                          propagParCryst, propagParCryst, propagParDrift, 
#                          propagParCryst, propagParCryst, propagParDrift])                                        # FE slit + CRL + M1 + monoC1 + monoC2 + monoC3 + monoC4
        
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

#       delta       = 9.758814e-7 #@18.6 keV  5.011848e-6 #RID @ 8247 eV    ###  4.3712962E-06 #Refractive index decrement of Be at 8830 eV
#       attenLen    = 26828e-6    #AL@18.6keV  ## 5749E-06 = AL@ 8247 eV      ###  6946.13E-06 #[m] Attenuation length of Be at 8830 eV
#       geomApertH  = 2E-03 #[m] Geometrical aperture of 1D CRL in the H plane
#       geomApertV  = 2E-03 #[m] Geometrical aperture of 1D CRL in the V plane
#       rMin        = 273.5e-6 ### 1.0525E-03 ###  0.5E-03 #[m] radius at tip of parabola of CRL
#       nCRL        = 7 ###  3
#       wallThick = 50E-06 #[m] wall thickness of CRL

        # CASE E = 11200 eV
        delta        = 2.712180e-6 #    @11.2keV (SIREPO)
        attenLen     = 12542e-6    #[m] @11.2keV (SIREPO) 
        geomApertH   = 1.1E-03       #[m] Geometrical aperture of 1D CRL in the H plane
        geomApertV   = 1.1E-03       #[m] Geometrical aperture of 1D CRL in the 

        rMin         = 400.0e-6
        nCRL         = 3
        wallThick    = 50E-06 #[m] wall thickness of CRL
        
        ftheo        = rMin / 2 / delta / nCRL 
        print('theoretical CRL focal length = '+str(ftheo))

#       
       
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


    if BLname == 'I04':
        # ----------------------------------------------------------------------
        # - entry slits
        # - drift.1
        # - KB1 
        # - drift.2
        # - KB2
        # - image plane
        # ----------------------------------------------------------------------
        
        #
        # A) define the optical elements
        #
        optApert   = SRWLOptA('r', 'a', slitDX*1e-6, slitDY*1e-6) # aperture     
        optDrift_1 = SRWLOptD(13.6)                               # drift (17.6)
        optDrift_2 = SRWLOptD(2.2)                                # drift (1.1)
        optDrift_3 = SRWLOptD(6.9)                                # drift (7.9)
        
        #
        # KB1 - elliptical mirror deflecting in the Vertical Plane
        #
        tG, sG  = 0.3, 0.3  
        grazG   = 3e-3 # 1. * np.pi / 180; 
        #optKB1   = SRWLOptMirEl(_p=34.9, _q=9, _ang_graz=grazG, _r_sag=1.e+23,
        #                        _size_tang=1, _size_sag=1,
        #                        _nvx=0, _nvy=cos(grazG), _nvz=-sin(grazG), _tvx=0, _tvy=-sin(grazG))
        optKB1   = SRWLOptMirEl(_p=30.9, _q=9.1, _ang_graz=grazG, _r_sag=1.e+23,
                                _size_tang=1, _size_sag=1,
                                _nvx=0, _nvy=cos(grazG), _nvz=-sin(grazG), _tvx=0, _tvy=-sin(grazG))
        #
        # KB2 - elliptical mirror deflecting in the Horizontal Plane
        #
        #tG, sG  = 0.3, 0.3  
        #grazG   = 3e-3 # 1. * np.pi / 180; 
        #optKB2   = SRWLOptMirEl(_p=36., _q=7.9, _ang_graz=grazG, _r_sag=1.e+23,
        #                        _size_tang=1, _size_sag=1,
        #                        _nvx=cos(grazG), _nvy=0, _nvz=-sin(grazG), _tvx=-sin(grazG), _tvy=0)
        tG, sG  = 0.3, 0.3  
        grazG   = 3e-3 # 1. * np.pi / 180; 
        optKB2   = SRWLOptMirEl(_p=33.1, _q=6.9, _ang_graz=grazG, _r_sag=1.e+23,
                                _size_tang=1, _size_sag=1,
                                _nvx=cos(grazG), _nvy=0, _nvz=-sin(grazG), _tvx=-sin(grazG), _tvy=0)
        #
        # B) define the propagation parameters
        #
        propagParApert = [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1., 0, 0, 0]
        propagParKB    = [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1., 0, 0, 0]
        propagParDrift = [0, 0, 1., 1, 0, 1.0, 1.1, 1.0, 1., 0, 0, 0]
        

#        optBL = SRWLOptC([optApert, optDrift_1, optDrift_2], [propagParApert, propagParDrift, propagParDrift]) #
#        optBL = SRWLOptC([optApert, optDrift_1, optCRL, optDrift_2], [propagParApert, propagParDrift, propagParLens, propagParDrift]) #
       
        optBL = SRWLOptC([optApert, optDrift_1, optKB1, optDrift_2,optKB2, optDrift_3],
                         [propagParApert, propagParDrift, propagParKB,propagParDrift, propagParKB, propagParDrift])

 # optKB2, optDrift_3], 
  #        [propagParApert, propagParDrift, propagParKB, propagParDrift, propagParKB, propagParDrift]) #

#        optBL = SRWLOptC([optApert, optDrift_1],
#                         [propagParApert, propagParDrift]) #





        

    return optBL
