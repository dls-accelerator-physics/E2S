#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu 16 Aug 2018

@author: MA

description: external file describing a BL    
"""

import os
import sys
import subprocess
import numpy as np
import interpol

SRWLIB      = '/dls/physics/xph53246/source_to_beamline/SRWLIB/' # MA 12/03/2018 - repository created for pure SRWlib files 

sys.path.insert(0, SRWLIB)
from srwlib import *
from uti_plot import *

def optBL(slitDX, slitDY, Ephot):

    print "+ ------------------------------------- + "
    print "| this is Beam Line I13d - coherence      "
    print "| rev. 16/8/2018,  M. Apollonio DLS       " 
    print "+ ------------------------------------- + "

    verba = 0

    # --------
    # CRL lens
    # -------- CASE E = 11200 eV
    delta        = 2.712180e-6   #    @11.2keV (SIREPO)
    attenLen     = 12542e-6      #[m] @11.2keV (SIREPO) 
    geomApertH   = 1.1E-03       #[m] Geometrical aperture of 1D CRL in the H plane
    geomApertV   = 1.1E-03       #[m] Geometrical aperture of 1D CRL in the 
    
    rMin         = 351.0e-6
    nCRL         = 3
    wallThick    = 50E-06 #[m] wall thickness of CRL
    
    ftheo        = rMin / 2 / delta / nCRL 
    print('theoretical CRL focal length = '+str(ftheo))       
    optCRL = srwl_opt_setup_CRL(3, delta, attenLen, 1, geomApertH, geomApertV, rMin, nCRL, wallThick, 0, 0)
    
    
    # -----------------------------------------
    # Si(111) Crystal Parameters: fn of energy!
    # -----------------------------------------
    if Ephot>0:
        print " MONO - Si111 crystal parameters set for Ephot = {}".format(Ephot)
        crpar       = interpol.CrPar(Ephot/1e3, 'Si111')           # Faissal's function (note: Ephot eV --> keV, unit used in Si111.csv)
        dSpSi111    = 3.1355713563754857                           # Crystal reflecting planes d-spacing for Si(111) crystal
        psi0rSi111  = crpar[0]; psi0iSi111 = crpar[1]              # Real and imaginary parts of 0-th Fourier component of crystal polarizability
        psihrSi111  = crpar[2]; psihiSi111 = crpar[3]              # Real and imaginary parts of h-th Fourier component of crystal polarizability
        psihbrSi111 = psihrSi111; psihbiSi111 = psihiSi111         # Real and imaginary parts of -h-th Fourier component of crystal polarizability
        thickCryst  = 10.e-03                                      # Thickness of each crystal [m]
        angAsCryst  = 0                                            # Asymmetry angle of each crystal [rad]
    else: 
        Ephot =np.abs(Ephot) # if Ephot < 0 use fixed params 
        dSpSi111    = 3.1355713563754857                             # Crystal reflecting planes d-spacing for Si(111) crystal
        psi0rSi111  = -7.757245827e-6; psi0iSi111 = 9.506848329e-8   # Real and imaginary parts of 0-th Fourier component of crystal polarizability
        psihrSi111  = -4.095903022e-6; psihiSi111 = 6.637430983e-8   # Real and imaginary parts of h-th Fourier component of crystal polarizability
        psihbrSi111 = psihrSi111; psihbiSi111 = psihiSi111           # Real and imaginary parts of -h-th Fourier component of crystal polarizability
        thickCryst  = 10.e-03                                        # Thickness of each crystal [m]
        angAsCryst  = 0                                              # Asymmetry angle of each crystal [rad]
        
        
    dE_error =  0.0 # energy error on the mono
        
    # -----------------------------------------------------
    # monoC1 - 1st crystal of the I13(coh) monochromator -- for reference PI = 3.1415926535897932384626433832795 
    # -----------------------------------------------------
    optCr_1 = SRWLOptCryst(_d_sp=dSpSi111, _psi0r=psi0rSi111,
                           _psi0i=psi0iSi111, _psi_hr=psihrSi111, _psi_hi=psihiSi111, _psi_hbr=psihbrSi111, 
                           _psi_hbi=psihbiSi111,_tc=thickCryst, _ang_as=angAsCryst)
    #Find appropriate orientation of the 1st crystal and the corresponding output beam frame (in the incident beam frame):
    orientDataCr1 = optCr_1.find_orient(Ephot+dE_error, np.pi/2 ) # (GsnBm.avgPhotEn) # ,0): deflect in the v-plane (default) / 3.1415926535897932384626433832795/2): deflect in the h-plane 
    orientCr1 = orientDataCr1[0] #1st crystal orientation
    tCr1 = orientCr1[0]; nCr1 = orientCr1[2] # Tangential and Normal vectors to crystal surface
    if verba == 1:
        print('   1st crystal orientation:'); 
        print('   t=', tCr1, 's=', orientCr1[1], 'n=', nCr1)

    #Set crystal orientation:
    optCr_1.set_orient(nCr1[0], nCr1[1], nCr1[2], tCr1[0], tCr1[1])
    orientOutFrCr1 = orientDataCr1[1] #Orientation (coordinates of base vectors) of the output beam frame 
    rxCr1 = orientOutFrCr1[0]; ryCr1 = orientOutFrCr1[1]; rzCr1 = orientOutFrCr1[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
    if verba == 11:
        print('   1st crystal output beam frame:'); print('   ex=', rxCr1, 'ey=', ryCr1, 'ez=', rzCr1)
    TrM = [rxCr1, ryCr1, rzCr1] #Input/Output beam transformation matrix (for debugging)
    if verba == 1:
        uti_math.matr_print(TrM)
        print('   Beam frame transformation matrix (from the begining of opt. scheme to output of current element):')
        
    
    # -----------------------------------------------------
    # monoC2 - 2nd crystal of the I13(coh) monochromator
    # -----------------------------------------------------
    optCr_2 = SRWLOptCryst(_d_sp=dSpSi111, _psi0r=psi0rSi111,
                           _psi0i=psi0iSi111, _psi_hr=psihrSi111, _psi_hi=psihiSi111, _psi_hbr=psihbrSi111, _psi_hbi=psihbiSi111,_tc=thickCryst, _ang_as=angAsCryst)
    #Find appropriate orientation of the 2nd crystal and the corresponding output beam frame (in the incident beam frame):
    orientDataCr2 = optCr_2.find_orient(Ephot+dE_error, _ang_dif_pl= np.pi + np.pi/2) # (GsnBm.avgPhotEn) 
    orientCr2 = orientDataCr2[0] #2nd crystal orientation
    tCr2 = orientCr2[0]; nCr2 = orientCr2[2] # Tangential and Normal vectors to crystal surface
    if verba == 2:
        print('   2nd crystal orientation:'); print('   t=', tCr2, 's=', orientCr2[1], 'n=', nCr2)
    #Set crystal orientation:
    optCr_2.set_orient(nCr2[0], nCr2[1], nCr2[2], tCr2[0], tCr2[1])
    orientOutFrCr2 = orientDataCr2[1] #Orientation (coordinates of base vectors) of the output beam frame 
    rxCr2 = orientOutFrCr2[0]; ryCr2 = orientOutFrCr2[1]; rzCr2 = orientOutFrCr2[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
    if verba == 22:
        print('   2nd crystal output beam frame:'); print('   ex=', rxCr2, 'ey=', ryCr2, 'ez=', rzCr2)
    TrM = uti_math.matr_prod(TrM, [rxCr2, ryCr2, rzCr2]) #Input/Output beam transformation matrix (for debugging)
    if verba == 1:
        print('   Beam frame transformation matrix (from the begining of opt. scheme to output of current element):')
        uti_math.matr_print(TrM)
        print('   After the two crystals of DCM, the transformation matrix should be close to the unit matrix.')
    
    # -----------------------------------------------------
    # monoC3 - 3rd crystal of the I13(coh) monochromator
    # -----------------------------------------------------
    optCr_3 = SRWLOptCryst(_d_sp=dSpSi111, _psi0r=psi0rSi111,
                           _psi0i=psi0iSi111, _psi_hr=psihrSi111, _psi_hi=psihiSi111, _psi_hbr=psihbrSi111, _psi_hbi=psihbiSi111,_tc=thickCryst, _ang_as=angAsCryst)
    #Find appropriate orientation of the 3rd crystal and the corresponding output beam frame (in the incident beam frame):
    orientDataCr3 = optCr_3.find_orient(Ephot+dE_error, _ang_dif_pl= + np.pi/2) # (GsnBm.avgPhotEn) 
    orientCr3 = orientDataCr3[0] #3rd crystal orientation
    tCr3 = orientCr3[0]; nCr3 = orientCr3[2] # Tangential and Normal vectors to crystal surface
    if verba == 3:
        print('   3rd crystal orientation:'); print('   t=', tCr3, 's=', orientCr3[1], 'n=', nCr3)
    #Set crystal orientation:
    optCr_3.set_orient(nCr3[0], nCr3[1], nCr3[2], tCr3[0], tCr3[1])
    orientOutFrCr3 = orientDataCr3[1] #Orientation (coordinates of base vectors) of the output beam frame 
    rxCr3 = orientOutFrCr3[0]; ryCr3 = orientOutFrCr3[1]; rzCr3 = orientOutFrCr3[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
    if verba == 33:
        print('   3rd crystal output beam frame:'); print('   ex=', rxCr3, 'ey=', ryCr3, 'ez=', rzCr3)
    TrM = uti_math.matr_prod(TrM, [rxCr3, ryCr3, rzCr3]) #Input/Output beam transformation matrix (for debugging)
    if verba == 1:
        print('   Beam frame transformation matrix (from the begining of opt. scheme to output of current element):')
        uti_math.matr_print(TrM)
    
    # -----------------------------------------------------
    # monoC4 - 4th crystal of the I13(coh) monochromator
    # -----------------------------------------------------
    optCr_4 = SRWLOptCryst(_d_sp=dSpSi111, _psi0r=psi0rSi111,
                           _psi0i=psi0iSi111, _psi_hr=psihrSi111, _psi_hi=psihiSi111, _psi_hbr=psihbrSi111, _psi_hbi=psihbiSi111,_tc=thickCryst, _ang_as=angAsCryst)
    #Find appropriate orientation of the 4th crystal and the corresponding output beam frame (in the incident beam frame):
    orientDataCr4 = optCr_4.find_orient(Ephot+dE_error, _ang_dif_pl= np.pi + np.pi/2) # (GsnBm.avgPhotEn) 
    orientCr4 = orientDataCr4[0] #4th crystal orientation
    tCr4 = orientCr4[0]; nCr4 = orientCr4[2] # Tangential and Normal vectors to crystal surface
    if verba == 4:
        print('   4th crystal orientation:'); print('   t=', tCr4, 's=', orientCr4[1], 'n=', nCr4)
    #Set crystal orientation:
    optCr_4.set_orient(nCr4[0], nCr4[1], nCr4[2], tCr4[0], tCr4[1])
    orientOutFrCr4 = orientDataCr4[1] #Orientation (coordinates of base vectors) of the output beam frame 
    rxCr4 = orientOutFrCr4[0]; ryCr4 = orientOutFrCr4[1]; rzCr4 = orientOutFrCr4[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
    if verba == 44:
        print('   4th crystal output beam frame:'); print('   ex=', rxCr4, 'ey=', ryCr4, 'ez=', rzCr4)
    TrM = uti_math.matr_prod(TrM, [rxCr4, ryCr4, rzCr4]) #Input/Output beam transformation matrix (for debugging)
    if verba == 1:
        print('   Beam frame transformation matrix (from the begining of opt. scheme to output of current element):')
        uti_math.matr_print(TrM)
        print('   After the two crystals of DCM, the transformation matrix should be close to the unit matrix.')
    
    # -----------------------------------------------------
    # M1 - planar mirror deflecting in the Horizontal Plane
    # -----------------------------------------------------
    tG, sG  = 0.2, 0.2  
    grazG   = 10e-3  # 10mrad  
    optM_1  = SRWLOptMirPl(_size_tang=tG, _size_sag=sG, _ap_shape='r',
                           _nvx=cos(grazG), _nvy=0, _nvz=-sin(grazG), _tvx=sin(grazG), _tvy=0)
    
    # --------------------------------------------------------
    # KB1 - elliptical mirror deflecting in the Vertical Plane
    # --------------------------------------------------------
    tG, sG  = 0.3, 0.3  
    grazG   = 3e-3 # 1. * np.pi / 180 
    optKB_1   = SRWLOptMirEl(_p=30.9, _q=9.1, _ang_graz=grazG, _r_sag=1.e+23,
                             _size_tang=1, _size_sag=1,
                             _nvx=0, _nvy=cos(grazG), _nvz=-sin(grazG), _tvx=0, _tvy=-sin(grazG))
    
    # ----------------------------------------------------------
    # KB2 - elliptical mirror deflecting in the Horizontal Plane
    # ----------------------------------------------------------
    tG, sG  = 0.3, 0.3  
    grazG   = 3e-3 # 1. * np.pi / 180; 
    optKB_2   = SRWLOptMirEl(_p=33.1, _q=6.9, _ang_graz=grazG, _r_sag=1.e+23,
                             _size_tang=1, _size_sag=1,
                             _nvx=cos(grazG), _nvy=0, _nvz=-sin(grazG), _tvx=-sin(grazG), _tvy=0)
    
    # ------
    # Drifts
    # ------
    
    D1            = 9.999
    optDrift_1    = SRWLOptD(D1)                              # 10m - 1mm drift 
    D2            = 4.7                                       # 4.7m drift
    optDrift_2    = SRWLOptD(D2)
    D3            = 187.1-D2 
    optDrift_3    = SRWLOptD(D3)                              # 182.5m drift from CRL (position of mono) 
    D4            = 1
    optDrift_4    = SRWLOptD(D4)
    D5            = 10
    optDrift_5    = SRWLOptD(D5)
    DKB           = 2.2
    optDrift_KB   = SRWLOptD(DKB)
    DKBSam        = 6.9-.8     #  .1:46/40   0:42/38   -.1: 38/36   -.2:34/34     -.3:30/32 
                               # -.4:26/30 -.5:22/29   -.6:19/27    -.7:15/25     -.8:11/23 
                               # -.9:8/22  -1.:5.2/21 -1.1:4.9/19  -1.2:7.3/17.3 -1.3:10.6/16.6
                               # -1.4:   -2: 14.3/16.7 
                               # -0.8m selected, brute figures: sx=11um/sy=23.6um, however the central spot is ~10.5um/11.2um 
    optDrift_KB_Sam   = SRWLOptD(DKBSam)
    
    
    # ----------------------
    # Propagation Parameters 
    # ----------------------
    
    propagParApert =  [0, 0, 1., 0, 0, 1.05, 1.0, 1.05, 1., 0, 0, 0]
    propagParLens  =  [0, 0, 1., 0, 0, 1.05, 1.0, 1.05, 1., 0, 0, 0]
    propagParDrift =  [0, 0, 1., 1, 0, 1.2, 1.1, 1.2, 1., 0, 0, 0]
    propagParPM    =  [0, 0, 1., 1, 0, 1.05, 1.1, 1.05, 1., 0, 0, 0]
    propagParCryst =  [0, 0, 1., 1, 0, 1.05, 1.1, 1.05, 1., 0, 0, 0]
    propagParKB    =  [0, 0, 1., 0, 0, 1.05, 1.0, 1.05, 1., 0, 0, 0]
    
    # --------------------------------------------------------------------------
    # Lists of Optical Elements (oe) and Propagation Parameters (pp) are defined
    # --------------------------------------------------------------------------
    
    s = []  # position along the line from the entry aperture

    oe=[]; pp=[] 
    bl = []
    bl.append("-||- ")                                                                                                                                                   # the optical elements
    oe.append( SRWLOptA('r', 'a', slitDX*1e-6, slitDY*1e-6) )                                                                                                       # oe(1):  aperture
    pp.append(propagParApert)
    s.append(0)

    bl.append( " --- " )
    oe.append( optDrift_1 )                                                                                                                                         # oe(2):  drift
    pp.append(propagParDrift)
    s.append(oe[-1].L)

    bl.append(" -()- ")
    oe.append( optCRL )                                                                                                                                             # oe(3):  CRL       
    pp.append(propagParLens)
    s.append(0)

    bl.append(" --- ")
    oe.append( optDrift_2 )                                                                                                                                        # oe(4):  drift
    pp.append(propagParDrift)
    s.append(oe[-1].L)


    bl.append(" -/- ")
    oe.append( optM_1 )                                                                                                                                             # oe(5):  plane mirror
    pp.append(propagParPM)
    s.append(0)
    
    bl.append(" --- ")
    oe.append( optDrift_3 )                                                                                                                                         # oe(6):  drift
    pp.append(propagParDrift)
    s.append(oe[-1].L)


    bl.append(" -## ")
    oe.append( optCr_1 )                                                                                                                                            # oe(7):  mono-element1
    pp.append(propagParCryst)
    s.append(0)

    bl.append(" ##- ")
    oe.append( optCr_2 )                                                                                                                                            # oe(8):  mono-element2
    pp.append(propagParCryst)
    s.append(0)

    bl.append(" --- ")
    oe.append( optDrift_4 )                                                                                                                                         # oe(9):  drift
    pp.append(propagParDrift)
    s.append(oe[-1].L)

    bl.append(" -## ")
    oe.append( optCr_3 )                                                                                                                                            # oe(10): mono-element3
    pp.append(propagParCryst)
    s.append(0)

    bl.append(" ##- ")
    oe.append( optCr_4 )                                                                                                                                            # oe(11): mono-element4
    pp.append(propagParCryst)
    s.append(0)

    bl.append(" --- ")
    oe.append( optDrift_5 )                                                                                                                                         # oe(12): drift
    pp.append(propagParDrift)
    s.append(oe[-1].L)


    bl.append( " -)- ")    
    oe.append( optKB_1 )                                                                                                                                            # oe(13): KB1 
    pp.append(propagParKB)
    s.append(0)

    bl.append(" --- ")
    oe.append( optDrift_KB )                                                                                                                                        # oe(14): drift KB 
    pp.append(propagParDrift)
    s.append(oe[-1].L)

    bl.append(" -(- ")
    oe.append( optKB_2 )                                                                                                                                            # oe(15): KB2 
    pp.append(propagParKB)
    s.append(0)

    bl.append(" --- ")
    oe.append( optDrift_KB_Sam )                                                                                                                                    # oe(16): drift KB-sample  
    pp.append(propagParDrift)
    s.append(oe[-1].L)

    bl.append(" --|S  ")
    bl = ''.join(bl)
    print "{}".format(bl)
    cs = np.cumsum(s)
    print("  ".join(str(int(np.ceil(x))) for x in cs))

    oe_pp =[]; oe_pp.append(oe); oe_pp.append(pp)
        
    return oe_pp


