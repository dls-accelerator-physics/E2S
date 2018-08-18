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


SRWLIB      = '/dls/physics/xph53246/source_to_beamline/SRWLIB/' # MA 12/03/2018 - repository created for pure SRWlib files 

sys.path.insert(0, SRWLIB)
from srwlib import *
from uti_plot import *

def optBL(slitDX, slitDY, Ephot):

    print "+ ------------------------------------- + "
    print "| this is Beam Line I04                   "
    print "| rev. 16/8/2018,  M. Apollonio DLS       " 
    print "+ ------------------------------------- + "

    verba = 0 

    optApert      = SRWLOptA('r', 'a', slitDX*1e-6, slitDY*1e-6) # aperture
    
    # ------
    # Drifts
    # ------
    
    D1            = 13.6
    optDrift_1    = SRWLOptD(D1)                             
    D2            = 2.2                                      
    optDrift_KB    = SRWLOptD(D2)
    D3                 = 6.9 
    optDrift_KB_Sam    = SRWLOptD(D3)                              

    # --------------------------------------------------------
    # KB1 - elliptical mirror deflecting in the Vertical Plane
    # --------------------------------------------------------
    tG, sG  = 0.3, 0.3
    grazG   = 3e-3 # 1. * np.pi / 180;
    optKB_1 = SRWLOptMirEl(_p=30.9, _q=9.1, _ang_graz=grazG, _r_sag=1.e+23,
                                _size_tang=1, _size_sag=1,
                                _nvx=0, _nvy=cos(grazG), _nvz=-sin(grazG), _tvx=0, _tvy=-sin(grazG))

    # ----------------------------------------------------------
    # KB2 - elliptical mirror deflecting in the Horizontal Plane
    # ----------------------------------------------------------
    tG, sG  = 0.3, 0.3
    grazG   = 3e-3 # 1. * np.pi / 180;
    optKB_2 = SRWLOptMirEl(_p=33.1, _q=6.9, _ang_graz=grazG, _r_sag=1.e+23,
                            _size_tang=1, _size_sag=1,
                            _nvx=cos(grazG), _nvy=0, _nvz=-sin(grazG), _tvx=-sin(grazG), _tvy=0)
    
    # ----------------------
    # Propagation Parameters 
    # ----------------------
    
    propagParApert =  [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1., 0, 0, 0]
    propagParDrift =  [0, 0, 1., 1, 0, 1.0, 1.1, 1.0, 1., 0, 0, 0]
    propagParKB    =  [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1., 0, 0, 0]
    
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

    bl.append(" -(- ")
    oe.append( optKB_1 )                                                                                                                                             # oe(3):  CRL       
    pp.append(propagParKB)
    s.append(0)

    bl.append(" --- ")
    oe.append( optDrift_KB )                                                                                                                                        # oe(4):  drift
    pp.append(propagParDrift)
    s.append(oe[-1].L)


    bl.append(" -)- ")
    oe.append( optKB_2 )                                                                                                                                             # oe(5):  plane mirror
    pp.append(propagParKB)
    s.append(0)
    
    bl.append(" --- ")
    oe.append( optDrift_KB_Sam )                                                                                                                                         # oe(6):  drift
    pp.append(propagParDrift)
    s.append(oe[-1].L)

    bl.append(" --|S  ")
    bl = ''.join(bl)
    print "{}".format(bl)
    cs = np.cumsum(s)
    print("  ".join(str(int(np.ceil(x))) for x in cs))

    oe_pp =[]; oe_pp.append(oe); oe_pp.append(pp)
        
    return oe_pp


