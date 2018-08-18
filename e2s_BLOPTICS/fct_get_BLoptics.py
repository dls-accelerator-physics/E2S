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
# importing the external functions defining the beamlines 
import optBL_I13d_test
import optBL_I04

SRWLIB      = '/dls/physics/xph53246/source_to_beamline/SRWLIB/' # MA 12/03/2018 - repository created for pure SRWlib files 



sys.path.insert(0, SRWLIB)
from srwlib import *
from uti_plot import *

os.system(' pwd ')



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
#optBL = SRWLOptC([optApert, optLens, optDrift], [propagParApert, propagParLens, propagParDrift]) #"Beamline" - Container of Optical Elements (together with the corresponding wavefront propagation instructions)



def DefineBLOptics(BLname,slitDX,slitDY):

    if BLname == 'I100':
        oe_pp = optBL_I100 # returns a container of Optical Elements and Propagation Parameters
        
    if BLname == 'I04':
        Ephot = 13841 # to be changed .... 
        oe_pp = optBL_I04.optBL(slitDX, slitDY, Ephot)  # define the oe/pp container 
            
    if BLname == 'I13d_test':            
        Ephot = 11209 # to be changed .... 
        oe_pp = optBL_I13d_test.optBL(slitDX, slitDY, Ephot)  # define the oe/pp container 
                

    optBL = SRWLOptC(oe_pp[0],oe_pp[1])   # define the actual optical BL object 
        
    return optBL

