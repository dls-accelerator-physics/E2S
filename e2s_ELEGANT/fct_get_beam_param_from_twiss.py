#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri 13 Feb 2018

@author: MA

description: uses a twiss list as input and calculates the beam sizes/divergences in return    
"""

import os
import sys
import subprocess
import numpy as np

os.system(' pwd ')

def GetBeamParam(tw):
    
    cou = tw[8]
    sz  = tw[9]
    bx  = tw[0]
    ax  = tw[1]
    gx  = (1+ax**2)/bx
    by  = tw[2]
    ay  = tw[3]
    gy  = (1+ay**2)/by
    ex  = tw[6]
    ey  = ex*cou
    hx  = tw[4]
    hxp = tw[5]
    dE  = tw[7]

    sx  = np.sqrt( bx*ex + (hx*dE)**2 )
    sy  = np.sqrt( by*ey )
    sxp = np.sqrt( gx*ex + (hxp*dE)**2  )
    syp = np.sqrt( gy*ey )

# -----------------------------------------------------------
# MA - 07032018 
# use instead the covariance matrix for the beam envelope 
# in order to introduce mixed terms 
# NOTE the term hx * hxp * (dE**2):
# this is not present in SPECTRA, determining 
# a different spectrum for cases when a non-zero dispersion
# is present (e.g. VMX-minibeta+HVF @DI)
# ----------------------------------------------------------

    sigXX   =  ex * bx + (hx * dE)**2
    sigXpXp =  ex * gx + (hxp * dE)**2
    sigXXp  = -ex * ax +  hx  * hxp * (dE**2)
    sigYY   =  ey * by                                        
    sigYpYp =  ey * gy                                       
    sigYYp  = -ey * ay   

    beam = [sx, sy, sxp, syp, dE, sz]
    mom  = [sigXX, sigXXp, sigXpXp, sigYY, sigYYp, sigYpYp] 
    return beam, mom 
