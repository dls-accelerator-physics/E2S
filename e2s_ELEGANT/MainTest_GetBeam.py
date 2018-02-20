#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 14/2/2018

This function retrieves the TWISS parameters of the lattice
And produces the beam sizes/divergences 

INPUT: argument 1 :    <root>.twi, <root>.rf file produced by ELEGANT
       argument 2 :    s location in the lattice, in meters
       e.g.  python MainTest_GetBeam.py  ../e2s_LATTICES/VMX 100
OUTPUT: beam parameters from Twiss:
        sx = ...
        sy = ...
        ...

@author: MA
"""

import os
import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import sys
import subprocess



from fct_get_beam_param_from_twiss import GetBeamParam

print(os.getcwd())

betax   = 7.6832911
betay   = 4.2698328
alphax  = -0.24845980
alphay  = -0.5471616
etax    = 0
etaxp   = 0
ex0     = 1.4e-10
Sdelta0 = 0.0006701929
Cou     = 0.01
sz      = 0.001728486 

twiss = [betax,alphax,betay,alphay,etax,etaxp,ex0,Sdelta0,Cou,sz]

beam  = GetBeamParam(twiss)


print("                                   ")
print(" beam parameters from Twiss:")
print(" ----------------------------------")
print(" sx    :", beam[0])
print(" sy    :", beam[1])
print(" sxp   :", beam[2])
print(" syp   :", beam[3])
print(" dE    :", beam[4])
print(" sz    :", beam[5])
