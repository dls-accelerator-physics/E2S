#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 13:05:29 2018

This function retrieves the TWISS parameters of the lattice

INPUT: argument 1 :    *.twi twiss file produced by ELEGANT
       argument 2 :    s location in the lattice, in meters

OUTPUT: 

@author: FBT
"""

import os
import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import sys
import subprocess



from fct_get_twiss_param_at_s_location import DisplayTwiss
from fct_get_twiss_param_at_s_location import GetTwissList

print(os.getcwd())






DisplayTwiss("VMX-twi.twi",54)

s,sIndex,betax,alphax,betay,alphay,etax,etaxp,ex0,Sdelta0 = GetTwissList("VMX-twi.twi",54)
print("***************************************************************")
print("***************************************************************")
print("closest s to the requested location    : ",s)
print("index of this value in the twiss file  : ", sIndex)
print("***************************************************************")
print("***************************************************************")
print("                                   ")
print(" Twiss parameters at that location:")
print(" ----------------------------------")
print(" betax     :", betax)
print(" alphax    :", alphax)
print(" betay     :", betay)
print(" alphay    :", alphay)
print(" etax      :",  etax)
print(" etaxp     :", etaxp)

print("                  ")
print(" Global parameters:")
print(" ------------------")
print(" the x-emittance ex0 is       :", ex0)
print(" the energy spread Sdelta0    :", Sdelta0)
print("                  ")
print("                  ")
print(' data retrieved. process completed.')