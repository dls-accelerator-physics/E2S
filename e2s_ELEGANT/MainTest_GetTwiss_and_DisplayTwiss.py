#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 13:05:29 2018

This function retrieves the TWISS parameters of the lattice

INPUT: argument 1 :    <root>.twi, <root>.rf file produced by ELEGANT
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
from fct_get_rf_param                      import DisplayRF
from fct_get_rf_param                      import GetRF


print(os.getcwd())


filin = sys.argv[1]
spos  = sys.argv[2]


#DisplayTwiss("VMX-twi.twi",54)
DisplayTwiss(filin+'.twi',spos) 

#s,sIndex,betax,alphax,betay,alphay,etax,etaxp,ex0,Sdelta0 = GetTwissList("VMX-twi.twi",54)
s,sIndex,betax,alphax,betay,alphay,etax,etaxp,ex0,Sdelta0 = GetTwissList(filin+'.twi',spos)
Sz0 = GetRF(filin+'.rf')
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
print(" emix       :", ex0)
print(" dE/E       :", Sdelta0)
print("                  ")
print(" sigma_z(0) :", Sz0)
print(' data retrieved. process completed.')
