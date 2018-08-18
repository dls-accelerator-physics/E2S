#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 10:13:26 2018

@author: mfc33124
"""


import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import subprocess

import subprocess
import pandas as pd
import os
import csv
from shlex import split



os.system('cd ')


path=os.getcwd()



from subprocess import Popen, PIPE


def fct_getCrystalParam(energy):
    #from shlex import split
    #p1 = Popen("./fb_getX0h2.sh", stdout=PIPE)
    energyStr=str(energy)
    #cmdString=("./fb_getX0h3.sh %s" %energyStr)
    #os.system(cmdString)
    tl_temp=subprocess.check_output(['./fb_getX0h3.sh',energyStr],stderr= subprocess.STDOUT).decode('UTF-8')
    crystalParams = [float(n) for n in tl_temp.split()]
    return crystalParams

    #tempecho=subprocess.check_output(['./fb_getX0h2.sh'],stderr= subprocess.STDOUT).decode('UTF-8')
    #x=[float(i) for i in tempecho.split()]    
    #return x


# =============================================================================
# energy = 11.208
# crystalParams=fct_getCrystalParam(energy)
# 
# 
# # =============================================================================
# # psi0rSi111 = crystalParams[1] # -7.757245827e-6; 
# # psi0iSi111 = crystalParams[2] #9.506848329e-8 #Real and imaginary parts of 0-th Fourier component of crystal polarizability
# # psihrSi111 = crystalParams[3] #-4.095903022e-6; 
# # psihiSi111 = crystalParams[4]  #6.637430983e-8 #Real and imaginary parts of h-th Fourier component of crystal polarizability
# # psihbrSi111 = psihrSi111; 
# # psihbiSi111 = psihiSi111 #Real and imaginary parts of -h-th Fourier component of crystal polarizability
# # =============================================================================
# 
# 
# print(crystalParams)
# print('finished')
# =============================================================================

