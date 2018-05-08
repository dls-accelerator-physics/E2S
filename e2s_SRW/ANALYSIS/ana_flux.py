# -*- coding: utf-8 -*-
#############################################################################
# reading output file and plot flux curve ... MA 11/1/2018
#
#  python ana_flux.py  ex08_res_int1.dat
#
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
import sys
import os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np

# Open file
# FILIN = 'ex08_res_int1.dat' 
   
FILIN = sys.argv[1]

f = open(FILIN, 'r')

header = {}
#for i in range(0, 100):
flag = 1
i    = 0
while flag > 0:
    linea     = f.readline()
    if linea[0] is '#':
        header[i] = linea
        if 'Initial Photon Energy' in header[i]:
            lin     = header[i].strip()
            col     = lin.split()
            Emin = col[0]
            Emin = Emin[1:]
            emin = float(Emin)
        if 'Final Photon Energy' in header[i]:
            lin     = header[i].strip()
            col     = lin.split()
            Emax = col[0]
            Emax = Emax[1:]
            emax = float(Emax)
        if 'Number of points vs Photon Energy' in header[i]:
            lin     = header[i].strip()
            col     = lin.split()
            Epoints = col[0]
            Epoints = Epoints[1:]
            epoints = int(Epoints)
            
        i = i + 1
        
    else:
        flag = -1 
  


# Loop over lines and extract variables of interest
cnt = 0
data = []
columns = linea.split()
data.append(float(columns[0]))  # the first data.append is on the linea 

for line in f:
    line = line.strip()
    columns = line.split()
    
    
    data.append(float(columns[0]))
    #print(data[cnt])
    cnt=cnt+1

f.close()

Z = data

F=np.array(Z)

E = []
for iE in range(0,epoints):
    E.append(emin + (emax-emin)/epoints*(iE-1))


fig, ax = plt.subplots(nrows=1, ncols=1)
#cpf = ax.plot(E,np.log(F))
cpf = ax.plot(E,F)


plt.show()


