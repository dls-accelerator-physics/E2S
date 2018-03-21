# -*- coding: utf-8 -*-
#############################################################################
# reading output file
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
import sys
import os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np

# open file: FILIN = '../plotxy_scatter.dat' 
FILIN = sys.argv[1]

f = open(FILIN, 'r')
    
header = {}
#for i in range(0, 100):
flag = 1
i    = 0
X = [] # scatter X coordinate
Y = [] # scatter Y coordinate
W = [] # scatter weight  
while True:
    linea     = f.readline().strip()
    if linea == '':
        break
    elif linea[0] == '#':
        print('skip header ...')
    else:
        X.append(float(linea.split()[0]))
        Y.append(float(linea.split()[1]))
        W.append(float(linea.split()[2]))

f.close()

fig, ax =  plt.subplots(nrows=1, ncols=1)
cpf     = ax.scatter(X,Y,s=W,alpha=0.5) 


# Reversed Greys colourmap for filled contours
#fig, ax = plt.subplots(nrows=1, ncols=1)
#cpf = ax.contourf(X,Y,W, 40, cmap=cm.hot)
#plt.colorbar(cpf)

plt.show()

# https://scipython.com/book/chapter-7-matplotlib/examples/simple-surface-plots/
# fig, ax = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': '3d'})
# ax.plot_surface(X,Y,A, rstride=20, cstride=20,  cmap=cm.hot)


#if 1 == 0:
# Set the colours of the contours and labels so they're white where the
# contour fill is dark (Z < 0) and black where it's light (Z >= 0)
#    colours = ['w' if level<0 else 'k' for level in cpf.levels]
#    cp = ax.contour(X, Y, Aresc, 5, colors=colours)

#    ax.clabel(cp, fontsize=10, colors=colours)





