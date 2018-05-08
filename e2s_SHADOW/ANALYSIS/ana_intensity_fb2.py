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
from matplotlib.ticker import NullFormatter

# open file: FILIN = '../plotxy_scatter.dat' 
FILIN = '../plotxy_scatter.dat'#sys.argv[1]




#### 
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]

diam=30   # width of the uniform function


#### back to xph

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

print(len(X))
############################################################################fig, ax =  plt.subplots(nrows=1, ncols=1)
############################################################################cpf     = ax.scatter(X,Y,s=W,alpha=0.5) 



# Reversed Greys colourmap for filled contours
#fig, ax = plt.subplots(nrows=1, ncols=1)
#cpf = ax.contourf(X,Y,W, 40, cmap=cm.hot)
#plt.colorbar(cpf)

#################################################################################plt.show()

# https://scipython.com/book/chapter-7-matplotlib/examples/simple-surface-plots/
# fig, ax = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': '3d'})
# ax.plot_surface(X,Y,A, rstride=20, cstride=20,  cmap=cm.hot)


#if 1 == 0:
# Set the colours of the contours and labels so they're white where the
# contour fill is dark (Z < 0) and black where it's light (Z >= 0)
#    colours = ['w' if level<0 else 'k' for level in cpf.levels]
#    cp = ax.contour(X, Y, Aresc, 5, colors=colours)

#    ax.clabel(cp, fontsize=10, colors=colours)

# start with a rectangular Figure

Nrays  = len(X)
W      = 0.04 # FWHM in cm
t0x    = Nrays * binwidth / W

aveX = np.average(X)
stdX = np.std(X)
aveY = np.average(Y)
stdY = np.std(Y)

x=X
y=Y
nullfmt = NullFormatter()         # no labels
plt.figure(1, figsize=(8, 8))

axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)
print(axHistx)
# no labels
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

# the scatter plot:
axScatter.scatter(x, y,s=W,alpha=0.5)

# now determine nice limits by hand:
binwidth = 0.0008
xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
lim = (int(xymax/binwidth) + 1) * binwidth
print(lim)
axScatter.set_xlim((-lim, lim))
axScatter.set_ylim((-lim, lim))

bins = np.arange(-lim, lim + binwidth, binwidth)
print(len(bins))
axHistx.hist(x, bins=bins)
axHisty.hist(y, bins=bins, orientation='horizontal')
print(x)
axHistx.set_xlim(axScatter.get_xlim())
axHisty.set_ylim(axScatter.get_ylim())

plt.show()





