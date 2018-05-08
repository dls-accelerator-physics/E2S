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

# Open file
#FILIN = 'ex08_res_int_part_coh1_X-0.001_Y0.0002.dat'
#FILIN='ex08_res_int_part_coh1_X0.001_Y-0.0002.dat'

#FILIN='ex08_res_int_part_coh1_X0.0_Y0.0_sX50um_sX50um_5000e.dat'
#FILIN='ex08_res_int_part_coh1_X0.0_Y0.0_sX400um_sX400um_5000e.dat'
#FILIN='ex08_res_int_part_coh1_X0.0_Y0.0_sX1um_sX1um_5000e.dat'
# FILIN='ex08_res_int_part_coh1_X0.0_Y0.0_sX1um_sX1um_50e.dat'
# FILIN='ex08_res_int_part_coh1_X0.0_Y0.0_sX200um_sX200um_500e.dat'


# FILIN = 'ex08_res_int_part_coh1.dat'
# FILIN = 'ex08_res_int3.dat' 
  
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
        if 'Initial Horizontal Position' in header[i]:
            lin     = header[i].strip()
            col     = lin.split()
            Xmin = col[0]
            Xmin = Xmin[1:]
            xmin = float(Xmin)
        if 'Final Horizontal Position' in header[i]:
            lin     = header[i].strip()
            col     = lin.split()
            Xmax = col[0]
            Xmax = Xmax[1:]
            xmax = float(Xmax)
        if 'Number of points vs Horizontal Position' in header[i]:
            lin     = header[i].strip()
            col     = lin.split()
            Xpoints = col[0]
            Xpoints = Xpoints[1:]
            xpoints = int(Xpoints)
        if 'Initial Vertical Position' in header[i]:
            lin     = header[i].strip()
            col     = lin.split()
            Ymin = col[0]
            Ymin = Ymin[1:]
            ymin = float(Ymin)
        if 'Final Vertical Position' in header[i]:
            lin     = header[i].strip()
            col     = lin.split()
            Ymax = col[0]
            Ymax = Ymax[1:]
            ymax = float(Ymax)
        if 'Number of points vs Vertical Position' in header[i]:
            lin     = header[i].strip()
            col     = lin.split()
            Ypoints = col[0]
            Ypoints = Ypoints[1:]
            ypoints = int(Ypoints)
            
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
#    print(data[cnt])
    cnt=cnt+1

f.close()

Z = data
# Z = data[0:10000]

y = []
for iY in range(0,ypoints):
    y.append(ymin + (ymax-ymin)/ypoints*iY + (ymax-ymin)/ypoints/2)

x = [] 
for iX in range(0,xpoints):
    x.append(xmin + (xmax-xmin)/xpoints*iX + (xmax-xmin)/xpoints/2)

X, Y = np.meshgrid(x,y)

A = []
for iY in range(0,ypoints):
    A.append(Z[iY*xpoints:iY*xpoints+xpoints])



A=np.array(A)
Aresc = A # /np.max(A)



nullfmt = NullFormatter()         # no labels


# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.02

rect_contour = [left, bottom, width, height]
rect_histx   = [left, 1.02*bottom_h, width, 0.15]
rect_histy   = [1.06*left_h, bottom, 0.15, height]
rect_com     = [1.06*left_h,1.02*bottom_h,0.15,0.15]

# start with a rectangular Figure
plt.figure(1, figsize=(9,9))


axContour = plt.axes(rect_contour)
axX       = plt.axes(rect_histx)
axY       = plt.axes(rect_histy)
axHisty   = plt.axes(rect_histy)
axcomm    = plt.axes(rect_com)
axContour.set_xlim((-0.0008e6, 0.0008e6))
axContour.set_ylim((-0.0003e6, 0.0003e6))
axX.set_xlim((-0.0008e6, 0.0008e6))
axY.set_ylim((-0.0003e6, 0.0003e6))

Xresc = np.multiply(X,1e6)
Yresc = np.multiply(Y,1e6)
axContour.contourf(Xresc,Yresc,Aresc, 40, cmap=cm.hot)
xx = np.multiply(x,1e6) # X[len(X)/2]
yy = np.multiply(y,1e6) # Y[len(Y)/2]
a_x = Aresc[:][126]
Aresc_T = np.transpose(Aresc)
#a_x = Aresc[126][:]
axX.plot(xx,a_x)
a_y = Aresc_T[:][2970/2]
axY.plot(a_y,yy)

mu1   = np.sum(xx*a_x)/np.sum(a_x)
sigma1= np.sqrt(np.sum(a_x*(xx-mu1)**2)/np.sum(a_x))
mu2   = np.sum(yy*a_y)/np.sum(a_y)
sigma2= np.sqrt(np.sum(a_y*(yy-mu2)**2)/np.sum(a_y))
peakBrilliance = np.max(a_x)

textstr = '$B=%.2e$ \n$\mu_1=%.2f (um)$ \n$\sigma_1=%.2f (um)$ \n$\mu_2=%.2f (um)$ \n$\sigma_2=%.2f (um)$' % (peakBrilliance, mu1, sigma1, mu2, sigma2) # multiply x 1000 pour avo

axcomm.xaxis.set_major_formatter(nullfmt)
axcomm.yaxis.set_major_formatter(nullfmt)
axcomm.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off
        
axcomm.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    left='off',
    labelbottom='off') # labels along the bottom edge are off

axcomm.text(0.05, 0.95, textstr, transform=axcomm.transAxes, fontsize=9,
        verticalalignment='top')


plt.show()

# https://scipython.com/book/chapter-7-matplotlib/examples/a-simple-contour-plot/
# Reversed Greys colourmap for filled contours
fig, ax = plt.subplots(nrows=1, ncols=1)
cpf = ax.contourf(X,Y,Aresc, 40, cmap=cm.hot)
plt.colorbar(cpf)

# https://scipython.com/book/chapter-7-matplotlib/examples/simple-surface-plots/
# fig, ax = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': '3d'})
# ax.plot_surface(X,Y,A, rstride=20, cstride=20,  cmap=cm.hot)


if 1 == 0:
# Set the colours of the contours and labels so they're white where the
# contour fill is dark (Z < 0) and black where it's light (Z >= 0)
    colours = ['w' if level<0 else 'k' for level in cpf.levels]
    cp = ax.contour(X, Y, Aresc, 5, colors=colours)

    ax.clabel(cp, fontsize=10, colors=colours)



plt.show()


