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


SMALL_SIZE = 30*0.4
MEDIUM_SIZE = 40*0.4
BIGGER_SIZE = 30*0.4

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

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


print(str(len(x)))
print(str(len(y)))

A = []
for iY in range(0,ypoints):
    A.append(Z[iY*xpoints:iY*xpoints+xpoints])



A=np.array(A)
Aresc = A # /np.max(A)



nullfmt = NullFormatter()         # no labels


# definitions for the axes

# default
if 1 == 0:
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02
else:
    left, width = 0.15, 0.55
    bottom, height = 0.15, 0.55
    bottom_h = left_h = left + width + 0.02
    

rect_contour = [left, bottom, width, height]
rect_histx   = [left, 1.0*bottom_h, width, 0.17]
rect_histy   = [1.00*left_h, bottom, 0.17, height]
rect_com     = [1.00*left_h,1.00*bottom_h,0.17,0.17]

# start with a rectangular Figure
fsx = 9.0 # 17.8 * .65
fsy = 9.0
plt.figure(1, figsize=(fsx,fsy))


axContour = plt.axes(rect_contour)
axX       = plt.axes(rect_histx)
axY       = plt.axes(rect_histy)
axHisty   = plt.axes(rect_histy)
axcomm    = plt.axes(rect_com)

xmin = xmin*1e6
xmax = xmax*1e6
ymin = ymin*1e6
ymax = ymax*1e6
 
axContour.set_xlim((xmin, xmax))
axContour.set_ylim((ymin, ymax))
axX.set_xlim((xmin, xmax))
axY.set_ylim((ymin, ymax))


Xresc = np.multiply(X,1e6)
Yresc = np.multiply(Y,1e6)

peakBrilliance = np.max(Aresc)
peakCoordinate = np.unravel_index(Aresc.argmax(),Aresc.shape)
xline          = peakCoordinate[1]

yline          = peakCoordinate[0]


axContour.contourf(Xresc,Yresc,Aresc, 100, cmap=cm.jet)
axContour.grid(color='w',linestyle='--')
axContour.set_xlabel('X ($\mu$m)')
axContour.set_ylabel('Y ($\mu$m)')
dy = (ymax-ymin)/len(y)
ycross = ymin+yline*dy
dx = (xmax-xmin)/len(x)
xcross = xmin+xline*dx
print(xline,ymax,ymin,dy,len(y),ycross)
print(yline,xmax,xmin,dx,len(x),xcross)

axContour.plot([xmin,xmax],[ycross,ycross],'w--')
axContour.plot([xcross,xcross],[ymin,ymax],'w--')

print(len(x))
print(len(y))

xx = np.multiply(x,1e6) # X[len(X)/2]
yy = np.multiply(y,1e6) # Y[len(Y)/2]


a_x = Aresc[yline,:] 
#a_x = Aresc.sum(axis=0) # projecting on axis-0 (x) 
#a_x = a_x/len(y)        # averaged on n. of pixels (each represents flux/mm2)

Aresc_T = np.transpose(Aresc)
axX.plot(xx,a_x)
axX.set_title('lattice: DII 6-HMBA')

a_y = Aresc[:,len(x)/2-1]
a_y = Aresc[:,xline]

#a_y = Aresc.sum(axis=1) # projecting on axis-1 (y) 
#a_y = a_y/len(x)        # averaged on n. of pixels (each represents flux/mm2)



axY.plot(a_y,yy,'C3')

binintx = len(x)/2 #495
bininty = len(y)/2 #154
int1  = range(len(x)/2-1-binintx,len(x)/2-1+binintx)
int2  = range(len(y)/2-1-bininty,len(y)/2-1+bininty)
xx    = xx[int1]
a_x   = a_x[int1]
yy    = yy[int2]
a_y   = a_y[int2]

mu1   = np.sum(xx*a_x)/np.sum(a_x)
sigma1= np.sqrt(np.sum(a_x*(xx-mu1)**2)/np.sum(a_x))
mu2   = np.sum(yy*a_y)/np.sum(a_y)
sigma2= np.sqrt(np.sum(a_y*(yy-mu2)**2)/np.sum(a_y))
#peakBrilliance = np.max(a_x)
print('mu1= '+str(mu1)+' n1= '+str(np.sum(xx*a_x))+'  d1='+str(np.sum(a_x)))
print('mu2= '+str(mu2)+' n2= '+str(np.sum(yy*a_y))+'  d2='+str(np.sum(a_y)))
print('NUM2 = '+str(np.sum(a_y*(yy-mu2)**2))+'   DEN2 = '+str(np.sum(a_y))) 

textstr = '$I=%.2e$ \n$(flux/mm^2)$ \n$\sigma_x=%.2f (\mu$m) \n$\sigma_y=%.2f (\mu$m)' % (peakBrilliance, sigma1, sigma2) # 

print('I  = '+str(peakBrilliance)+ ' (ph/s/0.1pcBW/mm2')
print('sigma_y = '+str(sigma1)+' (um)')
print('sigma_x = '+str(sigma2)+' (um)')

axX.set_xticklabels([])

    
axY.set_yticklabels([])

  
axcomm.set_frame_on(False)

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

txt_x = 0.05
txt_y = 0.95
txt_fs = 13 # 25
axcomm.text(txt_x, 0.95, textstr, transform=axcomm.transAxes, fontsize=txt_fs,
        verticalalignment='top')


plt.show()




