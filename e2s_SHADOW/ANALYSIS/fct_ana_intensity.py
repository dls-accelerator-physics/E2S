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

def ana_intensity(FILIN):
# open file: FILIN = '../plotxy_scatter.dat' 
# FILIN = sys.argv[1]

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
            
    sigmaX = np.std(X)
    aveX   = np.mean(X) 
    sigmaY = np.std(Y)
    aveY   = np.mean(Y) 
    
            
    print('X = '+str(aveX*10000)+' +/- '+str(sigmaX*10000)+' (um)')
    print('Y = '+str(aveY*10000)+' +/- '+str(sigmaY*10000)+' (um)')
    
    f.close()
            
    fig, ax =  plt.subplots(nrows=1, ncols=1)
    ax.grid()
    cpf     = ax.scatter(X,Y,s=W,alpha=0.5) 
    ax.set_xlim([-0.11, 0.11])
    ax.set_ylim([-0.040, 0.015])
    plt.suptitle('RMIRR = 8.75 cm')
    
    os.system('echo '+str(aveX)+' '+str(aveY)+' '+str(sigmaX)+' '+str(sigmaY)+'  > objectives')

    return aveX, aveY, sigmaX, sigmaY     


def ana_intensity_penalty(FILIN):
    ## see https://stackoverflow.com/questions/17197492/root-mean-square-error-in-python
    def rmse(predictions, targets, sigma):
        differences = predictions - targets           #the DIFFERENCEs.
        differences_squared = differences ** 2        #the SQUAREs of ^
        #mean_of_differences_squared = differences_squared.mean()  #the MEAN of ^
        #rmse_val = np.sqrt(mean_of_differences_squared)           #ROOT of ^
        rmse_val = np.sqrt(np.sum(differences_squared)) # * sigma
        return rmse_val                                           #get the ^
    
    
    
    
#### 
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02
    
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    
    diam=30   # width of the uniform function

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
   
    sigmaX = np.std(X)
    aveX   = np.mean(X) 
    sigmaY = np.std(Y)
    aveY   = np.mean(Y) 
    
            
    print('X = '+str(aveX*10000)+' +/- '+str(sigmaX*10000)+' (um)')
    print('Y = '+str(aveY*10000)+' +/- '+str(sigmaY*10000)+' (um)')        
    
    f.close()
    
    aveX = np.average(X)
    stdX = np.std(X)
    aveY = np.average(Y)
    stdY = np.std(Y)
    
    x=X
    y=Y


    diam=0.04  # FWHM in cm

    nullfmt = NullFormatter()         # no labels
    
# definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02
    
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, 1.02*bottom_h, width, 0.2]
    rect_histy = [1.06*left_h, bottom, 0.2, height]
    rect_com=[1.06*left_h,1.02*bottom_h,0.2,0.2]
    
    
# start with a rectangular Figure
    plt.figure(1, figsize=(9, 9))
    
    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    axcomm=plt.axes(rect_com)

    # the scatter plot:
    axScatter.scatter(x, y,s=W,alpha=0.9)
#axScatter.kde(x, y)
    
# now determine nice limits by hand:
    binwidth = 0.0004#  0.0008
    xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
    lim = (int(xymax/binwidth) + 1) * binwidth
    
    axScatter.set_xlim((-lim, lim))
    axScatter.set_ylim((-lim, lim))
    
    bins = np.arange(-lim, lim + binwidth, binwidth)
    hh=axHistx.hist(x, bins=bins)
    hr=hh[0]
    vv=axHisty.hist(y, bins=bins, orientation='horizontal')
    
    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())
    
    mu1 = np.mean(hh[1])
    ave1= np.average(hh[0])
    sigma1 = np.std(hh[1])
    mu2 = np.mean(vv[1])
    ave2= np.mean(vv[0])
    sigma2 = np.std(vv[1])
#textstr = '$\mu_1=%.2f$\n$\mathrm{median}=%.2f$\n$\sigma_1=%.2f$\n\n$\mu_2=%.2f$\n$\mathrm{median2}=%.2f$\n$\sigma_2=%.2f$' % (mu1, median1, sigma1,mu2,median2,sigma2)
    textstr = '$\mu_1=%.2f$\n$\mathrm{aver. 1}=%.2f$\n$\sigma_1=%.2f$\n\n$\mu_2=%.2f$\n$\mathrm{aver. 2}=%.2f$\n$\sigma_2=%.2f$' % (0, 0, 0,0,0,0)
    
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
    
    axcomm.text(0.05, 0.95, textstr, transform=axcomm.transAxes, fontsize=14,
                verticalalignment='top')
    
    # plt.show()
    
    horsteps=diam/binwidth #number of horizontal steps
    verheight=len(X)/horsteps #vertical height needed to add up to N, nb of particles
    
    kkinf=hh[1]<mu1-diam/2
    kksup=hh[1]>mu1+diam/2
    allval=hh[0]>-1
    tt=hh
    
    
    
    histop=axHistx.hist(x, bins=bins)
    tt[0][allval]=verheight
    
    
    kkinfred=kkinf[:-1]
    kksupred=kksup[:-1]
    tt[0][kkinfred]=0
    tt[0][kksupred]=0
    ty=tt[0]
    
    tx=tt[1]
    tx=tx[:-1]
    
    plt.figure(2)
    plt.scatter(tx,ty)
    
    
    
    penaltyHorizontal=rmse(ty,histop[0],sigmaX)
    os.system('echo '+str(aveX)+' '+str(aveY)+' '+str(sigmaX)+
              ' '+str(sigmaY)+' '+str(penaltyHorizontal)+'  > objectives')
    
    print('Uniformity_400um = '+str(penaltyHorizontal))

    return aveX, aveY, sigmaX, sigmaY, penaltyHorizontal

 




