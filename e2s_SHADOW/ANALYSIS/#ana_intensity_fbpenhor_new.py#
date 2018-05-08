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
FILIN = 'plotxy_scatter.dat'#sys.argv[1]



## see https://stackoverflow.com/questions/17197492/root-mean-square-error-in-python
def rmse(predictions, targets):
    differences = predictions - targets                       #the DIFFERENCEs.
    differences_squared = differences ** 2                    #the SQUAREs of ^
    mean_of_differences_squared = differences_squared.mean()  #the MEAN of ^
    rmse_val = np.sqrt(mean_of_differences_squared)           #ROOT of ^
    return rmse_val                                           #get the ^


#   Note on Uniform_distribution function below:
#   verification of consistency:
#   binwidth
#   Out[16]: 0.0008
#
#   histop[1][4]-histop[1][3]
#   Out[17]: 0.0008000000000000021

def uniform_distribution(diam,binwidth,axHistx):
    horsteps=diam/binwidth #number of horizontal steps
    verheight=len(X)/horsteps #vertical height needed to add up to N, nb of particles
    verheightmax=verheight
    kkinf=hh[1]<0-diam/2
    kksup=hh[1]>0+diam/2
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
    plt.title('in the function')
    plt.xlabel('tx')
    plt.ylabel('ty')
    blob500=500
    #to complete
    return tx,ty,histop,verheightmax,horsteps

def retrieve_hist_function(axHistx):
    histfct=axHistx.hist(x, bins=bins)
    histfctx=histfct[1]
    histfctx=histfctx[:-1] #remove the last point, for the function
    histfcty=histfct[0]
    return histfctx,histfcty
    
    

    
def two_peaked_parabolic_distribution(axHistx,delta_depth,diam,binwidth):
    #first, get the histogram info
    print('we are in the function, the binwidth is : ', binwidth)
    histfct_2p=axHistx.hist(x, bins=bins)
    histfctx_2p=histfct_2p[1]
    histfctx_2p=histfctx_2p[:-1] #remove the last point, for the function
    histfcty_2p=histfct_2p[0]
    
    #then, construct the parabolic double peak
    
    #number of horizontal steps
    horsteps=diam/binwidth 
    #vertical height needed to add up to N, nb of particles : this is the reference height of the uniform function
    verheight=len(X)/horsteps 
    verheightmax=verheight
    #construction of the high- and low-cutoff
    kkinf=hh[1]<0-(diam/2)
    kksup=hh[1]>0+(diam/2)
   # allval_2p=hh[0]>-1
    
    #Now, we set up the coefficient of the Ax^2+B parabola - only a x^2 term because of the symmetry
    area_test=np.trapz(histfcty_2p,dx=binwidth)
    print('This is just to test the area is correct...:', area_test)
    allval_2p=np.square(histfctx_2p)

    
    #BLOB
    tt=hh    
    histop=axHistx.hist(x, bins=bins)
   # print('allval_2p :', allval_2p)
   # print('tt[0] :', tt[0])
    ##############################################tt[0]=allval_2p       
    print('length of kkinf : ', len(kkinf))
    print('length of kksup : ', len(kksup))
    kkinfred=kkinf[:-1] #get rid of the extralength
    kksupred=kksup[:-1]
    
    
    allval_2p[kkinfred]=0
    allval_2p[kksupred]=0
    print(allval_2p)
    tt[0][kkinfred]=0
    tt[0][kksupred]=0
    Coeff_B=verheightmax-delta_depth
    print('verheightmax : ' , verheightmax)
    print('delta_depth : ', delta_depth)
    print('Coeff_B :', Coeff_B)
    print('area_test :', area_test)
    print(' allval_2p --------------------------------------------------------------------------')
    print(allval_2p)
    print('-------------------------------------------------------------------------------------')
    Coeff_A=3*(delta_depth)/(0.02*0.02)
    print('Coeff_A', Coeff_A)
    ty_2p=Coeff_A*allval_2p+Coeff_B  
#    tx=tt[1]
#    tx=tx[:-1]  
    
    plt.figure(9)
    plt.scatter(histfctx_2p,ty_2p)
    plt.title('in the parabolic function')
    plt.xlabel('tx')
    plt.ylabel('ty')
    blob500=500
    #to complete
    area_2p=np.trapz(ty_2p,dx=binwidth)
    print('area_2p :', area_2p)
    print('test of T : ' , Coeff_A*0.02*0.02/3 + Coeff_B)
    return tx,ty_2p,histop,verheightmax,horsteps





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


aveX = np.average(X)
stdX = np.std(X)
aveY = np.average(Y)
stdY = np.std(Y)

x=X
y=Y


plt.show()



##################################################################################################################################################


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

# no labels
#axHistx.xaxis.set_major_formatter(nullfmt)
#axHisty.yaxis.set_major_formatter(nullfmt)
#x = 30*np.random.randn(10000)


# the scatter plot:
axScatter.scatter(x, y,s=W,alpha=0.9)
#axScatter.kde(x, y)

# now determine nice limits by hand:
binwidth = 0.001
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

plt.show()



tx,ty,histop,verheightmax,horsteps=uniform_distribution(diam,binwidth,axHistx)

area_under_uniform=np.trapz(ty,dx=binwidth)
print('area under uniform :',area_under_uniform)



penaltyHorizontal=rmse(ty,histop[0])

histfctx,histfcty=retrieve_hist_function(axHistx)

area_under_hist=np.trapz(histfcty,dx=binwidth)
print('area under histogram :',area_under_hist)

#histop=axHistx.hist(x, bins=bins)
delta_depth=25
two_peaked_parabolic_distribution(axHistx,delta_depth,diam,binwidth)









