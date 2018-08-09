#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 11:17:16 2018

@author: mfc33124
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

t7InitLatwithErr-000263
Created on Tue Feb 20 11:54:58 2018

@author: mfc33124


Functions this file contains

sliceTwiss : the slicing needed for sddsfluxcurve or sddsbrightness
tuningCurveCalc
brightnessCalc : brightness
plotXRProp : plot Xr-ray properties of interest



"""

import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import sys


#============================================================================================
def concatenate_list_data(list):
    result= ''
    for element in list:
        result += str(element)
    return result
#============================================================================================

#============================================================================================
def extract_slice_from_Twiss(input_twi,output_slice_twi,IDpos,IDpos_min,IDpos_max):
    min_s=IDpos-IDpos_min
    max_s=IDpos+IDpos_max
    min_s_str=str(min_s)
    max_s_str=str(max_s)
    argCommandList=['-filter=column,s,',min_s_str,',',max_s_str]
    print('the argument command is:')
    print(' In python 3 you can write " print(*argCommandList) "')
    print(argCommandList)    
    argCommandString=concatenate_list_data(argCommandList)
    print('new:')
    print(argCommandString)    
    filter_status=subprocess.check_output(['sddsprocess',input_twi,output_slice_twi,argCommandString],stderr= subprocess.STDOUT).decode('UTF-8')    
    if filter_status=='warning: no rows selected for page 1\n': # FBT: si on oublie d'ajouter le \n ca echoue
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('%%%%%%                          SLICING FAILED !!!!!!!!      %%%%%%%%%%%%%%%%%%')
        print('')
        print(' No points exists in the specified interval, please check s, IDpos_max, or IDpos_min')
        print('-------------------------------------------------------------------------------')
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        #raise SystemExit
        sys.exit()
    else:
        return output_slice_twi
#================================================================================================


#================================================================================================
def calc_tuning_curves(input_twi,output_sddsfluxcurve, IDpos,IDpos_min,IDpos_max,modeCal,harmonics,methodCal,Ib,Cou,uPeriod,uNbPeriod,uKmin,uKmax,uPoints,e2s_ELEGANT):
    output_slice_twi='output_tuning_curve_temp.twi'
    jojolapin= extract_slice_from_Twiss(input_twi,output_slice_twi,IDpos,IDpos_min,IDpos_max)
    ######################################################################################
    # preprocessing before sddsfluxcurve
    input_slice_sddsfluxcurve=output_slice_twi
    arg1=['-mode=',modeCal]
    arg2=['-harmonics=',str(harmonics)]
    arg3=['-method=',methodCal]
    arg4=['-electronBeam=current=',str(Ib),',coupling=',str(Cou)]
    arg5=['-undulator=period=',str(uPeriod),',numberOfPeriods=',str(uNbPeriod),',kmin=',str(uKmin),',kmax=',str(uKmax),',points=',str(uPoints)]
    arg1Str=concatenate_list_data(arg1)
    arg2Str=concatenate_list_data(arg2)
    arg3Str=concatenate_list_data(arg3)
    arg4Str=concatenate_list_data(arg4)
    arg5Str=concatenate_list_data(arg5)
    status_sddsfluxcurve=subprocess.check_output([e2s_ELEGANT+'/sddsfluxcurve',input_slice_sddsfluxcurve,output_sddsfluxcurve,arg1Str,arg2Str,arg3Str,arg4Str,arg5Str],stderr= subprocess.STDOUT).decode('UTF-8')
#==================================================================================================


#==================================================================================================
def calc_brightness(input_twi,output_brightness, IDpos,IDpos_min,IDpos_max,harmonics,Kmin,Kmax,KRangeNbPoints,Ib,totalLength,periodLength,coupling,e2s_ELEGANT):
    
    output_slice_twi='output_brightness_temp.twi'
    jojolapin= extract_slice_from_Twiss(input_twi,output_slice_twi,IDpos,IDpos_min,IDpos_max)
#########################################################################################
#### for sddsbrightness
    input_slice_sddsbrightness=output_slice_twi
    output_sddsbrightness=output_brightness    
    arg1=['-harmonics=',str(harmonics)]
    arg2=['-Krange=start=',str(Kmin),',end=',str(Kmax),',points=',str(KRangeNbPoints)]
    arg3=['-current=',str(Ib)]
    arg4=['-totalLength=',str(totalLength)]
    arg5=['-periodLength=',str(periodLength)]
    arg6=['-coupling=',str(coupling)]    
    arg1Str=concatenate_list_data(arg1)
    arg2Str=concatenate_list_data(arg2)
    arg3Str=concatenate_list_data(arg3)
    arg4Str=concatenate_list_data(arg4)
    arg5Str=concatenate_list_data(arg5)
    arg6Str=concatenate_list_data(arg6)   
    
    status_sddsbrightness=subprocess.check_output([e2s_ELEGANT+'/sddsbrightness',input_slice_sddsbrightness,output_sddsbrightness,arg1Str,arg2Str,arg3Str,arg4Str,arg5Str,arg6Str],stderr= subprocess.STDOUT).decode('UTF-8')
#====================================================================================================


def plot_tune_curves(fichier):    
    cmdstring =("sddsprintout -spreadsheet=csv -col=* %s temp_sdds_tuningCurves.csv" %fichier)
    os.system(cmdstring)      
    thisdir=os.getcwd()  
    df=pd.read_csv(os.path.join(thisdir,'temp_sdds_tuningCurves.csv'))
    FluxDensity = [col for col in df if col.startswith('FluxDensity')]
    photonEnergy= [col for col in df if col.startswith('photonEnergy')]    
    ax2=df.plot(kind='scatter',y=FluxDensity[0],x=photonEnergy[0],color='Red', marker=".",alpha=.5, logy=1)    
    for bbb in range(1,len(FluxDensity)):
        print(bbb)
        df.plot(kind='scatter',y=FluxDensity[bbb],x=photonEnergy[bbb],color='DarkBlue',marker=".",ax=ax2, alpha=.5, logy=1)
    ax2.set_xlabel('E_phot (keV)')
    ax2.set_ylabel('Flux Density (phot/s/mr^2/0.1%B.W.)')
    plt.show(ax2)
    return df,FluxDensity,photonEnergy

    


def plot_brightness(fichier2): 
    cmdstring =("sddsprintout -spreadsheet=csv -col=* %s temp_sdds_brightness.csv" %fichier2)
    os.system(cmdstring)      
    thisdir=os.getcwd()  
    df=pd.read_csv(os.path.join(thisdir,'temp_sdds_brightness.csv'))
    Brightness = [col for col in df if col.startswith('Brightness')]
    photonEnergy= [col for col in df if col.startswith('photonEnergy')]   
    ax22=df.plot(kind='scatter',y=Brightness[0],x=photonEnergy[0],color='Red', marker=".",alpha=.5, logy=1)    
    for bbb in range(1,len(Brightness)):
        print(bbb)
        df.plot(kind='scatter',y=Brightness[bbb],x=photonEnergy[bbb],color='DarkBlue',marker=".",ax=ax22, alpha=.5,logy=1)
    ax22.set_xlabel('E_phot (keV)')
    ax22.set_ylabel('B (phot/s/mm^2/mrad^2/0.1%B.W.)')
    plt.show(ax22)
        
    return df,Brightness,photonEnergy



def test_tcb():
    ############################################################################################
    ############################################################################################
    ############################################################################################
    
    # EXAMPLE OF PREPARATION TO CALL THE TUNING CURVE FUNCTION
    
    ############################################################################################
    # PART 1/2 : these are the parameters for sddsprocess
    #input_twi= 'VMX.twi'
    #IDpos=282.29401
    #IDpos_min=0.00
    #IDpos_max=0.001
    input_twi='DTBA_C1a_AA.twi'
    IDpos=282.298
    IDpos_min=0.00
    IDpos_max=0.001
    
    #input_twi='DTBA_C1a_AA.twi'
    #IDpos=282.298
    #IDpos_min=0.00
    #IDpos_max=0.001
    ############################################################################################
    ### PART 2/2: now, parameters for sddsfluxcurve
    #
    output_sddsfluxcurve='dmd777.sdds'
    modeCal='density'
    #harmonics=7
    methodCal='dejus'
    Ib=0.3
    Cou=0.1
    #uPeriod=0.025
    lam_und=0.025
    #uNbPeriod=108
    Np_und=108
    uKmin=1
    uKmax=1.93095
    uPoints=10
    harm_1st  = 1
    harm_last = 15
    #############################################################################################
    # Now, we can call the function:
    calc_tuning_curves(input_twi,output_sddsfluxcurve, IDpos,IDpos_min,IDpos_max,modeCal,harm_last,methodCal,Ib,Cou,lam_und,Np_und,uKmin,uKmax,uPoints)
    
    ############################################################################################
    ############################################################################################
    ############################################################################################
    
    
    
    ############################################################################################
    #### post-processing: plotting of the tuning curves:
    
    kk,FluxDensity,photonEnergy=plot_tune_curves(output_sddsfluxcurve)
    
    
    
    ############################################################################################
    ############################################################################################
    ############################################################################################
    
    # EXAMPLE OF PREPARATION TO CALL THE CALC_BRIGHTNESS FUNCTION
    
    ############################################################################################
    # PART 1/2 : these are the parameters for sddsprocess
    input_twi= 'VMX.twi'
    IDpos=282.29401
    IDpos_min=0.00
    IDpos_max=0.001
    Cou=0.003
    
    input_twi='DTBA_C1a_AA.twi'
    IDpos=282.298
    IDpos_min=0.00
    IDpos_max=0.001
    Cou=0.1
    ############################################################################################
    # PART 2/2: setting up the inputs for sddsbrightness
    output_brightness='temp999.sdds'
    #harmonics=7
    harm_last=15
    Kmin=1
    Kmax=1.9
    KRangeNbPoints=100
    Ib=0.3
    
    #periodLength=0.025
    lam_und   = 0.025
    #coupling=0.003
    
    #totalLength=2.7
    totalLength=Np_und * lam_und
    # Now, making the call
    calc_brightness(input_twi,output_brightness, IDpos,IDpos_min,IDpos_max,harm_last,Kmin,Kmax,KRangeNbPoints,Ib,totalLength,lam_und,Cou)
    ############################################################################################
    ############################################################################################
    ############################################################################################
    
    ############################################################################################
    #### post-processing: plotting of the tuning curves:
    
    kk,Brightness,photonEnergy=plot_brightness(output_brightness)

 

 
    
   
