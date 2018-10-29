# -*- coding: utf-8 -*-
#############################################################################
# E2S - python version
# v 0.00
# MA 02/02/2018 
#
# Basic structure of the code 
#
# - select lattice (e.g. VMX, 6HMBA ...)
# - select ID (e.g. I13d = coherence branch) 
# - call elegant
# - extract Twiss parameters at ID
# - run SRW/SHADOW with these parameters as input
# - plot the results
#   1) Flux at Entry Slit
#   2) intensity image (somewhere down the optical beamline)
#
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
import os
import sys
import numpy as np

# create paths to subdirectories 
CWD = os.getcwd()
e2s_LATTICES = CWD+'/e2s_LATTICES/'  # typo LATTICE --> LATTICES corrected
e2s_SRW      = CWD+'/e2s_SRW/'
e2s_ELEGANT  = CWD+'/e2s_ELEGANT/'
e2s_BLOPTICS = CWD+'/e2s_BLOPTICS/'
e2s_SHADOW   = CWD+'/e2s_SHADOW/'


### SRWLIB      = '/dls/physics/xph53246/source_to_beamline/SRW_Dev/env/work/SRW_PROJECT/MyBeamline/'
SRWLIB      = '/dls/physics/xph53246/source_to_beamline/SRWLIB/' # MA 12/03/2018 - repository created for pure SRWlib files 


sys.path.insert(0, SRWLIB)
from srwlib import *
from uti_plot import *


sys.path.insert(0, e2s_LATTICES)    # typo LATTICE --> LATTICES corrected
sys.path.insert(0, e2s_SRW)
sys.path.insert(0, e2s_ELEGANT)
sys.path.insert(0, e2s_BLOPTICS)
sys.path.insert(0, e2s_SHADOW+'/ANALYSIS/')

 
# now we can import the functions 
from fct_get_twiss_param_at_s_location import GetTwissList
from fct_get_rf_param                  import GetRF
from fct_get_SR_param                  import GetCirc
from fct_get_SR_param                  import DisplayCirc

from fct_get_beam_param_from_twiss     import GetBeamParam

from fct_ana_intensity                 import ana_intensity
from fct_ana_intensity                 import ana_intensity_penalty

import fct_photonsPhysics_full as pp  #FBT


def read_input(filin):
    INPUT_file = filin # e.g. 'E2S.input'
    infile     = open(INPUT_file,'r')


    variables =[]; values=[];
    for line in infile:
        variable, value = line.split('=')
        variable = variable.strip()  # remove leading/traling blanks
        value    = value.strip()
        variables.append(variable)
        values.append(value)

    infile.close()

    dict={}              # create a dictionary for easy access to variables 
 
    for i in range(0,len(variables)) :
        dict[variables[i]] = values[i]

    dict['INPUT_file'] = INPUT_file
    return dict

def e2s(dict):

#
# get few input parameters for verbose summary ...
#
    
#*********** SR (choose the synchrotron radiation generator)
    SynchRad = str(dict['SynchRad'])

#*********** ID
    By_und  = float(dict['By_und'])
    lam_und = float(dict['lam_und'])
    Np_und  = float(dict['Np_und'])
    K_und   = 0.9338 * By_und * lam_und * 100
    IDpos    = float(dict['IDpos'])
    IDname   = str(dict['IDname'])

    IDpos_min = float(dict['IDpos_min'])   #FBT
    IDpos_max = float(dict['IDpos_max'])   #FBT
    Kmin      = float(dict['Kmin'])            #FBT
    Kmax      = float(dict['Kmax'])            #FBT
    KRangeNbPoints = float(dict['KRangeNbPoints'])            #FBT
    TCPoints  = float(dict['TCPoints'])            #FBT
    harm_last = float(dict['harm_last'])            #FBT
    Brightness    = float(dict['Brightness'])
    TuningCurves  = float(dict['TuningCurves'])

#*********** BEAM-shift-tilt at source 
    delta_x   = float(dict['delta_x']) 
    delta_y   = float(dict['delta_y']) 
    delta_xp  = float(dict['delta_xp']) 
    delta_yp  = float(dict['delta_yp']) 
    
#*********** MACHINE Parameters
    Ee      = float(dict['Ee'])
    Ib      = float(dict['Ib'])
# Circ    = float(dict['Circ'])
    Nbunch  = float(dict['Nbunch'])
    Cou     = float(dict['Cou'])
    LATTICE  = str(dict['LATTICE'])

#********** Calculation Parameters
    calc_type = str(dict['calc_type'])
    if SynchRad == 'SRW': 
        calc_meth = str(dict['calc_meth'])  # 0 = manual / 1 = undulator / 2 = wiggler 
        Ncores    = float(dict['Ncores'])   # only meaningful for multi-e individual cluster calculations 
    elif SynchRad == 'SHADOW':
        sour_type = str(dict['sour_type'])
        SSOUR_OE1     = str(dict['SSOUR_OE1']) # distance source-ellipsoidal mirror OE1
        SIMAG_OE1     = str(dict['SIMAG_OE1']) # distance ellipsoidal mirror - image OE1 (initially 1e+12, this is a collimating mirror)
        RMIRR_OE7     = str(dict['RMIRR_OE7']) # radius of cyclindrical mirror after 4-bump mono (specific of I20SCA) OE7
        AXMAJ     = str(dict['AXMAJ']) # axis-major of the elliptical mirror (I20SCA) OE8
        AXMIN     = str(dict['AXMIN']) # axis-minor of the elliptical mirror (I20SCA) OE8
        SSOUR_OE8     = str(dict['SSOUR_OE8']) # distance source-ellipsoidal mirror OE8 (initially 1e+12, it recieves from the collimating mirror OE1)
        SIMAG_OE8     = str(dict['SIMAG_OE8']) # distance ellipsoidal mirror - image OE8 (initially 1e+12, this is a collimating mirror)
        T_IMAGE_OE10   = str(dict['T_IMAGE_OE10']) # image distance at sample OE10
        
# ******* Input file name
    INPUT_file = dict['INPUT_file']
    
#    if SynchRad == 'SRW':
#        print('you have selected SRW ...')
#*********** BeamLine Parameters
#        slitZ   = float(dict['slitZ'])
#        slitDX  = float(dict['slitDX'])
#        slitDY  = float(dict['slitDY'])
    Ephot_ini = float(dict['Ephot_ini'])
    Ephot_end = float(dict['Ephot_end'])
#********** SRW Calculation Parameters
#        outfil    = str(dict['outfil'])
#        meshXsta  = float(dict['meshXsta'])
#        meshXfin  = float(dict['meshXfin'])
#        meshYsta  = float(dict['meshYsta'])
#        meshYfin  = float(dict['meshYfin'])
#        meshEsta  = float(dict['meshEsta'])
#        meshEfin  = float(dict['meshEfin'])
    
#    elif SynchRad == 'SHADOW':
#        print('you have selected SHADOW ...')
#********** SHADOW Beamline Parameters
#        slitZ   = float(dict['slitZ'])
#        slitDX  = float(dict['slitDX'])
#        slitDY  = float(dict['slitDY'])
#********** SHADOW Calculation Parameters
#        outfil    = str(dict['outfil'])
#    else :
#        print('no synch-rad mode selected ...')
        
#LATTICE = 'DTBA_C1a_AA'
    LATdir  = 'e2s_LATTICES/'
    SRWdir  = 'e2s_SRW/'
    SHAdir  = 'e2s_SHADOW/'
    
    # ---------------------------------
    # elegant lattice type: LATTICE.lte 
    # ---------------------------------
    
    eLTE = LATTICE+'.lte'
    
    # ----------------------------------
    # elegant steering file: LATTICE.ele 
    # ----------------------------------
    
    here = os.getcwd() # memorize the TOP directory
    
    cmd  = here+'/'+LATdir  # cd to ELEgant directory 
    os.chdir(cmd)
    eELE = LATTICE+'.ele'
    
    # ----------------------------------
    # RUN elegant
    # ----------------------------------
    
    cmd  = 'elegant '+eELE
    os.system(cmd)
    
    eTWI   = LATTICE+'.twi'
    eRF    = LATTICE+'.rf'
    eMAG   = LATTICE+'.mag'
    
    # ----------------------------------
    # retrieve results from elegant run
    # ----------------------------------
    
#spos = 282.298
    s,sIndex,betax,alphax,betay,alphay,etax,etaxp,ex0,Sdelta0 = GetTwissList(eTWI,IDpos)
    Sz0  = GetRF(eRF)
    
    Circ = GetCirc(LATTICE)
    
    #cou  = 0.01
    beam, mom = GetBeamParam([betax,alphax,betay,alphay,etax,etaxp,ex0,Sdelta0,Cou,Sz0])
    

    print("***************************************************************")
    print("closest s to the requested location    : ",s)
    print("index of this value in the twiss file  : ", sIndex)
    print("***************************************************************")
    print("                                   ")
    print(" Twiss parameters at that location:")
    print(" ----------------------------------")
    print(" betax     :", betax, " (m)")
    print(" alphax    :", alphax)
    print(" betay     :", betay, " (m)")
    print(" alphay    :", alphay)
    print(" etax      :",  etax, " (m)")
    print(" etaxp     :", etaxp)

    print(" sx     :", beam[0], " (m)")
    print(" sy     :", beam[1], " (m)")
    print(" sxp    :", beam[2], " (m)")
    print(" syp    :", beam[3], " (m)")
    

    print("                                   ")
    print(" Shifts/Tilts at source:")
    print(" ----------------------------------")
    print(" delta_x   :", delta_x, " (m)")
    print(" delta_y   :", delta_x, " (m)")
    print(" delta_xp  :", delta_xp, " (rad)")
    print(" delta_xp  :", delta_yp, " (rad)")


    
    print("                  ")
    print(" Global parameters:")
    print(" ------------------")
    print(" Eb         :", Ee, " (GeV)" )
    print(" Ib         :", Ib, " (A)")
    print(" emix       :", ex0, " (m)")
    print(" dE/E       :", Sdelta0)
    print(" Circ       :", Circ, " (m)")
    print("                  ")
    print(" sigma_z(0) :", Sz0)
    
    print("                  ")
    print(" ID:")
    print(" ------------------")
    print(" ID name       :", IDname)
    print(" Np_und        :", Np_und)
    
    
    if SynchRad == 'SRW':
        # ------------------------------
        # create SRW.input file to steer 
        # SRW calculation
        # ------------------------------  
        tgt  = here+'/'+SRWdir
        cmd  = 'cp '+here+'/'+INPUT_file+'  /'+tgt+'/SRW.input'
        os.system(cmd)
        cmd  = tgt
        os.chdir(cmd)
        os.system('echo BLname   = '+str(IDname)+' >> SRW.input\n')
        os.system('echo Circ     = '+str(Circ)+' >> SRW.input\n')
        os.system('echo sig_z    = '+str(Sz0)+' >> SRW.input\n')
        os.system('echo dE       = '+str(Sdelta0)+' >> SRW.input\n')
        os.system('echo -----Beam Twiss/Size/Moments:  =  >> SRW.input\n')
        os.system('echo emi_x    = '+str(ex0)+' >> SRW.input\n')
        os.system('echo beta_x   = '+str(betax)+' >> SRW.input\n')
        os.system('echo alpha_x  = '+str(alphax)+' >> SRW.input\n')
        os.system('echo beta_y   = '+str(betay)+' >> SRW.input\n')
        os.system('echo alpha_y  = '+str(alphay)+' >> SRW.input\n')
        os.system('echo eta_x    = '+str(etax)+' >> SRW.input\n')
        os.system('echo eta_xp   = '+str(etaxp)+' >> SRW.input\n')
        os.system('echo sig_x    = '+str(beam[0])+' >> SRW.input\n')
        os.system('echo sig_y    = '+str(beam[1])+' >> SRW.input\n')
        os.system('echo sig_xp   = '+str(beam[2])+' >> SRW.input\n')
        os.system('echo sig_yp   = '+str(beam[3])+' >> SRW.input\n')
        os.system('echo sigXX    = '+str(mom[0])+'  >> SRW.input\n')
        os.system('echo sigXXp   = '+str(mom[1])+'  >> SRW.input\n')
        os.system('echo sigXpXp  = '+str(mom[2])+'  >> SRW.input\n')
        os.system('echo sigYY    = '+str(mom[3])+'  >> SRW.input\n')
        os.system('echo sigYYp   = '+str(mom[4])+'  >> SRW.input\n')
        os.system('echo sigYpYp  = '+str(mom[5])+'  >> SRW.input\n')
        os.system('echo calc_meth   = '+str(calc_meth)+' >> SRW.input\n') # 0 =manual / 1 =undulator / 2 =wiggler
        
        cmd  = here
        os.chdir(cmd)

        # ----------------------------------
        # Run SRW 
        # ----------------------------------
    
        cmd  = here+'/'+SRWdir  # cd to ELEgant directory 
        os.chdir(cmd)
    
#os.system('python SRW_I13d_individual_electrons.py SRW.input')
###### os.system('./submit_runbatch_Individual.sh')
        
        print("Calc Type is "+calc_type)
        if calc_type == 'multie_from_individual':
            print("SRW - INTENSITY CALCULATION - front calculation summing up many electrons (partially coherent)")
            os.system('/dls_sw/apps/python/anaconda/1.7.0/64/bin/python SRW_individualelectrons.py SRW.input')
            
        elif calc_type == 'multie':
            print("SRW - INTENSITY CALCULATION - multi-e mode (fully coherent)") 
        #os.system(' /dls_sw/apps/python/anaconda/1.7.0/64/bin/python SRW_I13d_intensity.py SRW.input')
            os.system(' /dls_sw/apps/python/anaconda/1.7.0/64/bin/python SRW_intensity.py SRW.input')
            
        elif calc_type == 'flux':
            print("SRW - FLUX CALCULATION - interactive mode ... ")
            os.system(' /dls_sw/apps/python/anaconda/1.7.0/64/bin/python SRW_flux.py SRW.input')
            
        elif calc_type == 'flux_cluster':
            print("SRW - FLUX CALCULATION - using the cluster ... ")
    #os.system(' /dls_sw/apps/python/anaconda/1.7.0/64/bin/python SRW_I13d_flux.py SRW.input')
        #### os.system(' /dls_sw/apps/python/anaconda/1.7.0/64/bin/python SRW_flux.py SRW.input')
            os.system('./submit_runbatch_Flux.sh')
            
        elif calc_type == 'multie_from_individual_cluster':
            print("SRW - INTENSITY CALCULATION - front calculation summing up many electrons - qsub on the cluster (partially coherent)") 
            print("creating the submit_runbatch_Individual.sh script ...", str(int(Ncores)))
            
            os.system('echo qsub -q ap-high.q -l redhat_release=rhel6 -V -pe openmpi '+str(int(Ncores))+' runbatch_Individual.sh  > qsub_Individual.sh\n')
#            os.system('./submit_runbatch_Individual.sh')
            os.system('source ./qsub_Individual.sh')
        
        else:
            print("No calc_type found ... ")
            print("Try again using on the following cases: ")
            print("1. flux ")
            print("2. multie     (coherent) ")
            print("3. multie_from_individual         (partially coherent) ")
            print("4. multie_from_individual_cluster (partially coherent) ")
            
            
            

        cmd = here
        os.chdir(cmd)

        print('Brightness is: ', Brightness)  
        print('TuningCurves is: ', TuningCurves)  
        print('Kmin is : ', Kmin)
        print('Kmax is : ',Kmax)
        print('KRangeNbPoints is : ', KRangeNbPoints)

        if Brightness ==1:
            
            # preparatory input for Brightness calculation 
            # input_twi='DTBA_C1a_AA.twi'
            input_twi         = e2s_LATTICES+'/'+LATTICE+'.twi'  #'DTBA_C1a_AA.twi'
            # output_brightness = 'temp999.sdds'
            output_brightness = 'out_bright.sdds'
            totalLength       = Np_und * lam_und

            # Now, making the call
            pp.calc_brightness(input_twi,output_brightness, IDpos,IDpos_min,IDpos_max,harm_last,Kmin,Kmax,KRangeNbPoints,Ib,totalLength,lam_und,Cou,e2s_ELEGANT)
            kk,Brightness,photonEnergy=pp.plot_brightness(output_brightness)

        if TuningCurves ==1:
                    
            # preparatory input for Tuning Curves calculation
            # input_twi='DTBA_C1a_AA.twi'
            input_twi            = e2s_LATTICES+'/'+LATTICE+'.twi'  #'DTBA_C1a_AA.twi'

            # output_sddsfluxcurve = 'dmd777.sdds'
            output_sddsfluxcurve = 'out_tc.sdds'
   
            modeCal              = 'density'
            methodCal            = 'dejus'

            # Now, we can call the function:
            pp.calc_tuning_curves(input_twi,output_sddsfluxcurve, IDpos,IDpos_min,IDpos_max,modeCal,harm_last,methodCal,Ib,Cou,lam_und,Np_und,Kmin,Kmax,TCPoints,e2s_ELEGANT)
            
            #### post-processing: plotting of the tuning curves:            
            kk,FluxDensity,photonEnergy=pp.plot_tune_curves(output_sddsfluxcurve)
            
    elif SynchRad == 'SHADOW':
        Nrays = float(dict['Nrays'])

        # ------------------------------
        # create SHA.input file to steer 
        # SHAdow calculation
        # ------------------------------  
        tgt  = here+'/'+SHAdir
        cmd  = 'cp '+here+'/'+INPUT_file+'  /'+tgt+'/SHA.input'
        os.system(cmd)
        cmd  = tgt
        os.chdir(cmd)
        #
        # a) create the source input file    --> correspond to wiggler_source.TEMPLATE    
        #
        
        if sour_type == 'wiggler':
            print('')
            os.system('echo epath > SHA_source.input\n')
            os.system('echo 1 >> SHA_source.input\n')
            os.system('echo '+str(int(np.ceil(Np_und)))+' >> SHA_source.input\n')
            os.system('echo 1 >> SHA_source.input\n')
            os.system('echo '+str(lam_und)+' >> SHA_source.input\n')
            os.system('echo '+str(K_und)+' >> SHA_source.input\n')
            os.system('echo '+str(Ee)+' >> SHA_source.input\n')
            os.system('echo 501 >> SHA_source.input\n')
            os.system('echo 1.0 >> SHA_source.input\n')
            os.system('echo xshwig.par >> SHA_source.input\n')
            os.system('echo xshwig.traj >> SHA_source.input\n')
            os.system('echo nphoton >> SHA_source.input\n')
            os.system('echo xshwig.traj >> SHA_source.input\n')
            os.system('echo 1 >> SHA_source.input\n')
            os.system('echo '+str(Ephot_ini)+' >> SHA_source.input\n')
            os.system('echo '+str(Ephot_end)+' >> SHA_source.input\n')
            os.system('echo xshwig.sha >> SHA_source.input\n')
            os.system('echo input_source >> SHA_source.input\n')
            os.system('echo 1 >> SHA_source.input\n')
            os.system('echo 0 >> SHA_source.input\n')
            os.system('echo '+str(int(Nrays))+' >> SHA_source.input\n')
            os.system('echo 3398755 >> SHA_source.input\n')
            os.system('echo 2 >> SHA_source.input\n')
            os.system('echo xsh_slit_tmp.dat >> SHA_source.input\n')
            os.system('echo 0 >> SHA_source.input\n')
            os.system('echo 1 >> SHA_source.input\n')
            os.system('echo xshwig.sha >> SHA_source.input\n')
            os.system('echo 100 >> SHA_source.input\n')
            os.system('echo 1 >> SHA_source.input\n')
            os.system('echo 1 >> SHA_source.input\n')
            os.system('echo 1 >> SHA_source.input\n')
            os.system('echo 1 >> SHA_source.input\n')
            os.system('echo '+str(beam[0]*100)+' >> SHA_source.input\n')
            os.system('echo '+str(beam[1]*100)+' >> SHA_source.input\n')
            os.system('echo '+str(ex0*100)+' >> SHA_source.input\n')
            os.system('echo 0.0 >> SHA_source.input\n')
            os.system('echo '+str(ex0*Cou*100)+' >> SHA_source.input\n')
            os.system('echo 0.0 >> SHA_source.input\n')
            os.system('echo 3 >> SHA_source.input\n')
            os.system('echo 1 >> SHA_source.input\n')
            os.system('echo 1 >> SHA_source.input\n')
            os.system('echo source >> SHA_source.input\n')
            os.system('echo systemfile >> SHA_source.input\n')
            os.system('echo exit >> SHA_source.input\n')
             
        elif sour_type =='undulator':
            print('')

        #
        # b) create the SHADOW trace file        
        #
            
        os.system('echo trace > SHA_trace.input\n')
        os.system('echo systemfile >> SHA_trace.input\n')
        os.system('echo 0 >> SHA_trace.input\n')
        os.system('echo exit >> SHA_trace.input\n')

        #
        # c) create the plotxy files (for further analysis)        
        #

        os.system('echo plotxy > SHA_oe10.input\n')
        os.system('echo star.10 >> SHA_oe10.input\n')
        os.system('echo 2 >> SHA_oe10.input\n')
        os.system('echo ''viva'' >> SHA_oe10.input\n')
        os.system('echo 1 >> SHA_oe10.input\n')
        os.system('echo 3 >> SHA_oe10.input\n')
        os.system('echo 0 >> SHA_oe10.input\n')
        os.system('echo 0 >> SHA_oe10.input\n')
        os.system('echo 0 >> SHA_oe10.input\n')
        os.system('echo 50 >> SHA_oe10.input\n')
        os.system('echo 50 >> SHA_oe10.input\n')
        os.system('echo exit >> SHA_oe10.input\n')

        os.system('echo beg_plotxy > SHA_beg.input\n')
        os.system('echo begin.dat >> SHA_beg.input\n')
        os.system('echo 2 >> SHA_beg.input\n')
        os.system('echo ''viva'' >> SHA_beg.input\n')
        os.system('echo 1 >> SHA_beg.input\n')
        os.system('echo 3 >> SHA_beg.input\n')
        os.system('echo 0 >> SHA_beg.input\n')
        os.system('echo 0 >> SHA_beg.input\n')
        os.system('echo 0 >> SHA_beg.input\n')
        os.system('echo 50 >> SHA_beg.input\n')
        os.system('echo 50 >> SHA_beg.input\n')
        os.system('echo exit >> SHA_beg.input\n')

        #
        # d) modify beamline parameters (action on OE's)
        #    for now only the radius of M3 (cyl-mirror)
        #    RMIRR
        
        with open('start.01','r') as input_file, open('_start.01','w') as output_file:
            for line in input_file:
                L = line.split()[0]
                if L == 'SSOUR':
                    output_file.write(line.split()[0]+' '+line.split()[1]+'     '+SSOUR_OE1+'\n')
                elif L == 'SIMAG':
                    output_file.write(line.split()[0]+' '+line.split()[1]+'     '+SIMAG_OE1+'\n')    
                else:
                    output_file.write(line)
        os.system('cp _start.01 start.01')
            
        with open('start.07','r') as input_file, open('_start.07','w') as output_file:
            for line in input_file:
                L = line.split()[0]
                if L == 'RMIRR':
                    output_file.write(line.split()[0]+' '+line.split()[1]+'     '+RMIRR_OE7+'\n')
                else:
                    output_file.write(line)
        os.system('cp _start.07 start.07')

#        with open('start.08','r') as input_file, open('_start.08','w') as output_file:
#            for line in input_file:
#                L = line.split()[0]
#                if L == 'AXMAJ':
#                    output_file.write(line.split()[0]+' '+line.split()[1]+'     '+AXMAJ+'\n')
#                elif L == 'AXMIN':
#                    output_file.write(line.split()[0]+' '+line.split()[1]+'     '+AXMIN+'\n')    
#                else:
#                    output_file.write(line)
#        os.system('cp _start.08 start.08')

        with open('start.08','r') as input_file, open('_start.08','w') as output_file:
            for line in input_file:
                L = line.split()[0]
                if L == 'SSOUR':
                    output_file.write(line.split()[0]+' '+line.split()[1]+'     '+SSOUR_OE8+'\n')
                elif L == 'SIMAG':
                    output_file.write(line.split()[0]+' '+line.split()[1]+'     '+SIMAG_OE8+'\n')    
                else:
                    output_file.write(line)
        os.system('cp _start.08 start.08')
                    
        with open('start.10','r') as input_file, open('_start.10','w') as output_file:
            for line in input_file:
                L = line.split()[0]
                if L == 'T_IMAGE':
                    output_file.write(line.split()[0]+' '+line.split()[1]+'     '+T_IMAGE_OE10+'\n')
                else:
                    output_file.write(line)
        os.system('cp _start.10 start.10')

        

        os.system('echo BLname   = '+str(IDname)+' >> SHA.input\n')
        os.system('echo Circ     = '+str(Circ)+' >> SHA.input\n')
        os.system('echo sig_z    = '+str(Sz0)+' >> SHA.input\n')
        os.system('echo dE       = '+str(Sdelta0)+' >> SHA.input\n')
        os.system('echo -----Beam Twiss/Size/Moments:  =  >> SHA.input\n')
        os.system('echo emi_x    = '+str(ex0)+' >> SHA.input\n')
        os.system('echo beta_x   = '+str(betax)+' >> SHA.input\n')
        os.system('echo alpha_x  = '+str(alphax)+' >> SHA.input\n')
        os.system('echo beta_y   = '+str(betay)+' >> SHA.input\n')
        os.system('echo alpha_y  = '+str(alphay)+' >> SHA.input\n')
        os.system('echo eta_x    = '+str(etax)+' >> SHA.input\n')
        os.system('echo eta_xp   = '+str(etaxp)+' >> SHA.input\n')
        os.system('echo sig_x    = '+str(beam[0])+' >> SHA.input\n')
        os.system('echo sig_y    = '+str(beam[1])+' >> SHA.input\n')
        os.system('echo sig_xp   = '+str(beam[2])+' >> SHA.input\n')
        os.system('echo sig_yp   = '+str(beam[3])+' >> SHA.input\n')
        os.system('echo sigXX    = '+str(mom[0])+'  >> SHA.input\n')
        os.system('echo sigXXp   = '+str(mom[1])+'  >> SHA.input\n')
        os.system('echo sigXpXp  = '+str(mom[2])+'  >> SHA.input\n')
        os.system('echo sigYY    = '+str(mom[3])+'  >> SHA.input\n')
        os.system('echo sigYYp   = '+str(mom[4])+'  >> SHA.input\n')
        os.system('echo sigYpYp  = '+str(mom[5])+'  >> SHA.input\n')
        #os.system('echo calc_meth   = '+str(calc_meth)+' >> SHA.input\n')
        
        cmd  = here
        os.chdir(cmd)
        
        # ----------------------------------
        # Run SHADOW
        # ----------------------------------
    
        cmd  = here+'/'+SHAdir  # cd to ELEgant directory 
        os.chdir(cmd)
  
        print("Calc Type is "+calc_type)
 
        if calc_type == 'multie':
            print("SHADOW - INTENSITY CALCULATION - multi-e mode") 

            os.system('/dls_sw/apps/xop/2.4//extensions/shadowvui/shadow3/shadow3 < SHA_source.input')
            os.system('/dls_sw/apps/xop/2.4//extensions/shadowvui/shadow3/shadow3 < SHA_trace.input')#I20_SCA_branch.inp
            
            os.system('/dls_sw/apps/xop/2.4//extensions/shadowvui/shadow3/shadow3 < SHA_oe10.input')
            #os.system('/dls_sw/apps/xop/2.4//extensions/shadowvui/shadow3/shadow3 < SHA_beg.input')
            
            #### aveX, aveY, sigmaX, sigmaY = ana_intensity('plotxy_scatter.dat')
            aveX, aveY, sigmaX, sigmaY, penaltyHorizontal = ana_intensity_penalty('plotxy_scatter.dat')
            #### aveX, aveY, sigmaX, sigmaY, penaltyHorizontal = ana_intensity_penalty_parab('plotxy_scatter.dat')

        cmd = here
        os.chdir(cmd)
        


def main():
    filin = sys.argv[1]
    print(filin)
    dict = read_input(filin)
    print(dict)
    e2s( dict )


if __name__ == '__main__':
    main()

    
