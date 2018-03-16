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
# - run SRW with these parameters as input
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
e2s_LATTICE = CWD+'/e2s_LATTICE/'
e2s_SRW     = CWD+'/e2s_SRW/'
e2s_ELEGANT = CWD+'/e2s_ELEGANT/'

SRWLIB      = '/dls/physics/xph53246/source_to_beamline/SRW_Dev/env/work/SRW_PROJECT/MyBeamline/'
###SRWLIB      = '/dls/physics/xph53246/source_to_beamline/SRWLIB/' # MA 12/03/2018 - repository created for pure SRWlib files 


sys.path.insert(0, SRWLIB)
from srwlib import *
from uti_plot import *


sys.path.insert(0, e2s_LATTICE)
sys.path.insert(0, e2s_SRW)
sys.path.insert(0, e2s_ELEGANT)

 
# now we can import the fucntions 
from fct_get_twiss_param_at_s_location import GetTwissList
from fct_get_rf_param                  import GetRF
from fct_get_SR_param                  import GetCirc
from fct_get_SR_param                  import DisplayCirc

from fct_get_beam_param_from_twiss     import GetBeamParam


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
    
#*********** ID
    By_und  = float(dict['By_und'])
    lam_und = float(dict['lam_und'])
    Np_und  = float(dict['Np_und'])
    K_und   = 0.9338 * By_und * lam_und * 100
    IDpos    = float(dict['IDpos'])
#*********** SR Parameters
    Ee      = float(dict['Ee'])
    Ib      = float(dict['Ib'])
# Circ    = float(dict['Circ'])
    Nbunch  = float(dict['Nbunch'])
    Cou     = float(dict['Cou'])
    
#*********** BeamLine Parameters
    slitZ   = float(dict['slitZ'])
    slitDX  = float(dict['slitDX'])
    slitDY  = float(dict['slitDY'])
    Ephot_ini = float(dict['Ephot_ini'])
    Ephot_end = float(dict['Ephot_end'])
    
#********** Calculation Parameters (SRW)
    outfil    = str(dict['outfil'])
    meshXsta  = float(dict['meshXsta'])
    meshXfin  = float(dict['meshXfin'])
    meshYsta  = float(dict['meshYsta'])
    meshYfin  = float(dict['meshYfin'])
    meshEsta  = float(dict['meshEsta'])
    meshEfin  = float(dict['meshEfin'])
    
    IDname   = str(dict['IDname'])
    LATTICE  = str(dict['LATTICE'])
    calc_type = str(dict['calc_type'])
    calc_meth = str(dict['calc_meth'])  # 0 = manual / 1 = undulator / 2 = wiggler 

# ******* Input file name
    INPUT_file = dict['INPUT_file']


#LATTICE = 'DTBA_C1a_AA'
    LATdir  = 'e2s_LATTICES/'
    SRWdir  = 'e2s_SRW/'
    
    
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
    # RUN elelgant
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
    
    print("                  ")
    print(" Global parameters:")
    print(" ------------------")
    print(" emix       :", ex0)
    print(" dE/E       :", Sdelta0)
    print(" Circ       :", Circ)
    print("                  ")
    print(" sigma_z(0) :", Sz0)
    
    print("                  ")
    print(" ID:")
    print(" ------------------")
    print(" ID name       :", IDname)
    print(" Np_und        :", Np_und)
    
    
    
    
    tgt  = here+'/'+SRWdir
    cmd  = 'cp '+here+'/'+INPUT_file+'  /'+tgt+'/SRW.input'
    os.system(cmd)
    cmd  = tgt
    os.chdir(cmd)
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
    if calc_type == 'individual':
        print("INTENSITY CALCULATION - individual front calculation")
        os.system('./submit_runbatch_Individual.sh')
    
    elif calc_type == 'multie':
        print("INTENSITY CALCULATION - multi-e mode") 
        #os.system(' /dls_sw/apps/python/anaconda/1.7.0/64/bin/python SRW_I13d_intensity.py SRW.input')
        os.system(' /dls_sw/apps/python/anaconda/1.7.0/64/bin/python SRW_intensity.py SRW.input')
    
    elif calc_type == 'flux':
        print("FLUX CALCULATION - interactive mode ... ")
        os.system(' /dls_sw/apps/python/anaconda/1.7.0/64/bin/python SRW_flux.py SRW.input')
    
    elif calc_type == 'flux_cluster':
        print("FLUX CALCULATION - using the cluster ... ")
    #os.system(' /dls_sw/apps/python/anaconda/1.7.0/64/bin/python SRW_I13d_flux.py SRW.input')
        #### os.system(' /dls_sw/apps/python/anaconda/1.7.0/64/bin/python SRW_flux.py SRW.input')
        os.system('./submit_runbatch_Flux.sh')
    
        
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

    
