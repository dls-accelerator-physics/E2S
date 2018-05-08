from __future__ import print_function #Python 2.7 compatibility
import os
import sys
import numpy as np
import scipy.optimize as opt
import subprocess


e2s_ROOT = '/dls/physics/xph53246/source_to_beamline/E2S/'
sys.path.insert(0, e2s_ROOT)

here = os.getcwd() # memorize the TOP directory

def run_E2S(x,*names):
    #
    # change parameters in steering file
    # x = input parameters
    # names = names of the input parameters
    #### sf = 'VMX__I20SCA__SHA.input'
    sf = 'DTBA_C1a_AA__I20SCA__SHA.input'

    
    os.chdir(e2s_ROOT)
    nsf = alter_steer_file(sf, x, names) # nsf = new steer file 
    #cmd = 'python E2S.py '+nsf
    #os.system(cmd)
    subprocess.check_output(['python','E2S.py',nsf],stderr= subprocess.STDOUT)
    os.chdir(here)
    f = open('/dls/physics/xph53246/source_to_beamline/E2S/e2s_SHADOW/objectives','r')
#    os.system('echo '+str(x[0])+' '+str(x[1])+' '+str(x[2])+' '+str(x[3])+
#              ' '+str(x[4])+' '+str(x[5])+'  > /dls/physics/xph53246/source_to_beamline/E2S/e2s_SHADOW/params')
    os.system('echo '+str(x[0])+' '+str(x[1])+' > /dls/physics/xph53246/source_to_beamline/E2S/e2s_SHADOW/params')
    aveX=[]
    aveY=[]
    sigmaX=[]
    sigmaY=[]
    penX  =[] # penaltyHorizontal
    while True:
        linea     = f.readline().strip()
        if linea == '':
            break
        else:
            aveX.append(float(linea.split()[0]))
            aveY.append(float(linea.split()[1]))
            sigmaX.append(float(linea.split()[2]))
            sigmaY.append(float(linea.split()[3]))
            penX.append(float(linea.split()[4]))
    

#    obj = sigmaX[0] # need to read the sigmaX 
    obj = penX[0]   
    return obj


def alter_steer_file(sf, params, names):
    # here you need to identify the line with params[0]
    # and substitute the value
    nsf = '_'+sf
#    with open(sf,'r') as input_file, open(nsf,'w') as output_file:
#            for line in input_file:
#                L = line.split()[0]
#                if L == names:
#                    output_file.write(line.split()[0]+' '+line.split()[1]+'     '+str(params[0])+'\n')
#                else:
#                    output_file.write(line)
    with open(sf,'r') as input_file, open(nsf,'w') as output_file:
            for line in input_file:
                L = line.split()[0]
                if L == names[0]:
                    output_file.write(line.split()[0]+' '+line.split()[1]+'     '+str(params[0])+'\n')
                elif L == names[1]:
                    output_file.write(line.split()[0]+' '+line.split()[1]+'     '+str(params[1])+'\n')
#                elif L == names[2]:
#                    output_file.write(line.split()[0]+' '+line.split()[1]+'     '+str(params[2])+'\n')
#                elif L == names[3]:
#                    output_file.write(line.split()[0]+' '+line.split()[1]+'     '+str(params[3])+'\n')
#                elif L == names[4]:
#                    output_file.write(line.split()[0]+' '+line.split()[1]+'     '+str(params[4])+'\n')
#                elif L == names[5]:
#                    output_file.write(line.split()[0]+' '+line.split()[1]+'     '+str(params[5])+'\n')
                       
                else:
                    output_file.write(line)  
    
    return nsf

def main():
    # 1) change parameters in the steering file
    # 2) run E2S.py 
    # 3) read the results of the run
    
    #run_E2S() # this is working 
    ####### val0 =  9.2 # value for RMIRR
    ####### arg0 = 'RMIRR'
    #### val0 = [ 8.75, 7600., 23., 240.]  ##5000., 15.] # RMIRR_OE7, AXMAJ_OE8, AXmin_OE8
    #### arg0 = 'RMIRR','AXMAJ','AXMIN','T_IMAGE'
    ### val0 = [ 2350., 2000., 8.7435897, 2000, 2300., 150.] # SSOUR_OE1, SIMAG_OE1, RMIRR_OE7, SSOUR_OE8, SSOUR_OE7, T_IMAGE_OE10
    ### arg0 = 'SSOUR_OE1','SIMAG_OE1','RMIRR_OE7','SSOUR_OE8','SIMAG_OE8','T_IMAGE_OE10'
    val0 = [ 8.7435897, 340.] # SSOUR_OE1, SIMAG_OE1, RMIRR_OE7, SSOUR_OE8, SSOUR_OE7, T_IMAGE_OE10
    arg0 = 'RMIRR_OE7','T_IMAGE_OE10' # OE1 and OE8 alter only the vertical compnent, as such they are useless for our goal

    objectives = opt.minimize(run_E2S, val0,args=arg0, method='Nelder-Mead') #,options={'fatol':0.001}) 
    print('OBJECTIVES '+str(objectives.fun))
    # objectives = run_E2S(params)
    

if __name__ == '__main__':
    main()
