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
#from srwlib import *
#from uti_plot import *
import os
import sys
import numpy as np

def twiss_at_location(fichier,locpoint):

    import glob
    import os
    import pandas as pd
    import matplotlib.pyplot as plt
    import sys
    import subprocess

    os.system(' pwd ')

    print('reading the position s argument...')
    #locpoint=float(sys.argv[2])

    print('position read...')


###############################################################################
    def closest(list, Number):
        aux = []
        for valor in list:
            aux.append(abs(Number-valor))

        return aux.index(min(aux))
###############################################################################


    #fichier=sys.argv[1]

    cmdstring =("sddsprintout -spreadsheet=csv -col=* %s temp_table01234567.csv" %fichier)
    os.system(cmdstring)


    thisdir=os.getcwd()

    df=pd.read_csv(os.path.join(thisdir,'temp_table01234567.csv'))
    dft=df[['s','betax','alphax','betay','alphay','etax','etaxp','ElementName','ElementType']]


    bb=closest(dft['s'],locpoint)


    print("closest s to the given location      : ",dft.loc[bb,'s'])
    print("index of this value in the twiss file: ", bb)
    print("betax     :", dft.loc[bb,'betax'])
    print("alphax    :", dft.loc[bb,'alphax'])
    print("betay     :", dft.loc[bb,'betay'])
    print("alphay    :", dft.loc[bb,'alphay'])
    print("etax     :", dft.loc[bb,'etax'])
    print("etaxp    :", dft.loc[bb,'etaxp'])
    print('data retrieved. process completed.')


    os.system('rm temp_table01234567.csv')
    return [dft.loc[bb,'betax'], dft.loc[bb,'alphax'], dft.loc[bb,'betay'], dft.loc[bb,'alphay'], dft.loc[bb,'etax'], dft.loc[bb,'etaxp']]



LATTICE = 'VMX' 
ELEdir  = 'e2s_LATTICES/'

# ---------------------------------
# elegant lattice type: LATTICE.lte 
# ---------------------------------

eLTE = LATTICE+'.lte'

# ----------------------------------
# elegant steering file: LATTICE.ele 
# ----------------------------------

here = os.getcwd() 
cmd  = here+'/'+ELEdir
os.chdir(cmd)
eELE = LATTICE+'.ele'

# ----------------------------------
# RUN elelgant
# ----------------------------------

cmd  = 'elegant '+eELE
os.system(cmd)

eTWI   = LATTICE+'.twi'

twiss  = twiss_at_location(eTWI, 100)
betax  = twiss[0]
alphax = twiss[1]
betay  = twiss[2]
alphay = twiss[3]
etax   = twiss[4]
etaxp  = twiss[5]


print ('Twiss = ', twiss)
cmd  = here
os.chdir(cmd)

# ----------------------------------
# 
# ----------------------------------


