#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri 02 Feb 2018

@author: FBT

example of call is provided in the file MainTest_GetTwiss_and_DisplayTwiss.py
    
"""



import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import sys
import subprocess

os.system(' pwd ')



##################################################################################################################
# FBT: the little piece of code below ("closest" function) is taken from 
# https://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value
##################################################################################################################
def closest(list, Number): 
    aux = []
    for valor in list:
        aux.append(abs(Number-valor))

    return aux.index(min(aux))
##################################################################################################################





###########################################################################################################################################
def DisplayTwiss(name_of_File,s_loc):
    print('reading the position s argument...')
    locpoint=float(s_loc) #cast to float because the argument type is string, even if we type a number    
    print('position read...')      
    fichier=name_of_File    
    cmdstring =("sddsprintout -spreadsheet=csv -col=* %s temp_table01234567.csv" %fichier)
    os.system(cmdstring)      
    thisdir=os.getcwd()  
    df=pd.read_csv(os.path.join(thisdir,'temp_table01234567.csv'))
    dft=df[['s','betax','alphax','betay','alphay','etax','etaxp','ElementName','ElementType']]   
    bb=closest(dft['s'],locpoint)  
    ex0_temp=subprocess.check_output(['sddsprintout',fichier, '-par=ex0', '-notitle','-nolabel'],stderr= subprocess.STDOUT).decode('UTF-8')
    ex0=float(ex0_temp)
    Sdelta0_temp=subprocess.check_output(['sddsprintout',fichier, '-par=Sdelta0', '-notitle','-nolabel'],stderr= subprocess.STDOUT).decode('UTF-8')
    Sdelta0=float(Sdelta0_temp)  
    pCentral_temp=subprocess.check_output(['sddsprintout',fichier, '-par=pCentral', '-notitle','-nolabel'],stderr= subprocess.STDOUT).decode('UTF-8') # added, to take into account scaling between source energy and lattice energy
    m_e = 0.5109989461 
    pCentral=float(pCentral_temp) * m_e / 1e3  
    print("***************************************************************")
    print("***************************************************************")
    print("closest s to the requested location    : ",dft.loc[bb,'s'])
    print("index of this value in the twiss file  : ", bb)
    print("***************************************************************")
    print("***************************************************************")
    print("                                   ")
    print(" Twiss parameters at that location:")
    print(" ----------------------------------")
    print(" betax     :", dft.loc[bb,'betax'])
    print(" alphax    :", dft.loc[bb,'alphax'])
    print(" betay     :", dft.loc[bb,'betay'])
    print(" alphay    :", dft.loc[bb,'alphay'])
    print(" etax      :", dft.loc[bb,'etax'])
    print(" etaxp     :", dft.loc[bb,'etaxp'])
    print("                  ")
    print(" Global parameters:")
    print(" ------------------")
    print(" the x-emittance ex0 is       :", ex0)
    print(" the energy spread Sdelta0    :", Sdelta0)
    print(" the e-momentum is            :", pCentral," (GeV/c)")
    print("                  ")
    print("                  ")
    print(' data retrieved. process completed.')
      
    os.system('rm temp_table01234567.csv')
###############################################################################################################################################








###############################################################################################################################################
def GetTwissList(name_of_File,s_loc):

    print('reading the position s argument...')
    locpoint=float(s_loc) #cast to float because the argument type is string, even if we type a number   
    print('position read...')
    fichier=name_of_File    
    cmdstring =("sddsprintout -spreadsheet=csv -col=* %s temp_table01234567.csv" %fichier)
    os.system(cmdstring)       
    thisdir=os.getcwd()  
    df=pd.read_csv(os.path.join(thisdir,'temp_table01234567.csv'))
    dft=df[['s','betax','alphax','betay','alphay','etax','etaxp','ElementName','ElementType']]
    ########################################################################################################################
    # FBT : please note that writing:
    #dft=df[['s','betax','alphax','betay','alphay','etax','etaxp','Sdelta0','ex0','ElementName','ElementType']]
    # is wrong because the emittance and the energy spread are encoded as SDDS parameters and not as column. Therefore they will
    # be extracted using python's subprocess method, below
    #########################################################################################################################
    
    bb=closest(dft['s'],locpoint)
    ex0_temp=subprocess.check_output(['sddsprintout',fichier, '-par=ex0', '-notitle','-nolabel'],stderr= subprocess.STDOUT).decode('UTF-8')
    ex0=float(ex0_temp)
    Sdelta0_temp=subprocess.check_output(['sddsprintout',fichier, '-par=Sdelta0', '-notitle','-nolabel'],stderr= subprocess.STDOUT).decode('UTF-8')
    Sdelta0=float(Sdelta0_temp)
    pCentral_temp=subprocess.check_output(['sddsprintout',fichier, '-par=pCentral', '-notitle','-nolabel'],stderr= subprocess.STDOUT).decode('UTF-8')# added, to take into account scaling between source energy and lattice energy
    m_e = 0.5109989461
    pCentral=float(pCentral_temp) * m_e / 1e3
    print(" >>>>>>>>>>>>>>>>>>>>>>>>>>>>> P (GeV) =",pCentral)
    return (dft.loc[bb,'s'], bb,dft.loc[bb,'betax'],dft.loc[bb,'alphax'],dft.loc[bb,'betay'],dft.loc[bb,'alphay'],dft.loc[bb,'etax'],dft.loc[bb,'etaxp'],ex0,Sdelta0,pCentral )  # pCentral added [MA 29/10/2018]
    print("finished the retrieval. You can access the data with [],e.eg. ld[0] for s, etc.")
    os.system('rm temp_table01234567.csv')
#############################################################################################################################################





















