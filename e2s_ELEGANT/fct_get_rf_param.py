#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri 13 Feb 2018

@author: MA

example of call is provided in the file MainTest_GetTwiss_and_DisplayTwiss.py
    
"""



import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import sys
import subprocess

os.system(' pwd ')




###########################################################################################################
def DisplayRF(name_of_File):
     
    fichier=name_of_File    
    cmdstring =("sddsprintout -spreadsheet=csv -col=* %s temp_table9876543210.csv" %fichier)
    os.system(cmdstring)      
    thisdir=os.getcwd()  
     
    Sz0_temp=subprocess.check_output(['sddsprintout',fichier, '-par=Sz0', '-notitle','-nolabel'],stderr= subprocess.STDOUT).decode('UTF-8')
    Sz0=float(Sz0_temp)

    print("***************************************************************")
    print("***************************************************************")
    print("sigma_z(0))    : ",Sz0)
    print("***************************************************************")
    print("***************************************************************")
  
      
    os.system('rm temp_table9876543210.csv')
############################################################################################################





def GetRF(name_of_File):

    fichier=name_of_File    
    cmdstring =("sddsprintout -spreadsheet=csv -col=* %s temp_table9876543210.csv" %fichier)
    os.system(cmdstring)      
    thisdir=os.getcwd()  
     
    Sz0_temp=subprocess.check_output(['sddsprintout',fichier, '-par=Sz0', '-notitle','-nolabel'],stderr= subprocess.STDOUT).decode('UTF-8')
    Sz0=float(Sz0_temp)

    return (Sz0 )  
    print("finished the retrieval. You can access the data with [],e.eg. ld[0] for s, etc.")
    os.system('rm temp9876543210.csv')
#############################################################################################################################################


















