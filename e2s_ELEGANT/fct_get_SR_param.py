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
from shlex import split                                                              

os.system(' pwd ')




# -----------------------------------------------------------------------------------------------------------
def DisplayRF(rootname_of_File):
     
    fichier=rootname_of_File+'.rf'    
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
# ------------------------------------------------------------------------------------------------------------




# ------------------------------------------------------------------------------------------------------------
def GetRF(rootname_of_File):

    fichier=rootname_of_File+'.rf'    
    cmdstring =("sddsprintout -spreadsheet=csv -col=* %s temp_table9876543210.csv" %fichier)
    os.system(cmdstring)      
    thisdir=os.getcwd()  
     
    Sz0_temp=subprocess.check_output(['sddsprintout',fichier, '-par=Sz0', '-notitle','-nolabel'],stderr= subprocess.STDOUT).decode('UTF-8')
    Sz0=float(Sz0_temp)

    return (Sz0 )  
    print("finished the retrieval. You can access the data with [],e.eg. ld[0] for s, etc.")
    os.system('rm temp9876543210.csv')
# ------------------------------------------------------------------------------------------------------------



# -----------------------------------------------------------------------------------------------------------
def DisplayCirc(rootname_of_File):
    fichier = rootname_of_File+'.mag'
    p1 = subprocess.Popen(split("sdds2stream %s -col=s" %fichier), stdout=subprocess.PIPE)
    Circ_temp = subprocess.check_output(split("tail -1"), stdin=p1.stdout).decode('UTF-8')
    
    try:
        Circ = float(Circ_temp)
    except ValueError:
        Circ = -3.1415926999



    print("***************************************************************")
    print("Circumference (m)    : ",Circ)
    print("***************************************************************")
# ------------------------------------------------------------------------------------------------------------




# ------------------------------------------------------------------------------------------------------------
def GetCirc(rootname_of_File):
    fichier = rootname_of_File+'.mag'
    p1 = subprocess.Popen(split("sdds2stream %s -col=s" %fichier), stdout=subprocess.PIPE)
    Circ_temp = subprocess.check_output(split("tail -1"), stdin=p1.stdout).decode('UTF-8')

    try:
        Circ = float(Circ_temp)
    except ValueError:
        Circ = -3.1415926999

        
    return (Circ)  
# ------------------------------------------------------------------------------------------------------------




















