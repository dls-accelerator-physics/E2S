#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri 02 Feb 2018

@author: mfc33124

how to use: retrieve_twiss.py nameofthetwissfile numberlocation
example: python retrieve_twiss.py VMX.twi 23
    
"""

import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import sys
import subprocess

def main():
    os.system(' pwd ')

    print('reading the position s argument...')
    locpoint=float(sys.argv[2])

    print('position read...')


###############################################################################
    def closest(list, Number):
        aux = []
        for valor in list:
            aux.append(abs(Number-valor))

            return aux.index(min(aux))
###############################################################################


        fichier=sys.argv[1]

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

        return [1,2,3]



















