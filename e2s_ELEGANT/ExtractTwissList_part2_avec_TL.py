
import subprocess
import pandas as pd
import os
import csv
from shlex import split


# THIS is understod by PYTHON 3.6, but it's not what I want to do:
#blob=subprocess.check_output(['sdds2stream', ll[1], '-col=s'],stderr= subprocess.STDOUT).decode('UTF-8')


##########################################################################################################
##########################################################################################################
### FBT : the below block is the python equivalent of the lattice length calculated in the plott script
#from subprocess import Popen, PIPE
#from shlex import split
#p1 = Popen(split("sdds2stream %s -col=s" %ll[1]), stdout=PIPE)
#p2 = subprocess.check_output(split("tail -1"), stdin=p1.stdout).decode('UTF-8')
#
#print(p2)
#p3=float(p2)
#p4=p3*100
##########################################################################################################
##########################################################################################################

dd=pd.read_csv('blob.csv')



d32=dd[['Twi_Files_List']]



dft2=pd.DataFrame(columns=[['TwissName','LatticeLength']])




k=0

for nombre in range(0, len(d32)): #len(d32) is the total number in the list of Twiss Files
    fichier=d32.iloc[nombre]['Twi_Files_List']
    
    #df2=pd.read_csv(fichier) no !! coz fichier est un SDDS file pas un csv
    from subprocess import Popen, PIPE
    #from shlex import split
    p1 = Popen(split("sdds2stream %s -col=s" %fichier), stdout=PIPE)
    p2 = subprocess.check_output(split("tail -1"), stdin=p1.stdout).decode('UTF-8') 
    print('       ')
    print('    ')
    print(fichier)
    print(p2)
    try:
        p3=float(p2)
    except ValueError:
        p3=-3.1415926999
            
    dft2.loc[nombre,'TwissName']=fichier
    dft2.loc[nombre,'LatticeLength']=p3
    print(nombre)




dft2.to_csv('LatticesLengthsList.csv')




dfcsv = pd.read_csv('LatticesLengthsList.csv',float_precision='round_trip')

dfTemp2=dfcsv['LatticeLength']>561.54
dfTemp3=dfcsv['LatticeLength']<561.545

dfTemp4=dfcsv[dfTemp2]
diad=dfTemp4[dfTemp3]

diad.to_csv('diadLengths.csv')