# -*- coding: utf-8 -*-
#############################################################################
# SRWLIB Example#8: Simulating partially-coherent UR focusing with a CRL
# v 0.07
#############################################################################

from __future__ import print_function #Python 2.7 compatibility

import os
import sys
import datetime


###SRWLIB      = '/dls/physics/xph53246/source_to_beamline/SRW_Dev/env/work/SRW_PROJECT/MyBeamline/'
SRWLIB      = '/dls/physics/students/sug89938/E2S/SRWLIB/' # MA 12/03/2018 - repository created for pure SRWlib files 


sys.path.insert(0, SRWLIB)
from srwlib import *
from uti_plot import *
import numpy as np


print('SRWLIB Python Example # 6:')

INPUT_file = sys.argv[1]
infile = open(INPUT_file, 'r')

variables =[]; values=[];
for line in infile:
    # Typical line: variable = value
    variable, value = line.split('=')
    variable = variable.strip()  # remove leading/traling blanks
    value    = value.strip()
    variables.append(variable)
    values.append(value)

infile.close()

dict={}              # create a dictionary for easy access to variables 
 
for i in range(0,len(variables)) :
    dict[variables[i]] = values[i]

#
#***********Undulator
By_und  = float(dict['By_und'])
lam_und = float(dict['lam_und'])
Np_und  = float(dict['Np_und'])
K_und   = 0.9337 * By_und * lam_und * 100
#***********Beam Parameters
sig_x   = float(dict['sig_x'])
sig_y   = float(dict['sig_y'])
sig_xp  = float(dict['sig_xp'])
sig_yp  = float(dict['sig_yp'])
sigXX   = float(dict['sigXX'])
sigXXp   = float(dict['sigXXp'])
sigXpXp   = float(dict['sigXpXp'])
sigYY   = float(dict['sigYY'])
sigYYp   = float(dict['sigYYp'])
sigYpYp   = float(dict['sigYpYp'])
Ee      = float(dict['Ee'])
Ib      = float(dict['Ib'])
sigEperE = float(dict['dE'])
#***********BeamLine Parameters
slitZ   = float(dict['slitZ'])
slitDX  = float(dict['slitDX'])
slitDY  = float(dict['slitDY'])
Ephot_ini = float(dict['Ephot_ini'])
Ephot_end = float(dict['Ephot_end'])
BLname    = dict['IDname']
#**********Machine Parameters
meshXsta  = float(dict['meshXsta'])
meshXfin  = float(dict['meshXfin'])
meshYsta  = float(dict['meshYsta'])
meshYfin  = float(dict['meshYfin'])
meshEsta  = float(dict['meshEsta'])
meshEfin  = float(dict['meshEfin'])
Nelectr   = float(dict['Nelectr'])
outfil    = dict['outfil']
lattice   = dict['LATTICE']
calc_meth = float(dict['calc_meth'])
harm_1st  = float(dict['harm_1st'])
harm_last = float(dict['harm_last'])
#***********Extra Undulator Defs for Flux calculation
#harmB = SRWLMagFldH() #magnetic field harmonic 
#harmB.n = 1        # harmonic number
#harmB.h_or_v = 'v' # magnetic field plane: horzontal ('h') or vertical ('v')
#harmB.B = By_und   # 0.687566  #magnetic field amplitude [T]
#und = SRWLMagFldU([harmB])
#und.per = lam_und # 0.025 #period length [m]
#und.nPer = Np_und # number of periods (will be rounded to integer)
#magFldCnt = SRWLMagFldC([und], array('d', [0]), array('d', [0]), array('d', [0])) #Container of all magnetic field elements

print('+ ----------------------------------------------------------- ')
#print('| OUTPUT FILE: '+strExDataFolderName+'/'+strFluxOutFileName )
print('+ ----------------------------------------------------------- ')
print('| UNDULATOR PARAMETERS')
print('| By (T)           = '+str(By_und))
print('| lambda_u (m)     = '+str(lam_und))
print('| Np_u             = '+str(Np_und))
print('| K_und            = '+str(K_und))
print('+ ----------------------------------------------------------- ')
print('| BEAMLINE PARAMETERS')
print('| Z_slit (m)       = '+str(slitZ))
print('| DX (mm)          = '+str(dict['slitDX']))
print('| DY (mm)          = '+str(dict['slitDY']))
print('| Ephot_ini (eV)   = '+str(dict['Ephot_ini']))
print('| Ephot_end (eV)   = '+str(dict['Ephot_end']))
print('+ ----------------------------------------------------------- ')
print('| BEAM PARAMETERS @ centre of undulator')
print('| Ee       (GeV)   = '+str(Ee))
print('| sigma_x  (um)    = '+str(sig_x*1e6))
print('| sigma_xp (urad)  = '+str(sig_xp*1e6))
print('| sigma_y  (um)    = '+str(sig_y*1e6))
print('| sigma_yp (urad)  = '+str(sig_yp*1e6))
print('| dp/p             = '+str(sigEperE))
print('+ ----------------------------------------------------------- ')
print('| MOMENTS @ centre of undulator')
print('| sigXX     = '+str(sigXX))
print('| sigXXp    = '+str(sigXXp))
print('| sigXpXp   = '+str(sigXpXp))
print('| sigYY     = '+str(sigYY))
print('| sigYYp    = '+str(sigYYp))
print('| sigYpYp   = '+str(sigYpYp))
print('+ ----------------------------------------------------------- ')
print('| MACHINE PARAMETERS')
print('| Ib       (mA)   = '+str(Ib*1000))
print('+ ----------------------------------------------------------- ')
print('|COMPUTATION PARAMETERS')
print('| Nelectr                   = '+str(Nelectr))
print('| mesh X start       (mm)   = '+str(meshXsta))
print('| mesh X fin         (mm)   = '+str(meshXfin))
print('| mesh Y start       (mm)   = '+str(meshYsta))
print('| mesh Y fin         (mm)   = '+str(meshYfin))
print('| mesh E start       (eV)   = '+str(meshEsta))
print('| mesh E fin         (eV)   = '+str(meshEfin))
print('| outfil                    = '+outfil)
print('+ ----------------------------------------------------------- ')

#****************************Input Parameters:
strExDataFolderName =  outfil # 'data_example_08' #example data sub-folder name
os.system('[ -d '+outfil+' ] && echo "output directory exists ..." || mkdir '+outfil)
timestamp = "{:%Y-%b-%d_%H:%M:%S}".format(datetime.datetime.now())
strPowOutFileName = 'power__'+outfil+'__'+lattice+'_'+timestamp+'.dat' #file name for output power density data
#strExDataFolderName = 'SRW_I04' #example data sub-folder name

#strPowOutFileName = 'powerdensity_I04.dat' #file name for output power density data


#***********Undulator
numPer = Np_und # 73 #Number of ID Periods (without counting for terminations
undPer = lam_und # 0.033  # lam_und #Period Length [m]
Bx   = 0 #Peak Horizontal field [T]
By   = By_und # 0.3545 #Peak Vertical field [T]
phBx = 0 #Initial Phase of the Horizontal field component
phBy = 0 #Initial Phase of the Vertical field component
sBx  = 1 #Symmetry of the Horizontal field component vs Longitudinal position
sBy  = 1 #Symmetry of the Vertical field component vs Longitudinal position
xcID = 0 #Transverse Coordinates of Undulator Center [m]
ycID = 0
zcID = 0#lam_und*Np_und/2*1.055 
# Longitudinal Coordinate of Undulator Center wit hrespect to Straight Section Center [m]
# my understanding: you need to calculate from a point outside the undulator
# e.g. SIREPO fixes a -1.2705m offset for an undulator of 2.409m which 
# hence: 2.409/2*1.055 = 1.2707

print(zcID)

harmNum = 1 # default=1 for an undulator 
und = SRWLMagFldU([SRWLMagFldH(harmNum, 'v', By, phBy, sBy, 1), SRWLMagFldH(harmNum, 'h', Bx, phBx, sBx, 1)], undPer, numPer) #Ellipsoidal Undulator
magFldCnt = SRWLMagFldC([und], array('d', [xcID]), array('d', [ycID]), array('d', [zcID])) #Container of all Field Elements


#***********Electron Beam
elecBeam = SRWLPartBeam()
elecBeam.Iavg = Ib # 0.1 #Average Current [A]
elecBeam.partStatMom1.x = 0.00 #Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on) [m]
elecBeam.partStatMom1.y = 0.00
elecBeam.partStatMom1.z = -0.5*undPer*(numPer + 4) #Initial Longitudinal Coordinate (set before the ID)
elecBeam.partStatMom1.xp = 0 #Initial Relative Transverse Velocities
elecBeam.partStatMom1.yp = 0
elecBeam.partStatMom1.gamma = Ee/0.51099890221e-03 #Relative Energy
#2nd order statistical moments
elecBeam.arStatMom2[0]  = sigXX    # (sig_x)**2 # (5*118.027e-06)**2 #<(x-x0)^2> 
elecBeam.arStatMom2[1]  = sigXXp
elecBeam.arStatMom2[2]  = sigXpXp  # (sig_xp)**2 #(27.3666e-06)**2 #<(x'-x'0)^2>
elecBeam.arStatMom2[3]  = sigYY    #(sig_y)**2  # (15.4091e-06)**2 #<(y-y0)^2>
elecBeam.arStatMom2[4]  = sigYYp
elecBeam.arStatMom2[5]  = sigYpYp  #(sig_yp)**2 # (2.90738e-06)**2 #<(y'-y'0)^2>
elecBeam.arStatMom2[10] = (sigEperE)**2 #<(E-E0)^2>/E0^2

#***********Precision Parameters for SR calculation
meth = calc_meth # 1 #SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
relPrec = 0.001 #relative precision
zStartInteg = 0 #longitudinal position to start integration (effective if < zEndInteg)
zEndInteg = 0 #longitudinal position to finish integration (effective if > zStartInteg)
npTraj = 20000 #Number of points for trajectory calculation 
useTermin = 0 #1 #Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
sampFactNxNyForProp = 0.25*2 #sampling factor for adjusting nx, ny (effective if > 0)
arPrecPar = [meth, relPrec, zStartInteg, zEndInteg, npTraj, useTermin, 0]



########################## power density ###########


arPrecP = [0]*5 #for power density
arPrecP[0] = 1.5 #precision factor
arPrecP[1] = 1 #power density computation method (1- "near field", 2- "far field")
arPrecP[2] = 0 #initial longitudinal position (effective if arPrecP[2] < arPrecP[3])
arPrecP[3] = 0 #final longitudinal position (effective if arPrecP[2] < arPrecP[3])
arPrecP[4] = 200000 #number of points for (intermediate) trajectory calculatio




#if float(dict['slitDX']) > 0:
####### with aperture ###########
stkP = SRWLStokes() #for spectral flux vs photon energy
stkP.mesh.zStart = float(dict['slitZ'])     # 12.9 # 30. #longitudinal position [m] at which UR has to be calculated
stkP.allocate(1,100,100)  # (19001, 1, 1) #numbers of points vs photon energy, horizontal and vertical positions
RangeX = float(dict['slitDX'])#*K_und/elecBeam.partStatMom1.gamma
RangeY = float(dict['slitDY'])#/elecBeam.partStatMom1.gamma
stkP.mesh.xStart = -RangeX/2e6 #initial horizontal position [m]
stkP.mesh.xFin   =  RangeX/2e6 #final horizontal position [m]
stkP.mesh.yStart = -RangeY/2e6 #initial vertical position [m]
stkP.mesh.yFin   =  RangeY/2e6 #final vertical position [m] 

meshInitPartCoh  = deepcopy(stkP.mesh)
stkP.partBeam = elecBeam


print('   Performing Power Density calculation (from field) ... ', end='')
srwl.CalcPowDenSR(stkP, elecBeam, 0, magFldCnt, arPrecP)
#srwl.CalcStokesUR(stkP, elecBeam, und, arPrecP)
#srwl.CalcPowDenSR(stkPx, eBeam, 0, magFldCnt, arPrecP)
#srwl.CalcPowDenSR(stkPy, eBeam, 0, magFldCnt, arPrecP)
powTot =uti_math.integ_ar_2d(stkP.arS, 1,[stkP.mesh.xStart, stkP.mesh.xFin, stkP.mesh.nx], [stkP.mesh.yStart, stkP.mesh.yFin, stkP.mesh.ny])*1.e+06
#powInAp =uti_math.integ_ar_2d(stkP.arS, 1,[stkP.mesh.xStart, stkP.mesh.xFin, stkP.mesh.nx], [stkP.mesh.yStart, stkP.mesh.yFin, stkP.mesh.ny],[float(dict['slitDX']), float(dict['slitDX']), stkP.mesh.nx], [float(dict['slitDY']), float(dict['slitDY']), stkP.mesh.ny])*1.e+06

srwl_uti_save_intens_ascii(stkP.arS, stkP.mesh, os.path.join(os.getcwd(), strExDataFolderName, strPowOutFileName), 0, ['', 'Horizontal Position', 'Vertical Position', 'Power Density'], _arUnits=['', 'm', 'm', 'W/mm^2'])

print('   Saving power density data to file ... '+strExDataFolderName+'/'+strPowOutFileName+'... ', end='')
print('Power ~total:', powTot, 'W')
#print('Power within work aperture:', powInAp, 'W')
print('done')




 
