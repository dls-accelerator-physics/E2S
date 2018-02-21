# -*- coding: utf-8 -*-
#############################################################################
# SRWLIB Example#8: Simulating partially-coherent UR focusing with a CRL
# v 0.07
#############################################################################

from __future__ import print_function #Python 2.7 compatibility

import os
import sys
import numpy as np
import datetime

SRWLIB      = '/dls/physics/xph53246/source_to_beamline/SRW_Dev/env/work/SRW_PROJECT/MyBeamline/'
sys.path.insert(0, SRWLIB)
from srwlib import *
from uti_plot import *

print('SRWLIB Python Example # 8:')
print('Simulating emission and propagation of undulator radiation (UR) wavefront through a simple optical scheme including CRL')
print('')
print('First, single-electron UR (on-axis spectrum and a wavefront at a fixed photon energy) is calculated and propagated through the optical scheme. ', end='')
print('After this, calculation of partially-coherent UR from entire electron beam is started as a loop over "macro-electrons", using "srwl_wfr_emit_prop_multi_e" function. ', end='')
print('This function can run either in "normal" sequential mode, or in parallel mode under "mpi4py".', end='')
print('For this, an MPI2 package and the "mpi4py" Python package have to be installed and configured, and this example has to be started e.g. as:')
print('    mpiexec -n 5 python SRWLIB_Example08.py')
print('For more information on parallel calculations under "mpi4py" please see documentation to the "mpi4py" and MPI.')
print('Note that the long-lasting partially-coherent UR calculation saves from time to time instant average intensity to an ASCII file, ', end='')
print('so the execution of the long loop over "macro-electrons" can be aborted after some time without the danger that all results will be lost.')
print('')

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
K_und   = 0.9338 * By_und * lam_und * 100
#***********Beam Parameters
sig_x   = float(dict['sig_x'])
sig_y   = float(dict['sig_y'])
sig_xp  = float(dict['sig_xp'])
sig_yp  = float(dict['sig_yp'])
Ee      = float(dict['Ee'])
Ib      = float(dict['Ib'])
sigEperE = float(dict['dE'])
#***********BeamLine Parameters
slitZ   = float(dict['slitZ'])
slitDX  = float(dict['slitDX'])
slitDY  = float(dict['slitDY'])
Ephot_ini = float(dict['Ephot_ini'])
Ephot_end = float(dict['Ephot_end'])
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
#***********Extra Undulator Defs for Flux calculation
harmB = SRWLMagFldH() #magnetic field harmonic 
harmB.n = 1        # harmonic number
harmB.h_or_v = 'v' # magnetic field plane: horzontal ('h') or vertical ('v')
harmB.B = By_und   # 0.687566  #magnetic field amplitude [T]
und = SRWLMagFldU([harmB])
und.per = lam_und # 0.025 #period length [m]
und.nPer = Np_und # number of periods (will be rounded to integer)
magFldCnt = SRWLMagFldC([und], array('d', [0]), array('d', [0]), array('d', [0])) #Container of all magnetic field elements

print('+ ----------------------------------------------------------- ')
#print('| OUTPUT FILE: '+strExDataFolderName+'/'+strFluxOutFileName )
print('+ ----------------------------------------------------------- ')
print('| UNDULATOR PARAMETERS')
print('| By (T)           = '+str(By_und))
print('| lambda_u (m)     = '+str(lam_und))
print('| Np_u             = '+str(Np_und))
print('| K_max            = '+str(K_und))
print('+ ----------------------------------------------------------- ')
print('| BEAMLINE PARAMETERS')
print('| Z_slit (m)       = '+str(slitZ))
print('| DX (um)          = '+str(dict['slitDX']))
print('| DY (um)          = '+str(dict['slitDY']))
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
print('| MACHINE PARAMETERS')
print('| Ib       (mA)   = '+str(Ib*1000))
print('+ ----------------------------------------------------------- ')
print('|COMPUTATION PARAMETERS')
print('| Nelectr                   = '+str(Nelectr))
print('| mesh X start       (um)   = '+str(meshXsta))
print('| mesh X fin         (um)   = '+str(meshXfin))
print('| mesh Y start       (um)   = '+str(meshYsta))
print('| mesh Y fin         (um)   = '+str(meshYfin))
print('| mesh E start       (eV)   = '+str(meshEsta))
print('| mesh E fin         (eV)   = '+str(meshEfin))
print('| outfil                    = '+outfil)
print('+ ----------------------------------------------------------- ')

#****************************Input Parameters:
strExDataFolderName =  outfil # 'data_example_08' #example data sub-folder name
os.system('[ -d '+outfil+' ] && echo "output directory exists ..." || mkdir '+outfil)
timestamp = "{:%Y-%b-%d_%H:%M:%S}".format(datetime.datetime.now())
strFluxOutFileName  =  'flux__'+outfil+'__'+lattice+'_'+timestamp+'.dat' #'ex08_flux.dat' #file name for output UR flux data

#***********Undulator
numPer = Np_und # 73 #Number of ID Periods (without counting for terminations
undPer = lam_und # 0.033  # lam_und #Period Length [m]
Bx = 0 #Peak Horizontal field [T]
By = By_und # 0.3545 #Peak Vertical field [T]
phBx = 0 #Initial Phase of the Horizontal field component
phBy = 0 #Initial Phase of the Vertical field component
sBx = 1 #Symmetry of the Horizontal field component vs Longitudinal position
sBy = -1 #Symmetry of the Vertical field component vs Longitudinal position
xcID = 0 #Transverse Coordinates of Undulator Center [m]
ycID = 0
zcID = 1.25 # 0 #Longitudinal Coordinate of Undulator Center wit hrespect to Straight Section Center [m]

und = SRWLMagFldU([SRWLMagFldH(1, 'v', By, phBy, sBy, 1), SRWLMagFldH(1, 'h', Bx, phBx, sBx, 1)], undPer, numPer) #Ellipsoidal Undulator
magFldCnt = SRWLMagFldC([und], array('d', [xcID]), array('d', [ycID]), array('d', [zcID])) #Container of all Field Elements


#***********Electron Beam
elecBeam = SRWLPartBeam()
elecBeam.Iavg = Ib # 0.1 #Average Current [A]
elecBeam.partStatMom1.x = 0.00 #Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on) [m]
elecBeam.partStatMom1.y = 0.00
elecBeam.partStatMom1.z = 0. #-0.5*undPer*(numPer + 4) #Initial Longitudinal Coordinate (set before the ID)
elecBeam.partStatMom1.xp = 0 #Initial Relative Transverse Velocities
elecBeam.partStatMom1.yp = 0
elecBeam.partStatMom1.gamma = Ee/0.51099890221e-03 #Relative Energy
#2nd order statistical moments
elecBeam.arStatMom2[0] = (sig_x)**2 # (5*118.027e-06)**2 #<(x-x0)^2>
elecBeam.arStatMom2[1] = 0
elecBeam.arStatMom2[2] = (sig_xp)**2 #(27.3666e-06)**2 #<(x'-x'0)^2>
elecBeam.arStatMom2[3] = (sig_y)**2  # (15.4091e-06)**2 #<(y-y0)^2>
elecBeam.arStatMom2[4] = 0
elecBeam.arStatMom2[5] = (sig_yp)**2 # (2.90738e-06)**2 #<(y'-y'0)^2>
elecBeam.arStatMom2[10] = (sigEperE)**2 #<(E-E0)^2>/E0^2

#***********Precision Parameters for SR calculation
meth = 1 #SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
relPrec = 0.01 #relative precision
zStartInteg = 0 #longitudinal position to start integration (effective if < zEndInteg)
zEndInteg = 0 #longitudinal position to finish integration (effective if > zStartInteg)
npTraj = 20000 #Number of points for trajectory calculation 
useTermin = 0 #1 #Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
sampFactNxNyForProp = 0.25 #sampling factor for adjusting nx, ny (effective if > 0)
arPrecPar = [meth, relPrec, zStartInteg, zEndInteg, npTraj, useTermin, 0]

########################## flux through slit ###########
arPrecF = [0]*5 #for spectral flux vs photon energy
arPrecF[0] = 1 #initial UR harmonic to take into account
arPrecF[1] = 21 #final UR harmonic to take into account
arPrecF[2] = 1.5 * 2#longitudinal integration precision parameter
arPrecF[3] = 1.5 * 2#azimuthal integration precision parameter
arPrecF[4] = 1 #calculate flux (1) or flux per unit surface (2)

stkF = SRWLStokes() #for spectral flux vs photon energy
stkF.mesh.zStart = float(dict['slitZ'])     # 12.9 # 30. #longitudinal position [m] at which UR has to be calculated
stkF.mesh.eStart = float(dict['Ephot_ini']) # 1000. #initial photon energy [eV]
stkF.mesh.eFin   = float(dict['Ephot_end']) # 30000. #final photon energy [eV]
stkF.allocate(int(stkF.mesh.eFin-stkF.mesh.eStart+1),1,1)  # (19001, 1, 1) #numbers of points vs photon energy, horizontal and vertical positions
stkF.mesh.xStart = -float(dict['slitDX'])/2.0/1.e6 # -0.0001 #initial horizontal position [m]
stkF.mesh.xFin   =  float(dict['slitDX'])/2.0/1.e6 #  0.0001 #final horizontal position [m]
stkF.mesh.yStart = -float(dict['slitDY'])/2.0/1.e6 # -0.00006 #initial vertical position [m]
stkF.mesh.yFin   =  float(dict['slitDY'])/2.0/1.e6 #  0.00006 #final vertical position [m]    

print('   Performing Spectral Flux (Stokes parameters) calculation ... ', end='')
srwl.CalcStokesUR(stkF, elecBeam, und, arPrecF)
print('   Saving intensity data to file ... ', end='')
srwl_uti_save_intens_ascii(stkF.arS, stkF.mesh, os.path.join(os.getcwd(), strExDataFolderName, strFluxOutFileName), 0, ['Photon Energy', '', '', 'Flux'], _arUnits=['eV', '', '', 'ph/s/.1%bw'])

print('done')
##########################


 
