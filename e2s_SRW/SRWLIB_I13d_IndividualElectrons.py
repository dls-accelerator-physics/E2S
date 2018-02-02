# -*- coding: utf-8 -*-
#############################################################################
# SRWLIB Example#8: Simulating partially-coherent UR focusing with a CRL
# v 0.07
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
from uti_plot import *
import os
import sys
import numpy as np


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


#infile = open('DI_mb_input.txt', 'r')
infile = open(INPUT_file, 'r')

variables =[]; values=[];
for line in infile:
    # Typical line: variable = value
    variable, value = line.split('=')
    variable = variable.strip()  # remove leading/traling blanks
    value    = value.strip()
    variables.append(variable)
    values.append(value)

#print(variables)
#print(values)
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
print('| COMPUTATION PARAMETERS')
print('| mesh X start       (um)   = '+str(meshXsta))
print('| mesh X fin         (um)   = '+str(meshXfin))
print('| mesh Y start       (um)   = '+str(meshYsta))
print('| mesh Y fin         (um)   = '+str(meshYfin))
print('| mesh E start       (eV)   = '+str(meshEsta))
print('| mesh E fin         (eV)   = '+str(meshEfin))
print('+ ----------------------------------------------------------- ')

#****************************Input Parameters:
strExDataFolderName = 'data_example_08' #example data sub-folder name
strFluxOutFileName  = 'ex08_flux.dat' #file name for output UR flux data
strTrajOutFileName = 'ex08_res_traj.dat' #file name for output trajectory data
strIntOutFileName1 = 'ex08_res_int1.dat' #file name for output SR intensity data
strIntOutFileName2 = 'ex08_res_int2.dat' #file name for output SR intensity data
strIntOutFileName3 = 'ex08_res_int3.dat' #file name for output SR intensity data
strIntOutFileNamePartCoh = 'ex08_res_int_part_coh1.dat' #file name for output SR intensity data
strOpTrFileName = 'ex08_res_op_tr_crl.dat' #file name for output optical transmission data
strOpPathDifFileName = 'ex08_res_op_path_dif_crl.dat' #file name for optical path difference data

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

#***********Initial Wavefront data placeholder
wfr2 = SRWLWfr() #For intensity distribution at fixed photon energy
wfr2.allocate(1, 101, 101) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr2.mesh.zStart = 36.25 + 1.25 #Longitudinal Position [m] from Center of Straight Section at which SR has to be calculated

wfr2.mesh.zStart = 12.9 + 1.25 #Longitudinal Position [m] from Center of Straight Section at which SR has to be calculated

wfr2.mesh.eStart =    meshEsta # 8830 #Initial Photon Energy [eV] <<< 1st harm peak for a 7GeV machine
wfr2.mesh.eFin   =    meshEfin # 8830 #Final Photon Energy [eV]
wfr2.mesh.xStart = -0.0015  # -0.0015 #Initial Horizontal Position [m]
wfr2.mesh.xFin =    0.0015  # 0.0015 #Final Horizontal Position [m]
wfr2.mesh.yStart = -0.0006  # -0.0006 #Initial Vertical Position [m]
wfr2.mesh.yFin   =  0.0006  # 0.0006 #Final Vertical Position [m]
meshInitPartCoh = deepcopy(wfr2.mesh)
wfr2.partBeam = elecBeam


#***********Optical Elements and Propagation Parameters
fx = 1e+23 #Focal Length in Horizontal plane
fy = 1e+23 # 19.0939 #Focal Length in Vertical plane
optLens = SRWLOptL(fx, fy) #Ideal Lens

delta = 4.3712962E-06 #Refractive index decrement of Be at 8830 eV
attenLen = 6946.13E-06 #[m] Attenuation length of Be at 8830 eV
geomApertF = 1E-03 #[m] Geometrical aparture of 1D CRL in the Focusing plane
geomApertNF = 1E-03 #[m] Geometrical aparture of 1D CRL in the plane where there is no focusing
rMin = 0.5E-03 #[m] radius at tip of parabola of CRL
nCRL = 3
wallThick = 50E-06 #[m] wall thickness of CRL

optCRL = srwl_opt_setup_CRL(2, delta, attenLen, 1, geomApertNF, geomApertF, rMin, nCRL, wallThick, 0, 0) #1D CRL
print('Saving CRL transmission data to files (for viewing/debugging)...', end='')
optTrIntCRL = optCRL.get_data(2, 3)
#srwl_uti_save_intens_ascii(optTrIntCRL, optCRL.mesh, os.path.join(os.getcwd(), strExDataFolderName, strOpTrFileName), 0, ['', 'Horizontal Position', 'Vertical Position', 'Intensity Transmission'], _arUnits=['', 'm', 'm', 'r.u.'])

optPathDifCRL = optCRL.get_data(3, 3)
#srwl_uti_save_intens_ascii(optPathDifCRL, optCRL.mesh, os.path.join(os.getcwd(), strExDataFolderName, strOpPathDifFileName), 0, ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'], _arUnits=['', 'm', 'm', 'm'])
#print('done')

optApert = SRWLOptA('r', 'a', geomApertNF, geomApertF) #Aperture

# optDrift = SRWLOptD(38.73) #Drift space
# optDrift = SRWLOptD(18.37) #Drift space
# optDrift = SRWLOptD(0.1) #Drift space

optDrift = SRWLOptD(0.1) #Drift space
optApert = SRWLOptA('r', 'a', slitDX*1e-6, slitDY*1e-6) #Aperture

#optDrift = SRWLOptD(12.9) #Drift space

# propagParApert = [0, 0, 1., 0, 0, 1.5, 1.0, 1.1, 8., 0, 0, 0]
propagParApert = [0, 0, 1., 0, 0, 1.5, 1.0, 1.1, 6., 0, 0, 0]
propagParLens =  [0, 0, 1., 0, 0, 1.0, 1.0, 1.0, 1., 0, 0, 0]
propagParDrift = [0, 0, 1., 1, 0, 1.0, 1.2, 1.0, 1., 0, 0, 0]

#Wavefront Propagation Parameters:
#[0]: Auto-Resize (1) or not (0) Before propagation
#[1]: Auto-Resize (1) or not (0) After propagation
#[2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
#[3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
#[4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
#[5]: Horizontal Range modification factor at Resizing (1. means no modification)
#[6]: Horizontal Resolution modification factor at Resizing
#[7]: Vertical Range modification factor at Resizing
#[8]: Vertical Resolution modification factor at Resizing
#[9]: Type of wavefront Shift before Resizing (not yet implemented)
#[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
#[11]: New Vertical wavefront Center position after Shift (not yet implemented)
#optBL = SRWLOptC([optApert, optCRL, optDrift], [propagParApert, propagParLens, propagParDrift]) #"Beamline" - Container of Optical Elements (together with the 
optBL = SRWLOptC([optApert, optDrift], [propagParApert, propagParDrift]) #"Beamline" - Container of Optical Elements (together with the corresponding wavefront propagation instructions)
#optBL = SRWLOptC([optApert, optLens, optDrift], [propagParApert, propagParLens, propagParDrift]) #"Beamline" - Container of Optical Elements (together with the corresponding wavefront propagation instructions)


print('Simulating Partially-Coherent Wavefront Propagation by summing-up contributions of SR from individual electrons (takes time)... ')
#nMacroElec = 50000 #total number of macro-electrons
nMacroElec = 10000 #total number of macro-electrons
nMacroElecAvgPerProc = 5 # 5 #number of macro-electrons / wavefront to average on worker processes before sending data to master (for parallel calculation only)
nMacroElecSavePer = 10 #intermediate data saving periodicity (in macro-electrons)
srCalcMeth = 1 #SR calculation method
srCalcPrec = 0.01 #SR calculation rel. accuracy
radStokesProp = srwl_wfr_emit_prop_multi_e(elecBeam, magFldCnt, meshInitPartCoh, srCalcMeth, srCalcPrec, nMacroElec, nMacroElecAvgPerProc, nMacroElecSavePer, os.path.join(os.getcwd(), strExDataFolderName, strIntOutFileNamePartCoh), sampFactNxNyForProp, optBL)
print('done')

plotMeshX = [1000*meshInitPartCoh.xStart, 1000*meshInitPartCoh.xFin, meshInitPartCoh.nx]
plotMeshY = [1000*meshInitPartCoh.yStart, 1000*meshInitPartCoh.yFin, meshInitPartCoh.ny]
#uti_plot2d(radStokesProp.arS, plotMeshX, plotMeshY, ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Power Density'])


#uti_plot_show() #show all graphs (blocks script execution; close all graph windows to proceed)
print('done')

 
