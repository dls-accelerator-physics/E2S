import pandas as pd
import os

def CrPar(Ephot, CrystName):

    crpar = []
    print ">>>>>>>>> inside interpol <<<<<<<<<< {}".format(os.getcwd())
    print " Ephot =  {}".format(Ephot)

    filin = '/dls/physics/xph53246/source_to_beamline/E2S/e2s_BLOPTICS/'+CrystName+'.csv'
    df1 = pd.read_csv(filin).set_index('E(keV)')
    print filin
    before, after = df1.loc[:Ephot].iloc[-1], df1.loc[Ephot:].iloc[0]
    sandwich = pd.DataFrame([before, pd.Series(name=Ephot), after])
    s = sandwich.interpolate()
    crpar.append(s.loc[Ephot].psi0rSi111)
    crpar.append(s.loc[Ephot].psi0iSi111)
    crpar.append(s.loc[Ephot].psihrSi111)
    crpar.append(s.loc[Ephot].psihiSi111)
    crpar.append(s.loc[Ephot].psihbrSi111)
    crpar.append(s.loc[Ephot].psihbiSi111)


    print " crpar {}".format(crpar)
    return crpar
