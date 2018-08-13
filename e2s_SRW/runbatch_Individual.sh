#
# memo-script to run SRWLIB_individualelectrons.py in BATCH mode 
# on the cluster. 
#
# it makes use of mpiexec as discussed with M.Furseman, and it is launched 
# by submit_runbatch_Individual.sh 
# MA 13/8/2018
#
##/dls_sw/apps/python/anaconda/1.7.0/64/bin/python SRW_individualelectrons.py SRW.input
mpiexec /dls_sw/apps/python/anaconda/1.7.0/64/bin/python SRW_individualelectrons.py SRW.input

