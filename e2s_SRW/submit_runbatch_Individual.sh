#
# memo-script to run SRWLIB_individualelectrons.py in BATCH mode 
# on the cluster. 
#
# submit_runbatch_Individual.sh uses openmpi with a defined number of nodes
# to start a batch session using runbatch_Individual.sh
# MA 13/8/2018
#
#qsub -q ap-medium.q -l redhat_release=rhel6 runbatch_Individual.sh
qsub -q ap-high.q -l redhat_release=rhel6 -V -pe openmpi 50 runbatch_Individual.sh

