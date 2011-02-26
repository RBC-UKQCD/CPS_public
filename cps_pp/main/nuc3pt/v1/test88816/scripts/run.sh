#!/bin/sh
#--- qsub option ---#
#MJS: -pc

#MJS: -proc 4
#MJS: -thread 1
#MJS: -time 10:00:00
#MJS: -mem 20mb

#MJS: -compiler intel
#MJS: -eo
#MJS: -cwd
#MJS: -rerun Y

#--- FTL command ---#
#FTL_STAT: detail
#FTL_NO_RANK: on

#BEFORE: cpsMPI.x, do_arg.vml, meas_arg.vml
#BEFORE_R: configurations
 
#AFTER_R: work/, results/
#FTL_SUFFIX: on

qmp_geom="-qmp-geom 1 1 1 1"
exec_file=./cpsMPI.x
DIR=.

cmd=mpirun $exec_file do_arg.vml meas_arg.vml $DIR 0 0 0 




