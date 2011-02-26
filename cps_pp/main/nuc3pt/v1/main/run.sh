#!/bin/sh
#--- qsub option ---#
#MJS: -pc
#MJS: -prj G10011
#MJS: -proc 4
#MJS: -thread 1
#MJS: -time 1:00:00
#MJS: -mem 240mb
#MJS: -compiler intel
#MJS: -eo
#MJS: -cwd
#MJS: -rerun Y

#--- FTL command ---#
#FTL_STAT: detail

#FTL_NO_RANK: on

#BEFORE: *.x, *.vml, run.sh
#BEFORE_R: configurations

#AFTER_R: work, results, configurations
#FTL_SUFFIX: on

exec_file=./cpsMPI.x
DIR=.

cmd=mpirun $exec_file do_arg.vml meas_arg.vml $DIR 0 0 2 -qmp-geom 1 1 1 1 








