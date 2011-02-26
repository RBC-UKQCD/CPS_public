#!/bin/sh 
#------ qsub option --------#
#MJS: -prj G10011
#MJS: -pc
#MJS: -proc 64
#MJS: -thread 1
#MJS: -compiler intel
####MJS: -parallel openmpi
#MJS: -mem 1024mb
#MJS: -time 24:00:00
#MJS: -eo
#MJS: -cwd
#------- FTL command -------#
#FTL_STAT: detail
#FTLDIR: $MJS_CWD
#FTL_NO_RANK: on
#BEFORE: *.vml, *.x
#BEFORE_R: configurations
#AFTER:0: * 
####AFTER: res
#------- Program execution -------#


JOBID=propgen

qmp_geom="-qmp-geom 2 2 4 4"

ROOT="."
exec_file=$ROOT/cpsMPI.x
exec_dir=$ROOT

subm= mpirun \
$exec_file do_arg.vml meas_arg.vml $exec_dir 0 0 2 \
$qmp_geom 

echo $subm

$subm



