#!/bin/sh 
#------ qsub option --------#
#MJS: -prj G10012
#MJS: -pc
#MJS: -proc 2
#MJS: -thread 1
#MJS: -compiler intel
####MJS: -parallel openmpi
#MJS: -mem 20mb
#MJS: -time 10:00:00
#MJS: -eo
#MJS: -cwd
#------- FTL command -------#
#FTL_STAT: detail
#FTLDIR: $MJS_CWD
#FTL_NO_RANK: on
#BEFORE_R: 0: vml_files, configurations
#AFTER: 0: *
####AFTER: res
#------- Program execution -------#


JOBID=propgen

qmp_geom="-qmp-geom 2 1 1 1"

#ROOT="/home/izubuchi/Work/CPS508/main/prop_gen/v1/testrun"
ROOT="."
exec_file=$ROOT/cpsMPI.x

exec_dir="$ROOT"

configs_temp_dir_qrun="$ROOT/configurations"
config_no=500
t_src=0
 
ensemble_label=IWASAKI_Nf2p1_32c64s64_b1.75AuxDet_mu0.001ms0.045
ensemble_id=141
#logs_configs_dir_qrun=$ROOT/logs_t_src_${t_src}
logs_configs_dir_qrun=$ROOT
plaq_gfix_info_dir=$logs_configs_dir_qrun

midprop_contractions_dir=$ROOT/midprop_contractions
configs_done_log_qrun=$logs_configs_dir_qrun/donelog

mkdir -p $logs_configs_dir_qrun
mkdir -p $plaq_gfix_info_dir
mkdir -p $midprop_contractions_dir
#mkdir -p $configs_done_log_qrun

# remembering the finished mass
touch $logs_configs_dir_qrun/${config_no}.mdone

subm= mpirun \
$exec_file 0 $exec_dir vml_files/do_arg.in dummy vml_files/mass_list.in \
$configs_temp_dir_qrun/ckpoint_lat.$config_no none 17 $t_src \
$ensemble_label $ensemble_id $config_no $logs_configs_dir_qrun \
$plaq_gfix_info_dir $midprop_contractions_dir $configs_done_log_qrun \
$qmp_geom 

echo $subm

$subm



