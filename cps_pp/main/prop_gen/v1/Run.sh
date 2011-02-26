#!/bin/sh 

JOBID=propgen
#num of OpenMP threads
n_omp=2
#num of MPI processes
n_mpi=1

qmp_geom="-qmp-geom 1 1 1 1"


exec_file="/home/izubuchi/Work/CPS508/main/prop_gen/v1/cpsMPI.x"

ROOT="/home/izubuchi/Work/CPS508/main/prop_gen/v1/testrun"

exec_dir="$ROOT"

configs_temp_dir_qrun="$ROOT/configurations"
config_no=500
t_src=0

max_cg_iter=10

ensemble_label=IWASAKI_Nf2p1_32c64s64_b1.75AuxDet_mu0.001ms0.045
ensemble_id=141
logs_configs_dir_qrun=$ROOT/logs_t_src_${t_src}
plaq_gfix_info_dir=$logs_configs_dir_qrun

midprop_contractions_dir=$ROOT/midprop_contractions
configs_done_log_qrun=$logs_configs_dir_qrun/donelog

mkdir -p $logs_configs_dir_qrun
mkdir -p $plaq_gfix_info_dir
mkdir -p $midprop_contractions_dir
#mkdir -p $configs_done_log_qrun

# remembering the finished mass
touch $logs_configs_dir_qrun/$config_no

export  OMP_NUM_THREADS=$n_omp
subm= mpirun -np $n_mpi -x OMP_NUM_THREADS \
$exec_file 0 $exec_dir vml_files/do_arg.in dummy vml_files/mass_list.in \
$configs_temp_dir_qrun/ckpoint_lat.$config_no none 17 $t_src $max_cg_iter \
$ensemble_label $ensemble_id $config_no $logs_configs_dir_qrun \
$plaq_gfix_info_dir $midprop_contractions_dir $configs_done_log_qrun \
$qmp_geom 

echo $subm

$subm



