#!/bin/bash
#This script takes a list of configuration numbers and adds them into the
#appropriate logfiles_t_src_?/logs_configs directory if not already there
#so that the user.qcsh script will now include these in its "to do list"
#of configurations.
#
#
#Arguments:
#
#1. filename of the text file containing the list of configurations numbers
#2. a configuration number for which you the shift of the lattice to be 0
#   (and will go in the logfiles_t_src_0/logs_configs directory)
#3. the spacing between configuration numbers that are actually used

if [ $# -ne 3 ]; then
	echo "Usage: add_configs configlist_filename config_src0 delta_config"
	exit 1
fi
listfile=$1
config_src0=$2
delta_config=$3
for config_no in $( cat $listfile ); do
    let displ=config_no-config_src0
    let tmp=displ%delta_config
    if [ $tmp -ne 0 ]; then
	echo "Error: Invalid config number $config_no"
	continue
    fi
    let tsrc=displ/delta_config
    let tsrc=tsrc%4
    let tsrc=tsrc+4
    let tsrc=tsrc%4
    let tsrc=tsrc*16
    touch logfiles_t_src_$tsrc/logs_configs/$config_no
done
