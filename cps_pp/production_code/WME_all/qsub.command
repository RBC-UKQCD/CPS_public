#!/bin/tcsh                                                                                                                                                  
#setenv TZ GMT+24                                                                                                                                            
set binary = ../binaries/BGQ.x.v1
#set binary = ../binaries/BGQ.x.v1
qsub -A LatticeQCD -q default -n $1 -t $2 --env BG_THREADLAYOUT=1:BG_PROCESSESPERNODE=1:BG_SHAREDMEMSIZE=32:PAMI_DEVICE=M:BG_ALLOW_CACHELINE_LOCKING=1:BG_L2LOCK_L1_ONLY=1:L1P_POLICY=std:HPM_GROUP=0:HPM_SCOPE=node:PAMID_CONTEXT_MAX=1:PAMID_EAGER=65535 --mode c1 $binary `pwd`  -qmp-geom native
