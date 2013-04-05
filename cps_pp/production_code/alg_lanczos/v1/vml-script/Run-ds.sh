#! /bin/bash
#PBS -S /bin/bash
#PBS -N cps-lanczos
#PBS -j oe
#PBS -o ./o${PBS_JOBID}
#PBS -m n
#PBS -q test_ds -A dwf-ps -l nodes=2,walltime=00:03:00

#ECHO=echo

cd ${PBS_O_WORKDIR}

mkdir -p eigen

# print identifying info for this job
echo "Job ${PBS_JOBNAME} submitted from ${PBS_O_HOST} started "`date`" jobid ${PBS_JOBID}"
       
MPI_DIR=/usr/local/mvapich-1.2rc1/

# determine number of cores per host
coresPerNode=`cat /proc/cpuinfo | grep -c processor`
# count the number of nodes listed in PBS_NODEFILE
nNodes=$[`cat ${PBS_NODEFILE} | wc --lines`]
(( nNodes= 0 + nNodes ))
(( nCores = nNodes * coresPerNode ))
echo "NODEFILE nNodes=$nNodes ($nCores cores):"
cat ${PBS_NODEFILE}
       
# Always use fcp or rsync to stage any large input files from the cluster head node
# to your job's control worker node.
# All worker nodes have attached disk storage in /scratch
# Copy below is commented since there are no files to transfer
#fcp -c rsh -p ds1.fnal.gov:$MY_DIR/*.vml /scratch
       
#The directory the job was submitted from is $PBS_O_WORKDIR.
echo ${PBS_O_WORKDIR}
       


# LMA_SHIFTs,  set 0 if you don't want to do LMA or AMA
#
# For example, on 16^3 x 32 lattice,  and LMA_SHIFTS="2 2 2 4" then  source points are
# (0,0,0,0), (8,0,0,0) ... (0,0,0,8), ... (8,8,8,24), 
#  in total 2x2x2x4=32 source points
#
LMA_SHIFTS="0 0 0 0"

# do_cg  : wanna do full-CG ?  0 .. No
do_cg=0

# do_eigv_read :   wanna do reading eigenvector files into cache  before do anything ?  0 .. No
do_eigv_read=0


# geometry of nodes
geom="2 2 4 4"


#####################################################################
#work-directory LMA_SHIFT[0] LMA_SHIFT[1] LMA_SHIFT[2] LMA_SHIFT[3] do_cg do_eigv_read

# -hostfile $PBS_NODEFILE \
#$ECHO $MPI_DIR/bin/mpirun -np  1 \


# Is controlling NUMA effective ?
source /usr/local/mvapich/etc/mvapich.conf.sh

NUMA_WRAP=/usr/local/mvapich/bin/numa_32_mv

$ECHO $MPI_DIR/bin/mpirun -np  $nCores \
$NUMA_WRAP  ./cps-lanczos.x  ./ $LMA_SHIFTS $do_cg $do_eigv_read \
-qmp-geom $geom  \
> ./out.${PBS_JOBID}


# Always use fcp to copy any large result files you want to keep back
# to the head node before exiting your script. The /scratch area on the
# workers is wiped clean between jobs.
#
# Copy below is commented since there are no files to transfer
#fcp -c rsh -p /scratch/vacpol* ${PBS_O_WORKDIR}/.
       
exit





