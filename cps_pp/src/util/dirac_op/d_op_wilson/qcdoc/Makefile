NORECURSE_DIR = qcdoc qcdsp noarch mpi CVS 
QCDOC_DIR = qcdoc
QCDOC_SERIAL_DIR = qcdoc
QCDSP_DIR = qcdsp
MPI_DIR = mpi
NOARCH_DIR = noarch
include ../../../../../Makefile.vars
#include ../../../../../Makefile.rules

ABI_FLAGS = -ffixed-%r14 \
	-ffixed-%r15 \
	-ffixed-%r16 \
        -ffixed-%r17 \
        -ffixed-%r18 \
        -ffixed-%r19 \
        -ffixed-%r20 \
        -ffixed-%r21 \
        -ffixed-%r22 \
        -ffixed-%r23 \
        -ffixed-%r24 \
        -ffixed-%r25 \
        -ffixed-%r26 \
        -ffixed-%r27 \
        -ffixed-%r28 \
        -ffixed-%r29 \
        -ffixed-%r30 \
        -ffixed-%fr14 \
        -ffixed-%fr15 \
        -ffixed-%fr16 \
        -ffixed-%fr17 \
        -ffixed-%fr18 \
        -ffixed-%fr19 \
        -ffixed-%fr20 \
        -ffixed-%fr21 \
        -ffixed-%fr22 \
        -ffixed-%fr23 \
        -ffixed-%fr24 \
        -ffixed-%fr25 \
        -ffixed-%fr26 \
        -ffixed-%fr27 \
        -ffixed-%fr28 \
        -ffixed-%fr29 \
        -ffixed-%fr30 \
        -ffixed-%fr31 

CFLAGS:= $(CFLAGS) $(ABI_FLAGS)
CXXFLAGS:= $(CXXFLAGS) $(ABI_FLAGS)

include ../../../../../Makefile.compile
