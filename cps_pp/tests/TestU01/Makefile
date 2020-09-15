#Will have to be edited manually 

TESTU01=/home/chulwoo/CPS/TestU01
CPS=/home/chulwoo/CPS/branches/c11_RNG/

SRCDIR = /home/chulwoo/CPS/branches/c11_RNG/cps_pp
BUILDDIR = /home/chulwoo/CPS/branches/c11_RNG/C11

INCLUDE_FLAGS = -I${BUILDDIR} -I${SRCDIR}/include  -I/home/chulwoo/SciDAC//scidac/install/qmp-CPS-OPENMPI/include    -I/home/chulwoo/SciDAC/scidac/install/qio-2.3.9-OPENMPI/include

LDFLAGS =  -lgmp -L/home/chulwoo/SciDAC//scidac/install/qmp-CPS-OPENMPI/lib -lqmp -L/home/chulwoo/SciDAC//scidac/install/qmp-CPS-OPENMPI/lib -lqmp -L/home/chulwoo/SciDAC/scidac/install/qio-2.3.9-OPENMPI/lib -lqio -llime -L${TESTU01}/lib  -lmylib -ltestu01

DFLAGS +=  -DUSE_SSE -DGMP -DUSE_QMP -DUSE_QIO  -DUSE_OMP -DSSE_TO_C -DC11 -DC11_RNG
.SUFFIXES: .c .C .o .x

clean:
	rm -f *.o *.x
.c.o:	
	mpicc -I ${TESTU01}/include -c $< 

.C.o:
	mpic++ -c -I ${TESTU01}/include ${INCLUDE_FLAGS} ${DFLAGS}  $< 

cps.x:	u_cps.o cps.o
	mpic++ $^ ${CPS}/C11/cps.a ${LDFLAGS} -o $@
