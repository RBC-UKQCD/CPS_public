include ../../../../Makefile.users

#Choose Lanczos and Dslash library
DFLAGS += -DUSE_BFM_A2A -DUSE_BFM_LANCZOS
#DFLAGS += -DUSE_GRID_A2A -DUSE_GRID_LANCZOS

#Choose A2A precision (use double unless you run out of memory in the A2A contraction part)
DFLAGS += -DA2A_PREC_DOUBLE
#DFLAGS += -DA2A_PREC_SINGLE

#Choose BLAS implementation
#BLAS=-lgslcblas 
BLAS=-lcblas -latlas

#Optionally include EIGEN for Lanczos improvement
INCLIST += -Iinclude/eigen3  #Change to your eigen include directory
DFLAGS += -DUSE_EIGEN


LDFLAGS +=  -lgsl ${BLAS} -lrt

