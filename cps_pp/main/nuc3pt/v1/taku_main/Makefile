#all: cpsMPI.x sse-wilson_dslash.s
all: cpsMPI.x

pull:
	git pull origin

push:
	git push origin


test: cpsMPI.x
	sh ./test.sh

SRCDIR=~meifeng/code/CPS508/
OBJDIR=~meifeng/code/CPS508/objs-debug

clean:
	/bin/rm -f *.o *.x

CC=~meifeng/code/CPS508/objs-debug/mpicc
CXX=~meifeng/code/CPS508/objs-debug/mpic++


QMPD=~izubuchi/Work/CPS508/objs-debug/qmp-2.3.1
QIOD=~izubuchi/Work/CPS508/objs-debug/qmp-2.3.1


DFLAGS =  -DUSE_QMP -DUSE_QIO

INCLUDES =  -I $(OBJDIR)/ -I $(SRCDIR)/include/ \
  -I/home/izubuchi/lib/qmp-2.3.1/include    \
  -I/home/izubuchi/lib/qio-2.3.8/include \
  $(DFLAGS)

# this is 55.85 % 
ICCOPT=  -openmp -g -O2 -xS  -fno-alias  -align -inline-max-total-size=10000 -ip  -i-static $(INCLUDES)

GCCOPT= -fopenmp -g -O3


INTEL_LIBS=

HEADERS =    sse-blk-dag1.h  sse-macros-dag0.h sse-macros-dag1.h \
sse-bnd-dag0.h sse-bnd-dag1.h  sse-blk-dag0.h  sse-inlines.h  sse-subs.h test-order.h 


main.o: main.C
	$(CXX) $(ICCOPT)  -c  $<

#dirc objects
DOBJDIR= $(SRCDIR)/main/sse_f_dwf_inv/
OBJS=  $(DOBJDIR)/sse-wilson_dslash.o  $(DOBJDIR)/wilson_init.o $(DOBJDIR)/wilson_end.o \
        $(DOBJDIR)/sse-dwf_dslash_4.o $(DOBJDIR)/dwf_dslash.o  $(DOBJDIR)/dwf_dslash_5_plus.o 

cpsMPI.x: main.o
	$(CXX) $(ICCOPT) main.o  \
        $(OBJDIR)/cps-e.a  $(INTEL_LIBS)  \
	$(OBJS) \
        /home/izubuchi/lib/qmp-2.3.1/lib/libqmp.a \
        /home/izubuchi/lib/qio-2.3.8/lib/libqio.a \
        /home/izubuchi/lib/qio-2.3.8/lib/liblime.a \
        -o cpsMPI.x
