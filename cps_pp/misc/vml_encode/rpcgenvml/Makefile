CC= gcc
CXX= g++
CSRC:= $(wildcard *.c)
CPPSRC:= 
COBJ=$(CSRC:.c=.o)
CPPOBJ=$(CPPSRC:.C=.o)

.SUFFIXES: .C .c .o

.c.o: $(CSRC)
	$(CC) $(CFLAGS) -c $< -o $@
.C.o: $(CPPSRC)
	$(CXX) $(CFLAGS) -c $< -o $@

$(COBJ): $(CSRC) 

$(CPPOBJ): $(CPPSRC)

all: rpcgenvml

rpcgenvml: $(COBJ)
	$(CC) $(COBJ) -o $@

clean: 
	- rm -f $(COBJ) rpcgenvml *~
