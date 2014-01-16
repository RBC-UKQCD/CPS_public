TOP=.
CONFIG=$(TOP)/config/ssed.conf

CC       = @echo "Compile: <(CC)>"
XCC      = "<COMPILER FOR MAIN and LIBS>"
CFLAGS   = "<(CFLAGS)>"
COPT     = "<(COPT)>"
CINCDIRS = "<(CINCDIRS)>"
OBJECT   = "<object name>"
SOURCE   = "<source name>"

include $(CONFIG)

sources = dwf.nw interface.nw

code.c = dwf.c \
         dwf-ssef.c \
         dwf-ssed.c \
         dwf-bluelightf.c \
         dwf-bluelightd.c \
         dwf-altivecf.c
code.h = dwf-tables.h

interface.h = dwf-ssef.h dwf-ssed.h \
              dwf-bluelightd.h dwf-bluelightf.h \
              dwf-altivecf.h

TESTS= qdpc xqmp bluelight ssef
subdirs = $(TESTS:%=tests/%)

.PHONY: all clean realclean tar pdf programs compile

f.tex=$(sources:%.nw=%.tex)
f.dvi=$(sources:%.nw=%.dvi)
f.pdf=$(sources:%.nw=%.pdf)

lib:
	ar rcs libdwf-sse.a dwf-sse*o

all: pdf programs

clean:
	$(RM) *.tex *.log *.dvi *.aux *.toc *.nwi *.out $(code.c:%.c=%.o) *~
	for d in $(subdirs); do make -C $$d $@; done

realclean: clean
	$(RM) *.pdf $(code.c) $(interface.h) $(code.h)
	for d in $(subdirs); do make -C $$d $@; done

tar:
	base=`basename $$PWD` ; \
	cd .. ; \
	tar -cvf - $$base | bzip2 -9 > $$base.tar.bz2

pdf: $(f.pdf)

programs: $(interface.h) $(code.c) $(code.h)

compile: $(OBJECT)

$(OBJECT): $(SOURCE) $(SOURCE:%.c=%.h) dwf.c
	$(CC) $(CFLAGS) $(COPT) $(CINCDIRS) -c -o $(OBJECT) $(SOURCE)

$(f.pdf): %.pdf: %.dvi
	dvipdfm $<

$(f.dvi): %.dvi: %.tex
	latex $<
	latex $<
	latex $<

$(f.tex): %.tex: %.nw
	noweave -index -latex -delay $< > $@

$(interface.h): %.h: interface.nw
	notangle $(TFLAGS) -R$@ -c $< > $@ \
	|| (rm -f $@ ; exit 1)

$(code.c): %.c: dwf.nw
	notangle $(TFLAGS) -R$@ -c $< > $@ \
	|| (rm -f $@ ; exit 1)

$(code.h): %.h: dwf.nw
	notangle $(TFLAGS) -R$@ -c $< > $@ \
	|| (rm -f $@ ; exit 1)
