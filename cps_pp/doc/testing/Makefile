#include<config.h>
CPS_START_NAMESPACE
PDFLATEX=pdflatex
LATEX = latex
DVIPS = dvips

SRC =  testing

all:	$(SRC).pdf $(SRC).html $(SRC).ps $(SRC).tar.gz

$(SRC).pdf: clean $(SRC).tex
	$(PDFLATEX) $(SRC)
	$(PDFLATEX) $(SRC)
	$(PDFLATEX) $(SRC)

$(SRC).ps: $(SRC).dvi
	$(DVIPS) -o $(SRC).ps $(SRC).dvi

$(SRC).dvi: clean $(SRC).tex
	$(LATEX) $(SRC)
	$(LATEX) $(SRC)
	$(LATEX) $(SRC)

$(SRC).html: $(SRC).pdf
	tth -e2 -a -w1 -i -V $(SRC)

$(SRC).tar.gz: $(SRC).tex
	tar -cvf $(SRC).tar *
	gzip $(SRC).tar

clean:
	rm -rf "#"*"#" *~ core *.toc *.tar.gz *.lof *.log *.aux *.dvi *.ps
CPS_END_NAMESPACE
