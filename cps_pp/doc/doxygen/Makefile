#include<config.h>
CPS_START_NAMESPACE

all: html/index.html

html/index.html: doxygen.cfg
	doxygen doxygen.cfg

clean:
	rm -rf html/*.html
CPS_END_NAMESPACE
