#!/bin/bash -e

# environment variable QOS should be set as the top of QCDOC OS tree
#INSTALL_DIR=/home/chulwoo/QCDOC/CPS/QCDOC
#./configure --build sparc-sun-solaris2.9 --host powerpc-gnu-elf --enable-target=qcdoc --enable-installdir=${INSTALL_DIR}

INSTALL_DIR=/home/chulwoo/QCDOC/CPS/QCDOCHOST0/
./configure --enable-optimise --enable-installdir=${INSTALL_DIR}

mkdir -p ${INSTALL_DIR}/src
cp -p config.h Makefile.rules Makefile.vars Makefile.compile Makefile tests/Makefile_common* ${INSTALL_DIR}
touch ${INSTALL_DIR}/src/Makefile_depend
