AC_INIT(./src/)
dnl AC_CONFIG_HEADERS([conf.h])
AC_CHECK_HEADERS([malloc.h])

dnl Get host and target information:
AC_CANONICAL_TARGET


#CK: makedepend is soooooo sllooooooooww. Disabling in favour of the cmd line version
#AC_CHECK_PROG(HAVE_MAKEDEPEND,makedepend,yes,no)
HAVE_MAKEDEPEND="no"
AC_SUBST(HAVE_MAKEDEPEND)
dnl Are we compiling with QMP library?
AC_ARG_ENABLE(qmp, 
	AC_HELP_STRING([--enable-qmp], 
		[Compile with SciDAC Communication (QMP) lib. (default=no)]), 
	,[enable_qmp="no"])
AC_MSG_CHECKING([if including QMP])
AC_MSG_RESULT($enable_qmp)

dnl Set correct compiler information
if test "$cross_compiling" = "yes"; then
   if test "$host" != "none-none-none"; then
      if test "${host_alias}" = "powerpc64-bgq-linux"; then
        AC_PATH_PROGS(CC, mpicc)
        AC_PATH_PROGS(CXX, mpicxx)
      else
        AC_PATH_PROGS(CC, ${host_alias}-gcc)
        AC_PATH_PROGS(CXX, ${host_alias}-g++)
      fi
      AC_PATH_PROGS(LD, ${host_alias}-ld)
      AC_PATH_PROGS(AR, ${host_alias}-ar)
      AC_PATH_PROGS(AS, ${host_alias}-as)
      if test -z "$CXX" || test -z "$LD" || test -z "$AR" || test -z "$AS"; then
        echo "Missing required cross-compilation tools"
        exit 1
      fi
   fi
else
dnl check if Intel compiler is available if running on 686
   if test "$host" = "i686-pc-linux-gnu" && \
   test "$build" = "i686-pc-linux-gnu"; then
	   CC=""
       AC_CHECK_PROGS(CC, icc, gcc)
       AC_CHECK_PROGS(CXX, icc, g++)
   elif test "$host" = "powerpc-ibm-aix5.2.0.0" ||\
   test "$host" = "powerpc64-unknown-linux-gnu" ;then
	   CC=""
       AC_CHECK_PROGS(CC, xlC, gcc)
       AC_CHECK_PROGS(CXX, xlC, g++)
   elif test "$enable_qmp" != "no"; then
   CC=`$enable_qmp/bin/qmp-config --cc`
   AC_CHECK_PROGS(CXX, mpicxx mpic++ mpiCC, g++)
     AC_PATH_PROGS(CXX, mpicxx mpic++ mpiCC, g++)
   else
     AC_PROG_CXX
   fi
   AC_CHECK_PROGS(LD, ld)
   AC_CHECK_PROGS(AR, ar)
   AC_CHECK_PROGS(AS, as)
   if test -z "$CXX" || test -z "$LD" || test -z "$AR" || test -z "$AS"; then
     echo "Missing required compilation tools"
     exit 1
   fi
fi
dnl path for the compiler
CC_PATH=`which "$CXX" |sed -e 's/bin\/powerpc-gnu-elf-g++//'`
AC_SUBST(CC_PATH)

dnl Check for maths library
dnl AC_CHECK_LIB(m, sin, ,echo "CPS needs the math library - bailing out"; exit 1)

dnl Set the relevant include files
topwd_srcdir=`pwd`
AC_SUBST(topwd_srcdir)

dnl This is the _default_ default
namespace_default="yes" 

dnl What target are we compiling for?
AC_ARG_ENABLE(target, 
	AC_HELP_STRING([--enable-target], 
		[What target are we compiling for? (default=noarch)]), 
	,[enable_target="noarch"])
AC_MSG_CHECKING([what target we are compiling for])
AC_MSG_RESULT($enable_target)

dnl Can add more cases as we need to
if test "$enable_target" = "bgq"; then
   ARCH="BGQ"
   ARCHDIR="BGQ_DIR"
   TESTING_PARALLEL="yes";
   DFLAGS="${DFLAGS} -DUSE_OMP -DVEC_INLINE"
   LDFLAGS="-fopenmp ${LDFLAGS}"
   CFLAGS="${CFLAGS} -fopenmp"
   CXXFLAGS="${CXXFLAGS} -fopenmp"

#CK: Add an option to use bfm for noarch machines
elif test "$enable_target" = "noarch-bfm"; then
   ARCH=NOARCH
   ARCHDIR="NOARCH_BFM_DIR"
   if test "$enable_qmp" != "no"; then
   ARCHDIR="NOARCH_BFM_QMP_DIR"
   CC=`$enable_qmp/bin/qmp-config --cc`
   AC_CHECK_PROGS(CXX, mpic++ mpiCC mpCC, g++)
   TESTING_PARALLEL="yes";
   else
   TESTING_PARALLEL="no";
   fi

   #openmp assumed
   DFLAGS="${DFLAGS} -DUSE_OMP -DVEC_INLINE"

else

   ARCH=NOARCH
   ARCHDIR="NOARCH_DIR"
   if test "$enable_qmp" != "no"; then
   CC=`$enable_qmp/bin/qmp-config --cc`
   AC_CHECK_PROGS(CXX, mpic++ mpiCC mpCC, g++)
   TESTING_PARALLEL="yes";
   else
   TESTING_PARALLEL="no";
   fi
fi


dnl Need to put in machine/compiler specific optimisation flags here
AC_ARG_ENABLE(openmp, 
	AC_HELP_STRING([--enable-openmp], 
	       [Compile with OpenMP. (default=yes)]),
	,[enable_openmp="yes"])
AC_MSG_CHECKING([if using OpenMP])
AC_MSG_RESULT($enable_openmp)
if test "$enable_openmp" = "yes"; then
   DFLAGS="${DFLAGS} -DUSE_OMP"
fi

dnl Need to put in machine/compiler specific optimisation flags here
AC_ARG_ENABLE(optimise, 
	AC_HELP_STRING([--enable-optimise], 
	       [Compile with maximum optimisation. (default=yes)]),
	,[enable_optimise="yes"])
AC_MSG_CHECKING([if optimising])
AC_MSG_RESULT($enable_optimise)
if test "$enable_optimise" = "yes"; then
   if test "$enable_target" = "bgl"; then
     CFLAGS="${CFLAGS} -O -Wall"
     CXXFLAGS="${CXXFLAGS} -O -Wall"
   elif test "$enable_target" = "bgq"; then
     CFLAGS="${CFLAGS} -g -O2 -Wall -fno-strict-aliasing "
     CXXFLAGS="${CXXFLAGS} -g -O2 -Wall -fno-strict-aliasing "
   else
     CFLAGS="${CFLAGS} -g -O3 -Wall -fno-strict-aliasing "
     CXXFLAGS="${CXXFLAGS} -g -O3 -Wall -fno-strict-aliasing "
   fi
fi

dnl Are we compiling with SSE?
AC_ARG_ENABLE(sse, 
	AC_HELP_STRING([--enable-sse], 
	[Compile with SSE. (default=no)]), 
	,[enable_sse="no"])
AC_MSG_CHECKING([if including SSE])
AC_MSG_RESULT($enable_sse)
if test "$enable_target" != "noarch"; then
   enable_sse="no"
fi
if test "$enable_qmp" == "no"; then
   echo "QMP needed to compile SSE dslash"
   enable_sse="no"
fi
if test "$enable_sse" != "no"; then
   CFLAGS="${CFLAGS} -msse -msse2 -msse3"
   CXXFLAGS="${CXXFLAGS} -msse -msse2 -msse3"
   DFLAGS="${DFLAGS} -DUSE_SSE"
fi

dnl Are we compiling with C11?
AC_ARG_ENABLE(c11, 
	AC_HELP_STRING([--enable-c11], 
	[Compile with C11 features. (default=yes)]), 
	,[enable_c11="yes"])
AC_MSG_CHECKING([if compile with c++11 ])
AC_MSG_RESULT($enable_c11)
if test "$enable_c11" != "no"; then
   CXXFLAGS="${CXXFLAGS} -std=c++11"
   DFLAGS="${DFLAGS} -DUSE_C11 " 
fi

dnl Are we compiling with C11 RNG?
AC_ARG_ENABLE(c11_rng, 
	AC_HELP_STRING([--enable-c11-rng=sitmo|mt19937|ranlux48], 
	[Compile with C11 RNG. (default=yes, sitmo if enabled)]), 
	,[enable_c11_rng="no"])
	AC_MSG_CHECKING([if compile with c++11 RNG])
	AC_MSG_RESULT($enable_c11_rng)
if test "$enable_c11_rng" == "mt19937"; then
   DFLAGS="${DFLAGS} -DUSE_C11_RNG -DUSE_C11_MT " 
elif test "$enable_c11_rng" == "ranlux48"; then
   DFLAGS="${DFLAGS} -DUSE_C11_RNG -DUSE_C11_RANLUX " 
elif test "$enable_c11_rng" != "no"; then
   DFLAGS="${DFLAGS} -DUSE_C11_RNG -DUSE_C11_SITMO " 
fi


dnl -Og flag should overide any previous optimisation declarations
AC_ARG_ENABLE(debug, 
	AC_HELP_STRING([--enable-debug], 
		[Compile with debug information. (default=no)]),  
	,[enable_debug="no"])
AC_MSG_CHECKING([if debugging])
AC_MSG_RESULT($enable_debug)
if test "$enable_debug" = "yes"; then
   CFLAGS="${CFLAGS} -g -Og"
   CXXFLAGS="${CXXFLAGS} -g -Og"
fi

dnl Are we compiling with cg-dwf?
AC_ARG_ENABLE(cg_dwf, 
	AC_HELP_STRING([--enable-cg-dwf], 
		[Compile with cg-dwf (default=no)]), 
	,[enable_cg_dwf="no"])
AC_MSG_CHECKING([if including cg-dwf])
AC_MSG_RESULT($enable_cg_dwf)
if test "$enable_cg_dwf" != "no"; then
	if test "$enable_cg_dwf" != "yes"; then
   		LDFLAGS="-L${enable_cg_dwf} -ldwf-sse ${LDFLAGS}"
   		INCLUDE_FLAGS="${INCLUDE_FLAGS} -I${enable_cg_dwf}"
	else
   		LDFLAGS="-L\${SRCDIR}/aux/cg-dwf-1.3.3 -ldwf-sse ${LDFLAGS}"
   		INCLUDE_FLAGS="${INCLUDE_FLAGS} -I\${SRCDIR}/aux/cg-dwf-1.3.3"
	fi
   	DFLAGS="${DFLAGS} -DUSE_CG_DWF"
fi

dnl Are we compiling with BLAS?
AC_ARG_ENABLE(blas, 
	AC_HELP_STRING([--enable-blas], 
		[Compile with BLAS (default=no)]), 
	,[enable_blas="no"])
AC_MSG_CHECKING([if including blas])
AC_MSG_RESULT($enable_blas)
if test "$enable_blas" != "no"; then
   	DFLAGS="${DFLAGS} -DUSE_BLAS"
fi

dnl Are we compiling with LAPACK?
AC_ARG_ENABLE(lapack, 
	AC_HELP_STRING([--enable-lapack], 
		[Compile with lapack (default=no)]), 
	,[enable_lapack="no"])
AC_MSG_CHECKING([if including lapack])
AC_MSG_RESULT($enable_lapack)
if test "$enable_lapack" != "no"; then
   	DFLAGS="${DFLAGS} -DUSE_BLAS -DUSE_LAPACK"
fi


dnl Are we compiling with gnu multiprecision library?
AC_ARG_ENABLE(gmp, 
	AC_HELP_STRING([--enable-gmp], 
		[Compile with GNU Multiprecision lib. (default=no)]), 
	,[enable_gmp="no"])
AC_MSG_CHECKING([if including GMP])
AC_MSG_RESULT($enable_gmp)
AC_ARG_ENABLE(mpfr, 
	AC_HELP_STRING([--enable-mpfr], 
		[Compile with GMP mpfr routines. (default=no)]), 
	,[enable_mpfr="no"])
AC_MSG_CHECKING([if including GMP mpfr])
AC_MSG_RESULT($enable_mpfr)
if test "$enable_gmp" != "no"; then
   DFLAGS="${DFLAGS} -DUSE_GMP"
   LDFLAGS="-lgmp ${LDFLAGS}"
  if test "$enable_gmp" != "yes"; then
   LDFLAGS="-L${enable_gmp}/lib  ${LDFLAGS}"
   INCLUDE_FLAGS="${INCLUDE_FLAGS} -I${enable_gmp}/include"
  fi
  if test "$enable_mpfr" != "no"; then
     LDFLAGS="-lmpfr ${LDFLAGS}"
     DFLAGS="${DFLAGS} -DUSE_MPFR"
    if test "$enable_mpfr" != "yes"; then
     LDFLAGS="-L${enable_mpfr}/lib ${LDFLAGS} "
     INCLUDE_FLAGS="${INCLUDE_FLAGS} -I${enable_mpfr}/include"
    fi
  fi
fi

dnl Are we compiling with QIO library?
AC_ARG_ENABLE(qio, 
	AC_HELP_STRING([--enable-qio], 
		[Compile with SciDAC I/O (QIO) lib. (default=no)]), 
	,[enable_qio="no"])
AC_MSG_CHECKING([if including QIO])
AC_MSG_RESULT($enable_qio)

if test "$enable_qmp" != "no"; then
   if test "$enable_qmp" != "yes"; then
     QMP_DIR=${enable_qmp}
   elif test "X${QMP_DIR}X" = "XX"; then
     echo "QMP_DIR environment variable not set"
     exit 1
   fi
   QMP_CFLAGS=`${QMP_DIR}/bin/qmp-config --cflags`
   QMP_LDFLAGS=`${QMP_DIR}/bin/qmp-config --ldflags`
   QMP_LIBS=`${QMP_DIR}/bin/qmp-config --libs`
   LDFLAGS="${QMP_LDFLAGS} ${QMP_LIBS} ${LDFLAGS}"
   INCLUDE_FLAGS="${INCLUDE_FLAGS} ${QMP_CFLAGS} "
   DFLAGS="${DFLAGS} -DUSE_QMP"
   TESTING_PARALLEL="yes";
fi

dnl Are we compiling with MDWF library?
AC_ARG_ENABLE(mdwf, 
	AC_HELP_STRING([--enable-mdwf], 
		[Compile with Mobius DWF lib. (default=no)]), 
	,[enable_mdwf="no"])
AC_MSG_CHECKING([if including MDWF])
AC_MSG_RESULT($enable_mdwf)
if test "$enable_mdwf" != "no"; then
   if test "X${MDWF_DIR}X" = "XX"; then
     MDWF_DIR=${enable_mdwf}
   fi
   if test "X${QMP_DIR}X" = "XX"; then
     echo "QMP_DIR environment variable not set"
     exit 1
   fi
   MDWF_CFLAGS=`${MDWF_DIR}/bin/mdwf-config --cflags`
   MDWF_LDFLAGS=`${MDWF_DIR}/bin/mdwf-config --ldflags`
   MDWF_LIBS=`${MDWF_DIR}/bin/mdwf-config --libs`
   LDFLAGS="${MDWF_LDFLAGS} ${MDWF_LIBS} ${LDFLAGS}"
   INCLUDE_FLAGS="${INCLUDE_FLAGS} ${MDWF_CFLAGS}"
   DFLAGS="${DFLAGS} -DUSE_MDWF"
fi

if test "$enable_qio" != "no"; then
   if test "X${QMP_DIR}X" = "XX"; then
     echo "QMP directory not set"
     exit 1
   fi
   QMP_CFLAGS=`${QMP_DIR}/bin/qmp-config --cflags`
   QMP_LDFLAGS=`${QMP_DIR}/bin/qmp-config --ldflags`
   QMP_LIBS=`${QMP_DIR}/bin/qmp-config --libs`
   LDFLAGS="${QMP_LDFLAGS} ${QMP_LIBS} ${LDFLAGS}"
   INCLUDE_FLAGS="${INCLUDE_FLAGS} ${QMP_FLAGS} "
if test "$enable_qio" != "yes"; then
   QIO_DIR=${enable_qio}
   LDFLAGS="-L${QIO_DIR}/lib ${LDFLAGS}"
   INCLUDE_FLAGS="${INCLUDE_FLAGS} -I${QIO_DIR}/include"
fi
   LDFLAGS="${LDFLAGS} -lqio -llime "
   DFLAGS="${DFLAGS} -DUSE_QIO"
fi

AC_ARG_ENABLE(xml, 
	AC_HELP_STRING([--enable-xml], 
		[Compile with XML lib. (default=no)]), 
	,[enable_xml="no"])
AC_MSG_CHECKING([if including XML])
AC_MSG_RESULT($enable_xml)

if test "$enable_xml" != "no"; then
   if test "$enable_xml" != "yes"; then
     XML_DIR=${enable_xml}
#   elif test "X${XML_DIR}X" = "XX"; then
#     echo "XML_DIR environment variable not set"
#     exit 1

   XML_CFLAGS=`${XML_DIR}/bin/xml2-config --cflags`
   XML_LDFLAGS=`${XML_DIR}/bin/xml2-config --libs`

   LDFLAGS="${XML_LDFLAGS} ${LDFLAGS}"
   INCLUDE_FLAGS="${INCLUDE_FLAGS} ${XML_CFLAGS}"
   else 
   LDFLAGS="${XML_LDFLAGS} ${LDFLAGS}"
   INCLUDE_FLAGS="${INCLUDE_FLAGS} ${XML_CFLAGS}"
   fi
fi

dnl Are we compiling with BFM library?
AC_ARG_ENABLE(bfm, 
	AC_HELP_STRING([--enable-bfm], 
		[Compile with BAGEL fermion sparse-matrix (BFM) lib. (default=no)]), 
	,[enable_bfm="no"])
AC_MSG_CHECKING([if including BFM])
AC_MSG_RESULT($enable_bfm)

if test "$enable_bfm" != "no"; then
   if test "X${BFM_DIR}X" = "XX"; then
dnl     echo "BFM_DIR environment variable not set"
     BFM_DIR=${enable_bfm}
   fi


   if test "$enable_bfm" != "yes"; then
      LDFLAGS="-L${BFM_DIR}/lib ${LDFLAGS} "
      INCLUDE_FLAGS="-I${BFM_DIR}/include ${INCLUDE_FLAGS} "
   else
      BFM_DIR=${QMP_DIR}
      LDFLAGS="-L${QMP_DIR}/lib ${LDFLAGS} "
      INCLUDE_FLAGS="-I${QMP_DIR}/include ${INCLUDE_FLAGS} "
   fi
      DFLAGS="${DFLAGS} -DUSE_BFM -DUSE_BFM_MINV -DUSE_BFM_TM -DUSE_HDCG"

   #CK: don't assume QIO
   if test "$enable_qio" != "no"; then
	LDFLAGS="-lqio ${LDFLAGS}"
   fi		    					 

dnl      LDFLAGS="${LDFLAGS} -lqdp -lXPathReader -lxmlWriter -lxml2 -lqio -lbfm -lgsl -lgslcblas -lqmp"
dnl -L/soft/libraries/essl/current/lib64 -lesslbg -L/soft/compilers/ibmcmp-nov2012/xlf/bg/14.1/lib64 -lxlf90_r -lxlfmath -lxl -ldl 
   if test "X{$enable_xml}X" = "XX"; then
    enable_xml="yes"
   fi
fi

if test "$enable_xml" != "no"; then
   if test "$enable_xml" != "yes"; then
     XML_DIR=${enable_xml}
   elif test "X${XML_DIR}X" = "XX"; then
     XML_DIR=${BFM_DIR}
   fi

   if test "${XML_DIR}" = "no"; then
     echo "XML_DIR environment variable not set"
     exit 1
   fi

   XML_CFLAGS=`${XML_DIR}/bin/xml2-config --cflags`
   XML_LDFLAGS=`${XML_DIR}/bin/xml2-config --libs`

   LDFLAGS="${XML_LDFLAGS} ${LDFLAGS}"
   INCLUDE_FLAGS="${INCLUDE_FLAGS} ${XML_CFLAGS}"
fi

if test "$TESTING_PARALLEL" = "yes"; then
   ENABLE_PARALLEL="define PARALLEL 1"
else
   ENABLE_PARALLEL="undef PARALLEL"
fi
dnl -------------------------windy---------------------------
dnl  Edited by hjkim. This configuration is need by QUDA&CUDA
dnl -------------------------windy---------------------------
dnl Are we compiling with QUDA GPU library?
dnl Check whether --enable-cuda was given.
AC_ARG_ENABLE(cuda, 
	AC_HELP_STRING([--enable-cuda], 
		[Location of CUDA installation. (default=no)]), 
	,[enable_cuda="no"])
AC_MSG_CHECKING([if including CUDA])
AC_MSG_RESULT($enable_cuda)

dnl Check whether --enable-quda was given.
AC_ARG_ENABLE(quda, 
	AC_HELP_STRING([--enable-quda], 
		[Compile with QUDA GPU lib. (default=no)]), 
	,[enable_quda="no"])
AC_MSG_CHECKING([if including QUDA])
AC_MSG_RESULT($enable_quda)

dnl QUDA was given first and CUDA was given later, 
dnl this is because library dependency.
if test "$enable_quda" != "no"; then
   if test "X${QUDA_DIR}X" = "XX"; then
dnl     echo "QUDA_DIR environment variable not set"
     QUDA_DIR=${enable_quda}
   fi
   CXXFLAGS="${CXXFLAGS} -std=c++14"
   INCLUDE_FLAGS="${INCLUDE_FLAGS} -I${QUDA_DIR}/include"
   LDFLAGS="-L${QUDA_DIR}/lib -lquda ${LDFLAGS}"
   DFLAGS="${DFLAGS} -DUSE_QUDA"
fi

if test "$enable_cuda" != "no"; then
   if test "X${CUDA_DIR}X" = "XX"; then
dnl     echo "CUDA_DIR environment variable not set"
     CUDA_DIR=${enable_cuda}
   fi
   INCLUDE_FLAGS="${INCLUDE_FLAGS} -I${CUDA_DIR}/include"
   LDFLAGS="${LDFLAGS} -L${CUDA_DIR}/lib64 -lcudart -lcuda -lcublas "
fi

dnl If --enable-cuda was given for --enable-quda case, program will stop.
if test "$enable_quda" != "no"; then
   if test "X${CUDA_DIR}X" = "XX"; then
     echo "CUDA_DIR environment variable not set (required for QUDA)"
     exit 1
   fi
fi
#-------------------------windy---------------------------	

dnl Are we compiling with fftw?

AC_SEARCH_LIBS([fftw_execute],[fftw3],
	[AC_SEARCH_LIBS([fftwf_execute],[fftw3f], [], [AC_MSG_ERROR("single precision fftw library not found")])]
	[AC_DEFINE([USE_FFTW],[1],[Define 1 if you have fftw3 library])]
	[enable_fftw=true])

#AC_SEARCH_LIBS([crc32],[z],
#	[AC_DEFINE([USE_ZLIB],[1],[Define 1 if you have zlib])]
#	[enable_zlib=true],
#	[AC_MSG_ERROR("zlib needed for crc32")]
#	)

#AC_ARG_ENABLE(fftw, 
#	AC_HELP_STRING([--enable-fftw], 
#		[Compile with FFTW (default=no)]), 
#	,[enable_fftw="no"])
#AC_MSG_CHECKING([if including fftw])
#AC_MSG_RESULT($enable_fftw)
#if test "$enable_fftw" != "no"; then
#	if test "$enable_fftw" != "yes"; then
#   		LDFLAGS="${LDFLAGS} -L${enable_fftw}/lib"
#   		INCLUDE_FLAGS="${INCLUDE_FLAGS} -I${enable_fftw}/include"
##	else
#   		LDFLAGS="${LDFLAGS} -L\${SRCDIR}/aux/GotoBLAS2 -lgoto2"
#   		INCLUDE_FLAGS="${INCLUDE_FLAGS} -I\${SRCDIR}/aux/GotoBLAS2"
#	fi
#   	DFLAGS="${DFLAGS} -DUSE_FFTW"
#fi

if test "$enable_bfm" != "no"; then
   HAVE_BFM="define HAVE_BFM 1"
else
   HAVE_BFM="undef HAVE_BFM"
fi

dnl Are we compiling with grid?
AC_ARG_ENABLE(grid, 
	AC_HELP_STRING([--enable-grid], 
		[Compile with Grid library (default=no)]), 
	,[enable_grid="no"])
AC_MSG_CHECKING([if including Grid])
AC_MSG_RESULT($enable_grid)
if test "$enable_grid" != "no"; then
	if test "$enable_grid" != "yes"; then
   		LDFLAGS="-L${enable_grid}/build/Grid -lGrid ${LDFLAGS}"
   		INCLUDE_FLAGS="${INCLUDE_FLAGS} -I${enable_grid}/ -I${enable_grid}/build/Grid "
	else
   		LDFLAGS="-lGrid ${LDFLAGS}"
	fi
   	CXXFLAGS="${CXXFLAGS} -std=c++11"
   	DFLAGS="${DFLAGS} -DUSE_GRID"
fi
	
dnl Shall we use the cps namespace?
AC_ARG_ENABLE(namespace, 
	AC_HELP_STRING([--enable-namespace],
		[Build the library inside a 'cps' namespace. (default=no for QCDOC, yes otherwise)]),
	,[enable_namespace=$namespace_default])
AC_MSG_CHECKING([whether we are using the namespace])
AC_MSG_RESULT($enable_namespace)
	
if test "$enable_namespace" = "yes"
then
        CPS_END_NAMESPACE="}" 
        USING_NAMESPACE_CPS="using namespace cps;"
        CPS_NAMESPACE="cps"
        CPS_START_NAMESPACE="namespace cps {"
fi

dnl Make visible all variables declared above
AC_SUBST(AR)
AC_SUBST(ARCH)
AC_SUBST(ARCHDIR)
AC_SUBST(AS)
AC_SUBST(ASFLAGS)
AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(DFLAGS)
AC_SUBST(INCLUDE_FLAGS)
AC_SUBST(LD)
AC_SUBST(LDFLAGS)
AC_SUBST(LIBS)
AC_SUBST(TESTING_PARALLEL)
AC_SUBST(ENABLE_PARALLEL)
AC_SUBST(QOS)
AC_SUBST(DO_GMP)
AC_SUBST(CPS_END_NAMESPACE)
AC_SUBST(CPS_START_NAMESPACE)
AC_SUBST(USING_NAMESPACE_CPS)
AC_SUBST(CPS_NAMESPACE)



dnl Single/double precision switch:
AC_ARG_ENABLE(double-prec, 
	AC_HELP_STRING([--enable-double-prec],
		[Build the library in double precision (default=yes)]),
	,[enable_double_precision="yes"])
AC_MSG_CHECKING([if double-prec is enabled])
AC_MSG_RESULT($enable_double_precision)
if test "$enable_double_precision" = "yes"; then
  GSUMPRECISION="double"
  LOCALPRECISION="double"
  PRECDIR="DOUBLE_DIR"
else
  GSUMPRECISION="double"
  LOCALPRECISION="float"
  PRECDIR="FLOAT_DIR"
fi
AC_SUBST(GSUMPRECISION)
AC_SUBST(LOCALPRECISION)
AC_SUBST(PRECDIR)


if test "X${SRCDIR}X" = "XX"; then
	SRCDIR='${SRCDIR0}'
fi
if test "X${BUILDDIR}X" = "XX"; then
	BUILDDIR='${BUILDDIR0}'
fi
AC_SUBST(SRCDIR)
AC_SUBST(BUILDDIR)

dnl configuring version number
echo "`git log -n 1 --format=format:"#define GITHASH \\"%H:%d\\"%n" HEAD`" > $srcdir/include/version.h
#git log -n 1 --format=format:GITHASH=\\\"%H:%d\\\"%n\\\" HEAD > $srcdir/GITHASH
.  $srcdir/VERSION
#.  $srcdir/GITHASH
echo CPS source: $srcdir
echo "CPS Version MAJOR: ${VERSION_MAJOR}"
echo "CPS Version MINOR: ${VERSION_MINOR}"
echo "CPS Version SUB:   ${VERSION_SUB}"
echo Version ${VERSION_STR}
echo "SRCDIR = " ${SRCDIR}
echo "BUILDDIR = " ${BUILDDIR}
AC_SUBST(VERSION_MAJOR)
AC_SUBST(VERSION_MINOR)
AC_SUBST(VERSION_SUB)
AC_SUBST(VERSION_STR)
AC_SUBST(GITHASH)
AC_SUBST(SRCDIR)
AC_SUBST(BUILDDIR)

AC_SUBST(HAVE_BFM)

AC_PROG_INSTALL

AC_CHECK_PROG(have_doxygen, doxygen, "yes")
AC_CHECK_TYPES(bool_t)
AC_CHECK_TYPES(enum_t)
AC_CHECK_TYPES([int8_t],,,[#include<stdint.h>])
AC_CHECK_TYPES(int16_t)
AC_CHECK_TYPES(int32_t)
AC_CHECK_TYPES(int64_t)
AC_CHECK_TYPES(quad_t)
AC_CHECK_TYPES(uint_t)
AC_CHECK_TYPES(uint8_t)
AC_CHECK_TYPES(uint16_t)
AC_CHECK_TYPES(uint32_t)
AC_CHECK_TYPES(uint64_t)
AC_CHECK_TYPES(u_quad_t)
AC_CHECK_TYPES(netobj)
#AC_CHECK_FUNCS(sync)
AC_CHECK_FUNCS(posix_memalign)
dnl AC_CHECK_HEADERS([qcdocos/scu_checksum.h])
#AC_CHECK_HEADERS([gmp.h])
AC_CHECK_HEADERS([gmp.h])
AC_CHECK_HEADERS([cblas.h])
AC_CHECK_HEADERS([mkl_cblas.h])
AC_CHECK_HEADERS([gsl_cblas.h])
AC_CONFIG_HEADERS([conf.h])


AC_OUTPUT(config.h \
	  precision.h \
	  Makefile \
	  Makefile.rules \
	  Makefile.users \
          tests/Makefile_common \
          tests/regression.pl \
	  include/input/Makefile\
	  )	  
