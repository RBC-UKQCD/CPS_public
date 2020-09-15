/********************************************************************************
*                                                                              *
*    LICENSED MATERIALS - PROPERTY OF IBM                                      *
*                                                                              *
*    "RESTRICTED MATERIALS OF IBM"                                             *
*                                                                              *
*    5765-042                                                                  *
*    5765-C42                                                                  *
*    5765-F82                                                                  *
*    5765-G17                                                                  *
*    5765-H25                                                         *
*                                                                              *
*    (C) COPYRIGHT IBM CORP. 1991, 2010.  All rights reserved.                 *
*                                                                              *
*    US Government Users Restricted Rights - Use, duplication or disclosure    *
*    restricted by GSA ADP Schedule Contract with IBM Corp.                    *
*                                                                              *
********************************************************************************
*                                                                              *
*  Program name - <essl.h> header file                                         *
*  Descriptive name - ESSL C,C++ language header file                          *
*                                                                              *
*  Function : This file must be included in any C or C++ file                  *
*             containing ESSL calls in order for the ESSL calls                *
*             to work as documented in the ESSL Guide and Reference.           *
*                                                                              *
*  Change activity -                                                           *
*                    ESSL Version 3.1.0.0                     12/1997          *
*                    ESSL Version 3.1.0.1(APAR PQ16381)       06/1998          *
*                    ESSL Version 3.1.1.0                     10/1998          *
*                    ESSL Version 3.1.2.0                     10/1999          *
*                    ESSL Version 3.2.0.0                     07/2000          *
*                    ESSL Version 3.1.2.1(APAR PQ42486)       10/2000          *
*                    ESSL Version 3.2.0.1(APAR PQ47660)       08/2001          *
*                    ESSL Version 3.3.0.0                     12/2001          *
*                    ESSL Version 4.1.0.0                     07/2003          *
*                    ESSL Version 4.2.0.0                     10/2004          *
*                    ESSL Version 4.2.2.0                     11/2005          *
*                    ESSL Version 4.3.0.0                     08/2007          *
*                    ESSL Version 4.4.0.0                     11/2008          *
*                    ISO header file fix.                     11/2008          *
*                    ESSL Version 5.1.0.0                     02/2010          *
*******************************************************************************/
#ifdef _ESV6464
    #define _ESVINT long
#else
    #define _ESVINT int
#endif
 
#ifndef __cplusplus
  /* Begin C Path */
 
  #ifndef __essl
   #define __essl 1
 
 
   /* Definition of complex data types */
 
   #ifndef  _CMPLX
    #define  _CMPLX 1
    #ifndef  _REIM
      #define _REIM   1
    #endif
    typedef union { struct { float  _re, _im; } _data; double _align; }  cmplx;
   #endif
 
   #ifndef  _DCMPLX
    #define  _DCMPLX 1
    #ifndef  _REIM
      #define _REIM   1
    #endif
    typedef union { struct { double _re, _im; } _data; double _align; } dcmplx;
   #endif
 
   #ifdef  _REIM
    #define RE(x)   ((x)._data._re)
    #define IM(x)   ((x)._data._im)
   #endif
 
   /* Complex arguments defined differently for C and C++ */
 
   #define _ESVCM  cmplx
   #define _ESVCOM dcmplx
 
 
   /* Scalar output arguments defined differently for C and C++ */
 
   #define _ESVI    _ESVINT *
   #define _ESVS    float *
   #define _ESVD    double *
   #define _ESVC    _ESVCM *
   #define _ESVZ    _ESVCOM *
 
  #endif        /* End __essl */
 
  /* End C Path */
 
#else
  /* Begin C++ Path */
 
  #ifndef __essl
   #define __essl 1
   /* Complex arguments defined differently for C and C++ */
 
   #ifdef __blrts
       #ifndef __linux
           /* Force BlueGene/L through the same path as other Linuxes */
           #define __linux 22
       #endif
   #endif
 
   #ifndef __linux
 
   /* Check if user had included new iostream; if so, we don't allow old complex
      to be included. */
 
   #if !defined(_COMPLEX_) && !defined(_ESV_COMPLEX_) && !defined(_IOSTREAM_)
     #include <complex.h>
     #define _ESVCM cmplx
     #define _ESVCOM complex
   #else
     #include <complex>
     using std::complex;
     #define _ESVCM complex<float>
     #define _ESVCOM complex<double>
   #endif
   #else
     #include <complex>
     using std::complex;
     #define _ESVCM complex<float>
     #define _ESVCOM complex<double>
   #endif
 
   /* Scalar output arguments defined differently for C and C++ */
 
   #ifndef _ESVCPTR
     #define _ESVI    _ESVINT &
     #define _ESVS    float &
     #define _ESVD    double &
     #define _ESVC    _ESVCM &
     #define _ESVZ    _ESVCOM &
   #else
     #define _ESVI    _ESVINT *
     #define _ESVS    float *
     #define _ESVD    double *
     #define _ESVC    _ESVCM *
     #define _ESVZ    _ESVCOM *
   #endif
 
   /* Definition of complex data types */
 
   #ifndef __linux
   #ifdef COMPLEXH
     #ifndef  _CMPLX
       #define  _CMPLX 1
       class cmplx
         {
          private:
            float _re,_im;
          public:
             cmplx() { _re = 0.0; _im = 0.0; }
             cmplx(float r, float i = 0.0) { _re = r; _im = i; }
             friend inline float sreal(const cmplx& a) { return a._re; }
             friend inline float simag(const cmplx& a) { return a._im; }
         };
     #endif
   #endif
   #endif
 
  #endif        /* End __essl */
 
  /* End C++ Path */
#endif
 
#ifdef __cplusplus
  extern "C" {
#endif
 
/*  Vector-Scalar Linear Algebra Subprograms  */
 
#define isamax    esvisamax
#define idamax    esvidamax
#define icamax    esvicamax
#define izamax    esvizamax
 
#define isamin    esvisamin
#define idamin    esvidamin
 
#define ismax     esvismax
#define idmax     esvidmax
 
#define ismin     esvismin
#define idmin     esvidmin
 
#define sasum     esvsasum
#define dasum     esvdasum
#define scasum    esvscasum
#define dzasum    esvdzasum
 
#define saxpy     esvsaxpy
#define daxpy     esvdaxpy
#define caxpy     esvcaxpy
#define zaxpy     esvzaxpy
 
#define scopy     esvscopy
#define dcopy     esvdcopy
#define ccopy     esvccopy
#define zcopy     esvzcopy
 
#define sdot      esvsdot
#define ddot      esvddot
 
#define cdotu     esvcdotu
#define zdotu     esvzdotu
#define cdotc     esvcdotc
#define zdotc     esvzdotc
 
#define snaxpy    esvsnaxpy
#define dnaxpy    esvdnaxpy
 
#define sndot     esvsndot
#define dndot     esvdndot
 
#define snrm2     esvsnrm2
#define dnrm2     esvdnrm2
#define scnrm2    esvscnrm2
#define dznrm2    esvdznrm2
 
#define snorm2    esvsnorm2
#define dnorm2    esvdnorm2
#define cnorm2    esvcnorm2
#define znorm2    esvznorm2
 
#define srotg     esvsrotg
#define drotg     esvdrotg
#define crotg     esvcrotg
#define zrotg     esvzrotg
 
#define srot      esvsrot
#define drot      esvdrot
#define crot      esvcrot
#define zrot      esvzrot
#define csrot     esvcsrot
#define zdrot     esvzdrot
 
#define sscal     esvsscal
#define dscal     esvdscal
#define cscal     esvcscal
#define zscal     esvzscal
#define csscal    esvcsscal
#define zdscal    esvzdscal
 
#define sswap     esvsswap
#define dswap     esvdswap
#define cswap     esvcswap
#define zswap     esvzswap
 
#define svea      esvsvea
#define dvea      esvdvea
#define cvea      esvcvea
#define zvea      esvzvea
 
#define sves      esvsves
#define dves      esvdves
#define cves      esvcves
#define zves      esvzves
 
#define svem      esvsvem
#define dvem      esvdvem
#define cvem      esvcvem
#define zvem      esvzvem
 
#define syax      esvsyax
#define dyax      esvdyax
#define cyax      esvcyax
#define zyax      esvzyax
#define csyax     esvcsyax
#define zdyax     esvzdyax
 
#define szaxpy    esvszaxpy
#define dzaxpy    esvdzaxpy
#define czaxpy    esvczaxpy
#define zzaxpy    esvzzaxpy
 
/*  Sparse Vector-Scalar Linear Algebra Subprograms  */
 
#define ssctr      esvssctr
#define dsctr      esvdsctr
#define csctr      esvcsctr
#define zsctr      esvzsctr
 
#define sgthr      esvsgthr
#define dgthr      esvdgthr
#define cgthr      esvcgthr
#define zgthr      esvzgthr
 
#define sgthrz     esvsgthrz
#define dgthrz     esvdgthrz
#define cgthrz     esvcgthrz
#define zgthrz     esvzgthrz
 
#define saxpyi     esvsaxpyi
#define daxpyi     esvdaxpyi
#define caxpyi     esvcaxpyi
#define zaxpyi     esvzaxpyi
 
#define sdoti      esvsdoti
#define ddoti      esvddoti
#define cdotui     esvcdotui
#define zdotui     esvzdotui
#define cdotci     esvcdotci
#define zdotci     esvzdotci
 
/*  Dense Matrix-Vector Subprograms  */
 
#define sgemv      esvsgemv
#define dgemv      esvdgemv
#define cgemv      esvcgemv
#define zgemv      esvzgemv
 
#define sgemx      esvsgemx
#define dgemx      esvdgemx
 
#define sgemtx     esvsgemtx
#define dgemtx     esvdgemtx
 
#define sger       esvsger1
#define dger       esvdger1
#define sger1      esvsger1
#define dger1      esvdger1
#define cgeru      esvcgeru
#define zgeru      esvzgeru
#define cgerc      esvcgerc
#define zgerc      esvzgerc
 
#define sspmv      esvsspmv
#define dspmv      esvdspmv
#define chpmv      esvchpmv
#define zhpmv      esvzhpmv
 
#define ssymv      esvssymv
#define dsymv      esvdsymv
#define chemv      esvchemv
#define zhemv      esvzhemv
 
#define sslmx      esvsslmx
#define dslmx      esvdslmx
 
#define sspr       esvsspr
#define dspr       esvdspr
#define chpr       esvchpr
#define zhpr       esvzhpr
 
#define ssyr       esvssyr
#define dsyr       esvdsyr
#define cher       esvcher
#define zher       esvzher
 
#define sslr1      esvsslr1
#define dslr1      esvdslr1
 
#define sspr2      esvsspr2
#define dspr2      esvdspr2
#define chpr2      esvchpr2
#define zhpr2      esvzhpr2
 
#define ssyr2      esvssyr2
#define dsyr2      esvdsyr2
#define cher2      esvcher2
#define zher2      esvzher2
 
#define sslr2      esvsslr2
#define dslr2      esvdslr2
 
#define sgbmv      esvsgbmv
#define dgbmv      esvdgbmv
#define cgbmv      esvcgbmv
#define zgbmv      esvzgbmv
 
#define ssbmv      esvssbmv
#define dsbmv      esvdsbmv
#define chbmv      esvchbmv
#define zhbmv      esvzhbmv
 
#define strmv      esvstrmv
#define dtrmv      esvdtrmv
#define ctrmv      esvctrmv
#define ztrmv      esvztrmv
 
#define stpmv      esvstpmv
#define dtpmv      esvdtpmv
#define ctpmv      esvctpmv
#define ztpmv      esvztpmv
 
#define stbmv      esvstbmv
#define dtbmv      esvdtbmv
#define ctbmv      esvctbmv
#define ztbmv      esvztbmv
 
/*  Sparse Matrix-Vector Linear Algebra Subprograms  */
 
#define dsmmx      esvdsmmx
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define dsmtm      esvdsmtm_er
#else
  #define dsmtm      esvdsmtm
#endif
 
#define dsdmx      esvdsdmx
 
/*  Matrix Operation Subroutines  */
 
#define sgeadd     esvsgeadd
#define dgeadd     esvdgeadd
#define cgeadd     esvcgeadd
#define zgeadd     esvzgeadd
 
#define sgesub     esvsgesub
#define dgesub     esvdgesub
#define cgesub     esvcgesub
#define zgesub     esvzgesub
 
#define sgemul     esvsgemul
#define dgemul     esvdgemul
#define cgemul     esvcgemul
#define zgemul     esvzgemul
 
#define dgemlp     esvdgemlp
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define sgemms     esvsgemms_er
  #define dgemms     esvdgemms_er
  #define cgemms     esvcgemms_er
  #define zgemms     esvzgemms_er
#else
  #define sgemms     esvsgemms
  #define dgemms     esvdgemms
  #define cgemms     esvcgemms
  #define zgemms     esvzgemms
#endif
 
#define sgemm      esvsgemm
#define dgemm      esvdgemm
#define cgemm      esvcgemm
#define zgemm      esvzgemm
 
#define ssymm      esvssymm
#define dsymm      esvdsymm
#define csymm      esvcsymm
#define zsymm      esvzsymm
#define chemm      esvchemm
#define zhemm      esvzhemm
 
#define strmm      esvstrmm
#define dtrmm      esvdtrmm
#define ctrmm      esvctrmm
#define ztrmm      esvztrmm
 
#define ssyrk      esvssyrk
#define dsyrk      esvdsyrk
#define csyrk      esvcsyrk
#define zsyrk      esvzsyrk
#define cherk      esvcherk
#define zherk      esvzherk
 
#define ssyr2k     esvssyr2k
#define dsyr2k     esvdsyr2k
#define csyr2k     esvcsyr2k
#define zsyr2k     esvzsyr2k
#define cher2k     esvcher2k
#define zher2k     esvzher2k
 
#define sgetmi     esvsgetmi
#define dgetmi     esvdgetmi
#define cgetmi     esvcgetmi
#define zgetmi     esvzgetmi
 
#define sgetmo     esvsgetmo
#define dgetmo     esvdgetmo
#define cgetmo     esvcgetmo
#define zgetmo     esvzgetmo
 
/*  Dense Linear Algebraic Equation Subroutines  */
 
#define sgef       esvsgef
#define dgef       esvdgef
#define cgef       esvcgef
#define zgef       esvzgef
 
#define slange     esvslange
#define dlange     esvdlange
#define clange     esvclange
#define zlange     esvzlange
 
#define sgecon     esvsgecon
#define dgecon     esvdgecon
#define cgecon     esvcgecon
#define zgecon     esvzgecon
 
#define sgetrf     esvsgetrf
#define dgetrf     esvdgetrf
#define cgetrf     esvcgetrf
#define zgetrf     esvzgetrf
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define dgefp      esvdgefp_er
#else
  #define dgefp      esvdgefp
#endif
 
#define sges       esvsges
#define dges       esvdges
#define cges       esvcges
#define zges       esvzges
 
#define sgesm      esvsgesm
#define dgesm      esvdgesm
#define cgesm      esvcgesm
#define zgesm      esvzgesm
 
#define sgetrs     esvsgetrs
#define dgetrs     esvdgetrs
#define cgetrs     esvcgetrs
#define zgetrs     esvzgetrs
 
#define sgesv      esvsgesv
#define dgesv      esvdgesv
#define cgesv      esvcgesv
#define zgesv      esvzgesv
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define sgefcd     esvsgefcd_er
  #define dgefcd     esvdgefcd_er
#else
  #define sgefcd     esvsgefcd
  #define dgefcd     esvdgefcd
#endif
 
#define sppf       esvsppf
#define dppf       esvdppf
 
#define spptrf     esvspptrf
#define dpptrf     esvdpptrf
#define cpptrf     esvcpptrf
#define zpptrf     esvzpptrf
 
#define slansp     esvslansp
#define dlansp     esvdlansp
#define clanhp     esvclanhp
#define zlanhp     esvzlanhp
 
#define sppcon     esvsppcon
#define dppcon     esvdppcon
#define cppcon     esvcppcon
#define zppcon     esvzppcon
 
#define spptri     esvspptri
#define dpptri     esvdpptri
#define cpptri     esvcpptri
#define zpptri     esvzpptri
 
#define spof       esvspof
#define dpof       esvdpof
#define cpof       esvcpof
#define zpof       esvzpof
 
#define spotrf     esvspotrf
#define dpotrf     esvdpotrf
#define cpotrf     esvcpotrf
#define zpotrf     esvzpotrf
 
#define slansy     esvslansy
#define dlansy     esvdlansy
#define clanhe     esvclanhe
#define zlanhe     esvzlanhe
 
#define spocon     esvspocon
#define dpocon     esvdpocon
#define cpocon     esvcpocon
#define zpocon     esvzpocon
 
#if defined(__ESVERR) || defined(_ESVERR)
#define dppfp      esvdppfp_er
#else
#define dppfp      esvdppfp
#endif
 
#define spps       esvspps
#define dpps       esvdpps
 
#define spptrs     esvspptrs
#define dpptrs     esvdpptrs
#define cpptrs     esvcpptrs
#define zpptrs     esvzpptrs
 
#define sppsv      esvsppsv
#define dppsv      esvdppsv
#define cppsv      esvcppsv
#define zppsv      esvzppsv
 
#define sposm      esvsposm
#define dposm      esvdposm
#define cposm      esvcposm
#define zposm      esvzposm
 
#define spotrs     esvspotrs
#define dpotrs     esvdpotrs
#define cpotrs     esvcpotrs
#define zpotrs     esvzpotrs
 
#define sposv      esvsposv
#define dposv      esvdposv
#define cposv      esvcposv
#define zposv      esvzposv
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define sppfcd     esvsppfcd_er
  #define dppfcd     esvdppfcd_er
  #define spofcd     esvspofcd_er
  #define dpofcd     esvdpofcd_er
#else
  #define sppfcd     esvsppfcd
  #define dppfcd     esvdppfcd
  #define spofcd     esvspofcd
  #define dpofcd     esvdpofcd
#endif
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define dbssv      esvdbssv_er
  #define dbstrf     esvdbstrf_er
#else
  #define dbssv      esvdbssv
  #define dbstrf     esvdbstrf
#endif
#define dbstrs     esvdbstrs
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define sgeicd     esvsgeicd_er
  #define dgeicd     esvdgeicd_er
#else
  #define sgeicd     esvsgeicd
  #define dgeicd     esvdgeicd
#endif
 
#define sgetri     esvsgetri
#define dgetri     esvdgetri
#define cgetri     esvcgetri
#define zgetri     esvzgetri
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define sppicd     esvsppicd_er
  #define dppicd     esvdppicd_er
 
  #define spoicd     esvspoicd_er
  #define dpoicd     esvdpoicd_er
#else
  #define sppicd     esvsppicd
  #define dppicd     esvdppicd
 
  #define spoicd     esvspoicd
  #define dpoicd     esvdpoicd
#endif
 
#define spotri     esvspotri
#define dpotri     esvdpotri
#define cpotri     esvcpotri
#define zpotri     esvzpotri
 
#define strsv      esvstrsv
#define dtrsv      esvdtrsv
#define ctrsv      esvctrsv
#define ztrsv      esvztrsv
 
#define stpsv      esvstpsv
#define dtpsv      esvdtpsv
#define ctpsv      esvctpsv
#define ztpsv      esvztpsv
 
#define strsm      esvstrsm
#define dtrsm      esvdtrsm
#define ctrsm      esvctrsm
#define ztrsm      esvztrsm
 
#define stri       esvstri
#define dtri       esvdtri
 
#define strtri     esvstrtri
#define dtrtri     esvdtrtri
#define ctrtri     esvctrtri
#define ztrtri     esvztrtri
 
#define stpi       esvstpi
#define dtpi       esvdtpi
 
#define stptri     esvstptri
#define dtptri     esvdtptri
#define ctptri     esvctptri
#define ztptri     esvztptri
 
/*  Banded Linear Algebraic Equation Subroutines  */
 
#define sgbf       esvsgbf
#define dgbf       esvdgbf
 
#define sgbs       esvsgbs
#define dgbs       esvdgbs
 
#define spbf       esvspbf
#define dpbf       esvdpbf
 
#define spbchf     esvspbchf
#define dpbchf     esvdpbchf
 
#define spbs       esvspbs
#define dpbs       esvdpbs
 
#define spbchs     esvspbchs
#define dpbchs     esvdpbchs
 
#define sgtf       esvsgtf
#define dgtf       esvdgtf
 
#define sgts       esvsgts
#define dgts       esvdgts
 
#define sgtnp      esvsgtnp
#define dgtnp      esvdgtnp
#define cgtnp      esvcgtnp
#define zgtnp      esvzgtnp
 
#define sgtnpf     esvsgtnpf
#define dgtnpf     esvdgtnpf
#define cgtnpf     esvcgtnpf
#define zgtnpf     esvzgtnpf
 
#define sgtnps     esvsgtnps
#define dgtnps     esvdgtnps
#define cgtnps     esvcgtnps
#define zgtnps     esvzgtnps
 
#define sptf       esvsptf
#define dptf       esvdptf
 
#define spts       esvspts
#define dpts       esvdpts
 
#define stbsv      esvstbsv
#define dtbsv      esvdtbsv
#define ctbsv      esvctbsv
#define ztbsv      esvztbsv
 
/* Sparse Linear Algebraic Equations Subroutines   */
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define dgsf       esvdgsf_er
  #define dgss       esvdgss_er
  #define dgkfs      esvdgkfs_er
  #define dgkfsp     esvdgkfsp_er
  #define dskfs      esvdskfs_er
  #define dskfsp     esvdskfsp_er
  #define dsris      esvdsris_er
  #define dsmcg      esvdsmcg_er
  #define dsdcg      esvdsdcg_er
  #define dsmgcg     esvdsmgcg_er
  #define dsdgcg     esvdsdgcg_er
#else
  #define dgsf       esvdgsf
  #define dgss       esvdgss
  #define dgkfs      esvdgkfs
  #define dgkfsp     esvdgkfsp
  #define dskfs      esvdskfs
  #define dskfsp     esvdskfsp
  #define dsris      esvdsris
  #define dsmcg      esvdsmcg
  #define dsdcg      esvdsdcg
  #define dsmgcg     esvdsmgcg
  #define dsdgcg     esvdsdgcg
#endif
 
/* Linear Least Square Subroutines  */
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define sgesvf     esvsgesvf_er
  #define dgesvf     esvdgesvf_er
#else
  #define sgesvf     esvsgesvf
  #define dgesvf     esvdgesvf
#endif
 
#define sgesvs     esvsgesvs
#define dgesvs     esvdgesvs
 
#define sgeqrf     esvsgeqrf
#define dgeqrf     esvdgeqrf
#define cgeqrf     esvcgeqrf
#define zgeqrf     esvzgeqrf
 
#define sgels      esvsgels
#define dgels      esvdgels
#define cgels      esvcgels
#define zgels      esvzgels
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define sgells     esvsgells_er
  #define dgells     esvdgells_er
#else
  #define sgells     esvsgells
  #define dgells     esvdgells
#endif
 
/* Eigensystem Analysis Subroutines   */
 
#define sgeevx     esvsgeevx
#define dgeevx     esvdgeevx
#define cgeevx     esvcgeevx
#define zgeevx     esvzgeevx
 
#define sspevx     esvsspevx
#define dspevx     esvdspevx
#define chpevx     esvchpevx
#define zhpevx     esvzhpevx
 
#define ssyevx     esvssyevx
#define dsyevx     esvdsyevx
#define cheevx     esvcheevx
#define zheevx     esvzheevx
 
#define dsygvx     esvdsygvx
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define sgeev      esvsgeev_er
  #define dgeev      esvdgeev_er
  #define cgeev      esvcgeev_er
  #define zgeev      esvzgeev_er
 
  #define sspev      esvsspev_er
  #define dspev      esvdspev_er
  #define chpev      esvchpev_er
  #define zhpev      esvzhpev_er
 
  #define sslev      esvsspev_er
  #define dslev      esvdspev_er
  #define chlev      esvchpev_er
  #define zhlev      esvzhpev_er
 
  #define sspsv      esvsspsv_er
  #define dspsv      esvdspsv_er
  #define chpsv      esvchpsv_er
  #define zhpsv      esvzhpsv_er
 
  #define sgegv      esvsgegv_er
  #define dgegv      esvdgegv_er
 
  #define ssygv      esvssygv_er
  #define dsygv      esvdsygv_er
 
#else
  #define sgeev      esvsgeev
  #define dgeev      esvdgeev
  #define cgeev      esvcgeev
  #define zgeev      esvzgeev
 
  #define sspev      esvsspev
  #define dspev      esvdspev
  #define chpev      esvchpev
  #define zhpev      esvzhpev
 
  #define sslev      esvsspev
  #define dslev      esvdspev
  #define chlev      esvchpev
  #define zhlev      esvzhpev
 
  #define sspsv      esvsspsv
  #define dspsv      esvdspsv
  #define chpsv      esvchpsv
  #define zhpsv      esvzhpsv
 
  #define sgegv      esvsgegv
  #define dgegv      esvdgegv
 
  #define ssygv      esvssygv
  #define dsygv      esvdsygv
 
#endif
 
/*  Fourier Transforms Subroutines   */
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define scftd      esvscftd_er
  #define dcftd      esvdcftd_er
 
  #define srcftd     esvsrcftd_er
  #define drcftd     esvdrcftd_er
 
  #define scrftd     esvscrftd_er
  #define dcrftd     esvdcrftd_er
 
  #define scft       esvscft_er
  #define dcft       esvdcft_er
  #define scftp      esvscftp_er
 
  #define srcft      esvsrcft_er
  #define drcft      esvdrcft_er
 
  #define scrft      esvscrft_er
  #define dcrft      esvdcrft_er
 
  #define scosf      esvscosf_er
  #define dcosf      esvdcosf_er
  #define scosft     esvscosft_er
 
  #define ssinf      esvssinf_er
  #define dsinf      esvdsinf_er
 
  #define scft2      esvscft2_er
  #define dcft2      esvdcft2_er
  #define scft2p     esvscft2p_er
 
  #define srcft2     esvsrcft2_er
  #define drcft2     esvdrcft2_er
 
  #define scrft2     esvscrft2_er
  #define dcrft2     esvdcrft2_er
 
  #define scft3      esvscft3_er
  #define dcft3      esvdcft3_er
  #define scft3p     esvscft3p_er
 
  #define srcft3     esvsrcft3_er
  #define drcft3     esvdrcft3_er
 
  #define scrft3     esvscrft3_er
  #define dcrft3     esvdcrft3_er
#else
  #define scftd      esvscftd
  #define dcftd      esvdcftd
 
  #define srcftd     esvsrcftd
  #define drcftd     esvdrcftd
 
  #define scrftd     esvscrftd
  #define dcrftd     esvdcrftd
 
  #define scft       esvscft
  #define dcft       esvdcft
  #define scftp      esvscftp
 
  #define srcft      esvsrcft
  #define drcft      esvdrcft
 
  #define scrft      esvscrft
  #define dcrft      esvdcrft
 
  #define scosf      esvscosf
  #define dcosf      esvdcosf
  #define scosft     esvscosft
 
  #define ssinf      esvssinf
  #define dsinf      esvdsinf
 
  #define scft2      esvscft2
  #define dcft2      esvdcft2
  #define scft2p     esvscft2p
 
  #define srcft2     esvsrcft2
  #define drcft2     esvdrcft2
 
  #define scrft2     esvscrft2
  #define dcrft2     esvdcrft2
 
  #define scft3      esvscft3
  #define dcft3      esvdcft3
  #define scft3p     esvscft3p
 
  #define srcft3     esvsrcft3
  #define drcft3     esvdrcft3
 
  #define scrft3     esvscrft3
  #define dcrft3     esvdcrft3
#endif
 
/*  Convulutions/Correlation Subroutines   */
 
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define scon       esvscon_er
  #define scor       esvscor_er
#else
  #define scon       esvscon
  #define scor       esvscor
#endif
 
#define scond      esvscond
#define scord      esvscord
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define sconf      esvsconf_er
  #define scorf      esvscorf_er
#else
  #define sconf      esvsconf
  #define scorf      esvscorf
#endif
 
#define sdcon      esvsdcon
#define ddcon      esvddcon
#define sdcor      esvsdcor
#define ddcor      esvddcor
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define sacor      esvsacor_er
 
  #define sacorf     esvsacorf_er
#else
  #define sacor      esvsacor
 
  #define sacorf     esvsacorf
#endif
 
/*  Related Computations Subroutines   */
 
#define spoly      esvspoly
#define dpoly      esvdpoly
 
#define sizc       esvsizc
#define dizc       esvdizc
 
#define strec      esvstrec
#define dtrec      esvdtrec
 
#define sqint      esvsqint
#define dqint      esvdqint
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define swlev      esvswlev_er
  #define dwlev      esvdwlev_er
  #define cwlev      esvcwlev_er
  #define zwlev      esvzwlev_er
#else
  #define swlev      esvswlev
  #define dwlev      esvdwlev
  #define cwlev      esvcwlev
  #define zwlev      esvzwlev
#endif
 
/*  Sorting and Searching Subroutines  */
 
#define isort      esvisort
#define ssort      esvssort
#define dsort      esvdsort
 
#define isortx     esvisortx
#define ssortx     esvssortx
#define dsortx     esvdsortx
 
#define isorts     esvisorts
#define ssorts     esvssorts
#define dsorts     esvdsorts
 
#define ibsrch     esvibsrch
#define sbsrch     esvsbsrch
#define dbsrch     esvdbsrch
 
#define issrch     esvissrch
#define sssrch     esvsssrch
#define dssrch     esvdssrch
 
/*  Interpolation Subroutines  */
 
#define spint      esvspint
#define dpint      esvdpint
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define stpint     esvstpint_er
  #define dtpint     esvdtpint_er
#else
  #define stpint     esvstpint
  #define dtpint     esvdtpint
#endif
 
#define scsint     esvscsint
#define dcsint     esvdcsint
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define scsin2     esvscsin2_er
  #define dcsin2     esvdcsin2_er
#else
  #define scsin2     esvscsin2
  #define dcsin2     esvdcsin2
#endif
 
/*  Numerical Quadrature Subroutines  */
 
#define sptnq     esvsptnq
#define dptnq     esvdptnq
 
#define sglnq     esvsglnq
#define dglnq     esvdglnq
 
#define sglnq2    esvsglnq2
#define dglnq2    esvdglnq2
 
#define sglgq     esvsglgq
#define dglgq     esvdglgq
 
#define sgraq     esvsgraq
#define dgraq     esvdgraq
 
#define sghmq     esvsghmq
#define dghmq     esvdghmq
 
/*  Random Number Generation Subroutines  */
 
#define surand    esvsurand
#define durand    esvdurand
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define snrand    esvsnrand_er
  #define dnrand    esvdnrand_er
#else
  #define snrand    esvsnrand
  #define dnrand    esvdnrand
#endif
 
#define surxor    esvsurxor
#define durxor    esvdurxor
 
/*  Utility Subroutines  */
 
#define einfo     esveinfo
#define errset_   errset
#define errsav_   errsav
#define errstr_   errstr
 
#define ivsset    esvivsset
#define ievops    esvievops
 
#define stride    esvstride
 
#define dsrsm     esvdsrsm
 
#if defined(__ESVERR) || defined(_ESVERR)
  #define dgktrn    esvdgktrn_er
 
  #define dsktrn    esvdsktrn_er
#else
  #define dgktrn    esvdgktrn
 
  #define dsktrn    esvdsktrn
#endif
 
/*  Vector-Scalar Linear Algebra Subprograms  */
 
_ESVINT    esvisamax(_ESVINT, const   float *, _ESVINT);
_ESVINT    esvidamax(_ESVINT, const  double *, _ESVINT);
_ESVINT    esvicamax(_ESVINT, const  _ESVCM *, _ESVINT);
_ESVINT    esvizamax(_ESVINT, const _ESVCOM *, _ESVINT);
 
_ESVINT    esvisamin(_ESVINT, const  float *, _ESVINT);
_ESVINT    esvidamin(_ESVINT, const double *, _ESVINT);
 
_ESVINT    esvismax(_ESVINT, const  float *, _ESVINT);
_ESVINT    esvidmax(_ESVINT, const double *, _ESVINT);
 
_ESVINT    esvismin(_ESVINT, const  float *, _ESVINT);
_ESVINT    esvidmin(_ESVINT, const double *, _ESVINT);
 
float   esvsasum(_ESVINT, const   float *, _ESVINT);
double  esvdasum(_ESVINT, const  double *, _ESVINT);
float  esvscasum(_ESVINT, const  _ESVCM *, _ESVINT);
double esvdzasum(_ESVINT, const _ESVCOM *, _ESVINT);
 
void   esvsaxpy(_ESVINT,   float,   float *, _ESVINT,   float *,  _ESVINT);
void   esvdaxpy(_ESVINT,  double,  double *, _ESVINT,  double *,  _ESVINT);
void   esvcaxpy(_ESVINT,  _ESVCM,  _ESVCM *, _ESVINT,  _ESVCM *,  _ESVINT);
void   esvzaxpy(_ESVINT, _ESVCOM, _ESVCOM *, _ESVINT, _ESVCOM *,  _ESVINT);
 
void   esvscopy(_ESVINT,   float *, _ESVINT,   float *, _ESVINT);
void   esvdcopy(_ESVINT,  double *, _ESVINT,  double *, _ESVINT);
void   esvccopy(_ESVINT,  _ESVCM *, _ESVINT,  _ESVCM *, _ESVINT);
void   esvzcopy(_ESVINT, _ESVCOM *, _ESVINT, _ESVCOM *, _ESVINT);
 
float  esvsdot(_ESVINT, const  float *, _ESVINT, const  float *, _ESVINT);
double esvddot(_ESVINT, const double *, _ESVINT, const double *, _ESVINT);
 
_ESVCM  esvcdotu(_ESVINT, const  _ESVCM *, _ESVINT, const  _ESVCM *, _ESVINT);
_ESVCOM esvzdotu(_ESVINT, const _ESVCOM *, _ESVINT, const _ESVCOM *, _ESVINT);
_ESVCM  esvcdotc(_ESVINT, const  _ESVCM *, _ESVINT, const  _ESVCM *, _ESVINT);
_ESVCOM esvzdotc(_ESVINT, const _ESVCOM *, _ESVINT, const _ESVCOM *, _ESVINT);
 
void   esvsnaxpy(_ESVINT, _ESVINT, const  float *, _ESVINT, const void *, _ESVINT, _ESVINT,
                 void *, _ESVINT, _ESVINT);
void   esvdnaxpy(_ESVINT, _ESVINT, const double *, _ESVINT, const void *, _ESVINT, _ESVINT,
                 void *, _ESVINT, _ESVINT );
 
void   esvsndot(_ESVINT, _ESVINT,  float *, _ESVINT, _ESVINT, const void *, _ESVINT, _ESVINT,
                const void *, _ESVINT, _ESVINT);
void   esvdndot(_ESVINT, _ESVINT, double *, _ESVINT, _ESVINT, const void *, _ESVINT, _ESVINT,
                const void *, _ESVINT, _ESVINT);
 
float   esvsnrm2(_ESVINT, const   float *, _ESVINT);
double  esvdnrm2(_ESVINT, const  double *, _ESVINT);
float  esvscnrm2(_ESVINT, const  _ESVCM *, _ESVINT);
double esvdznrm2(_ESVINT, const _ESVCOM *, _ESVINT);
 
float  esvsnorm2(_ESVINT, const   float *, _ESVINT);
double esvdnorm2(_ESVINT, const  double *, _ESVINT);
float  esvcnorm2(_ESVINT, const  _ESVCM *, _ESVINT);
double esvznorm2(_ESVINT, const _ESVCOM *, _ESVINT);
 
void   esvsrotg(_ESVS, _ESVS, _ESVS, _ESVS);
void   esvdrotg(_ESVD, _ESVD, _ESVD, _ESVD);
void   esvcrotg(_ESVC, _ESVC, _ESVS, _ESVC);
void   esvzrotg(_ESVZ, _ESVZ, _ESVD, _ESVZ);
 
void    esvsrot(_ESVINT,   float *, _ESVINT,   float *, _ESVINT,  float,   float);
void    esvdrot(_ESVINT,  double *, _ESVINT,  double *, _ESVINT, double,  double);
void    esvcrot(_ESVINT,  _ESVCM *, _ESVINT,  _ESVCM *, _ESVINT,  float,  _ESVCM);
void    esvzrot(_ESVINT, _ESVCOM *, _ESVINT, _ESVCOM *, _ESVINT, double, _ESVCOM);
void   esvcsrot(_ESVINT,  _ESVCM *, _ESVINT,  _ESVCM *, _ESVINT,  float,   float);
void   esvzdrot(_ESVINT, _ESVCOM *, _ESVINT, _ESVCOM *, _ESVINT, double,  double);
 
void    esvsscal(_ESVINT,   float,   float *, _ESVINT);
void    esvdscal(_ESVINT,  double,  double *, _ESVINT);
void    esvcscal(_ESVINT,  _ESVCM,  _ESVCM *, _ESVINT);
void    esvzscal(_ESVINT, _ESVCOM, _ESVCOM *, _ESVINT);
void   esvcsscal(_ESVINT,   float,  _ESVCM *, _ESVINT);
void   esvzdscal(_ESVINT,  double, _ESVCOM *, _ESVINT);
 
void   esvsswap(_ESVINT,   float *, _ESVINT,   float *, _ESVINT);
void   esvdswap(_ESVINT,  double *, _ESVINT,  double *, _ESVINT);
void   esvcswap(_ESVINT,  _ESVCM *, _ESVINT,  _ESVCM *, _ESVINT);
void   esvzswap(_ESVINT, _ESVCOM *, _ESVINT, _ESVCOM *, _ESVINT);
 
void   esvsvea(_ESVINT,   float *, _ESVINT,   float *, _ESVINT,   float *, _ESVINT);
void   esvdvea(_ESVINT,  double *, _ESVINT,  double *, _ESVINT,  double *, _ESVINT);
void   esvcvea(_ESVINT,  _ESVCM *, _ESVINT,  _ESVCM *, _ESVINT,  _ESVCM *, _ESVINT);
void   esvzvea(_ESVINT, _ESVCOM *, _ESVINT, _ESVCOM *, _ESVINT, _ESVCOM *, _ESVINT);
 
void   esvsves(_ESVINT,   float *, _ESVINT,   float *, _ESVINT,   float *, _ESVINT);
void   esvdves(_ESVINT,  double *, _ESVINT,  double *, _ESVINT,  double *, _ESVINT);
void   esvcves(_ESVINT,  _ESVCM *, _ESVINT,  _ESVCM *, _ESVINT,  _ESVCM *, _ESVINT);
void   esvzves(_ESVINT, _ESVCOM *, _ESVINT, _ESVCOM *, _ESVINT, _ESVCOM *, _ESVINT);
 
void   esvsvem(_ESVINT,   float *, _ESVINT,   float *, _ESVINT,   float *, _ESVINT);
void   esvdvem(_ESVINT,  double *, _ESVINT,  double *, _ESVINT,  double *, _ESVINT);
void   esvcvem(_ESVINT,  _ESVCM *, _ESVINT,  _ESVCM *, _ESVINT,  _ESVCM *, _ESVINT);
void   esvzvem(_ESVINT, _ESVCOM *, _ESVINT, _ESVCOM *, _ESVINT, _ESVCOM *, _ESVINT);
 
void    esvsyax(_ESVINT,   float,   float *, _ESVINT,   float *, _ESVINT);
void    esvdyax(_ESVINT,  double,  double *, _ESVINT,  double *, _ESVINT);
void    esvcyax(_ESVINT,  _ESVCM,  _ESVCM *, _ESVINT,  _ESVCM *, _ESVINT);
void    esvzyax(_ESVINT, _ESVCOM, _ESVCOM *, _ESVINT, _ESVCOM *, _ESVINT);
void   esvcsyax(_ESVINT,   float,  _ESVCM *, _ESVINT,  _ESVCM *, _ESVINT);
void   esvzdyax(_ESVINT,  double, _ESVCOM *, _ESVINT, _ESVCOM *, _ESVINT);
 
void   esvszaxpy(_ESVINT,   float,   float *, _ESVINT,   float *, _ESVINT,   float *, _ESVINT);
void   esvdzaxpy(_ESVINT,  double,  double *, _ESVINT,  double *, _ESVINT,  double *, _ESVINT);
void   esvczaxpy(_ESVINT,  _ESVCM,  _ESVCM *, _ESVINT,  _ESVCM *, _ESVINT,  _ESVCM *, _ESVINT);
void   esvzzaxpy(_ESVINT, _ESVCOM, _ESVCOM *, _ESVINT, _ESVCOM *, _ESVINT, _ESVCOM *, _ESVINT);
 
/*  Sparse Vector-Scalar Linear Algebra Subrprograms  */
 
void   esvssctr(_ESVINT, const   float *, const _ESVINT *,   float *);
void   esvdsctr(_ESVINT, const  double *, const _ESVINT *,  double *);
void   esvcsctr(_ESVINT, const  _ESVCM *, const _ESVINT *,  _ESVCM *);
void   esvzsctr(_ESVINT, const _ESVCOM *, const _ESVINT *, _ESVCOM *);
 
void   esvsgthr(_ESVINT, const   float *,   float *, const _ESVINT *);
void   esvdgthr(_ESVINT, const  double *,  double *, const _ESVINT *);
void   esvcgthr(_ESVINT, const  _ESVCM *,  _ESVCM *, const _ESVINT *);
void   esvzgthr(_ESVINT, const _ESVCOM *, _ESVCOM *, const _ESVINT *);
 
void   esvsgthrz(_ESVINT,   float *,   float *, const _ESVINT *);
void   esvdgthrz(_ESVINT,  double *,  double *, const _ESVINT *);
void   esvcgthrz(_ESVINT,  _ESVCM *,  _ESVCM *, const _ESVINT *);
void   esvzgthrz(_ESVINT, _ESVCOM *, _ESVCOM *, const _ESVINT *);
 
void   esvsaxpyi(_ESVINT,   float, const   float *, const _ESVINT *,   float *);
void   esvdaxpyi(_ESVINT,  double, const  double *, const _ESVINT *,  double *);
void   esvcaxpyi(_ESVINT,  _ESVCM, const  _ESVCM *, const _ESVINT *,  _ESVCM *);
void   esvzaxpyi(_ESVINT, _ESVCOM, const _ESVCOM *, const _ESVINT *, _ESVCOM *);
 
float  esvsdoti(_ESVINT, const  float *, const _ESVINT *, const  float *);
double esvddoti(_ESVINT, const double *, const _ESVINT *, const double *);
 
_ESVCM  esvcdotui(_ESVINT, const  _ESVCM *, const _ESVINT *, const  _ESVCM *);
_ESVCOM esvzdotui(_ESVINT, const _ESVCOM *, const _ESVINT *, const _ESVCOM *);
 
_ESVCM  esvcdotci(_ESVINT, const  _ESVCM *, const _ESVINT *, const  _ESVCM *);
_ESVCOM esvzdotci(_ESVINT, const _ESVCOM *, const _ESVINT *, const _ESVCOM *);
 
/*  Dense Matrix-Vector Subprograms  */
 
void   esvsgemv(const char *, _ESVINT, _ESVINT,   float, const void *, _ESVINT,
                const   float *, _ESVINT,   float,   float *, _ESVINT);
void   esvdgemv(const char *, _ESVINT, _ESVINT,  double, const void *, _ESVINT,
                const  double *, _ESVINT,  double,  double *, _ESVINT);
void   esvcgemv(const char *, _ESVINT, _ESVINT,  _ESVCM, const void *, _ESVINT,
                const  _ESVCM *, _ESVINT,  _ESVCM,  _ESVCM *, _ESVINT);
void   esvzgemv(const char *, _ESVINT, _ESVINT, _ESVCOM, const void *, _ESVINT,
                const _ESVCOM *, _ESVINT, _ESVCOM, _ESVCOM *, _ESVINT);
 
void   esvsgemx(_ESVINT, _ESVINT,  float, const void *, _ESVINT,
                const  float *, _ESVINT,  float *, _ESVINT);
void   esvdgemx( _ESVINT, _ESVINT, double,const void *, _ESVINT,
                const double *, _ESVINT, double *, _ESVINT);
 
void   esvsgemtx(_ESVINT, _ESVINT,  float, const void *, _ESVINT,
                 const  float *, _ESVINT,  float *, _ESVINT);
void   esvdgemtx(_ESVINT, _ESVINT, double, const void *, _ESVINT,
                 const double *, _ESVINT, double *, _ESVINT);
 
void   esvsger1(_ESVINT, _ESVINT,   float, const   float *, _ESVINT,
                const   float *, _ESVINT, void *, _ESVINT);
void   esvdger1(_ESVINT, _ESVINT,  double, const  double *, _ESVINT,
                const  double *, _ESVINT, void *, _ESVINT);
void   esvcgeru(_ESVINT, _ESVINT,  _ESVCM, const  _ESVCM *, _ESVINT,
                const  _ESVCM *, _ESVINT, void *, _ESVINT);
void   esvzgeru(_ESVINT, _ESVINT, _ESVCOM, const _ESVCOM *, _ESVINT,
                const _ESVCOM *, _ESVINT, void *, _ESVINT);
void   esvcgerc(_ESVINT, _ESVINT,  _ESVCM, const  _ESVCM *, _ESVINT,
                const  _ESVCM *, _ESVINT, void *, _ESVINT);
void   esvzgerc(_ESVINT, _ESVINT, _ESVCOM, const _ESVCOM *, _ESVINT,
                const _ESVCOM *, _ESVINT, void *, _ESVINT);
 
void   esvsspmv(const char *, _ESVINT,   float, const   float *,
                const   float *, _ESVINT,   float,   float *, _ESVINT);
void   esvdspmv(const char *, _ESVINT,  double, const  double *,
                const  double *, _ESVINT,  double,  double *, _ESVINT);
void   esvchpmv(const char *, _ESVINT,  _ESVCM, const  _ESVCM *,
                const  _ESVCM *, _ESVINT,  _ESVCM,  _ESVCM *, _ESVINT);
void   esvzhpmv(const char *, _ESVINT, _ESVCOM, const _ESVCOM *,
                const _ESVCOM *, _ESVINT, _ESVCOM, _ESVCOM *, _ESVINT);
 
void   esvssymv(const char *, _ESVINT,   float, const void *, _ESVINT,
                const   float *, _ESVINT,   float,   float *, _ESVINT);
void   esvdsymv(const char *, _ESVINT,  double, const void *, _ESVINT,
                const  double *, _ESVINT,  double,  double *, _ESVINT);
void   esvchemv(const char *, _ESVINT,  _ESVCM, const void *, _ESVINT,
                const  _ESVCM *, _ESVINT,  _ESVCM,  _ESVCM *, _ESVINT);
void   esvzhemv(const char *, _ESVINT, _ESVCOM, const void *, _ESVINT,
                const _ESVCOM *, _ESVINT, _ESVCOM, _ESVCOM *, _ESVINT);
 
void   esvsslmx(_ESVINT,  float, const  float *,
                const  float *, _ESVINT,  float *, _ESVINT);
void   esvdslmx(_ESVINT, double, const double *,
                const double *, _ESVINT, double *, _ESVINT);
 
void   esvsspr(const char *, _ESVINT,  float, const   float *, _ESVINT,   float *);
void   esvdspr(const char *, _ESVINT, double, const  double *, _ESVINT,  double *);
void   esvchpr(const char *, _ESVINT,  float, const  _ESVCM *, _ESVINT,  _ESVCM *);
void   esvzhpr(const char *, _ESVINT, double, const _ESVCOM *, _ESVINT, _ESVCOM *);
 
void   esvssyr(const char *, _ESVINT,  float, const   float *, _ESVINT, void *, _ESVINT);
void   esvdsyr(const char *, _ESVINT, double, const  double *, _ESVINT, void *, _ESVINT);
void   esvcher(const char *, _ESVINT,  float, const  _ESVCM *, _ESVINT, void *, _ESVINT);
void   esvzher(const char *, _ESVINT, double, const _ESVCOM *, _ESVINT, void *, _ESVINT);
 
void   esvsslr1(_ESVINT,  float, const  float *, _ESVINT,  float *);
void   esvdslr1(_ESVINT, double, const double *, _ESVINT, double *);
 
void   esvsspr2(const char *, _ESVINT,   float, const   float *, _ESVINT,
                const   float *, _ESVINT,   float *);
void   esvdspr2(const char *, _ESVINT,  double, const  double *, _ESVINT,
                const  double *, _ESVINT,  double *);
void   esvchpr2(const char *, _ESVINT,  _ESVCM, const  _ESVCM *, _ESVINT,
                const  _ESVCM *, _ESVINT,  _ESVCM *);
void   esvzhpr2(const char *, _ESVINT, _ESVCOM, const _ESVCOM *, _ESVINT,
                const _ESVCOM *, _ESVINT, _ESVCOM *);
 
void   esvssyr2(const char *, _ESVINT,   float, const   float *, _ESVINT,
                const   float *,  _ESVINT, void *, _ESVINT);
void   esvdsyr2(const char *, _ESVINT,  double, const  double *, _ESVINT,
                const  double *, _ESVINT, void *, _ESVINT);
void   esvcher2(const char *, _ESVINT,  _ESVCM, const  _ESVCM *, _ESVINT,
                const  _ESVCM *, _ESVINT, void *, _ESVINT);
void   esvzher2(const char *, _ESVINT, _ESVCOM, const _ESVCOM *, _ESVINT,
                const _ESVCOM *, _ESVINT, void *, _ESVINT);
 
void   esvsslr2(_ESVINT,  float, const  float *, _ESVINT,
                const  float *, _ESVINT,  float *);
void   esvdslr2(_ESVINT, double, const double *, _ESVINT,
                const double *, _ESVINT, double *);
 
void   esvsgbmv(const char *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,   float, const void *, _ESVINT,
                const   float *, _ESVINT,   float,   float *, _ESVINT);
void   esvdgbmv(const char *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,  double, const void *, _ESVINT,
                const  double *, _ESVINT,  double,  double *, _ESVINT);
void   esvcgbmv(const char *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,  _ESVCM, const void *, _ESVINT,
                const  _ESVCM *, _ESVINT,  _ESVCM,  _ESVCM *, _ESVINT);
void   esvzgbmv(const char *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVCOM, const void *, _ESVINT,
                const _ESVCOM *, _ESVINT, _ESVCOM, _ESVCOM *, _ESVINT);
 
void   esvssbmv(const char *, _ESVINT, _ESVINT,   float, const void *, _ESVINT,
                const   float *, _ESVINT,   float,   float *, _ESVINT);
void   esvdsbmv(const char *, _ESVINT, _ESVINT,  double, const void *, _ESVINT,
                const  double *, _ESVINT,  double,  double *, _ESVINT);
void   esvchbmv(const char *, _ESVINT, _ESVINT,  _ESVCM, const void *, _ESVINT,
                const  _ESVCM *, _ESVINT,  _ESVCM,  _ESVCM *, _ESVINT);
void   esvzhbmv(const char *, _ESVINT, _ESVINT, _ESVCOM, const void *, _ESVINT,
                const _ESVCOM *, _ESVINT, _ESVCOM, _ESVCOM *, _ESVINT);
 
void   esvstrmv(const char *, const char *, const char *, _ESVINT,
                const void *, _ESVINT,   float *, _ESVINT);
void   esvdtrmv(const char *, const char *, const char *, _ESVINT,
                const void *, _ESVINT,  double *, _ESVINT);
void   esvctrmv(const char *, const char *, const char *, _ESVINT,
                const void *, _ESVINT,  _ESVCM *, _ESVINT);
void   esvztrmv(const char *, const char *, const char *, _ESVINT,
                const void *, _ESVINT, _ESVCOM *, _ESVINT);
 
void   esvstpmv(const char *, const char *, const char *,  _ESVINT,
                const   float *,   float *,  _ESVINT);
void   esvdtpmv(const char *, const char *, const char *,  _ESVINT,
                const  double *,  double *,  _ESVINT);
void   esvctpmv(const char *, const char *, const char *,  _ESVINT,
                const  _ESVCM *,  _ESVCM *,  _ESVINT);
void   esvztpmv(const char *, const char *, const char *,  _ESVINT,
                const _ESVCOM *, _ESVCOM *,  _ESVINT);
 
void   esvstbmv(const char *, const char *, const char *, _ESVINT, _ESVINT,
                const void *, _ESVINT,   float *, _ESVINT);
void   esvdtbmv(const char *, const char *, const char *, _ESVINT, _ESVINT,
                const void *, _ESVINT,  double *, _ESVINT);
void   esvctbmv(const char *, const char *, const char *, _ESVINT, _ESVINT,
                const void *, _ESVINT,  _ESVCM *, _ESVINT);
void   esvztbmv(const char *, const char *, const char *, _ESVINT, _ESVINT,
                const void *, _ESVINT, _ESVCOM *, _ESVINT);
 
/*  Sparse Matrix-Vector Subprograms  */
 
void   esvdsmmx(_ESVINT, _ESVINT, const void *, const void *, _ESVINT,
               const double *, double *);
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvdsmtm_er(_ESVINT, _ESVINT, const void *, const void *,
                    _ESVINT, _ESVI, _ESVI, void *, void *,
                    _ESVINT, float *, _ESVI);
#else
 void   esvdsmtm   (_ESVINT, _ESVINT, const void *, const void *,
                    _ESVINT, _ESVI, _ESVI, void *, void *,
                    _ESVINT, float *, _ESVINT);
#endif
 
void   esvdsdmx(_ESVINT, _ESVINT, _ESVINT, const void *, _ESVINT,
                const char *, const _ESVINT *, const double *, double *);
 
/*  Matrix Operation Subroutines  */
 
void   esvsgeadd(void *, _ESVINT, const char *, void *, _ESVINT, const char *,
                 void *, _ESVINT, _ESVINT, _ESVINT);
void   esvdgeadd(void *, _ESVINT, const char *, void *, _ESVINT, const char *,
                 void *, _ESVINT, _ESVINT, _ESVINT);
void   esvcgeadd(void *, _ESVINT, const char *, void *, _ESVINT, const char *,
                 void *, _ESVINT, _ESVINT, _ESVINT);
void   esvzgeadd(void *, _ESVINT, const char *, void *, _ESVINT, const char *,
                 void *, _ESVINT, _ESVINT, _ESVINT);
 
void   esvsgesub(void *, _ESVINT, const char *, void *, _ESVINT, const char *,
                 void *, _ESVINT, _ESVINT, _ESVINT);
void   esvdgesub(void *, _ESVINT, const char *, void *, _ESVINT, const char *,
                 void *, _ESVINT, _ESVINT, _ESVINT);
void   esvcgesub(void *, _ESVINT, const char *, void *, _ESVINT, const char *,
                 void *, _ESVINT, _ESVINT, _ESVINT);
void   esvzgesub(void *, _ESVINT, const char *, void *, _ESVINT, const char *,
                 void *, _ESVINT, _ESVINT, _ESVINT);
 
void   esvsgemul(const void *, _ESVINT, const char *, const void *, _ESVINT,
                 const char *, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT);
void   esvdgemul(const void *, _ESVINT, const char *, const void *, _ESVINT,
                 const char *, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT);
void   esvcgemul(const void *, _ESVINT, const char *, const void *, _ESVINT,
                 const char *, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT);
void   esvzgemul(const void *, _ESVINT, const char *, const void *, _ESVINT,
                 const char *, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT);
 
void   esvdgemlp(const void *, _ESVINT, const char *, const void *, _ESVINT,
                 const char *, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT);
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvsgemms_er(const void *, _ESVINT, const char *, const void *, _ESVINT,
                     const char *, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,  float *, _ESVI);
 _ESVINT    esvdgemms_er(const void *, _ESVINT, const char *, const void *, _ESVINT,
                     const char *, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, double *, _ESVI);
 _ESVINT    esvcgemms_er(const void *, _ESVINT, const char *, const void *, _ESVINT,
                     const char *, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,  float *, _ESVI);
 _ESVINT    esvzgemms_er(const void *, _ESVINT, const char *, const void *, _ESVINT,
                     const char *, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, double *, _ESVI);
#else
 void   esvsgemms   (const void *, _ESVINT, const char *, const void *, _ESVINT,
                     const char *, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,  float *, _ESVINT);
 void   esvdgemms   (const void *, _ESVINT, const char *, const void *, _ESVINT,
                     const char *, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, double *, _ESVINT);
 void   esvcgemms   (const void *, _ESVINT, const char *, const void *, _ESVINT,
                     const char *, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,  float *, _ESVINT);
 void   esvzgemms   (const void *, _ESVINT, const char *, const void *, _ESVINT,
                     const char *, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, double *, _ESVINT);
#endif
 
void   esvsgemm(const char *, const char *, _ESVINT, _ESVINT, _ESVINT,   float,
                const void *, _ESVINT, const void *, _ESVINT,   float, void *, _ESVINT);
void   esvdgemm(const char *, const char *, _ESVINT, _ESVINT, _ESVINT,  double,
                const void *, _ESVINT, const void *, _ESVINT,  double, void *, _ESVINT);
void   esvcgemm(const char *, const char *, _ESVINT, _ESVINT, _ESVINT,  _ESVCM,
                const void *, _ESVINT, const void *, _ESVINT,  _ESVCM, void *, _ESVINT);
void   esvzgemm(const char *, const char *, _ESVINT, _ESVINT, _ESVINT, _ESVCOM,
                const void *, _ESVINT, const void *, _ESVINT, _ESVCOM, void *, _ESVINT);
 
void   esvssymm(const char *, const char *, _ESVINT, _ESVINT,   float,
                const void *, _ESVINT, const void *, _ESVINT,   float, void *, _ESVINT);
void   esvdsymm(const char *, const char *, _ESVINT, _ESVINT,  double,
                const void *, _ESVINT, const void *, _ESVINT,  double, void *, _ESVINT);
void   esvcsymm(const char *, const char *, _ESVINT, _ESVINT,  _ESVCM,
                const void *, _ESVINT, const void *, _ESVINT,  _ESVCM, void *, _ESVINT);
void   esvzsymm(const char *, const char *, _ESVINT, _ESVINT, _ESVCOM,
                const void *, _ESVINT, const void *, _ESVINT, _ESVCOM, void *, _ESVINT);
void   esvchemm(const char *, const char *, _ESVINT, _ESVINT,  _ESVCM,
                const void *, _ESVINT, const void *, _ESVINT,  _ESVCM, void *, _ESVINT);
void   esvzhemm(const char *, const char *, _ESVINT, _ESVINT, _ESVCOM,
                const void *, _ESVINT, const void *, _ESVINT, _ESVCOM, void *, _ESVINT);
 
void   esvstrmm(const char *, const char *, const char *, const char *,
                _ESVINT, _ESVINT,   float, const void*, _ESVINT, void *, _ESVINT);
void   esvdtrmm(const char *, const char *, const char *, const char *,
                _ESVINT, _ESVINT,  double, const void*, _ESVINT, void *, _ESVINT);
void   esvctrmm(const char *, const char *, const char *, const char *,
                _ESVINT, _ESVINT,  _ESVCM, const void*, _ESVINT, void *, _ESVINT);
void   esvztrmm(const char *, const char *, const char *, const char *,
                _ESVINT, _ESVINT, _ESVCOM, const void*, _ESVINT, void *, _ESVINT);
 
void   esvssyrk(const char *, const char *, _ESVINT, _ESVINT,   float,
                const void *, _ESVINT,   float, void *, _ESVINT);
void   esvdsyrk(const char *, const char *, _ESVINT, _ESVINT,  double,
                const void *, _ESVINT,  double, void *, _ESVINT);
void   esvcsyrk(const char *, const char *, _ESVINT, _ESVINT,  _ESVCM,
                const void *, _ESVINT,  _ESVCM, void *, _ESVINT);
void   esvzsyrk(const char *, const char *, _ESVINT, _ESVINT, _ESVCOM,
                const void *, _ESVINT, _ESVCOM, void *, _ESVINT);
void   esvcherk(const char *, const char *, _ESVINT, _ESVINT,   float,
                const void *, _ESVINT,   float, void *, _ESVINT);
void   esvzherk(const char *, const char *, _ESVINT, _ESVINT,  double,
                const void *, _ESVINT,  double, void *, _ESVINT);
 
void   esvssyr2k(const char *, const char *, _ESVINT, _ESVINT,   float,
                 const void *, _ESVINT, const void *, _ESVINT,   float, void *, _ESVINT);
void   esvdsyr2k(const char *, const char *, _ESVINT, _ESVINT,  double,
                 const void *, _ESVINT, const void *, _ESVINT,  double, void *, _ESVINT);
void   esvcsyr2k(const char *, const char *, _ESVINT, _ESVINT,  _ESVCM,
                 const void *, _ESVINT, const void *, _ESVINT,  _ESVCM, void *, _ESVINT);
void   esvzsyr2k(const char *, const char *, _ESVINT, _ESVINT, _ESVCOM,
                 const void *, _ESVINT, const void *, _ESVINT, _ESVCOM, void *, _ESVINT);
void   esvcher2k(const char *, const char *, _ESVINT, _ESVINT,  _ESVCM,
                 const void *, _ESVINT, const void *, _ESVINT,   float, void *, _ESVINT);
void   esvzher2k(const char *, const char *, _ESVINT, _ESVINT, _ESVCOM,
                 const void *, _ESVINT, const void *, _ESVINT,  double, void *, _ESVINT);
 
void   esvsgetmi(void *, _ESVINT, _ESVINT);
void   esvdgetmi(void *, _ESVINT, _ESVINT);
void   esvcgetmi(void *, _ESVINT, _ESVINT);
void   esvzgetmi(void *, _ESVINT, _ESVINT);
 
void   esvsgetmo(const void *, _ESVINT, _ESVINT, _ESVINT, void *, _ESVINT);
void   esvdgetmo(const void *, _ESVINT, _ESVINT, _ESVINT, void *, _ESVINT);
void   esvcgetmo(const void *, _ESVINT, _ESVINT, _ESVINT, void *, _ESVINT);
void   esvzgetmo(const void *, _ESVINT, _ESVINT, _ESVINT, void *, _ESVINT);
 
/*  Dense Linear Algebraic Equation Subroutines  */
 
_ESVINT    esvsgef(void *, _ESVINT, _ESVINT, _ESVINT *);
_ESVINT    esvdgef(void *, _ESVINT, _ESVINT, _ESVINT *);
_ESVINT    esvcgef(void *, _ESVINT, _ESVINT, _ESVINT *);
_ESVINT    esvzgef(void *, _ESVINT, _ESVINT, _ESVINT *);
 
float    esvslange(const char *, _ESVINT, _ESVINT, const float   *, _ESVINT, float  *);
double   esvdlange(const char *, _ESVINT, _ESVINT, const double  *, _ESVINT, double *);
float    esvclange(const char *, _ESVINT, _ESVINT, const _ESVCM  *, _ESVINT, float  *);
double   esvzlange(const char *, _ESVINT, _ESVINT, const _ESVCOM *, _ESVINT, double *);
 
void   esvsgecon(const char *, _ESVINT, const float   *, _ESVINT, float , float  *, float   *, _ESVINT *, _ESVI);
void   esvdgecon(const char *, _ESVINT, const double  *, _ESVINT, double, double *, double  *, _ESVINT *, _ESVI);
void   esvcgecon(const char *, _ESVINT, const _ESVCM  *, _ESVINT, float , float  *, _ESVCM  *, float   *, _ESVI);
void   esvzgecon(const char *, _ESVINT, const _ESVCOM *, _ESVINT, double, double *, _ESVCOM *, double  *, _ESVI);
 
void   esvsgetrf(_ESVINT, _ESVINT, void *, _ESVINT, _ESVINT *, _ESVI);
void   esvdgetrf(_ESVINT, _ESVINT, void *, _ESVINT, _ESVINT *, _ESVI);
void   esvcgetrf(_ESVINT, _ESVINT, void *, _ESVINT, _ESVINT *, _ESVI);
void   esvzgetrf(_ESVINT, _ESVINT, void *, _ESVINT, _ESVINT *, _ESVI);
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvdgefp_er(void *, _ESVINT, _ESVINT, _ESVINT *, double *, _ESVI);
#else
 _ESVINT    esvdgefp   (void *, _ESVINT, _ESVINT, _ESVINT *, double *, _ESVINT);
#endif
 
void   esvsges(const void *, _ESVINT, _ESVINT, const _ESVINT *,   float *, _ESVINT);
void   esvdges(const void *, _ESVINT, _ESVINT, const _ESVINT *,  double *, _ESVINT);
void   esvcges(const void *, _ESVINT, _ESVINT, const _ESVINT *,  _ESVCM *, _ESVINT);
void   esvzges(const void *, _ESVINT, _ESVINT, const _ESVINT *, _ESVCOM *, _ESVINT);
 
void   esvsgesm(const char *, const void *, _ESVINT, _ESVINT, const _ESVINT *,
                void *, _ESVINT, _ESVINT);
void   esvdgesm(const char *, const void *, _ESVINT, _ESVINT, const _ESVINT *,
                void *, _ESVINT, _ESVINT);
void   esvcgesm(const char *, const void *, _ESVINT, _ESVINT, const _ESVINT *,
                void *, _ESVINT, _ESVINT);
void   esvzgesm(const char *, const void *, _ESVINT, _ESVINT, const _ESVINT *,
                void *, _ESVINT, _ESVINT);
 
void   esvsgetrs(const char *, _ESVINT, _ESVINT,const void *, _ESVINT, const _ESVINT *,
                 void *, _ESVINT,_ESVI);
void   esvdgetrs(const char *, _ESVINT, _ESVINT,const void *, _ESVINT, const _ESVINT *,
                 void *, _ESVINT,_ESVI);
void   esvcgetrs(const char *, _ESVINT, _ESVINT,const void *, _ESVINT, const _ESVINT *,
                 void *, _ESVINT,_ESVI);
void   esvzgetrs(const char *, _ESVINT, _ESVINT,const void *, _ESVINT, const _ESVINT *,
                 void *, _ESVINT,_ESVI);
 
void   esvsgesv(_ESVINT, _ESVINT,void *, _ESVINT, _ESVINT *,void *, _ESVINT,_ESVI);
void   esvdgesv(_ESVINT, _ESVINT,void *, _ESVINT, _ESVINT *,void *, _ESVINT,_ESVI);
void   esvcgesv(_ESVINT, _ESVINT,void *, _ESVINT, _ESVINT *,void *, _ESVINT,_ESVI);
void   esvzgesv(_ESVINT, _ESVINT,void *, _ESVINT, _ESVINT *,void *, _ESVINT,_ESVI);
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvsgefcd_er(void *, _ESVINT, _ESVINT, _ESVINT *, _ESVINT, _ESVS,  float *,  float *,
                     _ESVI);
 _ESVINT    esvdgefcd_er(void *, _ESVINT, _ESVINT, _ESVINT *, _ESVINT, _ESVD, double *, double *,
                     _ESVI);
#else
 _ESVINT    esvsgefcd   (void *, _ESVINT, _ESVINT, _ESVINT *, _ESVINT, _ESVS,  float *,  float *,
                     _ESVINT);
 _ESVINT    esvdgefcd   (void *, _ESVINT, _ESVINT, _ESVINT *, _ESVINT, _ESVD, double *, double *,
                     _ESVINT);
#endif
 
_ESVINT    esvsppf( float *, _ESVINT, _ESVINT);
_ESVINT    esvdppf(double *, _ESVINT, _ESVINT);
 
void    esvspptrf(const char *, _ESVINT, float  *, _ESVI);
void    esvdpptrf(const char *, _ESVINT, double *, _ESVI);
void    esvcpptrf(const char *, _ESVINT, _ESVCM *, _ESVI);
void    esvzpptrf(const char *, _ESVINT, _ESVCOM *, _ESVI);
 
float   esvslansp(const char *, const char *, _ESVINT, const float   *, float  *);
double  esvdlansp(const char *, const char *, _ESVINT, const double  *, double *);
float   esvclanhp(const char *, const char *, _ESVINT, const _ESVCM  *, float  *);
double  esvzlanhp(const char *, const char *, _ESVINT, const _ESVCOM *, double *);
 
void   esvsppcon(const char *, _ESVINT, const float   *, float,  float  *, float   *, _ESVINT *, _ESVI);
void   esvdppcon(const char *, _ESVINT, const double  *, double, double *, double  *, _ESVINT *, _ESVI);
void   esvcppcon(const char *, _ESVINT, const _ESVCM  *, float,  float  *, _ESVCM  *, float   *, _ESVI);
void   esvzppcon(const char *, _ESVINT, const _ESVCOM *, double, double *, _ESVCOM *, double  *, _ESVI);
 
void    esvspptri(const char *, _ESVINT, float   *, _ESVI);
void    esvdpptri(const char *, _ESVINT, double  *, _ESVI);
void    esvcpptri(const char *, _ESVINT, _ESVCM  *, _ESVI);
void    esvzpptri(const char *, _ESVINT, _ESVCOM *, _ESVI);
 
_ESVINT    esvspof(const char *, void *, _ESVINT, _ESVINT);
_ESVINT    esvdpof(const char *, void *, _ESVINT, _ESVINT);
_ESVINT    esvcpof(const char *, void *, _ESVINT, _ESVINT);
_ESVINT    esvzpof(const char *, void *, _ESVINT, _ESVINT);
 
void   esvspotrf(const char *, _ESVINT, void *, _ESVINT, _ESVI);
void   esvdpotrf(const char *, _ESVINT, void *, _ESVINT, _ESVI);
void   esvcpotrf(const char *, _ESVINT, void *, _ESVINT, _ESVI);
void   esvzpotrf(const char *, _ESVINT, void *, _ESVINT, _ESVI);
 
float   esvslansy(const char *, const char *, _ESVINT, const void *, _ESVINT, float  *);
double  esvdlansy(const char *, const char *, _ESVINT, const void *, _ESVINT, double *);
float   esvclanhe(const char *, const char *, _ESVINT, const void *, _ESVINT, float  *);
double  esvzlanhe(const char *, const char *, _ESVINT, const void *, _ESVINT, double *);
 
void   esvspocon(const char *, _ESVINT, const void *, _ESVINT, const float,  float  *, float  *, _ESVINT *, _ESVI);
void   esvdpocon(const char *, _ESVINT, const void *, _ESVINT, const double, double *, double *, _ESVINT *, _ESVI);
void   esvcpocon(const char *, _ESVINT, const void *, _ESVINT, const float,  float  *, float  *, float   *, _ESVI);
void   esvzpocon(const char *, _ESVINT, const void *, _ESVINT, const double, double *, double *, double  *, _ESVI);
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvdppfp_er(double *, _ESVINT, double *, _ESVI);
#else
 _ESVINT    esvdppfp   (double *, _ESVINT, double *, _ESVINT);
#endif
 
void   esvspps(const  float *, _ESVINT,  float *, _ESVINT);
void   esvdpps(const double *, _ESVINT, double *, _ESVINT);
 
void  esvspptrs(const char *, _ESVINT, _ESVINT, const float  *, void *, _ESVINT, _ESVI);
void  esvdpptrs(const char *, _ESVINT, _ESVINT, const double *, void *, _ESVINT, _ESVI);
void  esvcpptrs(const char *, _ESVINT, _ESVINT, const _ESVCM *, void *, _ESVINT, _ESVI);
void  esvzpptrs(const char *, _ESVINT, _ESVINT, const _ESVCOM *, void *, _ESVINT, _ESVI);
 
void  esvsppsv(const char *, _ESVINT, _ESVINT, float  *, void *, _ESVINT, _ESVI);
void  esvdppsv(const char *, _ESVINT, _ESVINT, double *, void *, _ESVINT, _ESVI);
void  esvcppsv(const char *, _ESVINT, _ESVINT, _ESVCM *, void *, _ESVINT, _ESVI);
void  esvzppsv(const char *, _ESVINT, _ESVINT, _ESVCOM *, void *, _ESVINT, _ESVI);
 
void  esvsposm(const char *, const void *,  _ESVINT,  _ESVINT, void *,  _ESVINT,  _ESVINT);
void  esvdposm(const char *, const void *,  _ESVINT,  _ESVINT, void *,  _ESVINT,  _ESVINT);
void  esvcposm(const char *, const void *,  _ESVINT,  _ESVINT, void *,  _ESVINT,  _ESVINT);
void  esvzposm(const char *, const void *,  _ESVINT,  _ESVINT, void *,  _ESVINT,  _ESVINT);
 
void  esvspotrs(const char *, _ESVINT, _ESVINT, const void *, _ESVINT, void *, _ESVINT, _ESVI);
void  esvdpotrs(const char *, _ESVINT, _ESVINT, const void *, _ESVINT, void *, _ESVINT, _ESVI);
void  esvcpotrs(const char *, _ESVINT, _ESVINT, const void *, _ESVINT, void *, _ESVINT, _ESVI);
void  esvzpotrs(const char *, _ESVINT, _ESVINT, const void *, _ESVINT, void *, _ESVINT, _ESVI);
 
void  esvsposv(const char *, _ESVINT, _ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVI);
void  esvdposv(const char *, _ESVINT, _ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVI);
void  esvcposv(const char *, _ESVINT, _ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVI);
void  esvzposv(const char *, _ESVINT, _ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVI);
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvsppfcd_er( float *, _ESVINT, _ESVINT, _ESVS,  float *,  float *, _ESVI);
 _ESVINT    esvdppfcd_er(double *, _ESVINT, _ESVINT, _ESVD, double *, double *, _ESVI);
 
 _ESVINT    esvspofcd_er(const char *, void*, _ESVINT, _ESVINT, _ESVINT,
                     _ESVS,  float *,  float *, _ESVI);
 _ESVINT    esvdpofcd_er(const char *, void*, _ESVINT, _ESVINT, _ESVINT,
                     _ESVD, double *, double *, _ESVI);
#else
 _ESVINT    esvsppfcd   ( float *, _ESVINT, _ESVINT, _ESVS,  float *,  float *, _ESVINT);
 _ESVINT    esvdppfcd   (double *, _ESVINT, _ESVINT, _ESVD, double *, double *, _ESVINT);
 
 _ESVINT    esvspofcd   (const char *, void *, _ESVINT, _ESVINT, _ESVINT,
                     _ESVS,  float*,  float *, _ESVINT);
 _ESVINT    esvdpofcd   (const char *, void *, _ESVINT, _ESVINT, _ESVINT,
                     _ESVD, double*, double *, _ESVINT);
#endif
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvdbssv_er(const char *, _ESVINT, _ESVINT, double *, _ESVINT *, void *, _ESVINT,
                    _ESVI);
 
 _ESVINT    esvdbstrf_er(const char *, _ESVINT, double *, _ESVINT *, _ESVI);
#else
 void   esvdbssv   (const char *, _ESVINT, _ESVINT, double *, _ESVINT *, void *, _ESVINT,
                    _ESVI);
 
 void   esvdbstrf   (const char *, _ESVINT, double *, _ESVINT *, _ESVI);
#endif
 
void   esvdbstrs(const char *, _ESVINT, _ESVINT, const double *, const _ESVINT *,
                 void *, _ESVINT, _ESVI);
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvsgeicd_er(void *, _ESVINT, _ESVINT, _ESVINT, _ESVS,  float *,  float *, _ESVI);
 _ESVINT    esvdgeicd_er(void *, _ESVINT, _ESVINT, _ESVINT, _ESVD, double *, double *, _ESVI);
#else
 _ESVINT    esvsgeicd   (void *, _ESVINT, _ESVINT, _ESVINT, _ESVS,  float *,  float *, _ESVINT);
 _ESVINT    esvdgeicd   (void *, _ESVINT, _ESVINT, _ESVINT, _ESVD, double *, double *, _ESVINT);
#endif
 
void   esvsgetri(_ESVINT, void *, _ESVINT, const _ESVINT *, float   *, _ESVINT, _ESVI);
void   esvdgetri(_ESVINT, void *, _ESVINT, const _ESVINT *, double  *, _ESVINT, _ESVI);
void   esvcgetri(_ESVINT, void *, _ESVINT, const _ESVINT *, _ESVCM  *, _ESVINT, _ESVI);
void   esvzgetri(_ESVINT, void *, _ESVINT, const _ESVINT *, _ESVCOM *, _ESVINT, _ESVI);
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvsppicd_er( float *, _ESVINT, _ESVINT, _ESVS,  float *,  float *, _ESVI);
 _ESVINT    esvdppicd_er(double *, _ESVINT, _ESVINT, _ESVD, double *, double *, _ESVI);
 
 _ESVINT    esvspoicd_er(const char *, void *, _ESVINT, _ESVINT, _ESVINT,
                     _ESVS,  float *,  float *, _ESVI);
 _ESVINT    esvdpoicd_er(const char *, void *, _ESVINT, _ESVINT, _ESVINT,
                     _ESVD, double *, double *, _ESVI);
#else
 _ESVINT    esvsppicd   ( float *, _ESVINT, _ESVINT, _ESVS,  float *,  float *, _ESVINT);
 _ESVINT    esvdppicd   (double *, _ESVINT, _ESVINT, _ESVD, double *, double *, _ESVINT);
 
 _ESVINT    esvspoicd   (const char *, void *, _ESVINT, _ESVINT, _ESVINT,
                     _ESVS,  float *,  float *, _ESVINT);
 _ESVINT    esvdpoicd   (const char *, void *, _ESVINT, _ESVINT, _ESVINT,
                     _ESVD, double *, double *, _ESVINT);
#endif
 
void   esvspotri(const char *, _ESVINT, void *, _ESVINT, _ESVI);
void   esvdpotri(const char *, _ESVINT, void *, _ESVINT, _ESVI);
void   esvcpotri(const char *, _ESVINT, void *, _ESVINT, _ESVI);
void   esvzpotri(const char *, _ESVINT, void *, _ESVINT, _ESVI);
 
void   esvstrsv(const char *, const char *, const char *, _ESVINT,
                const void *, _ESVINT, void *, _ESVINT);
void   esvdtrsv(const char *, const char *, const char *, _ESVINT,
                const void *, _ESVINT, void *, _ESVINT);
void   esvctrsv(const char *, const char *, const char *, _ESVINT,
                const void *, _ESVINT, void *, _ESVINT);
void   esvztrsv(const char *, const char *, const char *, _ESVINT,
                const void *, _ESVINT, void *, _ESVINT);
 
void   esvstpsv(const char *, const char *, const char *, _ESVINT,
                const   float *,   float *, _ESVINT);
void   esvdtpsv(const char *, const char *, const char *, _ESVINT,
                const  double *,  double *, _ESVINT);
void   esvctpsv(const char *, const char *, const char *, _ESVINT,
                const  _ESVCM *,  _ESVCM *, _ESVINT);
void   esvztpsv(const char *, const char *, const char *, _ESVINT,
                const _ESVCOM *, _ESVCOM *, _ESVINT);
 
void   esvstrsm(const char *, const char *, const char *, const char *,
                _ESVINT, _ESVINT,   float, const void *, _ESVINT, void *, _ESVINT);
void   esvdtrsm(const char *, const char *, const char *, const char *,
                _ESVINT, _ESVINT,  double, const void *, _ESVINT, void *, _ESVINT);
void   esvctrsm(const char *, const char *, const char *, const char *,
                _ESVINT, _ESVINT,  _ESVCM, const void *, _ESVINT, void *, _ESVINT);
void   esvztrsm(const char *, const char *, const char *, const char *,
                _ESVINT, _ESVINT, _ESVCOM, const void *, _ESVINT, void *, _ESVINT);
 
_ESVINT    esvstri(const char *, const char *, void *, _ESVINT, _ESVINT);
_ESVINT    esvdtri(const char *, const char *, void *, _ESVINT, _ESVINT);
 
void   esvstrtri(const char *, const char *, _ESVINT, void *, _ESVINT, _ESVI);
void   esvdtrtri(const char *, const char *, _ESVINT, void *, _ESVINT, _ESVI);
void   esvctrtri(const char *, const char *, _ESVINT, void *, _ESVINT, _ESVI);
void   esvztrtri(const char *, const char *, _ESVINT, void *, _ESVINT, _ESVI);
 
_ESVINT    esvstpi(const char *, const char *,  float *, _ESVINT);
_ESVINT    esvdtpi(const char *, const char *, double *, _ESVINT);
 
void   esvstptri(const char *, const char *, _ESVINT,   float *, _ESVI);
void   esvdtptri(const char *, const char *, _ESVINT,  double *, _ESVI);
void   esvctptri(const char *, const char *, _ESVINT,  _ESVCM *, _ESVI);
void   esvztptri(const char *, const char *, _ESVINT, _ESVCOM *, _ESVI);
 
/*  Banded Linear Algebraic Equation Subroutines  */
 
_ESVINT    esvsgbf(void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT *);
_ESVINT    esvdgbf(void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT *);
 
void   esvsgbs(const void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, const _ESVINT *,  float *);
void   esvdgbs(const void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, const _ESVINT *, double *);
 
_ESVINT    esvspbf(void *, _ESVINT, _ESVINT, _ESVINT);
_ESVINT    esvdpbf(void *, _ESVINT, _ESVINT, _ESVINT);
 
_ESVINT    esvspbchf(void *, _ESVINT, _ESVINT, _ESVINT);
_ESVINT    esvdpbchf(void *, _ESVINT, _ESVINT, _ESVINT);
 
void   esvspbs(const void *, _ESVINT, _ESVINT, _ESVINT,  float *);
void   esvdpbs(const void *, _ESVINT, _ESVINT, _ESVINT, double *);
 
void   esvspbchs(const void *, _ESVINT, _ESVINT, _ESVINT,  float *);
void   esvdpbchs(const void *, _ESVINT, _ESVINT, _ESVINT, double *);
 
_ESVINT    esvsgtf(_ESVINT,  float *,  float *,  float *,  float *, _ESVINT *);
_ESVINT    esvdgtf(_ESVINT, double *, double *, double *, double *, _ESVINT *);
 
void   esvsgts(_ESVINT, const  float *, const  float *, const  float *,
               const  float *, const _ESVINT *,  float *);
void   esvdgts(_ESVINT, const double *, const double *, const double *,
               const double *, const _ESVINT *, double *);
 
void   esvsgtnp(_ESVINT,   float *,   float *,   float *,   float *);
void   esvdgtnp(_ESVINT,  double *,  double *,  double *,  double *);
void   esvcgtnp(_ESVINT,  _ESVCM *,  _ESVCM *,  _ESVCM *,  _ESVCM *);
void   esvzgtnp(_ESVINT, _ESVCOM *, _ESVCOM *, _ESVCOM *, _ESVCOM *);
 
void   esvsgtnpf(_ESVINT,   float *,   float *,   float *, _ESVINT);
void   esvdgtnpf(_ESVINT,  double *,  double *,  double *, _ESVINT);
void   esvcgtnpf(_ESVINT,  _ESVCM *,  _ESVCM *,  _ESVCM *, _ESVINT);
void   esvzgtnpf(_ESVINT, _ESVCOM *, _ESVCOM *, _ESVCOM *, _ESVINT);
 
void   esvsgtnps(_ESVINT, const   float *, const   float *,
                 const   float *,   float *);
void   esvdgtnps(_ESVINT, const  double *, const  double *,
                 const  double *,  double *);
void   esvcgtnps(_ESVINT, const  _ESVCM *, const  _ESVCM *,
                 const  _ESVCM *,  _ESVCM *);
void   esvzgtnps(_ESVINT, const _ESVCOM *, const _ESVCOM *,
                 const _ESVCOM *, _ESVCOM *);
 
void   esvsptf(_ESVINT,  float *,  float *,  _ESVINT);
void   esvdptf(_ESVINT, double *, double *,  _ESVINT);
 
void   esvspts(_ESVINT, const  float *, const  float *,  float *);
void   esvdpts(_ESVINT, const double *, const double *, double *);
 
void   esvstbsv(const char *, const char *, const char *, _ESVINT, _ESVINT,
                const void *, _ESVINT,   float *, _ESVINT);
void   esvdtbsv(const char *, const char *, const char *, _ESVINT, _ESVINT,
                const void *, _ESVINT,  double *, _ESVINT);
void   esvctbsv(const char *, const char *, const char *, _ESVINT, _ESVINT,
                const void *, _ESVINT,  _ESVCM *, _ESVINT);
void   esvztbsv(const char *, const char *, const char *, _ESVINT, _ESVINT,
                const void *, _ESVINT, _ESVCOM *, _ESVINT);
 
/*  Sparse Linear Algebraic Equation Subroutines  */
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvdgsf_er(_ESVINT, _ESVINT, _ESVINT, double *, _ESVINT *, _ESVINT *, _ESVINT, const _ESVINT *,
                   const double *, double *, double *, _ESVI);
 _ESVINT    esvdgss_er(_ESVINT, _ESVINT, const double *, _ESVINT *, const _ESVINT *,
                   _ESVINT , double *, double *, _ESVI);
 _ESVINT    esvdgkfs_er (_ESVINT, double *, _ESVINT, _ESVINT *, double *, _ESVINT, _ESVINT *,
                     _ESVINT *, double *, double *, _ESVI, void *, _ESVINT, _ESVINT);
 _ESVINT    esvdgkfsp_er(_ESVINT, double *, _ESVINT, _ESVINT *, double *, _ESVINT, _ESVINT *,
                     _ESVINT *, double *, double *, _ESVI, void *, _ESVINT, _ESVINT);
 _ESVINT    esvdskfs_er (_ESVINT, double *, _ESVINT, _ESVINT *, _ESVINT *, double *,
                     double *, _ESVI, void *, _ESVINT, _ESVINT);
 _ESVINT    esvdskfsp_er(_ESVINT, double *, _ESVINT, _ESVINT *, _ESVINT *, double *,
                     double *, _ESVI, void *, _ESVINT, _ESVINT);
 _ESVINT    esvdsris_er(const char *, const char *, _ESVINT, double *, _ESVINT *,
                    const _ESVINT *, const double *, double *, _ESVINT *,
                    double *, double *, _ESVI, double *, _ESVI);
 _ESVINT    esvdsmcg_er(_ESVINT, _ESVINT, const void *, const void *, _ESVINT,
                    const double *, double *, _ESVINT *, double *,
                    double *, _ESVI, double *, _ESVI);
 _ESVINT    esvdsdcg_er(_ESVINT, _ESVINT, _ESVINT, const void *, _ESVINT, const _ESVINT *,
                    const double *, double *, _ESVINT *, double *,
                    double *, _ESVI, double *, _ESVI);
 _ESVINT    esvdsmgcg_er(_ESVINT, _ESVINT, const void *, const void *, _ESVINT,
                     const double *, double *, _ESVINT *, double *,
                     double *, _ESVI, double *, _ESVI);
 _ESVINT    esvdsdgcg_er(_ESVINT, _ESVINT, const void *, _ESVINT, const _ESVINT *,
                     const double *, double *, _ESVINT *, double *,
                     double *, _ESVI, double *, _ESVI);
#else
 _ESVINT    esvdgsf   (_ESVINT, _ESVINT, _ESVINT, double *, _ESVINT *, _ESVINT *, _ESVINT, const _ESVINT *,
                   const double *, double *, double *,  _ESVINT);
 void   esvdgss   (_ESVINT, _ESVINT, const double *, _ESVINT *, const _ESVINT *,
                   _ESVINT, double *, double *,  _ESVINT);
 _ESVINT    esvdgkfs    (_ESVINT, double *, _ESVINT, _ESVINT *, double *, _ESVINT, _ESVINT *,
                     _ESVINT *, double *, double *, _ESVINT, void *, _ESVINT, _ESVINT);
 _ESVINT    esvdgkfsp   (_ESVINT, double *, _ESVINT, _ESVINT *, double *, _ESVINT, _ESVINT *,
                     _ESVINT *, double *, double *, _ESVINT, void *, _ESVINT, _ESVINT);
 _ESVINT    esvdskfs    (_ESVINT, double *, _ESVINT, _ESVINT *, _ESVINT *, double *,
                     double *, _ESVINT, void *, _ESVINT, _ESVINT);
 _ESVINT    esvdskfsp   (_ESVINT, double *, _ESVINT, _ESVINT *, _ESVINT *, double *,
                     double *, _ESVINT, void *, _ESVINT, _ESVINT);
 _ESVINT    esvdsris   (const char *, const char *, _ESVINT, double *, _ESVINT *,
                    const _ESVINT *, const double *, double *, _ESVINT *,
                    double *, double *, _ESVINT, double *, _ESVINT);
 _ESVINT    esvdsmcg   (_ESVINT, _ESVINT, const void *, const void *, _ESVINT,
                    const double *, double *, _ESVINT *, double *,
                    double *, _ESVINT, double *, _ESVINT);
 _ESVINT    esvdsdcg   (_ESVINT, _ESVINT, _ESVINT, const void *, _ESVINT, const _ESVINT *,
                    const double *, double *, _ESVINT *, double *,
                    double *, _ESVINT, double *, _ESVINT);
 _ESVINT    esvdsmgcg   (_ESVINT, _ESVINT, const void *, const void *, _ESVINT,
                     const double *, double *, _ESVINT *, double *,
                     double *, _ESVINT, double *, _ESVINT);
 _ESVINT    esvdsdgcg   (_ESVINT, _ESVINT, const void *, _ESVINT, const _ESVINT *,
                     const double *, double *, _ESVINT *, double *,
                     double *, _ESVINT, double *, _ESVINT);
#endif
 
/*  Linear Least Squares Subroutines  */
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvsgesvf_er(_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVINT,  float *,
                     _ESVINT, _ESVINT,  float *, _ESVI);
 _ESVINT    esvdgesvf_er(_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVINT, double *,
                     _ESVINT, _ESVINT, double *, _ESVI);
#else
 _ESVINT    esvsgesvf   (_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVINT,  float *,
                     _ESVINT, _ESVINT,  float *, _ESVINT);
 _ESVINT    esvdgesvf   (_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVINT, double *,
                     _ESVINT, _ESVINT, double *, _ESVINT);
#endif
 
void   esvsgesvs(const void *, _ESVINT, void *, _ESVINT, _ESVINT, const  float *,
                 void *, _ESVINT, _ESVINT, _ESVINT,  float);
void   esvdgesvs(const void *, _ESVINT, void *, _ESVINT, _ESVINT, const double *,
                 void *, _ESVINT, _ESVINT, _ESVINT, double);
 
void   esvdgeqrf(_ESVINT, _ESVINT, void *, _ESVINT, double  *, double  *, _ESVINT, _ESVI);
void   esvsgeqrf(_ESVINT, _ESVINT, void *, _ESVINT, float   *, float   *, _ESVINT, _ESVI);
void   esvzgeqrf(_ESVINT, _ESVINT, void *, _ESVINT, _ESVCOM *, _ESVCOM *, _ESVINT, _ESVI);
void   esvcgeqrf(_ESVINT, _ESVINT, void *, _ESVINT, _ESVCM  *, _ESVCM  *, _ESVINT, _ESVI);
 
void   esvdgels(const char *, _ESVINT, _ESVINT, _ESVINT, void *, _ESVINT,
                void *, _ESVINT, double *, _ESVINT, _ESVI);
void   esvsgels(const char *, _ESVINT, _ESVINT, _ESVINT, void *, _ESVINT,
                void *, _ESVINT, float *, _ESVINT, _ESVI);
void   esvzgels(const char *, _ESVINT, _ESVINT, _ESVINT, void *, _ESVINT,
                void *, _ESVINT, _ESVCOM *, _ESVINT, _ESVI);
void   esvcgels(const char *, _ESVINT, _ESVINT, _ESVINT, void *, _ESVINT,
                void *, _ESVINT, _ESVCM *, _ESVINT, _ESVI);
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvsgells_er(_ESVINT, void *, _ESVINT, void *, _ESVINT, void *, _ESVINT,
                      float *,  float, _ESVINT, _ESVINT, _ESVINT, _ESVI,  float *, _ESVI);
 _ESVINT    esvdgells_er(_ESVINT, void *, _ESVINT, void *, _ESVINT, void *, _ESVINT,
                     double *, double, _ESVINT, _ESVINT, _ESVINT, _ESVI, double *, _ESVI);
#else
 void   esvsgells   (_ESVINT, void *, _ESVINT, void *, _ESVINT, void *, _ESVINT,
                      float *,  float, _ESVINT, _ESVINT, _ESVINT, _ESVI,  float *, _ESVINT);
 void   esvdgells   (_ESVINT, void *, _ESVINT, void *, _ESVINT, void *, _ESVINT,
                     double *, double, _ESVINT, _ESVINT, _ESVINT, _ESVI, double *, _ESVINT);
#endif
 
/*  Eigensystem Analysis Subroutines  */
 
void   esvssyevx(const char *, const char *, const char *, _ESVINT, void *,
                 _ESVINT, float, float, _ESVINT, _ESVINT, float, _ESVINT *,
                 float *, void *, _ESVINT, float *, _ESVINT, _ESVINT *,
                 _ESVINT *, _ESVI);
 
void   esvdsyevx(const char *, const char *, const char *, _ESVINT, void *,
                 _ESVINT, double, double, _ESVINT, _ESVINT, double, _ESVINT *,
                 double *, void *, _ESVINT, double *, _ESVINT, _ESVINT *,
                 _ESVINT *, _ESVI);
 
void   esvcheevx(const char *, const char *, const char *, _ESVINT, void *,
                 _ESVINT, float, float, _ESVINT, _ESVINT, float, _ESVINT *,
                 float *, void *, _ESVINT, _ESVCM *, _ESVINT, float *,
                  _ESVINT *, _ESVINT *, _ESVI);
 
void   esvzheevx(const char *, const char *, const char *, _ESVINT, void *,
                 _ESVINT, double, double, _ESVINT, _ESVINT, double, _ESVINT *,
                 double *, void *, _ESVINT, _ESVCOM *, _ESVINT, double *,
                  _ESVINT *, _ESVINT *, _ESVI);
 
void   esvsgeevx(const char *, const char *, const char *, const char *, _ESVINT,
                 void *, _ESVINT, float *, float *, void *, _ESVINT, void *,
                 _ESVINT, _ESVINT *, _ESVINT *, float *, float *, float *,
                 float *, float *, _ESVINT, _ESVINT *, _ESVI );
 
void   esvdgeevx(const char *, const char *, const char *, const char *, _ESVINT,
                 void *, _ESVINT, double *, double *, void *, _ESVINT, void *,
                 _ESVINT, _ESVINT *, _ESVINT *, double *, double *, double *,
                 double *, double *, _ESVINT, _ESVINT *, _ESVI );
 
void   esvcgeevx(const char *, const char *, const char *, const char *, _ESVINT,
                 void *, _ESVINT, _ESVCM *, void *, _ESVINT, void *, _ESVINT,
                 _ESVINT*, _ESVINT*, float *, float*, float *, float *, _ESVCM *,
                 _ESVINT, float *, _ESVI );
 
void   esvzgeevx(const char *, const char *, const char *, const char *, _ESVINT,
                 void *, _ESVINT, _ESVCOM *, void *, _ESVINT, void *, _ESVINT,
                 _ESVINT*, _ESVINT*, double *, double*, double *, double *,
                 _ESVCOM *, _ESVINT, double *, _ESVI );
 
void   esvsspevx(const char *, const char *, const char *, _ESVINT, void *, float,
                 float, _ESVINT, _ESVINT, float, _ESVINT *, float *, void *,
                 _ESVINT, float *, _ESVINT *, _ESVINT *, _ESVI);
 
void   esvdspevx(const char *, const char *, const char *, _ESVINT, void *, double,
                 double, _ESVINT, _ESVINT, double, _ESVINT *, double *, void *,
                 _ESVINT, double *, _ESVINT *, _ESVINT *, _ESVI);
 
void   esvchpevx(const char *, const char *, const char *, _ESVINT, void *, float,
                 float, _ESVINT, _ESVINT, float, _ESVINT *, float *, void *,
                 _ESVINT, _ESVCM *, float *, _ESVINT *, _ESVINT *, _ESVI);
 
void   esvzhpevx(const char *, const char *, const char *, _ESVINT, void *, double,
                 double, _ESVINT, _ESVINT, double, _ESVINT *, double *, void *,
                 _ESVINT, _ESVCOM *, double *, _ESVINT *, _ESVINT *, _ESVI);
 
void   esvdsygvx  (_ESVINT, const char *, const char *, const char *,
                 _ESVINT, void *, _ESVINT, void *, _ESVINT, double,
                 double, _ESVINT, _ESVINT, double, _ESVINT *, double *,
                 void *, _ESVINT, double *, _ESVINT, _ESVINT *, _ESVINT *,
                 _ESVI);
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvsgeev_er(_ESVINT, void *, _ESVINT,  _ESVCM *, void *, _ESVINT,
                    const _ESVINT *, _ESVINT,  float *, _ESVI);
 _ESVINT    esvdgeev_er(_ESVINT, void *, _ESVINT, _ESVCOM *, void *, _ESVINT,
                    const _ESVINT *, _ESVINT, double *, _ESVI);
 _ESVINT    esvcgeev_er(_ESVINT, void *, _ESVINT,  _ESVCM *, void *, _ESVINT,
                    const _ESVINT *, _ESVINT,  float *, _ESVI);
 _ESVINT    esvzgeev_er(_ESVINT, void *, _ESVINT, _ESVCOM *, void *, _ESVINT,
                    const _ESVINT *, _ESVINT, double *, _ESVI);
#else
 _ESVINT    esvsgeev   (_ESVINT, void *, _ESVINT,  _ESVCM *, void *, _ESVINT,
                    const _ESVINT *, _ESVINT,  float *, _ESVINT);
 _ESVINT    esvdgeev   (_ESVINT, void *, _ESVINT, _ESVCOM *, void *, _ESVINT,
                    const _ESVINT *, _ESVINT, double *, _ESVINT);
 _ESVINT    esvcgeev   (_ESVINT, void *, _ESVINT,  _ESVCM *, void *, _ESVINT,
                    const _ESVINT *, _ESVINT,  float *, _ESVINT);
 _ESVINT    esvzgeev   (_ESVINT, void *, _ESVINT, _ESVCOM *, void *, _ESVINT,
                    const _ESVINT *, _ESVINT, double *, _ESVINT);
#endif
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvsspev_er(_ESVINT,   float *,  float *, void *, _ESVINT,
                    _ESVINT,  float *, _ESVI);
 _ESVINT    esvdspev_er(_ESVINT,  double *, double *, void *, _ESVINT,
                    _ESVINT, double *, _ESVI);
 _ESVINT    esvchpev_er(_ESVINT,  _ESVCM *,  float *, void *, _ESVINT,
                    _ESVINT,  float *, _ESVI);
 _ESVINT    esvzhpev_er(_ESVINT, _ESVCOM *, double *, void *, _ESVINT,
                    _ESVINT, double *, _ESVI);
#else
 _ESVINT    esvsspev   (_ESVINT,   float *,  float *, void *, _ESVINT,
                    _ESVINT,  float *, _ESVINT);
 _ESVINT    esvdspev   ( _ESVINT, double *, double *, void *, _ESVINT,
                    _ESVINT, double *, _ESVINT);
 _ESVINT    esvchpev   (_ESVINT,  _ESVCM *,  float *, void *, _ESVINT,
                    _ESVINT,  float *, _ESVINT);
 _ESVINT    esvzhpev   (_ESVINT, _ESVCOM *, double *, void *, _ESVINT,
                    _ESVINT, double *, _ESVINT);
#endif
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvsspsv_er(_ESVINT,   float *,  float *, void *, _ESVINT, _ESVINT,
                    _ESVINT,  float *, _ESVI);
 _ESVINT    esvdspsv_er(_ESVINT,  double *, double *, void *, _ESVINT, _ESVINT,
                    _ESVINT, double *, _ESVI);
 _ESVINT    esvchpsv_er(_ESVINT,  _ESVCM *,  float *, void *, _ESVINT, _ESVINT,
                    _ESVINT,  float *, _ESVI);
 _ESVINT    esvzhpsv_er(_ESVINT, _ESVCOM *, double *, void *, _ESVINT, _ESVINT,
                    _ESVINT, double *, _ESVI);
#else
 _ESVINT    esvsspsv   (_ESVINT,   float *,  float *, void *, _ESVINT, _ESVINT,
                    _ESVINT,  float*, _ESVINT);
 _ESVINT    esvdspsv   (_ESVINT,  double *, double *, void *, _ESVINT, _ESVINT,
                    _ESVINT, double*, _ESVINT);
 _ESVINT    esvchpsv   (_ESVINT,  _ESVCM *,  float *, void *, _ESVINT, _ESVINT,
                    _ESVINT,  float*, _ESVINT);
 _ESVINT    esvzhpsv   (_ESVINT, _ESVCOM *, double *, void *, _ESVINT, _ESVINT,
                    _ESVINT, double*, _ESVINT);
#endif
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvsgegv_er(_ESVINT, void *, _ESVINT, void *, _ESVINT,  _ESVCM *,  float *,
                    void *, _ESVINT, _ESVINT,  float *, _ESVI);
 _ESVINT    esvdgegv_er(_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVCOM *, double *,
                    void *, _ESVINT, _ESVINT, double *, _ESVI);
#else
 _ESVINT    esvsgegv   (_ESVINT, void *, _ESVINT, void *, _ESVINT,  _ESVCM *,  float *,
                    void *, _ESVINT, _ESVINT,  float *, _ESVINT);
 _ESVINT    esvdgegv   (_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVCOM *, double *,
                    void *, _ESVINT, _ESVINT, double *, _ESVINT);
#endif
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvssygv_er(_ESVINT, void *, _ESVINT, void *, _ESVINT,   float *,
                    void *, _ESVINT, _ESVINT,  float *, _ESVI);
 _ESVINT    esvdsygv_er(_ESVINT, void *, _ESVINT, void *, _ESVINT,  double *,
                    void *, _ESVINT, _ESVINT, double *, _ESVI);
#else
 _ESVINT    esvssygv   (_ESVINT, void *, _ESVINT, void *, _ESVINT,   float *,
                    void *, _ESVINT, _ESVINT,  float *, _ESVINT);
 _ESVINT    esvdsygv   (_ESVINT, void *, _ESVINT, void *, _ESVINT,  double *,
                    void *, _ESVINT, _ESVINT, double *, _ESVINT);
#endif
 
/*  Fourier Transform Subroutines  */
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvscftd_er( _ESVINT, _ESVINT, void *, _ESVINT *, _ESVINT, void *, _ESVINT *, _ESVINT,
                    _ESVINT *, _ESVINT, _ESVINT *, float, double *, _ESVI, double *, _ESVI);
 _ESVINT    esvdcftd_er( _ESVINT, _ESVINT, void *, _ESVINT *, _ESVINT, void *, _ESVINT *, _ESVINT,
                    _ESVINT *, _ESVINT, _ESVINT *, double, double *, _ESVI, double *, _ESVI);
 
 _ESVINT    esvsrcftd_er( _ESVINT, _ESVINT, void *, _ESVINT *, _ESVINT, void *, _ESVINT *, _ESVINT,
                     _ESVINT *, _ESVINT, _ESVINT *, float, double *, _ESVI, double *, _ESVI);
 _ESVINT    esvdrcftd_er( _ESVINT, _ESVINT, void *, _ESVINT *, _ESVINT, void *, _ESVINT *, _ESVINT,
                     _ESVINT *, _ESVINT, _ESVINT *, double, double *, _ESVI, double *, _ESVI);
 
 _ESVINT    esvscrftd_er( _ESVINT, _ESVINT, void *, _ESVINT *, _ESVINT, void *, _ESVINT *, _ESVINT,
                     _ESVINT *, _ESVINT, _ESVINT *, float, double *, _ESVI, double *, _ESVI);
 _ESVINT    esvdcrftd_er( _ESVINT, _ESVINT, void *, _ESVINT *, _ESVINT, void *, _ESVINT *, _ESVINT,
                     _ESVINT *, _ESVINT, _ESVINT *, double, double *, _ESVI, double *, _ESVI);
 
 _ESVINT    esvscft_er(_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVI,
                   _ESVINT, _ESVINT,  float, double *, _ESVI, double *, _ESVI);
 _ESVINT    esvdcft_er(_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVI,
                   _ESVINT, _ESVINT, double, double *, _ESVI, double *, _ESVI);
 
 _ESVINT    esvscftp_er(_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT *, _ESVINT,
                    _ESVINT, float, double *, _ESVI, double *, _ESVI);
 
 _ESVINT    esvsrcft_er(_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVI, _ESVINT, _ESVINT,  float,
                    double *, _ESVI, double *, _ESVI, double *, _ESVI);
 _ESVINT    esvdrcft_er(_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVI, _ESVINT, _ESVINT, double,
                    double *, _ESVI, double *, _ESVI);
 
 _ESVINT    esvscrft_er(_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVI, _ESVINT, _ESVINT,  float,
                    double *, _ESVI, double *, _ESVI, double *, _ESVI);
 _ESVINT    esvdcrft_er(_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVI, _ESVINT, _ESVINT, double,
                    double *, _ESVI, double *, _ESVI);
#else
 void   esvscftd  ( _ESVINT, _ESVINT, void *, _ESVINT *, _ESVINT, void *, _ESVINT *, _ESVINT,
                    _ESVINT *, _ESVINT, _ESVINT *, float, double *, _ESVINT, double *, _ESVINT);
 void   esvdcftd  ( _ESVINT, _ESVINT, void *, _ESVINT *, _ESVINT, void *, _ESVINT *, _ESVINT,
                    _ESVINT *, _ESVINT, _ESVINT *, double, double *, _ESVINT, double *, _ESVINT);
 
 void   esvsrcftd  ( _ESVINT, _ESVINT, void *, _ESVINT *, _ESVINT, void *, _ESVINT *, _ESVINT,
                     _ESVINT *, _ESVINT, _ESVINT *, float, double *, _ESVINT, double *, _ESVINT);
 void   esvdrcftd  ( _ESVINT, _ESVINT, void *, _ESVINT *, _ESVINT, void *, _ESVINT *, _ESVINT,
                     _ESVINT *, _ESVINT, _ESVINT *, double, double *, _ESVINT, double *, _ESVINT);
 
 void   esvscrftd  ( _ESVINT, _ESVINT, void *, _ESVINT *, _ESVINT, void *, _ESVINT *, _ESVINT,
                     _ESVINT *, _ESVINT, _ESVINT *, float, double *, _ESVINT, double *, _ESVINT);
 void   esvdcrftd  ( _ESVINT, _ESVINT, void *, _ESVINT *, _ESVINT, void *, _ESVINT *, _ESVINT,
                     _ESVINT *, _ESVINT, _ESVINT *, double, double *, _ESVINT, double *, _ESVINT);
 
 void   esvscft   (_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT,
                   _ESVINT, _ESVINT,  float, double *, _ESVINT, double *, _ESVINT);
 void   esvdcft   (_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT,
                   _ESVINT, _ESVINT, double, double *, _ESVINT, double *, _ESVINT);
 void   esvscftp   (_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,
                    _ESVINT, float, double *, _ESVINT, double *, _ESVINT);
 
 void   esvsrcft   (_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,  float,
                    double *, _ESVINT, double *, _ESVINT, double *, _ESVINT);
 void   esvdrcft   (_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, double,
                    double *, _ESVINT, double *, _ESVINT);
 
 void   esvscrft   (_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,  float,
                    double *, _ESVINT, double *, _ESVINT, double *, _ESVINT);
 void   esvdcrft   (_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, double,
                    double *, _ESVINT, double *, _ESVINT);
#endif
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvscosf_er(_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVI,
                    _ESVINT, float, double *, _ESVI, double *, _ESVI);
 _ESVINT    esvdcosf_er(_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVI,
                    _ESVINT, double, double *, _ESVI, double *, _ESVI);
 _ESVINT    esvscosft_er(_ESVINT, void *,  _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVI,
                     _ESVINT, float, double *, _ESVI, double *, _ESVI);
 
 _ESVINT    esvssinf_er(_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVI,
                    _ESVINT,  float, double *, _ESVI, double *, _ESVI);
 _ESVINT    esvdsinf_er(_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVI,
                    _ESVINT, double, double *, _ESVI, double *, _ESVI);
#else
 void   esvscosf   (_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT,
                    _ESVINT,  float, double *, _ESVINT, double *, _ESVINT);
 void   esvdcosf   (_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT,
                    _ESVINT, double, double *, _ESVINT, double *, _ESVINT);
 void   esvscosft   (_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT,
                     _ESVINT,  float, double *, _ESVINT, double *, _ESVINT);
 
 void   esvssinf   (_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT,
                    _ESVINT,  float, double *, _ESVINT, double *, _ESVINT);
 void   esvdsinf   (_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT,
                    _ESVINT, double, double *, _ESVINT, double *, _ESVINT);
#endif
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvscft2_er(_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVI,
                    _ESVI, _ESVINT,  float, double *, _ESVI, double *, _ESVI);
 _ESVINT    esvdcft2_er(_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVI,
                    _ESVI, _ESVINT, double, double *, _ESVI, double *, _ESVI);
 _ESVINT    esvscft2p_er(_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT *,
                     _ESVINT *, _ESVINT, float, double *, _ESVI, double *, _ESVI);
 
 _ESVINT    esvsrcft2_er(_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVI, _ESVI,
                     _ESVINT,  float, double *, _ESVI, double *, _ESVI,
                     double *, _ESVI);
 _ESVINT    esvdrcft2_er(_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVI, _ESVI,
                     _ESVINT, double, double *, _ESVI, double *, _ESVI);
 
 _ESVINT    esvscrft2_er(_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVI, _ESVI,
                     _ESVINT,  float, double *, _ESVI, double *, _ESVI,
                     double *, _ESVI);
 _ESVINT    esvdcrft2_er(_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVI, _ESVI,
                     _ESVINT, double, double *, _ESVI, double *, _ESVI);
 
 _ESVINT    esvscft3_er(void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVI, _ESVI,
                    _ESVI, _ESVINT,  float, double *, _ESVI);
 _ESVINT    esvdcft3_er(void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVI, _ESVI,
                    _ESVI, _ESVINT, double, double *, _ESVI);
 _ESVINT    esvscft3p_er(void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT *, _ESVINT *,
                     _ESVINT *, _ESVINT, float, double *, _ESVI);
 
 _ESVINT    esvsrcft3_er(void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVI, _ESVI,
                     _ESVI, _ESVINT,  float, double *, _ESVI);
 _ESVINT    esvdrcft3_er(void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVI, _ESVI,
                     _ESVI, _ESVINT, double, double *, _ESVI);
 
 _ESVINT    esvscrft3_er(void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVI, _ESVI,
                     _ESVI, _ESVINT,  float, double *, _ESVI);
 _ESVINT    esvdcrft3_er(void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVI, _ESVI,
                     _ESVI, _ESVINT, double, double *, _ESVI);
#else
 void   esvscft2   (_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT,
                    _ESVINT, _ESVINT,  float, double *, _ESVINT, double *, _ESVINT);
 void   esvdcft2   (_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT,
                    _ESVINT, _ESVINT, double, double *, _ESVINT, double *, _ESVINT);
 void   esvscft2p   (_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,
                     _ESVINT, float, double *, _ESVINT, double *, _ESVINT);
 
 void   esvsrcft2   (_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT,
                     _ESVINT,  float, double *, _ESVINT, double *, _ESVINT, double *, _ESVINT);
 void   esvdrcft2   (_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT,
                     _ESVINT, double, double *, _ESVINT, double *, _ESVINT);
 
 void   esvscrft2   (_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT,
                     _ESVINT,  float, double *, _ESVINT, double *, _ESVINT, double *, _ESVINT);
 void   esvdcrft2   (_ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT,
                     _ESVINT, double, double *, _ESVINT, double *, _ESVINT);
 
 void   esvscft3   (void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,
                    _ESVINT, _ESVINT,  float, double *, _ESVINT);
 void   esvdcft3   (void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,
                    _ESVINT, _ESVINT, double, double *, _ESVINT);
 void   esvscft3p   (void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,
                     _ESVINT, _ESVINT, float, double *, _ESVINT);
 
 void   esvsrcft3   (void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,
                     _ESVINT, _ESVINT, float,  double *, _ESVINT);
 void   esvdrcft3   (void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,
                     _ESVINT, _ESVINT, double, double *, _ESVINT);
 
 void   esvscrft3   (void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,
                     _ESVINT, _ESVINT,  float, double *, _ESVINT);
 void   esvdcrft3   (void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT,
                     _ESVINT, _ESVINT, double, double *, _ESVINT);
#endif
 
/*  Convolutions/Correlations Subroutines  */
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvscon_er(_ESVINT, const float *, _ESVINT, const void *, _ESVINT, _ESVINT,
                   void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, double *,
                   _ESVI, double *, _ESVI);
 _ESVINT    esvscor_er(_ESVINT, const float *, _ESVINT, const void *, _ESVINT, _ESVINT,
                   void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, double *,
                   _ESVI, double *, _ESVI);
#else
 void   esvscon   (_ESVINT, const float *, _ESVINT, const void *, _ESVINT, _ESVINT,
                   void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, double *,
                   _ESVINT, double *, _ESVINT);
 void   esvscor   (_ESVINT, const float *, _ESVINT, const void *, _ESVINT, _ESVINT,
                   void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, double *,
                   _ESVINT, double *, _ESVINT);
#endif
 
void   esvscond(const float *, _ESVINT, const void *, _ESVINT, void *,
                _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT);
void   esvscord(const float *, _ESVINT, const void *, _ESVINT, void *,
                _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT);
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvsconf_er(_ESVINT, const float *, _ESVINT, void *, _ESVINT, _ESVINT,
                    void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT,
                    double *, _ESVI, double *, _ESVI);
 _ESVINT    esvscorf_er(_ESVINT, const float *, _ESVINT, void *, _ESVINT, _ESVINT,
                    void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT,
                    double *, _ESVI, double *, _ESVI);
#else
 void   esvsconf   (_ESVINT, const float *, _ESVINT, void *, _ESVINT, _ESVINT,
                    void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT,
                    double *, _ESVINT, double *, _ESVINT);
 void   esvscorf   (_ESVINT, const float *, _ESVINT, void *, _ESVINT, _ESVINT,
                    void *, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT, _ESVINT,
                    double *, _ESVINT, double *, _ESVINT);
#endif
 
void   esvsdcon(const  float *, _ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVINT,
                _ESVINT, _ESVINT, _ESVINT, _ESVINT);
void   esvddcon(const double *, _ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVINT,
                _ESVINT, _ESVINT, _ESVINT, _ESVINT);
void   esvsdcor(const  float *, _ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVINT,
                _ESVINT, _ESVINT, _ESVINT, _ESVINT);
void   esvddcor(const double *, _ESVINT, void *, _ESVINT, void *, _ESVINT, _ESVINT,
                _ESVINT, _ESVINT, _ESVINT, _ESVINT);
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvsacor_er(_ESVINT, const void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT,
                    _ESVINT, _ESVINT, _ESVINT, double *, _ESVI, double *, _ESVI);
 
 _ESVINT    esvsacorf_er(_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT,
                     _ESVINT, _ESVINT, double *, _ESVI, double *, _ESVI);
#else
 void   esvsacor   (_ESVINT, const void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT,
                    _ESVINT, _ESVINT, _ESVINT, double *, _ESVINT, double *, _ESVINT);
 
 void   esvsacorf   (_ESVINT, void *, _ESVINT, _ESVINT, void *, _ESVINT, _ESVINT, _ESVINT,
                     _ESVINT, _ESVINT, double *, _ESVINT, double *, _ESVINT);
#endif
 
/*  Related Computations Subroutines  */
 
void   esvspoly( const float *, _ESVINT, _ESVINT,const  float *, _ESVINT,  float *, _ESVINT, _ESVINT);
void   esvdpoly(const double *, _ESVINT, _ESVINT,const double *, _ESVINT, double *, _ESVINT, _ESVINT);
 
void   esvsizc( const float *, _ESVINT, _ESVINT, _ESVINT, _ESVINT *);
void   esvdizc(const double *, _ESVINT, _ESVINT, _ESVINT, _ESVINT *);
 
void   esvstrec( float, const  float *, _ESVINT, const  float *, _ESVINT,  float *, _ESVINT, _ESVINT, _ESVINT);
void   esvdtrec(double, const double *, _ESVINT, const double *, _ESVINT, double *, _ESVINT, _ESVINT, _ESVINT);
 
_ESVINT    esvsqint( float,  float,  float, const  float *, _ESVINT,  _ESVINT, const  float *, _ESVINT,  float *, _ESVINT, _ESVINT);
_ESVINT    esvdqint(double, double, double, const double *, _ESVINT,  _ESVINT, const double *, _ESVINT, double *, _ESVINT, _ESVINT);
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvswlev_er(const   float *, _ESVINT, const   float *, _ESVINT,   float *, _ESVINT, _ESVINT,  double *, _ESVI);
 _ESVINT    esvdwlev_er(const  double *, _ESVINT, const  double *, _ESVINT,  double *, _ESVINT, _ESVINT,  double *, _ESVI);
 _ESVINT    esvcwlev_er(const  _ESVCM *, _ESVINT, const  _ESVCM *, _ESVINT,  _ESVCM *, _ESVINT, _ESVINT, _ESVCOM *, _ESVI);
 _ESVINT    esvzwlev_er(const _ESVCOM *, _ESVINT, const _ESVCOM *, _ESVINT, _ESVCOM *, _ESVINT, _ESVINT, _ESVCOM *, _ESVI);
#else
 void   esvswlev   (const   float *, _ESVINT, const   float *, _ESVINT,   float *, _ESVINT, _ESVINT,  double *, _ESVINT);
 void   esvdwlev   (const  double *, _ESVINT, const  double *, _ESVINT,  double *, _ESVINT, _ESVINT,  double *, _ESVINT);
 void   esvcwlev   (const  _ESVCM *, _ESVINT, const  _ESVCM *, _ESVINT,  _ESVCM *, _ESVINT, _ESVINT, _ESVCOM *, _ESVINT);
 void   esvzwlev   (const _ESVCOM *, _ESVINT, const _ESVCOM *, _ESVINT, _ESVCOM *, _ESVINT, _ESVINT, _ESVCOM *, _ESVINT);
#endif
 
/*  Sorting and Searching Subroutines  */
 
void   esvisort(   _ESVINT *, _ESVINT, _ESVINT);
void   esvssort( float *, _ESVINT, _ESVINT);
void   esvdsort(double *, _ESVINT, _ESVINT);
 
void   esvisortx(   _ESVINT *, _ESVINT, _ESVINT, _ESVINT *);
void   esvssortx( float *, _ESVINT, _ESVINT, _ESVINT *);
void   esvdsortx(double *, _ESVINT, _ESVINT, _ESVINT *);
 
void   esvisorts(   _ESVINT *, _ESVINT, _ESVINT, _ESVINT *,    _ESVINT *, _ESVINT);
void   esvssorts( float *, _ESVINT, _ESVINT, _ESVINT *,  float *, _ESVINT);
void   esvdsorts(double *, _ESVINT, _ESVINT, _ESVINT *, double *, _ESVINT);
 
void   esvibsrch(const    _ESVINT *, _ESVINT, _ESVINT, const    _ESVINT *,
                 _ESVINT, _ESVINT, _ESVINT *, _ESVINT *, _ESVINT);
void   esvsbsrch(const  float *, _ESVINT, _ESVINT, const  float *,
                 _ESVINT, _ESVINT, _ESVINT *, _ESVINT *, _ESVINT);
void   esvdbsrch(const double *, _ESVINT, _ESVINT, const double *,
                 _ESVINT, _ESVINT, _ESVINT *, _ESVINT *, _ESVINT);
 
void   esvissrch(const    _ESVINT *, _ESVINT, _ESVINT, const    _ESVINT *,
                 _ESVINT, _ESVINT, _ESVINT, _ESVINT *);
void   esvsssrch(const  float *, _ESVINT, _ESVINT, const  float *,
                 _ESVINT, _ESVINT, _ESVINT, _ESVINT *);
void   esvdssrch(const double *, _ESVINT, _ESVINT, const double *,
                 _ESVINT, _ESVINT, _ESVINT, _ESVINT *);
 
/*  Interpolation Subroutines  */
 
void   esvspint(const  float *, const  float *, _ESVINT,  float *, _ESVI,
                const  float *,  float *, _ESVINT);
void   esvdpint(const double *, const double *, _ESVINT, double *, _ESVI,
                const double *, double *, _ESVINT);
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvstpint_er(const  float *, const  float *, _ESVINT, _ESVINT,
                     const  float *,  float *, _ESVINT,  float *, _ESVI);
 _ESVINT    esvdtpint_er(const double *, const double *, _ESVINT, _ESVINT,
                     const double *, double *, _ESVINT, double *, _ESVI);
#else
 void   esvstpint   (const  float *, const  float *, _ESVINT, _ESVINT,
                     const  float *,  float *, _ESVINT,  float *, _ESVINT);
 void   esvdtpint   (const double *, const double *, _ESVINT, _ESVINT,
                     const double *, double *, _ESVINT, double *, _ESVINT);
#endif
 
void   esvscsint(const  float *, const  float *, void *, _ESVINT, _ESVI,
                 const  float *,  float *, _ESVINT);
void   esvdcsint(const double *, const double *, void *, _ESVINT, _ESVI,
                 const double *, double *, _ESVINT);
 
#if defined(__ESVERR) || defined(_ESVERR)
_ESVINT    esvscsin2_er(const  float *, const  float *, const void *,
                    _ESVINT, _ESVINT, _ESVINT, const  float *, const  float *,
                    _ESVINT, _ESVINT, void *, _ESVINT,  float *, _ESVI);
_ESVINT    esvdcsin2_er(const double *, const double *, const void *,
                    _ESVINT, _ESVINT, _ESVINT, const double *, const double *,
                    _ESVINT, _ESVINT, void *, _ESVINT, double *, _ESVI);
#else
void   esvscsin2   (const  float *, const  float *, const void *,
                    _ESVINT, _ESVINT, _ESVINT, const  float *, const  float *,
                    _ESVINT, _ESVINT, void *, _ESVINT,  float *, _ESVINT);
void   esvdcsin2   (const double *, const double *, const void *,
                    _ESVINT, _ESVINT, _ESVINT, const double *, const double *,
                    _ESVINT, _ESVINT, void *, _ESVINT, double *, _ESVINT);
#endif
 
/*  Numerical Quadrature Subroutines  */
 
void   esvsptnq(const  float *, const  float *, _ESVINT, _ESVS, _ESVS);
void   esvdptnq(const double *, const double *, _ESVINT, _ESVD, _ESVD);
 
float  esvsglnq(void (*)(const  float *,  float *, const _ESVINT *),
                 float,  float, _ESVINT);
double esvdglnq(void (*)(const double *, double *, const _ESVINT *),
                double, double, _ESVINT);
 
float  esvsglnq2(void (*)(const  float *, const _ESVINT *, const  float *,
                 const _ESVINT *,  float *, const _ESVINT *),  float,  float,
                 _ESVINT,  float,  float, _ESVINT, const void *, _ESVINT);
double esvdglnq2(void (*)(const double *, const _ESVINT *, const double *,
                 const _ESVINT *, double *, const _ESVINT *), double, double,
                 _ESVINT, double, double, _ESVINT, const void *, _ESVINT);
 
float  esvsglgq(void (*)(const  float *,  float *, const _ESVINT *),
                 float,  float, _ESVINT);
double esvdglgq(void (*)(const double *, double *, const _ESVINT *),
                double, double, _ESVINT);
 
float  esvsgraq(void (*)(const  float *,  float *, const _ESVINT *),
                 float,  float, _ESVINT);
double esvdgraq(void (*)(const double *, double *, const _ESVINT *),
                double, double, _ESVINT);
 
float  esvsghmq(void (*)(const  float *,  float *, const _ESVINT *),
                 float,  float, _ESVINT);
double esvdghmq(void (*)(const double *, double *, const _ESVINT *),
                double, double, _ESVINT);
 
/*  Random Number Generation Subroutines  */
 
void   esvsurand(_ESVD, _ESVINT,  float *);
void   esvdurand(_ESVD, _ESVINT, double *);
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvsnrand_er(_ESVD, _ESVINT,  float *,  float *, _ESVI);
 _ESVINT    esvdnrand_er(_ESVD, _ESVINT, double *, double *, _ESVI);
#else
 void   esvsnrand(_ESVD, _ESVINT,  float *,  float *, _ESVINT);
 void   esvdnrand(_ESVD, _ESVINT, double *, double *, _ESVINT);
#endif
 
void esvsurxor(_ESVI, _ESVINT,  float *,  float *);
void esvdurxor(_ESVI, _ESVINT, double *, double *);
 
/*  Utility Subroutines  */
 
void   esveinfo(_ESVINT, _ESVI, _ESVI);
void   esvivsset(_ESVINT);
void   esvievops(_ESVINT);
_ESVINT    iessl(void);
void   esvstride(_ESVINT, _ESVINT, _ESVI, const char *, _ESVINT);
void   esvdsrsm(_ESVINT, const double *, const _ESVINT *, const _ESVINT *,
                _ESVINT, _ESVI, void *, void *, _ESVINT);
 
#if defined(__ESVERR) || defined(_ESVERR)
 _ESVINT    esvdgktrn_er(_ESVINT, double *, _ESVINT, _ESVINT *, double *, _ESVINT, _ESVINT *,
                     const _ESVINT *, double *, _ESVI);
 _ESVINT    esvdsktrn_er(_ESVINT, double *, _ESVINT, _ESVINT *, const _ESVINT *, double *, _ESVI);
#else
 void   esvdgktrn   (_ESVINT, double *, _ESVINT, _ESVINT *, double *, _ESVINT, _ESVINT *,
                     const _ESVINT *, double *, _ESVINT);
 void   esvdsktrn   (_ESVINT, double *, _ESVINT, _ESVINT *, const _ESVINT *, double *, _ESVINT);
#endif
 
#ifdef __cplusplus
}
/* extern "FORTRAN" { */
/* /\* These prototypes are for C++ only *\/ */
/*   void   errset(const _ESVINT &, const _ESVINT &, const _ESVINT &, const _ESVINT &, */
/*                 const void *, const _ESVINT &); */
/*   void   errsav(const _ESVINT &, void *); */
/*   void   errstr(const _ESVINT &, const void *); */
/* } */
#endif
