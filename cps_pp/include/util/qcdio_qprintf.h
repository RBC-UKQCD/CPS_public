#include<config.h>
CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*! The QCD I/O Interface:

  $Id: qcdio_qprintf.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $

  A.N.Jackson: ajackson@epcc.ed.ac.uk                       */
/*----------------------------------------------------------*/

#ifndef INCLUDED_QCDIO_QPRINTF_H_
#define INCLUDED_QCDIO_QPRINTF_H_

CPS_END_NAMESPACE
#include <stdio.h>
#include<config.h>
#include<util/data_types.h>
CPS_START_NAMESPACE


#ifdef PARALLEL
#define printf qprintf
#define fprintf qfprintf

  /* ------------------------------------------------------
  	Basic standard I/O:
    Separated from the rest of qcdio.h to avoid troublesome
    header file interdependancies.
    ------------------------------------------------------ */

  //! Reimplementation of printf that prints on only the zeroth node
  int qprintf( const char *format, ... );

  //! Reimplementation of fprintf that prints on only the zeroth node
  int qfprintf( FILE *stream, const char *format, ... );
  
  //! Reimplementation of printf that prints on all nodes
  int qprintf_all( const char *format, ... );

  //! Reimplementation of printf that prints on all nodes with node ID prefix
  int qprintf_allid( const char *format, ... );

  //! Reimplementation of fprintf that prints on all nodes
  int qfprintf_all( FILE *stream, const char *format, ... );

#endif

#endif
CPS_END_NAMESPACE
