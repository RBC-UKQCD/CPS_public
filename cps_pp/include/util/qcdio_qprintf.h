#include<config.h>
CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*!\file
  \brief Basic IO function declarations,

  Separated from the rest of qcdio.h to avoid troublesome
  header file interdependancies.

  $Id: qcdio_qprintf.h,v 1.2 2003-07-24 16:53:53 zs Exp $  
*/
/*  A.N.Jackson: ajackson@epcc.ed.ac.uk                       */
/*----------------------------------------------------------*/

#ifndef INCLUDED_QCDIO_QPRINTF_H_
#define INCLUDED_QCDIO_QPRINTF_H_

CPS_END_NAMESPACE
#include <stdio.h>
#include <util/data_types.h>
CPS_START_NAMESPACE


#ifdef PARALLEL
#define printf qprintf           //!<Parallel reimplementation of printf
#define fprintf qfprintf         //!<Parallel reimplementation of fprintf

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
