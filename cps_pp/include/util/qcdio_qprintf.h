#include<config.h>

/*----------------------------------------------------------*/
/*!\file
  \brief Basic IO function declarations,

  Separated from the rest of qcdio.h to avoid troublesome
  header file interdependancies.

  $Id: qcdio_qprintf.h,v 1.4 2004-06-02 11:03:19 zs Exp $  
*/
/*  A.N.Jackson: ajackson@epcc.ed.ac.uk                       */
/*----------------------------------------------------------*/

#ifndef INCLUDED_QCDIO_QPRINTF_H_
#define INCLUDED_QCDIO_QPRINTF_H_


#include <stdio.h>


#if TARGET == cpsMPI

CPS_BEGIN_NAMESPACE

namespace MPISCU{
  //! Reimplementation of printf that prints from only a single node
  int printf( const char *format, ... );

  //! Reimplementation of fprintf that prints from only a single node
  int fprintf( FILE *stream, const char *format, ... );
  
  //! Reimplementation of printf that prints on all nodes
  int printf_all( const char *format, ... );

  //! Reimplementation of printf that prints on all nodes with node ID prefix
  int printf_allid( const char *format, ... );

  //! Reimplementation of fprintf that prints on all nodes
  int fprintf_all( FILE *stream, const char *format, ... );
}


CPS_END_NAMESPACE

#endif
#endif

