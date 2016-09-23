#ifndef INCLUDED_QCDIO_QPRINTF_H_
#define INCLUDED_QCDIO_QPRINTF_H_

#include<config.h>

/*----------------------------------------------------------*/
/*!\file
  \brief Basic IO function declarations,

  Separated from the rest of qcdio.h to avoid troublesome
  header file interdependancies.

  $Id: qcdio_qprintf.h,v 1.5.464.1 2012/07/09 16:29:19 yinnht Exp $  
*/
/*  A.N.Jackson: ajackson@epcc.ed.ac.uk                       */
/*----------------------------------------------------------*/

#include <stdio.h>
#include <stdarg.h>


#if TARGET == cpsMPI

//namespace cps{
CPS_START_NAMESPACE

namespace MPISCU{
  //! Reimplementation of printf that prints from only a single node.
  int printf( const char *format, ... );

  //! Reimplementation of fprintf that prints from only a single node.
  int fprintf( FILE *stream, const char *format, ... );
  
  //! Reimplementation of printf that prints on all nodes.
  int printf_all( const char *format, ... );

  //! Reimplementation of printf that prints on all nodes with node ID prefix.
  int printf_allid( const char *format, ... );

  //! Reimplementation of fprintf that prints on all nodes.
  int fprintf_all( FILE *stream, const char *format, ... );

  //! Reimplementation of vprintf that prints from only a single node.
  int vprintf(const char*,  va_list);

  //! Reimplementation of vfprintf that prints from only a single node.
  int vfprintf(FILE*, const char*,  va_list);
  
  

} // namespace MPISCU


//} // namespace cps
CPS_END_NAMESPACE

#endif
#endif

