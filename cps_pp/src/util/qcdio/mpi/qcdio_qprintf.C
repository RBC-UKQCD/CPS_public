#include<config.h>

/*----------------------------------------------------------*/
/*!\file
  \brief Redefinitions of stdio functions for MPI. 

  $Id: qcdio_qprintf.C,v 1.3 2004-06-04 21:14:14 chulwoo Exp $
*/

/*  -----------------------------------------------------------
   CVS keywords
 
   $Author: chulwoo $ 
   $Date: 2004-06-04 21:14:14 $
   $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/qcdio/mpi/qcdio_qprintf.C,v 1.3 2004-06-04 21:14:14 chulwoo Exp $
   $Id: qcdio_qprintf.C,v 1.3 2004-06-04 21:14:14 chulwoo Exp $
   $Name: not supported by cvs2svn $
   $Locker:  $
   $RCSfile: qcdio_qprintf.C,v $
   $Revision: 1.3 $
   $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/qcdio/mpi/qcdio_qprintf.C,v $
   $State: Exp $  
 ----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <comms/sysfunc.h>
#include <util/gjp.h>

namespace cps{

namespace MPISCU{
  /*!
    Works just like printf.
    \param format The format string
    \param ... Optional arguments to the format string 
   */
  int printf( const char *format, ... ) {
        va_list ap;
        int nbytes;
	
	if( UniqueID() == 1 ) {
  	    va_start(ap, format);
  	    nbytes = vprintf( format, ap );
	    va_end(ap);
	} else {
	    nbytes = 0;
	}

        return nbytes;
  }

  /*!
    Works just like fprintf.
    \param stream The file handle to which to print
    \param format The format string
    \param ... Optional arguments to the format string 
   */
  int fprintf( FILE *stream, const char *format, ... ) {
        va_list ap;
        int nbytes;

	if( UniqueID() == 1 ) {
	    va_start(ap, format);
	    nbytes = vfprintf( stream, format, ap );
	    va_end(ap);
	} else {
	    nbytes = 0;
	}

        return nbytes;
  }
  
  /*!
    Works just like printf.
    \param format The format string
    \param ... Optional arguments to the format string 
   */
  int printf_all( const char *format, ... ) {
        va_list ap;
        int nbytes;
	
	va_start(ap, format);
	nbytes = vprintf( format, ap );
	va_end(ap);

        return nbytes;
  }
 
  /*!
    Works just like printf, but prefixes the message with the node ID number
    and the cartesian coordinates of the node in the node grid.
    \param format The format string
    \param ... Optional arguments to the format string 
   */
  int printf_allid( const char *format, ... ) {
        va_list ap;
        int nbytes;
	char* newformat = new char[strlen(format)+100];
	
	sprintf(newformat,"%i [%i,%i,%i,%i] %s",
			UniqueID(),
			GJP.XnodeCoor(),
			GJP.YnodeCoor(),
			GJP.ZnodeCoor(),
			GJP.TnodeCoor(),
			format );
	
	va_start(ap, newformat);
	nbytes = vprintf( newformat, ap );
	va_end(ap);

        return nbytes;
  }


  /*!
    Works just like fprintf.
    \param stream The file handle to which to print    
    \param format The format string
    \param ... Optional arguments to the format string 
   */

  int fprintf_all( FILE *stream, const char *format, ... ) {
        va_list ap;
        int nbytes;

	va_start(ap, format);
	nbytes = vfprintf( stream, format, ap );
	va_end(ap);

        return nbytes;
  }

} //namespace MPISCU
 
}//CPS_END_NAMESPACE 


