#include<config.h>

/*----------------------------------------------------------*/
/*!\file
  \brief Redefinitions of stdio functions for MPI. 

  $Id: qcdio_qprintf.C,v 1.6 2008-02-08 18:35:08 chulwoo Exp $
*/

/*  -----------------------------------------------------------
   CVS keywords
 
   $Author: chulwoo $ 
   $Date: 2008-02-08 18:35:08 $
   $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/qcdio/mpi/qcdio_qprintf.C,v 1.6 2008-02-08 18:35:08 chulwoo Exp $
   $Id: qcdio_qprintf.C,v 1.6 2008-02-08 18:35:08 chulwoo Exp $
   $Name: not supported by cvs2svn $
   $Locker:  $
   $RCSfile: qcdio_qprintf.C,v $
   $Revision: 1.6 $
   $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/qcdio/mpi/qcdio_qprintf.C,v $
   $State: Exp $  
 ----------------------------------------------------------------------*/

#include <util/qcdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <comms/sysfunc_cps.h>
#include <util/gjp.h>

namespace cps{

    namespace MPISCU{
	
	/*!
	  Works  like printf.
	  \param format The format string
	  \param ... Optional arguments to the format string 
	*/
	int printf( const char *format, ... ) {
	
	    if( UniqueID()!=1) return 0;

	    va_list ap;
	    int nbytes;
	    va_start(ap, format);
	    nbytes = ::vprintf( format, ap );
	    va_end(ap);

	    return nbytes;
	}

	/*!
	  Works  like fprintf.
	  \param stream The file handle to which to print
	  \param format The format string
	  \param ... Optional arguments to the format string 
	*/
	int fprintf( FILE *stream, const char *format, ... ) {

	    if( UniqueID() != 1 ) return 0;

    	    va_list ap;
	    int nbytes;
	    va_start(ap, format);
	    nbytes = ::vfprintf( stream, format, ap );
	    va_end(ap);
	    
	    return nbytes;
	}
  
	/*!
	  Works just  like printf.
	  \param format The format string
	  \param ... Optional arguments to the format string 
	*/
	int printf_all( const char *format, ... ) {
	    va_list ap;
	    int nbytes;
	
	    va_start(ap, format);
	    nbytes = ::vprintf( format, ap );
	    va_end(ap);

	    return nbytes;
	}
 
	/*!
	  Works  like printf, but prefixes the message with the node ID number
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
	    nbytes = ::vprintf( newformat, ap );
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
	    nbytes = ::vfprintf( stream, format, ap );
	    va_end(ap);

	    return nbytes;
	}

	/*!
	  Works like vprintf.
	  \param format The format string
	  \param arg The variable argument list.
	*/
	int vprintf(const char* format, va_list arg){
	    if( UniqueID() != 1 ) return 0;
	    return ::vprintf(format, arg);
	}

	/*!
	  Works like vprintf.
	  \param format The format string
	  \param arg The variable argument list.
	*/
	int vfprintf(FILE *f, const char* format, va_list arg){
	    if( UniqueID() != 1 ) return 0;
	    return ::vfprintf(f, format, arg);
	}
	
    } //namespace MPISCU
 
}//CPS_END_NAMESPACE 


