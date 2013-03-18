#include<config.h>
/*----------------------------------------------------------*/
/*!\file
  \brief Prototypes of gauge configuration IO functions.

  $Id: qcdio.h,v 1.12 2013-03-18 19:33:13 chulwoo Exp $
*/
/*2  A.N.Jackson: ajackson@epcc.ed.ac.uk                      
  -----------------------------------------------------------
   CVS keywords
 
   $Author: chulwoo $
   $Date: 2013-03-18 19:33:13 $
   $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/qcdio.h,v 1.12 2013-03-18 19:33:13 chulwoo Exp $
   $Id: qcdio.h,v 1.12 2013-03-18 19:33:13 chulwoo Exp $
   $Name: not supported by cvs2svn $
   $Locker:  $
   $RCSfile: qcdio.h,v $
   $Revision: 1.12 $
   $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/qcdio.h,v $
   $State: Exp $  */
/*----------------------------------------------------------*/

#ifndef INCLUDED_QCDIO_H_
#define INCLUDED_QCDIO_H_

#include <stdio.h>
#include <stdarg.h>
#include <util/data_types.h>
#include <util/lattice.h>
#include <util/qcdio_qprintf.h>
#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>

CPS_START_NAMESPACE

#ifndef GAUGE_CONF_PREC
/*! The default precision at which gauge configurations are stored in files. */
#define GAUGE_CONF_PREC (sizeof(Float))
#endif

#ifndef SWAP_BYTE_ORDER
/*! When loading gauge configurations, should the byte order be swapped before
 * putting the results into memory?  To swap the bytes by default, set this
 to 1, or set it to 0 to leave the byte order alone by default. */
#define SWAP_BYTE_ORDER 0
#endif

#ifndef TRANSPOSE_THE_MATRICES
/*! When the SU(3) matrices have been loaded, the row vs. column convention
  must be the same as for the code which generated the configuration.  To
  transpose the matrix by default, set this to be 1, otherwise use 0. */
#define TRANSPOSE_THE_MATRICES 1
#endif

  /* ------------------------------------------------------
  	Gauge configuration I/O:
    ------------------------------------------------------ */
  
  //! Routine to define whether the loaded data should be normalized.
  void qcdio_set_normalize( int );

  //! Routine for loading a gauge configuration
  void qload_gauge( char* fprefix, Lattice& lat,
             int prec = GAUGE_CONF_PREC,
	     int swap = SWAP_BYTE_ORDER,
	     int transp = TRANSPOSE_THE_MATRICES );

  //! Routine for saving the current gauge configuration
  void qsave_gauge( char* fprefix, Lattice& lat,
             int prec = GAUGE_CONF_PREC,
             int swap = SWAP_BYTE_ORDER,
	     int transp = TRANSPOSE_THE_MATRICES );

// File IO

  /*!\defgroup std_file_io Functions for doing file IO on a parallel machine
    @{ */

//! Type of IO
enum FileIoType{
    ZERO_ONLY,          /*!< Write from a single node only. */
    ADD_ID,              /*!< Each node writes to a seperate file which
			  has the node ID number appended to its name. */
    ALL_NODES
};

//! Opens a file
/*!
  This works like \a fopen in the C standard library, except that on a
  parallel machine we can decide whether just one node does the IO, or
  whether each node performs IO to a seperate file which has the node
  number appended to its name.
  
  \param type What sort of IO to do.
  \param filename the main name of the file
  \param mode the IO mode (as for \a fopen)
  \return The handle of the closed file, or NULL on failure.
  \pre In the mode where all nodes open a file, the length of the final
  filename, \e i.e. \a filename with a '.' and the node number appended,
  cannot exceed 200 characters.
  \post The file is opened by all nodes or just one. In the latter case
  a dummy  file handle is returned.
*/
FILE *Fopen( FileIoType type, const char *filename, const char *mode);

//! Closes a file
/*!
  This works like \a fclose in the C standard library
  
  \param type Ignored.
  \param stream The file handle to close.
  \return Normally, \c EOF or zero on error, but if the file was opened by
  a single node and this function is called by a different node (in which
  case \a stream will be the dummy file handle), then 1 is returned.
*/
int Fclose( FileIoType type, FILE *stream);

//! Prints to a file
/*!
  This works like \a fprint in the C standard library, except that on a
  parallel machine we can decide whether just one node does the IO, or
  whether each node performs IO to a seperate file which has the node
  number appended to its name.
  
  \param type Ignored
  \param stream The file handle
  \param format Format string - as in the C standard library \a fprintf
  \param ... Parameters, as in the C standard library \a fprintf
  \return Normally, the number of bytes written, but if the file was opened in
  ZERO_ONLY mode and this function is called by a different mode (in which
  case \a stream will be the dummy file handle), then 1 is returned.
*/

size_t Fwrite( const void *ptr, size_t size, size_t n, FILE *stream);
size_t Fread(void *ptr, size_t size, size_t n, FILE *stream);
int Fflush(FILE *stream);

//! Read & Write from a file
/*!
  This works like \a fread & fwrite in the C standard library, except that on a
  parallel machine we can decide whether just one node does the IO, or
  whether each node performs IO to a seperate file which has the node
  number appended to its name.
  
  \param stream The file handle
  \return Normally, the number of elements written, but if the file was opened 
  in ZERO_ONLY mode and this function is called by a different mode (in which
  case \a stream will be the dummy file handle), then n is returned from 
  Fwrite().
*/


int Fprintf( FileIoType type, FILE *stream, const char *format,...);

//! Prints a variable-length argument list to a file
/*!
  This works like \a vfprint in the C standard library, except that on a
  parallel machine we can decide whether just one node does the IO, or
  whether each node performs IO to a seperate file which has the node
  number appended to its name.
  
  \param type Ignored
  \param stream The file handle
  \param format Format string - as in the C standard library \a fprintf
  \param ap Parameter liss, as in the C standard library \a fprintf
  \return Normally, the number of bytes written, but if the file was opened in
  ZERO_ONLY mode and this function is called by a different mode  (in which
  case \a stream will be the dummy file handle), then 1 is returned.
*/
int Vfprintf( FileIoType type, FILE *stream, const char *format, va_list ap);

//! Opens a file
/*!
  This works like \a fopen in the C standard library, except that on a
  parallel machine just one node does the IO
  
  \param filename the main name of the file
  \param mode the IO mode (as for \a fopen)
  \post Only one node opens the file. A dummy handle is returned  on other
  nodes.
*/
inline FILE *Fopen( const char *filename, const char *mode)
    { return Fopen( ZERO_ONLY, filename, mode);}

//! Closes a file
/*!
  This works like \a fclose in the C standard library
  
  \param stream The file handle to close.
  \return Normally, \c EOF or zero on error, but if the file was opened by
  a single node and this function is called by a different node (in which
  case \a stream will be the dummy file handle), then 1 is returned.
*/
inline int Fclose( FILE *stream){return Fclose(ZERO_ONLY,stream);}

//! Prints to a file
/*!
  This works like \a fprint in the C standard library, except that on a
  parallel machine we can decide whether just one node does the IO, or
  whether each node performs IO to a seperate file which has the node
  number appended to its name.
  
  \param stream The file handle
  \param format Format string - as in the C standard library \a fprintf
  \param ... Parameters, as in the C standard library \a fprintf
  \return Normally, the number of bytes written, but if the file was opened in
  ZERO_ONLY mode and this function is called by a different mode (in which
  case \a stream will be the dummy file handle), then 1 is returned.
*/

int Fprintf( FILE *stream, const char *format,...);

//! Prints a variable-length argument list to a file
/*!
  This works like \a vfprint in the C standard library, except that on a
  parallel machine we can decide whether just one node does the IO, or
  whether each node performs IO to a seperate file which has the node
  number appended to its name.
  
  \param stream The file handle
  \param format Format string - as in the C standard library \a fprintf
  \param ap Parameter liss, as in the C standard library \a fprintf
  \return Normally, the number of bytes written, but if the file was opened in
  ZERO_ONLY mode and this function is called by a different mode, then 1
  is returned.
*/
inline int Vfprintf( FILE *stream, const char *format, va_list ap)
    { return Vfprintf(ZERO_ONLY,stream, format, ap);}

CPS_END_NAMESPACE
#endif

/*! @} */
