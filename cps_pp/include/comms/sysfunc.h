#include<config.h>

// If not MPI then must be either QCDOC or QCDSP
#if TARGET != cpsMPI
#include <sysfunc.h>
#else

CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*!\file
  \brief  Declarations for the MPI implementation of the QCDSP/QCDOC communications  layer.
  
  $Id: sysfunc.h,v 1.6 2004-02-06 12:20:43 zs Exp $
*/
/*----------------------------------------------------------------------
  The Sysfunc Comms Interface: sysfunc.h

  Declarations for the MPI implementation of the QCDSP SCU
  comms-layer.

  A.N.Jackson: ajackson@epcc.ed.ac.uk               
  -----------------------------------------------------------
  CVS keywords
 
  $Author: zs $
  $Date: 2004-02-06 12:20:43 $
  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/sysfunc.h,v 1.6 2004-02-06 12:20:43 zs Exp $
  $Id: sysfunc.h,v 1.6 2004-02-06 12:20:43 zs Exp $
  $Name: not supported by cvs2svn $
  $Locker:  $
  $RCSfile: sysfunc.h,v $
  $Revision: 1.6 $
  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/sysfunc.h,v $
  $State: Exp $  */
/*----------------------------------------------------------*/

CPS_END_NAMESPACE
CPS_START_NAMESPACE

#ifndef INCLUDED_SYSFUNC_H
#define INCLUDED_SYSFUNC_H

CPS_END_NAMESPACE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <comms/scu_dir_arg.h>
#include <comms/mpi_requests.h>
CPS_START_NAMESPACE

//! Number of dimensions.
#ifndef NDIM                     /* If NDIM has not been specified: */
#define NDIM 4                   
#endif

/* Macros for the required environment variable names: */

//! Name of the environment variable with the parallel execution parameters.
/*! This environment variable might define the parameters directly or
  name a file which does.
*/
#ifndef COMMS_ENVVAR
#define COMMS_ENVVAR "COMMS_DEF" 
#endif

//! Default name for the file containing the parallel execution parameters.
#ifndef COMMS_DEFFILE
#define COMMS_DEFFILE "commsMPI.def"
#endif

/*! Max size of temporary strings */
#define STRING_MAX_LEN  10000 

/* Note that this interface cannot be extern C'd because it uses
  overloaded subroutines */
//extern "C" {

//! The global comms-layer initialization flag.
extern int commsMPI_init;

//--------------------------------------------------------------------
/* Public interface subroutines for the comms: */
//--------------------------------------------------------------------
/*! \defgroup comms  Communications objects and functions 
  @{ */
 
//! Gets an ID which is unique for each node.
int UniqueID();

/*  Functions to return physics four dimensional coordinates of the
  node.  The physics coordinate axes are labeled T, X, Y, Z*/

int CoorT();  //!< Gets the grid coordinate of this node in the T direction.
int CoorX();  //!< Gets the grid coordinate of this node in the X direction.
int CoorY();  //!< Gets the grid coordinate of this node in the Y direction.
int CoorZ();  //!< Gets the grid coordinate of this node in the Z direction.

int SizeT(); //!< Gets the size of the grid  in the T direction.
int SizeX(); //!< Gets the size of the grid  in the X direction.
int SizeY(); //!< Gets the size of the grid  in the Y direction.
int SizeZ(); //!< Gets the size of the grid  in the Z direction.

//! Returns the total number of nodes in the processor grid.
int NumNodes();

/* Random number seeds that on QCDSP were loaded at boot time by the
  QOS.  N.B. Note that the MPI implementation of these functions are
  not currently compliant with the function comments */
unsigned int Seed();   //!< Gets a RNG seed.
unsigned int SeedS();  //!< Gets a RNG seed.
unsigned int SeedT();  //!< Gets a RNG seed.
unsigned int SeedST(); //!< Gets a RNG seed.

//! A barrier function.
unsigned int sync();

//! Gets the direction used internally by the comms layer.
int SCURemap( SCUDir dir );

/* The following are the primary functions for generic SCU transfers: */

//! Generic single communication.
void SCUTrans( SCUDirArg * arg );

//! Generic multiple communication.
void SCUTrans( SCUDirArg ** arg, int n );

//! Does a number of similar communications.
void SCUTrans( SCUDirArg * arg, unsigned int * offset, int n );

//! Initialise a data transfer.
void SCUSetDMA( SCUDirArg * arg );

//! Initialise multiple data transfers.
void SCUSetDMA( SCUDirArg ** arg, int n );

//! Performs a previously set-up data transfer.
void SCUTransAddr( SCUDirArg * arg );

//! Performs multiple previously set-up data transfers.
void SCUTransAddr( SCUDirArg ** arg, int n );

//! A communications barrier function,
void SCUTransComplete(void);


/*-------------------------------------------------------------------------*/
/*             Useful extra subroutines, part of the MPI version:          */
/*-------------------------------------------------------------------------*/
  
//! Initialises the MPI communications layer.
void SCUCommsInit( void );

//! Performs a clean exit from the MPI communications layer.
void SCUCommsFinalize( void );

//! Computes a global sum directly using MPI 
void SCUGlobalSum(Type_tag t, size_t tsize, int n, void *ivec, void *ovec );

//! Reports an error.
void SCURaiseError( char* errstr );

//! Reports an error.
void SCURaiseError( const char* errstring );

/*! @} */

/*-------------------------------------------------------------------------*/
/*              Implementation-specific internal subroutines:              */
/*              If this were a class, these would be private.              */
/*-------------------------------------------------------------------------*/

//! An implementation-specific internal subroutine
void SCUTrans_mpi( void* addr, MPI_Datatype mpi_dt, SCUDir dir, SCUXR sendrx );

//! An implementation-specific internal subroutine
void MPIParseCommsParam(void);

//! An implementation-specific internal subroutine
char *MPICommsStringTokenizer(char* str, const char* tokens, char** tok_pos );

//! An implementation-specific internal subroutine
MPI_Datatype SCUMPITypeConv( Type_tag t, size_t tsize );

//! An implementation-specific internal subroutine
unsigned int SCUReadSeedFile( void );


//}// End of extern "C".
#endif



CPS_END_NAMESPACE
#endif
