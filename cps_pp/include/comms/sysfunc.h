#include<config.h>
CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*! The Sysfunc Comms Interface: sysfunc.h

  Declarations for the MPI implementation of the QCDSP SCU
  comms-layer.

  A.N.Jackson: ajackson@epcc.ed.ac.uk               
  -----------------------------------------------------------
  CVS keywords
 
  $Author: mcneile $
  $Date: 2003-06-22 13:34:52 $
  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/sysfunc.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
  $Id: sysfunc.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
  $Name: not supported by cvs2svn $
  $Locker:  $
  $RCSfile: sysfunc.h,v $
  $Revision: 1.1.1.1 $
  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/sysfunc.h,v $
  $State: Exp $  */
/*----------------------------------------------------------*/

CPS_END_NAMESPACE
#include<config.h>
CPS_START_NAMESPACE

/* Allow the MPI stuff to be switched out, thus avoiding compiler
   errors (for the time being). */
#ifdef INCLUDE_MPI_SCU

#ifndef INCLUDED_SYSFUNC_H
#define INCLUDED_SYSFUNC_H

CPS_END_NAMESPACE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include<config.h>
#include<comms/comms/scu_dir_arg.h>
#include<comms/mpi_requests.h>
CPS_START_NAMESPACE

#ifndef NDIM                     /* If NDIM has not been specified: */
#define NDIM 4                   /*!< Number of dimensions */ 
#endif

/* Macros for the required environment variable names: */
/*! Name of the environment variable holding (or pointing to) 
  the simulation option definitions */
#ifndef COMMS_ENVVAR
#define COMMS_ENVVAR "COMMS_DEF" 
#endif

/*! Default name for the comms definition file. */
#ifndef COMMS_DEFFILE
#define COMMS_DEFFILE "commsMPI.def"
#endif

/*! Max size of temporary strings */
#define STRING_MAX_LEN  10000 

/*! Note that this interface cannot be extern C'd because it uses
  overloaded subroutines */
//extern "C" {

//! Declaration for the global comms-layer initialization flag.
extern int commsMPI_init;

//--------------------------------------------------------------------
/*! Public interface subroutines for the comms: */
//--------------------------------------------------------------------
//! An ID which is unique for each node.
int UniqueID();

/*!  Functions to return physics four dimensional coordinates of the
  node.  The physics coordinate axes are labeled T, X, Y, Z*/
int CoorT();
int CoorX();
int CoorY();
int CoorZ();

int SizeT();
int SizeX();
int SizeY();
int SizeZ();

/*! Returns the total number of nodes in the PE grid */
int NumNodes();

/*! Random number seeds that on QCDSP were loaded at boot time by the
  QOS.  N.B. Note that the MPI implementation of these functions are
  not currently compliant with the function comments */
unsigned int Seed();   /*<! Seed is different for each node and is
			 changed every time the machine is reset.*/
unsigned int SeedS();  /*<! SeedS is the same for each node
			 (spatially fixed, hence the S), but changes
			 in time*/
unsigned int SeedT();  /*<! SeedT is different for each node, but is
			 fixed in time (the T), so it is unchanged by
			 a reset.*/
unsigned int SeedST(); /*<! SeedST is the same for each node
			 (spatially fixed, hence the S), and the same
			 after every reset (fixed time, hence T).*/

/*! sync() is a function which blocks further code execution until all
  nodes in the machine have begun executing the code in the sync()
  routine.  Originally returned the idle-time, but now returns 0. */
unsigned int sync();

/*! SCURemap: On QCDSP, this function returned the explicit wire
  number (0 - 7) of the physics direction given by dir. In the MPI
  version, this returns the internal direction from the cartesian
  communicator which corresponds to the given physics direction. */
int SCURemap( SCUDir dir );

//! The following are the primary functions for generic SCU transfers:

void SCUTrans( SCUDirArg * arg );
/*!< This performs a single transfer, specified by arg. */

void SCUTrans( SCUDirArg ** arg, int n );
/*!< This function performs n transfers, as specified by the array of
  ACUDirArg objects in arg */

void SCUTrans( SCUDirArg * arg, unsigned int * offset, int n );
/*!<  This function does multiple transfers (n of them) for a specified
  direction.  All transfers have the same block, stride and number
  of blocks, but different addresses.  The address field of the
  arg object is the base address.  Each transfer is started at
  a specified offest relative to the base. */

/*! SCUSetDMA: Used to set up the block, stride and number of blocks,
  but no transfers are done. The transfer is instead inisialised using
  SCUTransAddr. */
void SCUSetDMA( SCUDirArg * arg );
void SCUSetDMA( SCUDirArg ** arg, int n );

/*! SCUTransAddr: This function also perform SCU transfers, but the
  existing block, stride and number of blocks, as defined by a call to
  SCUSetDMA.  Therefore, this command _must_ be preceded by a call to
  SCUSetDMA. The base-addresses of the data are taken from the
  SCUTransAddr arguments. */
void SCUTransAddr( SCUDirArg * arg );
void SCUTransAddr( SCUDirArg ** arg, int n );

/*!  SCUComplete() only returns when all transfers on this PE are completed. */
void SCUTransComplete(void);

/*-------------------------------------------------------------------------*/
/*             Useful extra subroutines, part of the MPI version:          */
/*-------------------------------------------------------------------------*/
//! Controls the comms initialisation, parses comm-params and calls MPI_Init:
void SCUCommsInit( void );

//! Controls the comms finalization; clean exit via MPI_Finalize:
void SCUCommsFinalize( void );

//! Perform a global sum directly using MPI */
void SCUGlobalSum(Type_tag t, size_t tsize, int n, void *ivec, void *ovec );

//! Wrapper for the comms-error reporting mechanism. */
void SCURaiseError( char* errstr );
//! Extra error wrapper to deal with string literals. */
void SCURaiseError( const char* errstring );


/*-------------------------------------------------------------------------*/
/*              Implementation-specific internal subroutines:              */
/*              If this were a class, these would be private.              */
/*-------------------------------------------------------------------------*/
//! Basic MPI transfer call, on which all others are based:
void SCUTrans_mpi( void* addr, MPI_Datatype mpi_dt, SCUDir dir, SCUXR sendrx );

//! Parses the comms parameters:
void MPIParseCommsParam(void);

//! String tokenizer, coded here to ensure portability:
char *MPICommsStringTokenizer(char* str, const char* tokens, char** tok_pos );

//! This performs on-the-fly type+size to MPI_Datatype conversion.
MPI_Datatype SCUMPITypeConv( Type_tag t, size_t tsize );

//! Reads a seed for every PE from a file specified during intialisation:
unsigned int SCUReadSeedFile( void );


//}// End of extern "C".
#endif

#endif /* INCLUDE_MPI_SCU */
CPS_END_NAMESPACE
