#include<config.h>

// If not MPI then must be either QCDOC or QCDSP
#if TARGET != cpsMPI
#include <sysfunc.h>
#else

/*----------------------------------------------------------*/
/*!\file
  \brief  Declarations for the MPI implementation of the QCDSP/QCDOC communications  layer.
  
  $Id: sysfunc.h,v 1.9 2004-08-18 11:57:36 zs Exp $
*/

/*----------------------------------------------------------------------
  The Sysfunc Comms Interface: sysfunc.h

  Declarations for the MPI implementation of the QCDSP SCU
  comms-layer.

  A.N.Jackson: ajackson@epcc.ed.ac.uk               
  -----------------------------------------------------------
  CVS keywords
 
  $Author: zs $
  $Date: 2004-08-18 11:57:36 $
  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/sysfunc.h,v 1.9 2004-08-18 11:57:36 zs Exp $
  $Id: sysfunc.h,v 1.9 2004-08-18 11:57:36 zs Exp $
  $Name: not supported by cvs2svn $
  $Locker:  $
  $RCSfile: sysfunc.h,v $
  $Revision: 1.9 $
  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/sysfunc.h,v $
  $State: Exp $  */
/*----------------------------------------------------------*/

#ifndef INCLUDED_SYSFUNC_H
#define INCLUDED_SYSFUNC_H

#include <comms/scu_dir_arg.h>
#include <comms/mpi_requests.h>

CPS_START_NAMESPACE



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
int CoorS();  //!< Gets the grid coordinate of this node in the S direction.

int SizeT(); //!< Gets the size of the grid  in the T direction.
int SizeX(); //!< Gets the size of the grid  in the X direction.
int SizeY(); //!< Gets the size of the grid  in the Y direction.
int SizeZ(); //!< Gets the size of the grid  in the Z direction.
int SizeS(); //!< Gets the size of the grid  in the S direction.

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


namespace MPISCU {

//! The global comms-layer initialization flag.
    extern bool Is_Initialised;


/* Useful extra subroutines, part of the MPI version:          */


//! Communicates the processor grid dimensions to the MPI-SCU layer.
    void set_pe_grid(int x, int y, int z, int t, int s=1);

//! Initialises the MPI communications layer.
    void CommsInit();

//! Performs a clean exit from the MPI communications layer.
    void CommsFinalize();

//! Computes a global sum directly using MPI 
    void GlobalSum(Type_tag t, size_t tsize, int n, void *ivec, void *ovec );

//! Reports an error.
    void RaiseError( char* errstr );

//! Reports an error.
    void RaiseError( const char* errstring );

/*! @} */

/*-------------------------------------------------------------------------*/
/*              Implementation-specific internal subroutines:              */
/*              If this were a class, these would be private.              */
/*-------------------------------------------------------------------------*/

//! An implementation-specific internal subroutine
    void Trans( void* addr, MPI_Datatype mpi_dt, SCUDir dir, SCUXR sendrx );

//! An implementation-specific internal subroutine
    void ParseCommsParam( char *);

//! An implementation-specific internal subroutine
    char *CommsStringTokenizer(char* str, const char* tokens, char** tok_pos );

//! An implementation-specific internal subroutine
    MPI_Datatype MPITypeConv( Type_tag t, size_t tsize );

//! An implementation-specific internal subroutine
    unsigned int ReadSeedFile(  );


} //namespace MPISCU

CPS_END_NAMESPACE

#endif //INCLUDED_SYSFUNC_H

#endif // TARGET




