#include<config.h>

CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*!\file
  \brief  Declarations for the serial emulation of the QCDSP/QCDOC communications  layer.
  
  $Id: sysfunc_noarch.h,v 1.9 2008-02-08 18:35:05 chulwoo Exp $
*/
/*----------------------------------------------------------------------
  The Sysfunc Comms Interface: sysfunc_cps.h

  Declarations for the noarch (fake) comms-layer

  Chulwoo Jung
  -----------------------------------------------------------
  CVS keywords
 
  $Author: chulwoo $
  $Date: 2008-02-08 18:35:05 $
  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/sysfunc_noarch.h,v 1.9 2008-02-08 18:35:05 chulwoo Exp $
  $Id: sysfunc_noarch.h,v 1.9 2008-02-08 18:35:05 chulwoo Exp $
  $Name: not supported by cvs2svn $
  $Locker:  $
  $RCSfile: sysfunc_noarch.h,v $
  $Revision: 1.9 $
  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/sysfunc_noarch.h,v $
  $State: Exp $  */
/*----------------------------------------------------------*/

CPS_END_NAMESPACE
CPS_START_NAMESPACE

#ifndef NOARCH_SYSFUNC_H
#define NOARCH_SYSFUNC_H





//--------------------------------------------------------------------
/* Public interface subroutines for the comms: */
//--------------------------------------------------------------------
/*! \defgroup comms  Communications objects and functions 
  @{ */
 
//! Gets an ID which is unique for each node.
inline int UniqueID(){return 0;}

/*  Functions to return physics four dimensional coordinates of the
  node.  The physics coordinate axes are labeled T, X, Y, Z*/

inline int CoorT(){return 0;}  //!< Gets the grid coordinate of this node in the T direction.
inline int CoorX(){return 0;}  //!< Gets the grid coordinate of this node in the X direction.
inline int CoorY(){return 0;}  //!< Gets the grid coordinate of this node in the Y direction.
inline int CoorZ(){return 0;}  //!< Gets the grid coordinate of this node in the Z direction.
inline int CoorS(){return 0;}  //!< Gets the grid coordinate of this node in the Z direction.
inline int CoorW(){return 0;}  //!< Gets the grid coordinate of this node in the Z direction.

inline int SizeT(){return 1;} //!< Gets the size of the grid  in the T direction.
inline int SizeX(){return 1;} //!< Gets the size of the grid  in the X direction.
inline int SizeY(){return 1;} //!< Gets the size of the grid  in the Y direction.
inline int SizeZ(){return 1;} //!< Gets the size of the grid  in the Z direction.
inline int SizeS(){return 1;} //!< Gets the size of the grid  in the Z direction.
inline int SizeW(){return 1;} //!< Gets the size of the grid  in the Z direction.

//! Returns the total number of nodes in the processor grid.
inline int NumNodes(){return 1;}

// Random number seeds that on QCDSP were loaded at boot time by the QOS.

static const unsigned int SERIAL_SEED = 112319;

inline unsigned int Seed(){return SERIAL_SEED;}   //!< Gets a RNG seed.
inline unsigned int SeedS(){return SERIAL_SEED;}  //!< Gets a RNG seed.
inline unsigned int SeedT(){return SERIAL_SEED;}  //!< Gets a RNG seed.
inline unsigned int SeedST(){return SERIAL_SEED;} //!< Gets a RNG seed.

#ifndef HAVE_SYNC
//! A barrier function.
inline unsigned int sync(){return 1;}
#endif

// done up to here yet.
#if 0 
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

#endif


#endif


CPS_END_NAMESPACE




