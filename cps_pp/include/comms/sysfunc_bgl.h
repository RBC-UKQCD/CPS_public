#include<config.h>

#ifndef INCLUDED_SYSFUNC_BGL_H
#define INCLUDED_SYSFUNC_BGL_H
/*----------------------------------------------------------*/
/*!\file
  \brief  Declarations for the MPI implementation of the QCDSP/QCDOC communications  layer.
*/
/*----------------------------------------------------------*/

/*----------------------------------------------------------*/
/*  
  The Sysfunc Comms Interface: sysfunc_cps.h

  Declarations for the MPI implementation of the QCDSP SCU
  comms-layer.
*/
/*----------------------------------------------------------*/

/* Allow the MPI stuff to be switched out, thus avoiding compiler
   errors (for the time being). */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <comms/scu_enum.h>
//#include <comms/scu_dir_arg.h>

CPS_START_NAMESPACE

//! Number of dimensions.
#ifndef NDIM                     /* If NDIM has not been specified: */
#define NDIM 4                   
#endif


//--------------------------------------------------------------------
/* Public interface subroutines for the comms: */
//--------------------------------------------------------------------
/*! \defgroup comms  Communications objects and functions 
  @{ */
 
/* TorusX, TorusY, TorusZ, TorusT return 1 if the corresponding direction
   is a torus and return 0 if it is a mesh */
int TorusT(); //! 1 if the dir is a torus, 0 if it is not
int TorusX(); //! 1 if the dir is a torus, 0 if it is not
int TorusY(); //! 1 if the dir is a torus, 0 if it is not
int TorusZ(); //! 1 if the dir is a torus, 0 if it is not

//! Gets an ID which is unique for each node.
int UniqueID();

/*  Functions to return physics four dimensional coordinates of the
  node.  The physics coordinate axes are labeled T, X, Y, Z*/

int CoorT();  //!< Gets the grid coordinate of this node in the T direction.
int CoorX();  //!< Gets the grid coordinate of this node in the X direction.
int CoorY();  //!< Gets the grid coordinate of this node in the Y direction.
int CoorZ();  //!< Gets the grid coordinate of this node in the Z direction.
inline int CoorS() {return 0;}
inline int CoorW() {return 0;}

int SizeT(); //!< Gets the size of the grid  in the T direction.
int SizeX(); //!< Gets the size of the grid  in the X direction.
int SizeY(); //!< Gets the size of the grid  in the Y direction.
int SizeZ(); //!< Gets the size of the grid  in the Z direction.
inline int SizeS() {return 1;}
inline int SizeW() {return 1;}

//! Returns the total number of nodes in the processor grid.
int NumNodes();

/* Random number seeds that on QCDSP were loaded at boot time by the
  QOS.  N.B. Note that the MPI implementation of these functions are
  not currently compliant with the function comments */
unsigned int Seed();   //!< Gets a RNG seed.
unsigned int SeedS();  //!< Gets a RNG seed.
unsigned int SeedT();  //!< Gets a RNG seed.
unsigned int SeedST(); //!< Gets a RNG seed.

#ifndef HAVE_SYNC
//! A barrier function.
unsigned int sync();
#endif

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
#endif



//}// End of extern "C".
CPS_END_NAMESPACE
#endif


