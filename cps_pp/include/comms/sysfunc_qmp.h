#include<config.h>

/*----------------------------------------------------------*/
/*!\file
  \brief  Declarations for the QMP implementation of the QCDSP/QCDOC communications  layer.
  
*/
/*----------------------------------------------------------------------
  The Sysfunc Comms Interface: sysfunc_qmp.h

  Declarations for the QMP implementation of the QCDSP/QCDOC SCU
  comms-layer.

  M. Cheng michaelc@phys.columbia.edu               
  -----------------------------------------------------------*/

/*----------------------------------------------------------*/
#include <comms/scu_enum.h>
CPS_START_NAMESPACE

#ifndef SYSFUNC_QMP_H
#define SYSFUNC_QMP_H






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
int CoorZ();  //!< Gets the grid coordinate of this node in the Z direction
int CoorS();  //!< Gets the grid coordinate of this node in the S direction.
int CoorW();  //!< Gets the grid coordinate of this node in the S direction.

int SizeT(); //!< Gets the size of the grid  in the T direction.
int SizeX(); //!< Gets the size of the grid  in the X direction.
int SizeY(); //!< Gets the size of the grid  in the Y direction.
int SizeZ(); //!< Gets the size of the grid  in the Z direction.
int SizeS(); //!< Gets the size of the grid  in the S direction.
int SizeW(); //!< Gets the size of the grid  in the S direction.

//! Returns the total number of nodes in the processor grid.
int NumNodes();

unsigned int Seed();   //!< Gets a RNG seed.
unsigned int SeedS();  //!< Gets a RNG seed.
unsigned int SeedT();  //!< Gets a RNG seed.
unsigned int SeedST(); //!< Gets a RNG seed.

#ifndef HAVE_SYNC
//! A barrier function.
unsigned int sync();
#endif

//! Gets the direction used internally by the comms layer.
int SCURemap( SCUDir dir );

/* The following are the primary functions for generic SCU transfers: */

#if 0
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

namespace QMPSCU {
  void init_qmp();
  void init_qmp(int * argc, char *** argv);
  void destroy_qmp();
}
#endif

CPS_END_NAMESPACE




