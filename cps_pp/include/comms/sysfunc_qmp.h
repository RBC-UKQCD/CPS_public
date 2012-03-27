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

#ifndef SYSFUNC_QMP_H
#define SYSFUNC_QMP_H
#ifndef UNIFORM_SEED_NO_COMMS
#include<qmp.h>
#endif
#include <comms/scu_enum.h>
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

//! A barrier function.
//inline void sync(){QMP_barrier();}
#ifndef HAVE_SYNC
unsigned int sync();
#endif

//! Gets the direction used internally by the comms layer.
//int SCURemap( SCUDir dir );


namespace QMPSCU {
  void init_qmp();
  void init_qmp(int * argc, char *** argv);
  void destroy_qmp();
}

CPS_END_NAMESPACE
#endif




