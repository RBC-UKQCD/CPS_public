#include<config.h>

CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*!\file
  \brief  Declarations for the serial emulation of the QCDSP/QCDOC communications  layer.
  
  $Id: sysfunc_noarch.h,v 1.9 2008/02/08 18:35:05 chulwoo Exp $
*/
/*----------------------------------------------------------------------
  The Sysfunc Comms Interface: sysfunc_cps.h

  Declarations for the noarch (fake) comms-layer

  Chulwoo Jung
  -----------------------------------------------------------
  CVS keywords
 
  $Author: chulwoo $
  $Date: 2008/02/08 18:35:05 $
  $Header: /space/cvs/cps/cps++/include/comms/sysfunc_noarch.h,v 1.9 2008/02/08 18:35:05 chulwoo Exp $
  $Id: sysfunc_noarch.h,v 1.9 2008/02/08 18:35:05 chulwoo Exp $
  $Name: v5_0_16_hantao_io_test_v7 $
  $Locker:  $
  $RCSfile: sysfunc_noarch.h,v $
  $Revision: 1.9 $
  $Source: /space/cvs/cps/cps++/include/comms/sysfunc_noarch.h,v $
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
inline unsigned int Barrier(){return 1;}

inline void broadcast( void *, size_t size){}
#endif


CPS_END_NAMESPACE




