#include<config.h>
CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*! The MPI comms direction and flag enums:

  This is very closely based on the original scu_enum.h from QCDSP.
  The only difference is the order of the direction enum `SCUDir'.
  This order reflects the way in which the MPI implementation
  identifies the `opposite direction', i.e. if were looking in
  direction dir:
    opposite_dir = (dir+NDIM)%(2*NDIM)
  where NDIM is the number of dimensions.  If this causes 
  problems, it can always be changed at a later stage.

  A.N.Jackson: ajackson@epcc.ed.ac.uk                       
  -----------------------------------------------------------
  CVS keywords
 
  $Author: mcneile $
  $Date: 2003-06-22 13:34:52 $
  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/scu_enum.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
  $Id: scu_enum.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
  $Name: not supported by cvs2svn $
  $Locker:  $
  $RCSfile: scu_enum.h,v $
  $Revision: 1.1.1.1 $
  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/scu_enum.h,v $
  $State: Exp $  */
/*----------------------------------------------------------*/

CPS_END_NAMESPACE
#include<config.h>
CPS_START_NAMESPACE

/* Allow the MPI stuff to be switched out, thus avoiding compiler
   errors (for the time being). */
#ifdef INCLUDE_MPI_SCU

#ifndef INCLUDED_SCU_ENUM_H
#define INCLUDED_SCU_ENUM_H


/*--------------------------------------------------------------------
//  The following enum's define the symbolic names for the physics
//  directions and axes.  The order here MUST not be changed, unless the
//  order of the English translations is changed in
//  $os/qio_qker/static_const/scu_sys.C.
//
//  The SCU system calls make use of relations like dir/2 = axis
//  for efficiency.
//--------------------------------------------------------------------
*/

//--------------------------------------------------------------------
//  Map physics directions to wires
//--------------------------------------------------------------------
/*! Like QCDSP, MPI implementation uses +t,-t,+x,-x,+y,-y,+z,-z.  Must go + then -, but the order of t, x, y and z is arbitrary. */
enum SCUDir { 
  SCU_TP, 
  SCU_TM, 
  SCU_XP, 
  SCU_XM, 
  SCU_YP, 
  SCU_YM, 
  SCU_ZP, 
  SCU_ZM,
  SCU_NoDir = -1 
};


//--------------------------------------------------------------------
//!  Label axes
//--------------------------------------------------------------------
enum SCUAxis { SCU_T, SCU_X, SCU_Y, SCU_Z, SCU_NoAxis = -1 };


//--------------------------------------------------------------------
//!  Define send and receive for SCU calls
//--------------------------------------------------------------------
enum SCUXR { SCU_REC, SCU_SEND = 8, SCU_NoXR = -1 };

#endif

#endif /* INCLUDE_MPI_SCU */


CPS_END_NAMESPACE
