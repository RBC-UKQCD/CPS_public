#include <config.h>
CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*!\file
  \brief  The MPI communication direction and flag enumerations:

  $Id: scu_enum.h,v 1.8 2003-10-23 13:38:59 zs Exp $
*/
/*---------------------------------------------------------------
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
 
  $Author: zs $
  $Date: 2003-10-23 13:38:59 $
  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/scu_enum.h,v 1.8 2003-10-23 13:38:59 zs Exp $
  $Id: scu_enum.h,v 1.8 2003-10-23 13:38:59 zs Exp $
  $Name: not supported by cvs2svn $
  $Locker:  $
  $RCSfile: scu_enum.h,v $
  $Revision: 1.8 $
  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/scu_enum.h,v $
  $State: Exp $  */
/*----------------------------------------------------------*/

CPS_END_NAMESPACE
CPS_START_NAMESPACE

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
//! Definition of the physical directions
/*!
  The order must alternate positive and negative directions,
  but the order of t, x, y and z is arbitrary.
*/
//--------------------------------------------------------------------
enum SCUDir { 
  SCU_TP,                      /*!< +t */
  SCU_TM,                      /*!< -t */
  SCU_XP,                      /*!< +x */  
  SCU_XM,                      /*!< -x */ 
  SCU_YP,                      /*!< +y */  
  SCU_YM,                      /*!< -y */  
  SCU_ZP,                      /*!< +z */  
  SCU_ZM,                      /*!< -z */ 
  SCU_NoDir=-1          /*!< Null direction */  
};




//--------------------------------------------------------------------
//! Definition of the physical axes
//--------------------------------------------------------------------
enum SCUAxis {
    SCU_T,          /*!< t axis */
    SCU_X,          /*!< x axis */
    SCU_Y,          /*!< y axis */
    SCU_Z,          /*!< z axis */
    SCU_NoAxis = -1 /*!< Dummy axis for serial code */
};


//--------------------------------------------------------------------
//!  Flags denoting 'send' and 'receive' in communications routines
//--------------------------------------------------------------------
enum SCUXR {
    SCU_REC,         /*!< Receive */
    SCU_SEND = 8,    /*!< Send */
    SCU_NoXR = -1    /*!< Dummy flag for serial code */
};




#endif

CPS_END_NAMESPACE



