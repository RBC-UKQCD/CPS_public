#include <config.h>

/*----------------------------------------------------------*/
/*!\file
  \brief  The MPI communication direction and flag enumerations:

  $Id: scu_enum.h,v 1.10 2004-08-18 11:57:36 zs Exp $
*/
/*---------------------------------------------------------------
  CVS keywords
 
  $Author: zs $
  $Date: 2004-08-18 11:57:36 $
  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/scu_enum.h,v 1.10 2004-08-18 11:57:36 zs Exp $
  $Id: scu_enum.h,v 1.10 2004-08-18 11:57:36 zs Exp $
  $Name: not supported by cvs2svn $
  $Locker:  $
  $RCSfile: scu_enum.h,v $
  $Revision: 1.10 $
  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/scu_enum.h,v $
  $State: Exp $  */
/*----------------------------------------------------------*/


#ifndef INCLUDED_SCU_ENUM_H
#define INCLUDED_SCU_ENUM_H

CPS_START_NAMESPACE

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
  SCU_SP,                      /*!< +s */  
  SCU_SM,                      /*!< -s */ 
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
    SCU_S,          /*!< s axis */    
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


CPS_END_NAMESPACE

#endif





