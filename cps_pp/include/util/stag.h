#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:37 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/stag.h,v 1.5 2004-08-18 11:57:37 zs Exp $
//  $Id: stag.h,v 1.5 2004-08-18 11:57:37 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: stag.h,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/stag.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/****************************************************************************/
/*                                                                          */
/* stag.h                                                                   */
/*                                                                          */
/* C header file for the staggered library.                                 */
/*                                                                          */
/****************************************************************************/

#ifndef INCLUDED_STAG_LIB_H
#define INCLUDED_STAG_LIB_H

//------------------------------------------------------------------
// The staggered dslash operator
// dir_flag is flag which takes value 0 when all direction contribute to D,
// 1 - when only the special anisotropic direction contributes to D,
// 2 - when all  except the special anisotropic direction. 
//------------------------------------------------------------------
extern "C" void stag_dirac(IFloat *f_out, IFloat *f_in, int cb, int dag , 
		      int dir_flag = 0);

//------------------------------------------------------------------
// Initialize all global variables and address tables needed by
// staggered dirac() function.
//------------------------------------------------------------------
extern "C" void stag_dirac_init(const void *gauge_field_addr);
extern "C" void stag_dirac_init_g();

//------------------------------------------------------------------
// Destroy all address tables needed by staggered dirac() function
//------------------------------------------------------------------
extern "C" void stag_destroy_dirac_buf();
extern "C" void stag_destroy_dirac_buf_g();

#endif

CPS_END_NAMESPACE
