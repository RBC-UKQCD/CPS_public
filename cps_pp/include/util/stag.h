#include<config.h>
#include <util/lattice.h>
CPS_START_NAMESPACE
/*!\file
  \brief Declaration of utility routines for the staggered fermions Dirac operator.

  $Id: stag.h,v 1.8 2010-07-26 17:39:53 chulwoo Exp $
*/
//--------------------------------------------------------------------
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2010-07-26 17:39:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/stag.h,v 1.8 2010-07-26 17:39:53 chulwoo Exp $
//  $Id: stag.h,v 1.8 2010-07-26 17:39:53 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: stag.h,v $
//  $Revision: 1.8 $
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
extern "C"{
void stag_dirac(Vector *f_out, Vector *f_in, int cb, int dag , 
		      int dir_flag = 0);

//------------------------------------------------------------------
// Initialize all global variables and address tables needed by
// staggered dirac() function.
//------------------------------------------------------------------
void stag_dirac_init(const void *gauge_field_addr);
void stag_dirac_init_g();
void stag_dirac_init_with_lat (Lattice& lat);

//------------------------------------------------------------------
// Destroy all address tables needed by staggered dirac() function
//------------------------------------------------------------------
void stag_destroy_dirac_buf();
void stag_destroy_dirac_buf_g();
}

#endif

CPS_END_NAMESPACE
