#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief Utility routines for the asqtad fermions Dirac operator

  $Id: asqtad.h,v 1.5 2004-09-02 16:59:57 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-09-02 16:59:57 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/asqtad.h,v 1.5 2004-09-02 16:59:57 zs Exp $
//  $Id: asqtad.h,v 1.5 2004-09-02 16:59:57 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: asqtad.h,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/asqtad.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/****************************************************************************/
/*                                                                          */
/* asqtad.h                                                                   */
/*                                                                          */
/* C header file for the improved staggered (Asqtad) library.     */
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
extern "C" void asqtad_dirac(IFloat *f_out, IFloat *f_in, int cb, int dag , 
		      int dir_flag = 0);

//------------------------------------------------------------------
// Initialize all global variables and address tables needed by
// staggered dirac() function.
//------------------------------------------------------------------
extern "C" void asqtad_dirac_init(const void *gauge_field_addr);
extern "C" void asqtad_dirac_init_g ();

//------------------------------------------------------------------
// Destroy all address tables needed by staggered dirac() function
//------------------------------------------------------------------
extern "C" void asqtad_destroy_dirac_buf();
extern "C" void asqtad_destroy_dirac_buf_g();

#endif

CPS_END_NAMESPACE
