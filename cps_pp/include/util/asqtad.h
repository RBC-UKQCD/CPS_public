#ifndef ASQTAD_H
#define ASQTAD_H
#include<config.h>
/*!\file
  \brief Utility routines for the asqtad fermions Dirac operator

  $Id: asqtad.h,v 1.8 2007-06-06 16:06:22 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2007-06-06 16:06:22 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/asqtad.h,v 1.8 2007-06-06 16:06:22 chulwoo Exp $
//  $Id: asqtad.h,v 1.8 2007-06-06 16:06:22 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: asqtad.h,v $
//  $Revision: 1.8 $
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


#if TARGET == QCDOC
#include <util/asqtad_int.h>
#endif

CPS_START_NAMESPACE
#ifdef __cplusplus
extern "C"{
#endif
//------------------------------------------------------------------
// The staggered dslash operator
// dir_flag is flag which takes value 0 when all direction contribute to D,
// 1 - when only the special anisotropic direction contributes to D,
// 2 - when all  except the special anisotropic direction. 
//------------------------------------------------------------------
#if TARGET == QCDOC
 inline void asqtad_dirac(IFloat *f_out, IFloat *f_in, int cb, int dag , 
		      int dir_flag = 0){
  asqd_p->dirac(f_out,f_in,cb,dag);
}
#else
 void asqtad_dirac(IFloat *f_out, IFloat *f_in, int cb, int dag , 
		      int dir_flag = 0);
#endif

//------------------------------------------------------------------
// Initialize all global variables and address tables needed by
// staggered dirac() function.
//------------------------------------------------------------------
void asqtad_dirac_init(Fasqtad *lat);
void asqtad_dirac_init_g (IFloat *frm_p);

//------------------------------------------------------------------
// Destroy all address tables needed by staggered dirac() function
//------------------------------------------------------------------
void asqtad_destroy_dirac_buf();
void asqtad_destroy_dirac_buf_g();
#ifdef __cplusplus
}
#endif


CPS_END_NAMESPACE
#endif
