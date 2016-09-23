#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  gamma_5 multiplication routine

  Used by derivatives of the FwilsonTypes class.

  $Id: gamma_5.C,v 1.1 2004/12/06 22:05:23 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004/12/06 22:05:23 $
//  $Header: /space/cvs/cps/cps++/src/util/lattice/f_wilson_types/qcdoc/gamma_5.C,v 1.1 2004/12/06 22:05:23 chulwoo Exp $
//  $Id: gamma_5.C,v 1.1 2004/12/06 22:05:23 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $Revision: 1.1 $
//  $Source: /space/cvs/cps/cps++/src/util/lattice/f_wilson_types/qcdoc/gamma_5.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
// gamma_5(IFloat *v_out, IFloat *v_in, int num_sites):
//
// v_out = Gamma5 * v_in. Gamme5 is in the chiral basis
//
//          [ 1  0  0  0]
// Gamma5 = [ 0  1  0  0]
//          [ 0  0 -1  0]
//          [ 0  0  0 -1]
//
// num_sites is the number of sites. It is assumed
// that each site has 24 components.
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/data_types.h>
CPS_START_NAMESPACE

void gamma_5(IFloat *v_out, IFloat *v_in, int num_sites) ;

void gamma_5(IFloat *v_out, IFloat *v_in, int num_sites) 
{
  IFloat *p_out = v_out ;
  IFloat *p_in = v_in ;

  int half_site_size = 12 ;

  for (int site=0; site<num_sites; ++site) {
    int comp;
    for (comp=0; comp<half_site_size; ++comp) {
      *p_out++ = *p_in++ ;
    }
    for (comp=0; comp<half_site_size; ++comp) {
      *p_out++ = -*p_in++ ;
    }
  }
}

CPS_END_NAMESPACE
