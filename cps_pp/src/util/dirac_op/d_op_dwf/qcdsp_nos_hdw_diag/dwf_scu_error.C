#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdsp_nos_hdw_diag/dwf_scu_error.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include <time.h>
#include<util/wilson.h>
#include<util/error.h>
CPS_START_NAMESPACE


extern "C" void dwf_scu_error(void)
{
  // Set the value of the clock
  wfm_scu_diag[1] = clock();

  // Exit
  ERR.Hardware(" ", "dwf_scu_error", 
  "Failed to complete transfer.\n\t19 words of diagnostic info in address %x\n\tSee wilson.h for definitions.\n", 
  wfm_scu_diag);
}
CPS_END_NAMESPACE
