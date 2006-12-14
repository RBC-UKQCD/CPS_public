#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/bgl/dwf_end.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
// 10/16/97
//
// dwf_end:
//
// This routine frees any memory that was allocated by dwf_init
//
// WARNING:
//
// This set of routines will work only if the node sublattices have
// even number of sites in each direction.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/dwf.h>
#include<util/smalloc.h>
#include<util/verbose.h>
CPS_START_NAMESPACE


void dwf_end( Dwf *dwf_p)
{
  char *cname = " ";
  char *fname = "dwf_end(Dwf*)";
  VRB.Func(cname,fname);

  //------------------------------------------------------------------
  // Free memory of the 12 word communications buffer needed
  // for the spread-out case.
  //------------------------------------------------------------------
  VRB.Sfree(cname,fname, "comm_buf", dwf_p->comm_buf);
  sfree(dwf_p->comm_buf);

  //----------------------------------------------------------------
  // Free temporary femrion fields 2 and 1
  //----------------------------------------------------------------
  VRB.Sfree(cname,fname, "frm_tmp2", dwf_p->frm_tmp2);
  sfree(dwf_p->frm_tmp2);

  VRB.Sfree(cname,fname, "frm_tmp1", dwf_p->frm_tmp1);
  sfree(dwf_p->frm_tmp1);

  //----------------------------------------------------------------
  // Un-initialize the wilson library. Memory is set free here.
  //----------------------------------------------------------------
  wilson_end(dwf_p->wilson_p);


}
CPS_END_NAMESPACE
