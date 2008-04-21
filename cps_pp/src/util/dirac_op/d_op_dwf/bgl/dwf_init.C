#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/bgl/dwf_init.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
// 11/26/97
//
// dwf_int:
//
// This routine performs all initializations needed before dwf
// func are called. It sets the addressing related arrays and 
// reserves memory for the needed temporary buffers. It only needs
// to be called once at the begining of the program (or after a 
// dwf_end call) before any number of calls to dwf funcs are made.
//
// WARNING:
//
// This set of routines will work only if the node sublattices have
// even number of sites in each direction.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/dwf.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
CPS_START_NAMESPACE



void dwf_init(Dwf *dwf_p)
{
  char *cname = " ";
  char *fname = "dwf_init(Dwf*)";
  VRB.Func(cname,fname);

//------------------------------------------------------------------
// Do initializations before the wilson library can be used
// Initialization involve memory allocation.
//------------------------------------------------------------------
  static Wilson wilson_struct;
  dwf_p->wilson_p = &wilson_struct;
  wilson_init(dwf_p->wilson_p);

//------------------------------------------------------------------
// Allocate memory for two temporary fermion checkerboard fields  
//------------------------------------------------------------------
  int f_size = 24 * GJP.VolNodeSites() * GJP.SnodeSites() / 2; 

  dwf_p->frm_tmp1 = (IFloat *) smalloc(f_size*sizeof(IFloat));
  if(dwf_p->frm_tmp1 == 0)
    ERR.Pointer(cname,fname, "frm_tmp1");
  VRB.Smalloc(cname,fname,
	      "frm_tmp1", dwf_p->frm_tmp1, f_size*sizeof(IFloat));

  dwf_p->frm_tmp2 = (IFloat *) smalloc(f_size*sizeof(IFloat));
  if(dwf_p->frm_tmp2 == 0)
    ERR.Pointer(cname,fname, "frm_tmp2");
  VRB.Smalloc(cname,fname,
	      "frm_tmp2", dwf_p->frm_tmp2, f_size*sizeof(IFloat));

//------------------------------------------------------------------
// Allocate memory for a 12 word communications buffer needed
// for the spread-out case.
//------------------------------------------------------------------
  dwf_p->comm_buf = (IFloat *) smalloc(cname,fname, "comm_buf", f_size * sizeof(IFloat));


//------------------------------------------------------------------
// Set the dwf coefficients
//------------------------------------------------------------------
  dwf_p->vol_4d = GJP.VolNodeSites();
  dwf_p->ls = GJP.SnodeSites();
  dwf_p->dwf_kappa = 
    1.0 / ( 2 * ( 4 + GJP.DwfA5Inv() - GJP.DwfHeight() ) );

}


CPS_END_NAMESPACE
