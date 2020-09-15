#include<config.h>
#ifdef USE_SSE
#include "../sse/dwf_init_sse.C"
#else
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /space/cvs/cps/cps++/src/util/dirac_op/d_op_dwf/noarch/dwf_init.C,v $
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
  size_t f_size = 24 * GJP.VolNodeSites() * GJP.SnodeSites() / 2; 
//  if(!UniqueID()) printf("Allocating temp ferms. Base ferm size is %d\n",f_size);//DEBUG
  if(GJP.Gparity()){
    f_size*=2;
    if(!UniqueID()) printf("G-parity active, ferm size is doubled to %d\n",f_size);//DEBUG
  }
  if(!UniqueID()) fflush(stdout);//DEBUG

  dwf_p->frm_tmp1 = (IFloat *) smalloc(cname,fname,"frm_tmp1",f_size*sizeof(IFloat));
//  if(dwf_p->frm_tmp1 == 0)
//    ERR.Pointer(cname,fname, "frm_tmp1");
//  VRB.Smalloc(cname,fname,
//	      "frm_tmp1", dwf_p->frm_tmp1, f_size*sizeof(IFloat));

  dwf_p->frm_tmp2 = (IFloat *) smalloc(cname,fname,"frm_tmp2",f_size*sizeof(IFloat));
//  if(dwf_p->frm_tmp2 == 0)
//    ERR.Pointer(cname,fname, "frm_tmp2");
//  VRB.Smalloc(cname,fname,
//	      "frm_tmp2", dwf_p->frm_tmp2, f_size*sizeof(IFloat));
  dwf_p->frm_tmp3 = (IFloat *) smalloc(cname,fname,"frm_tmp3",f_size*sizeof(IFloat));
  VRB.Debug(cname,fname,"frm_tmp1 frm_tmp2 frm_tmp3= %p %p %p\n",
	dwf_p->frm_tmp1, dwf_p->frm_tmp2, dwf_p->frm_tmp3);

//------------------------------------------------------------------
// Allocate memory for a 12 word communications buffer needed
// for the spread-out case.
//------------------------------------------------------------------
  dwf_p->comm_buf = (IFloat *) smalloc(12 * sizeof(IFloat));
  if(dwf_p->comm_buf == 0)
    ERR.Pointer(cname,fname, "comm_buf");
  VRB.Smalloc(cname,fname,
	      "comm_buf", dwf_p->comm_buf, 12*sizeof(IFloat));


//------------------------------------------------------------------
// Set the dwf coefficients
//------------------------------------------------------------------
  dwf_p->vol_4d = GJP.VolNodeSites();
  dwf_p->ls = GJP.SnodeSites();
  dwf_p->dwf_kappa = 
    1.0 / ( 2 * ( 4 + GJP.DwfA5Inv() - GJP.DwfHeight() ) );

}


CPS_END_NAMESPACE
#endif
