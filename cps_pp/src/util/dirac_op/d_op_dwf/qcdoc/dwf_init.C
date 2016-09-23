#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2009/03/23 19:13:32 $
//  $Header: /space/cvs/cps/cps++/src/util/dirac_op/d_op_dwf/qcdoc/dwf_init.C,v 1.10 2009/03/23 19:13:32 chulwoo Exp $
//  $Id: dwf_init.C,v 1.10 2009/03/23 19:13:32 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: dwf_init.C,v $
//  $Revision: 1.10 $
//  $Source: /space/cvs/cps/cps++/src/util/dirac_op/d_op_dwf/qcdoc/dwf_init.C,v $
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
#include<util/wfm.h>
CPS_START_NAMESPACE


static WilsonArg wil;
#include <stdio.h>

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

  wil.local_latt[0] = GJP.XnodeSites();
  wil.local_latt[1] = GJP.YnodeSites();
  wil.local_latt[2] = GJP.ZnodeSites();
  wil.local_latt[3] = GJP.TnodeSites();

  if ( GJP.Xnodes() > 1 )  wil.local_comm[0] = 0;
  else  wil.local_comm[0] = 1;
  if ( GJP.Ynodes() > 1 )  wil.local_comm[1] = 0;
  else  wil.local_comm[1] = 1;
  if ( GJP.Znodes() > 1 )  wil.local_comm[2] = 0;
  else  wil.local_comm[2] = 1;
  if ( GJP.Tnodes() > 1 )  wil.local_comm[3] = 0;
  else  wil.local_comm[3] = 1;

  wil.local_comm[0] = 0;
  wil.local_comm[1] = 0;
  wil.local_comm[2] = 0;
  wil.local_comm[3] = 0;

  wfm_vec_init(&wil);
//~~ 
//~~ twisted mass fermions:  added second argument WilsonArg *wil
//~~ for transfer of *spinor_tmp from WilsonArg to Wilson
//~~ 
//~~  wilson_compat_init(dwf_p->wilson_p);
  wilson_compat_init(dwf_p->wilson_p, &wil);
//~~


//------------------------------------------------------------------
// Allocate memory for two temporary fermion checkerboard fields  
//------------------------------------------------------------------
  int f_size = 24 * GJP.VolNodeSites() * GJP.SnodeSites() / 2; 

//  dwf_p->frm_tmp1 = (IFloat *) smalloc(f_size*sizeof(IFloat));
  dwf_p->frm_tmp1 = (IFloat *) fmalloc(f_size*sizeof(IFloat));
  if(dwf_p->frm_tmp1 == 0)
    ERR.Pointer(cname,fname, "frm_tmp1");
  VRB.Smalloc(cname,fname,
	      "frm_tmp1", dwf_p->frm_tmp1, f_size*sizeof(IFloat));

//  for(int i=0;i<f_size;i++) (dwf_p->frm_tmp1)[i]=0.0;

//  dwf_p->frm_tmp2 = (IFloat *) smalloc(f_size*sizeof(IFloat));
  dwf_p->frm_tmp2 = (IFloat *) fmalloc(f_size*sizeof(IFloat));
  if(dwf_p->frm_tmp2 == 0)
    ERR.Pointer(cname,fname, "frm_tmp2");
  VRB.Smalloc(cname,fname,
	      "frm_tmp2", dwf_p->frm_tmp2, f_size*sizeof(IFloat));
//  for(int i=0;i<f_size;i++) (dwf_p->frm_tmp2)[i]=0.0;


//------------------------------------------------------------------
// Set the dwf coefficients
//------------------------------------------------------------------
  dwf_p->vol_4d = GJP.VolNodeSites();
  dwf_p->ls = GJP.SnodeSites();
  dwf_p->dwf_kappa = 
    1.0 / ( 2 * ( 4 + GJP.DwfA5Inv() - GJP.DwfHeight() ) );

//------------------------------------------------------------------
// Allocate memory for a communications buffer needed
// for the spread-out case.
//------------------------------------------------------------------
  int ls_stride = 12*dwf_p->vol_4d;
//  dwf_p->comm_buf = (IFloat *) fmalloc(ls_stride*sizeof(IFloat));
  int ls = GJP.Snodes();
  if (ls >1)
  dwf_p->comm_buf = (IFloat *) qalloc(QFAST|QNONCACHE,ls_stride*sizeof(IFloat));
  else
  dwf_p->comm_buf = (IFloat *) smalloc(12*sizeof(IFloat));
  if(dwf_p->comm_buf == 0)
    ERR.Pointer(cname,fname, "comm_buf");
  VRB.Smalloc(cname,fname,
	      "comm_buf", dwf_p->comm_buf, 12*sizeof(IFloat));

  if (ls >1){
  dwf_p->PlusArg[0] = new SCUDirArgIR;
    (dwf_p->PlusArg[0]) ->Init (dwf_p->comm_buf,SCU_SP,SCU_SEND, ls_stride*sizeof(IFloat),1,0,IR_14);
  dwf_p->PlusArg[1] = new SCUDirArgIR;
    (dwf_p->PlusArg[1])->Init(dwf_p->comm_buf,SCU_SM,SCU_REC, ls_stride*sizeof(IFloat),1,0,IR_14);
	dwf_p->Plus = new SCUDirArgMulti;
    (dwf_p->Plus)->Init(dwf_p->PlusArg,2);

  dwf_p->MinusArg[0] = new SCUDirArgIR;
    (dwf_p->MinusArg[0]) ->Init (dwf_p->comm_buf ,SCU_SM,SCU_SEND,ls_stride*sizeof(IFloat),1,0,IR_15);
  dwf_p->MinusArg[1] = new SCUDirArgIR;
    (dwf_p->MinusArg[1])->Init(dwf_p->comm_buf,SCU_SP,SCU_REC, ls_stride*sizeof(IFloat),1,0,IR_15);
	dwf_p->Minus = new SCUDirArgMulti;
    (dwf_p->Minus)->Init(dwf_p->MinusArg,2);
  }

}
void dwf_end( Dwf *dwf_p)
{
  char *cname = " ";
  char *fname = "dwf_end(Dwf*)";
  VRB.Func(cname,fname);
int ls = GJP.Snodes();
  //------------------------------------------------------------------
  // Free memory of the 12 word communications buffer needed
  // for the spread-out case.
  //------------------------------------------------------------------
 
  VRB.Sfree(cname,fname, "comm_buf", dwf_p->comm_buf);
  ffree(dwf_p->comm_buf);
 
  //----------------------------------------------------------------
  // Free temporary femrion fields 2 and 1
  //----------------------------------------------------------------
  VRB.Sfree(cname,fname, "frm_tmp2", dwf_p->frm_tmp2);
  sfree(dwf_p->frm_tmp2);

  VRB.Sfree(cname,fname, "frm_tmp1", dwf_p->frm_tmp1);
  sfree(dwf_p->frm_tmp1);

 if(ls>1){
  for(int i = 0;i<2;i++){
    delete (dwf_p->PlusArg[i]);
    delete (dwf_p->MinusArg[i]);
  }
  delete (dwf_p->Plus);
  delete (dwf_p->Minus);
 }

  //----------------------------------------------------------------
  // Un-initialize the wilson library. Memory is set free here.
  //----------------------------------------------------------------
  wfm_vec_end(&wil);
  wilson_compat_end(dwf_p->wilson_p);

}


CPS_END_NAMESPACE
