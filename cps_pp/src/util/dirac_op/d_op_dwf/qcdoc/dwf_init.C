#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-07-01 21:22:36 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdoc/dwf_init.C,v 1.2 2004-07-01 21:22:36 chulwoo Exp $
//  $Id: dwf_init.C,v 1.2 2004-07-01 21:22:36 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.1.4.1  2004/06/09 04:31:54  chulwoo
//  *** empty log message ***
//
//  Revision 1.1.2.1  2004/05/20 14:36:32  pab
//  Files for three optimised dirac operators.
//  Added some patches to the ReadLattice too.
//
//  Revision 1.1.1.1  2003/06/22 13:34:46  mcneile
//  This is the cleaned up version of the Columbia Physics System.
//  The directory structure has been changed.
//  The include paths have been updated.
//
//
//  Revision 1.4  2001/08/16 10:50:17  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:12:39  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:05  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: dwf_init.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdoc/dwf_init.C,v $
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
  wilson_compat_init(dwf_p->wilson_p);

//------------------------------------------------------------------
// Allocate memory for two temporary fermion checkerboard fields  
//------------------------------------------------------------------
  int f_size = 24 * GJP.VolNodeSites() * GJP.SnodeSites() / 2; 

  dwf_p->frm_tmp1 = (IFloat *) smalloc(f_size*sizeof(IFloat));
  if(dwf_p->frm_tmp1 == 0)
    ERR.Pointer(cname,fname, "frm_tmp1");
  VRB.Smalloc(cname,fname,
	      "frm_tmp1", dwf_p->frm_tmp1, f_size*sizeof(IFloat));

//  for(int i=0;i<f_size;i++) (dwf_p->frm_tmp1)[i]=0.0;

  dwf_p->frm_tmp2 = (IFloat *) smalloc(f_size*sizeof(IFloat));
  if(dwf_p->frm_tmp2 == 0)
    ERR.Pointer(cname,fname, "frm_tmp2");
  VRB.Smalloc(cname,fname,
	      "frm_tmp2", dwf_p->frm_tmp2, f_size*sizeof(IFloat));
//  for(int i=0;i<f_size;i++) (dwf_p->frm_tmp2)[i]=0.0;

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

  int ls = dwf_p->ls;

}
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
  wfm_vec_end(&wil);
  wilson_compat_end(dwf_p->wilson_p);

}


CPS_END_NAMESPACE
