#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdsp_nos_hdw_diag/dwf_end.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: dwf_end.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:12:42  anj
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
//  $RCSfile: dwf_end.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdsp_nos_hdw_diag/dwf_end.C,v $
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
