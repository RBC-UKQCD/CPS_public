#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.

  $Id: wilson_end.C,v 1.2 2003-07-24 16:53:54 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:54 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/noarch/wilson_end.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Id: wilson_end.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:12:47  anj
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
//  Revision 1.2  2001/05/25 06:16:06  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: wilson_end.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/noarch/wilson_end.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/wilson.h>
#include <util/smalloc.h>
#include <util/verbose.h>
CPS_START_NAMESPACE

  
void wilson_end( Wilson *wilson_p)
{
  char *cname = " ";
  char *fname = "wilson_end(Wilson*)";
  VRB.Func(cname,fname);

  VRB.Sfree(cname,fname, "af[1]", wilson_p->af[1]);
  sfree(wilson_p->af[1]);

  VRB.Sfree(cname,fname, "af[0]", wilson_p->af[0]);
  sfree(wilson_p->af[0]);

  VRB.Sfree(cname,fname, "ptr", wilson_p->ptr);
  sfree(wilson_p->ptr);

}

CPS_END_NAMESPACE
