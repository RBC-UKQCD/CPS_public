#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-01-13 20:39:01 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_s_spect/aots_s.C,v 1.3 2004-01-13 20:39:01 chulwoo Exp $
//  $Id: aots_s.C,v 1.3 2004-01-13 20:39:01 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2.10.1  2003/11/06 00:17:43  cwj
//  *** empty log message ***
//
//  Revision 1.1.1.1  2003/11/04 05:04:58  chulwoo
//
//  starting again
//
//
//  Revision 1.2  2003/07/24 16:53:53  zs
//  Addition of documentation via doxygen:
//  doxygen-parsable comment blocks added to many source files;
//  New target in makefile and consequent alterations to configure.in;
//  New directories and files under the doc directory.
//
//  Revision 1.2  2001/06/19 18:11:29  anj
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
//  Revision 1.2  2001/05/25 06:16:00  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: aots_s.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_s_spect/aots_s.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// aots_s.C
CPS_END_NAMESPACE
#include <alg/aots_s.h>
#include <util/gjp.h>
#include <util/verbose.h>
CPS_START_NAMESPACE

//-----------------------------------------------------------------
// static member initialization
//-----------------------------------------------------------------
char Aots::cname[] = "Aots";

//-----------------------------------------------------------------
// CTOR
//-----------------------------------------------------------------
Aots::Aots(int s, int e, int stride)
: start(s), end(e), step(stride), current (s)
{
  char *fname = "Aots(int, int, int)";
  VRB.Func(cname, fname);

  if ( (step <= 0) || (start < 0) || (end < 0) || (start > end) || 
       (start >= GJP.TnodeSites() * GJP.Tnodes()) || 
       (end >= GJP.TnodeSites() * GJP.Tnodes()))
    ERR.General(cname, fname, "Invalid Aots Parameters\n");
}

CPS_END_NAMESPACE
