#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdsp_nos/dwf_init.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: dwf_init.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.6  2001/08/16 10:50:18  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.4  2001/07/03 17:01:01  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.3  2001/06/21 15:40:17  anj
//  Updated the _TARTAN ifdefs, using PARALLEL instead (where appropriate).Anj
//
//  Revision 1.2  2001/06/19 18:12:40  anj
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
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdsp_nos/dwf_init.C,v $
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
#include<config.h>
CPS_START_NAMESPACE
#ifdef PARALLEL
CPS_END_NAMESPACE
#include <sysfunc.h>
CPS_START_NAMESPACE
#endif
CPS_END_NAMESPACE
#include<util/dwf.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
CPS_START_NAMESPACE

int dwfso_wire_map[2];


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

#ifdef PARALLEL
//------------------------------------------------------------------
// Set the communication wires 
//------------------------------------------------------------------
  dwfso_wire_map[0] = SCURemap( gjp_scu_dir[8] );
  dwfso_wire_map[1] = SCURemap( gjp_scu_dir[9] );
#endif

}






CPS_END_NAMESPACE
