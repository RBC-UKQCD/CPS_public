#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  gamma_5 multiplication routine

  Used by derivatives of the FwilsonTypes class.

  $Id: gamma_5.C,v 1.2 2003-07-24 16:53:54 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:54 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_wilson_types/noarch/gamma_5.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Id: gamma_5.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:36  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:25  anj
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
//  Revision 1.2  2001/05/25 06:16:09  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: gamma_5.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_wilson_types/noarch/gamma_5.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
// gamma_5(IFloat *v_out, IFloat *v_in, int num_sites):
//
// v_out = Gamma5 * v_in. Gamme5 is in the chiral basis
//
//          [ 1  0  0  0]
// Gamma5 = [ 0  1  0  0]
//          [ 0  0 -1  0]
//          [ 0  0  0 -1]
//
// num_sites is the number of sites. It is assumed
// that each site has 24 components.
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/data_types.h>
CPS_START_NAMESPACE

void gamma_5(IFloat *v_out, IFloat *v_in, int num_sites) ;

void gamma_5(IFloat *v_out, IFloat *v_in, int num_sites) 
{
  IFloat *p_out = v_out ;
  IFloat *p_in = v_in ;

  int half_site_size = 12 ;

  for (int site=0; site<num_sites; ++site) {
    int comp;
    for (comp=0; comp<half_site_size; ++comp) {
      *p_out++ = *p_in++ ;
    }
    for (comp=0; comp<half_site_size; ++comp) {
      *p_out++ = -*p_in++ ;
    }
  }
}

CPS_END_NAMESPACE
