#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_lcl_nos/wfm_comm_forward_l.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: wfm_comm_forward_l.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:24  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:12:56  anj
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
//  Revision 1.2  2001/05/25 06:16:07  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: wfm_comm_forward_l.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_lcl_nos/wfm_comm_forward_l.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include<util/wilson.h>
CPS_START_NAMESPACE

extern "C" void wfm_comm_forward_l(IFloat *ab0, IFloat *ab1, IFloat *ab2, IFloat *ab3,
			     Wilson *wilson_p)
{
   int   mu, i, j;
   IFloat *send_ad, *receive_ad;
   IFloat *ab[4];

   ab[0]=ab0;
   ab[1]=ab1;
   ab[2]=ab2;
   ab[3]=ab3;

   for(mu=0; mu<ND; ++mu)
   {
      send_ad = wilson_p->comm_offset[mu] + ab[mu];
      receive_ad = 0 + ab[mu];
      for(i=0; i<wilson_p->comm_numblk[mu]; ++i)
      { 
	 for(j=0; j<wilson_p->comm_blklen[mu]; ++j)
	    *receive_ad++ = *send_ad++;

	 send_ad = send_ad + wilson_p->comm_stride[mu] - 1;
	 receive_ad = receive_ad + wilson_p->comm_stride[mu] - 1;
      }
   }
}
CPS_END_NAMESPACE
