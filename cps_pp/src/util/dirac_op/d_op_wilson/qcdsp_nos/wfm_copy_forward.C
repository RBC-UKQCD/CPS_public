#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:09 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_copy_forward.C,v 1.2 2004-06-04 21:14:09 chulwoo Exp $
//  $Id: wfm_copy_forward.C,v 1.2 2004-06-04 21:14:09 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_copy_forward.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include<util/wilson.h>
CPS_START_NAMESPACE

extern int gjp_local_axis[];

extern "C" void wfm_copy_forward(IFloat *ab0, IFloat *ab1, IFloat *ab2, IFloat *ab3,
			         Wilson *wilson_p)
{
   int   mu, i, j;
   IFloat *send_ad, *receive_ad;
   IFloat *ab[4];

   ab[0]=ab0;
   ab[1]=ab1;
   ab[2]=ab2;
   ab[3]=ab3;

   for(mu=0; mu<ND; ++mu) {
     if(gjp_local_axis[mu] == 1) {
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
}
CPS_END_NAMESPACE
