#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_copy_backward.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include<util/wilson.h>
CPS_START_NAMESPACE

extern int gjp_local_axis[];

extern "C" void wfm_copy_backward(IFloat *af0, IFloat *af1, IFloat *af2, IFloat *af3, 
				  Wilson *wilson_p)
{
   int   mu, i, j;
   IFloat *receive_ad, *send_ad;
   IFloat *af[4];

   af[0] = af0;
   af[1] = af1;
   af[2] = af2;
   af[3] = af3;

   for(mu=0; mu<ND; ++mu)
   {
     if(gjp_local_axis[mu] == 1) {
       receive_ad = wilson_p->comm_offset[mu] + af[mu];
       send_ad = 0 + af[mu];
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
