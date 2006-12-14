#include<config.h>
CPS_START_NAMESPACE
/*--------------------------------------------------------------------*/

CPS_END_NAMESPACE
#include<util/wilson.h>
#include <comms/scu.h>
#if TARGET_BGL == 1
#include <bgl_sys/bgl_sys_all.h>
#endif
CPS_START_NAMESPACE

void wfm_comm_forward(IFloat *ab0, IFloat *ab1, IFloat *ab2, IFloat *ab3,
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
      for(i=0; i<wilson_p->comm_numblk[mu]; ++i) { 

	for(j=0; j<wilson_p->comm_blklen[mu]; ++j) {
	  getMinusData(receive_ad, send_ad, BLOCK, mu);
	  send_ad    = send_ad    + BLOCK;
	  receive_ad = receive_ad + BLOCK;
	}
	 
	send_ad = send_ad + wilson_p->comm_stride[mu] - 1;
	receive_ad = receive_ad + wilson_p->comm_stride[mu] - 1;
      }
   }
}
CPS_END_NAMESPACE
