#include<config.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include<util/wilson.h>
#include <comms/scu.h>
#if TARGET_BGL == 1
#include <bgl_sys/bgl_sys_all.h>
#endif
CPS_START_NAMESPACE

void wfm_comm_backward(IFloat *af0, IFloat *af1, IFloat *af2, IFloat *af3, 
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
      receive_ad = wilson_p->comm_offset[mu] + af[mu];
      send_ad = 0 + af[mu];
      for(i=0; i<wilson_p->comm_numblk[mu]; ++i) { 

	for(j=0; j<wilson_p->comm_blklen[mu]; ++j) {
	  getPlusData(receive_ad, send_ad, BLOCK, mu);
	  send_ad    = send_ad    + BLOCK;
	  receive_ad = receive_ad + BLOCK;
	}

	send_ad = send_ad + wilson_p->comm_stride[mu] - 1;
	receive_ad = receive_ad + wilson_p->comm_stride[mu] - 1;
      }
   }
}
CPS_END_NAMESPACE
