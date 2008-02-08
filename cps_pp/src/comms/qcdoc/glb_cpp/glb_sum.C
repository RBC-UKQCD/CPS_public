#include<config.h>
#include<util/qcdio.h>
#include<qalloc.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief Definition of glb_sum routine.
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdoc/glb_cpp/glb_sum.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
// glb_sum
//
// Sum over all nodes 
// {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()}
//--------------------------------------------------------------
CPS_END_NAMESPACE
#include<comms/glb.h>
#include<comms/scu.h>
#include<util/gjp.h>
#include<comms/double64.h>
#include <comms/sysfunc_cps.h>
#include <comms/glb_sum_internal.h>
CPS_START_NAMESPACE



static Double64 *transmit_buf = NULL;
static Double64 *receive_buf = NULL;
static Double64 *gsum_buf = NULL;

static int output = 0;


//----------------------------------------------------------------------
/*!
  \param float_p The number to be summed.
  \post The number pointed to by \a float_p is summed over all nodes
  and that sum is written back to \a float_p, which is identical on all nodes.

  \ingroup comms
*/
//---------------------------------------------------------------------- 
void glb_sum(Float * float_p)
{
  glb_sum_internal2(float_p,4);
}

void glb_sum_gimp(Float * float_p)
{
  int NP[4] = {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()};
  static int counter = 0;
  if (transmit_buf == NULL) 
      transmit_buf = (Double64 *)qalloc(QFAST|QNONCACHE,sizeof(Double64));
  if (receive_buf == NULL) 
      receive_buf = (Double64 *)qalloc(QFAST|QNONCACHE,sizeof(Double64));
  if (gsum_buf == NULL) 
      gsum_buf = (Double64 *)qalloc(QFAST|QNONCACHE,sizeof(Double64));
  if (output) printf("glb_sum cpp %d before = %e ", counter, (double)*float_p);
  *gsum_buf = (Double64)*float_p;

  for(int i = 0; i < 4; ++i) {

      *transmit_buf = *gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
	SCUDirArg send(transmit_buf, gjp_scu_dir[2*i], SCU_SEND, sizeof(Double64) );
	SCUDirArg rcv(receive_buf, gjp_scu_dir[2*i+1], SCU_REC, sizeof(Double64) );

	send.StartTrans();
	rcv.StartTrans();
	send.TransComplete();
	rcv.TransComplete();

        *gsum_buf += *receive_buf;
        *transmit_buf = *receive_buf;
      }
  }
  *float_p = (Float)*gsum_buf;
  if (output)   printf("after = %e\n", (double)*float_p);
  counter++;
}

CPS_END_NAMESPACE
