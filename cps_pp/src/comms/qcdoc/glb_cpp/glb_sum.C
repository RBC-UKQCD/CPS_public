#include<config.h>
#include<stdio.h>
#include<qalloc.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief Definition of glb_sum routine.

  $Id: glb_sum.C,v 1.3 2004-04-27 03:51:17 cwj Exp $ 
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: cwj $
//  $Date: 2004-04-27 03:51:17 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdoc/glb_cpp/glb_sum.C,v 1.3 2004-04-27 03:51:17 cwj Exp $
//  $Id: glb_sum.C,v 1.3 2004-04-27 03:51:17 cwj Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: glb_sum.C,v $
//  $Revision: 1.3 $
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
#include <comms/sysfunc.h>
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
  int NP[4] = {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()};
  static int counter = 0;
  if (transmit_buf == NULL) 
      transmit_buf = (Double64 *)qalloc(QFAST|QNONCACHE,sizeof(Double64));
  if (receive_buf == NULL) 
      receive_buf = (Double64 *)qalloc(QFAST|QNONCACHE,sizeof(Double64));
  if (gsum_buf == NULL) 
      gsum_buf = (Double64 *)qalloc(QFAST|QNONCACHE,sizeof(Double64));

  if (output) printf("glb_sum %d before = %e ", counter, (double)*float_p);
  *gsum_buf = (Double64)*float_p;

  for(int i = 0; i < 4; ++i) {

      *transmit_buf = *gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
	SCUDirArg send(transmit_buf, gjp_scu_dir[2*i], SCU_SEND, sizeof(Double64) );
	SCUDirArg rcv(receive_buf, gjp_scu_dir[2*i+1], SCU_REC, sizeof(Double64) );

//	printf("sent=%e received %e\n",transmit_buf,receive_buf);
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
