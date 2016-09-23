#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief Definition of glb_sum routine.

  $Id: glb_sum.C,v 1.4 2008/02/08 18:35:06 chulwoo Exp $ 
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008/02/08 18:35:06 $
//  $Header: /space/cvs/cps/cps++/src/comms/mpi/glb/glb_sum.C,v 1.4 2008/02/08 18:35:06 chulwoo Exp $
//  $Id: glb_sum.C,v 1.4 2008/02/08 18:35:06 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: glb_sum.C,v $
//  $Revision: 1.4 $
//  $Source: /space/cvs/cps/cps++/src/comms/mpi/glb/glb_sum.C,v $
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
CPS_START_NAMESPACE



static Double64 transmit_buf;
static Double64 receive_buf;
static Double64 gsum_buf;


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

  gsum_buf = *float_p;

  for(int i = 0; i < 4; ++i) {

      transmit_buf = gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
// If we are using the non-Tartan setup then GLOBALSUM_TYPE will be properly
// defined as a float or double native type:
#ifdef GLOBALSUM_TYPE
	SCUDirArg send(&transmit_buf, gjp_scu_dir[2*i], SCU_SEND, (sizeof(Double64)/COMMS_DATASIZE) );
	SCUDirArg rcv(&receive_buf, gjp_scu_dir[2*i+1], SCU_REC, (sizeof(Double64)/COMMS_DATASIZE) );
// If not, then sizeof(Double64) will break because it returns the pointer
// size (?, 4 bytes instead of 8??).
// FIXME:  There should be a more elegant solution to this problem.
#else
	SCUDirArg send(&transmit_buf, gjp_scu_dir[2*i], SCU_SEND, 2 );
	SCUDirArg rcv(&receive_buf, gjp_scu_dir[2*i+1], SCU_REC, 2 );
#endif

	SCUTrans(&send);
	SCUTrans(&rcv);

	SCUTransComplete();

        gsum_buf += receive_buf;
        transmit_buf = receive_buf;
      }
  }
  *float_p = gsum_buf;
}
CPS_END_NAMESPACE
