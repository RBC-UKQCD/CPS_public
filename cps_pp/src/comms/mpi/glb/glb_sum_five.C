#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_sum_five routine.

  $Id: glb_sum_five.C,v 1.4 2008/02/08 18:35:06 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008/02/08 18:35:06 $
//  $Header: /space/cvs/cps/cps++/src/comms/mpi/glb/glb_sum_five.C,v 1.4 2008/02/08 18:35:06 chulwoo Exp $
//  $Id: glb_sum_five.C,v 1.4 2008/02/08 18:35:06 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: glb_sum_five.C,v $
//  $Revision: 1.4 $
//  $Source: /space/cvs/cps/cps++/src/comms/mpi/glb/glb_sum_five.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
//  glb_sum_five
//
// Sum over all nodes of the "virtual" 5-dimensional volume.
// {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes(), GJP.Snodes}
// Relevant for spread-out DWF (GJP.s_nodes not 1) only.
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
  This routine need only be used by domain-wall fermion code where
  the 5th dimension is parallelised.
  
  \param float_p The number to be summed.
  \post The number pointed to by \a float_p is summed over all nodes
  and that sum is written back to \a float_p, which is identical on all nodes.

  \ingroup comms
*/
//---------------------------------------------------------------------- 

void glb_sum_five(Float * float_p)
{
  int NP[5] = {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes(), GJP.Snodes()};


  // Sum over the "virtual" 5-dimensional mesh
  //------------------------------------------------------------
  gsum_buf = *float_p;

  int i;
  for(i = 0; i < 5; ++i) {

      transmit_buf = gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
	SCUDirArg send(&transmit_buf, gjp_scu_dir[2*i], SCU_SEND, 2);
	SCUDirArg rcv(&receive_buf, gjp_scu_dir[2*i+1], SCU_REC, 2);

	SCUTrans(&send);
	SCUTrans(&rcv);

	SCUTransComplete();

        gsum_buf += receive_buf;
        transmit_buf = receive_buf;
      }
  }

  // Broadcast the result of node (0,0,0,0,0)
  //------------------------------------------------------------
  if(GJP.XnodeCoor() != 0 || 
     GJP.YnodeCoor() != 0 || 
     GJP.ZnodeCoor() != 0 || 
     GJP.TnodeCoor() != 0 || 
     GJP.SnodeCoor() != 0 ) {
    gsum_buf = 0;
  }

  for(i = 0; i < 5; ++i) {
    
      transmit_buf = gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
	SCUDirArg send(&transmit_buf, gjp_scu_dir[2*i], SCU_SEND, 2);
	SCUDirArg rcv(&receive_buf, gjp_scu_dir[2*i+1], SCU_REC, 2);

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
