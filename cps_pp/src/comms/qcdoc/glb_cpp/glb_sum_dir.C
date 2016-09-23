#include<config.h>
#include<qalloc.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_sum_dir routine.

*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /space/cvs/cps/cps++/src/comms/qcdoc/glb_cpp/glb_sum_dir.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
// glb_sum_dir
//
// Sum over all nodes along a direction
// (0,1,2,3,4) <-> (x,y,z,t,s)
//--------------------------------------------------------------
CPS_END_NAMESPACE
#include<comms/glb.h>
#include<comms/scu.h>
#include<util/gjp.h>
#include<comms/double64.h>
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE


static Double64 *transmit_buf = NULL;
static Double64 *receive_buf = NULL;
static Double64 *gsum_buf = NULL;

//----------------------------------------------------------------------
/*!
  \param float_p The number to be summed.
  \param dir The direction in which to sum; one of {0, 1, 2, 3, 4},
  corresponding to {x, y, z, t, s}.
  \post The number pointed to by \a float_p is summed over all nodes along the
  \a dir direction, \e i.e. over each strip of nodes where the grid
  coordinates in all other directions are constant.
  That sum is written back to \a float_p, which is identical on all nodes in
  this strip.

  \ingroup comms
*/
//---------------------------------------------------------------------- 

void glb_sum_dir(Float * float_p, int dir)
{
  int NP[5] = {GJP.Xnodes(), 
	       GJP.Ynodes(), 
	       GJP.Znodes(), 
	       GJP.Tnodes(), 
	       GJP.Snodes()};

  int COOR[5] = {GJP.XnodeCoor(), 
		 GJP.YnodeCoor(), 
		 GJP.ZnodeCoor(), 
		 GJP.TnodeCoor(), 
		 GJP.SnodeCoor()}; 


  if (transmit_buf == NULL) 
      transmit_buf = (Double64 *)qalloc(QFAST|QNONCACHE,sizeof(Double64));
  if (receive_buf == NULL) 
      receive_buf = (Double64 *)qalloc(QFAST|QNONCACHE,sizeof(Double64));
  if (gsum_buf == NULL) 
      gsum_buf = (Double64 *)qalloc(QFAST|QNONCACHE,sizeof(Double64));

  // Sum along dir
  //--------------------------------------------------------------
  *gsum_buf = *float_p;

  *transmit_buf = *gsum_buf;

  int itmp;
  for (itmp = 1; itmp < NP[dir]; itmp++) {
    SCUDirArg send(transmit_buf, gjp_scu_dir[2*dir], SCU_SEND, sizeof(Double64));
    SCUDirArg rcv(receive_buf, gjp_scu_dir[2*dir+1], SCU_REC, sizeof(Double64));
    send.StartTrans();
    rcv.StartTrans();
    send.TransComplete();
    rcv.TransComplete();

    *gsum_buf += *receive_buf;
    *transmit_buf = *receive_buf;
  }


  // Broadcast the result of node with dir coordinate == 0
  //--------------------------------------------------------------

  if(COOR[dir] != 0) {
    *gsum_buf = 0;
  }
    
  *transmit_buf = *gsum_buf;
  
  for (itmp = 1; itmp < NP[dir]; itmp++) {
    SCUDirArg send(transmit_buf, gjp_scu_dir[2*dir], SCU_SEND, sizeof(Double64));
    SCUDirArg rcv(receive_buf, gjp_scu_dir[2*dir+1], SCU_REC, sizeof(Double64));
    send.StartTrans();
    rcv.StartTrans();
    send.TransComplete();
    rcv.TransComplete();
    
    *gsum_buf += *receive_buf;
    *transmit_buf = *receive_buf;
  }


  *float_p = *gsum_buf;
//  fprintf(stdout,"glb_sum_dir():NP[%d]=%d COOR[%d]=%d sum=%0.16e \n",dir,NP[dir],dir,COOR[dir],float_p);
}


CPS_END_NAMESPACE
