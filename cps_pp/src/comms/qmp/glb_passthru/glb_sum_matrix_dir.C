#include<config.h>
#include<qalloc.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_sum_matrix_dir routine.
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qmp/glb_passthru/glb_sum_matrix_dir.C,v $
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
#include<util/vector.h>
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE



static Matrix *transmit_buf = NULL;
static Matrix *receive_buf = NULL;
static Matrix *gsum_buf = NULL;

//----------------------------------------------------------------------
/*!
  \param matrix_p The Matrix to be summed.
  \param dir The direction in which to sum; one of {0, 1, 2, 3, 4},
  corresponding to {x, y, z, t, s}.
  \post The matrix pointed to by \a float_p is summed over all nodes along the
  \a dir direction, \e i.e. over each strip of nodes where the grid
  coordinates in all other directions are constant.
  The matrix sum is written back to \a matrix_p, which is identical on all
  nodes in this strip.

  \ingroup comms
*/
//---------------------------------------------------------------------- 

void glb_sum_matrix_dir(Matrix * matrix_p, int dir)
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
      transmit_buf = (Matrix *)qalloc(QFAST|QNONCACHE,sizeof(Matrix));
  if (receive_buf == NULL) 
      receive_buf = (Matrix *)qalloc(QFAST|QNONCACHE,sizeof(Matrix));
  if (gsum_buf == NULL) 
      gsum_buf = (Matrix *)qalloc(QFAST|QNONCACHE,sizeof(Matrix));

  // Sum along dir
  //--------------------------------------------------------------
  *gsum_buf = *matrix_p;

  *transmit_buf = *gsum_buf;
  int blocksize=sizeof(Matrix);
  int itmp;
  for (itmp = 1; itmp < NP[dir]; itmp++) {
    SCUDirArg send(transmit_buf, gjp_scu_dir[2*dir], SCU_SEND, blocksize);
    SCUDirArg rcv(receive_buf, gjp_scu_dir[2*dir+1], SCU_REC, blocksize);

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
    gsum_buf->ZeroMatrix();
  }
    
  *transmit_buf = *gsum_buf;
  
  for (itmp = 1; itmp < NP[dir]; itmp++) {
    SCUDirArg send(transmit_buf, gjp_scu_dir[2*dir], SCU_SEND, blocksize);
    SCUDirArg rcv(receive_buf, gjp_scu_dir[2*dir+1], SCU_REC, blocksize);
    
    send.StartTrans();
    rcv.StartTrans();
    send.TransComplete();
    rcv.TransComplete();
    
    *gsum_buf += *receive_buf;
    *transmit_buf = *receive_buf;
  }


  *matrix_p = *gsum_buf;
}


CPS_END_NAMESPACE
