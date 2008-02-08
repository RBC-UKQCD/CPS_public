#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_sum_matrix_dir routine.

  $Id: glb_sum_matrix_dir.C,v 1.4 2008-02-08 18:35:06 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008-02-08 18:35:06 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/mpi/glb/glb_sum_matrix_dir.C,v 1.4 2008-02-08 18:35:06 chulwoo Exp $
//  $Id: glb_sum_matrix_dir.C,v 1.4 2008-02-08 18:35:06 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: glb_sum_matrix_dir.C,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/mpi/glb/glb_sum_matrix_dir.C,v $
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



static Matrix transmit_buf;
static Matrix receive_buf;
static Matrix gsum_buf;

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

  // Sum along dir
  //--------------------------------------------------------------
  gsum_buf = *matrix_p;

  transmit_buf = gsum_buf;
  int blocksize=sizeof(Float)*sizeof(Matrix);
  int itmp;
  for (itmp = 1; itmp < NP[dir]; itmp++) {
    SCUDirArg send(&transmit_buf, gjp_scu_dir[2*dir], SCU_SEND, blocksize);
    SCUDirArg rcv(&receive_buf, gjp_scu_dir[2*dir+1], SCU_REC, blocksize);

    SCUTrans(&send);
    SCUTrans(&rcv);

    SCUTransComplete();

    gsum_buf += receive_buf;
    transmit_buf = receive_buf;
  }


  // Broadcast the result of node with dir coordinate == 0
  //--------------------------------------------------------------

  if(COOR[dir] != 0) {
    gsum_buf.ZeroMatrix();
  }
    
  transmit_buf = gsum_buf;
  
  for (itmp = 1; itmp < NP[dir]; itmp++) {
    SCUDirArg send(&transmit_buf, gjp_scu_dir[2*dir], SCU_SEND, blocksize);
    SCUDirArg rcv(&receive_buf, gjp_scu_dir[2*dir+1], SCU_REC, blocksize);
    
    SCUTrans(&send);
    SCUTrans(&rcv);
    
    SCUTransComplete();
    
    gsum_buf += receive_buf;
    transmit_buf = receive_buf;
  }


  *matrix_p = gsum_buf;
}


CPS_END_NAMESPACE
