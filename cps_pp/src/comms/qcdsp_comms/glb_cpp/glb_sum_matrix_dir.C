#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_sum_matrix_dir routine.

  $Id: glb_sum_matrix_dir.C,v 1.2 2003-07-24 16:53:54 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:54 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/glb_cpp/glb_sum_matrix_dir.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Id: glb_sum_matrix_dir.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.1.1.1  2003/06/22 13:34:47  mcneile
//  This is the cleaned up version of the Columbia Physics System.
//  The directory structure has been changed.
//  The include paths have been updated.
//
//
//  Revision 1.2  2001/06/19 18:11:51  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:01  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: glb_sum_matrix_dir.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/glb_cpp/glb_sum_matrix_dir.C,v $
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
#include <comms/sysfunc.h>
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
