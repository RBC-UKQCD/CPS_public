#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/glb_cpp_25MHz_nos/glb_sum_matrix_dir.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: glb_sum_matrix_dir.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:49:50  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:11:55  anj
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
//  Revision 1.2  2001/05/25 06:16:02  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: glb_sum_matrix_dir.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/glb_cpp_25MHz_nos/glb_sum_matrix_dir.C,v $
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
#include <sysfunc.h>
CPS_START_NAMESPACE



static Matrix transmit_buf;
static Matrix receive_buf;
static Matrix gsum_buf;
static IFloat *send_buf ;
static IFloat *rcv_buf ;

static volatile unsigned* dsp_scu_base0x10 =
        (volatile unsigned* )(DSP_SCU_BASE + 0x10);



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

  send_buf = (IFloat *) &transmit_buf;
  rcv_buf = (IFloat *) &receive_buf;

  int blocksize=sizeof(Float)*sizeof(Matrix);
  int itmp;
  for (itmp = 1; itmp < NP[dir]; itmp++) {
    bsm(send_buf, blocksize, 0, 1, gjp_scu_wire_map[2*dir], TRANSMIT);
    bsm(rcv_buf, blocksize, 0, 1, gjp_scu_wire_map[2*dir+1], RECEIVE);
    while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*dir+1]) ) ;
    while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*dir]) ) ;

 
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
     bsm(send_buf, blocksize, 0, 1, gjp_scu_wire_map[2*dir], TRANSMIT);
    bsm(rcv_buf, blocksize, 0, 1, gjp_scu_wire_map[2*dir+1], RECEIVE);
    while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*dir+1]) ) ;
    while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*dir]) ) ;
 
    gsum_buf += receive_buf;
    transmit_buf = receive_buf;
  }


  *matrix_p = gsum_buf;
}


CPS_END_NAMESPACE
