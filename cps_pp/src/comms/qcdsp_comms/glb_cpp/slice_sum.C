#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/glb_cpp/slice_sum.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: slice_sum.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.8  2002/03/11 22:26:36  anj
//  This should now be the correct, fully merged code from our two versions. Anj
//
//  Revision 1.5.2.1  2002/03/08 16:36:02  anj
//  Checking in the Columbia code branch on tag Columbia4_1_1_test-branch, to be
//  merged with the UKQCD head branch shortly.  Anj
//
//  Revision 1.5  2001/08/16 12:54:02  anj
//  Some fixes follosin the float-> IFloat change, mostly of the (variable
//  anme) IFloat_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.4  2001/08/16 10:49:48  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
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
//  $RCSfile: slice_sum.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/glb_cpp/slice_sum.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//====================================================================
//*  SUI 3/27/97
//*  slice_sum.C
//*  last modified 11/7/97
//====================================================================


CPS_END_NAMESPACE
#include<comms/glb.h>
#include<util/smalloc.h>
#include<util/gjp.h>
#include <sysfunc.h>
CPS_START_NAMESPACE

const SCUDir pos_dir[] = { SCU_XP, SCU_YP, SCU_ZP, SCU_TP };
const SCUDir neg_dir[] = { SCU_XM, SCU_YM, SCU_ZM, SCU_TM };


//-------------------------------------------------------------------
//* sum over a slice(hyperplane) which is orthogonal to the direction
//* "dir". There are "blcklength" summations to be done:
//* float_p[i] = sum_over_nodes_of_this_slice(float_p[i])
//-------------------------------------------------------------------

void slice_sum(Float * float_p, int blcklength, int dir)
{
  char *cname = "slice_sum";
  char *fname = "slice_sum(*float_p, int, int)";

  int NP[4] = { GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes() };
  const int MAX=1023;

  if (blcklength > MAX)
    ERR.General(cname, fname, "blcklength (%d) too big > MAX (%d) \n",blcklength, MAX);
		    
  IFloat *transmit_buf = (IFloat *)smalloc(blcklength*sizeof(IFloat));

  // added by manke
  if(transmit_buf == 0)
    ERR.Pointer(cname,fname, "transmit_buf");
  VRB.Smalloc(cname,fname, "transmit_buf", transmit_buf, blcklength*sizeof(float));
  // end added

  IFloat *receive_buf  = (IFloat *)smalloc(blcklength*sizeof(IFloat));

  // added by manke
  if(receive_buf == 0)
    ERR.Pointer(cname,fname, "receive_buf");
  VRB.Smalloc(cname,fname, "receive_buf", receive_buf, blcklength*sizeof(float));
  // end added

  IFloat *transmit_buf_p = transmit_buf;
  IFloat *receive_buf_p = receive_buf;
  
  IFloat *free_buffer_p = transmit_buf_p; 
	// buffer to be used as next receive 
  int count;		  
	// loop index where 0 <= count < blcklength 
  int itmp;		  
	// loop index with 1<= itmp < NP[i] 
  int i;

  for(i = 0; i < 4; ++i) {

      if(i == dir) continue;

      //--------------------------------------------------------------
      // address of buffer of to be sent (data on this node) 
      //--------------------------------------------------------------
      transmit_buf_p = (IFloat *)float_p; 

      SCUDirArg send(transmit_buf_p, pos_dir[i], SCU_SEND, blcklength);
      SCUDirArg recv(receive_buf_p, neg_dir[i], SCU_REC, blcklength);

      //--------------------------------------------------------------
      // tranmit & receive NP[i] - 1 times in snd_dir[i] direction
      //--------------------------------------------------------------
      for ( itmp = 1; itmp < NP[i]; itmp++) {

	 //-----------------------------------------------------------
         // do SCU transfers
	 //-----------------------------------------------------------
         SCUTrans(&send);
         SCUTrans(&recv);
	 SCUTransComplete();

	 //-----------------------------------------------------------
         // accumulate received data	
	 //-----------------------------------------------------------
   	 for (count = 0; count < blcklength; ++count) 
           float_p[count] += receive_buf_p[count];

	 //-----------------------------------------------------------
         // the received data will be sent out     
	 // the free buffer will be used to receive      
	 //-----------------------------------------------------------
	 send.Addr(transmit_buf_p = receive_buf_p);
	 recv.Addr(receive_buf_p = free_buffer_p);

	 //-----------------------------------------------------------
	 // transmit_buf WILL be free buffer NEXT round of transmit
	 //-----------------------------------------------------------
	 free_buffer_p = transmit_buf_p;
      }
  }
  sfree(transmit_buf);
  sfree(receive_buf);
}


CPS_END_NAMESPACE
