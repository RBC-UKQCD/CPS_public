#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief Definition of slice_sum routine

  $Id: slice_sum.C,v 1.4 2008-02-08 18:35:06 chulwoo Exp $
 */
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008-02-08 18:35:06 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/mpi/glb/slice_sum.C,v 1.4 2008-02-08 18:35:06 chulwoo Exp $
//  $Id: slice_sum.C,v 1.4 2008-02-08 18:35:06 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: slice_sum.C,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/mpi/glb/slice_sum.C,v $
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
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE

const SCUDir pos_dir[] = { SCU_XP, SCU_YP, SCU_ZP, SCU_TP };
const SCUDir neg_dir[] = { SCU_XM, SCU_YM, SCU_ZM, SCU_TM };


//-------------------------------------------------------------------
//* sum over a slice(hyperplane) which is orthogonal to the direction
//* "dir". There are "blcklength" summations to be done:
//* float_p[i] = sum_over_nodes_of_this_slice(float_p[i])
//-------------------------------------------------------------------

//-------------------------------------------------------------------
/*!
  The vector pointed to by \a float_p is summed over each 3-dimensional
  hyperplane of nodes which is perpendiculat to the \a dir direction

  \param float_p The number to be summed.
  \param blcklength The number of floating point numbers in the vector.
  \param dir The normal direction defining the hyperplane; one of {0, 1, 2, 3} 
  corresponding to {x, y, z, t}.
  \post The vector sum is written back to \a float_p, which is identical on
  all nodes in this hyperplane.

  \ingroup comms
*/
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
  VRB.Smalloc(cname,fname, "transmit_buf", transmit_buf, blcklength*sizeof(IFloat));
  // end added

  IFloat *receive_buf  = (IFloat *)smalloc(blcklength*sizeof(IFloat));

  // added by manke
  if(receive_buf == 0)
    ERR.Pointer(cname,fname, "receive_buf");
  VRB.Smalloc(cname,fname, "receive_buf", receive_buf, blcklength*sizeof(IFloat));
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
