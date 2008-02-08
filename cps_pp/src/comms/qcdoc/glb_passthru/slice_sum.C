#include<config.h>
#include<qalloc.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief Definition of slice_sum routine
 */
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdoc/glb_passthru/slice_sum.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

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

  IFloat *transmit_buf = (IFloat *)qalloc(QFAST|QNONCACHE,blcklength*sizeof(IFloat));

  // added by manke
  if(transmit_buf == 0)
    ERR.Pointer(cname,fname, "transmit_buf");
  VRB.Smalloc(cname,fname, "transmit_buf", transmit_buf, blcklength*sizeof(IFloat));
  // end added

  IFloat *receive_buf  = (IFloat *)qalloc(QFAST|QNONCACHE,blcklength*sizeof(IFloat));

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
  int i,j;

  for(i = 0; i < 4; ++i) 
      if( (i != dir) && (NP[i]>1)){

      //--------------------------------------------------------------
      // address of buffer of to be sent (data on this node) 
      //--------------------------------------------------------------
      for(j = 0;j<blcklength;j++) transmit_buf_p[j] = float_p[j];

      SCUDirArg send(transmit_buf_p, pos_dir[i], SCU_SEND, blcklength*sizeof(IFloat));
      SCUDirArg rcv(receive_buf_p, neg_dir[i], SCU_REC, blcklength*sizeof(IFloat));

      //--------------------------------------------------------------
      // tranmit & receive NP[i] - 1 times in snd_dir[i] direction
      //--------------------------------------------------------------
      for ( itmp = 1; itmp < NP[i]; itmp++) {

	 //-----------------------------------------------------------
         // do SCU transfers
	 //-----------------------------------------------------------
    send.StartTrans();
    rcv.StartTrans();
    send.TransComplete();
    rcv.TransComplete();

	 //-----------------------------------------------------------
         // accumulate received data	
	 //-----------------------------------------------------------
 //        printf("float_p[%d]=%e\n",blcklength-1,float_p[blcklength-1]);
   	 for (count = 0; count < blcklength; ++count) 
           float_p[count] += receive_buf_p[count];
//         printf("float_p[%d](after)=%e\n",blcklength-1,float_p[blcklength-1]);

	 //-----------------------------------------------------------
         // the received data will be sent out     
	 // the free buffer will be used to receive      
	 //-----------------------------------------------------------
         transmit_buf_p = receive_buf_p;
         receive_buf_p = free_buffer_p;
	 send.Addr(transmit_buf_p );
	 rcv.Addr(receive_buf_p );

	 //-----------------------------------------------------------
	 // transmit_buf WILL be free buffer NEXT round of transmit
	 //-----------------------------------------------------------
	 free_buffer_p = transmit_buf_p;
      }
  }
  qfree(transmit_buf);
  qfree(receive_buf);
}


CPS_END_NAMESPACE
