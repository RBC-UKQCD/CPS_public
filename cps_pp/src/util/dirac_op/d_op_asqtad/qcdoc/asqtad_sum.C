#include "asq_data_types.h"
#include "asqtad_int.h"

typedef double Double64;

static Double64 *transmit_buf = NULL;
static Double64 *receive_buf = NULL;
static Double64 *gsum_buf = NULL;

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

static int initted=0;
static SCUDirArgIR *Send[5];
static SCUDirArgIR *Recv[5];
void AsqD::Sum(Float * float_p)
{
//  static int NP[5] ={0,0,0,0,0};
//  static int coor[5] ={0,0,0,0,0};

  int ndir = 4;
  if (!initted){
      int max=NP[0];
      for (int i = 1;i<ndir;i++)
	if (max <NP[i]) max = NP[i];
      transmit_buf = (Double64 *)qalloc(QFAST|QNONCACHE,sizeof(Double64)*2);
      if(!transmit_buf) PointerErr(cname, "Sum", "transmit_buf");
      receive_buf = transmit_buf+1;
      gsum_buf = (Double64 *)qalloc(QFAST,sizeof(Double64)*max);
      if(!gsum_buf) PointerErr(cname, "Sum", "gsum_buf");
      for(int i = 0;i<ndir;i++)
      if (NP[i]>1){
      Send[i] = new SCUDirArgIR(transmit_buf, scudir[i+4], SCU_SEND, sizeof(Double64));
      Recv[i] = new SCUDirArgIR(receive_buf, scudir[i], SCU_REC, sizeof(Double64));
      }
  }
  initted = 1;

  // Sum over the "virtual" 5-dimensional mesh
  //------------------------------------------------------------
//  gsum_buf[0] = (Double64)*float_p;

  Double64 tmp_sum = (Double64)*float_p;
  
  for(int i = 0; i < ndir; ++i) 
  if (NP[i] >1) {
      int pos = coor[i];
//printf("coor[%d]=%d\n",i,coor[i]);
      *transmit_buf = gsum_buf[pos]= tmp_sum;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
	pos = (pos+1)%NP[i];
	Send[i]->StartTrans(); Recv[i]->StartTrans();
	Send[i]->TransComplete(); Recv[i]->TransComplete();

        gsum_buf[pos] = *receive_buf;
        *transmit_buf = *receive_buf;
      }
      tmp_sum = gsum_buf[0];
      for (int itmp = 1; itmp < NP[i]; itmp++) {
	    tmp_sum += gsum_buf[itmp];
      }
  }
  *float_p = (Float)tmp_sum;

}

