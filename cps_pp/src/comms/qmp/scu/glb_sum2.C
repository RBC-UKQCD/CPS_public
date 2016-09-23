#include<config.h>
#ifdef USE_QMP
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_sum_five routine.

*/
//--------------------------------------------------------------------
//  CVS keywords
//  $Source: /space/cvs/cps/cps++/src/comms/qmp/scu/glb_sum2.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
//  glb_sum_internal2
//
//  An internal routine for glb_sum(), glb_sum_five()
//--------------------------------------------------------------


CPS_END_NAMESPACE
#include<comms/glb.h>
#include<comms/scu.h>
#include<util/gjp.h>
#include<util/checksum.h>
#include<util/data_shift.h>
#include<comms/double64.h>
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE

union DoubleBytes {
    Float dblval;
    char byte[8];
};

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
void glb_sum_internal2(Float * float_p,int ndir)
{
  static int NP[5] ={0,0,0,0,0};
  static int coor[5] ={0,0,0,0,0};
//  static SCUDirArgIR *Send[5];
//  static SCUDirArgIR *Recv[5];

#ifndef  UNIFORM_SEED_NO_COMMS
  static QMP_mem_t *trans_mem, *gsum_mem;
  static QMP_msgmem_t msgmem2[10];
  static QMP_msghandle_t msghandle2[10];
  static QMP_msghandle_t sndrcv_msghandle2[5];
#endif

  if (!initted){
      NP[0] = GJP.Xnodes();
      int max=NP[0];
      NP[1] = GJP.Ynodes();
      NP[2] = GJP.Znodes();
      NP[3] = GJP.Tnodes();
      NP[4] = GJP.Snodes();
      for (int i = 1;i<5;i++)
	if (max <NP[i]) max = NP[i];
#ifndef  UNIFORM_SEED_NO_COMMS
      trans_mem = QMP_allocate_memory(sizeof(Double64)*2);
      transmit_buf = (Double64 *)QMP_get_memory_pointer(trans_mem);
      gsum_mem = QMP_allocate_memory(sizeof(Double64)*max);
      gsum_buf = (Double64 *)QMP_get_memory_pointer(gsum_mem);
      if(!UniqueID())
        printf("transmit_buf=%p gsum_buf=%p\n",transmit_buf,gsum_buf);
      receive_buf = transmit_buf+1;

      for(int i = 0;i<5;i++)
      if (NP[i]>1){
      msgmem2[2*i] = QMP_declare_msgmem(transmit_buf, sizeof(Double64));
      msghandle2[2*i] = QMP_declare_send_relative(msgmem2[2*i], i, -1, 0);
      msgmem2[2*i+1] = QMP_declare_msgmem(receive_buf, sizeof(Double64));
      msghandle2[2*i+1] = QMP_declare_receive_relative(msgmem2[2*i+1], i, +1, 0);
//      sndrcv_msghandle2[i] = QMP_declare_multiple(&msghandle2[2*i], 2);
      }
#endif
  }
  initted = 1;


  // Sum over the "virtual" 5-dimensional mesh
  //------------------------------------------------------------
//  gsum_buf[0] = (Double64)*float_p;

  Double64 tmp_sum = (Double64)*float_p;

  // Save checksum of local floating points
  //---------------------------------------------------
  unsigned long sum_sum = local_checksum(float_p,1);
 CSM.AccumulateCsum(CSUM_GLB_LOC,sum_sum);



#ifndef  UNIFORM_SEED_NO_COMMS
  for(int i = 0; i < ndir; ++i) 
  if (NP[i] >1) {
      int coor = (GJP.NodeCoor(i)-GDS.Origin(i)+NP[i])%NP[i];
//fprintf(stderr,"coor[%d]=%d\n",i,coor);
      *transmit_buf = gsum_buf[coor]= tmp_sum;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
	coor = (coor+1)%NP[i];
//	QMP_start(sndrcv_msghandle2[i]);
    QMP_start(msghandle2[2*i]);
    QMP_start(msghandle2[2*i+1]);
//	QMP_status_t status = QMP_wait(sndrcv_msghandle2[i]);
    QMP_status_t status = QMP_wait (msghandle2[2*i]);
	if (status != QMP_SUCCESS)
	  QMP_error("Error in glb_sum_internal2: %s\n", QMP_error_string(status));
    status = QMP_wait (msghandle2[2*i+1]);
	if (status != QMP_SUCCESS)
	  QMP_error("Error in glb_sum_internal2: %s\n", QMP_error_string(status));

        gsum_buf[coor] = *receive_buf;
        *transmit_buf = *receive_buf;
      }
      tmp_sum = gsum_buf[0];
      for (int itmp = 1; itmp < NP[i]; itmp++) {
	    tmp_sum += gsum_buf[itmp];
      }
  }
  *float_p = (Float)tmp_sum;
#endif


  // accumulate final global sum checksum
  //------------------------------------
  sum_sum = local_checksum(float_p,1);
  CSM.AccumulateCsum(CSUM_GLB_SUM,sum_sum);  
}

static int initted_u=0;

static unsigned long long *transmit_buf_u = NULL;
static unsigned long long *receive_buf_u = NULL;
static unsigned long long *gsum_buf_u = NULL;

void glb_sum_internal2(unsigned int *uint_p, int ndir, int sum_flag) {
  static int NP[5] = {0,0,0,0,0};
  static int coor[5] = {0,0,0,0,0};
#ifndef  UNIFORM_SEED_NO_COMMS
  static QMP_mem_t *trans_mem, *gsum_mem;
  static QMP_msgmem_t msgmem[10];
  static QMP_msghandle_t msghandle[10];
  static QMP_msghandle_t sndrcv_msghandle[5];
#endif

  if (!initted_u) {
	NP[0] = GJP.Xnodes();
	int max = NP[0];
	NP[1] = GJP.Ynodes();
	NP[2] = GJP.Znodes();
	NP[3] = GJP.Tnodes();
	NP[4] = GJP.Snodes();
	for (int i = 1;i<5;i++)
	  if (max<NP[i]) max = NP[i];
#ifndef  UNIFORM_SEED_NO_COMMS
      trans_mem = QMP_allocate_memory(sizeof(unsigned long long)*2);
      transmit_buf_u = (unsigned long long *)QMP_get_memory_pointer(trans_mem);
      gsum_mem = QMP_allocate_memory(sizeof(unsigned long long)*max);
      gsum_buf_u = (unsigned long long *)QMP_get_memory_pointer(gsum_mem);
	receive_buf_u = transmit_buf_u+1;

	for(int i=0; i<5; i++)
      if (NP[i]>1) {

      msgmem[2*i] = QMP_declare_msgmem(transmit_buf_u, sizeof(unsigned long long));
      msghandle[2*i] = QMP_declare_send_relative(msgmem[2*i], i, -1, 0);
      msgmem[2*i+1] = QMP_declare_msgmem(receive_buf_u, sizeof(unsigned long long));
      msghandle[2*i+1] = QMP_declare_receive_relative(msgmem[2*i+1], i, +1, 0);
      sndrcv_msghandle[i] = QMP_declare_multiple(&msghandle[2*i], 2);
      }
#endif
  }
  initted_u = 1;
  
  // Sum over the "virtual" 5-dimensional mesh
  //------------------------------------------------------------
  //  gsum_buf_u[0] = (unsigned long long)*float_p;
  
  unsigned int tmp_sum = *uint_p;
  
#ifndef  UNIFORM_SEED_NO_COMMS
  for(int i=0; i<ndir; ++i) 
	if (NP[i] > 1) {
      int coor = GJP.NodeCoor(i);
	  //printf("coor[%d]=%d\n",i,coor);
      *transmit_buf_u = gsum_buf_u[coor] = (unsigned long long)tmp_sum;
	  
      for (int itmp = 1; itmp < NP[i]; itmp++) {
		coor = (coor+1)%NP[i];
		QMP_start(sndrcv_msghandle[i]);
		QMP_status_t status = QMP_wait(sndrcv_msghandle[i]);
		if (status != QMP_SUCCESS)
		  QMP_error("Error in glb_sum_internal2: %s\n", QMP_error_string(status));
		
        gsum_buf_u[coor] = *receive_buf_u;
        *transmit_buf_u = *receive_buf_u;
      }
      tmp_sum = (unsigned int)gsum_buf_u[0];
      for (int itmp = 1; itmp < NP[i]; itmp++) {
	if (sum_flag) tmp_sum += (unsigned int)gsum_buf_u[itmp];
	else  tmp_sum ^= (unsigned int)gsum_buf_u[itmp];
      }
	}

  *uint_p = tmp_sum;
#endif
}


CPS_END_NAMESPACE
#endif
