#include<config.h>
#include<util/qcdio.h>
//#include<qalloc.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief Definition of glb_sum routine.
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qmp/glb_cpp/glb_sum.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
// glb_sum
//
// Sum over all nodes 
// {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()}
//--------------------------------------------------------------
CPS_END_NAMESPACE
#include<comms/glb.h>
#include<comms/scu.h>
#include<util/gjp.h>
#include<comms/double64.h>
#include <comms/sysfunc.h>
#include <comms/glb_sum_internal.h>
CPS_START_NAMESPACE
#ifndef USE_QMP
#define USE_QMP
#endif


static Double64 *transmit_buf = NULL;
static Double64 *receive_buf = NULL;
static Double64 *gsum_buf = NULL;

static int output = 0;


//----------------------------------------------------------------------
/*!
  \param float_p The number to be summed.
  \post The number pointed to by \a float_p is summed over all nodes
  and that sum is written back to \a float_p, which is identical on all nodes.

  \ingroup comms
*/
//---------------------------------------------------------------------- 
void glb_sum(Float * float_p)
{
#ifdef USE_QMP
  QMP_sum_double(float_p);
#else
  glb_sum_internal2(float_p,4);
#endif
}

void glb_sum_gimp(Float * float_p)
{
  int NP[4] = {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()};
  static int counter = 0;

  #ifdef USE_QMP
  QMP_sum_double(float_p);
  #else
  if (transmit_buf == NULL) 
      transmit_buf = (Double64 *)qalloc(QFAST|QNONCACHE,sizeof(Double64));
  if (receive_buf == NULL) 
      receive_buf = (Double64 *)qalloc(QFAST|QNONCACHE,sizeof(Double64));
  if (gsum_buf == NULL) 
      gsum_buf = (Double64 *)qalloc(QFAST|QNONCACHE,sizeof(Double64));
  if (output) printf("glb_sum cpp %d before = %e ", counter, (double)*float_p);
  *gsum_buf = (Double64)*float_p;

  #ifdef USE_QMP
  QMP_msgmem_t msgmem[2];
  QMP_msghandle_t msghandle[2];
  QMP_msghandle_t sndrcv;
  #endif

  for(int i = 0; i < 4; ++i) {

      *transmit_buf = *gsum_buf;
      #ifdef USE_QMP
      msgmem[0] = QMP_declare_msgmem((void *)transmit_buf, sizeof(Double64));
      msghandle[0] = QMP_declare_send_relative(msgmem[0], i, 1, 0);
      msgmem[1] = QMP_declare_msgmem((void *)receive_buf, sizeof(Double64));
      msghandle[1] = QMP_declare_receive_relative(msgmem[1], i, -1, 0);
      sndrcv = QMP_declare_multiple(msghandle, 2);
      #endif

      for (int itmp = 1; itmp < NP[i]; itmp++) {
	#ifndef USE_QMP
	SCUDirArg send(transmit_buf, gjp_scu_dir[2*i], SCU_SEND, sizeof(Double64) );
	SCUDirArg rcv(receive_buf, gjp_scu_dir[2*i+1], SCU_REC, sizeof(Double64) );

	send.StartTrans();
	rcv.StartTrans();
	send.TransComplete();
	rcv.TransComplete();
	#else
	QMP_start(sndrcv);
	QMP_status_t status = QMP_wait(sndrcv);
	if (status != QMP_SUCCESS)
	  QMP_error("Communication error in glb_sum_gimp:%s\n", QMP_error_string(status));
	#endif

        *gsum_buf += *receive_buf;
        *transmit_buf = *receive_buf;
      }
      #ifdef USE_QMP
      QMP_free_msghandle(sndrcv);
      QMP_free_msgmem(msgmem[0]);
      QMP_free_msgmem(msgmem[1]);
      #endif
  }
  *float_p = (Float)*gsum_buf;
  if (output)   printf("after = %e\n", (double)*float_p);
  #endif
  counter++;
}

CPS_END_NAMESPACE
