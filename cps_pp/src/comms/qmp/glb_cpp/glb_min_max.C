#include<config.h>
#ifdef USE_QMP
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_min and glb_max routine.
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /space/cvs/cps/cps++/src/comms/qmp/glb_cpp/glb_min_max.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/*
 *  glb_min_max.C (C++ version)
 */


CPS_END_NAMESPACE
#include<comms/glb.h>
#include<comms/scu.h>
#include<util/gjp.h>
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE
#ifndef USE_QMP
#define USE_QMP
#endif

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))


// const SCUDir dir[] = { SCU_XP, SCU_XM, SCU_YP, SCU_YM, SCU_ZP, SCU_ZM, SCU_TP, SCU_TM };



static Float *transmit_buf = NULL;
static Float *receive_buf = NULL;
static Float *gsum_buf = NULL;


//----------------------------------------------------------------------
/*!
  The number pointed to by \a float_p on all nodes is compared.
  \param float_p The number to be compared.
  \post The maximum is found and is written back to \a float_p, which is
  then identical on all nodes.

  \ingroup comms
*/
//---------------------------------------------------------------------- 

void glb_max(Float * float_p)
{
#if 1
#ifndef UNIFORM_SEED_TESTING
  QMP_max_double(float_p);
#endif
#else
  int NP[4] = {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()};

  if (transmit_buf == NULL) 
      transmit_buf = (Float *)qalloc(QFAST|QNONCACHE,sizeof(Float));
  if (receive_buf == NULL) 
      receive_buf = (Float *)qalloc(QFAST|QNONCACHE,sizeof(Float));
  if (gsum_buf == NULL) 
      gsum_buf = (Float *)qalloc(QFAST|QNONCACHE,sizeof(Float));
  *gsum_buf = *float_p;

  for(int i = 0; i < 4; ++i) {

      *transmit_buf = *gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {

	QMP_msgmem_t msgmem[2];
	QMP_msghandle_t msghandle[2];
	QMP_msghandle_t sndrcv;
	msgmem[0] = QMP_declare_msgmem((void *)transmit_buf, sizeof(Float));
	msghandle[0] = QMP_declare_send_relative(msgmem[0], i, 1, 0);
       	msgmem[1] = QMP_declare_msgmem((void *)receive_buf, sizeof(Float));
	msghandle[1] = QMP_declare_receive_relative(msgmem[1], i, -1, 0);
	sndrcv = QMP_declare_multiple(msghandle, 2);

	QMP_start(sndrcv);
	QMP_status_t status = QMP_wait(sndrcv);
	if (status != QMP_SUCCESS)
	  QMP_error("Communication error in glb_min:%s\n", QMP_error_string(status));
	QMP_free_msghandle(sndrcv);
	QMP_free_msgmem(msgmem[0]);
	QMP_free_msgmem(msgmem[1]);

        *gsum_buf = max(*gsum_buf, *receive_buf) ;
        *transmit_buf = *receive_buf;
      }
  }
  *float_p = *gsum_buf;
#endif
}


//----------------------------------------------------------------------
/*!
  The number pointed to by \a float_p on all nodes is compared.
  \param float_p The number to be compared.
  \post The minimum is found and is written back to \a float_p, which is
  then identical on all nodes.

  \ingroup comms
*/
//---------------------------------------------------------------------- 
void glb_min(Float * float_p)
{
#if 1
#ifndef UNIFORM_SEED_TESTING
  QMP_min_double(float_p);
#endif
#else
  int NP[4] = {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()};
  if (transmit_buf == NULL) 
      transmit_buf = (Float *)qalloc(QFAST|QNONCACHE,sizeof(Float));
  if (receive_buf == NULL) 
      receive_buf = (Float *)qalloc(QFAST|QNONCACHE,sizeof(Float));
  if (gsum_buf == NULL) 
      gsum_buf = (Float *)qalloc(QFAST|QNONCACHE,sizeof(Float));

  *gsum_buf = *float_p;

  for(int i = 0; i < 4; ++i) {

      *transmit_buf = *gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
	QMP_msgmem_t msgmem[2];
	QMP_msghandle_t msghandle[2];
	QMP_msghandle_t sndrcv;
	msgmem[0] = QMP_declare_msgmem((void *)transmit_buf, sizeof(Float));
	msghandle[0] = QMP_declare_send_relative(msgmem[0], i, 1, 0);
       	msgmem[1] = QMP_declare_msgmem((void *)receive_buf, sizeof(Float));
	msghandle[1] = QMP_declare_receive_relative(msgmem[1], i, -1, 0);
	sndrcv = QMP_declare_multiple(msghandle, 2);

	QMP_start(sndrcv);
	QMP_status_t status = QMP_wait(sndrcv);
	if (status != QMP_SUCCESS)
	  QMP_error("Communication error in glb_max:%s\n", QMP_error_string(status));
	QMP_free_msghandle(sndrcv);
	QMP_free_msgmem(msgmem[0]);
	QMP_free_msgmem(msgmem[1]);


        *gsum_buf = min(*gsum_buf, *receive_buf) ;
        *transmit_buf = *receive_buf;
      }
  }
  *float_p = *gsum_buf;
#endif
}

CPS_END_NAMESPACE
#endif
