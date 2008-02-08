#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_min and glb_max routine.
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/bgl/glb/glb_min_max.C,v $
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
#include<util/qcdio.h>
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))




//----------------------------------------------------------------------
/*!
  The number pointed to by \a float_p on all nodes is compared.
  \param float_p The number to be compared.
  \post The maximum is found and is written back to \a float_p, which is
  then identical on all nodes.

  \ingroup comms
*/
//---------------------------------------------------------------------- 

static int initted=0;
static FILE *fp;

inline void glb_minmax(Float * float_p, int n_dir, int minmax)
{
  int NP[4] = {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()};

//  if (!initted)
//    fp = Fopen(ADD_ID,"glb_minmax","a");

  initted=1;

   Float transmit_buf;
   Float receive_buf;
   Float gsum_buf;
  gsum_buf = *float_p;

  for(int i = 0; i < n_dir; ++i) {

      transmit_buf = gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
        getPlusData(&receive_buf,&transmit_buf,1,i);

        if (minmax)
        gsum_buf = max(gsum_buf, receive_buf) ;
        else
        gsum_buf = min(gsum_buf, receive_buf) ;
        transmit_buf = receive_buf;
      }
  }
  *float_p = gsum_buf;
//  Fprintf(fp," %d %0.20e\n",minmax,*float_p);
//  printf("glb_minmax %d: %d %0.20e\n",UniqueID(),minmax,*float_p);
  
}

void glb_min(Float * float_p){
  glb_minmax(float_p,4,0);
}

void glb_max(Float * float_p){
  glb_minmax(float_p,4,1);
}

#if 0
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
	SCUDirArg send(transmit_buf, dir[2*i], SCU_SEND, sizeof(Float));
	SCUDirArg rcv(receive_buf, dir[2*i+1], SCU_REC, sizeof(Float));

	send.StartTrans();
	rcv.StartTrans();
	send.TransComplete();
	rcv.TransComplete();

        *gsum_buf = min(*gsum_buf, *receive_buf) ;
        *transmit_buf = *receive_buf;
      }
  }
  *float_p = *gsum_buf;
}
#endif

CPS_END_NAMESPACE
