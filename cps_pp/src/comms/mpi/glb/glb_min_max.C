#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_min and glb_max routine.

  $Id: glb_min_max.C,v 1.4 2008-02-08 18:35:06 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008-02-08 18:35:06 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/mpi/glb/glb_min_max.C,v 1.4 2008-02-08 18:35:06 chulwoo Exp $
//  $Id: glb_min_max.C,v 1.4 2008-02-08 18:35:06 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: glb_min_max.C,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/mpi/glb/glb_min_max.C,v $
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

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))


const SCUDir dir[] = { SCU_XP, SCU_XM, SCU_YP, SCU_YM,
                       SCU_ZP, SCU_ZM, SCU_TP, SCU_TM };



static Float transmit_buf;
static Float receive_buf;
static Float gsum_buf;


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
  int NP[4] = {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()};

  gsum_buf = *float_p;

  for(int i = 0; i < 4; ++i) {

      transmit_buf = gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
	SCUDirArg send(&transmit_buf, dir[2*i], SCU_SEND, sizeof(Float));
	SCUDirArg rcv(&receive_buf, dir[2*i+1], SCU_REC, sizeof(Float));

	SCUTrans(&send);
	SCUTrans(&rcv);

	SCUTransComplete();

        gsum_buf = max(gsum_buf, receive_buf) ;
        transmit_buf = receive_buf;
      }
  }
  *float_p = gsum_buf;
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
  int NP[4] = {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()};

  gsum_buf = *float_p;

  for(int i = 0; i < 4; ++i) {

      transmit_buf = gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
	SCUDirArg send(&transmit_buf, dir[2*i], SCU_SEND, sizeof(Float));
	SCUDirArg rcv(&receive_buf, dir[2*i+1], SCU_REC, sizeof(Float));

	SCUTrans(&send);
	SCUTrans(&rcv);

	SCUTransComplete();

        gsum_buf = min(gsum_buf, receive_buf) ;
        transmit_buf = receive_buf;
      }
  }
  *float_p = gsum_buf;
}

CPS_END_NAMESPACE
