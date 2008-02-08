#include<config.h>
#include<qalloc.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_min and glb_max routine.
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdoc/glb_passthru/glb_min_max.C,v $
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
#include <qcdocos/gsum64.h>
CPS_START_NAMESPACE

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))


const SCUDir dir[] = { SCU_XP, SCU_XM, SCU_YP, SCU_YM,
                       SCU_ZP, SCU_ZM, SCU_TP, SCU_TM };



//static Float *transmit_buf = NULL;
//static Float *receive_buf = NULL;
static Gsum64 gminmax;
static SCUAxis gsum_axis[]={SCU_X,SCU_Y,SCU_Z,SCU_T};
static int initted=0;
//static int counter = 0;


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
  if (!initted){ gminmax.Init(gsum_axis,4); initted=1;};
  *float_p = (Float) gminmax.Max((double)*float_p);
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
  if (!initted){ gminmax.Init(gsum_axis,4); initted=1;};
  *float_p = (Float) gminmax.Min((double)*float_p);
}

CPS_END_NAMESPACE
