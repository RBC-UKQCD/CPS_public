#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief Definition of glb_sum routine.
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qmp/glb_passthru/glb_sum_five.C,v $
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
#include <qcdocos/gsum64.h>
CPS_START_NAMESPACE



static Double64 transmit_buf;
static Double64 receive_buf;
static Double64 gsum_buf;
static Gsum64 gsum;
static SCUAxis gsum_axis[]={SCU_X,SCU_Y,SCU_Z,SCU_T,SCU_S};
static int initted=0;


//----------------------------------------------------------------------
/*!
  \param float_p The number to be summed.
  \post The number pointed to by \a float_p is summed over all nodes
  and that sum is written back to \a float_p, which is identical on all nodes.

  \ingroup comms
*/
//---------------------------------------------------------------------- 
void glb_sum_five(Float * float_p)
{
//  static int called=0;
//  printf("glb_sum_five: before =%e\n",*float_p);
  if (!initted){ gsum.Init(gsum_axis,5); initted=1;}
  *float_p = (Float) gsum.Sum((double)*float_p);
//  printf("glb_sum_five: after =%e\n",*float_p);
//  called++;
//  if(called>=10) exit(-43);
}
CPS_END_NAMESPACE
