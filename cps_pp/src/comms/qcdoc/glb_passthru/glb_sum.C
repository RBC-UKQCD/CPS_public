#include<config.h>
#include<stdio.h>
#include<qcdoc_align.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief Definition of glb_sum routine.

  $Id: glb_sum.C,v 1.4 2004-04-27 03:51:17 cwj Exp $ 
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: cwj $
//  $Date: 2004-04-27 03:51:17 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdoc/glb_passthru/glb_sum.C,v 1.4 2004-04-27 03:51:17 cwj Exp $
//  $Id: glb_sum.C,v 1.4 2004-04-27 03:51:17 cwj Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: glb_sum.C,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdoc/glb_passthru/glb_sum.C,v $
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



static Gsum64 gsum;
static SCUAxis gsum_axis[]={SCU_X,SCU_Y,SCU_Z,SCU_T};
static int initted=0;
static int counter = 0;


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
   if (!initted){ gsum.Init(gsum_axis,4); initted=1;};
//  printf("glb_sum %d before = %e ",counter,(double)*float_p);
  *float_p = (Float) gsum.Sum((double)*float_p);
//  printf("after = %e\n",(double)*float_p);
//  counter ++;
}
CPS_END_NAMESPACE
