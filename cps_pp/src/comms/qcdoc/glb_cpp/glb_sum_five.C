#include<config.h>
#include<qalloc.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_sum_five routine.

*/
//--------------------------------------------------------------------
//  CVS keywords
//  $Source: /space/cvs/cps/cps++/src/comms/qcdoc/glb_cpp/glb_sum_five.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
//  glb_sum_five
//
// Sum over all nodes of the "virtual" 5-dimensional volume.
// {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes(), GJP.Snodes}
// Relevant for spread-out DWF (GJP.s_nodes not 1) only.
//--------------------------------------------------------------


CPS_END_NAMESPACE
#include<comms/glb.h>
#include<comms/scu.h>
#include<util/gjp.h>
#include<comms/double64.h>
#include <comms/sysfunc_cps.h>
#include <comms/glb_sum_internal.h>
CPS_START_NAMESPACE


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

void glb_sum_five(Float * float_p)
{
	glb_sum_internal2(float_p,5);
}

CPS_END_NAMESPACE
