#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_sum_dir routine.

*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /space/cvs/cps/cps++/src/comms/bgl/glb/glb_sum_dir.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
// glb_sum_dir
//
// Sum over all nodes along a direction
// (0,1,2,3,4) <-> (x,y,z,t,s)
//--------------------------------------------------------------
CPS_END_NAMESPACE
#include<comms/glb.h>
#include<comms/scu.h>
#include<util/gjp.h>
#include<comms/double64.h>
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE

void glb_sum_internal (Float * float_p, int dir,int len);


static Double64 *transmit_buf = NULL;
static Double64 *receive_buf = NULL;
static Double64 *gsum_buf = NULL;

//----------------------------------------------------------------------
/*!
  \param float_p The number to be summed.
  \param dir The direction in which to sum; one of {0, 1, 2, 3, 4},
  corresponding to {x, y, z, t, s}.
  \post The number pointed to by \a float_p is summed over all nodes along the
  \a dir direction, \e i.e. over each strip of nodes where the grid
  coordinates in all other directions are constant.
  That sum is written back to \a float_p, which is identical on all nodes in
  this strip.

  \ingroup comms
*/
//---------------------------------------------------------------------- 

void glb_sum_dir(Float * float_p, int dir)
{
  glb_sum_internal(float_p,dir,1);
}


CPS_END_NAMESPACE
