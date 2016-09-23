#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_sum_dir routine.
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /space/cvs/cps/cps++/src/comms/qmp/glb_passthru/glb_sum_dir.C,v $
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
#include <qcdocos/gsum64.h>
CPS_START_NAMESPACE


static Double64 transmit_buf;
static Double64 receive_buf;
static Double64 gsum_buf;
static Gsum64 gsum;
static SCUAxis gsum_axis[]={SCU_X,SCU_Y,SCU_Z,SCU_T,SCU_S,SCU_W};
static int initted=0;

//----------------------------------------------------------------------
/*!
  \param float_p The number to be summed.
  \param dir The direction in which to sum; one of {0, 1, 2, 3, 4,5},
  corresponding to {x, y, z, t, s,w}.
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
	if(!initted){ gsum.Init(&gsum_axis[dir],1); initted=1;}
	*float_p = (Float) gsum.Sum( (double) *float_p);
}


CPS_END_NAMESPACE
