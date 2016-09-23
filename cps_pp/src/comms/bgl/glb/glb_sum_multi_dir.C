#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_sum_multi_dir routine.
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /space/cvs/cps/cps++/src/comms/bgl/glb/glb_sum_multi_dir.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
// glb_sum_dir
// sum over multiple IFloating point data with single precision
// Sum over all nodes along a direction
// (0,1,2,3,4) <-> (x,y,z,t,s)
//--------------------------------------------------------------
CPS_END_NAMESPACE
#include<comms/glb.h>
#include<comms/scu.h>
#include<util/gjp.h>
#include<comms/double64.h>
#include <comms/sysfunc_cps.h>
#include <util/lat_data.h>
#include <comms/glb_sum_internal.h>
CPS_START_NAMESPACE

#define MAX_NUM_WORDS 40  //set to larger than size of double precision 3x3 matrix
static Double64 *transmit_buf = NULL;
static Double64 *receive_buf = NULL;
static Double64 *gsum_buf = NULL;

//----------------------------------------------------------------------
/*!
  \param float_p The number to be summed.
  \param dir The direction in which to sum; one of {0, 1, 2, 3, 4},
  corresponding to {x, y, z, t, s}.
  \param len The number of floating point numbers in the vector.
  \post The vector pointed to by \a float_p is summed over all nodes along the
  \a dir direction, \e i.e. over each strip of nodes where the grid
  coordinates in all other directions are constant.
  That sum vector is written back to \a float_p, which is identical on all
  nodes in this strip.

  \ingroup comms
*/
//---------------------------------------------------------------------- 

void glb_sum_multi_dir(LatData &lat, const int dir){
    glb_sum_multi_dir(lat.Field(),dir,(const int)lat.Size());
}

void glb_sum_multi_dir(const Float * float_p, const int dir, const int len)
#if 1
{
//   printf("glb_sum_multi_dir(%p, %d, %d)\n",float_p,dir,len);
	int len2 = len ;
	Float *tmp = (Float *)float_p;
	while( len2 > MAX_BUF){
		glb_sum_internal(tmp,dir,MAX_BUF);
		len2 -=MAX_BUF;
		tmp += MAX_BUF;
	}
	if (len2>0)
		glb_sum_internal(tmp,dir,len2);
  //      printf("done\n");
}
#else
#endif


CPS_END_NAMESPACE
