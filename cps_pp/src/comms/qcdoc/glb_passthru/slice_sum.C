#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief Definition of slice_sum routine

  $Id: slice_sum.C,v 1.3 2004-06-04 21:14:01 chulwoo Exp $
 */
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:01 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdoc/glb_passthru/slice_sum.C,v 1.3 2004-06-04 21:14:01 chulwoo Exp $
//  $Id: slice_sum.C,v 1.3 2004-06-04 21:14:01 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdoc/glb_passthru/slice_sum.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//====================================================================
//*  SUI 3/27/97
//*  slice_sum.C
//*  last modified 11/7/97
//====================================================================


CPS_END_NAMESPACE
#include<comms/glb.h>
#include<util/smalloc.h>
#include<util/gjp.h>
#include <comms/sysfunc.h>
#include <qcdocos/gsum64.h>
CPS_START_NAMESPACE

static SCUAxis gsum_axis[]= {SCU_X,SCU_Y,SCU_Z,SCU_T};
static Gsum64 gsum;


//-------------------------------------------------------------------
//* sum over a slice(hyperplane) which is orthogonal to the direction
//* "dir". There are "blcklength" summations to be done:
//* float_p[i] = sum_over_nodes_of_this_slice(float_p[i])
//-------------------------------------------------------------------

//-------------------------------------------------------------------
/*!
  The vector pointed to by \a float_p is summed over each 3-dimensional
  hyperplane of nodes which is perpendiculat to the \a dir direction

  \param float_p The number to be summed.
  \param blcklength The number of floating point numbers in the vector.
  \param dir The normal direction defining the hyperplane; one of {0, 1, 2, 3} 
  corresponding to {x, y, z, t}.
  \post The vector sum is written back to \a float_p, which is identical on
  all nodes in this hyperplane.

  \ingroup comms
*/
//-------------------------------------------------------------------
void slice_sum(Float * float_p, int blcklength, int dir)
{
  char *cname = "slice_sum";
  char *fname = "slice_sum(*float_p, int, int)";

  const int MAX=1023;
  SCUAxis axes[3];
  int i,j =0 ;
  for(i=0;i<4;i++)
    if(i!=dir){
      axes[j]=gsum_axis[i]; j++;
    } 
  gsum.Init(axes,3);

  if (blcklength > MAX)
    ERR.General(cname, fname, "blcklength (%d) too big > MAX (%d) \n",blcklength, MAX);
  for(i=0;i<blcklength;i++){
    *(float_p+i) = (Float) gsum.Sum( (double) *(float_p+i));
  }
		    
}


CPS_END_NAMESPACE
