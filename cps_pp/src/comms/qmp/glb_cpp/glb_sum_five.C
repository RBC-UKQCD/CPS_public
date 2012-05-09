#include<config.h>
#ifdef USE_QMP
//#include<qalloc.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_sum_five routine.

*/
//--------------------------------------------------------------------
//  CVS keywords
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qmp/glb_cpp/glb_sum_five.C,v $
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
#include<util/data_shift.h>
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
  static int initted=0;
#ifdef UNIFORM_SEED_TESTING
  double temp=0.;
  unsigned long *temp_i = (unsigned long *)&temp;
  if (initted==0 && !UniqueID()) fprintf(stderr,"USING UNIFORM_SEED glb_sum_five() %x %x\n",*temp_i, *(temp_i+1));
  if(!CoorX() &&!CoorY()&&!CoorZ()&&!CoorT())
    temp=*float_p;
  glb_sum_internal2(&temp,4);
//  QMP_sum_double(&temp);
  if (*float_p != temp){
      fprintf(stderr, "glb_sum_five %d  : Node %d : Oops I did it again: me (%d,%d,%d,%d,%d) %0.14e != Node 0: %0.14e\n",
        initted, UniqueID(), CoorX(), CoorY(), CoorZ(), CoorT(), CoorS(),
               *float_p,  temp );
//    exit(-30);
	    float *temp = NULL; *temp=1.;
  }
  glb_sum_internal2(float_p,5);
//  QMP_sum_double(float_p);
#else
  Float shifted = *float_p;
  QMP_sum_double(float_p);
#if 0
  if (initted==0 && !UniqueID()) fprintf(stderr,"shifting glb_sum_five (1,1,1,0)\n");
  GDS.Set(1,1,1,0);
  GDS.Shift(&shifted,sizeof(Float));
  QMP_sum_double(&shifted);
  if(*float_p != shifted){
      fprintf(stderr, "glb_sum_five %d: Node %d : shifted sum doesn't agree : me (%d,%d,%d,%d,%d) %0.16e != shfted: %0.16e\n",
        initted,   UniqueID(), CoorX(), CoorY(), CoorZ(), CoorT(), CoorS(),
               *float_p,  shifted );
    exit(-30);
  }
  if (initted==0 && !UniqueID()) fprintf(stderr,"shifting glb_sum_five (1,1,1,0) done \n");
#endif
#endif
  initted++;
}

CPS_END_NAMESPACE
#endif
