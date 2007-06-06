#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief Routines used by the AlgOverRelax class methods:

  $Id: alg_overrelax_kern.C,v 1.4 2007-06-06 16:06:22 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2007-06-06 16:06:22 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_overrelax/noarch/alg_overrelax_kern.C,v 1.4 2007-06-06 16:06:22 chulwoo Exp $
//  $Id: alg_overrelax_kern.C,v 1.4 2007-06-06 16:06:22 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_overrelax_kern.C,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_overrelax/noarch/alg_overrelax_kern.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// alg_overrelax.C
//
// AlgOverRelax is derived from Alg and is relevant to the 
// gauge overrelaxation algorithm. The type of gauge action is
// determined by the argument to the constructor.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdlib.h>	// exit()
#include <util/qcdio.h>
#include <math.h>
#include<alg/alg_overrelax.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/random.h>
#include<util/smalloc.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/algebra.h>
CPS_START_NAMESPACE





//static int	NHITS = 10;
//static Float	SMALL = 0.2;
//static Float	EPSILON[MATRIXSIZE];
static Float	MTEMP1[MATRIXSIZE];
static Float	u_sigma[MATRIXSIZE];

//*******************************************************
//     Here is a Metropolis update kernel...  		*
//*******************************************************
void AlgOverRelax::kernel(  Float *sigma, Float *U)
{

 
  // These variables are associated with the actual update. 
  Float  old_action;
  Float  new_action;
  Float  accept_probability;

  
  // remember - you cannot use multi-hit here, it breaks detailed balance
  //generate G0xGxG0
  m_multiply3( u_sigma, U, sigma );
  m_multiply2r( u_sigma, sigma );
  m_conjugate(u_sigma);
  //yes I know this is slow
  
  
  // Now begin the update process. 
  
  m_multiply3( MTEMP1, sigma, U );
  old_action = fGamma * (9. - .5 * (MTEMP1[0]+MTEMP1[4]+MTEMP1[8]) );
    
  m_multiply3( MTEMP1, u_sigma, U );
  
  new_action = fGamma * (9. - .5 * (MTEMP1[0]+MTEMP1[4]+MTEMP1[8]) );
  
  accept_probability = exp( old_action - new_action );
  
  
  if(fabs(LRG.Urand()) < accept_probability ) {
    m_equal( U, u_sigma );
  }

  
}             /* End HIT_CNT loop for multiple hits.           */




CPS_END_NAMESPACE
