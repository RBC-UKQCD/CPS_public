#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief   Methods for the Random Number Generator classes.

  $Id: random_asm.C,v 1.5 2004-07-02 14:13:43 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-07-02 14:13:43 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/random/noarch/random_asm.C,v 1.5 2004-07-02 14:13:43 chulwoo Exp $
//  $Id: random_asm.C,v 1.5 2004-07-02 14:13:43 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: random_asm.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/random/noarch/random_asm.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

//---------------------------------------------------------------
//  This is the routine ran3 from Numerical Recipes in C 
//---------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdio.h>
#include <util/random.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <math.h>
CPS_START_NAMESPACE


/*!
  \pre ::Reset must have already been called.
  \return A random number from a uniform distribution over (0,1)
 */
IFloat RandomGenerator::Rand(void)
{
    //---------------------------------------------------------
    // Start generating random numbers here
    //---------------------------------------------------------
    if ( ++inext == state_size ) inext = 0 ;
    if ( ++inextp == state_size ) inextp = 0 ;
    int mj = ma[inext] - ma[inextp] ;
    if ( mj < 0 ) mj += MBIG ;
    ma[inext] = mj ;
//    printf("mj=%d\n",mj);
    return mj*FAC ;
}




/*!
  \pre ::Reset must have already been called.
  \return A random number from a gaussian distribution with zero mean.
*/
IFloat GaussianRandomGenerator::Rand()
{
    if(iset == 0) {	// We don't have an extra deviate handy
	IFloat v1, v2, rsq;
        do {
	    v1 = 2.0 * RandomGenerator::Rand() - 1.0;
	    v2 = 2.0 * RandomGenerator::Rand() - 1.0;
//	    printf("v1=%e v2=%e\n",v1,v2);
	    rsq = v1*v1 + v2*v2;
	} while(rsq >= 1.0 || rsq == 0);
	    // pick 2 uniform numbers in the square extending from
	    // -1 to 1 in each direction, see if they are in the
	    // unit circle, and try again if they are not.

	IFloat fac = sqrt(-2.0 * sigma2 * log(rsq)/rsq);

	gset = v1 * fac;   iset = 1;
	return v2 * fac;

    } else {
    	iset = 0;
	return gset;
    }
}

/*!
  \pre ::Reset must have already been called.
  \return A random number from a uniform distribution.
*/
IFloat UniformRandomGenerator::Rand()
{
    return A + ((B - A) * RandomGenerator::Rand());
}



CPS_END_NAMESPACE
