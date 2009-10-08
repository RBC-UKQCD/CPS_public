#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief   Methods for the Random Number Generator classes.

  $Id: random_asm.C,v 1.17 2009-10-08 15:00:24 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2009-10-08 15:00:24 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/random/noarch/random_asm.C,v 1.17 2009-10-08 15:00:24 chulwoo Exp $
//  $Id: random_asm.C,v 1.17 2009-10-08 15:00:24 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: random_asm.C,v $
//  $Revision: 1.17 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/random/noarch/random_asm.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

//---------------------------------------------------------------
//  This is the routine ran3 from Numerical Recipes in C 
//---------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/qcdio.h>
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
    static int called = 0;
    char *fname="Rand()";
    char *cname="RandomGenerator";
    //---------------------------------------------------------
    // Start generating random numbers here
    //---------------------------------------------------------
//    VRB.Result(cname,fname,"inext=%d inextp=%d\n",inext,inextp);
    if ( ++inext == state_size ) inext = 0 ;
    if ( ++inextp == state_size ) inextp = 0 ;
    long mj = ma[inext] - ma[inextp] ;
    if ( mj < 0 ) mj += MBIG ;
    ma[inext] = mj ;
#if 0
    if (!UniqueID() && called %100==0)
    printf("mj = %d mj*FAC=%0.16e\n",mj,mj*FAC);
    called++;
#endif
    return mj*FAC ;
}




/*!
  \pre ::Reset must have already been called.
  \return A random number from a gaussian distribution with zero mean.
*/
IFloat GaussianRandomGenerator::Rand(int noexit)
{
    char *cname = "GaussianRandomGenerator";
    char *fname = "Rand()";
//    VRB.Result(cname,fname,"noexit=%d iset=%d\n",noexit,iset);
    if(iset == 0) {	// We don't have an extra deviate handy
        int num_try = 1;
	IFloat v1, v2, rsq;
        do {
//	    VRB.Result(cname,fname,"v1 = 2.0 * RandomGenerator::Rand() - 1.0;\n");
	    v1 = 2.0 * RandomGenerator::Rand() - 1.0;
//	    VRB.Result(cname,fname,"v2 = 2.0 * RandomGenerator::Rand() - 1.0;\n");
	    v2 = 2.0 * RandomGenerator::Rand() - 1.0;
            if ((num_try %1000)==0) 
	      VRB.Result(cname,fname,"num_try=%d v1=%e v2=%e\n",cname,fname,num_try,v1,v2);
	    rsq = v1*v1 + v2*v2;
            num_try++;
	} while((num_try<10000) &&(rsq >= 1.0 || rsq == 0) );
        if (num_try >9999) {
	  if(noexit){
	      fprintf(stderr,"%s::%s: failed after 10000 tries (corrupted RNG?), returning ridiculous numbers (1e+10)\n",cname,fname);
             gset=1e+10; iset=1;
             return 1e+10;
          }
    	  else 
	  ERR.General(cname,fname,"num_try=%d rsq>1.0",num_try);
 	}
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

/*!
  Generates a uniform random number from a distribution determined by
  the function parameters.
  \param high The upper distribution bound 
  \param low The lower distribution bound 
  \return A random number from a uniform distribution.
  \note The user must ensure that the lower bound is not greater than the
  upper bound.
*/
IFloat UniformRandomGenerator::Rand(Float high, Float low)
{
    return low + (high - low) * RandomGenerator::Rand();
}



CPS_END_NAMESPACE
