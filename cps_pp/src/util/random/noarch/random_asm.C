#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief   Methods for the Random Number Generator classes.

  $Id: random_asm.C,v 1.2 2003-07-24 16:53:54 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:54 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/random/noarch/random_asm.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Id: random_asm.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:38  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:35  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:10  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: random_asm.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/random/noarch/random_asm.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//  random.C
//---------------------------------------------------------------
//  This is the routine from Numerical Recipes in C PP.283
//	ran3
//---------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/random.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <math.h>
CPS_START_NAMESPACE



// Constants
const int   MBIG  = 1000000000;
const int   MZ    = 0;
const IFloat FAC   = 1.0e-9;			// 1.0/MBIG




/*!
  \pre ::Reset must have already been called.
  \return A random number from a uniform distribution over (0,1)
 */
IFloat RandomGenerator::Rand(void)
{
    //---------------------------------------------------------
    // Start generating random numbers here
    //---------------------------------------------------------
    if ( ++inext == 55 ) inext = 0 ;
    if ( ++inextp == 55 ) inextp = 0 ;
    int mj = ma[inext] - ma[inextp] ;
    if ( mj < MZ ) mj += MBIG ;
    ma[inext] = mj ;
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
