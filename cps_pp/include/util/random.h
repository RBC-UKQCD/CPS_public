#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of RNG classes.

  $Id: random.h,v 1.3 2004-01-13 20:38:57 chulwoo Exp $
 */
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-01-13 20:38:57 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/random.h,v 1.3 2004-01-13 20:38:57 chulwoo Exp $
//  $Id: random.h,v 1.3 2004-01-13 20:38:57 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2.10.1  2003/11/25 22:53:17  cwj
//  *** empty log message ***
//
//  Revision 1.2  2003/07/24 16:53:53  zs
//  Addition of documentation via doxygen:
//  doxygen-parsable comment blocks added to many source files;
//  New target in makefile and consequent alterations to configure.in;
//  New directories and files under the doc directory.
//
//  Revision 1.6  2002/12/04 17:16:27  zs
//  Merged the new 2^4 RNG into the code.
//  This new RNG is implemented in the LatRanGen class.
//  The following algorithm and utility classes are affected:
//
//  AlgEig                  Fdwf
//  AlgGheatBath            Fstag
//  AlgHmd                  GlobalJobParameter
//  AlgNoise                Lattice
//  AlgPbp                  Matrix
//  AlgThreept              RandomGenerator
//                          Vector
//
//  Revision 1.5  2001/08/16 10:50:30  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:18  anj
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
//  Revision 1.2  2001/05/25 06:16:09  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: random.h,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/random.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


#ifndef INCLUDED_RANDOM_H
#define INCLUDED_RANDOM_H              //!< Prevent multiple inclusion

CPS_END_NAMESPACE
#include <util/data_types.h>
#include <util/vector.h>
#include <util/enum.h>
CPS_START_NAMESPACE
//---------------------------------------------------------------
//! A random number generator generating uniform random numbers in (0,1)
/*!
  This uses the Fibonacci RNG routine (ran3) from Numerical Recipes in C
  pp.283.
*/
//---------------------------------------------------------------
class RandomGenerator {
  private:
    int ma[55];	// The value 55(range ma[0...54])
    				// is special and should not be
				// modified.
    int inext;
    int inextp;

  public:

    RandomGenerator();

    //! Gets a random number 
    virtual IFloat Rand(void);

    //! Seeds the RNG.
    void Reset(int seed = SERIAL_SEED);

    //! Stores the RNG state.
    void StoreSeeds(unsigned *to);

    //! Loads the RNG state.
    void RestoreSeeds(unsigned *from);

};


inline RandomGenerator::RandomGenerator() {}



//---------------------------------------------------------------
//! A random number generator generating uniform random numbers.
/*!
  The range is defined in the constructor.
  
  This should be inherited by UGrandomGenerator class rather than used
  directly.
*/
//---------------------------------------------------------------
class UniformRandomGenerator
: public virtual RandomGenerator
{
  private:
    IFloat A;
    IFloat B;

  public:
    UniformRandomGenerator(IFloat high_limit = 0.5,
			IFloat low_limit = -0.5);

//! Sets the interval over which the uniform distribution is defined.
/*!
  The default range is (-0.5, 0.5).
  \param high_limit the upper bound of the distribution range
  \param lower_limit the lower bound of the distribution range
*/
    void SetInterval(IFloat high_limit, IFloat low_limit)
                {   A = low_limit;  B = high_limit; }
    IFloat Rand(void);
};
 
/*!
  The default range is (-0.5, 0.5).
  \param high_limit the upper bound of the distribution range
  \param lower_limit the lower bound of the distribution range
*/
inline UniformRandomGenerator
::UniformRandomGenerator(IFloat high_limit, IFloat low_limit)
  : RandomGenerator(), A(low_limit), B(high_limit) {}
 


//---------------------------------------------------------------
// Gaussian Random Generator. Generate Exp(-x^2/(2*sigma2))
// 
//! A random number generator generating gaussian random numbers.
/*!
  The mean of the distribution is 0; the variance is defined in the
  constructor. It is based on the Box-Muller algorithm (see "Numerical
  Recipes in C" pp.289).

  This should be inherited by UGrandomGenerator class rather than used
  directly.
*/ 
//---------------------------------------------------------------

class GaussianRandomGenerator
: public virtual RandomGenerator
{
  private:
    IFloat sigma2;
    int iset;			 // flag
    IFloat gset;			 // saved random number

  public:
    GaussianRandomGenerator(IFloat s2 = 1.0); 	// s2 is the sigma^2
    //! Sets the variance of the distribution.
/*!
  The default variance is 1.0.
  \param s2 the variance of the gaussian distribution.
*/
    void SetSigma(IFloat s2) {sigma2 = s2; }
    IFloat Rand(void);
};
 
/*!
  The default variance is 1.0.
  \param s2 the variance of the gaussian distribution.
*/
inline GaussianRandomGenerator
::GaussianRandomGenerator(IFloat s2 /* = 1.0 */)
  : RandomGenerator(), sigma2(s2), iset(0) {}


//---------------------------------------------------------------
//! The random number generator for a single 2^4 hypercube in the lattice.
/*!
  For each 2^4 hypercube there is a uniform RNG and a gaussian RNG.
  They should be accessed through the LatRanGen class rather than directly.
*/
//  LatRanGen possesses an array of these generators,
//  one for each 2^4 hypercube.
//---------------------------------------------------------------
class UGrandomGenerator
: public UniformRandomGenerator,
    public GaussianRandomGenerator
{
  public:
    UGrandomGenerator();

    //! This should not be used.
    IFloat Rand(void); // This will return an error message - not valid option

    //! Get a gaussian random number
    IFloat Grand(void) { return GaussianRandomGenerator::Rand(); }

    //! Get a uniform random number.
    IFloat Urand(void) { return UniformRandomGenerator::Rand(); }

//! Sets the interval over which the uniform distribution is defined.
/*!
  The default range is (-0.5, 0.5).
  \param high the upper bound of the distribution range
  \param lower the lower bound of the distribution range
*/
    void SetInterval(IFloat high, IFloat low)
	{UniformRandomGenerator::SetInterval(high, low); }

//! Sets the variance of the distribution.
/*!
  The default variance is 1.0 (and the mean is zero).
  \param s2 the variance of the gaussian distribution.
*/
    void SetSigma(IFloat s2) { GaussianRandomGenerator::SetSigma(s2);  }
};

inline UGrandomGenerator
::UGrandomGenerator()
  :  UniformRandomGenerator(), GaussianRandomGenerator() {}

//---------------------------------------------------------------
//! The lattice random number generator.
/*!
  This class contains a uniform and a gaussian RNG for each 2^4 hypercube
  on the lattice.
  To ensure cross-platform reproducibility, these RNGs should be used, not the
  ones defined in the classes from which this inherits.
*/
//---------------------------------------------------------------
class LatRanGen
{
  private:
    static int n_rgen;  // Gives the number of generators (and hypercubes)
    static int rgen_pos;// ID of the generator being used
    static int can[4];  // Needed for Canonical Assignment of Sites
    static int hx[4];   // int (GJP.(X,Y,Z,T)nodeSites / 2)
    static int is_initialized; // = 0 when LatRanGen is not initialized
                                // = 1 when LatRanGen is initialized
    static UGrandomGenerator *ugran;
  public:
    LatRanGen();
    ~LatRanGen() {}
    void Initialize();  // Identical to the Constructor

    //! Get a uniform random number.
    IFloat Urand(void);

    //! Get a gaussian random number
    IFloat Grand(void);

    //! Get a uniform random number which is the same on all nodes.
    IFloat Lrand(void); 

    //! Sets the variance of the distribution.
    void SetSigma(IFloat sigma);

    //! Sets the interval over which the uniform distribution is defined.
    void SetInterval(IFloat high, IFloat low);

    //! Specifies which hypercube RNG to use.
    void AssignGenerator(int x, int y, int z, int t);
    //! Specifies which hypercube RNG to use.
    void AssignGenerator(int * coor);
    //! Specifies which hypercube RNG to use.
    void AssignGenerator(int i);
};

/*! An instance of the LatRanGen class, named LRG, should be
  created at the highest scope (outside main). This external declaration
  allows control of and access to the random number generation.
*/
extern LatRanGen LRG;

#endif


CPS_END_NAMESPACE
