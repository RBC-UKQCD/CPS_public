#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of RNG classes.

  $Id: random.h,v 1.13 2004-08-17 03:33:08 chulwoo Exp $
 */
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-08-17 03:33:08 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/random.h,v 1.13 2004-08-17 03:33:08 chulwoo Exp $
//  $Id: random.h,v 1.13 2004-08-17 03:33:08 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: random.h,v $
//  $Revision: 1.13 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/random.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


#ifndef INCLUDED_RANDOM_H
#define INCLUDED_RANDOM_H              //!< Prevent multiple inclusion

CPS_END_NAMESPACE
#include <util/data_types.h>
#include <util/enum.h>
#include <util/smalloc.h>
CPS_START_NAMESPACE
//---------------------------------------------------------------
//! A random number generator generating uniform random numbers in (0,1)
/*!
  This uses the Fibonacci RNG routine (ran3) from Numerical Recipes in C
*/
//---------------------------------------------------------------
const int MBIG  = 1000000000;
const IFloat FAC = 1.0E-09;			// 1.0/MBIG
class RandomGenerator {
  private:

    static const int state_size = 55;
    int ma[state_size];	// The value 55(range ma[0...54])
    				// is special and should not be
				// modified.
    int inext;
    int inextp;
    
  public:

    RandomGenerator() {
//      FAC = 1e-9;
    }

    //! Gets a random number 
//    virtual IFloat Rand();
    IFloat Rand();

    //! Seeds the RNG.
    void Reset(int seed);

    //! Stores the RNG state.
    void StoreSeeds(unsigned int *to) const;

    //! Loads the RNG state.
    void RestoreSeeds(const unsigned int *from);

    //! Size of the RNG state.
    int StateSize() const;
    
};






//---------------------------------------------------------------
//! A random number generator generating uniform random numbers.
/*!
  The range is defined in the constructor.
  
  This should be inherited by UGrandomGenerator class rather than used
  directly.
*/
//---------------------------------------------------------------
class UniformRandomGenerator: public virtual RandomGenerator
{
  private:
    static IFloat A;
    static IFloat B;

  public:

/*!
  The default range is (-0.5, 0.5).
  \param high_limit the upper bound of the distribution range
  \param lower_limit the lower bound of the distribution range
*/
    UniformRandomGenerator(IFloat high_limit = 0.5, IFloat low_limit = -0.5):
//	RandomGenerator(), A(low_limit), B(high_limit) {}
	RandomGenerator() {}

//! Sets the interval over which the uniform distribution is defined.
/*!
  The default range is (-0.5, 0.5).
  \param high_limit the upper bound of the distribution range
  \param lower_limit the lower bound of the distribution range
*/
    void SetInterval(IFloat high_limit, IFloat low_limit){
	A = low_limit;
	B = high_limit;
    }

    IFloat Rand();
    IFloat Rand(Float high, Float low);
};
 
 


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

class GaussianRandomGenerator : public virtual RandomGenerator
{
  private:
    static IFloat sigma2;
    int iset;			 // flag
    IFloat gset;			 // saved random number

  public:
/*!
  The default variance is 1.0.
  \param s2 the variance of the gaussian distribution.
*/
    GaussianRandomGenerator(IFloat s2 = 1.0):
	RandomGenerator(), iset(0) {SetSigma(s2);}    

//! Sets the variance of the distribution.
/*!
  \param s2 the variance of the gaussian distribution.
*/
    void SetSigma(IFloat s2) {
	sigma2 = s2;
    }

    IFloat Rand();
};
 


//---------------------------------------------------------------
//! The random number generator for a single 2^4 hypercube in the lattice.
/*!
  For each 2^4 hypercube there is a uniform RNG and a gaussian RNG.
  They should be accessed through the LatRanGen class rather than directly.
*/
//  LatRanGen possesses an array of these generators,
//  one for each 2^4 hypercube.
//---------------------------------------------------------------
class UGrandomGenerator:
public UniformRandomGenerator, public GaussianRandomGenerator
{
  public:
    UGrandomGenerator():
	UniformRandomGenerator(),
	GaussianRandomGenerator() {};

    //! This should not be used.
//    IFloat Rand(); // This will return an error message - not valid option

    //! Get a gaussian random number
    IFloat Grand() { return GaussianRandomGenerator::Rand(); }

    //! Get a uniform random number.
    IFloat Urand() { return UniformRandomGenerator::Rand(); }
    IFloat Urand(Float high, Float low) 
	{ return UniformRandomGenerator::Rand(high,low); }

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

    static const int default_seed = 112319;
    
    int n_rgen;  // Gives the number of generators (and hypercubes)
    int rgen_pos;// ID of the generator being used
    int can[5];  // Needed for Canonical Assignment of Sites
    int hx[5];   // int (GJP.(X,Y,Z,T,S)nodeSites / 2)
    int is_initialized; // = 0 when LatRanGen is not initialized
                                // = 1 when LatRanGen is initialized
    UGrandomGenerator *ugran;
    int n_rgen_4d ;  // CJ: Gives the number of 4d generators (and hypercubes)
    int rgen_pos_4d;// CJ: ID of the 4D generator being used
    UGrandomGenerator *ugran_4d; // CJ: 4D RNG for gauge field

    char *cname;
    
  public:
    LatRanGen();
    ~LatRanGen();
    void Initialize();  

    //! Get a uniform random number.
    IFloat Urand(FermionFieldDimension frm_dim=FIVE_D);

    //! Get a gaussian random number
    IFloat Grand(FermionFieldDimension frm_dim=FIVE_D);

    //! Get a uniform random number which is the same on all nodes.
    IFloat Lrand(); 

    //! Sets the variance of the distribution.
    void SetSigma(IFloat sigma);

    //! Sets the interval over which the uniform distribution is defined.
    void SetInterval(IFloat high, IFloat low);

    //! Specifies which hypercube RNG to use.
    void AssignGenerator(int x, int y, int z, int t,int s = 0);
    //! Specifies which hypercube RNG to use.
    void AssignGenerator(const int * coor);
    //! Specifies which hypercube RNG to use.
    void AssignGenerator(int i);

    //! Size of the RNG state (per hypercube).
    int StateSize() const;


    //! Assign the  state to a selected RNG.
    void SetState(const unsigned int *,
        FermionFieldDimension frm_dim =FIVE_D );

    //! Assign the  state of all RNGs.      
    void SetStates(unsigned int **, 
        FermionFieldDimension frm_dim =FIVE_D );

    //! Get the total number of states
    int NStates() const;
    
    //! Retrieve the state of a single RNG
    void GetState(unsigned int *, 
        FermionFieldDimension frm_dim = FIVE_D ) const;

    //! Retrieve the state of all RNGs
    void GetStates(unsigned int **,
        FermionFieldDimension frm_dim = FIVE_D ) const;
    
	
};

/*! An instance of the LatRanGen class, named LRG, should be
  created at the highest scope (outside main). This external declaration
  allows control of and access to the random number generation.
*/
extern LatRanGen LRG;

#endif


CPS_END_NAMESPACE
