#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/random.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: random.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
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
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/random.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// random.h
// Header file for the implementation of a random number generator
//
//  Implemented in "random_asm.asm" and "random.C"
//

//////////////////////////////////////////////////////////////////
//								//
//		Fibonacci Random number generator		//
//								//
//////////////////////////////////////////////////////////////////


#ifndef INCLUDED_RANDOM_H
#define INCLUDED_RANDOM_H

CPS_END_NAMESPACE
#include<util/data_types.h>
#include<util/vector.h>
#include<util/enum.h>
CPS_START_NAMESPACE
//---------------------------------------------------------------
// Generate uniform random number in (0,1)
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
    virtual IFloat Rand(void);
    void Reset(int seed = SERIAL_SEED);
	// This function should be called before any random
	// generator is used.

    void StoreSeeds(unsigned *to);
    void RestoreSeeds(unsigned *from);
    	// to[] anf from[] are buffers of size 55+2
};


inline RandomGenerator::RandomGenerator() {}



//---------------------------------------------------------------
// Generate uniform random number in (a,b)
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
    void SetInterval(IFloat high_limit, IFloat low_limit)
                {   A = low_limit;  B = high_limit; }
    IFloat Rand(void);
};
 

inline UniformRandomGenerator
::UniformRandomGenerator(IFloat high_limit, IFloat low_limit)
  : RandomGenerator(), A(low_limit), B(high_limit) {}
 


//---------------------------------------------------------------
// Gaussian Random Generator. Generate Exp(-x^2/(2*sigma2))
// see "Numerical Recipes in C" pp.289
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
    void SetSigma(IFloat s2) {sigma2 = s2; }
    IFloat Rand(void);
};
 

inline GaussianRandomGenerator
::GaussianRandomGenerator(IFloat s2 /* = 1.0 */)
  : RandomGenerator(), sigma2(s2), iset(0) {}


//---------------------------------------------------------------
//  This is a random number generator for a single 2^4 hypercube
//  in the lattice.  LatRanGen possesses an array of these generators,
//  one for each 2^4 hypercube.
//---------------------------------------------------------------
class UGrandomGenerator
: public UniformRandomGenerator,
    public GaussianRandomGenerator
{
  public:
    UGrandomGenerator();
    IFloat Rand(void); // This will return an error message - not valid option
    IFloat Grand(void) { return GaussianRandomGenerator::Rand(); }
    IFloat Urand(void) { return UniformRandomGenerator::Rand(); }
    void SetInterval(IFloat high, IFloat low)
                {UniformRandomGenerator::SetInterval(high, low); }
    void SetSigma(IFloat s2) { GaussianRandomGenerator::SetSigma(s2);  }
};

inline UGrandomGenerator
::UGrandomGenerator()
  :  UniformRandomGenerator(), GaussianRandomGenerator() {}

//---------------------------------------------------------------
// Possesses a unique random number generator for each 2^4
// hypercube in the lattice.
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
    IFloat Urand(void);
    IFloat Grand(void);
    IFloat Lrand(void); // Returns same random number on all nodes
    void SetSigma(IFloat sigma);
    void SetInterval(IFloat high, IFloat low);
    void AssignGenerator(int x, int y, int z, int t);
    void AssignGenerator(int * coor);
    void AssignGenerator(int i);
};

extern LatRanGen LRG;

#endif

CPS_END_NAMESPACE
