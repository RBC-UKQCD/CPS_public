#ifndef INCLUDED_ALG_REMEZ_H
#define INCLUDED_ALG_REMEZ_H

#include<config.h>

#ifdef GMP        // If GMP is defined 

#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/common_arg.h>
#include <alg/bigfloat.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// alg_remez.h
//
// AlgRemez is relevant to the Rational Hybrid Molecular Dynamics
// algorithm.  The Remez algorithm is used for generating the nth root
// rational approximation.
//
// Note this class can only be used if
// the gnu multiprecision (GNU MP) library is present.
//
//------------------------------------------------------------------

#define JMAX 1000 //Maximum number of iterations of Newton's approximation

class AlgRemez
{
 private:
  char *cname;

  // The approximation parameters
  bigfloat *param, *roots, *poles;
  bigfloat norm;

  // The numerator and denominator degree (n=d)
  int n, d;
  
  // The bounds of the approximation
  bigfloat apstrt, apwidt, apend;

  // the numerator and denominator of the power we are approximating
  unsigned long power_num; 
  unsigned long power_den;

  // Flag to determine whether the arrays have been allocated
  int alloc;

  // Variables used to calculate the approximation
  int nd1, iter;
  bigfloat *xx, *mm, *step;
  bigfloat delta, spread, tolerance;

  // The number of equations we must solve at each iteration (n+d+1)
  int neq;

  // The precision of the GNU MP library
  long prec;

  // Initial values of maximal and minmal errors
  void initialGuess();

  // Solve the equations
  void equations();

  // Search for error maxima and minima
  void search(bigfloat *step); 

  // Initialise step sizes
  void stpini(bigfloat *step);

  // Calculate the roots of the approximation
  int root();

  // Evaluate the polynomial
  bigfloat polyEval(bigfloat x, bigfloat *poly, long size);

  // Evaluate the differential of the polynomial
  bigfloat polyDiff(bigfloat x, bigfloat *poly, long size);

  // Newton's method to calculate roots
  bigfloat rtnewt(bigfloat *poly, long i, bigfloat x1, bigfloat x2, bigfloat xacc);

  // Evaluate the partial fraction expansion of the rational function
  // with res roots and poles poles.  Result is overwritten on input
  // arrays.
  void pfe(bigfloat *res, bigfloat* poles, bigfloat norm);

  // Evaluate the rational form P(x)/Q(x) using coefficients from the
  // solution vector param
  bigfloat approx(bigfloat x);

  // Calculate function required for the approximation
  bigfloat func(bigfloat x);

  // Compute size and sign of the approximation error at x
  bigfloat getErr(bigfloat x, int *sign);

  // Solve the system AX=B
  int simq(bigfloat *A, bigfloat *B, bigfloat *X, int n);

  // Free memory and reallocate as necessary
  void allocate(int num_degree, int den_degree);
 public:
  
  // Constructor
  AlgRemez(Float lower, Float upper, long prec);

  // Destructor
  virtual ~AlgRemez();

  // Reset the bounds of the approximation
  void setBounds(Float lower, Float upper);

  // Generate the rational approximation x^(pnum/pden)
  Float generateApprox(int num_degree, int den_degree, unsigned long power_num, unsigned long power_den);
  Float generateApprox(int degree, unsigned long power_num, unsigned long power_den);

  // Return the partial fraction expansion of the approximation x^(pnum/pden)
  int getPFE(Float *res, Float *pole, Float *norm);

  // Return the partial fraction expansion of the approximation x^(-pnum/pden)
  int getIPFE(Float *res, Float *pole, Float *norm);

};
CPS_END_NAMESPACE

#else             // If not defined GMP


#include <util/error.h>
#include <util/data_types.h>
CPS_START_NAMESPACE


// Dummy class for case when gmp is not present

class AlgRemez
{
 private:
  char *cname;

 public:
  
  AlgRemez(Float lower, Float upper, long prec) {
    cname = "AlgRemez";
    char *fname = "AlgRemez(Float, Float, long)";
    ERR.General(cname,fname,"AlgRemez cannot be instantiated without GMP installed\n");
  }
  ~AlgRemez() {;}
  void setBounds(Float lower, Float upper) {;}
  Float generateApprox(int num_degree, int den_degree, unsigned long power_num, unsigned long power_den);
  Float generateApprox(int degree, unsigned long power_num, unsigned long power_den)
    {return 0;}
  int getPFE(Float *res, Float *pole, Float *norm) {return 0;}
  int getIPFE(Float *res, Float *pole, Float *norm) {return 0;}

};
CPS_END_NAMESPACE

#endif  // Ifdef GMP

#endif  // Include guard



