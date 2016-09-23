#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definitions of the AlgPbp class.

  $Id: alg_pbp.h,v 1.3 2004/08/18 11:57:35 zs Exp $
*/
//---------------------------------------------------------------------------

#ifndef INCLUDED_ALG_PBP_H
#define INCLUDED_ALG_PBP_H           //!< Prevent multiple inclusion.

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/pbp_arg.h>
CPS_START_NAMESPACE

//---------------------------------------------------------------------------
//! Class for quark condensate calculation.
/*! 
  The condensate \f$ \bar\psi\psi \f$ (and similar things) are computed
  stochastically using the Conjugate Gradient algorithm. They can be computed
  for a number of different fermion masses.

  It is normalized so that for large values of the
  mass, \f$ \bar\psi\psi \f$ =  1/mass for any fermion type.

 This normalization results to the following small mass
 behavior for a trivial background gauge field with periodic
 boundary conditions:
- Staggered    = 16 / ( Volume * mass )
- Wilson       =  1 / ( Volume * mass )
- Domain wall  =  1 / ( Volume * mass )

  The algorithm runs roughly as follows:
-#  Create a random gaussian vector R
-#  Solve M psi = R for psi where M is the fermion matrix. The initial guess
  for psi is a vector with every complex component equal to 1.
-#  Compute the normalised real part of the dot product <R,psi>.
-#  For Wilson type quarks, compute the normalised real part of the dot
product <R, gamma_5 psi>; for staggered quarks, compute something else.

\ingroup alg 
 */
//---------------------------------------------------------------------------
class AlgPbp : public Alg
{
 private:
    char *cname;

    PbpArg *alg_pbp_arg;
        // The argument structure for the pbp algorithm
 
    int f_size;       
        // Node checkerboard size of the fermion field

    Vector *src;
        // The source vector

    Vector *sol;
        // The solution vector

 public:
    AlgPbp(Lattice & latt, CommonArg *c_arg, PbpArg *arg);

    virtual ~AlgPbp();

    void run(void);

    void runPointSource(int x, int y, int z, int t);
};

#endif





CPS_END_NAMESPACE
