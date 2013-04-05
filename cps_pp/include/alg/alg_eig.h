#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief Definitions of the AlgEig class.

  $Id: alg_eig.h,v 1.5 2013-04-05 17:46:30 chulwoo Exp $
*/
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_EIG_H
#define INCLUDED_ALG_EIG_H           //!< Prevent multiple inclusion.

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/eig_arg.h>
CPS_START_NAMESPACE

//! A class implementing a Ritz eigenvalue solver for the fermion matrix.
/*! \ingroup alg */
class AlgEig : public Alg
{
 private:
    char *cname;

    EigArg *alg_eig_arg;
        // The argument structure for the eig algorithm
 
    int Ncb;       
        // Number of checkerboards for fermion field (1 or 2)

    Vector **eigenv;
        // The eigenvectors (initial and final)

    Float *lambda;
        // The eigenvalues (final)

    Float *chirality;
        // The chirality of the eigenvalues (final)

    int *valid_eig;
        // Whether the eigenvalues are valid or not (final)

    int n_masses;
    // The number of masses in the loop

 public:
    AlgEig(Lattice & latt, CommonArg *c_arg, EigArg *arg);

    virtual ~AlgEig();

    void run(void);
    void run(Float **lambda, Vector** in_eigv=0);
    void runLanczos();//Float *eval, Vector **evec);
};

#endif





CPS_END_NAMESPACE
