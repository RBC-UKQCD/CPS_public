#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// alg_eig.h
//
// Header file for the AlgEig class.
//
// AlgEig is derived from Alg and is relevant to the 
// Ritz eigenvector solver. The type of fermion is
// determined by the argument to the constructor.
//
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_EIG_H
#define INCLUDED_ALG_EIG_H

CPS_END_NAMESPACE
#include<util/lattice.h>
#include<util/smalloc.h>
#include<util/pmalloc.h>
#include<alg/alg_base.h>
#include<alg/common_arg.h>
#include<alg/eig_arg.h>
CPS_START_NAMESPACE


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

 public:
    AlgEig(Lattice & latt, CommonArg *c_arg, EigArg *arg);

    virtual ~AlgEig();

    void run(void);
};

#endif




CPS_END_NAMESPACE
