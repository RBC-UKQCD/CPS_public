#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// alg_pbp.h
//
// Header file for the AlgPbp class.
//
// AlgPbp is derived from Alg and is relevant to the 
// stochastic measurement of PsiBar Psi using the
// Conjugate Gradient algorithm. The type of fermion is
// determined by the argument to the constructor.
//
// PsiBarPsi is normalized so that for large values of the
// PbpArg.mass  PsiBarPsi =  1 / mass for any fermion type.
// This normalization results to the following small mass
// behavior for a trivial background gauge field with periodic
// boundary conditions:
// Staggered = 16 / ( Volume * mass )
// Wilson    =  1 / ( Volume * mass )
// Dwf       =  1 / ( Volume * mass )
//
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_PBP_H
#define INCLUDED_ALG_PBP_H

CPS_END_NAMESPACE
#include<util/lattice.h>
#include<util/smalloc.h>
#include<util/pmalloc.h>
#include<alg/alg_base.h>
#include<alg/common_arg.h>
#include<alg/pbp_arg.h>
CPS_START_NAMESPACE


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
