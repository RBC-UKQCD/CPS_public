#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// alg_ghb.h
//
// Header file for the AlgGheatBath class.
//
// AlgGheatBath is derived from Alg and is relevant to the 
// gauge heat bath algorithm. The type of gauge action is
// determined by the argument to the constructor.
//
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_GHB_H
#define INCLUDED_ALG_GHB_H

CPS_END_NAMESPACE
#include<util/lattice.h>
#include<util/smalloc.h>
#include<util/pmalloc.h>
#include<alg/alg_base.h>
#include<alg/common_arg.h>
#include<alg/ghb_arg.h>
CPS_START_NAMESPACE

class AlgGheatBath : public Alg
{
 private:
    char *cname;

    GhbArg *alg_ghb_arg;
        // The argument structure for the GheatBath algorithm

    void relocate();
    void preserve_seed();
    void UpdateLink(Matrix* mp, const Matrix & stap);

 public:
    AlgGheatBath(Lattice & latt, CommonArg *c_arg, GhbArg *arg);

    virtual ~AlgGheatBath();

    void run(void);
    void NoCheckerBoardRun();
    void NodeCheckerBoardRun();
};

#endif




CPS_END_NAMESPACE
