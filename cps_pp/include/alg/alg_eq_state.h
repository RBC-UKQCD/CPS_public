#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// alg_eq_state.h
//
// Header file for the AlgEqState class.
//
// AlgEqState is derived from Alg and it measures the sum, normalized
// by volume, of the plaquette on particular hyperplane(s).
//
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_EQ_STATE_H
#define INCLUDED_ALG_EQ_STATE_H

CPS_END_NAMESPACE
#include<util/lattice.h>
CPS_START_NAMESPACE
// Reduntant  #include "../../util/include/smalloc.h"
// Reduntant  #include "../../util/include/pmalloc.h"
CPS_END_NAMESPACE
#include<alg/alg_base.h>
#include<alg/common_arg.h>
#include<alg/eq_state_arg.h>
CPS_START_NAMESPACE


class AlgEqState : public Alg
{
 private:
    char *cname;

    EqStateArg *alg_eq_state_arg;
        // The argument structure for AlgEqState
 
    Float norm_fac;       
        // normalization factor

 public:
    AlgEqState(Lattice & latt, CommonArg *c_arg, EqStateArg *arg);

    virtual ~AlgEqState();

    void run(void);
};



#endif




CPS_END_NAMESPACE
