#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// alg_inst.h
//
// Header file for the AlgInst class.
//
// AlgInst is derived from Alg and is relevant to the 
// generation of an instanton configuration.
//
//----- Begin Adrian
// ???
//
//      This operator creates a single instanton at
//      the center of the lattice.  It is aware of the
//      overall size of the lattice, and uses this, to
//      place the instanton at the actual center of the
//      whole lattice, regardless of the number of processors.
//        The instanton requires two arguments.
//              0       'Instanton Size'
//              1       average value of delta-nu.  If this
//                       argument is +1 or -1 then the move
//                       is forced, if it is 0, then the
//                       move is in a random direction and
//                       there is an accept reject step.
//              2       For squashed instantons, this is
//                        the value of "r_max"
//      The return value is the change in the action of the 
//      configuration as a result of the operation.
//
//        The actual 'instanton' configuration is the solution
//      outlined in Coleman(p297):
//              g1 = (1/r) * (x_4 + i x_dot_sigma )
//              f(r) = r^2/(r^2/rho^2)
//              A_mu = f(r) g1 d_mu [g1]^-1
//
//----- End Adrian
//
//------------------------------------------------------------------

#ifndef INCLUDED_INST_PBP_H
#define INCLUDED_INST_PBP_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/inst_arg.h>
CPS_START_NAMESPACE


class AlgInst : public Alg
{
 private:
    char *cname;

    InstArg *alg_inst_arg;
        // The argument structure for the inst algorithm

 public:
    AlgInst(Lattice & latt, CommonArg *c_arg, InstArg *arg);

    virtual ~AlgInst();

    void run(void);
};

#endif





CPS_END_NAMESPACE
