#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// alg_mom.h
//
// Header file for the AlgMom class.
//
// AlgMom it calculates the phase factor for each
// lattice site given a number of momenta and the source parameters
//
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_MOM_H
#define INCLUDED_ALG_MOM_H

CPS_END_NAMESPACE
#include <math.h>    // for cos and sin
#include <alg/mom_arg.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
CPS_START_NAMESPACE


class AlgMom
{
 private:
    char *cname;
    MomArg *alg_mom_arg; // argument structure for momentum states

    Float PI;
    int dir,i,j,k;      // propagation direction and 3 orthogonal
    int no_of_mom;      // number of momenta to be calculated
    int deg;            // calculate degenerate momenta separately/together

    int nx[4];                // local lattice extent
    int glb_L[4];             // global lattice extent
    int glb_sour_center[4];   // global source location
    Complex *mom_fact;

 public:
    AlgMom(CommonArg *c_arg, MomArg *arg);
    virtual ~AlgMom();
    void run(void);
        
    // returns the complex phase factor for momentum "imom" at site "s"
    Complex fact(int imom, int *s);
};



#endif





CPS_END_NAMESPACE
