#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// alg_plaq.h
//
// Header file for the AlgPlaq class.
//
// AlgPlaq is derived from Alg and it measures the average
// value of the plaquette. 
//
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_PLAQ_H
#define INCLUDED_ALG_PLAQ_H

CPS_END_NAMESPACE
#include<util/lattice.h>
#include<util/smalloc.h>
#include<util/pmalloc.h>
#include<alg/alg_base.h>
#include<alg/common_arg.h>
#include<alg/no_arg.h>
CPS_START_NAMESPACE


class AlgPlaq : public Alg
{
 private:
    char *cname;

    NoArg *alg_plaq_arg;
        // The argument structure for the plaquette
 
    Float norm_fac;       
        // normalization factor

 public:
    AlgPlaq(Lattice & latt, CommonArg *c_arg, NoArg *arg);

    virtual ~AlgPlaq();

    void run(void);
};



#endif




CPS_END_NAMESPACE
