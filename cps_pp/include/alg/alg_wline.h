#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// alg_wline.h
//
// Header file for the AlgWline class.
//
// AlgWline is derived from Alg and it measures the average
// value of the Wilson line for each direction. 
//
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_WLINE_H
#define INCLUDED_ALG_WLINE_H

CPS_END_NAMESPACE
#include<util/lattice.h>
#include<util/smalloc.h>
#include<util/pmalloc.h>
#include<alg/alg_base.h>
#include<alg/common_arg.h>
#include<alg/no_arg.h>
CPS_START_NAMESPACE


class AlgWline : public Alg
{
 private:
    char *cname;

    NoArg *alg_wline_arg;
        // The argument structure for the plaquette
 
 public:
    AlgWline(Lattice & latt, CommonArg *c_arg, NoArg *arg);

    virtual ~AlgWline();

    void run(void);
};



#endif




CPS_END_NAMESPACE
