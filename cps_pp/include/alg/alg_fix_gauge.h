#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// alg_fix_gauge.h
//
// Header file for the AlgFixGauge class.
//
// AlgFixGauge is derived from Alg and is relevant to the 
// gauge fixing algorithms. The type of glue and fermion is
// determined by the argument to the constructor.
//
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_FIX_GAUGE_H
#define INCLUDED_ALG_FIX_GAUGE_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/fix_gauge_arg.h>
CPS_START_NAMESPACE


class AlgFixGauge : public Alg
{
 private:
    char *cname;

    FixGaugeArg *alg_fix_gauge_arg;
        // The argument structure for the AlgFixGauge algorithm
 

 public:
    AlgFixGauge(Lattice & latt, CommonArg *c_arg, FixGaugeArg *arg);

    virtual ~AlgFixGauge();

    void run(void);  
       // Allocates memory and constructs the gauge fixing matrices

    void free(void);
       // Free the memory of the gauge fixing matrices.
};




#endif





CPS_END_NAMESPACE
