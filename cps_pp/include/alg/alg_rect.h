#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------------------
//
// alg_rect.h
//
// Header file for the AlgRect class.
//
// AlgRect is derived from Alg and it measures the average
// value of the rectangle, the 1x2 planar Wilson loop.
//
//------------------------------------------------------------------------------

#ifndef INCLUDED_ALG_RECT_H
#define INCLUDED_ALG_RECT_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/no_arg.h>
CPS_START_NAMESPACE

class AlgRect : public Alg
{
 private:
    char *cname;

    NoArg *alg_rect_arg;
        // The argument structure for the rectangle
 
    Float norm_fac;       
        // normalization factor

 public:
    AlgRect(Lattice & latt, CommonArg *c_arg, NoArg *arg);

    virtual ~AlgRect();

    void run(void);
};

#endif

CPS_END_NAMESPACE
