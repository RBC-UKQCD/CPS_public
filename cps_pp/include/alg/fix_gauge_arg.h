#include<config.h>
CPS_START_NAMESPACE
//  fix_gauge_arg.h

#ifndef INCLUDED_FIX_GAUGE_ARG_H
#define INCLUDED_FIX_GAUGE_ARG_H

CPS_END_NAMESPACE
#include<util/enum.h>
#include<util/vector.h>
CPS_START_NAMESPACE

struct FixGaugeArg {
  FixGaugeType fix_gauge_kind;   // The kind of gauge fixing
  int hyperplane_start;          // The full lattice coordinate of the first 
                                 // hyperplane of gauge fixing matrices.
  int hyperplane_step;           // The coordinate step between hyperplanes.
  int hyperplane_num;            // The number of hyperplanes.
  Float stop_cond;               // The stopping condition.
  int max_iter_num;              // Maximum number of iterations.
                                 // It exits if reached. If is set to 0 then 
                                 // the number of interations is not checked
                                 // and will continue untill the stopping 
                                 // condition is satisfied.
};


#endif
CPS_END_NAMESPACE
