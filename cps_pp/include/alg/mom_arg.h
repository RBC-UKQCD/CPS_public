#include<config.h>
CPS_START_NAMESPACE
//  mom_arg.h

#ifndef INCLUDED_MOM_ARG_H
#define INCLUDED_MOM_ARG_H

struct MomArg {
  int no_of_momenta; // number of different momenta    
  int deg;           // control flag: average over degenerate momenta on/off
  int dir;           // propagation direction
  int src_begin[4];  // source 
  int src_end[4];

};


#endif
CPS_END_NAMESPACE
