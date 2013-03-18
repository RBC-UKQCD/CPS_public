// -*- mode:c++; c-basic-offset:4 -*-
#ifndef INCLUDED_RUN_2PION_H_KL3
#define INCLUDED_RUN_2PION_H_KL3

#include <vector>
#include <alg/qpropw.h>
#include "prop_container.h"

// \pi^+ \pi^+ 2 pion scattering, D and C diagrams only.
// 
// no separation 
//
// pd: momentum in the d quark propagator
void run_2pionDC(const AllProp &uprop,
                 const AllProp &dprop,
                 const char fn_t[], int traj,
                 PROP_TYPE ptype, const int pd[3]);

#endif
