// -*- mode:c++; c-basic-offset:4 -*-
#ifndef INCLUDED_RUN_KL3_H_KL3
#define INCLUDED_RUN_KL3_H_KL3

#include <vector>
#include <alg/qpropw.h>
#include "prop_container.h"

// compute kl3 correlation functions, source type is
// determined implicitly by the propagator.
//
// This version computes all possible K-pi separations
void run_kl3(const AllProp &sprop,
             const AllProp &lprop,
             const AllProp &ltwst,
             const char fn_t[], int traj,
             PROP_TYPE ptype);

#endif
