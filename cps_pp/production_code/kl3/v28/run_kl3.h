// -*- mode:c++; c-basic-offset:4 -*-
#ifndef INCLUDED_RUN_KL3_H_KL3
#define INCLUDED_RUN_KL3_H_KL3

#include <vector>
#include <string>
#include <alg/qpropw.h>
#include "prop_container.h"
#include "run_wsnk.h"

CPS_START_NAMESPACE

// compute kl3 correlation functions, source type is
// determined implicitly by the propagator.
//
// This version computes all possible K-pi separations
void run_kl3(const AllProp &sprop,
             const AllProp &lprop,
             const AllProp &ltwst,
             const std::string &fn,
             PROP_TYPE ptype);

CPS_END_NAMESPACE

#endif
