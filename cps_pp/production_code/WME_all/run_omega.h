// -*- mode:c++; c-basic-offset:4 -*-
#ifndef INCLUDED_RUN_OMEGA_H_KL3
#define INCLUDED_RUN_OMEGA_H_KL3

#include <vector>
#include <string>
#include <alg/qpropw.h>
#include "prop_container.h"
#include "my_util.h"

CPS_START_NAMESPACE

// compute omega correlation functions point sink, source type is
// determined implicitly by the propagator.
void run_omega_pt(const AllProp &prop,
                  Operator op,
                  const std::string &fn,
                  PROP_TYPE ptype);

CPS_END_NAMESPACE

#endif
