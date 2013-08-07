// -*- mode:c++; c-basic-offset:4 -*-
#ifndef INCLUDED_RUN_2PION_H_KL3
#define INCLUDED_RUN_2PION_H_KL3

#include <vector>
#include <string>
#include <alg/qpropw.h>
#include "prop_container.h"
#include "run_wsnk.h"

CPS_START_NAMESPACE

// 1. 2 pion contractions, D and C diagram only (for I=2 pi pi
// scattering), no separation.
//
// 2. Assuming cosine source with antiperiodic boundary condition, use
// momentum sink.
//
// 3. pd: momentum in the d quark propagator.
//
// FIXME: there may be one or more factors of 2 in the sink side.
void run_2pionDC(const AllProp &uprop,
                 const AllProp &dprop,
                 const std::string &fn,
                 PROP_TYPE ptype, const int pd[3]);

CPS_END_NAMESPACE

#endif
