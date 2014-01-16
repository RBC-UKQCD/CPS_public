// -*- mode:c++; c-basic-offset:4 -*-
#ifndef INCLUDED_RUN_BK_H_KL3
#define INCLUDED_RUN_BK_H_KL3

#include <vector>
#include <string>
#include <alg/qpropw.h>
#include "my_util.h"
#include "prop_container.h"

CPS_START_NAMESPACE

// compute Bk correlation functions, source type is
// determined implicitly by the propagator.
void run_bk(const AllProp &lpropA,
            const AllProp &spropA,
            const AllProp &lpropB,
            const AllProp &spropB,
            const std::string &fn,
            PROP_TYPE ptype);

CPS_END_NAMESPACE

#endif
