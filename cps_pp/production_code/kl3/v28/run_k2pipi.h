// -*- mode:c++; c-basic-offset:4 -*-
#ifndef INCLUDED_RUN_K2PIPI_H_KL3
#define INCLUDED_RUN_K2PIPI_H_KL3

#include <vector>
#include <string>
#include <alg/qpropw.h>
#include "prop_container.h"
#include "run_wsnk.h"

CPS_START_NAMESPACE

// run I=2 K to pi pi.
void run_k2pipi(const AllProp &sprop,
                const AllProp &uprop,
                const AllProp &dprop,
                const std::string &fn,
                PROP_TYPE ptype);

CPS_END_NAMESPACE

#endif
