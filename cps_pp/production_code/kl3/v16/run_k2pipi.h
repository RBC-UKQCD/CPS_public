// -*- mode:c++; c-basic-offset:4 -*-
#ifndef INCLUDED_RUN_K2PIPI_H_KL3
#define INCLUDED_RUN_K2PIPI_H_KL3

#include <vector>
#include <alg/qpropw.h>
#include "prop_container.h"

// run I=2 K to pi pi.
void run_k2pipi(const AllProp &sprop,
                const AllProp &uprop,
                const AllProp &dprop,
                const char fn_t[], int traj,
                PROP_TYPE ptype);

#endif
