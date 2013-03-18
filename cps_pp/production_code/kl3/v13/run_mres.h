// -*- mode:c++; c-basic-offset:4 -*-
#ifndef INCLUDED_RUN_MRES_H_KL3
#define INCLUDED_RUN_MRES_H_KL3

#include <vector>
#include <alg/qpropw.h>
#include "prop_container.h"

void run_mres(const cps::QPropW &qp, unsigned tsrc, const char fn[]);
void run_za(const cps::QPropW &qp, cps::Float mass, unsigned tsrc, const char fn[]);

#endif
