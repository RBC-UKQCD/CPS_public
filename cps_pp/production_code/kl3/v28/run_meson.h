// -*- mode:c++; c-basic-offset:4 -*-
#ifndef INCLUDED_RUN_MESON_H_KL3
#define INCLUDED_RUN_MESON_H_KL3

#include <vector>
#include <string>
#include <alg/qpropw.h>
#include "prop_container.h"
#include "my_util.h"
#include "run_wsnk.h"

CPS_START_NAMESPACE

// compute meson correlation functions point sink, source type is
// determined implicitly by the propagator.
void run_meson_pt(const AllProp &propA,
                  const AllProp &propB,
                  Operator snk_op,
                  Operator src_op,
                  const std::string &fn,
                  PROP_TYPE ptype);

// compute meson correlation functions wall sink, source type is
// determined implicitly by the propagator.
void run_meson_wall(const AllProp &propA,
                    const AllProp &propB,
                    Operator snk_op,
                    Operator src_op,
                    const std::string &fn,
                    PROP_TYPE ptype);

// compute meson correlation functions wall sink, source type is
// determined implicitly by the propagator.
//
// Disconnected diagram
void run_meson_disc(const AllProp &propL,
                    const AllProp &propR,
                    Operator l_op,
                    Operator r_op,
                    const std::string &fn,
                    PROP_TYPE ptype);

CPS_END_NAMESPACE

#endif
