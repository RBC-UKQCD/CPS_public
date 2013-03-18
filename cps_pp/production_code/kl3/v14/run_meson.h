// -*- mode:c++; c-basic-offset:4 -*-
#ifndef INCLUDED_RUN_MESON_H_KL3
#define INCLUDED_RUN_MESON_H_KL3

#include <vector>
#include <alg/qpropw.h>
#include "prop_container.h"

//! IMPORTANT: The function apply_op_src() and apply_op_snk() depend
//! on the specific order the following values are organized. DON'T
//! CHANGE THE ORDERING OF THE FOLLOWING ENUMERATE VALUES.
enum Operator {
    GAMMA_0, GAMMA_1, GAMMA_2, GAMMA_3, GAMMA_5,
    ID,
    GAMMA_05, GAMMA_15, GAMMA_25, GAMMA_35,
    GAMMA_50, GAMMA_51, GAMMA_52, GAMMA_53,
};

// compute meson correlation functions point sink, source type is
// determined implicitly by the propagator.
void run_meson_pt(const AllProp &propA,
                  const AllProp &propB,
                  Operator snk_op,
                  Operator src_op,
                  const char fn_t[], int traj,
                  PROP_TYPE ptype);

// compute meson correlation functions wall sink, source type is
// determined implicitly by the propagator.
void run_meson_wall(const AllProp &propA,
                    const AllProp &propB,
                    Operator snk_op,
                    Operator src_op,
                    const char fn_t[], int traj,
                    PROP_TYPE ptype);

// compute meson correlation functions wall sink, source type is
// determined implicitly by the propagator.
//
// Disconnected diagram
void run_meson_dis_wall(const AllProp &propL,
                        const AllProp &propR,
                        Operator l_op,
                        Operator r_op,
                        const char fn_t[], int traj,
                        PROP_TYPE ptype);

#endif
