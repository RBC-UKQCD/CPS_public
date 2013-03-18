// -*- mode:c++; c-basic-offset:4 -*-
#ifndef INCLUDED_RUN_PION_H_KL3
#define INCLUDED_RUN_PION_H_KL3

#include <vector>
#include <alg/qpropw.h>
#include "prop_container.h"

// compute meson correlation functions point sink, source type is
// determined implicitly by the propagator.
void run_meson_pt(const AllProp &propA,
                  const AllProp &propB,
                  const char fn_t[], int traj);

// compute meson correlation functions wall sink, source type is
// determined implicitly by the propagator.
void run_meson_wall(const AllProp &propA,
                    const AllProp &propB,
                    const char fn_t[], int traj);

// compute meson correlation functions wall sink, source type is
// determined implicitly by the propagator.
//
// Disconnected diagram
void run_meson_dis_wall(const AllProp &propL,
                        const AllProp &propR,
                        const char fn_t[], int traj);

// compute kl3 correlation functions, source type is
// determined implicitly by the propagator.
//
// This version computes all possible K-pi separations
void run_kl3(const AllProp &sprop,
             const AllProp &lprop,
             const AllProp &ltwst,
             const char fn_t[], int traj);

#endif
