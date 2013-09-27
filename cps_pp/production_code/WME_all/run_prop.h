// -*- mode:c++; c-basic-offset:4 -*-
#ifndef INCLUDED_RUN_PROP_H_KL3
#define INCLUDED_RUN_PROP_H_KL3

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <string>
#include <vector>
#include <cassert>

//#include <util/lattice.h>
#include <alg/qpropw_arg.h>
#include <alg/array_arg.h>
#include <alg/eigcg_arg.h>

#include "prop_container.h"

namespace cps {
    class Lattice;
};

void run_wall_prop(const char *pname,
		   AllProp *prop_e,
                   AllProp *prop,
                   cps::IntArray &eloc, 
                   cps::Lattice &lat,
                   cps::QPropWArg &qp_arg,
                   cps::EigCGArg *eigcg_arg,
                   int traj,
                   bool do_mres);

void run_mom_prop(const char *pname,
		  AllProp *prop_e,
                  AllProp *prop,
                  cps::IntArray &eloc,
                  cps::Lattice &lat,
                  cps::QPropWArg &qp_arg,
                  cps::EigCGArg *eigcg_arg,
                  int traj,
                  const int mom[3]);

#endif
