#include<config.h>
CPS_START_NAMESPACE
/*  threept_arg.h */

/*  The structure type ThreePtArg holds the parameters specific to
    meson three point functions for Wilson type fermion's  */

#ifndef INCLUDED_3PT_ARG_H
#define INCLUDED_3PT_ARG_H

CPS_END_NAMESPACE
#include<alg/cg_arg.h>
#include<util/vector.h>
CPS_START_NAMESPACE

struct ThreePtArg {
  /* ??? */

  CgArg cg;		/* The conjugate gradient argument for 
			   the quark propagator inversion */
  int seed;		/* seed for random source */

  int t_src;		// src time for quarks in I graphs

  int t_Op;		// the operator time slice

  int t_Op_2;		// the 2nd operator time slice

  int t_sink;		// sink time for spectator quark in I graphs

  int num_masses;	// number of masses to do

  Float mass[20];	// the list of masses

};

#endif /* !INCLUDED_3PT_ARG_H */
CPS_END_NAMESPACE
