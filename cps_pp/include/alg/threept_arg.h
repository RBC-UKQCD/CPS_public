/*  threept_arg.h */

/*  The structure type ThreePtArg holds the parameters specific to
    meson three point functions for Wilson type fermion's  */

#ifndef INCLUDED_3PT_ARG_H
#define INCLUDED_3PT_ARG_H

#include <alg/cg_arg.h>
#include <util/vector.h>

struct ThreePtArg {
  /* ??? */

  CgArg cg;		/* The conjugate gradient argument for 
			   the quark propagator inversion */
  int seed;		/* seed for random source */

  int tOpStart;		// the first operator time slice

  int tOpEnd;		// the last operator time slice

  int t_src[4];	        // source times for quarks

  int num_masses;	// number of masses to do

  int num_charm_masses;// number of masses to do

  int num_hits;         // number of random sources to do

  Float mass[10];	// the list of light masses

  Float charm_mass[10];	// the list of charm masses
  
  unsigned lat_check_sum; // lattice checksum for prop tagging.
};

#endif /* !INCLUDED_3PT_ARG_H */
