#include<config.h>
CPS_START_NAMESPACE
/*  hmd_arg.h  */

/*  The type definition for the structure used to hold the arguments
    for the hybrid molecular dynamics code. This structure has entries
    that accomodate both the Phi and the R algorithms. Some entries are 
    relevent to one but not the other. An entry that is not relevant to
    Phi or R is ignored in the algorithm. 
*/


#ifndef INCLUDED_HMD_ARG_H
#define INCLUDED_HMD_ARG_H

CPS_END_NAMESPACE
#include<util/vector.h>
CPS_START_NAMESPACE


enum ReunitarizeType { REUNITARIZE_NO    = 0,
		       REUNITARIZE_YES   = 1 };

enum MetropolisType { METROPOLIS_NO    = 0,
		      METROPOLIS_YES   = 1 };



#define MAX_HMD_MASSES  4           /* The maximum number of dynamical 
		 	               masses. */


struct HmdArg {
  /* ??? */

  int n_frm_masses;                 /* The number of dynamical fermion 
				    masses. It must be <= MAX_HMD_MASSES. 
				    To each mass corresponds 
				    Lattice::ExactFlavors() flavors
				    of dynamical fermions. */
				
  Float frm_mass[MAX_HMD_MASSES];   /* The array of masses for each "set" of
				    dynamical fermions. */

  int frm_flavors[MAX_HMD_MASSES];  /* The array of flavors for each "set" of
				    dynamical fermions. This is relevant
				    to the R algorithm only. Bosonic fields
                                    are simulated by setting frm_flavors
                                    to a negative flavor number. */

  int n_bsn_masses;                 /* The number of dynamical boson 
				    masses. It must be <= MAX_HMD_MASSES. 
				    To each mass corresponds 
				    Lattice::ExactFlavors() flavors
				    of dynamical bosons. Relevant to the
                                    Phi algorithm only. */
				
  Float bsn_mass[MAX_HMD_MASSES];   /* The array of masses for each "set" of
				    dynamical bosons. Relevant to the Phi
                                    algorithm only.*/

  int max_num_iter[MAX_HMD_MASSES]; /* The maximum number of iterations
				    of the conjugate gradient to do when 
				    evaluating the effects of the 
				    corresponding "set" of dynamical 
				    fermions/bosons */
                                   

  Float stop_rsd[MAX_HMD_MASSES];   /* The residual for the stopping condition
				    for the conjugate gradient for the 
				    corresponding set of dynamical 
				    fermions/fermions. */


  int steps_per_traj;		    /* The number of steps per trajectory
				    (steps between momentum updates). */

  Float step_size;		    /* The step size. */

  MetropolisType metropolis;        /* Metropolis accep-reject step YES/NO.
                                    This is relevant to the Phi
				    algorithm only. */ 

  ReunitarizeType reunitarize;      /* Reunitarize YES/NO at the end of the
                                    trajectory (for HmcPhi before 
                                    metropolis step) */

};


#endif /* !INCLUDED_HMD_ARG_H */


CPS_END_NAMESPACE
