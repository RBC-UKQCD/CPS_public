#include<config.h>
CPS_START_NAMESPACE
/*  eig_arg.h */

/*  The structure type EigArg holds the parameters specific to
    the subroutine eig(). */

#ifndef INCLUDED_EIG_ARG_H
#define INCLUDED_EIG_ARG_H

CPS_END_NAMESPACE
#include<util/vector.h>
CPS_START_NAMESPACE

struct EigArg {

  Float Mass_init;	// The initial mass to use 
  Float Mass_final;	// The final mass to use
  Float Mass_step;      // The step size in mass

  int N_eig;		// The number of eigenvectors
  int Kalk_Sim;         // Switch to use Jacobi rotations (Kalkreuter-Simma alg.)
  int MaxCG;		// The maximum number of Ritz iterations to do.
  Float RsdR_a;         // Absolute residual
  Float RsdR_r;         // Relative residual
  Float Rsdlam;         // Residual eigenvalue accuracy
  Float Cv_fact;        // Convergence factor 
  int N_min;            // Minimum number of ritz iterations
  int N_max;            // Maximum number of ritz iterations
  int N_KS_max;         // Maximum Number of KS iterations
  int n_renorm;         // Renormalize this often in Ritz
  int ProjApsiP;        // Flag whether to reorthonalize every Ritz iter.

  enum RitzMatType RitzMatOper; // Which operator to determine eigenvalues of

  int print_hsum;       // Flag to indicate whether to print 1-D eigenvec.
  int hsum_dir;         // propagation direction
                                // 0->x, 1->y, 2->z, 3->t

  Float mass;		// The mass to use in the eigenvector solver
};

#endif /* !INCLUDED_EIG_ARG_H */
CPS_END_NAMESPACE
