#include<config.h>
CPS_START_NAMESPACE
/*------------------------------------------------------------------*/
/*!\file
  \brief  Definitions of the EigArg structure.
  
  $Id: eig_arg.x,v 1.2 2004-12-11 20:57:32 chulwoo Exp $
*/
/*------------------------------------------------------------------*/

#ifndef INCLUDED_EIG_ARG_H
#define INCLUDED_EIG_ARG_H           /*!< Prevent multiple inclusion.*/

CPS_END_NAMESPACE
#include <util/vector.h>
CPS_START_NAMESPACE

/*! A structure holding the parameters relevant to the eigenvalue measurement.*/
/*! \ingroup algargs */

struct EigArg {

    Float Mass_init;	/*!< The initial mass to use  */
    Float Mass_final;	/*!< The final mass to use */
    Float Mass_step;      /*!< The step size in mass */
    

    int N_eig;		/*!< The number of eigenvectors/values to calculate */
    int Kalk_Sim;       /*!< Use the Kalkreuter-Simma algorithm. */
    int MaxCG;		/*!<  */
    Float RsdR_a;         /*!< Absolute residual.  */
    Float RsdR_r;         /*!< Relative residual.  */
    Float Rsdlam;         /*!< Residual eigenvalue accuracy.  */
    Float Cv_fact;        /*!< Convergence factor. */
    int N_min;            /*!< Minimum number of Ritz iterations.  */
    int N_max;            /*!< Maximum number of Ritz iterations.  */
    int N_KS_max;         /*!< Maximum Number of Kalkreuter-Simma iterations.
			     */
    int n_renorm;         /*!< How frequently to renormalize in Ritz. */
    int ProjApsiP;        /*!< Whether to reorthogonalize every Ritz iteration.*/

    enum RitzMatType RitzMatOper; /*!< Which operator to determine the
				    eigenvalues/vectors of */

    int print_hsum;       /*!<
			    The eigenvectors are sliced into 3-dim lattice
			    slices perpendicular some direction and the real
			    part of the square norm is calculated on each
			    slice. If this variable is true then these
			    results are printed.
			  */

    int hsum_dir;         /*!< The direction perpendicular to which the
			    eigenvectors are sliced when the sliced square
			    norm is calculated.
			   */

    Float mass;		/*!< The mass to use in the eigenvector solver */
    
};

#endif /* !INCLUDED_EIG_ARG_H */

CPS_END_NAMESPACE
