/*------------------------------------------------------------------*/
/*!\file
  \brief  Definitions of the EigArg structure.
  
  $Id: eig_arg.x,v 1.8 2013-04-05 17:46:30 chulwoo Exp $
*/
/*------------------------------------------------------------------*/

/*! A structure holding the parameters relevant to the eigenvalue measurement.*/
/*! \ingroup algargs */

class EigArg {

    PatternType   pattern_kind;	/*!< Specifies the pattern used
                                to obtain the mass values. For
                                each mass value the eigenvalue is measured.*/
 
    Float Mass_init;	/*!< The initial mass to use  */
    Float Mass_final;	/*!< The final mass to use */
    Float Mass_step;      /*!< The step size in mass */
    
    Float Mass<>;   /*!< Relevant to the ::ARRAY pattern:
				The array of mass values for which 
                                the eigenvalue is measured. */
    int n_masses; 

    int N_eig;		/*!< The number of eigenvectors/values to calculate */
    /* CLAUDIO: added number of accurate eigenvalues to compute */
    int N_eigacc;	/*!< The number of accurate eigenvectors/values to calculate */
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

      // the following two variables should be added, as in QCDSP -- Sam
    string fname<>;   // should be moved from common_arg to here -- Sam
    int ncorr;

    memfun void resize(u_int nmass);
    memfun EigArg();
    memfun ~EigArg();
};


