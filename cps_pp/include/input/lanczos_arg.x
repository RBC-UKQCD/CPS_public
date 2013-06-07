/*------------------------------------------------------------------*/
/*!\file
  \brief  Definitions of the LanczosArg structure.
  
*/
/*------------------------------------------------------------------*/

/*! A structure holding the parameters relevant to the Lanczos eigenvalue measurement.*/
/*! \ingroup algargs */

class  MatrixPolynomialArg {
   int Npol; // degree of the polynomial

   Float params<>;

   Pointer tmp1;
   Pointer tmp2; // the pointer to the temporaly vectors
};

class LanczosArg {
    enum RitzMatType RitzMat_lanczos; /*!< Which operator to be used
                                      in the lanczos process */
    enum RitzMatType RitzMat_convcheck; /*!< Which operator to determine the
				      eigenvalues/vectos of 
                                      in the convergence check */

    Float mass; // mass of the operator

    int nk_lanczos_vectors;  //desired number of Lanczos eigenvectors
    int np_lanczos_vectors;  // extra Lanczos eigenvectors used in Implicit Restart
    Float eigen_shift;     // shift for target eigenvalues, only for DWF, for now 
    Float stop_residual;     /*!< Absolute residual.  */
    int maxiters;            /*!< max number of restartings  */
    int save; /* save eig vecs or not (always save evals) */
    int conv_check; /* do the convergence check every conv_check iters */

    string results<>;  // the file name for ascii output file (number of iteration and whatnot)
    string file<>;  // the file name for eigen vector/values
    MatrixPolynomialArg matpoly_arg; // (will-be-casted) pointer to the MatrixPolynomialArg, used for the Polynomial filtering

 
};


