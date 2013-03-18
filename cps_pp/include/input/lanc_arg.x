/*------------------------------------------------------------------*/
/*!\file
  \brief  Definitions of the LancArg structure.
  
  $Id: lanc_arg.x,v 1.2 2013-03-18 19:33:13 chulwoo Exp $
*/
/*------------------------------------------------------------------*/

/*! A structure holding the parameters relevant to the eigenvalue measurement with IRL.*/
/*! \ingroup algargs */

/* Added by jasper to specify Lanczos eigen solver type */
enum EigenType {
  D         = 0, 
  DDAG      = 1, 
  G5D       = 2, 
  DDAGD     = 3}; 

class LancArg {

    Float mass;	         	       /*!< The mass to use in the eigenvector solver */
    Float stop_rsd;              /*!< Target stopping residual.  */
    Float qr_rsd;                /*!< Targe residual for the QR algorithm.  */
    enum EigenType EigenOper;    /*!< Which operator to determine the eigenvalues/vectors of */
    bool precon;                 /*!< Preconditioned = 1 not = 0 */
    int N_get;		               /*!< The number of eigenvectors/values to get */
    int N_use;		               /*!< The number of eigenvectors/values used in calculation */
    int N_true_get;		           /*!< The number of eigenvectors/values you actually get */
    int ch_ord;                  /*!< Order of Chebyshev polynomial  */
    Float ch_alpha;              /*!< Spectral radius of G5D. Spetrum between alpha and beta are suppressed */
    Float ch_beta;               /*!< Upperbound for the wanted eigenvalue of G5D.  */
    bool ch_sh;                  /*!< Do the internal Chebyshev shift = 1 not = 0.*/
    Float ch_mu;                 /*!< Internal Chebyshev shift size */
    bool lock;                   /*!< Lock the eigenvector using locking transformation = 1 not = 0 */
    int maxits;                  /*!< Maximum iterations for implicit restart*/
 
    string fname<>;              /*!< Output file   */

    //memfun LancArg();
    //memfun ~LancArg();
};


