#include<config.h>
/*!\file
  \brief  Definition of the Dirac operator classes: DiracOp, DiracOpStagTypes.

  $Id: dirac_op.h,v 1.40 2013-06-07 19:26:34 chulwoo Exp $
*/

#ifndef INCLUDED_DIRAC_OP_H
#define INCLUDED_DIRAC_OP_H

#include <util/enum.h>
#include <util/lattice.h>
#include <util/vector.h>
#include <alg/cg_arg.h>
#include <alg/eig_arg.h>
#include <stdint.h>
#include <alg/lanczos_arg.h>

//#include <alg/matrixpolynomial_arg.h>
#include <util/error.h>

#ifdef USE_QUDA
#include <alg/quda_arg.h>
#endif

CPS_START_NAMESPACE

typedef unsigned long long IntFlopCounter;


//------------------------------------------------------------------
//! A class representing operations on the Dirac operator.
/*!
  The constructor changes the storage order of the gauge field to an order
  appropriate to the type of fermion action. This change persists
  throughout the lifetime of the object.

  Only one instance of this class is allowed to be in existence at any time.
  
  This is an abstract base class, so the details specific to the various
  types of fermion action are defined in the derived classes.
 */
//------------------------------------------------------------------
class DiracOp
{
 private:
  const char *cname;                     // Class name.


 protected:
  Lattice& lat;                    //!< The lattice..
  Vector *f_out;                   //!< Pointer aliasing vector \a f_field_out.
  Vector *f_in;                    //!< Pointer aliasing vector \a f_field_out.
  CgArg *dirac_arg;                //!< Solver parameters.
  CnvFrmType cnv_frm;              //!< Field storage order conversion flag.
  Matrix *gauge_field;             //!< Pointer to the gauge field.
  
  //! Global sum of a floating point number.
  virtual void DiracOpGlbSum(Float *float_p);					
     // The global sum used by InvCg. At the base class
     // level it is just the usual global sum. It is overloaded
     // by any definitions by the derived classes. This is needed
     // for example by DiracOpDwf where the global sum has different
     // meaning for s_nodes = 1 and s_nodes > 1.


 public:
  static int scope_lock;
  //!< lock that prevents more than one DiracOp object to be on scope at any time.

  // PAB... uniform performance counting
  static IntFlopCounter CGflops; //!< Flops counter.
  static unsigned int CGcount;

  DiracOp(Lattice& latt,           // Lattice object.
	  Vector *f_field_out,     // Output fermion field ptr.
	  Vector *f_field_in,      // Input fermion field ptr.
	  CgArg *arg,              // Argument structure
	  CnvFrmType convert); // Fermion conversion flag

  virtual ~DiracOp();

  //eigen CG //by Qi Liu
  int InvEigCg(Vector *sol, Vector *src, Float *true_res, const int nev, const int m, Vector **V, const int vec_len, Float *M, float **U, Rcomplex *invH, const int def_len, const Float *restart, const int restart_len);


  // Inverse matrix from eigen mode, for Low Mode Average
  int InvLowModeApprox(Vector* out,
		       Vector* in, 
		       char* fname_eigenv, 
		       int neig,
		       Float* true_res );

  // lowmode projection onto src
  void InvLowModeProj( Vector* in, 
		       char* fname_eigenv, 
		       int neig );

  //! The matrix inversion used in the molecular dynamics algorithms.
  int InvCgShift(Vector *out,
	    Vector *in,
	    Float src_norm_sq,
	    Float *true_res,
           Float *shift	);

  int InvCg(Vector *out,
	    Vector *in,
	    Float src_norm_sq,
	    Float *true_res);


#if 0
	{
		InvCgShift(out,in,src_norm_sq,true_res,NULL);
	    }
#endif
     // The Conjugate Gradient inverter. 
     // Source is *in, initial guess and solution is *out.
     // src_norm_sq is the norm of the source squared 
     // (dot product of source * source). If src_norm_sq=0 
     // it is calculated inside InvCg. 
     // If true_res !=0 the value of the true residual is returned
     // in true_res.
     // *true_res = |src - MatPcDagMatPc * sol| / |src|

  //! The matrix inversion used in the molecular dynamics algorithms.  
  int InvCg(Vector *out, Vector *in, Float src_norm_sq);
      // Same as original but with true_res=0

  //! The matrix inversion used in the molecular dynamics algorithms.
  virtual int InvCg(Vector *out, Vector *in, Float *true_res);
     // Same as original but with src_norm_sq=0.0

  //! The matrix inversion used in the molecular dynamics algorithms.
  int InvCg(Vector *out, Vector *in);
  // Same as original but with s
  
  //! The matrix inversion used in the molecular dynamics algorithms.
  int InvCg(Float src_norm_sq, Float *true_res);
     // Same as original but with 
     // in=f_in, out=f_out

  //! The matrix inversion used in the molecular dynamics algorithms.
  int InvCg(Float src_norm_sq);
     // Same as original but with 
     // in=f_in, out=f_out, true_res=0

  //! The matrix inversion used in the molecular dynamics algorithms.
  virtual int InvCg(Float *true_res);
     // Same as original but with 
     // in=f_in, out=f_out, src_norm_sq=0.0

  //! The matrix inversion used in the molecular dynamics algorithms.
  int InvCg(void);
  // Same as original but with 
  // in=f_in, out=f_out, src_norm_sq=0.0, true_res=0
  
  //! Ritz eigensolver.
  int Ritz(Vector **psi_all, int N_eig, Float &lambda,
	   Float RsdR_a, Float RsdR_r, Float Rsdlam, Float Cutl_zero,
	   int n_renorm, int Kalk_Sim, int N_min, int N_max, Float Cv_fact,
	   int MaxCG, int ProjApsiP);
  // Ritz minimizer for eigenvectors.
  // The matrix that is used is RitzMat.
  // RitzMat must be hermitian and positive def.

  
  //! Implicitly restarted Lanczos (TB Dec. 2010)
  //  with Polynomial accerelation+shift+ many (TIZB Aug. 2011)
  int ImpResLanczos( Vector **V, //Lanczos vectors, eigenvectors of RitzMat on return
		     Float *alpha, // eigenvalues
		     LanczosArg* lanczos_arg
		     ); // arguments for eigen solver algorithm
  //
  // In this routine two different RitzMat will be used :
  // 1.  lanczos_arg->RitzMat_lanczos, which is used
  //     in the kernel of lanczos process,
  //      e.g. MATPCDAG_MATPC
  // 2.  lanczos_arg->RitzMat_convcheck, which is used for check of convergence,
  //     e.g. MATPC_HERM
  //
  // two matrixces has to be commutable, so that they shares
  // the same set of eigenvectors.
  //
  // RitzMat_convcheck will be used for the eigenvalue computation
  //  i.e. if RitzMat_convcheck==MATPC_HERM, then the eigenvalue file will be
  //  MATPC_HERM's eigenvalue not the squared one.
  //
  // Use RitzMat_convcheck=MATPC_HERM to save computation (by factor of ~two).
  // if MATPC_HERM is not implemented yet for your femion
  // you could either implmented (not hard) 
  // or just uncomment the following to use the
  // squared matrix (easiest), MATPCDAG_MATPC
  //

  //! the fundamental lanczos steps (TB Dec. 2010, TIZB Aug. 2011)
  void lanczos(int k, int m, int f_size,
	       Vector *Apsi, Vector *r, Vector **V, Float *alpha, Float* beta,
	       MatrixPolynomialArg* cheby_arg,
	       RitzMatType RitMat_lanczos );


  
  //! Vector orthogonalisation
  void GramSchm(Vector **psi, int Npsi, Vector **vec, int Nvec, int f_size);
  void GramSchm(Vector *psi, Vector **vec, int Nvec, int f_size);

  //! Jacobi diagonalisation of a matrix.
  int Jacobi(Vector **psi, int N_eig, Float *lambda, 
	     Complex *off_diag, 
	     Float Toler, int N_max);
  // Jacobi diagonalizer of a single site matrix.

  //! Multishift CG invertor used in RHMC.
  //
  virtual int MInvCG(Vector **out, Vector *in, Float in_norm, Float *shift, 
	     int Nshift, int isz, Float *RsdCG, 
	     MultiShiftSolveType type, Float *alpha);

  //! Chronological initial guess for the solver. 
  void MinResExt(Vector *psi, Vector *phi, Vector **psi_old, 
		 Vector **vm, int degree);

// Pure virtual functions
//------------------------------------------------------------------

//! Resets the parameters used in the solver (and elsewhere).
/*!
  The quark mass is set from the set parameters used in the solver
  and made available to all methods in this class.
  \param arg The new set of solver parameters.
 */
  virtual void DiracArg(CgArg *arg) = 0;
     // It sets the dirac_arg pointer to arg and initializes
     // the relevant parameters (kappa, m^2, ...).
  virtual CgArg * const DiracArg() {return dirac_arg;}

  //! Multiplication by the square of the odd-even preconditioned fermion matrix.
  /*!
    A vector, defined on a lattice sites of a single parity, is multiplied by
    \f$ M^\dagger M \f$ where \e M is a single parity of the odd-even
    preconditioned fermion matrix.
    \param out The result
    \param in The vector to be multiplied.
    \param dot_prd Whether or not to compute the real part of the
    local dot product (\e i.e. with no global sum)
    of the vectors \a out and \a in. If \a dot_prod is initially
    non-zero then this is computed and the result placed here.
  */
  virtual void MatPcDagMatPc(Vector *out, Vector *in, Float *dot_prd=0) = 0;
 
  //! The derivative part of the fermion matrix.
  /*!
    A vector defined on lattice sites of a single parity is multiplied by
    the derivative part of the fermion matrix.
    \param out The resulting vector.
    \param in The vector to be multiplied.
    \param cb The parity on which the vector \a in is defined.
    \param dag Whether to multiply by the hermitian conjugate matrix:
    Should be set to 1 to multiply by the hermitian conjugate, 0 otherwise.
   */
  virtual void Dslash(Vector *out, 
		      Vector *in,
		      ChkbType cb, 
		      DagType dag)=0;

  //! Fermion matrix inversion.
  /*!
    Solves <em> A x = b </em> for \e x, where \a A is the
    fermion matrix. The vectors are defined on the whole lattice,
    not just on sites of a single parity.
    \param out The initial guess of the solution vector \e x.
    \param in The source vector \e b.
    \param true_res Whether or not to report the true residual. The true
    residual will be written here if this is non-zero.
    \param prs_in Whether or not the source vector is allowed to be
    overwritten, thereby saving memory. For staggered fermions it is
    preserved regardless of the value of \a prs_in. 
    \return The number of solver iterations.
    \post \a out contains the solution vector \e x.
    \post \a true_res contains the true residual, if this was non-zero to start with.
    The residual is  <em>  |b - A x| / |b| </em>
  */
  virtual int MatInv(Vector *out, 
		     Vector *in, 
		     Float *true_res, 
		     PreserveType prs_in = PRESERVE_YES) = 0;

  //! Fermion matrix inversion.
  /*!
    Solves <em> A x = b </em> for \e x, where \e A is the
    fermion matrix. The vectors are defined on the whole lattice,
    not just on sites of a single parity.
    \param out The initial guess of the solution vector \e x.
    \param in The source vector \e b.
    \param prs_in Whether or not the source vector is allowed to be
    overwritten, thereby saving memory. For staggered fermions it is
    preserved regardless of the value of \a prs_in. 
    \return The number of solver iterations.
    \post \a out contains the solution vector \e x.
  */

  virtual int eig_MatInv(Vector **V, const int vec_len, Float *M, const int nev, const int m, float **U, Rcomplex *invH, const int def_len, const Float *restart, const int restart_len, Vector *out, Vector *in, Float *true_res, PreserveType prs_in = PRESERVE_YES)
  {
	  const char *fname="eig_MatInv()";
	  ERR.General(cname, fname,"only implemented for dwf currently, so not pure virtual\n");
	  return 0;
  } 
  //eigCG inversion //by Qi Liu

  virtual int MatInv(Vector *out, Vector *in, 
		     PreserveType prs_in = PRESERVE_YES) = 0;
     // Same as original but true_res=0.

  //! Fermion matrix inversion.
  /*!
    Solves <em> A x = x </em> for \a f_out,
    where \a A is the fermion matrix.
    The vectors \a f_field_in and \a f_field_out, defined in the constructor,
    are used for \e b and \e x respectively

    \pre \a f_field_in and \a f_field_out, defined in the
    constructor, must point to spin-colour vectors defined on the
    whole lattice, not just on sites of a single parity.
    \pre \a f_field_out contains the initial guess for \e x.
    
    \param true_res Whether or not to report the true residual. The true
    residual will be written here if this is non-zero.
    \param prs_in Whether or not the source vector is allowed to be
    overwritten, thereby saving memory. For staggered fermions it is
    preserved regardless of the value of \a prs_in. 
    \return The number of solver iterations.
    \post \a f_field_out contains the solution vector \e x.
    \post \a true_res contains the true residual, if this was non-zero to start with.
    The residual is  <em>  |b - A x| / |b| </em>
  */
  virtual int MatInv(Float *true_res, 
		     PreserveType prs_in = PRESERVE_YES) = 0;
     // Same as original but in = f_in and out = f_out.

  //! Fermion matrix inversion.
  /*!
    Solves <em> A x = b </em> for \a x,
    where \a A is the fermion matrix.
    The vectors \a f_field_in and \a f_field_out, defined in the constructor,
    are used for \e b and \e x respectively
    
    \pre \a f_field_in and \a f_field_out, defined in the
    constructor, must point to spin-colour vectors defined on the
    whole lattice, not just on sites of a single parity.
    \pre \a f_field_out contains the initial guess for \e x.

    \param prs_in Whether or not the source vector is allowed to be
    overwritten, thereby saving memory. For staggered fermions it is
    preserved regardless of the value of \a prs_in. 
    \return The number of solver iterations.
    \post \a f_field_out contains the solution vector \e x.
  */
  virtual int MatInv(PreserveType prs_in = PRESERVE_YES) = 0;
     // Same as original but in = f_in, out = f_out, true_res=0.


//! Computes eigenvectors and eigenvalues.
/*!
  Uses a Conjugate Gradient based Ritz functional minimisation method
  to compute the \e n lowest eigenvectors and eigenvalues of some fermion
  matrix operator.
 */
   virtual int RitzEig(Vector **eigenv, Float lambda[], int valid_eig[], 
		      EigArg *eig_arg) = 0;
     // The eigenvector solver using Ritz minimization and Jacobi reconstruction
     // (using the Kalkreuter-Simma algorithm) of the eigenvectors of RitzEigMat.
     // The matrix that is used is RitzEigMat.
     // RitzEigMat is not necessarily to be positive def. (but must be hermitian).
     // For RitzEigJac, RitzMat = RitzEigMat*RitzEigMat when using KS algorithm.
     // Can also use non KS algorithm to return eigenvalues only of 
     // RitzEigMat = RitzMat.

   //! Multiplication by a fermion matrix operator.
   /*!
     Multiplies a vector by the hermitian matrix defined for the
     eigenvalue calculations. Exactly what this is and whether the operation
     acts on vectors defined on the full lattice or just one parity depends
     on the fermion action used.
     \param out The resulting vector.
     \param in The vector to be multiplied.
   */
  virtual void RitzEigMat(Vector *out, Vector *in) = 0;

   //! Multiplication by a fermion matrix operator.
   /*!
     Multiplies a vector by something like
     the hermitian matrix defined (in the structure
     \a arg used in the constructor) for the
     eigenvalue calculations. Exactly what this is and whether the operation
     acts on vectors defined on the full lattice or just one parity depends
     on the fermion action used.
     \param out The resulting vector.
     \param in The vector to be multiplied.
   */
  virtual void RitzMat(Vector *out, Vector *in) = 0;

  //TIZB  PolynomialAccerelation
  virtual void RitzMat(Vector *out, Vector *in,
		       MatrixPolynomialArg* cheby_arg);

//! The size of a fermion field.
/*!
  \return The size on each node, in terms of floating point numbers,
  of the fermion field on which the matrix defined (in the structure
  \a arg used in the constructor) for the eigenvalue calculations operates.
*/
  virtual int RitzLatSize() = 0;

};


//------------------------------------------------------------------
//
// DiracOpStagTypes class.
//
//! A class describing the Dirac operator for all sorts of staggered fermions.
/*!
  This is an abstract base class from which the staggered fermion Dirac
  operator classes are derived.
*/
//------------------------------------------------------------------
class DiracOpStagTypes : public DiracOp
{
 private:
  const char *cname;    // Class name.

 public:
  DiracOpStagTypes(Lattice& latt,            // Lattice object.
		   Vector *f_field_out,      // Output fermion field ptr.
		   Vector *f_field_in,       // Input fermion field ptr.
		   CgArg *arg,               // Argument structure
		   CnvFrmType convert);  // Fermion conversion flag

  virtual ~DiracOpStagTypes();

  virtual int RitzEig(Vector **eigenv, Float lambda[], int valid_eig[], 
	              EigArg *eig_arg);
     // The eigenvector solver using Ritz minimization and Jacobi reconstruction
     // (using the Kalkreuter-Simma algorithm) of the eigenvectors of RitzEigMat.
     // The matrix that is used is RitzEigMat.
     // RitzEigMat is not necessarily to be positive def. (but must be hermitian).
     // For RitzEigJac, RitzMat = RitzEigMat*RitzEigMat when using KS algorithm.
     // Can also use non KS algorithm to return eigenvalues only of 
     // RitzEigMat = RitzMat.

  virtual void RitzMat(Vector *out, Vector *in)
  {ERR.NotImplemented(cname,"RitzMat(V*,V*)");}
     // RitzMat is the fermion matrix used in Ritz
     // RitzMat works on the full lattice or half lattice
     // The in, out fields are defined on the full or half lattice.

//  virtual void RitzMat(Vector *out, Vector *in,
 //                      MatrixPolynomialArg* cheby_arg);
    // TIZB


  virtual int RitzLatSize();

};

//------------------------------------------------------------------
//
// DiracOpStag is derived from DiracOpStagTypes and is the front
// end for all Dirac operators associated with Staggered fermions.

//! A class describing the Dirac operator for staggered fermions.
/*!
  The staggered fermion action is
  \f[
  M_{xy} =
  m_0 - \sum_mu e^{\sum_{i=0}^{\mu-1}x_i} 
  [ U^\dagger_\mu(x) \delta_{x\,y+\mu} - U_\mu(x-\mu) \delta_{x\,y-\mu} ]
  \f]

  \e N.B  The phases are implemented in the gauge field when it is
  converted to staggered (::STAG) storage order.

  The constructor changes the storage order of the gauge field to ::STAG order
  This change persists throughout the lifetime of the object.

  Only one instance of this class is allowed to be in existence at any time.
 */
//------------------------------------------------------------------
class DiracOpStag : public DiracOpStagTypes
{
 private:
  const char *cname;         // Class name.

  Float mass_rs;       // rescaled mass
  Float mass_sq;       // = mass^2

  int f_size_cb;       //The node checkerbrd. size of the ferm. field

  Vector *frm_tmp;     // Temporary fermion field

 public:
  DiracOpStag(Lattice& latt,            // Lattice object.
	      Vector *f_field_out,      // Output fermion field ptr.
	      Vector *f_field_in,       // Input fermion field ptr.
	      CgArg *arg,               // Argument structure
	      CnvFrmType convert);  // Fermion conversion flag

  virtual ~DiracOpStag();

  void DiracArg(CgArg *arg);
     // It sets the dirac_arg pointer to arg and initializes
     // mass_sq = 4 * mass^2

  void MatPcDagMatPc(Vector *out, Vector *in, Float *dot_prd=0);
     // MatPcDagMatPc is the fermion matrix that appears in the HMC 
     // evolution. It is a Hermitian matrix where M is
     // the Dirac Operator matrix.        
     // MatPcDagMatPc connects only even-->even or odd-->odd sites.
     // The in, out fields are defined on a checkerboard.
     // If dot_prd is not 0 then the dot product (on node)
     // <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.

  void Dslash(Vector *out, 
	      Vector *in,
	      ChkbType cb, 
	      DagType dag);
     // Dslash is the derivative part of the fermion matrix. 
     // Dslash conects only odd-->even or even-->odd sites.
     // The in, out fields are defined on a checkerboard.
     // cb refers to the checkerboard of the in field.

  //! The derivative part of the fermion matrix.
  void Dslash(Vector *out, 
	      Vector *in,
	      ChkbType cb, 
	      DagType dag,
	      int dir_flag);
     // Dslash is the derivative part of the fermion matrix. 
     // Dslash conects only odd-->even or even-->odd sites.
     // The in, out fields are defined on a checkerboard.
     // cb refers to the checkerboard of the in field.
     // dir_flag is flag which takes value 0 when all direction contribute 
     // 1 - when only the special anisotropic direction contributes to D,
     // 2 - when all  except the special anisotropic direction.

  int MatInv(Vector *out, 
	     Vector *in, 
	     Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
     // The inverse of the Dirac Operator (D+m)
     // using Conjugate gradient.
     // Assume: the vector in contains both even and odd src.
     //		even part is the 1st part.
     // Return: the vector out contains both even and odd solutions.
     // 	the even solution is the 1st part.
     // If true_res !=0 the value of the true residual is returned
     // in true_res.
     // *true_res = |src - MatPcDagMatPc * sol| / |src|
     // prs_in is not used. The source in is always preserved.
     // The function returns the total number of CG iterations.

  int MatInv(Vector *out, 
	     Vector *in,
	     PreserveType prs_in = PRESERVE_YES);
     // Same as original but true_res=0.

  int MatInv(Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
     // Same as original but in = f_in and out = f_out.

  int MatInv(PreserveType prs_in = PRESERVE_YES);
     // Same as original but in = f_in, out = f_out, true_res=0.

  void RitzEigMat(Vector *out, Vector *in);
     // RitzEigMat is the fermion matrix used in RitzEig
     // RitzEigMat works on the full lattice or half lattice
     // The in, out fields are defined on the full or half lattice.

  //void RitzMat(Vector *out, Vector *in);
     // RitzMat is the fermion matrix used in Ritz
     // RitzMat works on the full lattice or half lattice
     // The in, out fields are defined on the full or half lattice.

  //virtual void RitzMat(Vector *out, Vector *in,
	//	       MatrixPolynomialArg* cheby_arg);
    // TIZB

};

//------------------------------------------------------------------
//
// DiracOpAsqtad is derived from DiracOpStagTypes and is the front
// end for all Dirac operators associated with Staggered fermions.

//! A class describing the Dirac operator for Asqtad improved staggered fermions.
/*!
  The Asqtad fermion action is complicated but is described in the User
  Guide.
  
  The constructor changes the storage order of the gauge field to ::STAG order
  This change persists throughout the lifetime of the object.

  Only one instance of this class is allowed to be in existence at any time.
 */

//------------------------------------------------------------------
class DiracOpAsqtad : public DiracOpStagTypes
{
 private:
  const char *cname;         // Class name.

  Float mass_rs;       // rescaled mass
  Float mass_sq;       // = mass^2

  int f_size_cb;       //The node checkerbrd. size of the ferm. field

  Vector *frm_tmp;     // Temporary fermion field


 public:
  DiracOpAsqtad(Lattice& latt,            // Lattice object.
	      Vector *f_field_out,      // Output fermion field ptr.
	      Vector *f_field_in,       // Input fermion field ptr.
	      CgArg *arg,               // Argument structure
	      CnvFrmType convert);  // Fermion conversion flag

  virtual ~DiracOpAsqtad();

  void DiracArg(CgArg *arg);
     // It sets the dirac_arg pointer to arg and initializes
     // mass_sq = 4 * mass^2

  void MatPcDagMatPc(Vector *out, Vector *in, Float *dot_prd=0);
     // MatPcDagMatPc is the fermion matrix that appears in the HMC 
     // evolution. It is a Hermitian matrix where M is
     // the Dirac Operator matrix.        
     // MatPcDagMatPc connects only even-->even or odd-->odd sites.
     // The in, out fields are defined on a checkerboard.
     // If dot_prd is not 0 then the dot product (on node)
     // <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.

  void Dslash(Vector *out, 
	      Vector *in,
	      ChkbType cb, 
	      DagType dag);
     // Dslash is the derivative part of the fermion matrix. 
     // Dslash conects only odd-->even or even-->odd sites.
     // The in, out fields are defined on a checkerboard.
     // cb refers to the checkerboard of the in field.

  //! The derivative part of the fermion matrix.
  void Dslash(Vector *out, 
	      Vector *in,
	      ChkbType cb, 
	      DagType dag,
	      int dir_flag);
  
     // Dslash is the derivative part of the fermion matrix. 
     // Dslash conects only odd-->even or even-->odd sites.
     // The in, out fields are defined on a checkerboard.
     // cb refers to the checkerboard of the in field.
     // dir_flag is flag which takes value 0 when all direction contribute 
     // 1 - when only the special anisotropic direction contributes to D,
     // 2 - when all  except the special anisotropic direction.

  int MatInv(Vector *out, 
	     Vector *in, 
	     Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
     // The inverse of the Dirac Operator (D+m)
     // using Conjugate gradient.
     // Assume: the vector in contains both even and odd src.
     //		even part is the 1st part.
     // Return: the vector out contains both even and odd solutions.
     // 	the even solution is the 1st part.
     // If true_res !=0 the value of the true residual is returned
     // in true_res.
     // *true_res = |src - MatPcDagMatPc * sol| / |src|
     // prs_in is not used. The source in is always preserved.
     // The function returns the total number of CG iterations.

  int MatInv(Vector *out, 
	     Vector *in,
	     PreserveType prs_in = PRESERVE_YES);
     // Same as original but true_res=0.

  int MatInv(Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
     // Same as original but in = f_in and out = f_out.

  int MatInv(PreserveType prs_in = PRESERVE_YES);
     // Same as original but in = f_in, out = f_out, true_res=0.

  void RitzEigMat(Vector *out, Vector *in);
     // RitzEigMat is the fermion matrix used in RitzEig
     // RitzEigMat works on the full lattice or half lattice
     // The in, out fields are defined on the full or half lattice.

  //void RitzMat(Vector *out, Vector *in);
     // RitzMat is the fermion matrix used in Ritz
     // RitzMat works on the full lattice or half lattice
     // The in, out fields are defined on the full or half lattice.

  //virtual void RitzMat(Vector *out, Vector *in,
//		       MatrixPolynomialArg* cheby_arg)
 // {ERR.NotImplemented(cname,"RitzMat(V*,V*,MatrixPolynomialArg*");}
  // TIZB

};


//------------------------------------------------------------------
//
// DiracOpP4 is derived from DiracOpStagTypes and is the front
// end for all Dirac operators associated with Staggered fermions.

//! A class describing the Dirac operator for P4 staggered fermions.
/*!
  The P4 action includes the naive one-link terms in the fermion actions,
  but also incorporates a "knight's move" three-link term.

  The motivation is to improve the rotational invariance of the Dirac
  propogator.  The p4 action seems to do the best job of reducing
  leading order lattice artifacts in the dispersion relation.  The relative
  coefficient of the naive staggered term and the "knight's move" term have 
  also been calculated to one-loop order in perturbation theory.

  Only one instance of this class is allowed to be in existence at any time.
 */

//------------------------------------------------------------------
class DiracOpP4 : public DiracOpStagTypes
{
 private:
  const char *cname;         // Class name.

  Float mass_rs;       // rescaled mass
  Float mass_sq;       // = mass^2

  int f_size_cb;       //The node checkerbrd. size of the ferm. field

  Vector *frm_tmp;     // Temporary fermion field


 public:
  DiracOpP4(Lattice& latt,            // Lattice object.
	      Vector *f_field_out,      // Output fermion field ptr.
	      Vector *f_field_in,       // Input fermion field ptr.
	      CgArg *arg,               // Argument structure
	      CnvFrmType convert);  // Fermion conversion flag

  virtual ~DiracOpP4();

  void DiracArg(CgArg *arg);
     // It sets the dirac_arg pointer to arg and initializes
     // mass_sq = 4 * mass^2

  void MatPcDagMatPc(Vector *out, Vector *in, Float *dot_prd=0);
     // MatPcDagMatPc is the fermion matrix that appears in the HMC 
     // evolution. It is a Hermitian matrix where M is
     // the Dirac Operator matrix.        
     // MatPcDagMatPc connects only even-->even or odd-->odd sites.
     // The in, out fields are defined on a checkerboard.
     // If dot_prd is not 0 then the dot product (on node)
     // <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.

  void Dslash(Vector *out, 
	      Vector *in,
	      ChkbType cb, 
	      DagType dag);
     // Dslash is the derivative part of the fermion matrix. 
     // Dslash conects only odd-->even or even-->odd sites.
     // The in, out fields are defined on a checkerboard.
     // cb refers to the checkerboard of the in field.

  //! The derivative part of the fermion matrix.
  /* void Dslash(Vector *out, 
	      Vector *in,
	      ChkbType cb, 
	      DagType dag,
	      int dir_flag);*/
     // Dslash is the derivative part of the fermion matrix. 
     // Dslash conects only odd-->even or even-->odd sites.
     // The in, out fields are defined on a checkerboard.
     // cb refers to the checkerboard of the in field.
     // dir_flag is flag which takes value 0 when all direction contribute 
     // 1 - when only the special anisotropic direction contributes to D,
     // 2 - when all  except the special anisotropic direction.

  void dMdmu(Vector *out, 
	      Vector *in,
	      ChkbType cb, 
	      DagType dag,
	      int order);
     // dMDmu is the derivative of the fermion matrix with
     // respect to the chemical potential. 
     // dMdmu conects only odd-->even or even-->odd sites.
     // The in, out fields are defined on a checkerboard.
     // cb refers to the checkerboard of the in field.
     // order refers to the order of the derivative.

  int MatInv(Vector *out, 
	     Vector *in, 
	     Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
     // The inverse of the Dirac Operator (D+m)
     // using Conjugate gradient.
     // Assume: the vector in contains both even and odd src.
     //		even part is the 1st part.
     // Return: the vector out contains both even and odd solutions.
     // 	the even solution is the 1st part.
     // If true_res !=0 the value of the true residual is returned
     // in true_res.
     // *true_res = |src - MatPcDagMatPc * sol| / |src|
     // prs_in is not used. The source in is always preserved.
     // The function returns the total number of CG iterations.

  int MatInv(Vector *out, 
	     Vector *in,
	     PreserveType prs_in = PRESERVE_YES);
     // Same as original but true_res=0.

  int MatInv(Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
     // Same as original but in = f_in and out = f_out.

  int MatInv(PreserveType prs_in = PRESERVE_YES);
     // Same as original but in = f_in, out = f_out, true_res=0.

  void RitzEigMat(Vector *out, Vector *in);
     // RitzEigMat is the fermion matrix used in RitzEig
     // RitzEigMat works on the full lattice or half lattice
     // The in, out fields are defined on the full or half lattice.

  //void RitzMat(Vector *out, Vector *in);
     // RitzMat is the fermion matrix used in Ritz
     // RitzMat works on the full lattice or half lattice
     // The in, out fields are defined on the full or half lattice.

  //v////oid RitzMat(Vector *out, Vector *in,
//		       MatrixPolynomialArg* cheby_arg)
 // {ERR.NotImplemented(cname,"RitzMat(V*,V*,MatrixPolynomialArg*");}
  // TIZB


};


//------------------------------------------------------------------
//
// DiracOpWilsonTypes class.

//! A class describing the Dirac operator for all sorts of Wilson fermions.
/*! This is an abstract base class from which the Wilson Dirac operator
  classes are derived.

  The general Wilson type of fermion matrix is \f$ M = A - \kappa D \f$
  where the Wilson hopping matrix \e D is
  \f[
  D_{xy} = \sum_{\mu}(1-\gamma_{\mu})U_\mu(x)\delta_{y\,x+\mu}
  +(1+\gamma_{\mu})U_\mu^{\dagger}(x-\mu)\delta_{y\,x-\mu}
  \f]w

  Note that \e D connects sites of opposite parity.
  
  For Wilson fermions \e A is just a unit matrix.

  For clover improved fermions   \e A is the clover matrix
  \f[
  A(x) = 1 - \frac{1}{2} \kappa c_\mathrm{sw}   \sum_{\mu<\nu}
  [\gamma_\mu, \gamma_\nu] F_{\mu\nu}(x) 
  \f]
  where 
  the field strength tensor \f$ F_{\mu\nu} \f$ is the matrix
  \f[
  F_{\mu\nu}(x) = \frac{1}{4}\frac{1}{2}
                (P_{\mu\nu}(x) - P^\dagger_{\mu\nu}(x))
  \f]
  in which
\f[
	P_{\mu\nu}(x)  =
    U_\mu(x) U_\nu(x+mu) U_\mu(x+\nu)^\dagger U_\nu(x)^\dagger +  
    U_\nu(x) U_\mu(x-\mu+\nu)^\dagger U_\nu(x-\mu)^\dagger U_\mu(x-\mu) +
 \f]\f[    
    U_\mu(x-\mu)^\dagger U_\nu(x-\mu-\nu)^\dagger U_\mu(x-\mu-\nu) U_\nu(x-\nu) + 
    U_\nu(x-\nu)^\dagger U_\mu(x-\nu) U_\nu(x+\mu-\nu) U_\mu(x)^\dagger
\f]
  
  The odd-even preconditioned fermion matrix is
  \f$ A - \kappa^2 D A^{-1} D \f$  on odd parity lattice sites
  and \e A on even parity sites.

  Note that for Wilson fermions the even sites are completely decoupled.  
  
The following representation of the Euclidean gamma matrices is probably used:
\f[
\gamma_t = \pmatrix{
0&0&1&0 \cr 0&0&0&1 \cr1&0&0&0 \cr 0&1&0&0 \cr 
}

\gamma_x = \pmatrix{
0&0&0&i \cr 0&0&i&0 \cr 0&-i&0&0 \cr -i&0&0&0 \cr
}
\f]
\f[
\gamma_y = \pmatrix{
0&0&0&-1 \cr 0&0&1&0 \cr 0&1&0&0 \cr -1&0&0&0 \cr
}

\gamma_z = \pmatrix{
0&0&i&0 \cr 0&0&0&-i \cr -i&0&0&0 \cr 0&i&0&0 \cr
}
\f]

  The constructor changes the storage order of the gauge field to ::WILSON
  order This change persists throughout the lifetime of the object.

  Only one instance of this class is allowed to be in existence at any time.
*/
//------------------------------------------------------------------
class DiracOpWilsonTypes : public DiracOp
{
 private:
  const char *cname;    // Class name.

  protected:

  //! Multiply a vector by  gamma_5.
  void MultGamma(Vector *out, const Vector *in, int gamma_num, int nodevol);

 public:
  DiracOpWilsonTypes(Lattice& latt,            // Lattice object.
		     Vector *f_field_out,      // Output fermion field ptr.
		     Vector *f_field_in,       // Input fermion field ptr.
		     CgArg *arg,               // Argument structure
		     CnvFrmType convert);  // Fermion conversion flag

  virtual ~DiracOpWilsonTypes();

  
  virtual void MatHerm(Vector *out, Vector *in) = 0;
  //!< Multiplication by the hermitian fermion matrix.
/*!<
  The fermion matrix is left multiplied by gamma_5 to make it hermitian.
  The vectors are defined on the whole lattice.
  \param out The resulting vector.
  \param in The vector to be multiplied.
*/


  virtual void Mat(Vector *out, Vector *in) = 0;
  //!< Multiplication by the full fermion matrix.
/*!<
  The vectors are defined on the whole lattice.
  \param out The resulting vector.
  \param in The vector to be multiplied.
*/

  virtual void MatDag(Vector *out, Vector *in) = 0;
  //!< Multiplication by the hermitian conjugate of the full fermion matrix.  
/*!<
  The vectors are defined on the whole lattice.
  \param out The resulting vector.
  \param in The vector to be multiplied.
*/

  //! Multiplication by the square of the fermion matrix.
  virtual void MatDagMat(Vector *out, Vector *in);

  virtual void MatPc(Vector *out, Vector *in) = 0;

  virtual int RitzEig(Vector **eigenv, Float lambda[], int valid_eig[], 
		      EigArg *eig_arg);
     // The eigenvector solver using Ritz minimization and Jacobi reconstruction
     // (using the Kalkreuter-Simma algorithm) of the eigenvectors of RitzEigMat.
     // The matrix that is used is RitzEigMat.
     // RitzEigMat is not necessarily to be positive def. (but must be hermitian).
     // For RitzEigJac, RitzMat = RitzEigMat*RitzEigMat when using KS algorithm.
     // Can also use non KS algorithm to return eigenvalues only of 
     // RitzEigMat = RitzMat.

  virtual void RitzEigMat(Vector *out, Vector *in);
     // RitzEigMat is the fermion matrix used in RitzEig
     // RitzEigMat works on the full lattice or half lattice
     // The in, out fields are defined on the full or half lattice.

  virtual void RitzMat(Vector *out, Vector *in);
//  {ERR.NotImplemented(cname,"RitzMat(V*,V*)");}
     // RitzMat is the fermion matrix used in Ritz
     // RitzMat works on the full lattice or half lattice
     // The in, out fields are defined on the full or half lattice.

  virtual int RitzLatSize();
     // RitzLatSize returns the size of a fermion on a node
     // It uses the RitzMatType flag to determine the operator
     // to use and relevant checkerboard sizes.

  //! The matrix inversion used in to solve M out = in
  int BiCGstab(Vector *out,
	       Vector *in,
	       Float src_norm_sq,
	       int n,
	       Float *true_res);
     // The Conjugate Gradient inverter. 
     // Source is *in, initial guess and solution is *out.
     // src_norm_sq is the norm of the source squared 
     // (dot product of source * source). If src_norm_sq=0 
     // it is calculated inside InvCg. 
     // n determines which BiCGstab(n) algorithm we do
     // If true_res !=0 the value of the true residual is returned
     // in true_res.
     // *true_res = |src - MatPcDagMatPc * sol| / |src|
  
};



//-----------------------------------------------------------------
//! A class describing the Dirac operator for Wilson fermions.
/*!
  See the description of the DiracOpWilsonTypes class for the definition
  of the clover fermion matrix.

  The constructor changes the storage order of the gauge field to ::WILSON
  order This change persists throughout the lifetime of the object.

  Only one instance of this class is allowed to be in existence at any time. 
*/
//------------------------------------------------------------------
class DiracOpWilson : public DiracOpWilsonTypes
{
 private:
  const char *cname;    // Class name.

 public:

  void *wilson_lib_arg;  // pointer to an argument structure related
                         // to the wilson library.
  Float kappa;           // Wilson kappa = 1 /[2 (4 +mass) ]

  DiracOpWilson(Lattice& latt,            // Lattice object.
		Vector *f_field_out,      // Output fermion field ptr.
		Vector *f_field_in,       // Input fermion field ptr.
		CgArg *arg,               // Argument structure
		CnvFrmType convert);  // Fermion conversion flag

  virtual ~DiracOpWilson();

  void DiracArg(CgArg *arg);
     // It sets the dirac_arg pointer to arg and initializes
     // kappa.

  void MatPcDagMatPc(Vector *out, Vector *in, Float *dot_prd=0);
     // MatPcDagMatPc is the fermion matrix that appears in the HMC 
     // evolution. It is a Hermitian matrix where M is
     // the even/odd preconditioned Dirac Operator matrix.        
     // MatPcDagMatPc connects only even-->even or odd-->odd sites.
     // The in, out fields are defined on a checkerboard.
     // If dot_prd is not 0 then the dot product (on node)
     // <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.

  void Dslash(Vector *out, 
		      Vector *in,
		      ChkbType cb, 
		      DagType dag);
     // Dslash is the derivative part of the fermion matrix. 
     // Dslash conects only odd-->even or even-->odd sites.
     // The in, out fields are defined on a checkerboard.
     // cb refers to the checkerboard of the in field.

  //! Multiplication by the odd-even preconditioned fermion matrix.
  void MatPc(Vector *out, Vector *in);

  //! Multiplication by the hermitian conjugate of the odd-even preconditioned fermion matrix.
  void MatPcDag(Vector *out, Vector *in);

  int MatInv(Vector *out, 
	     Vector *in, 
	     Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
     // The inverse of the unconditioned Dirac Operator 
     // using Conjugate gradient. source is *in, initial
     // guess and solution is *out.
     // If true_res !=0 the value of the true residual is returned
     // in true_res.
     // *true_res = |src - MatPcDagMatPc * sol| / |src| 
     // prs_in is used to specify if the source
     // in should be preserved or not. If not the memory usage
     // is less by half the size of a fermion vector.
     // The function returns the total number of CG iterations.

  int MatInv(Vector *out, 
	     Vector *in,
	     PreserveType prs_in = PRESERVE_YES);
     // Same as original but true_res=0.

  int MatInv(Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
     // Same as original but in = f_in and out = f_out.

  int MatInv(PreserveType prs_in = PRESERVE_YES);
     // Same as original but in = f_in, out = f_out, true_res=0.

  void MatHerm(Vector *out, Vector *in);

  void Mat(Vector *out, Vector *in);

  void MatDag(Vector *out, Vector *in);

  void CalcHmdForceVecs(Vector *chi) ;
  //!< Computes vectors used in the HMD pseudofermionic force term.
  
    // GRF
    // chi is the solution to MatPcInv.  The user passes two full size
    // CANONICAL fermion vectors with conversion enabled to the
    // constructor.  Using chi, the function fills these vectors;
    // the result may be used to compute the HMD fermion force.

#ifdef USE_QUDA
  int QudaInvert(Vector *out, Vector *in, Float *true_res, int mat_type);
#endif


};

//-----------------------------------------------------------------

//-----------------------------------------------------------------
//! A class describing the Dirac operator for Naive fermions.
//------------------------------------------------------------------
class DiracOpNaive : public DiracOpWilsonTypes
{
 private:
  char *cname;    // Class name.

 public:

  void *wilson_lib_arg;  // pointer to an argument structure related
                         // to the wilson library.
  Float kappa;           // Wilson kappa = 1 /[2 (4 +mass) ]

  DiracOpNaive(Lattice& latt,            // Lattice object.
		Vector *f_field_out,      // Output fermion field ptr.
		Vector *f_field_in,       // Input fermion field ptr.
		CgArg *arg,               // Argument structure
		CnvFrmType convert);  // Fermion conversion flag

  virtual ~DiracOpNaive();

  void DiracArg(CgArg *arg);
     // It sets the dirac_arg pointer to arg and initializes
     // kappa.

  void MatPcDagMatPc(Vector *out, Vector *in, Float *dot_prd=0);
     // MatPcDagMatPc is the fermion matrix that appears in the HMC 
     // evolution. It is a Hermitian matrix where M is
     // the even/odd preconditioned Dirac Operator matrix.        
     // MatPcDagMatPc connects only even-->even or odd-->odd sites.
     // The in, out fields are defined on a checkerboard.
     // If dot_prd is not 0 then the dot product (on node)
     // <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.

  void Dslash(Vector *out, 
		      Vector *in,
		      ChkbType cb, 
		      DagType dag);
     // Dslash is the derivative part of the fermion matrix. 
     // Dslash conects only odd-->even or even-->odd sites.
     // The in, out fields are defined on a checkerboard.
     // cb refers to the checkerboard of the in field.

  //! Multiplication by the odd-even preconditioned fermion matrix.
  void MatPc(Vector *out, Vector *in);

  //! Multiplication by the hermitian conjugate of the odd-even preconditioned fermion matrix.
  void MatPcDag(Vector *out, Vector *in);

  int MatInv(Vector *out, 
	     Vector *in, 
	     Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
     // The inverse of the unconditioned Dirac Operator 
     // using Conjugate gradient. source is *in, initial
     // guess and solution is *out.
     // If true_res !=0 the value of the true residual is returned
     // in true_res.
     // *true_res = |src - MatPcDagMatPc * sol| / |src| 
     // prs_in is used to specify if the source
     // in should be preserved or not. If not the memory usage
     // is less by half the size of a fermion vector.
     // The function returns the total number of CG iterations.

  int MatInv(Vector *out, 
	     Vector *in,
	     PreserveType prs_in = PRESERVE_YES);
     // Same as original but true_res=0.

  int MatInv(Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
     // Same as original but in = f_in and out = f_out.

  int MatInv(PreserveType prs_in = PRESERVE_YES);
     // Same as original but in = f_in, out = f_out, true_res=0.

  void MatHerm(Vector *out, Vector *in);

  void Mat(Vector *out, Vector *in);

  void MatDag(Vector *out, Vector *in);

  void CalcHmdForceVecs(Vector *chi) ;
  //!< Computes vectors used in the HMD pseudofermionic force term.
  
    // GRF
    // chi is the solution to MatPcInv.  The user passes two full size
    // CANONICAL fermion vectors with conversion enabled to the
    // constructor.  Using chi, the function fills these vectors;
    // the result may be used to compute the HMD fermion force.

};

//-----------------------------------------------------------------
//! ~~ A class describing the Dirac operator for twisted-mass Wilson fermions.
/*!
  The constructor changes the storage order of the gauge field to ::WILSON
  order This change persists throughout the lifetime of the object.

  Only one instance of this class is allowed to be in existence at any time. 
*/
//------------------------------------------------------------------
class DiracOpWilsonTm : public DiracOpWilson
{
 private:
  const char *cname;    // Class name.

  Float epsilon;	 // fractional imaginary mass 
  Float ctheta, stheta;	 // gamma_5(theta) parameters

 public:

//
// Functions DiracOpWilson to CalcHmdForceVecs are as in DiracOpWilson 
// with twisted mass term - gamma_5(theta) added where necessary
//

  DiracOpWilsonTm(Lattice& latt,            // Lattice object.
		Vector *f_field_out,      // Output fermion field ptr.
		Vector *f_field_in,       // Input fermion field ptr.
		CgArg *arg,               // Argument structure
		CnvFrmType convert);  // Fermion conversion flag

  virtual ~DiracOpWilsonTm();

  void DiracArg(CgArg *arg);
     // It sets the dirac_arg pointer to arg and initializes
     // kappa.

  void MatPcDagMatPc(Vector *out, Vector *in, Float *dot_prd=0);

  void Dslash(Vector *out, 
		      Vector *in,
		      ChkbType cb, 
		      DagType dag);

  void Dslash_tm(Vector *out, 
		      Vector *in,
		      ChkbType cb, 
		      DagType dag);
//
  void MatPc(Vector *out, Vector *in);

  void MatPcDag(Vector *out, Vector *in);

  void CalcHmdForceVecs(Vector *chi) ;

  void CalcBsnForceVecs(Vector *chi, Vector *eta) ;
 
//
// Function MatInv to MatDag are not (yet) implemented because
// no calls of these functions have been found.  Correct twisted 
// mass modification depends on calling context.

#ifdef USE_BFM_TM
  int InvCg(Vector *out, Vector *in,
	    Float src_norm_sq, Float *true_res);
//  int InvCg(Vector *out, Vector *in, Float *true_res);
  int InvCg(Float *true_res){
 InvCg(f_out,f_in,0,true_res);
 }
#endif
//

  int MatInv(Vector *out, 
	     Vector *in, 
	     Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
     // The inverse of the unconditioned Dirac Operator 
     // using Conjugate gradient.
     
  int MatInv(Vector *out, 
	     Vector *in,
	     PreserveType prs_in = PRESERVE_YES);
     
  int MatInv(Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
     
  int MatInv(PreserveType prs_in = PRESERVE_YES);
     
  void MatHerm(Vector *out, Vector *in);

  void Mat(Vector *out, Vector *in);

  void MatDag(Vector *out, Vector *in);


  //  virtual void RitzMat(Vector *out, Vector *in);
     // RitzMat is the fermion matrix used in Ritz
     // RitzMat works on the full lattice or half lattice
     // The in, out fields are defined on the full or half lattice.

  
};

//------------------------------------------------------------------
// Forward struct declaration for class DiracOpClover
//------------------------------------------------------------------
struct Clover;
class Fclover;

//------------------------------------------------------------------
//! A class describing the Dirac operator for clover-improved Wilson fermions.
/*!
  See the description of the DiracOpWilsonTypes class for the definition
  of the clover fermion matrix.

  The constructor changes the storage order of the gauge field to ::WILSON
  order This change persists throughout the lifetime of the object.

  Only one instance of this class is allowed to be in existence at any time.
*/
//------------------------------------------------------------------
class DiracOpClover : public DiracOpWilsonTypes
{
 private:
  const char *cname;            // Class name.
  Clover *clover_lib_arg; // pointer to an argument structure related
                          // to the clover library.
  Float kappa;    /* For isotropic lattices:
                   *     kappa = 1 / [ 2 (4 +mass) ]
                   * For anisotropic lattices:
                   *     kappa = velocity / { 2 [xi (1 +mass) + 3 velocity]} 
                   */
  Float omega;    /* For isotropic lattices:
                   *     omega_t = omega_s = kappa * Csw
                   * For anisotropic lattices:
                   *     omega_t = Csw_t vel^2 / {xi 2 [xi (1 +mass) + 3 vel]} 
                   *     omega_s = Csw_s       / {   2 [xi (1 +mass) + 3 vel]} 
                   */
 
  Float omega_xi; // see comments for omega
 
 public:
  DiracOpClover(Lattice& latt,            // Lattice object.
		Vector *f_field_out,      // Output fermion field ptr.
		Vector *f_field_in,       // Input fermion field ptr.
		CgArg *arg,               // Argument structure
		CnvFrmType convert);  // Fermion conversion flag

  virtual ~DiracOpClover();

  // Auxiliary member functions
  //----------------------------------------------------------------------
  const Matrix& GetLink(const int *site, int dir) const;
  //!< Gets a reference to a link.
  // The site[] is [x,y,z,t] and could be off-node.
  // If off-node, the link is stored in a static buffer and will be
  // overwritten by the next off-node link, so BE CAREFUL.

  void SiteFuv(Matrix &Fuv, const int *site, int mu, int nu) const;
  //!< Computes the gauge field strength tensor for a specified site and orientation.
  // calculate Fuv, on return Fuv follows the matrix convention
  // of IFloat[row][col][2], instead of IFloat[col][row][2] used in
  // Wilson order.

  void SiteCloverMat(const int *site, IFloat *mat_72_IFloats) const;
  //!< Computes the clover matrix at a lattice site.
  // calculate the local 12x12 clover matrix, and stored it as
  // two blocks of 6x6 compressed (that is, only the lower diagonal
  // part is stored) hermitian matrices.

  void CloverMatChkb(ChkbType chkb, int inverse = 0) const;
  //!< Computes the clover matrix or its inverse at all sites with a given parity.  

  //! Multiplication by the odd-even preconditioned fermion matrix or its hermitian conjugate.
  void MatPcDagOrNot(Vector *out, const Vector *in, int dag, Float *dot=NULL) const;

  // FOR THE PURPOSE OF DEBUGGING ONLY
  //----------------------------------------------------------------------
//! For  debugging  purposes only, apparently
  void MatDagOrNotDbg(Vector *out, const Vector *in, int dag, int direct) const;
  // out = MatDag in    OR  out = Mat in

//! For  debugging  purposes only, apparently
  void MatDagMatDbg(Vector *out, Vector *in, Float *dot_prd=0, int direct=0);
  // out = MatDagMat in
  // If direct = 0,   use MatPc to implement
  //           = 1,   calculate directly

  // More functions
  //----------------------------------------------------------------------
  void DiracArg(CgArg *arg);
     // It sets the dirac_arg pointer to arg and initializes
     // kappa

  void MatPcDagMatPc(Vector *out, Vector *in, Float *dot_prd=0);
     // MatPcDagMatPc is the Hermitian matrix M^dag M, where M is
     // the even/odd preconditioned Dirac Operator matrix.        
     // MatPcDagMatPc connects only odd-->odd sites.
     // The in, out fields are defined on the odd checkerboard.
     // If dot_prd is not 0 then the dot product (on node)
     // <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.

  void Dslash(Vector *out, 
		      Vector *in,
		      ChkbType cb, 
		      DagType dag);
     // Not implemented

  //! Multiplication by the odd-even preconditioned fermion matrix.
/*!
  The vectors are defined on odd parity lattice sites and the multiplication
  is by the odd parity part of the matrix.
  \param out The resulting vector.
  \param in The vector to be multiplied.
*/
  void MatPc(Vector *out, Vector *in) {
    MatPcDagOrNot(out, in, 0);
  }

  //! Multiplication by the  hermitian conjugate odd-even preconditioned fermion matrix.
/*!
  The vectors are defined on odd parity lattice sites and the multiplication
  is by the odd parity part of the matrix.
  \param out The resulting vector.
  \param in The vector to be multiplied.
*/
  void MatPcDag(Vector *out, Vector *in) {
    MatPcDagOrNot(out, in, 1);
  }
     // MatPcDag(Vector *out, Vector *in) :
     // MatPcDag is the dagger of the preconditioned fermion matrix.  
     // MatPcDag connects only odd-->odd sites.
     // The in, out fields are defined on the odd checkerboard.

  int MatInv(Vector *out, 
	     Vector *in, 
	     Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
     // It calculates out where A * out = in and
     // A is the fermion matrix (Dirac operator) with no 
     // preconditioning. The preconditioned matrix is inverted 
     // and from the result the non-preconditioned out is 
     // calculated . The inversion of the odd checkerboard
     // piece is done with the conjugate gradient algorithm
     // while the inversion of the even checkerboard is done
     // using standard explicit hermitian matrix inversion of the
     // clover matrix. cg_arg is the structure that contains
     // all the control parameters for the CG, in is the
     // fermion field source vector, out should be set to be
     // the initial guess and on return is the solution.
     // in and out are defined on the whole latice.
     // If true_res !=0 the value of the true residual of the CG and
     // is returned in true_res.
     // *true_res = |src - MatPcDagMatPc * sol| / |src|
     // prs_in is used to specify if the source
     // in should be preserved or not. If not the memory usage
     // is less by the size of a fermion vector.
     // The function returns the total number of CG iterations.

  int MatInv(Vector *out, 
	     Vector *in,
	     PreserveType prs_in = PRESERVE_YES);
     // Same as original but true_res=0.

  int MatInv(Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
     // Same as original but in = f_in and out = f_out.

  int MatInv(PreserveType prs_in = PRESERVE_YES);
     // Same as original but in = f_in, out = f_out, true_res=0.

  int MatEvlInv(Vector *out, Vector *in, Float *true_res);
  //!< The matrix inversion used in the molecular dynamics evolution.
     // It calculates out where A * out = in and
     // A is the preconditioned fermion matrix that appears
     // in the HMC evolution (even/odd  preconditioning 
     // of [Dirac^dag Dirac]. The inversion of the odd checkerboard
     // piece is done with the conjugate gradient algorithm
     // while the inversion of the even checkerboard is done
     // using standard explicit hermitian matrix inversion of the
     // clover matrix. cg_arg is the structure that contains
     // all the control parameters for the CG, in is the
     // fermion field source vector, out should be set to be
     // the initial guess and on return is the solution.
     // in and out are defined on the whole latice.
     // If true_res !=0 the value of the true residual of the CG and
     // is returned in true_res.
     // *true_res = |src - MatPcDagMatPc * sol| / |src|
     // The function returns the total number of CG iterations.

  int MatEvlInv(Vector *out, Vector *in);
  //!< The matrix inversion used in the molecular dynamics evolution.
     // Same as original but true_res=0.

  int MatEvlInv(Float *true_res);
    //!< The matrix inversion used in the molecular dynamics evolution.
     // Same as original but in = f_in and out = f_out.

  int MatEvlInv(void);
  //!< The matrix inversion used in the molecular dynamics evolution.  
     // Same as original but in = in, f_out = f_out, true_res=0.

  void MatHerm(Vector *out, Vector *in);
     // MatHerm is the hermitian version of Mat.
     // MatHerm works on the full lattice.
     // The in, out fields are defined on the ful.

  void Mat(Vector *out, Vector *in);
  //!< Not implemented.
     // Mat is the unpreconditioned fermion matrix.  
     // Mat works on the full lattice
     // The in, out fields are defined on the full lattice.

  void MatDag(Vector *out, Vector *in);
  //!< Not implemented.  
     // MatDag is the dagger of the unpreconditioned fermion matrix. 
     // MatDag works on the full lattice
     // The in, out fields are defined on the full lattice.

  void CalcHmdForceVecs(Vector *chi) ;
  //!< Computes vectors used in the HMD pseudofermionic force term.
  
    // Lingling changed the prvious defintion of the fermion fields by
    // some constant factors. chi is the solution to MatPcInv on the  
    // odd sites.  The user passes two full size  
    // CANONICAL fermion vectors with conversion enabled to the
    // constructor.  Using chi, the function fills these vectors;
    // the result may be used to compute the HMD fermion force.

};
CPS_END_NAMESPACE

#include<util/dirac_op/d_op_dwf.h>

CPS_START_NAMESPACE

//------------------------------------------------------------------
// A class describing the Dirac operator for Mobius domain-wall fermions.
//------------------------------------------------------------------
class DiracOpMdwf : public DiracOpWilsonTypes
{
 private:
  const char *cname;    // Class name.

  Float epsilon;
  int max_num_iter;
  int use_single_precision;

  struct QOP_MDWF_State * mdwf_state;
  struct QOP_MDWF_Parameters * mdwf_param;
  // gauge field
  union {
    struct QOP_F3_MDWF_Gauge * f;
    struct QOP_D3_MDWF_Gauge * d;
  } mdwf_gauge_ptr;
  // pointer to the lattice instance
  // as a consequence, you can't destroy the lattice instance while
  // this class is in use.
  Lattice * latt_ptr;
 protected:
  void DiracOpGlbSum(Float *float_p);					
     // The global sum used by InvCg. If s_nodes = 1
     // it is the usual global sum. If s_nodes > 1 it
     // is the 5-dimensional globals sum glb_sum_five.


 public:
  DiracOpMdwf(Lattice& latt,            // Lattice object.
              MdwfArg * mdwf_arg_p);

  virtual ~DiracOpMdwf();

  void MatPcDagMatPc(Vector *out, Vector *in, Float *dot_prd=0);
     // MatPcDagMatPc is the fermion matrix that appears in the HMC 
     // evolution. It is a Hermitian matrix.
     // The in, out fields are defined on the checkerboard lattice.
     // If dot_prd is not 0 then the dot product (on node)
     // <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.

  void MatPcDagMatPcShift(Vector *out, Vector *in, Float *dot_prd=0);
  // MatPcDagMatPc is the fermion matrix that appears in the HMC 
  // evolution. It is a Hermitian matrix.
  // The in, out fields are defined on the checkerboard lattice.
  // If dot_prd is not 0 then the dot product (on node)
  // <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.
  //
  // When dirac_arg->eigen_shift is non-zero, it shift the spectrum of matrix:
  //    MatPcDagMatPc = H^2  ->  (H-shift)(H-shift)
  // where H = Gamma_5 MatPc
  //
  // For other fermions, one could also implement similar shifts.
  // For wilson, H = gamma_5 MatPC .
  //------------------------------------------------------------------
  
  
  void Dslash(Vector *out, 
		      Vector *in,
		      ChkbType cb, 
		      DagType dag);
     // Dslash is the derivative part of the fermion matrix. 
     // The in, out fields are defined on the checkerboard lattice
     // cb = 0/1 <--> even/odd checkerboard of in field.
     // dag = 0/1 <--> Dslash/Dslash^dagger is calculated.

  //! Multiplication by the odd-even preconditioned fermion matrix.
  void MatPc(Vector *out, Vector *in);
     // MatPc is the fermion matrix.  
     // The in, out fields are defined on the checkerboard lattice.

  //! Multiplication by the  hermitian conjugate odd-even preconditioned fermion matrix.
  void MatPcDag(Vector *out, Vector *in);
     // MatPcDag is the dagger of the fermion matrix. 
     // The in, out fields are defined on the checkerboard lattice.

  int MatInv(Vector *out, 
	     Vector *in, 
	     Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
     // The inverse of the unconditioned Dirac Operator 
     // using Conjugate gradient.  source is *in, initial
     // guess and solution is *out.
     // If true_res !=0 the value of the true residual is returned
     // in true_res.
     // *true_res = |src - MatPcDagMatPc * sol| / |src|
     // prs_in is used to specify if the source
     // in should be preserved or not. If not the memory usage
     // is less by half the size of a fermion vector.
     // The function returns the total number of CG iterations.

  int MatInv(Vector *out, 
	     Vector *in,
	     PreserveType prs_in = PRESERVE_YES);
     // Same as original but true_res=0.

  int MatInv(Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
     // Same as original but in = f_in and out = f_out.

  int MatInv(PreserveType prs_in = PRESERVE_YES);
     // Same as original but in = f_in, out = f_out, true_res=0.
 
  void MatHerm(Vector *out, Vector *in);
     // MatHerm is the hermitian version of Mat.
     // MatHerm works on the full lattice.
     // The in, out fields are defined on the ful.

  void Mat5doe(Vector *out, Vector *in);
  void Mat(Vector *out, Vector *in);
     // Mat is the unpreconditioned fermion matrix.  
     // Mat works on the full lattice
     // The in, out fields are defined on the full lattice.

  void MatDag(Vector *out, Vector *in);
     // MatDag is the dagger of the unpreconditioned fermion matrix. 
     // MatDag works on the full lattice
     // The in, out fields are defined on the full lattice.

  //! Multiplication by the square of the fermion matrix.
  virtual void MatDagMat(Vector *out, Vector *in);

  void CalcHmdForceVecs(Vector *chi);
  //!< Computes vectors used in the HMD pseudofermionic force term.
    // GRF
    // chi is the solution to MatPcInv.  The user passes two full size
    // CANONICAL fermion vectors with conversion enabled to the
    // constructor.  Using chi, the function fills these vectors;
    // the result may be used to compute the HMD fermion force.

  void Reflex(Vector *out, Vector *in);
  //!< Not implemented
    // Reflexion in s operator, needed for the hermitian version 
    // of the dirac operator in the Ritz solver.

  void DiracArg(CgArg *arg);
     // It sets the dirac_arg pointer to arg and initializes
     // the relevant parameters (kappa, m^2, ...).

//  virtual void RitzMat(Vector *out, Vector *in);
     // RitzMat is the fermion matrix used in Ritz
     // RitzMat works on the full lattice or half lattice
     // The in, out fields are defined on the full or half lattice.

//  virtual void RitzMat(Vector *out, Vector *in,
//		       MatrixPolynomialArg* cheby_arg);
  
  void MatPcHerm(Vector *out, Vector *in);
     // MatPcHerm is the hermitian version of MatPc.
     // MatHerm works on the 4d-eo half-parity lattice.
     // The in, out fields are defined on the even-even parity.



};


//------------------------------------------------------------------
//! A class describing the Dirac operator for Mobius Wilson fermions.
/*!
  See the description of the DiracOpWilsonTypes class for the definition
  of the Wilson fermion matrix.

  The constructor changes the storage order of the gauge field to ::WILSON
  order This change persists throughout the lifetime of the object.

  Only one instance of this class is allowed to be in existence at any time.
*/
//------------------------------------------------------------------
class DiracOpMobius : public DiracOpWilsonTypes
{
 private:
  char *cname;    // Class name.

  void *mobius_lib_arg;     // pointer to an argument structure related
                         // to the mobius library.

  Float mass;            // Mobius mass (couples left-right components)


  protected:
  void DiracOpGlbSum(Float *float_p);					
     // The global sum used by InvCg. If s_nodes = 1
     // it is the usual global sum. If s_nodes > 1 it
     // is the 5-dimensional globals sum glb_sum_five.


 public:
  DiracOpMobius(Lattice& latt,            // Lattice object.
		Vector *f_field_out,      // Output fermion field ptr.
		Vector *f_field_in,       // Input fermion field ptr.
		CgArg *arg,               // Argument structure
		CnvFrmType convert);  // Fermion conversion flag
  
  virtual ~DiracOpMobius();

  void DiracArg(CgArg *arg);
     // It sets the dirac_arg pointer to arg and initializes
     // kappa

  void MatPcDagMatPc(Vector *out, Vector *in, Float *dot_prd=0);
     // MatPcDagMatPc is the fermion matrix that appears in the HMC 
     // evolution. It is a Hermitian matrix.
     // The in, out fields are defined on the checkerboard lattice.
     // If dot_prd is not 0 then the dot product (on node)
     // <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.

  void MatPcDagMatPcShift(Vector *out, Vector *in, Float *dot_prd=0);
  // MatPcDagMatPc is the fermion matrix that appears in the HMC 
  // evolution. It is a Hermitian matrix.
  // The in, out fields are defined on the checkerboard lattice.
  // If dot_prd is not 0 then the dot product (on node)
  // <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.
  //
  // When dirac_arg->eigen_shift is non-zero, it shift the spectrum of matrix:
  //    MatPcDagMatPc = H^2  ->  (H-shift)(H-shift)
  // where H = Gamma_5 MatPc
  //
  // For other fermions, one could also implement similar shifts.
  // For wilson, H = gamma_5 MatPC .
  //------------------------------------------------------------------
  
  
  void Dslash(Vector *out, 
	      Vector *in,
	      ChkbType cb, 
	      DagType dag);
     // Dslash is the derivative part of the fermion matrix. 
     // The in, out fields are defined on the checkerboard lattice
     // cb = 0/1 <--> even/odd checkerboard of in field.
     // dag = 0/1 <--> Dslash/Dslash^dagger is calculated.

  //! Multiplication by the odd-even preconditioned fermion matrix.
  void MatPc(Vector *out, Vector *in);
     // MatPc is the fermion matrix.  
     // The in, out fields are defined on the checkerboard lattice.

  //! Multiplication by the  hermitian conjugate odd-even preconditioned fermion matrix.
  void MatPcDag(Vector *out, Vector *in);
     // MatPcDag is the dagger of the fermion matrix. 
     // The in, out fields are defined on the checkerboard lattice.

  int MatInv(Vector *out, 
	     Vector *in, 
	     Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
     // The inverse of the unconditioned Dirac Operator 
     // using Conjugate gradient.  source is *in, initial
     // guess and solution is *out.
     // If true_res !=0 the value of the true residual is returned
     // in true_res.
     // *true_res = |src - MatPcDagMatPc * sol| / |src|
     // prs_in is used to specify if the source
     // in should be preserved or not. If not the memory usage
     // is less by half the size of a fermion vector.
     // The function returns the total number of CG iterations.

  int MatInv(Vector *out, 
	     Vector *in,
	     PreserveType prs_in = PRESERVE_YES);
  // Same as original but true_res=0.
  
  int MatInv(Float *true_res,
	     PreserveType prs_in = PRESERVE_YES);
  // Same as original but in = f_in and out = f_out.
  
  int MatInv(PreserveType prs_in = PRESERVE_YES);
  // Same as original but in = f_in, out = f_out, true_res=0.
  
  void MatHerm(Vector *out, Vector *in);
  // MatHerm is the hermitian version of Mat.
  // MatHerm works on the full lattice.
  // The in, out fields are defined on the ful.
  
  void Mat5doe(Vector *out, Vector *in);
  void Mat(Vector *out, Vector *in);
  // Mat is the unpreconditioned fermion matrix.  
  // Mat works on the full lattice
  // The in, out fields are defined on the full lattice.
  
  void MatDag(Vector *out, Vector *in);
  // MatDag is the dagger of the unpreconditioned fermion matrix. 
  // MatDag works on the full lattice
  // The in, out fields are defined on the full lattice.
  
  void CalcHmdForceVecs(Vector *chi) ;
  //!< Computes vectors used in the HMD pseudofermionic force term.
  // GRF
  // chi is the solution to MatPcInv.  The user passes two full size
  // CANONICAL fermion vectors with conversion enabled to the
  // constructor.  Using chi, the function fills these vectors;
  // the result may be used to compute the HMD fermion force.
  
  void Reflex(Vector *out, Vector *in);
  //!< Not implemented
  // Reflexion in s operator, needed for the hermitian version 
  // of the dirac operator in the Ritz solver.
  
  virtual void RitzMat(Vector *out, Vector *in);
  // RitzMat is the fermion matrix used in Ritz
  // RitzMat works on the full lattice or half lattice
  // The in, out fields are defined on the full or half lattice.
  
  virtual void RitzMat(Vector *out, Vector *in,
		       MatrixPolynomialArg* cheby_arg);
  
  void MatPcHerm(Vector *out, Vector *in);
  // MatPcHerm is the hermitian version of MatPc.
  // MatHerm works on the 4d-eo half-parity lattice.
  // The in, out fields are defined on the even-even parity.

  // Hantao's operator used to construct Hermitian mobius Dirac op
  void Dminus(Vector *out, Vector *in);

};



CPS_END_NAMESPACE
#endif




