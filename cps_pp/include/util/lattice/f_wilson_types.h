#ifndef INCLUDED_F_WILSON_TYPES_H
#define INCLUDED_F_WILSON_TYPES_H           //!< Prevent multiple inclusion

CPS_START_NAMESPACE

//------------------------------------------------------------------
//! A class containing methods relevant to all Wilson type fermion actions.
//------------------------------------------------------------------
class FwilsonTypes : public virtual Lattice
{
 private:
    const char *cname;    // Class name.
    
 public:
// protected:
    void (*sproj_tr[8])(IFloat *f, 
			IFloat *v, 
			IFloat *w, 
			int num_blk, 
			int v_stride,
			int w_stride) ;
    //!< Array of pointers to external functions.
    /*!<
      These functions compute 
      \f$ f_{ij} = Tr_{spins}[ (1 \pm \gamma_\mu) v_i w^\dagger_j ] \f$
      for spin-colour vectors \a v and \a w where \a i and \a j are
      colour indices.
    */
    // Array with entries that point to 12 non-member functions.
    // These functions are called as follows:
    // Sigmaproj_tr[SigmaprojType mu_nu]();
    // For the various SprojTypes see enum.h
    // These functions return a color matrix in f constructed from
    // the spinors v, w using: 
    // f_(i,j) = Tr_spin[ (1 +/- gamma_mu) v_i w^dag_j ]
    //
    // num_blk is the number of spinors v, w. The routines 
    // accumulate the sum over spinors in f.
    //
    // v_stride and w_stride are the number of Floats between spinors
    // (a stride = 0 means that the spinors are consecutive in memory)

    void (*Sigmaproj_tr[12])(IFloat *f, 
			IFloat *v, 
			IFloat *w, 
			int num_blk, 
			int v_stride,
			int w_stride) ;
    //!< Array of pointers to external functions.
    /*!<
      These functions compute 
    \f$ f_{ij} = \frac{1}{2} Tr_{spins}[ \sigma_{\mu\nu} v_i w^\dagger_j ] \f$
      for spin-colour vectors \a v and \a w where \a i and \a j
      are colour indices.
    */
    // Array with entries that point to 12 non-member functions.
    // These functions are called as follows:
    // Sigmaproj_tr[SigmaprojType mu_nu]();
    // For the various SigmaprojTypes see enum.h 
    // These functions return a color matrix in f constructed from
    // the spinors v, w using: 
    // f_(i,j) = 1/2 Tr_spin[ Sigma_{mu,nu} v_i w^dag_j ]
    //
    // num_blk is the number of spinors v, w. The routines 
    // accumulate the sum over spinors in f.
    //
    // v_stride and w_stride are the number of Floats between spinors
    // (a stride = 0 means that the spinors are consecutive in memory)

 public:

    FwilsonTypes();

    virtual ~FwilsonTypes();

    //! Multiplication of a lattice spin-colour vector by gamma_5.
    void Gamma5(Vector *v_out, Vector *v_in, int num_sites);

    void Fconvert(Vector*, StrOrdType, StrOrdType, int cb=2);

    Float FhamiltonNode( Vector*,  Vector*) ;

    int FsiteOffsetChkb(const int*) const;

    int FsiteOffset(const int*) const;

    int FsiteSize() const;
	
    virtual FclassType Fclass() const = 0;
    
    int SpinComponents() const;

    int ExactFlavors() const;
    
    //!< Method to ensure bosonic force works (does nothing for Wilson
    //!< theories.
    void BforceVector(Vector *in, CgArg *cg_arg);
};


//------------------------------------------------------------------
//! A class implementing Wilson fermions.
/*!
  \ingroup factions
*/
//------------------------------------------------------------------
class Fwilson : public virtual FwilsonTypes
{
 private:
    const char *cname;    // Class name.
    
 public:

    Fwilson();

    virtual ~Fwilson();

    int FsiteSize() const;

    FclassType Fclass() const; 

    int FchkbEvl() const;
	// returns 1 => The fermion fields in the evolution
        //      or the CG that inverts the evolution matrix
	//      are defined on a single checkerboard (half the 
	//      lattice).

    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm = CNV_FRM_YES);
        // It calculates f_out where A * f_out = f_in and
        // A is the preconditioned fermion matrix that appears
        // in the HMC evolution (even/odd  preconditioning 
        // of [Dirac^dag Dirac]. The inversion is done
	// with the conjugate gradient. cg_arg is the structure
        // that contains all the control parameters, f_in is the
        // fermion field source vector, f_out should be set to be
        // the initial guess and on return is the solution.
	// f_in and f_out are defined on a checkerboard.
        // If true_res !=0 the value of the true residual is returned
        // in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|
	// The function returns the total number of CG iterations.

    int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
		    int Nshift, int isz, CgArg **cg_arg, 
		    CnvFrmType cnv_frm, MultiShiftSolveType type, 
		    Float *alpha, Vector **f_out_d);
    //!< The matrix inversion used in the molecular dynamics algorithms.
    /*!<
      Solves \f$ (M^\dagger M + shift) f_{out} = f_{in} \f$ for \f$ f_{out} \f$,
      where \a M is the (possibly odd-even preconditioned) fermionic matrix.

      \param f_out The solution vectors.
      \param f_in The source vector
      \param shift The shifts of the fermion matrix.
      \param Nshift The number of shifts
      \param isz The smallest shift (required by MInvCG)
      \param cg_arg The solver parameters
      \param cnv_frm Whether the lattice fields need to be converted to
      to a new storage order appropriate for the type of fermion action.
      If this is ::CNV_FRM_NO, then just the gauge field is converted.
      If this is ::CNV_FRM_YES, then the fields \a f_out and \a f_in
      are also converted: This assumes they are initially in the same order as
      the gauge field. Fields that are converted are restored to their original
      order upon exit of this method. \e N.B. If the fields are already in the
      suitable order, then specifying ::CNV_FRM_YES here has no effect.
      \param type The type of multimass inverter.
      If type == MULTI, then regular multishift inversion is performed, each solution
      stored separately.  If type == SINGLE, the there is a single solution vector, and
      each solution is summed to this vector with amount alpha.
      \param alpha The contribution of each shifted solution to the total solution vector
      \param f_out_d ?
      \return The number of solver iterations.
      \post \a f_out contains the solution vector.
    */

    void FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
		     Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm);

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		Float *true_res,
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
        // It calculates f_out where A * f_out = f_in and
        // A is the fermion matrix (Dirac operator). The inversion
	// is done with the conjugate gradient. cg_arg is the 
        // structure that contains all the control parameters, f_in 
        // is the fermion field source vector, f_out should be set 
        // to be the initial guess and on return is the solution.
	// f_in and f_out are defined on the whole lattice.
        // If true_res !=0 the value of the true residual is returned
        // in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|
        // cnv_frm is used to specify if f_in should be converted 
        // from canonical to fermion order and f_out from fermion 
        // to canonical. 
        // prs_f_in is used to specify if the source
        // f_in should be preserved or not. If not the memory usage
        // is less by half the size of a fermion vector.
	// The function returns the total number of CG iterations.

    int FeigSolv(Vector **f_eigenv, Float *lambda,
		 LanczosArg *eig_arg, 
		 CnvFrmType cnv_frm = CNV_FRM_YES){ return 0; };
    int FeigSolv(Vector **f_eigenv, Float *lambda, 
		 Float *chirality, int *valid_eig,
		 Float **hsum,
		 EigArg *eig_arg, 
		 CnvFrmType cnv_frm = CNV_FRM_YES);
        // It finds the eigenvectors and eigenvalues of A where
        // A is the fermion matrix (Dirac operator). The solution
	// uses Ritz minimization. eig_arg is the 
        // structure that contains all the control parameters, f_eigenv
        // are the fermion field source vectors which should be
        // defined initially, lambda are the eigenvalues returned 
        // on solution. f_eigenv is defined on the whole lattice.
	// The function returns the total number of Ritz iterations.

    Float SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
		 Float mass, DagType dag);
	// It sets the pseudofermion field phi from frm1, frm2.

    ForceArg EvolveMomFforce(Matrix *mom, Vector *frm, 
				 Float mass, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the fermion force. 

    ForceArg EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
			  Float mass, Float step_size);
        // It evolve the canonical momentum mom  by step_size
        // using the bosonic quotient force.

    ForceArg RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
			      int isz, Float *alpha, Float mass, Float dt,
			      Vector **sol_d, ForceMeasure measure);

    Float BhamiltonNode(Vector *boson, Float mass);
        // The boson Hamiltonian of the node sublattice.

};


//------------------------------------------------------------------
//! A class implementing Wilson fermions.
/*!
  \ingroup factions
*/
//------------------------------------------------------------------
class Fnaive : public virtual FwilsonTypes
{
 private:
    char *cname;    // Class name.
    
 public:

    Fnaive();

    virtual ~Fnaive() ;

    int FsiteSize() const;

    FclassType Fclass() const; 

    int FchkbEvl() const;
	// returns 1 => The fermion fields in the evolution
        //      or the CG that inverts the evolution matrix
	//      are defined on a single checkerboard (half the 
	//      lattice).

    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm = CNV_FRM_YES)
		{ ERR.NotImplemented(cname,"FmatEvlInv()");
			return -1;}
        // It calculates f_out where A * f_out = f_in and
        // A is the preconditioned fermion matrix that appears
        // in the HMC evolution (even/odd  preconditioning 
        // of [Dirac^dag Dirac]. The inversion is done
	// with the conjugate gradient. cg_arg is the structure
        // that contains all the control parameters, f_in is the
        // fermion field source vector, f_out should be set to be
        // the initial guess and on return is the solution.
	// f_in and f_out are defined on a checkerboard.
        // If true_res !=0 the value of the true residual is returned
        // in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|
	// The function returns the total number of CG iterations.

    int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
		    int Nshift, int isz, CgArg **cg_arg, 
		    CnvFrmType cnv_frm, MultiShiftSolveType type, 
		    Float *alpha, Vector **f_out_d)
		{ ERR.NotImplemented(cname,"FmatEvlMInv()");
			return -1;}
    //!< The matrix inversion used in the molecular dynamics algorithms.
    /*!<
      Solves \f$ (M^\dagger M + shift) f_{out} = f_{in} \f$ for \f$ f_{out} \f$,
      where \a M is the (possibly odd-even preconditioned) fermionic matrix.

      \param f_out The solution vectors.
      \param f_in The source vector
      \param shift The shifts of the fermion matrix.
      \param Nshift The number of shifts
      \param isz The smallest shift (required by MInvCG)
      \param cg_arg The solver parameters
      \param cnv_frm Whether the lattice fields need to be converted to
      to a new storage order appropriate for the type of fermion action.
      If this is ::CNV_FRM_NO, then just the gauge field is converted.
      If this is ::CNV_FRM_YES, then the fields \a f_out and \a f_in
      are also converted: This assumes they are initially in the same order as
      the gauge field. Fields that are converted are restored to their original
      order upon exit of this method. \e N.B. If the fields are already in the
      suitable order, then specifying ::CNV_FRM_YES here has no effect.
      \param type The type of multimass inverter.
      If type == MULTI, then regular multishift inversion is performed, each solution
      stored separately.  If type == SINGLE, the there is a single solution vector, and
      each solution is summed to this vector with amount alpha.
      \param alpha The contribution of each shifted solution to the total solution vector
      \param f_out_d ?
      \return The number of solver iterations.
      \post \a f_out contains the solution vector.
    */

    void FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
		    Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm){};

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		Float *true_res,
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
        // It calculates f_out where A * f_out = f_in and
        // A is the fermion matrix (Dirac operator). The inversion
	// is done with the conjugate gradient. cg_arg is the 
        // structure that contains all the control parameters, f_in 
        // is the fermion field source vector, f_out should be set 
        // to be the initial guess and on return is the solution.
	// f_in and f_out are defined on the whole lattice.
        // If true_res !=0 the value of the true residual is returned
        // in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|
        // cnv_frm is used to specify if f_in should be converted 
        // from canonical to fermion order and f_out from fermion 
        // to canonical. 
        // prs_f_in is used to specify if the source
        // f_in should be preserved or not. If not the memory usage
        // is less by half the size of a fermion vector.
	// The function returns the total number of CG iterations.

    int FeigSolv(Vector **f_eigenv, Float *lambda,
		 LanczosArg *eig_arg, 
		 CnvFrmType cnv_frm = CNV_FRM_YES){ return 0; };
    int FeigSolv(Vector **f_eigenv, Float *lambda, 
		 Float *chirality, int *valid_eig,
		 Float **hsum,
		 EigArg *eig_arg, 
		 CnvFrmType cnv_frm = CNV_FRM_YES);
        // It finds the eigenvectors and eigenvalues of A where
        // A is the fermion matrix (Dirac operator). The solution
	// uses Ritz minimization. eig_arg is the 
        // structure that contains all the control parameters, f_eigenv
        // are the fermion field source vectors which should be
        // defined initially, lambda are the eigenvalues returned 
        // on solution. f_eigenv is defined on the whole lattice.
	// The function returns the total number of Ritz iterations.

    Float SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
		 Float mass, DagType dag);
	// It sets the pseudofermion field phi from frm1, frm2.

    ForceArg EvolveMomFforce(Matrix *mom, Vector *frm, 
			     Float mass, Float step_size)
		{ ERR.NotImplemented(cname,"EvolveMomFforce()");
		ForceArg null_force; return null_force ;}
        // It evolves the canonical momentum mom by step_size
        // using the fermion force. 

    ForceArg EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
			     Float mass, Float step_size)
		{ ERR.NotImplemented(cname,"EvolveMomFforce()");
		ForceArg null_force; return null_force ;}
        // It evolve the canonical momentum mom  by step_size
        // using the bosonic quotient force.

    ForceArg RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
			      int isz, Float *alpha, Float mass, Float dt,
				  Vector **sol_d, ForceMeasure measure)
		{ ERR.NotImplemented(cname,"RHMC_EvolveMomFforce()");
		ForceArg null_force; return null_force ;}

    Float BhamiltonNode(Vector *boson, Float mass)
		{ ERR.NotImplemented(cname,"BhamiltonNode()");
		return 0.;}
        // The boson Hamiltonian of the node sublattice.

};

//------------------------------------------------------------------
//! A class implementing twisted mass Wilson fermions.
/*!
  \ingroup factions
*/
//------------------------------------------------------------------
class FwilsonTm : public virtual Fwilson
{
 private:
    const char *cname;    // Class name.
    
 public:

    FwilsonTm();

    virtual ~FwilsonTm();

    FclassType Fclass() const; 

   //
   // ~~ modified in f_wilsonTm to create wilsonTm fermions 
   // 
    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm = CNV_FRM_YES);

   //
   //~~ the following functions are versions with the epsilon parameter  
   //~~ for twisted mass Wilson fermions; all implemented here
   //
    Float SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
		 Float mass, Float epsilon, DagType dag);

   // ~~ evolves the  momentum mom using the fermion force.
    ForceArg EvolveMomFforce(Matrix *mom, Vector *frm, 
				 Float mass, Float epsilon, Float step_size);

   // ~~ evolves the  momentum mom using the boson force.
    ForceArg EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
				  Float mass, Float epsilon, Float step_size);

    Float BhamiltonNode(Vector *boson, Float mass, Float epsilon);

    //
    //~~ the following functions are "normal" versions without the
    //~~ epsilon parameter; should never be called by wilsonTm fermions
    //~~ and are implemented here as errors
    // 
    Float SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
		 Float mass, DagType dag);
 
    ForceArg EvolveMomFforce(Matrix *mom, Vector *frm, 
				 Float mass, Float step_size);
    ForceArg EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
				  Float mass, Float step_size);
    Float BhamiltonNode(Vector *boson, Float mass);
};


//------------------------------------------------------------------
//! A class implementing clover improved Wilson fermions.
/*!
  \ingroup factions
*/
//------------------------------------------------------------------
class Fclover : public virtual FwilsonTypes
{
 private:
    const char *cname;    // Class name.

    void EvolveMomFforceSupp(Matrix *mom, Vector *v1, Vector *v2,
		 Vector *v3, Vector *v4, Float mass, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the clover contribution of fermion force 

 public:

    Fclover();
        // Among other things the constructor allocates
        // memory for the even/odd checkerpoard clover
        // matrices. aux0_ptr of the base class is set
        // to the pointer of the even checkerboard matrices
        // and aux1_ptr to the odd.

    virtual ~Fclover();

    FclassType Fclass() const; 
        // It returns the type of fermion class

    int FchkbEvl() const;
        // returns 0 => The fermion fields in the evolution
        // are defined on ODD-EVEN checkerboard (whole
        // lattice).

    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm = CNV_FRM_YES);
        // It calculates f_out where A * f_out = f_in and
        // A is the preconditioned fermion matrix that appears
        // in the HMC evolution (even/odd  preconditioning 
        // of [Dirac^dag Dirac]. The inversion of the odd checkerboard
        // piece is done with the conjugate gradient algorithm
        // while the inversion of the even checkerboard is done
        // using standard explicit hermitian matrix inversion of the
        // clover matrix. cg_arg is the structure that contains
        // all the control parameters for the CG, f_in is the
        // fermion field source vector, f_out should be set to be
        // the initial guess and on return is the solution.
        // f_in and f_out are defined on the whole latice.
        // If true_res !=0 the value of the true residual of the CG and
        // is returned in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|
        // The function returns the total number of CG iterations.

    int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
		    int Nshift, int isz, CgArg **cg_arg, 
		    CnvFrmType cnv_frm, MultiShiftSolveType type, Float *alpha,
		    Vector **f_out_d);
    //!< The matrix inversion used in the molecular dynamics algorithms.
    /*!<
      Solves \f$ (M^\dagger M + shift) f_{out} = f_{in} \f$ for \f$ f_{out} \f$,
      where \a M is the (possibly odd-even preconditioned) fermionic matrix.

      \param f_out The solution vectors.
      \param f_in The source vector
      \param shift The shifts of the fermion matrix.
      \param Nshift The number of shifts
      \param isz The smallest shift (required by MInvCG)
      \param cg_arg The solver parameters
      \param cnv_frm Whether the lattice fields need to be converted to
      to a new storage order appropriate for the type of fermion action.
      If this is ::CNV_FRM_NO, then just the gauge field is converted.
      If this is ::CNV_FRM_YES, then the fields \a f_out and \a f_in
      are also converted: This assumes they are initially in the same order as
      the gauge field. Fields that are converted are restored to their original
      order upon exit of this method. \e N.B. If the fields are already in the
      suitable order, then specifying ::CNV_FRM_YES here has no effect.
      \param type The type of multimass inverter.
      If type == MULTI, then regular multishift inversion is performed, each solution
      stored separately.  If type == SINGLE, the there is a single solution vector, and
      each solution is summed to this vector with amount alpha.
      \param alpha The contribution of each shifted solution to the total solution vector
      \param f_out_d ?
      \return The number of solver iterations.
      \post \a f_out contains the solution vector.
    */

    void FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
		     Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm);

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		Float *true_res,
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
    // It calculates f_out where A * f_out = f_in and
    // A is the fermion matrix (Dirac operator) with no 
    // preconditioning. The preconditioned matrix is inverted 
    // and from the result the non-preconditioned f_out is 
    // calculated . The inversion of the odd checkerboard
    // piece is done with the conjugate gradient algorithm
    // while the inversion of the even checkerboard is done
    // using standard explicit hermitian matrix inversion of the
    // clover matrix. cg_arg is the structure that contains
    // all the control parameters for the CG, f_in is the
    // fermion field source vector, f_out should be set to be
    // the initial guess and on return is the solution.
    // f_in and f_out are defined on the whole latice.
    // If true_res !=0 the value of the true residual of the CG and
    // is returned in true_res.
    // *true_res = |src - MatPcDagMatPc * sol| / |src|
    // cnv_frm is used to specify if f_in should be converted 
    // from canonical to fermion order and f_out from fermion 
    // to canonical. 
    // prs_f_in is used to specify if the source
    // f_in should be preserved or not. If not the memory usage
    // is less by the size of one fermion vector.
    // The function returns the total number of CG iterations.


    int FeigSolv(Vector **f_eigenv, Float *lambda,
		 LanczosArg *eig_arg, 
		 CnvFrmType cnv_frm = CNV_FRM_YES){ return 0; };
    int FeigSolv(Vector **f_eigenv, Float *lambda,
		 Float *chirality, int *valid_eig,
		 Float **hsum,
		 EigArg *eig_arg, 
		 CnvFrmType cnv_frm = CNV_FRM_YES);
        // It finds the eigenvectors and eigenvalues of A where
        // A is the fermion matrix (Dirac operator). The solution
	// uses Ritz minimization. eig_arg is the 
        // structure that contains all the control parameters, f_eigenv
        // are the fermion field source vectors which should be
        // defined initially, lambda are the eigenvalues returned 
        // on solution. f_eigenv is defined on the whole lattice.
	// The function returns the total number of Ritz iterations.

    Float SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
		 Float mass, DagType dag);
	// It sets the pseudofermion field phi from frm1, frm2.

    ForceArg EvolveMomFforce(Matrix *mom, Vector *frm, 
				 Float mass, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the fermion force.

    ForceArg EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
			  Float mass, Float step_size);
        // It evolve the canonical momentum mom  by step_size
        // using the bosonic quotient force.

    ForceArg RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
			      int isz, Float *alpha, Float mass, Float dt,
			      Vector **sol_d, ForceMeasure measure);

    Float BhamiltonNode(Vector *boson, Float mass);
    // The boson Hamiltonian of the node sublattice
    
};


//------------------------------------------------------------------
//! A class implementing  domain wall fermions.
/*!
  \ingroup factions
*/
//------------------------------------------------------------------
class FdwfBase : public virtual FwilsonTypes
{
 private:
    const char *cname;    // Class name.
    
 public:

    FdwfBase(void);

    virtual ~FdwfBase(void);

    FclassType Fclass() const;
        // It returns the type of fermion class

    int FsiteOffsetChkb(const int *x) const;
        // Sets the offsets for the fermion fields on a 
        // checkerboard. The fermion field storage order
        // is not the canonical one but it is particular
        // to the Dwf fermion type. x[i] is the 
        // ith coordinate where i = {0,1,2,3,4} = {x,y,z,t,s}.

#if 0
    int FsiteOffset(const int *x) const;
        // Sets the offsets for the fermion fields on a 
        // checkerboard. The fermion field storage order
        // is the canonical one. X[I] is the
        // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
#endif

    int FsiteSize() const;
        // Returns the number of fermion field 
        // components (including real/imaginary) on a
        // site of the 4-D lattice.

    int FchkbEvl() const;
        // Returns 0 => If no checkerboard is used for the evolution
        //      or the CG that inverts the evolution matrix.

    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm = CNV_FRM_YES);
        // It calculates f_out where A * f_out = f_in and
        // A is the preconditioned fermion matrix that appears
        // in the HMC evolution (even/odd preconditioning 
        // of [Dirac^dag Dirac]). The inversion is done
	// with the conjugate gradient. cg_arg is the structure
        // that contains all the control parameters, f_in is the
        // fermion field source vector, f_out should be set to be
        // the initial guess and on return is the solution.
	// f_in and f_out are defined on a checkerboard.
        // If true_res !=0 the value of the true residual is returned
        // in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|
	// The function returns the total number of CG iterations.
    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   CnvFrmType cnv_frm = CNV_FRM_YES);

    int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
		    int Nshift, int isz, CgArg **cg_arg, 
		    CnvFrmType cnv_frm, MultiShiftSolveType type, Float *alpha,
		    Vector **f_out_d);

    void FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
		     Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm);
    
    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		Float *true_res,
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
        // It calculates f_out where A * f_out = f_in and
        // A is the fermion matrix (Dirac operator). The inversion
	// is done with the conjugate gradient. cg_arg is the 
        // structure that contains all the control parameters, f_in 
        // is the fermion field source vector, f_out should be set 
        // to be the initial guess and on return is the solution.
	// f_in and f_out are defined on the whole lattice.
        // If true_res !=0 the value of the true residual is returned
        // in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|
        // cnv_frm is used to specify if f_in should be converted 
        // from canonical to fermion order and f_out from fermion 
        // to canonical. 
        // prs_f_in is used to specify if the source
        // f_in should be preserved or not. If not the memory usage
        // is less by half the size of a fermion vector.
	// The function returns the total number of CG iterations.

    int FmatInvMobius(Vector * f_out,
                      Vector * f_in,
                      CgArg * cg_arg_dwf,
                      MdwfArg * mdwf_arg,
                      Float * true_res,
                      CnvFrmType cnv_frm,
                      PreserveType prs_f_in);
    // FmatInvMobius: same as FmatInv, except that we use mobius DWF
    // formalism to speed up the CG inversion (via constructing initial guess).
    // ======================================================================
    // n_restart: How many restarts we perform
      
    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
	
    int eig_FmatInv(Vector **V, const int vec_len, Float *M, const int nev, const int m, float **U, Rcomplex *invH, const int def_len, const Float *restart,const int restart_len,
			Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		Float *true_res,
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
	
    void Ffour2five(Vector *five, Vector *four, int s_u, int s_l, int Ncb=2);
    //!< Transforms a 4-dimensional fermion field into a 5-dimensional field.
    /* The 5d field is zero */
    // The 5d field is zero
    // except for the upper two components (right chirality)
    // at s = s_u which are equal to the ones of the 4d field
    // and the lower two components (left chirality) 
    // at s_l, which are equal to the ones of the 4d field
    // For spread-out DWF s_u, s_l refer to the global
    // s coordinate i.e. their range is from 
    // 0 to [GJP.Snodes() * GJP.SnodeSites() - 1]

    void Ffive2four(Vector *four, Vector *five, int s_u, int s_l, int Ncb=2);
    //!< Transforms a 5-dimensional fermion field into a 4-dimensional field.
    //The 4d field has
        // the upper two components (right chirality) equal to the
        // ones of the 5d field at s = s_u and the lower two 
        // components (left chirality) equal to the
        // ones of the 5d field at s = s_l, where s is the 
        // coordinate in the 5th direction.
        // For spread-out DWF s_u, s_l refer to the global
        // s coordinate i.e. their range is from 
        // 0 to [GJP.Snodes() * GJP.SnodeSites() - 1]
        // The same 4D field is generarted in all s node slices.

    void Fsolfour2five(Vector *sol_5d, Vector *sol_4d, Vector *src_5d, CgArg *cg_arg);
    // Recover the 5D solution from the 4D solution, without solve the equation again.

    int FeigSolv(Vector **f_eigenv, Float *lambda,
		 Float *chirality, int *valid_eig,
		 Float **hsum,
		 EigArg *eig_arg, 
		 CnvFrmType cnv_frm = CNV_FRM_YES);
        // It finds the eigenvectors and eigenvalues of A where
        // A is the fermion matrix (Dirac operator). The solution
	// uses Ritz minimization. eig_arg is the 
        // structure that contains all the control parameters, f_eigenv
        // are the fermion field source vectors which should be
        // defined initially, lambda are the eigenvalues returned 
        // on solution. f_eigenv is defined on the whole lattice.
	// The function returns the total number of Ritz iterations.
    int FeigSolv(Vector **f_eigenv, Float *lambda,
		 LanczosArg *eig_arg, 
		 CnvFrmType cnv_frm = CNV_FRM_YES);
        // It finds the eigenvectors and eigenvalues of A where
        // A is the fermion matrix (Dirac operator). The solution
	// uses implicitly started Lanczos. eig_arg is the 
        // structure that contains all the control parameters, f_eigenv
        // are the fermion field source vectors which should be
        // defined initially, lambda are the eigenvalues returned 
        // on solution. f_eigenv is defined on the whole lattice.
	// The function returns the total number of iterations.

    Float SetPhi(Vector *phi, Vector *frm1, Vector *frm2,	       
		 Float mass, DagType dag);
	// It sets the pseudofermion field phi from frm1, frm2.
	
    ForceArg EvolveMomFforce(Matrix *mom, Vector *frm, 
			  Float mass, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the fermion force.

    ForceArg EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
			  Float mass, Float step_size);
        // It evolve the canonical momentum mom  by step_size
        // using the bosonic quotient force.

    ForceArg EvolveMomFforceInt(Matrix *mom, Vector *v1, Vector *v2,
			  Float mass, Float step_size);

    ForceArg RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
			      int isz, Float *alpha, Float mass, Float dt,
			      Vector **sol_d, ForceMeasure measure);

    Float FhamiltonNode( Vector *phi,  Vector *chi) ;
        // The fermion Hamiltonian of the node sublattice.
        // chi must be the solution of Cg with source phi.	       

    void Fconvert(Vector *f_field,
			  StrOrdType to,
		  StrOrdType from, int cb=2);
        // Convert fermion field f_field from -> to

    void SpinProject(Vector * out, Vector *in, int s_size, int type);
    //--------------------------------------------------------------------
    // void SpinProject():
    //
    // Does a spin projection along s direction, specifically (the
    // matrices appear below apply in s direction, the example are given
    // for Ls=5):
    // 
    // for type == 0:
    //
    //       [ P- P+ 0  0  0  ]
    //       [ 0  P- P+ 0  0  ]
    // out = [ 0  0  P- P+ 0  ] in
    //       [ 0  0  0  P- P+ ]
    //       [ P+ 0  0  0  P- ]
    //
    // for type == 1:
    //
    //       [ P- 0  0  0  P+ ]
    //       [ P+ P- 0  0  0  ]
    // out = [ 0  P+ P- 0  0  ] in
    //       [ 0  0  P+ P- 0  ]
    //       [ 0  0  0  P+ P- ]
    //
    // in and out are assumed to be in CANONICAL storage order.
    // Since this function may be used by Mobius fermions, the
    // parameter s_size specifies the length in s direction for both
    // in and out vectors.
    //--------------------------------------------------------------------

    Float BhamiltonNode(Vector *boson, Float mass);
        // The boson Hamiltonian of the node sublattice

    void Freflex (Vector *out, Vector *in);
    //!< Does something really cool.
       // Reflexion in s operator, needed for the hermitian version 
       // of the dirac operator in the Ritz solver.
};

//------------------------------------------------------------------
//! A class implementing  domain wall fermions.
/*!
  This adds nothing to the base class.
  \ingroup factions
*/
//------------------------------------------------------------------
class Fdwf : public FdwfBase {
 private:
    const char *cname;    // Class name.
    
 public:

    Fdwf(void);
    ~Fdwf(void);
    ForceArg EvolveMomFforce(Matrix *mom, Vector *frm, 
				 Float mass, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the fermion force.
};

class Fmdwf : public virtual Lattice {
 private:
  char * cname;

 public:
  Fmdwf(void);
  ~Fmdwf(void);

  // added for FeigSolv. Since we need a CgArg anyway, and we can't change FeigSolv.
  //  Fmdwf(const CgArg * cg_arg);

  FclassType Fclass() const;
  // It returns the type of fermion class
  
  //! Multiplication of a lattice spin-colour vector by gamma_5.
  void Gamma5(Vector *v_out, Vector *v_in, int num_sites);

  int FsiteOffsetChkb(const int *x) const;
  // Sets the offsets for the fermion fields on a 
  // checkerboard. The fermion field storage order
  // is not the canonical one but it is particular
  // to the Dwf fermion type. x[i] is the 
  // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
  
  int FsiteOffset(const int *x) const;
  // Sets the offsets for the fermion fields on a 
  // checkerboard. The fermion field storage order
  // is the canonical one. X[I] is the
  // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
  
  int FsiteSize() const;
  // Returns the number of fermion field 
  // components (including real/imaginary) on a
  // site of the 4-D lattice.
  
  int FchkbEvl() const;
  // Returns 0 => If no checkerboard is used for the evolution
  //      or the CG that inverts the evolution matrix.
  
  int FmatEvlInv(Vector *f_out, Vector *f_in, 
                 CgArg *cg_arg, 
                 Float *true_res,
                 CnvFrmType cnv_frm = CNV_FRM_YES);
  // It calculates f_out where A * f_out = f_in and
  // A is the preconditioned fermion matrix that appears
  // in the HMC evolution (even/odd preconditioning 
  // of [Dirac^dag Dirac]). The inversion is done
  // with the conjugate gradient. cg_arg is the structure
  // that contains all the control parameters, f_in is the
  // fermion field source vector, f_out should be set to be
  // the initial guess and on return is the solution.
  // f_in and f_out are defined on a checkerboard.
  // If true_res !=0 the value of the true residual is returned
  // in true_res.
  // *true_res = |src - MatPcDagMatPc * sol| / |src|
  // The function returns the total number of CG iterations.
  int FmatEvlInv(Vector *f_out, Vector *f_in, 
                 CgArg *cg_arg, 
                 CnvFrmType cnv_frm = CNV_FRM_YES);
  
  int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
                  int Nshift, int isz, CgArg **cg_arg, 
                  CnvFrmType cnv_frm, MultiShiftSolveType type, Float *alpha,
                  Vector **f_out_d);
  
  void FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
                  Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm);
  
  int FmatInv(Vector *f_out, Vector *f_in, 
              CgArg *cg_arg, 
              Float *true_res,
              CnvFrmType cnv_frm = CNV_FRM_YES,
              PreserveType prs_f_in = PRESERVE_YES);
  // It calculates f_out where A * f_out = f_in and
  // A is the fermion matrix (Dirac operator). The inversion
  // is done with the conjugate gradient. cg_arg is the 
  // structure that contains all the control parameters, f_in 
  // is the fermion field source vector, f_out should be set 
  // to be the initial guess and on return is the solution.
  // f_in and f_out are defined on the whole lattice.
  // If true_res !=0 the value of the true residual is returned
  // in true_res.
  // *true_res = |src - MatPcDagMatPc * sol| / |src|
  // cnv_frm is used to specify if f_in should be converted 
  // from canonical to fermion order and f_out from fermion 
  // to canonical. 
  // prs_f_in is used to specify if the source
  // f_in should be preserved or not. If not the memory usage
  // is less by half the size of a fermion vector.
  // The function returns the total number of CG iterations.
  int FmatInv(Vector *f_out, Vector *f_in, 
              CgArg *cg_arg, 
              CnvFrmType cnv_frm = CNV_FRM_YES,
              PreserveType prs_f_in = PRESERVE_YES);
  
  // FmatInvMobius: same as FmatInv, except that we use mobius DWF
  // formalism to speed up the CG inversion (via constructing initial guess).
  // n_restart: How many restarts we perform
  int FmatInvMobius(Vector * f_out,
                    Vector * f_in,
                    MdwfArg * mob_l,
                    MdwfArg * mob_s,
                    Float * true_res,
                    CnvFrmType cnv_frm,
                    PreserveType prs_f_in,
                    int n_restart, Float rsd_vec[]);
  
  void Ffour2five(Vector *five, Vector *four, int s_u, int s_l, int Ncb=2);
  //!< Transforms a 4-dimensional fermion field into a 5-dimensional field.
  /* The 5d field is zero */
  // The 5d field is zero
  // except for the upper two components (right chirality)
  // at s = s_u which are equal to the ones of the 4d field
  // and the lower two components (left chirality) 
  // at s_l, which are equal to the ones of the 4d field
  // For spread-out DWF s_u, s_l refer to the global
  // s coordinate i.e. their range is from 
  // 0 to [GJP.Snodes() * GJP.SnodeSites() - 1]
  
  void Ffive2four(Vector *four, Vector *five, int s_u, int s_l, int Ncb=2);
  //!< Transforms a 5-dimensional fermion field into a 4-dimensional field.
  //The 4d field has
  // the upper two components (right chirality) equal to the
  // ones of the 5d field at s = s_u and the lower two 
  // components (left chirality) equal to the
  // ones of the 5d field at s = s_l, where s is the 
  // coordinate in the 5th direction.
  // For spread-out DWF s_u, s_l refer to the global
  // s coordinate i.e. their range is from 
  // 0 to [GJP.Snodes() * GJP.SnodeSites() - 1]
  // The same 4D field is generarted in all s node slices.
  
  int FeigSolv(Vector **f_eigenv, Float *lambda,
               Float *chirality, int *valid_eig,
               Float **hsum,
               EigArg *eig_arg, 
               CnvFrmType cnv_frm = CNV_FRM_YES);
  // It finds the eigenvectors and eigenvalues of A where
  // A is the fermion matrix (Dirac operator). The solution
  // uses Ritz minimization. eig_arg is the 
  // structure that contains all the control parameters, f_eigenv
  // are the fermion field source vectors which should be
  // defined initially, lambda are the eigenvalues returned 
  // on solution. f_eigenv is defined on the whole lattice.
  // The function returns the total number of Ritz iterations.
  
  Float SetPhi(Vector *phi, Vector *frm1, Vector *frm2,	       
               Float mass, DagType dag);
  // It sets the pseudofermion field phi from frm1, frm2.
  
  ForceArg EvolveMomFforce(Matrix *mom, Vector *frm, 
                           Float mass, Float step_size);
  // It evolves the canonical momentum mom by step_size
  // using the fermion force.
  
  ForceArg EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
                           Float mass, Float step_size);
  // It evolve the canonical momentum mom  by step_size
  // using the bosonic quotient force.
  
  ForceArg RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
                                int isz, Float *alpha, Float mass, Float dt,
                                Vector **sol_d, ForceMeasure measure);
  
  Float FhamiltonNode( Vector *phi,  Vector *chi) ;
  // The fermion Hamiltonian of the node sublattice.
  // chi must be the solution of Cg with source phi.	       
  
  void Fconvert(Vector *f_field,
                StrOrdType to,
                StrOrdType from);
  // Convert fermion field f_field from -> to
  
  Float BhamiltonNode(Vector *boson, Float mass);
  // The boson Hamiltonian of the node sublattice
  
  void SpinProject(Vector * out, Vector *in, int s_size, int type);
  //--------------------------------------------------------------------
  // void SpinProject():
  //
  // Does a spin projection along s direction, specifically (the
  // matrices appear below apply in s direction, the example are given
  // for Ls=5):
  // 
  // for type == 0:
  //
  //       [ P- P+ 0  0  0  ]
  //       [ 0  P- P+ 0  0  ]
  // out = [ 0  0  P- P+ 0  ] in
  //       [ 0  0  0  P- P+ ]
  //       [ P+ 0  0  0  P- ]
  //
  // for type == 1:
  //
  //       [ P- 0  0  0  P+ ]
  //       [ P+ P- 0  0  0  ]
  // out = [ 0  P+ P- 0  0  ] in
  //       [ 0  0  P+ P- 0  ]
  //       [ 0  0  0  P+ P- ]
  //
  // in and out are assumed to be in CANONICAL storage order.
  //--------------------------------------------------------------------

  void Freflex (Vector *out, Vector *in);
  //!< Does something really cool.
  // Reflexion in s operator, needed for the hermitian version 
  // of the dirac operator in the Ritz solver.

  int SpinComponents() const;

  int ExactFlavors() const;
    
  //!< Method to ensure bosonic force works (does nothing for Wilson
  //!< theories.
  void BforceVector(Vector *in, CgArg *cg_arg);
};

class Fmobius : public FdwfBase {
 private:
    char *cname;    // Class name.
    
 public:

    Fmobius(void);
    ~Fmobius(void);

    FclassType Fclass(void) const;

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		Float *true_res,
		CnvFrmType cnv_frm,
		PreserveType prs_f_in);

    int FeigSolv(Vector **f_eigenv, Float *lambda,
		 Float *chirality, int *valid_eig,
		 Float **hsum,
		 EigArg *eig_arg, 
		 CnvFrmType cnv_frm);

    int FeigSolv(Vector **f_eigenv, Float *lambda,
		 LanczosArg *eig_arg, 
		 CnvFrmType cnv_frm);
};

CPS_END_NAMESPACE
#endif


