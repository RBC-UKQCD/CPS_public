#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/dirac_op.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: dirac_op.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.5  2001/08/16 12:54:30  anj
//  Some fixes follosin the float-> IFloat change, mostly of the (variable
//  anme) IFloat_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.4  2001/08/16 10:50:29  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:16  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:09  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: dirac_op.h,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/dirac_op.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// dirac_op.h
//
// Header file for the DiracOp base class. Each type of frmion has 
// a derived class. DiracOp is the front end for all Dirac operators 
// associated with a specific kind of  fermions.
//
//------------------------------------------------------------------


#ifndef INCLUDED_DIRAC_OP_H
#define INCLUDED_DIRAC_OP_H

CPS_END_NAMESPACE
#include<util/lattice.h>
#include<util/vector.h>
#include<alg/cg_arg.h>
#include<alg/eig_arg.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
//
// DiracOp base class.
//
//------------------------------------------------------------------
class DiracOp
{
 private:
  char *cname;                     // Class name.

  static int scope_lock;           // lock that forbids more than
                                   // one DiracOp object to be on
                                   // scope at any time.

 protected:
  Lattice& lat;                    // Lattice object.
  Vector *f_out;                   // Output fermion field pointer.
  Vector *f_in;                    // Input fermion field pointer.
  CgArg *dirac_arg;                // Argument structure
  CnvFrmType cnv_frm;              // Fermion conversion flag
  Matrix *gauge_field;             // pointer to the gauge field
  
  virtual void DiracOpGlbSum(Float *float_p);					
     // The global sum used by InvCg. At the base class
     // level it is just the usual global sum. It is overloaded
     // by any definitions by the derived classes. This is needed
     // for example by DiracOpDwf where the global sum has different
     // meaning for s_nodes = 1 and s_nodes > 1.


 public:
  DiracOp(Lattice& latt,           // Lattice object.
	  Vector *f_field_out,     // Output fermion field ptr.
	  Vector *f_field_in,      // Input fermion field ptr.
	  CgArg *arg,              // Argument structure
	  CnvFrmType cnv_frm_flg); // Fermion conversion flag

  virtual ~DiracOp();


  int InvCg(Vector *out,
	    Vector *in,
	    Float src_norm_sq,
	    Float *true_res);
     // The Conjugate Gradient inverter. 
     // Source is *in, initial guess and solution is *out.
     // src_norm_sq is the norm of the source squared 
     // (dot product of source * source). If src_norm_sq=0 
     // it is calculated inside InvCg. 
     // If true_res !=0 the value of the true residual is returned
     // in true_res.
     // *true_res = |src - MatPcDagMatPc * sol| / |src|

  int InvCg(Vector *out, Vector *in, Float src_norm_sq);
      // Same as original but with true_res=0

  int InvCg(Vector *out, Vector *in, Float *true_res);
     // Same as original but with src_norm_sq=0.0

  int InvCg(Vector *out, Vector *in);
     // Same as original but with src_norm_sq=0.0, true_res=0

  int InvCg(Float src_norm_sq, Float *true_res);
     // Same as original but with 
     // in=f_in, out=f_out

  int InvCg(Float src_norm_sq);
     // Same as original but with 
     // in=f_in, out=f_out, true_res=0

  int InvCg(Float *true_res);
     // Same as original but with 
     // in=f_in, out=f_out, src_norm_sq=0.0

  int InvCg(void);
     // Same as original but with 
     // in=f_in, out=f_out, src_norm_sq=0.0, true_res=0

  int Ritz(Vector **psi_all, int N_eig, Float &lambda,
	   Float RsdR_a, Float RsdR_r, Float Rsdlam, Float Cutl_zero,
	   int n_renorm, int Kalk_Sim, int N_min, int N_max, Float Cv_fact,
	   int MaxCG, int ProjApsiP);
     // Ritz minimizer for eigenvectors.
     // The matrix that is used is RitzMat.
     // RitzMat must be hermitian and positive def.

  int Jacobi(Vector **psi, int N_eig, Float lambda[], 
	     Complex off_diag[], 
	     Float Toler, int N_max);
     // Jacobi diagonalizer of a single site matrix.

// Pure virtual functions
//------------------------------------------------------------------
  virtual void DiracArg(CgArg *arg) = 0;
     // It sets the dirac_arg pointer to arg and initializes
     // the relevant parameters (kappa, m^2, ...).

  virtual void MatPcDagMatPc(Vector *out, Vector *in, Float *dot_prd=0) = 0;
     // MatPcDagMatPc is the fermion matrix that appears in the HMC 
     // evolution. It is a Hermitian matrix where M is
     // the preconditioned (if relevant) Dirac Operator matrix.        
     // MatPcDagMatPc connects only even-->even or odd-->odd sites.
     // The in, out fields are defined on a checkerboard.
     // If dot_prd is not 0 then the dot product (on node)
     // <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.
 
  virtual void Dslash(Vector *out, 
		      Vector *in,
		      ChkbType cb, 
		      DagType dag)=0;
     // Dslash is the derivative part of the fermion matrix. 
     // Dslash conects only odd-->even or even-->odd sites.
     // The in, out fields are defined on a checkerboard.
     // cb refers to the checkerboard of the in field.

  virtual int MatInv(Vector *out, 
		     Vector *in, 
		     Float *true_res, 
		     PreserveType prs_in = PRESERVE_YES) = 0;
     // The inverse of the unconditioned (if relevant) Dirac Operator 
     // using Conjugate gradient. source is *in, initial
     // guess and solution is *out.
     // If true_res !=0 the value of the true residual is returned
     // in true_res.
     // *true_res = |src - MatPcDagMatPc * sol| / |src|
     // prs_in is used to specify if the source
     // in should be preserved or not. If not the memory usage
     // is less by the size of one fermion vector or by the size
     // of one checkerboard fermion vector (half a fermion vector).
     // For staggered fermions in is preserved regardles of
     // the value of prs_in. 
     // The function returns the total number of CG iterations.

  virtual int MatInv(Vector *out, Vector *in, 
		     PreserveType prs_in = PRESERVE_YES) = 0;
     // Same as original but true_res=0.

  virtual int MatInv(Float *true_res, 
		     PreserveType prs_in = PRESERVE_YES) = 0;
     // Same as original but in = f_in and out = f_out.

  virtual int MatInv(PreserveType prs_in = PRESERVE_YES) = 0;
     // Same as original but in = f_in, out = f_out, true_res=0.

  virtual int RitzEig(Vector **eigenv, Float lambda[], int valid_eig[], 
		      EigArg *eig_arg) = 0;
     // The eigenvector solver using Ritz minimization and Jacobi reconstruction
     // (using the Kalkreuter-Simma algorithm) of the eigenvectors of RitzEigMat.
     // The matrix that is used is RitzEigMat.
     // RitzEigMat is not necessarily to be positive def. (but must be hermitian).
     // For RitzEigJac, RitzMat = RitzEigMat*RitzEigMat when using KS algorithm.
     // Can also use non KS algorithm to return eigenvalues only of 
     // RitzEigMat = RitzMat.

  virtual void RitzEigMat(Vector *out, Vector *in) = 0;
     // RitzEigMat is the fermion matrix used in RitzEig
     // RitzEigMat works on the full lattice or half lattice
     // The in, out fields are defined on the full or half lattice.

  virtual void RitzMat(Vector *out, Vector *in) = 0;
     // RitzMat is the fermion matrix used in Ritz
     // RitzMat works on the full lattice or half lattice
     // The in, out fields are defined on the full or half lattice.

  virtual int RitzLatSize() = 0;
     // RitzLatSize returns the number of checkerboards of a fermion.
     // It uses the RitzMatType flag to determine the operator
     // to use and relevant checkerboard sizes.
};


//------------------------------------------------------------------
//
// DiracOpStagTypes class.
//
// Is derived from DiracOp and is relevant to
// all DiracOp classes with Staggered type fermions.
// These classes are derived from DiracOpStagTypes
//
//------------------------------------------------------------------
class DiracOpStagTypes : public DiracOp
{
 private:
  char *cname;    // Class name.

 public:
  DiracOpStagTypes(Lattice& latt,            // Lattice object.
		   Vector *f_field_out,      // Output fermion field ptr.
		   Vector *f_field_in,       // Input fermion field ptr.
		   CgArg *arg,               // Argument structure
		   CnvFrmType cnv_frm_flg);  // Fermion conversion flag

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

  virtual int RitzLatSize();
     // RitzLatSize returns the size of a fermion on a node
     // It uses the RitzMatType flag to determine the operator
     // to use and relevant checkerboard sizes.
};


//------------------------------------------------------------------
//
// DiracOpStag is derived from DiracOpStagTypes and is the front
// end for all Dirac operators associated with Staggered fermions.
//
//------------------------------------------------------------------
class DiracOpStag : public DiracOpStagTypes
{
 private:
  char *cname;         // Class name.

  Float mass_rs;       // rescaled mass
  Float mass_sq;       // = mass^2

  int f_size_cb;       //The node checkerbrd. size of the ferm. field

  Vector *frm_tmp;     // Temporary fermion field


 public:
  DiracOpStag(Lattice& latt,            // Lattice object.
	      Vector *f_field_out,      // Output fermion field ptr.
	      Vector *f_field_in,       // Input fermion field ptr.
	      CgArg *arg,               // Argument structure
	      CnvFrmType cnv_frm_flg);  // Fermion conversion flag

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

  void RitzMat(Vector *out, Vector *in);
     // RitzMat is the fermion matrix used in Ritz
     // RitzMat works on the full lattice or half lattice
     // The in, out fields are defined on the full or half lattice.
};


//------------------------------------------------------------------
//
// DiracOpWilsonTypes class.
//
// Is derived from DiracOp and is relevant to
// all DiracOp classes with Wilson type fermions 
// (e.g DiracOpWilson, DiracOpClover, DiracOpDwf, ...). 
//  These classes are derived from DiracOpWilsonTypes
//
//------------------------------------------------------------------
class DiracOpWilsonTypes : public DiracOp
{
 private:
  char *cname;    // Class name.

 protected:
  void MultGamma(Vector *out, const Vector *in, int gamma_num, int nodevol);
     // Multiply vector on left by 
     // gamma_1^n1 * gamma_2^n2 * gamma_3^n3 * gamma_4^n4
     // Where  (n1,n2,n3,n4) = bit wise decomposition of gamma_num
     // Vector has nodevol number of sites to run over
     // out = gamma_1^n1*gamma_2^n2*gamma_3^n3*gamma_4^n4 * in

 public:
  DiracOpWilsonTypes(Lattice& latt,            // Lattice object.
		     Vector *f_field_out,      // Output fermion field ptr.
		     Vector *f_field_in,       // Input fermion field ptr.
		     CgArg *arg,               // Argument structure
		     CnvFrmType cnv_frm_flg);  // Fermion conversion flag

  virtual ~DiracOpWilsonTypes();

  virtual void MatHerm(Vector *out, Vector *in) = 0;
     // MatHerm is the hermitian version of Mat.
     // MatHerm works on the full lattice.
     // The in, out fields are defined on the ful.

  virtual void Mat(Vector *out, Vector *in) = 0;
     // Mat is the unpreconditioned fermion matrix.  
     // Mat works on the full lattice
     // The in, out fields are defined on the full lattice.

  virtual void MatDag(Vector *out, Vector *in) = 0;
     // MatDag is the dagger of the unpreconditioned fermion matrix. 
     // MatDag works on the full lattice
     // The in, out fields are defined on the full lattice.

  virtual void MatDagMat(Vector *out, Vector *in);
     // MatDagMat is the MatDag*Mat.
     // MatDagMat works on the full lattice.
     // The in, out fields are defined on the full.

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
     // RitzMat is the fermion matrix used in Ritz
     // RitzMat works on the full lattice or half lattice
     // The in, out fields are defined on the full or half lattice.

  virtual int RitzLatSize();
     // RitzLatSize returns the size of a fermion on a node
     // It uses the RitzMatType flag to determine the operator
     // to use and relevant checkerboard sizes.
};



//------------------------------------------------------------------
//
// DiracOpWilson is derived from DiracOp and is the frot end for 
// all Dirac operators associated with Wilson fermions.
//
//------------------------------------------------------------------
class DiracOpWilson : public DiracOpWilsonTypes
{
 private:
  char *cname;    // Class name.

  void *wilson_lib_arg;  // pointer to an argument structure related
                         // to the wilson library.

  Float kappa;           // Wilson kappa = 1 /[2 (4 +mass) ]


 public:
  DiracOpWilson(Lattice& latt,            // Lattice object.
		Vector *f_field_out,      // Output fermion field ptr.
		Vector *f_field_in,       // Input fermion field ptr.
		CgArg *arg,               // Argument structure
		CnvFrmType cnv_frm_flg);  // Fermion conversion flag

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

  void MatPc(Vector *out, Vector *in);
     // MatPc is the preconditioned fermion matrix.  
     // MatPc connects only odd-->odd sites.
     // The in, out fields are defined on the odd checkerboard.

  void MatPcDag(Vector *out, Vector *in);
     // MatPcDag is the dagger of the preconditioned fermion matrix. 
     // MatPcDag connects only odd-->odd sites.
     // The in, out fields are defined on the odd checkerboard.

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
     // MatHerm is the hermitian version of Mat.
     // MatHerm works on the full lattice.
     // The in, out fields are defined on the ful.

  void Mat(Vector *out, Vector *in);
     // Mat is the unpreconditioned fermion matrix.  
     // Mat works on the full lattice
     // The in, out fields are defined on the full lattice.

  void MatDag(Vector *out, Vector *in);
     // MatDag is the dagger of the unpreconditioned fermion matrix. 
     // MatDag works on the full lattice
     // The in, out fields are defined on the full lattice.

  void CalcHmdForceVecs(Vector *chi) ;
    // GRF
    // chi is the solution to MatPcInv.  The user passes two full size
    // CANONICAL fermion vectors with conversion enabled to the
    // constructor.  Using chi, the function fills these vectors;
    // the result may be used to compute the HMD fermion force.
};


//------------------------------------------------------------------
//
// DiracOpClover is derived from DiracOp and is the frot end for 
// all Dirac operators associated with Clover fermions.
//
//------------------------------------------------------------------

//------------------------------------------------------------------
// Forward struct declaration for class DiracOpClover
//------------------------------------------------------------------
struct Clover;
class Fclover;

class DiracOpClover : public DiracOpWilsonTypes
{
 private:
  char *cname;            // Class name.
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
		CnvFrmType cnv_frm_flg);  // Fermion conversion flag

  virtual ~DiracOpClover();

  // Auxiliary member functions
  //----------------------------------------------------------------------
  const Matrix& GetLink(const int *site, int dir) const;
  // returns a reference to a link in Wilson order.
  // The site[] is [x,y,z,t] and could be off-node.
  // If off-node, the link is stored in a static buffer and will be
  // overwritten by the next off-node link, so BE CAREFUL.

  void SiteFuv(Matrix &Fuv, const int *site, int mu, int nu) const;
  // calculate Fuv, on return Fuv follows the matrix convention
  // of IFloat[row][col][2], instead of IFloat[col][row][2] used in
  // Wilson order.

  void SiteCloverMat(const int *site, IFloat *mat_72_IFloats) const;
  // calculate the local 12x12 clover matrix, and stored it as
  // two blocks of 6x6 compressed (that is, only the lower diagonal
  // part is stored) hermitian matrices.

  void CloverMatChkb(ChkbType chkb, int inverse = 0) const;
  // calculate the clover matrices for all even, or all odd sites, and
  // invert them if specified.

  void MatPcDagOrNot(Vector *out, const Vector *in, int dag) const;
  // out = MatPcDag in   OR  out = MatPc in

  // FOR THE PURPOSE OF DEBUGGING ONLY
  //----------------------------------------------------------------------
  void MatDagOrNotDbg(Vector *out, const Vector *in, int dag, int direct) const;
  // out = MatDag in    OR  out = Mat in

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

  void MatPc(Vector *out, Vector *in) {
    MatPcDagOrNot(out, in, 0);
  }
     // MatPc(Vector *out, Vector *in) :
     // MatPc is the preconditioned fermion matrix.  
     // MatPc connects only odd-->odd sites.
     // The in, out fields are defined on the odd checkerboard.

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
     // Same as original but true_res=0.

  int MatEvlInv(Float *true_res);
     // Same as original but in = f_in and out = f_out.

  int MatEvlInv(void);
     // Same as original but in = in, f_out = f_out, true_res=0.

  void MatHerm(Vector *out, Vector *in);
     // MatHerm is the hermitian version of Mat.
     // MatHerm works on the full lattice.
     // The in, out fields are defined on the ful.

  void Mat(Vector *out, Vector *in);
     // Mat is the unpreconditioned fermion matrix.  
     // Mat works on the full lattice
     // The in, out fields are defined on the full lattice.

  void MatDag(Vector *out, Vector *in);
     // MatDag is the dagger of the unpreconditioned fermion matrix. 
     // MatDag works on the full lattice
     // The in, out fields are defined on the full lattice.

  void CalcHmdForceVecs(Vector *chi) ;
    // Lingling changed the prvious defintion of the fermion fields by
    // some constant factors. chi is the solution to MatPcInv on the  
    // odd sites.  The user passes two full size  
    // CANONICAL fermion vectors with conversion enabled to the
    // constructor.  Using chi, the function fills these vectors;
    // the result may be used to compute the HMD fermion force.

};


//------------------------------------------------------------------
//
// DiracOpDwf is derived from DiracOp and is the frot end for 
// all Dirac operators associated with Dwf fermions.
//
//------------------------------------------------------------------
class DiracOpDwf : public DiracOpWilsonTypes
{
 private:
  char *cname;    // Class name.

  void *dwf_lib_arg;     // pointer to an argument structure related
                         // to the dwf library.

  Float mass;            // Dwf mass (couples left-right components)


 protected:
  void DiracOpGlbSum(Float *float_p);					
     // The global sum used by InvCg. If s_nodes = 1
     // it is the usual global sum. If s_nodes > 1 it
     // is the 5-dimensional globals sum glb_sum_five.


 public:
  DiracOpDwf(Lattice& latt,            // Lattice object.
	     Vector *f_field_out,      // Output fermion field ptr.
	     Vector *f_field_in,       // Input fermion field ptr.
	     CgArg *arg,               // Argument structure
	     CnvFrmType cnv_frm_flg);  // Fermion conversion flag

  virtual ~DiracOpDwf();

  void DiracArg(CgArg *arg);
     // It sets the dirac_arg pointer to arg and initializes
     // kappa

  void MatPcDagMatPc(Vector *out, Vector *in, Float *dot_prd=0);
     // MatPcDagMatPc is the fermion matrix that appears in the HMC 
     // evolution. It is a Hermitian matrix.
     // The in, out fields are defined on the checkerboard lattice.
     // If dot_prd is not 0 then the dot product (on node)
     // <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.

  void Dslash(Vector *out, 
		      Vector *in,
		      ChkbType cb, 
		      DagType dag);
     // Dslash is the derivative part of the fermion matrix. 
     // The in, out fields are defined on the checkerboard lattice
     // cb = 0/1 <--> even/odd checkerboard of in field.
     // dag = 0/1 <--> Dslash/Dslash^dagger is calculated.

  void MatPc(Vector *out, Vector *in);
     // MatPc is the fermion matrix.  
     // The in, out fields are defined on the checkerboard lattice.

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

  void Mat(Vector *out, Vector *in);
     // Mat is the unpreconditioned fermion matrix.  
     // Mat works on the full lattice
     // The in, out fields are defined on the full lattice.

  void MatDag(Vector *out, Vector *in);
     // MatDag is the dagger of the unpreconditioned fermion matrix. 
     // MatDag works on the full lattice
     // The in, out fields are defined on the full lattice.

  void CalcHmdForceVecs(Vector *chi) ;
    // GRF
    // chi is the solution to MatPcInv.  The user passes two full size
    // CANONICAL fermion vectors with conversion enabled to the
    // constructor.  Using chi, the function fills these vectors;
    // the result may be used to compute the HMD fermion force.

  void Reflex(Vector *out, Vector *in);
    // Reflexion in s operator, needed for the hermitian version 
    // of the dirac operator in the Ritz solver.
};

#endif



CPS_END_NAMESPACE
