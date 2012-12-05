#include<config.h>
/*!\file
  \brief  Definition of the Dirac operator classes: DiracOp, DiracOpStagTypes.

  $Id: d_op_dwf.h,v 1.3 2012-12-05 16:39:19 chulwoo Exp $
*/

#ifndef INCLUDED_D_OP_DWF_H
#define INCLUDED_D_OP_DWF_H

#include <util/dirac_op.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
//! A class describing the Dirac operator for domain-wall Wilson fermions.
/*!
  See the description of the DiracOpWilsonTypes class for the definition
  of the Wilson fermion matrix.

  The constructor changes the storage order of the gauge field to ::WILSON
  order This change persists throughout the lifetime of the object.

  Only one instance of this class is allowed to be in existence at any time.
*/
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
	     CnvFrmType convert);  // Fermion conversion flag

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

  //! Multiplication by the odd-even preconditioned fermion matrix.
  void MatPc(Vector *out, Vector *in);
     // MatPc is the fermion matrix.  
     // The in, out fields are defined on the checkerboard lattice.

  //! Multiplication by the  hermitian conjugate odd-even preconditioned fermion matrix.
  void MatPcDag(Vector *out, Vector *in);
     // MatPcDag is the dagger of the fermion matrix. 
     // The in, out fields are defined on the checkerboard lattice.

#ifdef USE_BFM
//  int InvCg(Vector *out,
//	    Vector *in,
//	    Float src_norm_sq,
//	    Float *true_res);
  int InvCg(Vector *out, Vector *in, Float *true_res);
  int InvCg(Float *true_res){
   InvCg(f_out,f_in,true_res);
  }
#ifdef USE_BFM_MINV
  int MInvCG(Vector **out, Vector *in, Float in_norm, Float *shift, 
	     int Nshift, int isz, Float *RsdCG, 
	     MultiShiftSolveType type, Float *alpha);
#endif

#endif

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
   int eig_MatInv(Vector **V, const int vec_len, Float *M, const int nev, const int m, float **U, Rcomplex *invH, const int def_len, const Float *restart, const int restart_len,Vector *out, Vector *in, Float *true_res, PreserveType prs_in = PRESERVE_YES);
   //eigCG inversion //by Qi Liu
   //


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
#ifdef USE_QUDA
  int QudaInvert(Vector *out, Vector *in, Float *true_res, int mat_type);
#endif

};
CPS_END_NAMESPACE

#endif
