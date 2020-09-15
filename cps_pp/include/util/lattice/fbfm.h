#ifndef INCLUDED_FBFM_H__
#define INCLUDED_FBFM_H__

#include<config.h>

#ifdef USE_BFM
#include <util/lattice/bfm_evo.h>
#include <util/lattice/bfm_mixed_solver.h>
#include <util/lattice/bfm_mixed_solver_multi.h>
#include <util/lattice/eff_overlap.h>
# endif

#include <util/lattice.h>

CPS_START_NAMESPACE

#ifdef USE_BFM
class Fbfm : public virtual Lattice,public virtual FwilsonTypes {
 public:
  //CK: Modified to allow for multiple sets of arguments for different fermion types
//  static int current_arg_idx;  //current array index of bfmarg within arg array. Switch is performed either manually or automatically in LatticeFactory
    //static bfmarg bfm_arg;

    //CK: Modified to allow for multiple sets of arguments for different fermion types
//    static int current_arg_idx;  //current array index of bfmarg within arg array. Switch is performed either manually or automatically in LatticeFactory
//    static bfmarg bfm_args[2]; //currently setup to allow 2 different choices corresponding to F_CLASS_BFM and F_CLASS_BFM_TYPE2, can be extended in principle
//    static int nthreads[2];
  
    static std::map<Float, bfmarg> arg_map;
    static Float default_key_mass;
    static Float current_key_mass;
	Float key_mass;
//BfmSolver solver;
	static BfmSolver CurrentSolver(){
     return arg_map.at(current_key_mass).solver;
   };
	BfmSolver Solver(Float mass){
     return arg_map.at(mass).solver;
   };

    

    // Lets the user specify what MADWF parameters should be used for each key mass during measurements.
    // Leave empty to use regular CG.
  static bool use_mixed_solver; // set true to use mixed solver for regular matrix inversion. For multi-shift solves this activates the use of the MultiShiftCGcontroller, otherwise they 
  // are just performed in double precision
    
  bfm_evo<double> bd;
  bfm_evo<float> bf;
  static std::map<Float, MADWFParams> madwf_arg_map;
//    bfm_evo<double> kernel;

 private:
  const char *cname;

    bool bfm_initted;
  // These are eigenvectors/eigenvalues obtained from Rudy's Lanczos
  // code. Use them for deflation.
  multi1d<bfm_fermion> *evec;
  multi1d<double> *evald;
  multi1d<float> *evalf;
  int ecnt;

 public:
  Fbfm(void);
  virtual ~Fbfm(void);

  template<typename EVAL_TYPE>
    void set_deflation(multi1d<Fermion_t[2]> *_evec,
                       multi1d<EVAL_TYPE> *_eval,
                       int _ecnt)
    {
      evec = _evec;

      evald = NULL;
      evalf = NULL;
      if(sizeof(EVAL_TYPE) == sizeof(double)) {
	evald = (multi1d<double> *)_eval;
      } else {
	evalf = (multi1d<float> *)_eval;
      }
      ecnt = _ecnt;
    }

  void unset_deflation(void) {
    evec = NULL;
    evald = NULL;
    evalf = NULL;
    ecnt = 0;
  }

  void CalcHmdForceVecsBilinear(Float *v1, Float *v2,
				Vector *phi1, Vector *phi2,
				Float mass);

#if 1
  ForceArg EvolveMomFforceBaseThreaded(Matrix *mom,
				       Vector *phi1,
				       Vector *phi2,
				       Float mass,
				       Float coef);
  // It evolves the canonical Momemtum mom:
  // mom += coef * (phi1^\dag e_i(M) \phi2 + \phi2^\dag e_i(M^\dag) \phi1)
  // note: this function does not exist in the base Lattice class.
  //CK: This function is not used, so I have not modified it for WilsonTM
#endif

  ForceArg EvolveMomFforceBase(Matrix *mom,
			       Vector *phi1,
			       Vector *phi2,
			       Float mass,
			       Float coef);
  // It evolves the canonical Momemtum mom:
  // mom += coef * (phi1^\dag e_i(M) \phi2 + \phi2^\dag e_i(M^\dag) \phi1)
  // note: this function does not exist in the base Lattice class.
  //CK: Version for non-WilsonTm quarks. Will throw an error if used for WilsonTm

    FclassType Fclass()const {
        return F_CLASS_BFM;
    }

  // It returns the type of fermion class
  
  //! Multiplication of a lattice spin-colour vector by gamma_5.
  void Gamma5(Vector *v_out, Vector *v_in, int num_sites);

  int FsiteOffsetChkb(const int *x) const;
  // Sets the offsets for the fermion fields on a 
  // checkerboard. The fermion field storage order
  // is not the canonical one but it is particular
  // to the Dwf fermion type. x[i] is the 
  // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
  
#if 0
  int FsiteOffset(const int *x) const;
  // Sets the offsets for the fermion fields on a 
  // checkerboard. The fermion field storage order
  // is the canonical one. X[I] is the
  // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
#endif
  
    int FsiteSize() const {
	const char* fname = "FsiteSize()";
	if (arg_map.count(current_key_mass) == 0) {
	    ERR.General(cname, fname, "No entry for current key mass %e in arg_map!\n", current_key_mass);
	    return 0;
	} else {
	    int Ls = arg_map.at(current_key_mass).Ls;
	    int ret = 24 * Ls;
	    //printf("FsiteSize() using current_key_mass = %e -> Ls = %d -> site size = %d!\n", current_key_mass, Ls, ret);
	    return ret;
	}
    }
  // Returns the number of fermion field 
  // components (including real/imaginary) on a
  // site of the 4-D lattice.
  
  int FchkbEvl() const {
    if((bd.solver == HtCayleyEOFA) || (bd.solver == HmCayleyEOFA)){ return 0; }
    else{ return 1; }
  }
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
		 CnvFrmType cnv_frm = CNV_FRM_YES)
  {
    return FmatEvlInv(f_out, f_in, cg_arg, NULL, cnv_frm);
  }
  
  int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
		  int Nshift, int isz, CgArg **cg_arg, 
		  CnvFrmType cnv_frm, MultiShiftSolveType type, Float *alpha,
		  Vector **f_out_d);
  
  void FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
		  Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm);
  
  int FmatInv(Vector *f_out, Vector *f_in, 
	      CgArg *cg_arg, 
	      Float *true_res,
	      CnvFrmType cnv_frm ,
	      PreserveType prs_f_in ,int if_dminus);
  int FmatInv(Vector *f_out, Vector *f_in, 
	      CgArg *cg_arg, 
	      Float *true_res,
	      CnvFrmType cnv_frm = CNV_FRM_YES,
	      PreserveType prs_f_in = PRESERVE_YES){
	FmatInv(f_out,f_in, cg_arg,true_res,cnv_frm,prs_f_in,1);
}

  int FmatInvTest(Vector *f_out, Vector *f_in, 
	      CgArg *cg_arg, 
	      Float *true_res,
	      CnvFrmType cnv_frm = CNV_FRM_YES,
	      PreserveType prs_f_in = PRESERVE_YES){
	FmatInv(f_out,f_in, cg_arg,true_res,cnv_frm,prs_f_in,0);
}

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
	      PreserveType prs_f_in = PRESERVE_YES)
  {
    return FmatInv(f_out, f_in, cg_arg, NULL, cnv_frm, prs_f_in);
  }
  
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

//  void MatPc(Vector *out, Vector *in, Float mass, Float epsilon, DagType dag);

  void MatPc(Vector *out, Vector *in, Float mass, DagType dag);

  ForceArg EvolveMomFforce(Matrix *mom, Vector *frm,
			   Float mass, Float step_size);
  // It evolves the canonical momentum mom by step_size
  // using the fermion force.
  
  ForceArg EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
			   Float mass, Float step_size) {
    return EvolveMomFforceBase(mom, phi, eta, mass, -step_size);
  }

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
  
  void Freflex (Vector *out, Vector *in);
  //!< Does something really cool.
  // Reflexion in s operator, needed for the hermitian version 
  // of the dirac operator in the Ritz solver.

  int SpinComponents() const {
    return 4;
  }

  int ExactFlavors() const {
    if((bd.solver == HtCayleyEOFA) || (bd.solver == HmCayleyEOFA)){ return 1; }
    else{ return 2; }
  }
    
  //!< Method to ensure bosonic force works (does nothing for Wilson
  //!< theories.
  void BforceVector(Vector *in, CgArg *cg_arg);

  // !< Special for Mobius fermions, applies the D_- 5D matrix to an
  // !< unpreconditioned fermion vector.
  //
  // !< The following gives an example of D_- with Ls = 4:
  //       [ -D_-^1 0      0      0      ]
  //       [ 0      -D_-^2 0      0      ]
  // D_- = [ 0      0      -D_-^3 0      ]
  //       [ 0      0      0      -D_-^4 ]
  //
  // !< where D_-^s = c[s] D_W - 1, D_W is the 4D Wilson Dirac operator.
  void Dminus(Vector *out, Vector *in);

  //!< Toggle boundary condition
  //
  //!< Note: Agent classes which needs to import gauge field to
  //!external libraries need to overwrite this function.
#ifndef NO_BFM_BC
  void BondCond(){
    Lattice::BondCond();
    if (bfm_initted) ImportGauge();
  }
#endif



  void ImportGauge();

    void SetBfmArg(Float key_mass);

    void SetMassArg(Float mass){
	VRB.Result(cname,"SetMassArg(F)","called\n");
	SetBfmArg(mass);
    }

    void SetMass(Float mass) {
	const char *fname="SetMass(Float)";
    if(!bfm_initted) ERR.General(cname,fname,"Fbfm not initted\n");
        if(bd.mass != mass) {
            bd.mass = mass;
            bd.GeneralisedFiveDimEnd();
            bd.GeneralisedFiveDimInit();
        }
//        if(use_mixed_solver && bf.mass != mass) {
        if( bf.mass != mass) {
            bf.mass = mass;
            bf.GeneralisedFiveDimEnd();
            bf.GeneralisedFiveDimInit();
        }
//        if( kernel.mass != mass) {
//            kernel.mass = mass;
//            kernel.GeneralisedFiveDimEnd();
//            kernel.GeneralisedFiveDimInit();
//        }
}
  //CK: Added WilsonTm twist parameter epsilon. Use -12345 as a default for non WilsonTm fermions. Add a catch if using WilsonTm and this value of epsilon is passed in
  void SetMass(Float mass, Float epsilon) {
	const char *fname="SetMass(Float,Float)";
    if(!bfm_initted) ERR.General(cname,fname,"Fbfm not initted\n");
    if(epsilon == -12345 && CurrentSolver() == WilsonTM) ERR.General("Fbfm","SetMass(Float,Float)","Must specify epsilon for twisted mass Wilson fermions"); 
    if(epsilon != -12345 && CurrentSolver() != WilsonTM) ERR.General("Fbfm","SetMass(Float,Float)","Must specify epsilon for twisted mass Wilson fermions"); 

    if(bd.mass != mass || bd.twistedmass != epsilon) {
      bd.mass = mass;
      bd.twistedmass = epsilon;
      bd.GeneralisedFiveDimEnd();
      bd.GeneralisedFiveDimInit();
    }
    if(use_mixed_solver && (bf.mass != mass || bf.twistedmass != epsilon) ) {
      bf.mass = mass;
      bf.twistedmass = epsilon;
      bf.GeneralisedFiveDimEnd();
      bf.GeneralisedFiveDimInit();
    }
  }
        void Fdslash(Vector *f_out, Vector *f_in, CgArg *cg_arg,
                    CnvFrmType cnv_frm, int dir_flag);

    void SaveFmatEvlInvProblem(const char* path, Vector *f_in, CgArg *cg_arg);
    void LoadFmatEvlInvProblem(const char* path, int save_number, Vector *f_in, CgArg *cg_arg);

    //-------------------------------------------------------------------------------------------
    // DJM: additions for EOFA

    bool EOFA_solver(){ return ((bd.solver == HtCayleyEOFA) || (bd.solver == HmCayleyEOFA)) ? true : false; }

    void SetEOFAParams(Float m1, Float m2, Float m3, Float a, int pm) override
    {
      bool reinit_flag(false);
      Float alpha = Fbfm::bfm_args[current_arg_idx].mobius_scale;

      if(bd.mass         != m1   ){ bd.mass         = m1;    reinit_flag = true; }
      if(bd.mq1          != m1   ){ bd.mq1          = m1;    reinit_flag = true; }
      if(bd.mq2          != m2   ){ bd.mq2          = m2;    reinit_flag = true; }
      if(bd.mq3          != m3   ){ bd.mq3          = m3;    reinit_flag = true; }
      if(bd.a            != a    ){ bd.a            = a;     reinit_flag = true; }
      if(bd.pm           != pm   ){ bd.pm           = pm;    reinit_flag = true; }
      if(bd.mobius_scale != alpha){ bd.mobius_scale = alpha; reinit_flag = true; }
      
      if(bf.mass         != m1   ){ bf.mass         = m1;    reinit_flag = true; }
      if(bf.mq1          != m1   ){ bf.mq1          = m1;    reinit_flag = true; }
      if(bf.mq2          != m2   ){ bf.mq2          = m2;    reinit_flag = true; }
      if(bf.mq3          != m3   ){ bf.mq3          = m3;    reinit_flag = true; }
      if(bf.a            != a    ){ bf.a            = a;     reinit_flag = true; }
      if(bf.pm           != pm   ){ bf.pm           = pm;    reinit_flag = true; }
      if(bf.mobius_scale != alpha){ bf.mobius_scale = alpha; reinit_flag = true; }

      if(reinit_flag) {
        bd.GeneralisedFiveDimEnd();
        bd.GeneralisedFiveDimInit();
        if(use_mixed_solver) {
          bf.GeneralisedFiveDimEnd();
          bf.GeneralisedFiveDimInit();
        }
      }
    }

    void ChiralProj(Vector* f_out, const Vector* f_in, int pm) override;
    Float kt(Float m1, Float m2) override;
    void g5R5(Vector* f_out, Vector* f_in) override;
    void Omega(Vector* f_out, Vector* f_in, int pm, int dag) override;
    void Delta(Vector* f_out, Vector* f_in, Float m1, Float m2, int pm) override;
    void Dtilde(Vector* f_out, Vector* f_in, Float m) override;
    void DtildeInv(Vector* f_out, Vector* f_in, Float m) override;
    void HeofapaD_proj(Vector* f_out, Vector* f_in, Float m1,
        Float m2, Float m3, Float a, int pm) override;
    int FmatEvlMeofa(Vector* f_out, Vector* f_in, CgArg* cg_arg,
        Float m1, Float m2, Float m3, Float a, int pm, 
        Vector* prec_soln = nullptr) override;
};

class GnoneFbfm
    : public virtual Lattice,
      public virtual Gnone,
      public virtual Fbfm {
private:
    const char *cname;
public:
    GnoneFbfm(void);
    virtual ~GnoneFbfm();
};

class GimprRectFbfm
    : public virtual Lattice,
      public virtual GimprRect,
      public virtual Fbfm {
private:
    const char *cname;
public:
    GimprRectFbfm(void);
    virtual ~GimprRectFbfm();
};

#endif

CPS_END_NAMESPACE

#endif
