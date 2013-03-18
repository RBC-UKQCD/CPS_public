#ifndef INCLUDED_FBFM_H__
#define INCLUDED_FBFM_H__

#include<config.h>

#ifdef USE_BFM
#include <util/lattice/bfm_evo.h>
#include <util/lattice/bfm_mixed_solver.h>
#endif

#include <util/lattice.h>

CPS_START_NAMESPACE

#ifdef USE_BFM
class Fbfm : public virtual Lattice {
public:
    // have to do this since lattice factory does not accept any input
    // parameters.
    static bfmarg bfm_arg;

    // set true to use single precision BFM object.
    static bool use_mixed_solver;

    bfm_evo<double> bd;
    bfm_evo<float> bf;
private:
    const char *cname;

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

    ForceArg EvolveMomFforceBaseThreaded(Matrix *mom,
                                         Vector *phi1,
                                         Vector *phi2,
                                         Float mass,
                                         Float coef);
    // It evolves the canonical Momemtum mom:
    // mom += coef * (phi1^\dag e_i(M) \phi2 + \phi2^\dag e_i(M^\dag) \phi1)
    // note: this function does not exist in the base Lattice class.

    ForceArg EvolveMomFforceBase(Matrix *mom,
                                 Vector *phi1,
                                 Vector *phi2,
                                 Float mass,
                                 Float coef);
    // It evolves the canonical Momemtum mom:
    // mom += coef * (phi1^\dag e_i(M) \phi2 + \phi2^\dag e_i(M^\dag) \phi1)
    // note: this function does not exist in the base Lattice class.

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
  
    int FsiteOffset(const int *x) const;
    // Sets the offsets for the fermion fields on a 
    // checkerboard. The fermion field storage order
    // is the canonical one. X[I] is the
    // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
  
    int FsiteSize() const {
        return 24 * Fbfm::bfm_arg.Ls;
    }
    // Returns the number of fermion field 
    // components (including real/imaginary) on a
    // site of the 4-D lattice.
  
    int FchkbEvl() const {
        return 1;
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
  
    void MatPc(Vector *out, Vector *in, Float mass, DagType dag);

    Float SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
                 Float mass, DagType dag);
    // It sets the pseudofermion field phi from frm1, frm2.
  
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
        return 2;
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
    virtual void BondCond();

    void ImportGauge();

    void SetMass(Float mass) {
        if(bd.mass != mass) {
            bd.mass = mass;
            bd.GeneralisedFiveDimEnd();
            bd.GeneralisedFiveDimInit();
        }
        if(use_mixed_solver && bf.mass != mass) {
            bf.mass = mass;
            bf.GeneralisedFiveDimEnd();
            bf.GeneralisedFiveDimInit();
        }
    }
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
