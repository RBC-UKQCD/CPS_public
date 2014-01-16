// -*- c-basic-offset: 4 -*-
#include<config.h>
#include<math.h>

#ifdef USE_BFM

#include <util/lattice/bfm_evo.h>
#include <util/lattice/bfm_eigcg.h>
//#include <util/lattice/bfm_hdcg.h>
#include <util/lattice/fbfm.h>
#include <util/wilson.h>
#include <util/verbose.h>
#include <util/gjp.h>
#include <util/error.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <util/enum_func.h>
#include <util/sproj_tr.h>
#include <util/time_cps.h>
#include <util/lattice/fforce_wilson_type.h>

#include<omp.h>
#include<bfm_hdcg_wrapper.h>
#include<BfmMultiGrid.h>

#if 0
class HDCGInstance{
	public:
	static BfmMultiGridParams  Params;
	static HDCG_wrapper  * _instance ;
	static HDCG_wrapper *getInstance(){return _instance;} 
	static HDCG_wrapper *setInstance(HDCG_wrapper *_new){_instance = _new;} 
	static void free(){
		if (_instance){ 
			_instance->HDCG_end();
			delete _instance;
		}
		_instance=NULL;
	}
};
#endif

#include <util/lattice/hdcg_controller.h>



HDCGInstance hdcg_instance; // to invoke constructor with defaults
BfmMultiGridParams HDCGInstance::Params;
HDCG_wrapper  *HDCGInstance:: _instance=NULL;



CPS_START_NAMESPACE

bfmarg Fbfm::bfm_arg;
bool Fbfm::use_mixed_solver = false;

// NOTE:
//
// 1. Initialize QDP++ and the static copy Fbfm::bfm_arg before
// using this class.
//
// 2. This class acts like a DiracOp class, while it is in scope the
// gauge field has the boundary condition on.
Fbfm::Fbfm(void):cname("Fbfm")
{
    const char *fname = "Fbfm()";
    VRB.Func(cname,fname);

    if(GJP.Snodes() != 1) {
        ERR.NotImplemented(cname, fname);
    }
    if(sizeof(Float) == sizeof(float)) {
        ERR.NotImplemented(cname, fname);
    }

    bd.init(bfm_arg);

    if(use_mixed_solver) {
//    if(1) {
        bd.comm_end();
        bf.init(bfm_arg);
        bf.comm_end();
        bd.comm_init();
    }

    // call our own version to import gauge field.
    Fbfm::BondCond();

    evec = NULL;
    evald = NULL;
    evalf = NULL;
    ecnt = 0;
}

Fbfm::~Fbfm(void)
{
    const char *fname = "~Fbfm()";
    VRB.Func(cname,fname);
    // we call base version just to revert the change, no need to
    // import to BFM in a destructor.
    Lattice::BondCond();
    bd.end();
    if(use_mixed_solver) {
//    if(1) {
        bf.end();
    }
}

// This function differs from the original CalcHmdForceVecsBilinear()
// in that it stores v1 and v2 in (color, spin, s, x, y, z, t) order
// to facilitate force evaluation.
void Fbfm::CalcHmdForceVecsBilinear(Float *v1,
                                    Float *v2,
                                    Vector *phi1,
                                    Vector *phi2,
                                    Float mass)
{
    Fermion_t pi[2] = {bd.allocFermion(), bd.allocFermion()};
    Fermion_t po[4] = {bd.allocFermion(), bd.allocFermion(),
                       bd.allocFermion(), bd.allocFermion()};

    SetMass(mass);
    bd.cps_impexcbFermion((Float *)phi1, pi[0], 1, 1);
    bd.cps_impexcbFermion((Float *)phi2, pi[1], 1, 1);

#pragma omp parallel
    {
        bd.calcMDForceVecs(po + 0, po + 2, pi[0], pi[1]);
    }

    bd.cps_impexFermion_s(v1, po + 0, 0);
    bd.cps_impexFermion_s(v2, po + 2, 0);

    bd.freeFermion(pi[0]);
    bd.freeFermion(pi[1]);
    bd.freeFermion(po[0]);
    bd.freeFermion(po[1]);
    bd.freeFermion(po[2]);
    bd.freeFermion(po[3]);
}

ForceArg Fbfm::EvolveMomFforceBaseThreaded(Matrix *mom,
                                           Vector *phi1, Vector *phi2,
                                           Float mass, Float coef)
{
    const char *fname = "EvolveMomFforceBaseThreaded()";

    Float dtime = -dclock();

    Fermion_t in[2] = {bd.allocFermion(), bd.allocFermion()};

    SetMass(mass);

    bd.cps_impexcbFermion((Float *)phi1, in[0], 1, 1);
    bd.cps_impexcbFermion((Float *)phi2, in[1], 1, 1);

    Float *gauge = (Float *)(this->GaugeField());
#pragma omp parallel
    {
        bd.compute_force((Float *)mom, gauge, in[0], in[1], coef);
    }

    bd.freeFermion(in[0]);
    bd.freeFermion(in[1]);

    dtime += dclock();

    VRB.Result(cname, fname, "takes %17.10e seconds\n", dtime);
    return ForceArg();
}

// It evolves the canonical Momemtum mom:
// mom += coef * (phi1^\dag e_i(M) \phi2 + \phi2^\dag e_i(M^\dag) \phi1)
//
// NOTE:
//
// 1. This function does not exist in the base Lattice class.
//
// 2. The 2 auxiliary vectors v1 and v2 calculated by
// CalcHmdForceVecsBilinear must be in (reim, color, spin, s, x, y, z,
// t) order.
//
// 3. For BFM M is M = M_oo - M_oe M^{-1}_ee M_eo
ForceArg Fbfm::EvolveMomFforceBase(Matrix *mom,
                                   Vector *phi1,
                                   Vector *phi2,
                                   Float mass,
                                   Float coef)
{
    const char *fname = "EvolveMomFforceBase()";

#if 0
    return EvolveMomFforceBaseThreaded(mom, phi1, phi2, mass, coef);
#endif

    long f_size = (long)SPINOR_SIZE * GJP.VolNodeSites() * Fbfm::bfm_arg.Ls;
    Float *v1 = (Float *)smalloc(cname, fname, "v1", sizeof(Float) * f_size);
    Float *v2 = (Float *)smalloc(cname, fname, "v2", sizeof(Float) * f_size);

    CalcHmdForceVecsBilinear(v1, v2, phi1, phi2, mass);

    FforceWilsonType cal_force(mom, this->GaugeField(),
                               v1, v2, Fbfm::bfm_arg.Ls, coef);
    ForceArg ret = cal_force.run();

    sfree(cname, fname, "v1", v1);
    sfree(cname, fname, "v2", v2);
    return ret;
}

//------------------------------------------------------------------
//! Multiplication of a lattice spin-colour vector by gamma_5.
//------------------------------------------------------------------
void Fbfm::Gamma5(Vector *v_out, Vector *v_in, int num_sites)
{
    Float *p_out = (Float *)v_out;
    Float *p_in  = (Float *)v_in;

    int half_site_size = 12 ;
    for (int site = 0; site < num_sites; ++site) {

        for(int comp = 0; comp < half_site_size; ++comp) {
            *p_out++ = *p_in++ ;
        }
        for(int comp = 0; comp < half_site_size; ++comp) {
            *p_out++ = -*p_in++ ;
        }
    }
}

// Sets the offsets for the fermion fields on a 
// checkerboard. The fermion field storage order
// is not the canonical one but it is particular
// to the Dwf fermion type. x[i] is the 
// ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
int Fbfm::FsiteOffsetChkb(const int *x) const
{
    const char *fname = "FsiteOffsetChkb()";
    ERR.NotImplemented(cname, fname);
}

// Sets the offsets for the fermion fields on a 
// checkerboard. The fermion field storage order
// is the canonical one. X[I] is the
// ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
int Fbfm::FsiteOffset(const int *x) const
{
    const char *fname = "FsiteOffset()";
    ERR.NotImplemented(cname, fname);
}

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
int Fbfm::FmatEvlInv(Vector *f_out, Vector *f_in, 
                     CgArg *cg_arg, 
                     Float *true_res,
                     CnvFrmType cnv_frm)
{
    const char *fname = "FmatEvlInv()";

    if(cg_arg == NULL)
        ERR.Pointer(cname, fname, "cg_arg");

    Fermion_t in  = bd.allocFermion();
    Fermion_t out = bd.allocFermion();

    SetMass(cg_arg->mass);
    bd.residual = cg_arg->stop_rsd;
    bd.max_iter = bf.max_iter = cg_arg->max_num_iter;
    // FIXME: pass single precision rsd in a reasonable way.
    bf.residual = 1e-5;

    bd.cps_impexcbFermion((Float *)f_in , in,  1, 1);
    bd.cps_impexcbFermion((Float *)f_out, out, 1, 1);

    int iter = -1;
#pragma omp parallel
    {
        iter =
            use_mixed_solver 
            ? mixed_cg::threaded_cg_mixed_MdagM(out, in, bd, bf, 5)
            : bd.CGNE_prec_MdagM(out, in);
    }

    bd.cps_impexcbFermion((Float *)f_out, out, 0, 1);

    bd.freeFermion(in);
    bd.freeFermion(out);

    return iter;
}

int Fbfm::FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
                      int Nshift, int isz, CgArg **cg_arg, 
                      CnvFrmType cnv_frm, MultiShiftSolveType type, Float *alpha,
                      Vector **f_out_d)
{
    const char *fname = "FmatEvlMInv(V*,V*,F*, ...)";
  
    if(isz != 0) {
        ERR.General(cname, fname, "Non-zero isz is not implemented.\n");
    }

    Fermion_t *sol_multi = new Fermion_t[Nshift];
    double *ones = new double[Nshift];
    double *mresidual = new double[Nshift];
    for(int i = 0; i < Nshift; ++i) {
        sol_multi[i] = bd.allocFermion();
        ones[i] = 1.0;
        mresidual[i] = cg_arg[i]->stop_rsd;
    }

    // source
    Fermion_t src = bd.allocFermion();
    bd.cps_impexcbFermion((Float *)f_in, src, 1, 1);

    SetMass(cg_arg[0]->mass);
    bd.residual = cg_arg[0]->stop_rsd;
    bd.max_iter = cg_arg[0]->max_num_iter;

    int iter;
#pragma omp parallel
    {
        iter = bd.CGNE_prec_MdagM_multi_shift(sol_multi, src, shift, ones, Nshift, mresidual, 0);
    }

    if(type == SINGLE) {
        // FIXME
        int f_size_cb = GJP.VolNodeSites() * SPINOR_SIZE * Fbfm::bfm_arg.Ls / 2;
        Vector *t = (Vector *)smalloc(cname, fname, "t", sizeof(Float) * f_size_cb);

        for(int i = 0; i < Nshift; ++i) {
            bd.cps_impexcbFermion((Float *)t, sol_multi[i], 0, 1);
            f_out[0]->FTimesV1PlusV2(alpha[i], t, f_out[0], f_size_cb);
        }
        sfree(cname, fname, "t", t);
    } else {
        for(int i = 0; i < Nshift; ++i) {
            bd.cps_impexcbFermion((Float *)f_out[i], sol_multi[i], 0, 1);
        }
    }

    bd.freeFermion(src);
    for(int i = 0; i < Nshift; ++i) {
        bd.freeFermion(sol_multi[i]);
    }

    delete[] sol_multi;
    delete[] ones;
    delete[] mresidual;

    return iter;
}

void Fbfm::FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
                      Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm)
{
    const char *fname = "FminResExt(V*, V*, V**, ...)";

    int f_size_cb = GJP.VolNodeSites() * SPINOR_SIZE * Fbfm::bfm_arg.Ls / 2;

    // does nothing other than setting sol to zero
    sol->VecZero(f_size_cb);
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
int Fbfm::FmatInv(Vector *f_out, Vector *f_in,
                  CgArg *cg_arg,
                  Float *true_res,
                  CnvFrmType cnv_frm,
                  PreserveType prs_f_in)
{
    const char *fname = "FmatInv()";
    VRB.Func(cname,fname);

    if(cg_arg == NULL)
        ERR.Pointer(cname, fname, "cg_arg");
    int threads =   omp_get_max_threads();
if (cg_arg->Inverter == HDCG){
    if(!use_mixed_solver)
	ERR.General(cname,fname,"Fbfm::use_mixed_solver should be set true to use HDCG\n");
//    HDCGController<Float> *control = HDCGController<Float>::getInstance();
	HDCG_wrapper *control = HDCGInstance::getInstance();
    if (!control){
	bfmActionParams BAP_;
	BAP_.M5 = bfm_arg.M5;
	BAP_.mass = cg_arg->mass;
	BAP_.twistedmass=0;
	BAP_.Csw=0;
	BAP_.solver=bfm_arg.solver;
	BAP_.mobius_scale=bfm_arg.mobius_scale;
	BAP_.zolo_hi=bfm_arg.zolo_hi;
	BAP_.zolo_lo=bfm_arg.zolo_lo;
	BAP_.Ls=bfm_arg.Ls;
	BAP_.precon_5d=bfm_arg.precon_5d;

	BAP_.solveMobiusDminus=1;

//	int _Ns = 20;
//	int block[5]={4,4,4,4,Ls};
//	int quad[4]={4,4,4,4};
//	control = HDCGController<Float>::setInstance(Ls,_Ns,block,quad);
	control = new HDCG_wrapper;
	HDCGInstance::setInstance(control);
	HDCGInstance::Params.SubspaceRationalMass=cg_arg->mass;
//	HDCGInstance::Params.PreconditionerKrylovResidual=1e-4;
//	HDCGInstance::Params.PreconditionerKrylovIterMax=8;
//	HDCGInstance::Params.PreconditionerKrylovShift=1.0;
	control->HDCG_init( HDCGInstance::Params, BAP_); 
    	Float *gauge = (Float *)(this->GaugeField());
	control->HDCG_gauge_import_cps<Float>(gauge);
    	SetMass(cg_arg->mass);
	control->HDCG_set_mass(cg_arg->mass);
    	VRB.Result(cname,fname,"HDCG first called with nthreads=%d. Initialzing Ldop\n",threads);
	control->HDCG_subspace_init();
	control->HDCG_subspace_compute(0);
	control->HDCG_subspace_refine();
	control->HDCG_subspace_compute(0);
//    }
//    BfmMultiGrid<Float> *hdcg = control->ldop_d;
//    if (!hdcg){
//	control->setHDCG(bd,bf); 
//	hdcg = control->getHDCG();
    }
////    hdcg->InnerKrylovIterMax  = 8;
//    hdcg->SubspaceSurfaceDepth= 256;
//    hdcg->InnerKrylovShift    = 1.0;
//    hdcg->PcgShift       = 1.0;
//    hdcg->PcgType        = PcgAdef2f;
//    filecg=file+".PcgADef2f_1e-2";
//    bd.InverterLoggingBegin("BfmMultiGrid.log");
}

    Fermion_t in[2]  = {bd.allocFermion(), bd.allocFermion()};
    Fermion_t out[2] = {bd.allocFermion(), bd.allocFermion()};

    SetMass(cg_arg->mass);
    bd.residual = cg_arg->stop_rsd;
    bd.max_iter = bf.max_iter = cg_arg->max_num_iter;
    // FIXME: pass single precision rsd in a reasonable way.
    bf.residual = 1e-5;

    // deal with Mobius Dminus
    if(bd.solver == HmCayleyTanh) {
        bd.cps_impexFermion((Float *)f_in , out,  1);
#pragma omp parallel
        {
            bd.G5D_Dminus(out, in, 0);
        }
    } else {
        bd.cps_impexFermion((Float *)f_in , in,  1);
    }

    bd.cps_impexFermion((Float *)f_out, out, 1);

    int iter = -1;
if((cg_arg->Inverter == HDCG)) {
	HDCG_wrapper *control = HDCGInstance::getInstance();
	control->HDCG_set_mass(cg_arg->mass);
        control->HDCG_invert(out, in, cg_arg->stop_rsd, cg_arg->max_num_iter);
}
else {
#pragma omp parallel
    {
        if(use_mixed_solver && (cg_arg->Inverter != HDCG)) {
            iter = mixed_cg::threaded_cg_mixed_M(out, in, bd, bf, 5, cg_arg->Inverter, evec, evalf, ecnt);
        } else {
            switch(cg_arg->Inverter) {
            case CG:
                if(evec && evald && ecnt) {
                    iter = bd.CGNE_M(out, in, *evec, *evald);
                } else {
                    iter = bd.CGNE_M(out, in);
                }
                break;
            case EIGCG:
                iter = bd.EIG_CGNE_M(out, in);
                break;
            case HDCG:
                if(bd.isBoss()) {
                    printf("%s::%s: HDCG implemented outside threaded region. Shouldn't have reached this line!\n", cname, fname);
                }
                break;
            default:
                if(bd.isBoss()) {
                    printf("%s::%s: Not implemented\n", cname, fname);
                }
                exit(-1);
                break;
            }
        }
    }
}
if (cg_arg->Inverter == HDCG){
//    HDCGController<Float> *control = HDCGController::getInstance();
//    BfmMultiGrid<Float> *hdcg = control->getHDCG();
    bd.InverterLoggingEnd();
}

    bd.cps_impexFermion((Float *)f_out, out, 0);

    bd.freeFermion(in[0]);
    bd.freeFermion(in[1]);
    bd.freeFermion(out[0]);
    bd.freeFermion(out[1]);

    return iter;
}

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
void Fbfm::Ffour2five(Vector *five, Vector *four, int s_u, int s_l, int Ncb)
{
    const char *fname = "Ffour2five(V*, V*, ...)";

    // Note: we don't allow splitting in s direction
    if(GJP.Snodes() != 1) {
        ERR.NotImplemented(cname, fname);
    }
    // what does Ncb do?
    if(Ncb != 2) {
        ERR.NotImplemented(cname, fname);
    }

    Float *f5d = (Float *)five;
    Float *f4d = (Float *)four;

    const int size_4d = GJP.VolNodeSites() * SPINOR_SIZE;
    const int size_5d = size_4d * Fbfm::bfm_arg.Ls;

    // zero 5D vector
#pragma omp parallel for
    for(int i=0; i< size_5d; ++i) {
        f5d[i]  = 0.0;
    }

    Float *f4du = f4d;
    Float *f4dl = f4d + 12;
    Float *f5du = f5d + s_u * size_4d;
    Float *f5dl = f5d + s_l * size_4d + 12;

#pragma omp parallel for
    for(int x = 0; x < size_4d; x += SPINOR_SIZE) {
        memcpy(f5du + x, f4du + x, sizeof(Float) * 12);
        memcpy(f5dl + x, f4dl + x, sizeof(Float) * 12);
    }
}

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
void Fbfm::Ffive2four(Vector *four, Vector *five, int s_u, int s_l, int Ncb)
{
    const char *fname = "Ffive2four(V*,V*,i,i)";

    // Note: we don't allow splitting in s direction
    if(GJP.Snodes() != 1) {
        ERR.NotImplemented(cname, fname);
    }
    // what does Ncb do?
    if(Ncb != 2) {
        ERR.NotImplemented(cname, fname);
    }

    Float *f5d = (Float *)five;
    Float *f4d = (Float *)four;

    const int size_4d = GJP.VolNodeSites() * SPINOR_SIZE;

    // zero 4D vector
#pragma omp parallel for
    for(int i=0; i< size_4d; ++i) {
        f4d[i]  = 0.0;
    }

    Float *f4du = f4d;
    Float *f4dl = f4d + 12;
    Float *f5du = f5d + s_u * size_4d;
    Float *f5dl = f5d + s_l * size_4d + 12;

#pragma omp parallel for
    for(int x = 0; x < size_4d; x += SPINOR_SIZE) {
        memcpy(f4du + x, f5du + x, sizeof(Float) * 12);
        memcpy(f4dl + x, f5dl + x, sizeof(Float) * 12);
    }
}

// It finds the eigenvectors and eigenvalues of A where
// A is the fermion matrix (Dirac operator). The solution
// uses Ritz minimization. eig_arg is the 
// structure that contains all the control parameters, f_eigenv
// are the fermion field source vectors which should be
// defined initially, lambda are the eigenvalues returned 
// on solution. f_eigenv is defined on the whole lattice.
// The function returns the total number of Ritz iterations.
int Fbfm::FeigSolv(Vector **f_eigenv, Float *lambda,
                   Float *chirality, int *valid_eig,
                   Float **hsum,
                   EigArg *eig_arg, 
                   CnvFrmType cnv_frm)
{
    const char *fname = "FeigSolv(EigArg*,V*,F*,CnvFrmType)";

    // only 1 eigenvalue can be computed now.
    if(eig_arg->N_eig != 1) {
        ERR.NotImplemented(cname, fname);
    }
    if(eig_arg->RitzMatOper != MATPCDAG_MATPC &&
       eig_arg->RitzMatOper != NEG_MATPCDAG_MATPC) {
        ERR.NotImplemented(cname, fname);
    }
    
    SetMass(eig_arg->mass);
    bd.residual = eig_arg->Rsdlam;
    bd.max_iter = eig_arg->MaxCG;

    VRB.Result(cname, fname, "residual = %17.10e max_iter = %d mass = %17.10e\n",
               bd.residual, bd.max_iter, bd.mass);

    Fermion_t in = bd.allocFermion();
    bd.cps_impexcbFermion((Float *)f_eigenv[0], in, 1, 1);

#pragma omp parallel
    {
        lambda[0] = bd.ritz(in, eig_arg->RitzMatOper == MATPCDAG_MATPC);
    }

    bd.cps_impexcbFermion((Float *)f_eigenv[0], in, 0, 1);

    // correct the eigenvalue for a dumb convention problem.
    if(eig_arg->RitzMatOper == NEG_MATPCDAG_MATPC) lambda[0] = -lambda[0];

    valid_eig[0] = 1;
    bd.freeFermion(in);

    return 0;
}

// It sets the pseudofermion field phi from frm1, frm2.
Float Fbfm::SetPhi(Vector *phi, Vector *frm1, Vector *frm2,	       
                   Float mass, DagType dag)
{
    const char *fname = "SetPhi(V*,V*,V*,F)";

    if (phi == 0)
        ERR.Pointer(cname,fname,"phi") ;

    if (frm1 == 0)
        ERR.Pointer(cname,fname,"frm1") ;

    MatPc(phi, frm1, mass, dag);
    return FhamiltonNode(frm1, frm1);
}

void Fbfm::MatPc(Vector *out, Vector *in, Float mass, DagType dag)
{
    const char *fname = "MatPc()";

    Fermion_t i = bd.allocFermion();
    Fermion_t o = bd.allocFermion();
    Fermion_t t = bd.allocFermion();

    SetMass(mass);

    bd.cps_impexcbFermion((Float *)in , i, 1, 1);
#pragma omp parallel
    {
        bd.Mprec(i, o, t, dag == DAG_YES, 0);
    }
    bd.cps_impexcbFermion((Float *)out, o, 0, 1);

    bd.freeFermion(i);
    bd.freeFermion(o);
    bd.freeFermion(t);
}

// It evolves the canonical momentum mom by step_size
// using the fermion force.
ForceArg Fbfm::EvolveMomFforce(Matrix *mom, Vector *frm, 
                               Float mass, Float step_size)
{
    const char *fname = "EvolveMomFforce()";
  
    const int f_size_4d = SPINOR_SIZE * GJP.VolNodeSites();
    const int f_size_cb = f_size_4d * Fbfm::bfm_arg.Ls / 2;
  
    Vector *tmp = (Vector *)smalloc(cname, fname, "tmp", sizeof(Float)*f_size_cb);
    MatPc(tmp, frm, mass, DAG_NO);

    ForceArg f_arg = EvolveMomFforceBase(mom, tmp, frm, mass, step_size);
    sfree(cname, fname, "tmp", tmp);

    return f_arg;
}

ForceArg Fbfm::RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
                                    int isz, Float *alpha, Float mass, Float dt,
                                    Vector **sol_d, ForceMeasure force_measure)
{
    const char *fname = "RHMC_EvolveMomFforce()";
    char *force_label=NULL;

    int g_size = GJP.VolNodeSites() * GsiteSize();

    Matrix *mom_tmp;

    if (force_measure == FORCE_MEASURE_YES) {
        mom_tmp = (Matrix*)smalloc(g_size * sizeof(Float),cname, fname, "mom_tmp");
        ((Vector*)mom_tmp) -> VecZero(g_size);
        force_label = new char[100];
    } else {
        mom_tmp = mom;
    }

    for (int i=0; i<degree; i++) {
        ForceArg Fdt = EvolveMomFforce(mom_tmp, sol[i], mass, alpha[i]*dt);

        if (force_measure == FORCE_MEASURE_YES) {
            sprintf(force_label, "Rational, mass = %e, pole = %d:", mass, i+isz);
            Fdt.print(dt, force_label);
        }
    }

    ForceArg ret;

    // If measuring the force, need to measure and then sum to mom
    if (force_measure == FORCE_MEASURE_YES) {
        ret.measure(mom_tmp);
        ret.glb_reduce();

        fTimesV1PlusV2((IFloat*)mom, 1.0, (IFloat*)mom_tmp, (IFloat*)mom, g_size);

        delete[] force_label;
        sfree(mom_tmp, cname, fname, "mom_tmp");
    }

    return ret;
}

// The fermion Hamiltonian of the node sublattice.
// chi must be the solution of Cg with source phi.
// copied from FdwfBase
Float Fbfm::FhamiltonNode(Vector *phi, Vector *chi)
{
    const char *fname = "FhamiltonNode(V*, V*)";

    if (phi == 0) ERR.Pointer(cname, fname, "phi");
    if (chi == 0) ERR.Pointer(cname, fname, "chi");

    int f_size = GJP.VolNodeSites() * FsiteSize() / 2;

    // Sum accross s nodes is not necessary for MDWF since the library
    // does not allow lattice splitting in s direction.
    return phi->ReDotProductNode(chi, f_size);
}

// Convert fermion field f_field from -> to
// Moved to fbfm.h by CJ
#if 0
void Fbfm::Fconvert(Vector *f_field,
                    StrOrdType to,
                    StrOrdType from)
{
    const char *fname = "Fconvert()";
    VRB.Func(cname,fname);

    // nothing needs to be done
//    ERR.NotImplemented(cname, fname);
}
#endif


// The boson Hamiltonian of the node sublattice
Float Fbfm::BhamiltonNode(Vector *boson, Float mass)
{
    const char *fname = "BhamiltonNode()";
    ERR.NotImplemented(cname, fname);
}

// Reflexion in s operator, needed for the hermitian version 
// of the dirac operator in the Ritz solver.
void Fbfm::Freflex(Vector *out, Vector *in)
{
    const char *fname = "Freflex(V*,V*)";
    ERR.NotImplemented(cname, fname);
}

//!< Method to ensure bosonic force works (does nothing for Wilson
//!< theories.
void Fbfm::BforceVector(Vector *in, CgArg *cg_arg)
{
    return;
}

// !< Special for Mobius fermions, applies the D_- 5D matrix to an
// !< unpreconditioned fermion vector.
//
// !< The following gives an example of D_- with Ls = 4:
//       [ D_-^1 0      0      0     ]
//       [ 0     D_-^2  0      0     ]
// D_- = [ 0     0      D_-^3  0     ]
//       [ 0     0      0      D_-^4 ]
//
// !< where D_-^s = 1 - c[s] D_W, D_W is the 4D Wilson Dirac operator.
void Fbfm::Dminus(Vector *out, Vector *in)
{
    const char *fname = "Dminus(V*, V*)";

    // should be very easy to implement ...
    ERR.NotImplemented(cname, fname);
}

void Fbfm::BondCond()
{
    Lattice::BondCond();
    ImportGauge();
}

void Fbfm::ImportGauge()
{
    Float *gauge = (Float *)(this->GaugeField());
    bd.cps_importGauge(gauge);
    if(use_mixed_solver) {
        bd.comm_end();
        bf.comm_init();
        bf.cps_importGauge(gauge);
        bf.comm_end();
        bd.comm_init();
    }
}

CPS_END_NAMESPACE

#endif
