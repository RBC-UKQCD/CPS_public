// -*- c-basic-offset: 4 -*-
#include<config.h>
#include<math.h>

#include<util/multi_cg_controller.h>
CPS_START_NAMESPACE
MultiShiftCGcontroller MultiShiftController;
CPS_END_NAMESPACE
#ifdef USE_BFM

#include <util/lattice/bfm_evo.h>
#include <util/lattice/bfm_eigcg.h>
//#include <util/lattice/bfm_hdcg.h>
#include <util/lattice/fbfm.h>
#include <util/wilson.h>
#include <util/verbose.h>
#include <util/gjp.h>
#include <util/error.h>
#include <util/pt.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <util/enum_func.h>
#include <util/sproj_tr.h>
#include <util/time_cps.h>
#include <util/lattice/fforce_wilson_type.h>
#include <util/timer.h>
#include <util/lattice/hdcg_controller.h>

#include<omp.h>

#include <util/qioarg.h>
#include <util/WriteLatticePar.h>
#include <util/ReadLatticePar.h>

// These have to be defined somewhere. Why not here?
BfmHDCGParams HDCGInstance::Params;
HDCG_wrapper  *HDCGInstance::_instance = NULL;

CPS_START_NAMESPACE

std::map<Float, bfmarg> Fbfm::arg_map;
Float Fbfm::current_key_mass = -1789.8;

std::map<Float, MADWFParams> Fbfm::madwf_arg_map;

bool Fbfm::use_mixed_solver = false;


// NOTE:
//
// 1. Initialize QDP++ and the static Fbfm::arg_map before
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

    Lattice::BondCond();

    bfm_inited = false;
    Float current_key_mass = -1789.8;

    evec = NULL;
    evald = NULL;
    evalf = NULL;
    ecnt = 0;
}

Fbfm::~Fbfm(void)
{
    const char *fname = "~Fbfm()";
    VRB.Result(cname, fname,"start");
    // we call base version just to revert the change, no need to
    // import to BFM in a destructor.
    Lattice::BondCond();
    VRB.Result(cname, fname,"BondCond");

    if (bfm_inited) {
	bd.end();
    VRB.Result(cname, fname,"bd.end()");
#if 0
	kernel.end();
    VRB.Result(cname, fname,"kernel.end()");
#endif
	if (use_mixed_solver) {
	    bf.end();
    VRB.Result(cname, fname,"bf.end()");
	}
    }
}

// automatically fills in some bfmarg fields
void AutofillBfmarg(bfmarg &arg)
{
    // Make sure some fields are filled in properly
    multi1d<int> sub_latt_size = QDP::Layout::subgridLattSize();
    arg.node_latt[0] = sub_latt_size[0];
    arg.node_latt[1] = sub_latt_size[1];
    arg.node_latt[2] = sub_latt_size[2];
    arg.node_latt[3] = sub_latt_size[3];

    multi1d<int> procs = QDP::Layout::logicalSize();
    arg.local_comm[0] = procs[0] > 1 ? 0 : 1;
    arg.local_comm[1] = procs[1] > 1 ? 0 : 1;
    arg.local_comm[2] = procs[2] > 1 ? 0 : 1;
    arg.local_comm[3] = procs[3] > 1 ? 0 : 1;

    arg.ncoor[0] = 0;
    arg.ncoor[1] = 0;
    arg.ncoor[2] = 0;
    arg.ncoor[3] = 0;

    arg.max_iter = 100000;
    arg.verbose = BfmMessage | BfmError;
}

void Fbfm::SetBfmArg(Float key_mass)
{
    const char* fname = "SetBfmArg(F)";

    if (arg_map.count(key_mass) == 0) {
	ERR.General(cname, fname, "No entry for key mass %e in arg_map!\n", key_mass);
    }

    VRB.Result(cname, fname, "SetBfmArg: (Re)initing BFM objects from key mass %e)\n", key_mass);

    bfmarg new_arg = arg_map.at(key_mass);
    bfmarg kernel_arg = new_arg;
    kernel_arg.solver=DWFKernel;
    kernel_arg.Ls=1;

    if (!bfm_inited) {
	AutofillBfmarg(new_arg);
//	AutofillBfmarg(kernel_arg);
 
	bd.init(new_arg);
	if (use_mixed_solver) {
	    bd.comm_end();
	    bf.init(new_arg);
	    bf.comm_end();
	    bd.comm_init();
	}
#if 0
	bd.comm_end();
	kernel.init(kernel_arg);
	kernel.comm_end();
	bd.comm_init();
#endif

	ImportGauge();
	VRB.Result(cname, fname, "inited BFM objects with new BFM arg: solver = %d, mass = %e, Ls = %d, mobius_scale = %e\n", bd.solver, bd.mass, bd.Ls, bd.mobius_scale);
    } else {
	if (key_mass == current_key_mass) {
	    VRB.Result(cname, fname, "Already inited from desired key mass %e\n", key_mass);
	    return; // already inited with desired params
	}

        bool bad_change = false;
	if (bd.solver != new_arg.solver || bd.CGdiagonalMee != new_arg.CGdiagonalMee) {
          bad_change = true;
        } else if (bd.solver != WilsonTM) {
          if (bd.mobius_scale != new_arg.mobius_scale
	      || bd.Ls != new_arg.Ls
	      || bd.precon_5d != new_arg.precon_5d) {
            bad_change = true;
          }
        }
        if (bad_change) {
	    ERR.General(cname, fname, "Can't change solver, mobius_scale, Ls, precon_5d, or CGdiagonalMee "
                "during lifetime of Fbfm object: must destroy and recreate lattice object "
                "(solver=%d->%d, mobius_scale=%e->%e, Ls=%d->%d, precon_5d=%d->%d, CGdiagonalMee=%d->%d).\n", 
                bd.solver, new_arg.solver, bd.mobius_scale, new_arg.mobius_scale, bd.Ls, new_arg.Ls, 
                bd.precon_5d, new_arg.precon_5d, bd.CGdiagonalMee, new_arg.CGdiagonalMee);
	}

	if (new_arg.solver == WilsonTM) {
	    bd.mass = new_arg.mass;
	    bf.mass = new_arg.mass;
	    bd.twistedmass = new_arg.twistedmass;
	    bf.twistedmass = new_arg.twistedmass;
	} else {
	    SetMass(new_arg.mass);
	}
	VRB.Result(cname, fname, "Just set new mass %e for solver = %d, Ls = %d\n", bd.mass, bd.solver, bd.Ls);
    }

    bfm_inited = true;
    current_key_mass = key_mass;
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
    SetBfmArg(mass);

    VRB.Result(cname, "CalcHmdForceVecsBilinear()", "bd.CGdiagonalMee = %d\n", bd.CGdiagonalMee);

    Fermion_t pi[2] = { bd.allocFermion(), bd.allocFermion() };
    Fermion_t po[4] = {bd.allocFermion(), bd.allocFermion(),
                       bd.allocFermion(), bd.allocFermion()};
    Fermion_t tmp = bd.allocFermion();

    bd.cps_impexcbFermion((Float *)phi1, pi[0], 1, 1);
    bd.cps_impexcbFermion((Float *)phi2, pi[1], 1, 1);

#pragma omp parallel
    {
	// For CGdiagonalMee == 2 there is an extra factor of
	// Moo^{-1} in front of phi2_o
	if (bd.CGdiagonalMee == 2) {
	    bd.MooeeInv(pi[1], tmp, DaggerNo);
	    bd.copy(pi[1], tmp);
	}

	// For CGdiagonalMee == 1 there is an extra factor of
	// Moo^{\dag-1} in front of phi1_o
	if (bd.CGdiagonalMee == 1) {
	    bd.MooeeInv(pi[0], tmp, DaggerYes);
	    bd.copy(pi[0], tmp);
	}

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
    bd.freeFermion(tmp);
}

ForceArg Fbfm::EvolveMomFforceBaseThreaded(Matrix *mom,
                                           Vector *phi1, Vector *phi2,
                                           Float mass, Float coef)
{
    const char *fname = "EvolveMomFforceBaseThreaded()";

    Float dtime = -dclock();

    SetBfmArg(mass);

    Fermion_t in[2] = { bd.allocFermion(), bd.allocFermion() };

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
    static Timer time(cname, fname);
    time.start(true);

    SetBfmArg(mass);

#if 0
    return EvolveMomFforceBaseThreaded(mom, phi1, phi2, mass, coef);
#endif

    long f_size = (long)SPINOR_SIZE * GJP.VolNodeSites() * bd.Ls;
    Float *v1 = (Float *)smalloc(cname, fname, "v1", sizeof(Float) * f_size);
    Float *v2 = (Float *)smalloc(cname, fname, "v2", sizeof(Float) * f_size);

    CalcHmdForceVecsBilinear(v1, v2, phi1, phi2, mass);

    FforceWilsonType cal_force(mom, this->GaugeField(),
	v1, v2, bd.Ls, coef);
    ForceArg ret = cal_force.run();

    sfree(cname, fname, "v1", v1);
    sfree(cname, fname, "v2", v2);

    time.stop(true);
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
    int iter = -1;

    static Timer timer(cname, fname);
    static std::map<Float, Timer*> timers;
    if (timers.count(cg_arg->mass) == 0) {
	char timer_mass_name[512];
	sprintf(timer_mass_name, "FmatEvlInv(mass=%0.4f)", cg_arg->mass);
	timers[cg_arg->mass] = new Timer(cname, timer_mass_name);
    }
    timer.start(true);
    timers[cg_arg->mass]->start(true);

    SetBfmArg(cg_arg->mass);

    VRB.Result(cname, fname, "target residual = %e\n", cg_arg->stop_rsd);

    if (cg_arg == NULL)
	ERR.Pointer(cname, fname, "cg_arg");

    Fermion_t in = bd.allocFermion();
    Fermion_t out = bd.allocFermion();

    bd.residual = cg_arg->stop_rsd;
    bd.max_iter = bf.max_iter = cg_arg->max_num_iter;
    // FIXME: pass single precision rsd in a reasonable way.
    bf.residual = 1e-5;

    bd.cps_impexcbFermion((Float *)f_in, in, 1, 1);
    bd.cps_impexcbFermion((Float *)f_out, out, 1, 1);

#pragma omp parallel
    {
	iter = use_mixed_solver ?
	    mixed_cg::threaded_cg_mixed_MdagM(out, in, bd, bf, 5) :
	    bd.CGNE_prec_MdagM(out, in);
    }

    bd.cps_impexcbFermion((Float *)f_out, out, 0, 1);

    bd.freeFermion(in);
    bd.freeFermion(out);

    timers[cg_arg->mass]->stop(true);
    timer.stop(true);

    return iter;
}

int Fbfm::FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift,
    int Nshift, int isz, CgArg **cg_arg,
    CnvFrmType cnv_frm, MultiShiftSolveType type, Float *alpha,
    Vector **f_out_d)
{
    const char *fname = "FmatEvlMInv(V*,V*,F*, ...)";

    if (isz != 0) {
	ERR.General(cname, fname, "Non-zero isz is not implemented.\n");
    }

    SetBfmArg(cg_arg[0]->mass);

    Fermion_t *sol_multi = new Fermion_t[Nshift];
    double *ones = new double[Nshift];
    double *mresidual = new double[Nshift];
    for (int i = 0; i < Nshift; ++i) {
	sol_multi[i] = bd.allocFermion();
	ones[i] = 1.0;
	mresidual[i] = cg_arg[i]->stop_rsd;
    }

    // source
    Fermion_t src = bd.allocFermion();
    bd.cps_impexcbFermion((Float *)f_in, src, 1, 1);

    bd.residual = cg_arg[0]->stop_rsd;
    bd.max_iter = cg_arg[0]->max_num_iter;

    int iter;
    if (use_mixed_solver && bd.solver != WilsonTM) {
	MultiShiftController.MInv(sol_multi, src, shift, Nshift, mresidual, ones, 0, bd, bf);
    } else {
#pragma omp parallel
	{
	    iter = bd.CGNE_prec_MdagM_multi_shift(sol_multi, src, shift, ones, Nshift, mresidual, 0);
	}
    }

    if (type == SINGLE) {
	// FIXME
	int f_size_cb = GJP.VolNodeSites() * SPINOR_SIZE * bd.Ls / 2;
	Vector *t = (Vector *)smalloc(cname, fname, "t", sizeof(Float) * f_size_cb);

	for (int i = 0; i < Nshift; ++i) {
	    bd.cps_impexcbFermion((Float *)t, sol_multi[i], 0, 1);
	    f_out[0]->FTimesV1PlusV2(alpha[i], t, f_out[0], f_size_cb);
	}
	sfree(cname, fname, "t", t);
    } else {
	for (int i = 0; i < Nshift; ++i) {
	    bd.cps_impexcbFermion((Float *)f_out[i], sol_multi[i], 0, 1);
	}
    }

    bd.freeFermion(src);
    for (int i = 0; i < Nshift; ++i) {
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

    SetBfmArg(cg_arg->mass);

    // does nothing other than setting sol to zero
    int f_size_cb = GJP.VolNodeSites() * SPINOR_SIZE * bd.Ls / 2;
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
    VRB.Func(cname, fname);

    if (cg_arg == NULL)
	ERR.Pointer(cname, fname, "cg_arg");
    int threads = omp_get_max_threads();

    SetBfmArg(cg_arg->mass);

    Fermion_t in[2] = { bd.allocFermion(), bd.allocFermion() };
    Fermion_t out[2] = { bd.allocFermion(), bd.allocFermion() };

    bd.residual = cg_arg->stop_rsd;
    bd.max_iter = bf.max_iter = cg_arg->max_num_iter;
    // FIXME: pass single precision rsd in a reasonable way.
    bf.residual = 1e-5;

    // deal with Mobius Dminus
    if (bd.solver == HmCayleyTanh) {
	bd.cps_impexFermion((Float *)f_in, out, 1);
#pragma omp parallel
	{
	    bd.G5D_Dminus(out, in, 0);
	}
    } else {
	bd.cps_impexFermion((Float *)f_in, in, 1);
    }

    bd.cps_impexFermion((Float *)f_out, out, 1);

    int iter = -1;

    if (madwf_arg_map.count(cg_arg->mass) > 0) {
	// MADWF inversion
	VRB.Result(cname, fname, "Using MADWF: Main Ls = %d, cheap approx Ls = %d.\n", bd.Ls, madwf_arg_map[cg_arg->mass].cheap_approx.Ls);

	iter = MADWF_CG_M(bd, bf, use_mixed_solver,
	    out, in, bd.mass, this->GaugeField(), cg_arg->stop_rsd, madwf_arg_map[cg_arg->mass], cg_arg->Inverter);
    } else if (cg_arg->Inverter == HDCG) {
	HDCG_wrapper *control = HDCGInstance::getInstance();
	assert(control != NULL);
	control->HDCG_set_mass(cg_arg->mass);
	control->HDCG_invert(out, in, cg_arg->stop_rsd, cg_arg->max_num_iter);
    } else {
	// no MADWF:
#pragma omp parallel
	{
	    if (use_mixed_solver) {
		iter = mixed_cg::threaded_cg_mixed_M(out, in, bd, bf, 5, cg_arg->Inverter, evec, evalf, ecnt);
	    } else {
		switch (cg_arg->Inverter) {
		    case CG:
			if (evec && evald && ecnt) {
			    iter = bd.CGNE_M(out, in, *evec, *evald);
			} else {
			    iter = bd.CGNE_M(out, in);
			}
			break;
		    case EIGCG:
			iter = bd.EIG_CGNE_M(out, in);
			break;
		    default:
			if (bd.isBoss()) {
			    printf("%s::%s: Not implemented\n", cname, fname);
			}
			exit(-1);
			break;
		}
	    }
	}
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
    VRB.Result(cname, fname, "Taking Ls from current_key_mass = %e!\n", current_key_mass);
    const int size_5d = size_4d * arg_map.at(current_key_mass).Ls; // current_key_mass must be set correctly!!!

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
    
    SetBfmArg(eig_arg->mass);
    bd.residual = eig_arg->Rsdlam;
    bd.max_iter = eig_arg->MaxCG;

    VRB.Result(cname, fname, "residual = %17.10e max_iter = %d mass = %17.10e\n",
               bd.residual, bd.max_iter, bd.mass);
#if 0
    if( eig_arg->RitzMatOper == MATPCDAG_MATPC) 
{
    Fermion_t x[2];
    x[0] = bd.allocFermion();
    x[1] = bd.allocFermion();
    LatVector random0(4*GJP.SnodeSites());
    LatVector random1(4*GJP.SnodeSites());
    RandGaussVector(random0.Vec(),0.5,1);
    RandGaussVector(random1.Vec(),0.5,1);
//    Float * f_tmp = (Float *)f_eigenv[0];
    bd.cps_impexcbFermion(random0.Field(), x[0], 1, 1);
//    f_tmp += (GJP.VolNodeSites()/2)*(4*3*2);// checkerboarded 4D volume
    bd.cps_impexcbFermion(random1.Field(), x[1], 1, 1);

#pragma omp parallel
    {
        lambda[0] = bd.simple_lanczos(x);
    }

    bd.freeFermion(x[0]);
    bd.freeFermion(x[1]);
}
#endif

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

    SetBfmArg(mass);

    MatPc(phi, frm1, mass, dag);
    Float ret = FhamiltonNode(frm1, frm1);
    return ret;
}

void Fbfm::MatPc(Vector *out, Vector *in, Float mass, DagType dag)
{
    const char *fname = "MatPc()";

    VRB.Result(cname, fname, "start MatPc: mass = %e\n", mass);

    SetBfmArg(mass);

    Fermion_t i = bd.allocFermion();
    Fermion_t o = bd.allocFermion();
    Fermion_t t = bd.allocFermion();

    bd.cps_impexcbFermion((Float *)in , i, 1, 1);
#pragma omp parallel
    {
        bd.Mprec(i, o, t, dag == DAG_YES, 0);
    }
    bd.cps_impexcbFermion((Float *)out, o, 0, 1);

    bd.freeFermion(i);
    bd.freeFermion(o);
    bd.freeFermion(t);

    VRB.Result(cname, fname, "end MatPc: mass = %e\n", mass);
}

// It evolves the canonical momentum mom by step_size
// using the fermion force.
ForceArg Fbfm::EvolveMomFforce(Matrix *mom, Vector *frm, 
                               Float mass, Float step_size)
{
    const char *fname = "EvolveMomFforce()";
  
    SetBfmArg(mass);

    const int f_size_4d = SPINOR_SIZE * GJP.VolNodeSites();
    const int f_size_cb = f_size_4d * bd.Ls / 2;
  
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
    static Timer time(cname, fname);
    time.start(true);

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

    time.stop(true);
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
                    StrOrdType from)
{
    const char *fname = "Fconvert()";

    // nothing needs to be done
    //ERR.NotImplemented(cname, fname);
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

#if 1
void Fbfm::ImportGauge()
{
    const char *fname="ImportGauge()";
    VRB.Result(cname,fname,"NEW VERSION with CPS parallel transport\n");
    LatMatrix One;
    LatMatrix LatDir[8];
    Matrix *mout[8],*min[8];
#if 1
  for(int i=0;i<One.Vol();i++){
    *(One.Mat(i))= 1.;
  }
#endif
    int dirs[8]={0,2,4,6,1,3,5,7};
  for(int i=0;i<8;i++){
    min[i] = One.Mat();
    mout[i] = LatDir[i].Mat();
  }
{
    ParTransGauge pt_g(*this);
//    pt_g.run(8,mout,min,dirs);
    pt_g.run(4,mout,min,dirs); //positive Dirs
    pt_g.run(4,mout+4,min+4,dirs+4); //positive Dirs
}
    bd.cps_importGauge_dir(LatDir[0].Field(),1); //Plus X
    bd.cps_importGauge_dir(LatDir[1].Field(),3); //Plus Y
    bd.cps_importGauge_dir(LatDir[2].Field(),5); //Plus Z
    bd.cps_importGauge_dir(LatDir[3].Field(),7); //Plus T
    bd.cps_importGauge_dir(LatDir[4].Field(),0); //Minus X
    bd.cps_importGauge_dir(LatDir[5].Field(),2); //Minus Y
    bd.cps_importGauge_dir(LatDir[6].Field(),4); //Minus Z
    bd.cps_importGauge_dir(LatDir[7].Field(),6); //Minus T
    if(use_mixed_solver) {
        bd.comm_end();
        bf.comm_init();
    bf.cps_importGauge_dir(LatDir[0].Field(),1); //Plus X
    bf.cps_importGauge_dir(LatDir[1].Field(),3); //Plus Y
    bf.cps_importGauge_dir(LatDir[2].Field(),5); //Plus Z
    bf.cps_importGauge_dir(LatDir[3].Field(),7); //Plus T
    bf.cps_importGauge_dir(LatDir[4].Field(),0); //Minus X
    bf.cps_importGauge_dir(LatDir[5].Field(),2); //Minus Y
    bf.cps_importGauge_dir(LatDir[6].Field(),4); //Minus Z
    bf.cps_importGauge_dir(LatDir[7].Field(),6); //Minus T
        bf.comm_end();
        bd.comm_init();
    }
}
#else
void Fbfm::ImportGauge()
{
    const char *fname="ImportGauge()";
    VRB.Result(cname,fname,"OLD VERSION with qpd++ parallel transport\n");
    Float *gauge = (Float *)(this->GaugeField());
    bd.cps_importGauge(gauge);
#if 0
if (0){
        bd.comm_end();
        kernel.comm_init();
        kernel.cps_importGauge(gauge);
        kernel.comm_end();
        bd.comm_init();
}
#endif
    if(use_mixed_solver) {
        bd.comm_end();
        bf.comm_init();
        bf.cps_importGauge(gauge);
        bf.comm_end();
        bd.comm_init();
    }
}
#endif

CPS_END_NAMESPACE

#endif
