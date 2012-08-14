// -*- c-basic-offset: 4 -*-
#include<config.h>
#include<math.h>

#ifdef USE_BFM

#include <util/lattice/bfm_evo.h>
#include <util/lattice/fbfm.h>
#include <util/wilson.h>
#include <util/verbose.h>
#include <util/gjp.h>
#include <util/error.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <util/enum_func.h>
#include <util/sproj_tr.h>
#include<util/time_cps.h>

#include<omp.h>

CPS_START_NAMESPACE

inline void compute_coord(int x[4], const int hl[4], const int low[4], int i)
{
    x[0] = i % hl[0] + low[0]; i /= hl[0];
    x[1] = i % hl[1] + low[1]; i /= hl[1];
    x[2] = i % hl[2] + low[2]; i /= hl[2];
    x[3] = i % hl[3] + low[3];
}

// ----------------------------------------------------------------
// static void BondCond: toggle boundary condition on/off for any
// gauge-like field. Based on code from
// src/util/dirac_op/d_op_base/comsrc/dirac_op_base.C
//
// u_base must be in CANONICAL order.
// ----------------------------------------------------------------
template<typename Float>
static void BondCond(Float *u_base)
{
    for(int mu = 0; mu < 4; ++mu) {
        if(GJP.NodeBc(mu) != BND_CND_APRD) continue;

        int low[4] = { 0, 0, 0, 0 };
        int high[4] = { GJP.XnodeSites(), GJP.YnodeSites(),
                        GJP.ZnodeSites(), GJP.TnodeSites() };
        low[mu] = high[mu] - 1;

        int hl[4] = { high[0] - low[0], high[1] - low[1],
                      high[2] - low[2], high[3] - low[3] };

        const int hl_sites = hl[0] * hl[1] * hl[2] * hl[3];

#pragma omp parallel for
        for(int i = 0; i < hl_sites; ++i) {
            int x[4];
            compute_coord(x, hl, low, i);

            int off = mu + 4 * (x[0] + high[0] *
                                (x[1] + high[1] *
                                 (x[2] + high[2] * x[3])));
            Float *m = u_base + off * 18;
            for(int j = 0; j < 18; ++j) {
                m[j] = -m[j];
            }
        }
    }
}

bfmarg Fbfm::bfm_arg;
bool Fbfm::use_mixed_solver = 0;

// NOTE: Initialize QDP++ before using this class!
Fbfm::Fbfm(void):cname("Fbfm")
{
    const char *fname = "Fbfm()";

    if(GJP.Snodes() != 1) {
        ERR.NotImplemented(cname, fname);
    }

    bevo.init(bfm_arg);

    Float *gauge = (Float *)(this->GaugeField());
    BondCond(gauge);
    bevo.cps_importGauge(gauge);
    BondCond(gauge);

    // Fill in the array of sproj_tr functions, used for evolution.
    sproj_tr[SPROJ_XM] = sprojTrXm;
    sproj_tr[SPROJ_YM] = sprojTrYm;
    sproj_tr[SPROJ_ZM] = sprojTrZm;
    sproj_tr[SPROJ_TM] = sprojTrTm;
    sproj_tr[SPROJ_XP] = sprojTrXp;
    sproj_tr[SPROJ_YP] = sprojTrYp;
    sproj_tr[SPROJ_ZP] = sprojTrZp;
    sproj_tr[SPROJ_TP] = sprojTrTp;

    lclx[0] = GJP.XnodeSites();
    lclx[1] = GJP.YnodeSites();
    lclx[2] = GJP.ZnodeSites();
    lclx[3] = GJP.TnodeSites();
    lclx[4] = GJP.SnodeSites();

    int vol_5d = lclx[0] * lclx[1] * lclx[2] * lclx[3] * lclx[4];
    surf_size_all = 0;
    for(int i = 0; i < 4; ++i) {
        surf_size[i] = SPINOR_SIZE * (vol_5d / lclx[i]);
        surf_size_all += surf_size[i];
    }

    // calculate offset of surface vectors v1 and v2
    surf_v1[0] = 0;
    surf_v2[0] = surf_size[0];
    for(int i = 1; i < 4; ++i) {
        surf_v1[i] = surf_v1[i-1] + surf_size[i-1] * 2;
        surf_v2[i] = surf_v1[i] + surf_size[i];
    }
}

Fbfm::~Fbfm(void)
{
    bevo.end();
}

// copy 3d surface data from v4d to v3d in mu direction, use this
// function to fill the buffer v4d before communication.
//
// If send_neg == true, then it copies data on x[mu] == 0 surface to
// v3d (sends data in negative direction), otherwise it copies data on
// x[mu] == size - 1 surface to v3d.
//
// !!!NOTE: v4d is assumed to be in sxyzt order, i.e., the s index
// changes fastest.
void Fbfm::CopySendFrmData(Float *v3d, Float *v4d, int mu, bool send_neg)
{
    int low[4] = { 0, 0, 0, 0 };
    int high[4] = { lclx[0], lclx[1], lclx[2], lclx[3] };
    low[mu] = send_neg ? 0 : lclx[mu] - 1;
    high[mu] = low[mu] + 1;

    int block_size = SPINOR_SIZE * lclx[4]; // s inner most

    const int hl[4] = {high[0] - low[0],
                       high[1] - low[1],
                       high[2] - low[2],
                       high[3] - low[3] };
    const int hl_sites = hl[0] * hl[1] * hl[2] * hl[3];

#pragma omp parallel for 
    for(int i = 0; i < hl_sites; ++i) {
        int x[4];
        compute_coord(x, hl, low, i);
        int off_4d = idx_4d(x, lclx);
        int off_3d = idx_4d_surf(x, lclx, mu);
        
        memcpy(v3d + off_3d * block_size,
               v4d + off_4d * block_size,
               sizeof(Float) * block_size);
    }
}

// just TrLessAntiHermMatrix() ...
static inline void trless_am(Float *p, Float coef)
{
    p[0] = p[8] = p[16] = 0.;

    Float tmp = 0.5*(p[2] - p[6]) * coef;
    p[2]=tmp; p[6] = -tmp;

    tmp = 0.5*(p[3] + p[7]) * coef;
    p[3]=tmp; p[7] = tmp;

    tmp = 0.5*(p[4] - p[12]) * coef;
    p[4]=tmp; p[12] = -tmp;

    tmp = 0.5*(p[5] + p[13]) * coef;
    p[5]=tmp; p[13] = tmp;

    tmp = 0.5*(p[10] - p[14]) * coef;
    p[10]=tmp; p[14] = -tmp;

    tmp = 0.5*(p[11] + p[15]) * coef;
    p[11]=tmp; p[15] = tmp;

    IFloat c = 1./3. * (p[1] + p[9] + p[17]);

    p[1] = (p[1] - c) * coef;
    p[9] = (p[9] - c) * coef;
    p[17] = (p[17] - c) * coef;
}

static void thread_work_partial(int nwork, int me, int nthreads,
                                int &mywork, int &myoff)
{
    int basework = nwork / nthreads;
    int backfill = nthreads - (nwork % nthreads);
    mywork = (nwork + me) / nthreads;
    myoff  = basework * me;
    if ( me > backfill ) 
        myoff += (me-backfill);
}

ForceArg Fbfm::EvolveMomFforceInternal(Matrix *mom,
                                       Float *v1, Float *v2, // only internal data will be used
                                       Float coef, int mu,
                                       int nthreads)
{
    int low[4] = { 0, 0, 0, 0 };
    int high[4] = { lclx[0], lclx[1], lclx[2], lclx[3] };
    --high[mu];
    const int hl[4] = {high[0] - low[0],
                       high[1] - low[1],
                       high[2] - low[2],
                       high[3] - low[3] };
    const int hl_sites = hl[0] * hl[1] * hl[2] * hl[3];

    Matrix *gauge = GaugeField();

    int block_size = SPINOR_SIZE * lclx[4];

    int me = omp_get_thread_num();
    int mywork, myoff;
    // some threads are used in communication
    thread_work_partial(hl_sites, me, nthreads, mywork, myoff);

    ForceArg f_arg(0, 0, 0);
    for(int i = 0; i < mywork; ++i) {
        int x[4];
        compute_coord(x, hl, low, i + myoff);
        int off_4d = idx_4d(x, lclx);
        int gid = mu + 4 * off_4d;
        int fid = block_size * off_4d;

        int y[4] = {x[0], x[1], x[2], x[3]};
        ++y[mu];
        int fidp = block_size * idx_4d(y, lclx);

        Matrix force;
        FforceSiteS(force, gauge[gid],
                    v2 + fid, v2 + fidp,
                    v1 + fid, v1 + fidp, mu);
        trless_am((Float *)&force, -coef);
        // force.TrLessAntiHermMatrix();
        // force *= -coef;
        *(mom + gid) += force;
        updateForce(f_arg, force);
    }

    return f_arg;
}

ForceArg Fbfm::EvolveMomFforceSurface(Matrix *mom,
                                      Float *v1, Float *v2, // internal data
                                      Float *v1_s, Float *v2_s, // surface data
                                      Float coef, int mu)
{
    int low[4] = { 0, 0, 0, 0 };
    int high[4] = { lclx[0], lclx[1], lclx[2], lclx[3] };
    low[mu] = lclx[mu] - 1;
    high[mu] = lclx[mu];
    const int hl[4] = {high[0] - low[0],
                       high[1] - low[1],
                       high[2] - low[2],
                       high[3] - low[3] };
    const int hl_sites = hl[0] * hl[1] * hl[2] * hl[3];

    Matrix *gauge = GaugeField();

    int block_size = SPINOR_SIZE * lclx[4];
    int sign = GJP.NodeBc(mu) == BND_CND_APRD ? -1.0 : 1.0;

    int nthreads = omp_get_num_threads();
    int me = omp_get_thread_num();
    int mywork, myoff;

    // here all threads participate
    thread_work_partial(hl_sites, me, nthreads, mywork, myoff);

    ForceArg f_arg(0, 0, 0);
    for(int i = 0; i < mywork; ++i) {
        int x[4];
        compute_coord(x, hl, low, i + myoff);

        int off_4d = idx_4d(x, lclx);
        int gid = mu + 4 * off_4d;
        int fid = block_size * off_4d;
        int fid_s = block_size * idx_4d_surf(x, lclx, mu);

        Matrix force;
        FforceSiteS(force, gauge[gid],
                    v2 + fid, v2_s + fid_s,
                    v1 + fid, v1_s + fid_s, mu);
        trless_am((Float *)&force, -coef * sign);
        // force.TrLessAntiHermMatrix();
        // force *= -coef * sign;
        *(mom + gid) += force;
        updateForce(f_arg, force);
    }

    return f_arg;
}

// Calculate fermion force on a specific site, also do the
// summation over s direction.
void Fbfm::FforceSiteS(Matrix& force, Matrix &gauge,
                       Float *v1, Float *v1p,
                       Float *v2, Float *v2p, int mu)
{
    Matrix t1, t2;

    sproj_tr[mu](   (Float *)&t1, v1p, v2, lclx[4], 0, 0);
    sproj_tr[mu+4]( (Float *)&t2, v2p, v1, lclx[4], 0, 0);
    
    t1 += t2;

    force.DotMEqual(gauge, t1);
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
    Fermion_t pi[2] = {bevo.allocFermion(), bevo.allocFermion()};
    Fermion_t po[4] = {bevo.allocFermion(), bevo.allocFermion(),
                       bevo.allocFermion(), bevo.allocFermion()};

    bevo.mass = mass;
    // reinitialize since we are using a new mass.
    bevo.GeneralisedFiveDimEnd();
    bevo.GeneralisedFiveDimInit();

    bevo.cps_impexcbFermion((Float *)phi1, pi[0], 1, 1);
    bevo.cps_impexcbFermion((Float *)phi2, pi[1], 1, 1);

#pragma omp parallel
    {
        bevo.calcMDForceVecs(po + 0, po + 2, pi[0], pi[1]);
    }

    bevo.cps_impexFermion_s(v1, po + 0, 0);
    bevo.cps_impexFermion_s(v2, po + 2, 0);

    bevo.freeFermion(pi[0]);
    bevo.freeFermion(pi[1]);
    bevo.freeFermion(po[0]);
    bevo.freeFermion(po[1]);
    bevo.freeFermion(po[2]);
    bevo.freeFermion(po[3]);
}

ForceArg Fbfm::EvolveMomFforceBaseThreaded(Matrix *mom,
                                           Vector *phi1, Vector *phi2,
                                           Float mass, Float coef)
{
    const char *fname = "EvolveMomFforceBaseThreaded()";

    Float dtime = -dclock();

    Fermion_t in[2] = {bevo.allocFermion(), bevo.allocFermion()};

    bevo.mass = mass;
    // reinitialize since we are using a new mass.
    bevo.GeneralisedFiveDimEnd();
    bevo.GeneralisedFiveDimInit();

    Float *gauge = (Float *)(this->GaugeField());
    BondCond(gauge);

    bevo.cps_impexcbFermion((Float *)phi1, in[0], 1, 1);
    bevo.cps_impexcbFermion((Float *)phi2, in[1], 1, 1);

#pragma omp parallel
    {
        bevo.compute_force((Float *)mom, gauge, in[0], in[1], coef);
    }

    bevo.freeFermion(in[0]);
    bevo.freeFermion(in[1]);

    BondCond(gauge);
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

    Float dtime1 = dclock();

    int f_size_4d = SPINOR_SIZE * GJP.VolNodeSites();
    int f_size = f_size_4d * GJP.SnodeSites();
  
    Float *v1 = (Float *)smalloc(cname, fname, "v1", sizeof(Float) * f_size);
    Float *v2 = (Float *)smalloc(cname, fname, "v2", sizeof(Float) * f_size);
    Float *sndbuf = (Float *)smalloc(cname, fname, "sndbuf", sizeof(Float) * surf_size_all * 2);
    Float *rcvbuf = (Float *)smalloc(cname, fname, "rcvbuf", sizeof(Float) * surf_size_all * 2);

    CalcHmdForceVecsBilinear(v1, v2, phi1, phi2, mass);

    Float dtime2 = dclock();

    for(int i = 0; i < 4; ++i) {
        CopySendFrmData(sndbuf + surf_v1[i], v1, i, true);
        CopySendFrmData(sndbuf + surf_v2[i], v2, i, true);
    }

    Float dtime3 = dclock();

    // single threaded comm
    for(int dir = 0; dir < 4; ++dir) {
        getPlusData(rcvbuf + surf_v1[dir], sndbuf + surf_v1[dir],
                    surf_size[dir] * 2, dir);
    }

    omp_set_num_threads(bfm_arg.threads);
    ForceArg ret;

#pragma omp parallel
    {
        int nthreads = omp_get_num_threads();
        // int me = omp_get_thread_num();
        ForceArg f_arg; // threaded

        // internal forces
        for(int i = 0; i < 4; ++i) {
            ForceArg t = EvolveMomFforceInternal(mom, v1, v2, coef, i, nthreads);
            f_arg.combine(t);
        }

// #pragma omp barrier

        for(int i = 0; i < 4; ++i) {
            ForceArg t = EvolveMomFforceSurface(mom, v1, v2,
                                                rcvbuf + surf_v1[i], // v1 surface
                                                rcvbuf + surf_v2[i], // v2 surface
                                                coef, i);
            f_arg.combine(t);
        }

#pragma omp critical
        {
            ret.combine(f_arg);
        }
    }

    Float dtime4 = dclock();

    sfree(cname, fname, "v1", v1);
    sfree(cname, fname, "v2", v2);
    sfree(cname, fname, "sndbuf", sndbuf);
    sfree(cname, fname, "rcvbuf", rcvbuf);

    VRB.Result(cname, fname, "cal aux vectors  : takes %17.10e seconds\n", dtime2 - dtime1);
    VRB.Result(cname, fname, "prepare for comm : takes %17.10e seconds\n", dtime3 - dtime2);
    VRB.Result(cname, fname, "comm/forces      : takes %17.10e seconds\n", dtime4 - dtime3);
    VRB.Result(cname, fname, "total            : takes %17.10e seconds\n", dtime4 - dtime1);

    glb_sum(&ret.L1);
    glb_sum(&ret.L2);
    glb_max(&ret.Linf);

    ret.unitarize(4 * GJP.VolSites());
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

//------------------------------------------------------------------
// returns the type of fermion class
//------------------------------------------------------------------
FclassType Fbfm::Fclass(void)const
{
    return F_CLASS_BFM;
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

// Returns the number of fermion field 
// components (including real/imaginary) on a
// site of the 4-D lattice.
int Fbfm::FsiteSize(void)const
{
    return 24 * GJP.SnodeSites();
}

// Returns 0 => If no checkerboard is used for the evolution
//      or the CG that inverts the evolution matrix.
int Fbfm::FchkbEvl(void)const
{
    return 1;
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
    const char *fname = "FmatEvlInv(V*, V*, CgArg *, ...)";

    if(cg_arg == NULL)
        ERR.Pointer(cname, fname, "cg_arg");

    if(use_mixed_solver) {
        return FmatEvlInvMixed(f_out, f_in, cg_arg, 1e-5,
                               cg_arg->max_num_iter,
                               5);
    }

    Fermion_t in  = bevo.allocFermion();
    Fermion_t out = bevo.allocFermion();

    bevo.mass = cg_arg->mass;
    bevo.residual = cg_arg->stop_rsd;
    bevo.max_iter = cg_arg->max_num_iter;
    // reinitialize since we are using a new mass.
    bevo.GeneralisedFiveDimEnd();
    bevo.GeneralisedFiveDimInit();

    bevo.cps_impexcbFermion((Float *)f_in , in,  1, 1);
    bevo.cps_impexcbFermion((Float *)f_out, out, 1, 1);

    int iter;
#pragma omp parallel
    {
        iter = bevo.CGNE_prec_MdagM(out, in);
    }

    bevo.cps_impexcbFermion((Float *)f_out, out, 0, 1);

    bevo.freeFermion(in);
    bevo.freeFermion(out);

    return iter;
}

int Fbfm::FmatEvlInvMixed(Vector *f_out, Vector *f_in, 
                          CgArg *cg_arg,
                          Float single_rsd,
                          int max_iter,
                          int max_cycle)
{
    bfm_arg.mass = cg_arg->mass;

    bfm_evo<float> bfm_f;
    bfm_f.init(bfm_arg);
    bfm_f.residual = single_rsd;
    bfm_f.max_iter = max_iter;

    Float *gauge = (Float *)(this->GaugeField());
    BondCond(gauge);
    bfm_f.cps_importGauge(gauge);
    BondCond(gauge);

    bfm_f.comm_end();
    bevo.comm_init();

    bevo.mass = cg_arg->mass;
    bevo.residual = cg_arg->stop_rsd;
    bevo.max_iter = cg_arg->max_num_iter;
    // reinitialize since we are using a new mass.
    bevo.GeneralisedFiveDimEnd();
    bevo.GeneralisedFiveDimInit();

    Fermion_t src = bevo.allocFermion();
    Fermion_t sol = bevo.allocFermion();

    bevo.cps_impexcbFermion((Float *)f_in , src, 1, 1);
    bevo.cps_impexcbFermion((Float *)f_out, sol, 1, 1);

    int iter = -1;
#pragma omp parallel
    {
        iter = mixed_cg::threaded_cg_mixed_MdagM(sol, src, bevo, bfm_f, max_cycle);

        // bevo.max_iter = 20;
        // iter = mixed_cg::cg_MdagM_single_precnd(sol, src, bevo, bfm_f);
        // bevo.max_iter = cg_arg->max_num_iter;
    }

    bevo.comm_end();
    bfm_f.comm_init();
    bfm_f.end();
    bevo.comm_init();

    bevo.cps_impexcbFermion((Float *)f_out, sol, 0, 1);

    bevo.freeFermion(src);
    bevo.freeFermion(sol);

    return iter;
}

int Fbfm::FmatEvlInv(Vector *f_out, Vector *f_in, 
                     CgArg *cg_arg, 
                     CnvFrmType cnv_frm)
{
    return FmatEvlInv(f_out, f_in, cg_arg, NULL, cnv_frm);
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
        sol_multi[i] = bevo.allocFermion();
        ones[i] = 1.0;
        mresidual[i] = cg_arg[i]->stop_rsd;
    }

    // source
    Fermion_t src = bevo.allocFermion();
    bevo.cps_impexcbFermion((Float *)f_in, src, 1, 1);

    // reinitialize since we are using a new mass.
    bevo.mass = cg_arg[0]->mass;
    bevo.residual = cg_arg[0]->stop_rsd;
    bevo.max_iter = cg_arg[0]->max_num_iter;
    bevo.GeneralisedFiveDimEnd();
    bevo.GeneralisedFiveDimInit();

    int iter;
#pragma omp parallel
    {
        iter = bevo.CGNE_prec_MdagM_multi_shift(sol_multi, src, shift, ones, Nshift, mresidual, 0);
    }

    if(type == SINGLE) {
        if(1) {
            // FIXME: Never use this in production code!
            int f_size_cb = GJP.VolNodeSites() * SPINOR_SIZE * GJP.SnodeSites() / 2;
            Vector *t = (Vector *)smalloc(cname, fname, "t", sizeof(Float) * f_size_cb);

            for(int i = 0; i < Nshift; ++i) {
                bevo.cps_impexcbFermion((Float *)t, sol_multi[i], 0, 1);
                f_out[0]->FTimesV1PlusV2(alpha[i], t, f_out[0], f_size_cb);
            }
            sfree(cname, fname, "t", t);
        } else {
            // Can't do this since axpy needs a threaded environment (and also
            // other bfm related problems, mainly how it uses alpha now).
            // bevo.cps_impexcbFermion((Float *)f_out[0], bevo.psi_i[0], 1, 1);
            // for(int i = 0; i < Nshift; ++i) {
            //     bevo.axpy(bevo.psi_i[0],
            //               bevo.psi_multi[i],
            //               bevo.psi_i[0],
            //               alpha[i]);
            // }
            // bevo.axpy(bevo.psi_i[0],
            //           bevo.psi_multi[0],
            //           bevo.psi_i[0],
            //           1.0);
            // bevo.cps_impexcbFermion((Float *)f_out[0], bevo.psi_i[0], 0, 1);
        }
    } else {
        for(int i = 0; i < Nshift; ++i) {
            bevo.cps_impexcbFermion((Float *)f_out[i], sol_multi[i], 0, 1);
        }
    }

    bevo.freeFermion(src);
    for(int i = 0; i < Nshift; ++i) {
        bevo.freeFermion(sol_multi[i]);
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

    int f_size_cb = GJP.VolNodeSites() * SPINOR_SIZE * GJP.SnodeSites() / 2;

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

    if(cg_arg == NULL)
        ERR.Pointer(cname, fname, "cg_arg");

    Fermion_t in[2]  = {bevo.allocFermion(), bevo.allocFermion()};
    Fermion_t out[2] = {bevo.allocFermion(), bevo.allocFermion()};

    bevo.mass = cg_arg->mass;
    bevo.residual = cg_arg->stop_rsd;
    bevo.max_iter = cg_arg->max_num_iter;
    // reinitialize since we are using a new mass.
    bevo.GeneralisedFiveDimEnd();
    bevo.GeneralisedFiveDimInit();

    bevo.cps_impexFermion((Float *)f_in , in,  1);
    bevo.cps_impexFermion((Float *)f_out, out, 1);

    int iter;
#pragma omp parallel
    {
        iter = bevo.prop_solve(out, in);
    }

    bevo.cps_impexFermion((Float *)f_out, out, 0);

    bevo.freeFermion(in[0]);
    bevo.freeFermion(in[1]);
    bevo.freeFermion(out[0]);
    bevo.freeFermion(out[1]);

    return iter;
}

int Fbfm::FmatInv(Vector *f_out, Vector *f_in, 
                  CgArg *cg_arg, 
                  CnvFrmType cnv_frm,
                  PreserveType prs_f_in)
{
    return FmatInv(f_out, f_in, cg_arg, NULL, cnv_frm, prs_f_in);
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
    VRB.Func(cname,fname);

    int x;
    int i;
    Float *field_4D;
    Float *field_5D;

    //------------------------------------------------------------------
    // Initializations
    //------------------------------------------------------------------
    int f_size = GJP.VolNodeSites() * FsiteSize()*Ncb/2;
    int ls = GJP.SnodeSites();
    int vol_4d = GJP.VolNodeSites()*Ncb/2;
    int ls_stride = 24 * vol_4d;
    int s_u_local = s_u % GJP.SnodeSites();
    int s_l_local = s_l % GJP.SnodeSites();
    int s_u_node = s_u / GJP.SnodeSites();
    int s_l_node = s_l / GJP.SnodeSites();


    //------------------------------------------------------------------
    // Set *five using the 4D field *four. 
    //------------------------------------------------------------------

    // Set all components of the 5D field to zero.
    //---------------------------------------------------------------
    field_5D  = (Float *) five;
    for(i=0; i<f_size; i++){
        field_5D[i]  = 0.0;
    }

    // Do the two upper spin components if s_u is in the node
    //---------------------------------------------------------------
    if( s_u_node == GJP.SnodeCoor() ){
        field_4D  = (Float *) four;
        field_5D  = (Float *) five;
        field_5D  = field_5D  + s_u_local * ls_stride;
        for(x=0; x<vol_4d; x++){
            for(i=0; i<12; i++){
                field_5D[i]  = field_4D[i];
            }
            field_4D  = field_4D  + 24;
            field_5D  = field_5D  + 24;
        }
    }

    // Do the two lower spin components if s_l is in the node
    //----------------------------------------------------------------
    if( s_l_node == GJP.SnodeCoor() ){
        field_4D  = (Float *) four;
        field_5D  = (Float *) five;
        field_4D  = field_4D  + 12;
        field_5D  = field_5D  + 12 + s_l_local * ls_stride;
        for(x=0; x<vol_4d; x++){
            for(i=0; i<12; i++){
                field_5D[i]  = field_4D[i];
            }
            field_4D  = field_4D  + 24;
            field_5D  = field_5D  + 24;
        }
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

    int x;
    int i;
    Float *field_4D;
    Float *field_5D;
    VRB.Func(cname,fname);

    //------------------------------------------------------------------
    // Initializations
    //------------------------------------------------------------------
    int ls = GJP.SnodeSites();
    int f_size = GJP.VolNodeSites() * FsiteSize()*Ncb / (ls*2);
    int vol_4d = GJP.VolNodeSites()*Ncb/2;
    int ls_stride = 24 * vol_4d;
    int s_u_local = s_u % GJP.SnodeSites();
    int s_l_local = s_l % GJP.SnodeSites();
    int s_u_node = s_u / GJP.SnodeSites();
    int s_l_node = s_l / GJP.SnodeSites();


    //------------------------------------------------------------------
    // Set *four using the 5D field *five. 
    //------------------------------------------------------------------

    // Set all components of the 4D field to zero.
    //---------------------------------------------------------------
    field_4D  = (Float *) four;
    for(i=0; i<f_size; i++){
        field_4D[i]  = 0.0;
    }

    // Do the two upper spin components if s_u is in the node
    //---------------------------------------------------------------
    if( s_u_node == GJP.SnodeCoor() ){
        field_4D = (Float *) four;
        field_5D = (Float *) five;
        field_5D = field_5D + s_u_local * ls_stride;
        for(x=0; x<vol_4d; x++){
            for(i=0; i<12; i++){
                field_4D[i] = field_5D[i];
            }
            field_4D = field_4D + 24;
            field_5D = field_5D + 24;
        }
    }

    // Do the two lower spin components if s_l is in the node
    //----------------------------------------------------------------
    if( s_l_node == GJP.SnodeCoor() ){
        field_4D = (Float *) four;
        field_5D = (Float *) five;
        field_4D = field_4D + 12;
        field_5D = field_5D + 12 + s_l_local * ls_stride;
        for(x=0; x<vol_4d; x++){
            for(i=0; i<12; i++){
                field_4D[i] = field_5D[i];
            }
            field_4D = field_4D + 24;
            field_5D = field_5D + 24;
        }
    }

    // Sum along s direction to get the same 4D field in all 
    // s node slices.
    //----------------------------------------------------------------
    if( GJP.Snodes() > 1) {
        Float sum;
        field_4D  = (Float *) four;
        for(i=0; i<f_size; i++){
            sum = field_4D[i];
            glb_sum_dir(&sum, 4);
            field_4D[i] = sum;    
        }
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
    
    bevo.residual = eig_arg->Rsdlam;
    bevo.max_iter = eig_arg->MaxCG;
    bevo.mass = eig_arg->mass;
    // reinitialize since we are using a new mass.
    bevo.GeneralisedFiveDimEnd();
    bevo.GeneralisedFiveDimInit();

    VRB.Result(cname, fname, "residual = %17.10e max_iter = %d mass = %17.10e\n",
               bevo.residual, bevo.max_iter, bevo.mass);

    Fermion_t in = bevo.allocFermion();
    bevo.cps_impexcbFermion((Float *)f_eigenv[0], in, 1, 1);

#pragma omp parallel
    {
        lambda[0] = bevo.ritz(in, eig_arg->RitzMatOper == MATPCDAG_MATPC);
    }

    bevo.cps_impexcbFermion((Float *)f_eigenv[0], in, 0, 1);

    // correct the eigenvalue for a dumb convention problem.
    if(eig_arg->RitzMatOper == NEG_MATPCDAG_MATPC) lambda[0] = -lambda[0];

    valid_eig[0] = 1;
    bevo.freeFermion(in);

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

    Fermion_t i = bevo.allocFermion();
    Fermion_t o = bevo.allocFermion();
    Fermion_t t = bevo.allocFermion();

    bevo.mass = mass;
    // reinitialize since we are using a new mass.
    bevo.GeneralisedFiveDimEnd();
    bevo.GeneralisedFiveDimInit();

    bevo.cps_impexcbFermion((Float *)in , i, 1, 1);
#pragma omp parallel
    {
        bevo.Mprec(i, o, t, dag == DAG_YES, 0);
    }
    bevo.cps_impexcbFermion((Float *)out, o, 0, 1);

    bevo.freeFermion(i);
    bevo.freeFermion(o);
    bevo.freeFermion(t);
}

// It evolves the canonical momentum mom by step_size
// using the fermion force.
ForceArg Fbfm::EvolveMomFforce(Matrix *mom, Vector *frm, 
                               Float mass, Float step_size)
{
    const char *fname = "EvolveMomFforce()";
  
    const int f_size_4d = SPINOR_SIZE * GJP.VolNodeSites();
    const int f_size_cb = f_size_4d * GJP.SnodeSites() / 2;
  
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

    Float L1 = 0.0;
    Float L2 = 0.0;
    Float Linf = 0.0;

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

    // If measuring the force, need to measure and then sum to mom
    if (force_measure == FORCE_MEASURE_YES) {
        for (int i=0; i<g_size/18; i++) {
            Float norm = (mom_tmp+i)->norm();
            Float tmp = sqrt(norm);
            L1 += tmp;
            L2 += norm;
            Linf = (tmp>Linf ? tmp : Linf);
        }
        glb_sum(&L1);
        glb_sum(&L2);
        glb_max(&Linf);

        L1 /= 4.0*GJP.VolSites();
        L2 /= 4.0*GJP.VolSites();

        fTimesV1PlusV2((IFloat*)mom, 1.0, (IFloat*)mom_tmp, (IFloat*)mom, g_size);

        delete[] force_label;
        sfree(mom_tmp, cname, fname, "mom_tmp");
    }

    return ForceArg(L1, sqrt(L2), Linf);
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
void Fbfm::Fconvert(Vector *f_field,
                    StrOrdType to,
                    StrOrdType from)
{
    const char *fname = "Fconvert()";

    // nothing needs to be done
    //ERR.NotImplemented(cname, fname);
}


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

int Fbfm::SpinComponents()const
{
    return 4;
}

int Fbfm::ExactFlavors()const
{
    return 2;
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

CPS_END_NAMESPACE

#endif
