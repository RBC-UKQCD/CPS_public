// -*- mode:c++;c-basic-offset:4 -*-

#include <util/gjp.h>
#include <util/error.h>
#include <util/vector.h>
#include <util/omp_wrapper.h>
#include <util/sproj_tr.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <alg/force_arg.h>
#include <util/wilson.h>
#include <util/lattice/fforce_wilson_type.h>


USING_NAMESPACE_CPS;

static void compute_coord(long x[4], const long hl[4], const long low[4], long i)
{
    x[0] = i % hl[0] + low[0]; i /= hl[0];
    x[1] = i % hl[1] + low[1]; i /= hl[1];
    x[2] = i % hl[2] + low[2]; i /= hl[2];
    x[3] = i % hl[3] + low[3];
}

static long idx_4d(const long x[4], const long lx[4])
{
    long ret = 0;
    for(int i = 3; i >= 0; --i) {
        ret = ret * lx[i] + x[i];
    }
    return ret;
}

static long idx_4d_surf(const long x[4], const long lx[4], int mu)
{
    long ret = 0;
    for(int i = 3; i >= 0; --i) {
        if(i == mu) continue;
        ret = ret * lx[i] + x[i];
    }
    return ret;
}

static void thread_work_partial(long nwork, int me, int nthreads,
                                long &mywork, long &myoff)
{
    long basework = nwork / nthreads;
    long backfill = nthreads - (nwork % nthreads);
    mywork = (nwork + me) / nthreads;
    myoff  = basework * me;
    if ( me > backfill ) 
        myoff += (me-backfill);
}

void updateForce(ForceArg *f_arg, const Matrix &m)
{
    Float a2 = m.norm();
    Float a = sqrt(a2);
    
    f_arg->L1 += a;
    f_arg->L2 += a2;
    f_arg->Linf = f_arg->Linf > a ? f_arg->Linf : a;
}

// Calculate fermion force on a specific site, also do the
// summation over s direction (if s > 1).
static void do_site_force(Matrix *force, const Matrix &gauge,
                          Float *v1, Float *v1p,
                          Float *v2, Float *v2p,
                          int mu, int ls)
{
    Matrix t1, t2;
    
    Float *t1f = (Float *)&t1;
    Float *t2f = (Float *)&t2;

    // sproj_tr[mu](   (Float *)&t1, v1p, v2, ls, 0, 0);
    // sproj_tr[mu+4]( (Float *)&t2, v2p, v1, ls, 0, 0);
    
    switch(mu) {
    case 0:
      sprojTrXm(t1f, v1p, v2, ls, 0, 0);
      sprojTrXp(t2f, v2p, v1, ls, 0, 0);
      break;
    case 1:
      sprojTrYm(t1f, v1p, v2, ls, 0, 0);
      sprojTrYp(t2f, v2p, v1, ls, 0, 0);
      break;
    case 2:
      sprojTrZm(t1f, v1p, v2, ls, 0, 0);
      sprojTrZp(t2f, v2p, v1, ls, 0, 0);
      break;
    default:
      sprojTrTm(t1f, v1p, v2, ls, 0, 0);
      sprojTrTp(t2f, v2p, v1, ls, 0, 0);
    }

    t1 += t2;
    force->DotMEqual(gauge, t1);
}

FforceWilsonType::FforceWilsonType(Matrix *momentum, Matrix *gauge_field,
                                   Float *vec1, Float *vec2, int Ls, Float dt)
{
    lcl[0] = GJP.XnodeSites();
    lcl[1] = GJP.YnodeSites();
    lcl[2] = GJP.ZnodeSites();
    lcl[3] = GJP.TnodeSites();
    // This has to be passed in as a parameter (say, we are using 4D
    // fermions together with 5D fermions, then GJP.SnodeSites() is
    // not reliable).
    lcl[4] = Ls;

    mom = momentum;
    gauge = gauge_field;
    v1 = vec1;
    v2 = vec2;
    coef = dt;

    long vol_5d = lcl[0] * lcl[1] * lcl[2] * lcl[3] * lcl[4];
    bufsize = 0;
    for(int mu = 0; mu < 4; ++mu) {
        surf_size[mu] = SPINOR_SIZE * (vol_5d / lcl[mu]);
        v1so[mu] = bufsize;
        v2so[mu] = bufsize + surf_size[mu];
        surf_size[mu] *= 2;
        bufsize += surf_size[mu];
    }

    sndbuf = new Float[bufsize];
    rcvbuf = new Float[bufsize];

    if(sndbuf == NULL || rcvbuf == NULL) {
        ERR.Pointer("FforceWilsonType", "FforceWilsonType", "snd/rcv buf");
    }
}

FforceWilsonType::~FforceWilsonType()
{
    delete[] sndbuf;
    delete[] rcvbuf;
}

// copy 3d surface data from v4d to v3d in mu direction, use this
// function to fill the buffer v4d before communication.
//
// NOTE:
//
// 1. v4d is assumed to be in sxyzt order, i.e., the s index
// changes fastest.
//
// 2. we always send in negative direction.
void FforceWilsonType::collect_surface(int mu)
{
    long l[4] = { 0, 0, 0, 0 };
    long h[4] = { lcl[0], lcl[1], lcl[2], lcl[3] };
    h[mu] = 1;

    const long block = SPINOR_SIZE * lcl[4];
    const long sites = h[0] * h[1] * h[2] * h[3];

    Float *v1s = sndbuf + v1so[mu];
    Float *v2s = sndbuf + v2so[mu];

#pragma omp parallel for 
    for(long i = 0; i < sites; ++i) {
        long x[4];
        compute_coord(x, h, l, i);
        long o4d = idx_4d(x, lcl);
        long o3d = idx_4d_surf(x, lcl, mu);
        
        memcpy(v1s + o3d * block, v1  + o4d * block,
               sizeof(Float) * block);
        memcpy(v2s + o3d * block, v2  + o4d * block,
               sizeof(Float) * block);
    }
}

void FforceWilsonType::comm(int mu)
{
    getPlusData(rcvbuf + v1so[mu], sndbuf + v1so[mu],
                surf_size[mu], mu);
}

ForceArg FforceWilsonType::do_internal(int mu, int nthreads)
{
    long low[4] = { 0, 0, 0, 0 };
    long high[4] = { lcl[0], lcl[1], lcl[2], lcl[3] };
    --high[mu];
    const long sites = high[0] * high[1] * high[2] * high[3];
    const long block = SPINOR_SIZE * lcl[4];

    int me = omp_get_thread_num();
    long mywork, myoff;
    // some threads are used in communication
    thread_work_partial(sites, me, nthreads, mywork, myoff);

    ForceArg ret;
    for(long i = 0; i < mywork; ++i) {
        long x[4];
        compute_coord(x, high, low, i + myoff);
        long off_4d = idx_4d(x, lcl);
        long gid = mu + 4 * off_4d;
        long fid = block * off_4d;

        long y[4] = {x[0], x[1], x[2], x[3]};
        ++y[mu];
        long fidp = block * idx_4d(y, lcl);

        Matrix force;
        do_site_force(&force, gauge[gid],
                      v2 + fid, v2 + fidp,
                      v1 + fid, v1 + fidp, mu, lcl[4]);

        force.TrLessAntiHermMatrix();
        force *= -coef;
        *(mom + gid) += force;
        updateForce(&ret, force);
    }

    return ret;
}

ForceArg FforceWilsonType::do_surface(int mu, int nthreads)
{
    long low[4] = { 0, 0, 0, 0 };
    long high[4] = { lcl[0], lcl[1], lcl[2], lcl[3] };
    low[mu] = lcl[mu] - 1;
    high[mu] = lcl[mu];
    const long hl[4] = {high[0] - low[0],
                        high[1] - low[1],
                        high[2] - low[2],
                        high[3] - low[3] };
    const long sites = hl[0] * hl[1] * hl[2] * hl[3];
    const long block = SPINOR_SIZE * lcl[4];

    int me = omp_get_thread_num();
    long mywork, myoff;
    thread_work_partial(sites, me, nthreads, mywork, myoff);

    Float *v1_s = rcvbuf + v1so[mu];
    Float *v2_s = rcvbuf + v2so[mu];

    ForceArg ret;
    for(long i = 0; i < mywork; ++i) {
        long x[4];
        compute_coord(x, hl, low, i + myoff);

        long off_4d = idx_4d(x, lcl);
        long gid = mu + 4 * off_4d;
        long fid = block * off_4d;
        long fid_s = block * idx_4d_surf(x, lcl, mu);

        Matrix force;
        do_site_force(&force, gauge[gid],
                      v2 + fid, v2_s + fid_s,
                      v1 + fid, v1_s + fid_s, mu, lcl[4]);

        force.TrLessAntiHermMatrix();
        force *= -coef;
        *(mom + gid) += force;
        updateForce(&ret, force);
    }
    return ret;
}

ForceArg FforceWilsonType::run()
{
    for(int mu = 0; mu < 4; ++mu) {
        collect_surface(mu);
    }

    // single threaded comm
    for(int mu = 0; mu < 4; ++mu) {
        comm(mu);
    }

    ForceArg ret;

#pragma omp parallel
    {
        int nthreads = omp_get_num_threads();
        ForceArg f_arg; // threaded

        for(int mu = 0; mu < 4; ++mu) {
            f_arg.combine(do_internal(mu, nthreads));
            f_arg.combine(do_surface(mu, nthreads));
        }

#pragma omp critical
        {
            ret.combine(f_arg);
        }
    }

    ret.glb_reduce();
    return ret;
}
