// -*- mode:c++; c-basic-offset:4 -*-
#include <cstring>

#include <qmp.h>
#include <util/gjp.h>
#include <util/rcomplex.h>
#include <util/qcdio.h>
#include <alg/qpropw.h>
#include <vector>
#include <cassert>

#include <omp.h>

#include "my_util.h"
#include "run_pion.h"

USING_NAMESPACE_CPS
using namespace std;

// compute mres, note: no need to use P+A here, we simply store mres
// computed from all propagators we have.
void run_mres(const cps::QPropW &qp,
              unsigned tsrc, const char fn[])
{
    if(qp.StoreMidprop() == 0) return; 

    const int t_size_glb = GJP.TnodeSites() * GJP.Tnodes();

    const int lcl[4] = {GJP.XnodeSites(), GJP.YnodeSites(),
                        GJP.ZnodeSites(), GJP.TnodeSites(),};
    const int lcl_vol = lcl[0] * lcl[1] * lcl[2] * lcl[3];
    const int shift = GJP.TnodeSites() * GJP.TnodeCoor();

    vector<Rcomplex> pion(t_size_glb, Rcomplex(0, 0));
    vector<Rcomplex> j5_q(t_size_glb, Rcomplex(0, 0));

#pragma omp parallel
    {
        // threaded results
        vector<Rcomplex> pion_tmp(t_size_glb, Rcomplex(0, 0));
        vector<Rcomplex> j5_q_tmp(t_size_glb, Rcomplex(0, 0));

#pragma omp for
        for(int i = 0; i < lcl_vol; ++i) {
            int x[4];
            compute_coord(x, lcl, i);
            int t_glb = x[3] + shift;
            
            // J5 contraction (pion)
            WilsonMatrix p[2]  = {qp[i], qp[i]};
            p[1].hconj();
            // J5q contraction (midplane)
            WilsonMatrix q[2]  = {qp(i), qp(i)};
            q[1].hconj();
            
            pion_tmp[t_glb] += Trace(p[0], p[1]);
            j5_q_tmp[t_glb] += Trace(q[0], q[1]);
        } // sites
#pragma omp critical
        for(int t = 0; t < t_size_glb; ++t) {
            pion[t] += pion_tmp[(t+tsrc)%t_size_glb];
            j5_q[t] += j5_q_tmp[(t+tsrc)%t_size_glb];
        } // critical, for
    }//omp

    // FIXME
    assert(GJP.Snodes() == 1);
    QMP_sum_double_array((double *)pion.data(), 2 * t_size_glb);
    QMP_sum_double_array((double *)j5_q.data(), 2 * t_size_glb);

    FILE *fp = Fopen(fn, "a+");
    for(unsigned t = 0; t < t_size_glb; ++t) {
        Fprintf(fp, "%3u %3u %17.10e %17.10e %17.10e %17.10e\n", tsrc, t,
                pion[t].real(), pion[t].imag(),
                j5_q[t].real(), j5_q[t].imag());
    }
    Fprintf(fp, "\n");
    Fclose(fp);
}

// Compute Z_A, axial current normalization factor. Note: no need to
// use P+A here, we simply store the results computed from all
// propagators we have.
//
// Does not work yet.
void run_za(const cps::QPropW &qp, Float mass, 
            unsigned tsrc, const char fn[])
{
    if(qp.StoreMidprop() == 0) return; 

    const int t_size_glb = GJP.TnodeSites() * GJP.Tnodes();

    const int lcl[4] = {GJP.XnodeSites(), GJP.YnodeSites(),
                        GJP.ZnodeSites(), GJP.TnodeSites(),};
    const int lcl_vol = lcl[0] * lcl[1] * lcl[2] * lcl[3];
    const int shift = GJP.TnodeSites() * GJP.TnodeCoor();

    vector<Rcomplex> a4(t_size_glb, Rcomplex(0, 0));
    vector<Rcomplex> a5(t_size_glb, Rcomplex(0, 0));

#pragma omp parallel
    {
        // threaded results
        vector<Rcomplex> a4_tmp(t_size_glb, Rcomplex(0, 0));
        vector<Rcomplex> a5_tmp(t_size_glb, Rcomplex(0, 0));

#pragma omp for
        for(int i = 0; i < lcl_vol; ++i) {
            int x[4];
            compute_coord(x, lcl, i);
            int t_glb = x[3] + shift;
            
            // a4t contraction
            WilsonMatrix a[2] = {qp[i].glV(3), qp[i]};
            a[1].hconj();

            // J5 contraction (pion)
            WilsonMatrix p[2]  = {qp[i], qp[i]};
            p[1].hconj();
            // J5q contraction (midplane)
            WilsonMatrix q[2]  = {qp(i), qp(i)};
            q[1].hconj();
            
            a4_tmp[t_glb] += -Trace(a[0], a[1]);

            a5_tmp[t_glb] += 2 * mass * Trace(p[0], p[1]);
            a5_tmp[t_glb] += 2 * Trace(q[0], q[1]);
        } // sites
#pragma omp critical
        for(int t = 0; t < t_size_glb; ++t) {
            a4[t] += a4_tmp[(t+tsrc)%t_size_glb];
            a5[t] += a5_tmp[(t+tsrc)%t_size_glb];
        } // critical, for
    }//omp

    // FIXME
    assert(GJP.Snodes() == 1);
    QMP_sum_double_array((double *)a4.data(), 2 * t_size_glb);
    QMP_sum_double_array((double *)a5.data(), 2 * t_size_glb);

    FILE *fp = Fopen(fn, "a+");
    for(unsigned t = 0; t < t_size_glb; ++t) {
        Fprintf(fp, "%3u %3u %17.10e %17.10e %17.10e %17.10e\n",
                tsrc, t,
                a4[t].real() - a4[(t + t_size_glb - 1)%t_size_glb].real(),
                a4[t].imag() - a4[(t + t_size_glb - 1)%t_size_glb].imag(),
                a5[t].real(),
                a5[t].imag());
    }
    Fprintf(fp, "\n");
    Fclose(fp);
}
