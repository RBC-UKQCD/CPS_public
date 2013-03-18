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
#include "run_wsnk.h"
#include "run_k2pipi.h"

USING_NAMESPACE_CPS
using namespace std;

// The 8 different types of contractions related to I=2 K to pi pi are
// defined using Qi's convention. However 1-8 are relabelled as 0-7.
static void cal_8types(vector<vector<Rcomplex> > *tmp,
                       const WilsonMatrix p[2],
                       int t_glb)
{
    for(int mu = 0; mu < 4; ++mu) {
        WilsonMatrix l[2] = {p[0].glL(mu), p[1].glL(mu)};
        WilsonMatrix r[2] = {p[0].glR(mu), p[1].glR(mu)};

        (*tmp)[0][t_glb] += l[0].Trace() * l[1].Trace(); // 1
        (*tmp)[2][t_glb] += l[0].Trace() * r[1].Trace(); // 3
        (*tmp)[4][t_glb] += Trace(l[0], l[1]);           // 5
        (*tmp)[6][t_glb] += Trace(r[0], l[1]);           // 7

        (*tmp)[1][t_glb] += (SpinTrace(l[0]) * SpinTrace(l[1])).Tr();  // 2
        (*tmp)[3][t_glb] += (SpinTrace(l[0]) * SpinTrace(r[1])).Tr();  // 4
        (*tmp)[5][t_glb] += Trace(ColorTrace(l[0]), ColorTrace(l[1])); // 6
        (*tmp)[7][t_glb] += Trace(ColorTrace(r[0]), ColorTrace(l[1])); // 8
    }
}

////////////////////////////////////////////////////////////////////////
// Compute I=2 K to pipi correlation functions
//
// Note: since all d quarks are connected to the 4 point operator, we
// don't need its momentum.
void run_k2pipi(const AllProp &sprop,
                const AllProp &uprop,
                const AllProp &dprop,
                const char fn_t[], int traj,
                PROP_TYPE ptype)
{
    const int t_scale = ptype == PROP_PA ? 2 : 1;
    const int t_size = GJP.TnodeSites() * GJP.Tnodes();
    const int t_size_ap = t_scale * t_size;
    const int lcl[4] = {GJP.XnodeSites(), GJP.YnodeSites(),
                        GJP.ZnodeSites(), GJP.TnodeSites(),};
    const int lcl_vol = lcl[0] * lcl[1] * lcl[2] * lcl[3];
    const int shift = GJP.TnodeSites() * GJP.TnodeCoor();

    vector<vector<WilsonMatrix> > uwsnk;
    run_wall_snk(&uwsnk, uprop, ptype);

    // truncate files
    char buf[1024];
    sprintf(buf, fn_t, ptype_str[ptype], traj);
    FILE *fp = Fopen(buf, "w");

    for(unsigned sep = 0; sep < t_size_ap; ++sep) {
        for(unsigned tk = 0; tk < t_size; ++tk) {
            const unsigned tp = (tk + sep) % t_size_ap;
            if(sprop.empty(tk, ptype)) continue;
            if(uwsnk[tp].empty()) continue;
            if(dprop.empty(tp, ptype)) continue;

            vector<vector<Rcomplex> > kpp(8, vector<Rcomplex>(t_size_ap, Rcomplex(0, 0)));

#pragma omp parallel
            {
                // threaded results
                vector<vector<Rcomplex> > tmp(8, vector<Rcomplex>(t_size_ap, Rcomplex(0, 0)));

#pragma omp for
                for(int i = 0; i < t_scale * lcl_vol; ++i) {
                    int x[4];
                    compute_coord_ap(x, lcl, i, t_size);

                    WilsonMatrix p[3];
                    p[1] = sprop(tk, ptype)(i, ptype);
                    p[1].hconj();
                    p[2].glV(p[1], -5);
                    p[2].gr(-5);
                    p[1] = uwsnk[tp][tk];
                    p[1].hconj();

                    p[0] = dprop(tp, ptype)(i, ptype) * p[1] * p[2]; // part 0

                    p[2] = uprop(tp, ptype)(i, ptype);
                    p[2].hconj();
                    p[2].gr(-5);
                    p[1] = dprop(tp, ptype)(i, ptype) * p[2]; // part 1

                    // p[2] is not used in the following function.
                    cal_8types(&tmp, p, x[3] + shift);
                } // sites

#pragma omp critical
                for(int op_id = 0; op_id < 8; ++op_id) {
                    for(int t = 0; t < t_size_ap; ++t) {
                        kpp[op_id][t] += tmp[op_id][(t+tk)%t_size_ap];
                    }
                } // critical, for
            }//omp

            // FIXME
            assert(GJP.Snodes() == 1);
            for(int op_id = 0; op_id < 8; ++op_id) {
                QMP_sum_double_array((double *)kpp[op_id].data(), 2 * t_size_ap);
            }

            for(unsigned t = 0; t < t_size_ap; ++t) {
                Fprintf(fp, "%3u %3u %3u", sep, tk, t);
                for(int op_id = 0; op_id < 8; ++op_id) {
                    Fprintf(fp, " %17.10e %17.10e",
                            kpp[op_id][t].real(),
                            kpp[op_id][t].imag());
                }
                Fprintf(fp, "\n");
            }
        } //tk
    } //sep

    Fclose(fp);
}
