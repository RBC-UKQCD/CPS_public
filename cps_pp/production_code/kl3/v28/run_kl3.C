// -*- mode:c++; c-basic-offset:4 -*-
#include <vector>
#include <string>
#include <alg/qpropw.h>
#include <util/qcdio.h>
#include "prop_container.h"
#include "run_wsnk.h"

CPS_START_NAMESPACE

////////////////////////////////////////////////////////////////////////
// Operators that we calculate for the kl3 contraction.
//
// Make sure op_str and apply_op() are consistent!
//
// We only do the twisting in x direction, so we put gamma_x and
// gamma_t here. The others are statistically zero.

// 5 contractions:
// 0: gamma_0 (x)
// 1: gamma_1 (y)
// 2: gamma_2 (z)
// 3: gamma_3 (t)
// 4: I       (identity matrix)
const int total_ops = 5;
static const char *op_str[total_ops] = {
    "0", "1", "2", "3", "I",
};

static const WilsonMatrix &apply_op(WilsonMatrix &out,
                                    const WilsonMatrix &in,
                                    int op_id)
{
    switch(op_id) {
    case 0:
    case 1:
    case 2:
    case 3:
        return out.glV(in, op_id);
    case 4:
        return out = in;
    default:
        ERR.General("", "apply_op", "Invalid op_id = %d\n", op_id);
        return out;
    }
}

// compute kl3 correlation functions, source type is
// determined implicitly by the propagator.
//
// This version computes all possible K-pi separations
void run_kl3(const AllProp &sprop,
             const AllProp &lprop,
             const AllProp &ltwst,
             const std::string &fn,
             PROP_TYPE ptype)
{
    const int t_scale = ptype == PROP_PA ? 2 : 1;
    const int t_size = GJP.TnodeSites() * GJP.Tnodes();
    const int t_size_ap = t_scale * t_size;
    const int lcl[4] = {GJP.XnodeSites(), GJP.YnodeSites(),
                        GJP.ZnodeSites(), GJP.TnodeSites(),};
    const int lcl_vol = lcl[0] * lcl[1] * lcl[2] * lcl[3];
    const int shift = GJP.TnodeSites() * GJP.TnodeCoor();

    std::vector<std::vector<WilsonMatrix> > lwsnk;
    run_wall_snk(&lwsnk, lprop, ptype);

    FILE *fp = Fopen(fn.c_str(), "w");

    for(unsigned sep = 0; sep < t_size_ap; ++sep) {
        for(unsigned tk = 0; tk < t_size; ++tk) {
            const unsigned tp = (tk + sep) % t_size_ap;
            if(sprop.empty(tk, ptype)) continue;
            if(lwsnk[tp].empty()) continue;
            if(ltwst.empty(tp, ptype)) continue;
            // 3 operators, see previous list
            std::vector<std::vector<Rcomplex> >
                kl3(total_ops, std::vector<Rcomplex>(t_size_ap, Rcomplex(0, 0)));

#pragma omp parallel
            {
                // threaded results
                std::vector<std::vector<Rcomplex> >
                    tmp(total_ops, std::vector<Rcomplex>(t_size_ap, Rcomplex(0, 0)));

#pragma omp for
                for(int i = 0; i < t_scale * lcl_vol; ++i) {
                    int x[4];
                    // Note: since i will never be larger than lcl_vol
                    // in P or A case, compute_coord_ap() suits all
                    // ptypes.
                    compute_coord_ap(x, lcl, i, t_size);
                    int t_glb = x[3] + shift;
            
                    WilsonMatrix p[3] = {0, ltwst(i, tp, ptype), 0};
                    p[0].glV(sprop(i, tk, ptype), -5);
                    p[2].glV(lwsnk[tp][tk], -5);
                    p[0].hconj();
                    p[2].hconj();

                    p[1] *= p[2] * p[0];

                    for(int op_id = 0; op_id < total_ops; ++op_id) {
                        apply_op(p[0], p[1], op_id);
                        tmp[op_id][t_glb] += p[0].Trace();
                    }
                } // sites

#pragma omp critical
                for(int op_id = 0; op_id < total_ops; ++op_id) {
                    for(int t = 0; t < t_size_ap; ++t) {
                        kl3[op_id][t] += tmp[op_id][(t+tk)%t_size_ap];
                    }
                } // critical, for
            }//omp

            // FIXME
            assert(GJP.Snodes() == 1);
            for(int op_id = 0; op_id < total_ops; ++op_id) {
                QMP_sum_double_array((double *)kl3[op_id].data(), 2 * t_size_ap);
            }

            for(unsigned t = 0; t < t_size_ap; ++t) {
                Fprintf(fp, "%3u %3u %3u", sep, tk, t);

                for(int op_id = 0; op_id < total_ops; ++op_id) {
                    Fprintf(fp, " %17.10e %17.10e",
                            real(kl3[op_id][t]),
                            imag(kl3[op_id][t]));
                }
                Fprintf(fp, "\n");
            }
        } //tk
    } //sep

    Fclose(fp);
}

CPS_END_NAMESPACE
