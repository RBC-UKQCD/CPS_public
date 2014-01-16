// -*- mode:c++; c-basic-offset:4 -*-
#include <vector>
#include <string>
#include <alg/qpropw.h>
#include <util/qcdio.h>
#include <util/time_cps.h>
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
    const char *fname = "run_kl3()";
    const int t_scale = ptype == PROP_PA ? 2 : 1;
    const int t_size = GJP.TnodeSites() * GJP.Tnodes();
    const int t_size_ap = t_scale * t_size;
    const int lcl[4] = {GJP.XnodeSites(), GJP.YnodeSites(),
                        GJP.ZnodeSites(), GJP.TnodeSites(),};
    const int lcl_vol = lcl[0] * lcl[1] * lcl[2] * lcl[3];
    const int shift = GJP.TnodeSites() * GJP.TnodeCoor();

    std::vector<std::vector<WilsonMatrix> > lwsnk;
    run_wall_snk(&lwsnk, lprop, ptype);

    const long total_size = (long)t_size_ap * t_size * t_size_ap * total_ops;
    std::vector<Rcomplex> kl3v(total_size, 0);
    Rcomplex *kl3 = kl3v.data();

    Float dtime0 = dclock();

#pragma omp parallel for
    for(int sep_tk = 0; sep_tk < t_size_ap * t_size; ++sep_tk) {
        // for(int sep = 0; sep < t_size_ap; ++sep) {
        //     for(int tk = 0; tk < t_size; ++tk) {
        const int sep = sep_tk / t_size;
        const int tk = sep_tk % t_size;
        const int tp = (tk + sep) % t_size_ap;
        const long loc = (long)sep_tk * t_size_ap * total_ops;
        if(sprop.empty(tk, ptype)) continue;
        if(lwsnk[tp].empty()) continue;
        if(ltwst.empty(tp, ptype)) continue;

        for(int i = 0; i < t_scale * lcl_vol; ++i) {
            int x[4];
            // Note: since i will never be larger than lcl_vol
            // in P or A case, compute_coord_ap() suits all
            // ptypes.
            compute_coord_ap(x, lcl, i, t_size);
            int t_glb = x[3] + shift;
            int t_delta = (t_glb + t_size_ap - tk) % t_size_ap;
            
            WilsonMatrix p[3] = {0, ltwst(i, tp, ptype), 0};
            p[0].glV(sprop(i, tk, ptype), -5);
            p[2].glV(lwsnk[tp][tk], -5);
            p[0].hconj();
            p[2].hconj();

            p[1] *= p[2] * p[0];

            for(int op_id = 0; op_id < total_ops; ++op_id) {
                apply_op(p[0], p[1], op_id);
                // tmp[op_id][t_glb] += p[0].Trace();
                kl3[loc + t_delta * total_ops + op_id] += p[0].Trace();
            }
        } // sites
    } //sep, tk

    // FIXME
    assert(GJP.Snodes() == 1);
    QMP_sum_double_array((double*)kl3, total_size * 2);

    Float dtime1 = dclock();

    // binary dump
    {
        const string fn_bin = fn + ".bin";
        FILE *fp = Fopen(fn_bin.c_str(), "w");
        Fwrite(kl3, total_size * sizeof(Rcomplex), 1, fp);
        Fclose(fp);
    }

    // output
    const int n_nodes = GJP.Xnodes() * GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();

    for(int sep = UniqueID(); sep < t_size_ap; sep += n_nodes) {
        const string fn_node = fn + tostring(".") + tostring(sep);
        FILE *fp = fopen(fn_node.c_str(), "w");

        for(int tk = 0; tk < t_size; ++tk) {
            const int tp = (tk + sep) % t_size_ap;
            const long loc = (long)(sep * t_size + tk) * t_size_ap * total_ops;

            if(sprop.empty(tk, ptype)) continue;
            if(lwsnk[tp].empty()) continue;
            if(ltwst.empty(tp, ptype)) continue;
            for(int t = 0; t < t_size_ap; ++t) {
                fprintf(fp, "%3d %3d %3d", sep, tk, t);

                for(int op_id = 0; op_id < total_ops; ++op_id) {
                    fprintf(fp, " %17.10e %17.10e",
                            real(kl3[loc + t * total_ops + op_id]),
                            imag(kl3[loc + t * total_ops + op_id]));
                }
                fprintf(fp, "\n");
            }
        }

        fclose(fp);
    }

    Float dtime2 = dclock();

    VRB.Result("", fname, "time1-0 = %17.10e seconds\n", dtime1  - dtime0);
    VRB.Result("", fname, "time2-1 = %17.10e seconds\n", dtime2  - dtime1);
    VRB.Result("", fname, "time2-0 = %17.10e seconds\n", dtime2  - dtime0);
}

CPS_END_NAMESPACE
