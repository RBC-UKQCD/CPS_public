// -*- mode:c++; c-basic-offset:4 -*-
#include <vector>
#include <string>
#include <alg/qpropw.h>
#include <util/qcdio.h>
#include <util/time_cps.h>
#include "my_util.h"
#include "prop_container.h"

CPS_START_NAMESPACE

// compute Bk correlation functions, source type is
// determined implicitly by the propagator.
void run_bk(const AllProp &lpropA,
            const AllProp &spropA,
            const AllProp &lpropB,
            const AllProp &spropB,
            const std::string &fn,
            PROP_TYPE ptype)
{
    const char *fname = "run_bk()";
    const int t_scale = ptype == PROP_PA ? 2 : 1;
    const int t_size = GJP.TnodeSites() * GJP.Tnodes();
    const int t_size_ap = t_scale * t_size;
    const int lcl[4] = {GJP.XnodeSites(), GJP.YnodeSites(),
                        GJP.ZnodeSites(), GJP.TnodeSites(),};
    const int lcl_vol = lcl[0] * lcl[1] * lcl[2] * lcl[3];
    const int shift = GJP.TnodeSites() * GJP.TnodeCoor();

    const int total_ops = 2;
    const long total_size = (long)t_size_ap * t_size * t_size_ap * total_ops;
    std::vector<Rcomplex> bkv(total_size, 0);
    Rcomplex *bk = bkv.data();

    Float dtime0 = dclock();

#pragma omp parallel for
    for(int sep_ka = 0; sep_ka < t_size_ap * t_size; ++sep_ka) {
        // for(int sep = 0; sep < t_size_ap; ++sep) {
        //     for(int ka = 0; ka < t_size; ++ka) {
        const int sep = sep_ka / t_size;
        const int ka = sep_ka % t_size;
        const int kb = (ka + sep) % t_size_ap;
        if(lpropA.empty(ka, ptype)) continue;
        if(spropA.empty(ka, ptype)) continue;
        if(lpropB.empty(kb, ptype)) continue;
        if(spropB.empty(kb, ptype)) continue;

        for(int i = 0; i < t_scale * lcl_vol; ++i) {
            int x[4];
            compute_coord_ap(x, lcl, i, t_size);
            int t_glb = x[3] + shift;
            int t_delta = (t_glb + t_size_ap - ka) % t_size_ap;
            long loc = ((long)sep_ka * t_size_ap + t_delta) * total_ops;

            // K & Kb
            WilsonMatrix K[2]  = {
                lpropA(i, ka, ptype), spropA(i, ka, ptype)
            };
            WilsonMatrix Kb[2] = {
                lpropB(i, kb, ptype), spropB(i, kb, ptype)
            };

            K[1].hconj();
            Kb[1].hconj();

            WilsonMatrix lines[2] = {K[0] * K[1], Kb[0] * Kb[1]};

            // type 0 = Tr(loop0) * Tr(loop1)
            // type 1 = Tr(loop0 * loop1)
            for(int mu = 0; mu < 4; ++mu) {
                WilsonMatrix tmp[2];
                // AA term
                tmp[0].glA(lines[0], mu);
                tmp[1].glA(lines[1], mu);

                bk[loc + 0] += tmp[0].Trace() * tmp[1].Trace();
                bk[loc + 1] += Trace(tmp[0], tmp[1]);

                // VV term
                tmp[0].glV(lines[0], mu);
                tmp[1].glV(lines[1], mu);

                bk[loc + 0] += tmp[0].Trace() * tmp[1].Trace();
                bk[loc + 1] += Trace(tmp[0], tmp[1]);
            } // mu
        } // sites
    } // sep, ka

    // FIXME
    assert(GJP.Snodes() == 1);
    QMP_sum_double_array((double *)bk, total_size * 2);

    Float dtime1 = dclock();

    // binary dump
    {
        const string fn_bin = fn + ".bin";
        FILE *fp = Fopen(fn_bin.c_str(), "w");
        Fwrite(bk, total_size * sizeof(Rcomplex), 1, fp);
        Fclose(fp);
    }

    // output
    const int n_nodes = GJP.Xnodes() * GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
    for(int sep = UniqueID(); sep < t_size_ap; sep += n_nodes) {
        const string fn_node = fn + tostring(".") + tostring(sep);
        FILE *fp = fopen(fn_node.c_str(), "w");

        for(int ka = 0; ka < t_size; ++ka) {
            const int kb = (ka + sep) % t_size_ap;
            const long loc = (long)(sep * t_size + ka) * t_size_ap * total_ops;

            if(lpropA.empty(ka, ptype)) continue;
            if(spropA.empty(ka, ptype)) continue;
            if(lpropB.empty(kb, ptype)) continue;
            if(spropB.empty(kb, ptype)) continue;

            for(int t = 0; t < t_size_ap; ++t) {
                fprintf(fp, "%3d %3d %3d", sep, ka, t);

                for(int op_id = 0; op_id < total_ops; ++op_id) {
                    fprintf(fp, " %17.10e %17.10e",
                            real(bk[loc + t * total_ops + op_id]),
                            imag(bk[loc + t * total_ops + op_id]));
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
