// -*- mode:c++; c-basic-offset:4 -*-
#include <vector>
#include <string>
#include <alg/qpropw.h>
#include <util/qcdio.h>
#include <util/time_cps.h>
#include "prop_container.h"
#include "run_wsnk.h"

CPS_START_NAMESPACE

// The 8 different types of contractions related to I=2 K to pi pi are
// defined using Qi's convention. However 1-8 are relabelled as 0-7.
//
// Contractions written in LR style (i.e., operators are
// gamma(i)(1+/-gamma(5))).
static inline void cal_8types_LR(Rcomplex *ret,
                                 const WilsonMatrix p[2],
                                 int t_glb)
{
    for(int mu = 0; mu < 4; ++mu) {
        WilsonMatrix l[2] = {p[0].glL(mu), p[1].glL(mu)};
        WilsonMatrix r[2] = {p[0].glR(mu), p[1].glR(mu)};

        ret[0] += l[0].Trace() * l[1].Trace(); // 1
        ret[2] += l[0].Trace() * r[1].Trace(); // 3
        ret[4] += Trace(l[0], l[1]);           // 5
        ret[6] += Trace(r[0], l[1]);           // 7

        ret[1] += (SpinTrace(l[0]) * SpinTrace(l[1])).Tr();  // 2
        ret[3] += (SpinTrace(l[0]) * SpinTrace(r[1])).Tr();  // 4
        ret[5] += Trace(ColorTrace(l[0]), ColorTrace(l[1])); // 6
        ret[7] += Trace(ColorTrace(r[0]), ColorTrace(l[1])); // 8
    }
}

// The 8 different types of contractions related to I=2 K to pi pi are
// defined using Qi's convention. However 1-8 are relabelled as 0-7.
//
// Contractions written in VA style (i.e., operators are gamma(i) and
// gamma(i)gamma(5)).
static inline void cal_8types_VA(Rcomplex *ret,
                                 const WilsonMatrix p[2],
                                 int t_glb)
{
    for(int mu = 0; mu < 4; ++mu) {
        WilsonMatrix v[2] = {p[0].glV(mu), p[1].glV(mu)};
        WilsonMatrix a[2] = {p[0].glA(mu), p[1].glA(mu)};

        ret[0] += v[0].Trace() * a[1].Trace(); // 1
        ret[2] += a[0].Trace() * v[1].Trace(); // 3
        ret[4] += Trace(v[0], a[1]);           // 5
        ret[6] += Trace(a[0], v[1]);           // 7

        ret[1] += (SpinTrace(v[0]) * SpinTrace(a[1])).Tr();  // 2
        ret[3] += (SpinTrace(a[0]) * SpinTrace(v[1])).Tr();  // 4
        ret[5] += Trace(ColorTrace(v[0]), ColorTrace(a[1])); // 6
        ret[7] += Trace(ColorTrace(a[0]), ColorTrace(v[1])); // 8
    }
}

// run I=2 K to pi pi.
void run_k2pipi(const AllProp &sprop,
                const AllProp &uprop,
                const AllProp &dprop,
                const std::string &fn,
                PROP_TYPE ptype)
{
    const char *fname = "run_k2pipi()";
    const int t_scale = ptype == PROP_PA ? 2 : 1;
    const int t_size = GJP.TnodeSites() * GJP.Tnodes();
    const int t_size_ap = t_scale * t_size;
    const int lcl[4] = {GJP.XnodeSites(), GJP.YnodeSites(),
                        GJP.ZnodeSites(), GJP.TnodeSites(),};
    const int lcl_vol = lcl[0] * lcl[1] * lcl[2] * lcl[3];
    const int shift = GJP.TnodeSites() * GJP.TnodeCoor();

    std::vector<std::vector<WilsonMatrix> > uwsnk;
    run_wall_snk(&uwsnk, uprop, ptype);

    const int total_ops = 8;
    const long total_size = (long)t_size_ap * t_size * t_size_ap * total_ops;
    std::vector<Rcomplex> kppv(total_size, 0);
    Rcomplex *kpp = kppv.data();

    Float dtime0 = dclock();

#pragma omp parallel for
    for(int sep_tk = 0; sep_tk < t_size_ap * t_size; ++sep_tk) {
        // for(int sep = 0; sep < t_size_ap; ++sep) {
        //     for(int tk = 0; tk < t_size; ++tk) {
        const int sep = sep_tk / t_size;
        const int tk = sep_tk % t_size;
        const int tp = (tk + sep) % t_size_ap;

        if(sprop.empty(tk, ptype)) continue;
        if(uwsnk[tp].empty()) continue;
        if(dprop.empty(tp, ptype)) continue;

        for(int i = 0; i < t_scale * lcl_vol; ++i) {
            int x[4];
            compute_coord_ap(x, lcl, i, t_size);
            int t_glb = x[3] + shift;
            int t_delta = (t_glb + t_size_ap - tk) % t_size_ap;
            const long loc = ((long)sep_tk * t_size_ap + t_delta) * total_ops;

            WilsonMatrix p[3];
            p[1] = sprop(i, tk, ptype);
            p[1].hconj();
            p[2].glV(p[1], -5);
            p[2].gr(-5);
            p[1] = uwsnk[tp][tk];
            p[1].hconj();

            p[0] = dprop(i, tp, ptype) * p[1] * p[2]; // part 0

            p[2] = uprop(i, tp, ptype);
            p[2].hconj();
            p[2].gr(-5);
            p[1] = dprop(i, tp, ptype) * p[2]; // part 1

            // p[2] is not used in the following function.
            cal_8types_VA(kpp + loc, p, x[3] + shift);
        } // sites
    } //sep_tk

    // FIXME
    assert(GJP.Snodes() == 1);
    QMP_sum_double_array((double *)kpp, total_size * 2);

    Float dtime1 = dclock();

    // binary dump
    {
        const string fn_bin = fn + ".bin";
        FILE *fp = Fopen(fn_bin.c_str(), "w");
        Fwrite(kpp, total_size * sizeof(Rcomplex), 1, fp);
        Fclose(fp);
    }

    //output
    const int n_nodes = GJP.Xnodes() * GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();

    for(int sep = UniqueID(); sep < t_size_ap; sep += n_nodes) {
        const string fn_node = fn + tostring(".") + tostring(sep);
        FILE *fp = fopen(fn_node.c_str(), "w");

        for(int tk = 0; tk < t_size; ++tk) {
            const int tp = (tk + sep) % t_size_ap;
            const long loc = (long)(sep * t_size + tk) * t_size_ap * total_ops;

            if(sprop.empty(tk, ptype)) continue;
            if(uwsnk[tp].empty()) continue;
            if(dprop.empty(tp, ptype)) continue;
            for(int t = 0; t < t_size_ap; ++t) {
                fprintf(fp, "%3d %3d %3d", sep, tk, t);

                for(int op_id = 0; op_id < total_ops; ++op_id) {
                    fprintf(fp, " %17.10e %17.10e",
                            real(kpp[loc + t * total_ops + op_id]),
                            imag(kpp[loc + t * total_ops + op_id]));
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
