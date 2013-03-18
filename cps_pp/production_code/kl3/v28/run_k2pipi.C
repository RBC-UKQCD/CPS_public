// -*- mode:c++; c-basic-offset:4 -*-
#include <vector>
#include <string>
#include <alg/qpropw.h>
#include <util/qcdio.h>
#include "prop_container.h"
#include "run_wsnk.h"

CPS_START_NAMESPACE

// The 8 different types of contractions related to I=2 K to pi pi are
// defined using Qi's convention. However 1-8 are relabelled as 0-7.
//
// Contractions written in LR style (i.e., operators are
// gamma(i)(1+/-gamma(5))).
static inline void cal_8types_LR(std::vector<std::vector<Rcomplex> > *tmp,
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

// The 8 different types of contractions related to I=2 K to pi pi are
// defined using Qi's convention. However 1-8 are relabelled as 0-7.
//
// Contractions written in VA style (i.e., operators are gamma(i) and
// gamma(i)gamma(5)).
static inline void cal_8types_VA(std::vector<std::vector<Rcomplex> > *tmp,
                                 const WilsonMatrix p[2],
                                 int t_glb)
{
    for(int mu = 0; mu < 4; ++mu) {
        WilsonMatrix v[2] = {p[0].glV(mu), p[1].glV(mu)};
        WilsonMatrix a[2] = {p[0].glA(mu), p[1].glA(mu)};

        (*tmp)[0][t_glb] += v[0].Trace() * a[1].Trace(); // 1
        (*tmp)[2][t_glb] += a[0].Trace() * v[1].Trace(); // 3
        (*tmp)[4][t_glb] += Trace(v[0], a[1]);           // 5
        (*tmp)[6][t_glb] += Trace(a[0], v[1]);           // 7

        (*tmp)[1][t_glb] += (SpinTrace(v[0]) * SpinTrace(a[1])).Tr();  // 2
        (*tmp)[3][t_glb] += (SpinTrace(a[0]) * SpinTrace(v[1])).Tr();  // 4
        (*tmp)[5][t_glb] += Trace(ColorTrace(v[0]), ColorTrace(a[1])); // 6
        (*tmp)[7][t_glb] += Trace(ColorTrace(a[0]), ColorTrace(v[1])); // 8
    }
}

// run I=2 K to pi pi.
void run_k2pipi(const AllProp &sprop,
                const AllProp &uprop,
                const AllProp &dprop,
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

    std::vector<std::vector<WilsonMatrix> > uwsnk;
    run_wall_snk(&uwsnk, uprop, ptype);

    FILE *fp = Fopen(fn.c_str(), "w");

    for(unsigned sep = 0; sep < t_size_ap; ++sep) {
        for(unsigned tk = 0; tk < t_size; ++tk) {
            const unsigned tp = (tk + sep) % t_size_ap;
            if(sprop.empty(tk, ptype)) continue;
            if(uwsnk[tp].empty()) continue;
            if(dprop.empty(tp, ptype)) continue;

            std::vector<std::vector<Rcomplex> >
                kpp(8, std::vector<Rcomplex>(t_size_ap, Rcomplex(0, 0)));

#pragma omp parallel
            {
                // threaded results
                std::vector<std::vector<Rcomplex> >
                    tmp(8, std::vector<Rcomplex>(t_size_ap, Rcomplex(0, 0)));

#pragma omp for
                for(int i = 0; i < t_scale * lcl_vol; ++i) {
                    int x[4];
                    compute_coord_ap(x, lcl, i, t_size);

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
                    cal_8types_VA(&tmp, p, x[3] + shift);
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

CPS_END_NAMESPACE
