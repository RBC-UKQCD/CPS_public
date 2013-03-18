// -*- mode:c++; c-basic-offset:4 -*-
#include <vector>
#include <string>
#include <alg/qpropw.h>
#include <util/qcdio.h>
#include "prop_container.h"
#include "my_util.h"
#include "run_wsnk.h"

CPS_START_NAMESPACE

static inline const WilsonMatrix &apply_op_src(WilsonMatrix &out,
                                               const WilsonMatrix &in,
                                               Operator src_op)
{
    WilsonMatrix tmp;
    switch(src_op) {
    case GAMMA_0:
    case GAMMA_1:
    case GAMMA_2:
    case GAMMA_3:
        return out.glA(in, src_op - GAMMA_0);
    case GAMMA_5:
        return out = in;
    case ID:
        return out.glV(in, -5);
    case GAMMA_05:
    case GAMMA_15:
    case GAMMA_25:
    case GAMMA_35:
        return out.glV(in, src_op - GAMMA_05);
    case GAMMA_50:
    case GAMMA_51:
    case GAMMA_52:
    case GAMMA_53:
        tmp.glV(in, src_op - GAMMA_50);
        return out = -1*tmp;
    default:
        ERR.General("", "apply_op_src()", "Invalid op_id = %d\n", src_op);
        return out;
    }
}

static inline const WilsonMatrix &apply_op_snk(WilsonMatrix &out,
                                               const WilsonMatrix &in,
                                               Operator snk_op)
{
    WilsonMatrix tmp;
    switch(snk_op) {
    case GAMMA_0:
    case GAMMA_1:
    case GAMMA_2:
    case GAMMA_3:
        tmp.glA(in, snk_op - GAMMA_0);
        return out = -1*tmp;
    case GAMMA_5:
        return out = in;
    case ID:
        return out.glV(in, -5);
    case GAMMA_05:
    case GAMMA_15:
    case GAMMA_25:
    case GAMMA_35:
        tmp.glV(in, snk_op - GAMMA_05);
        return out = -1*tmp;
    case GAMMA_50:
    case GAMMA_51:
    case GAMMA_52:
    case GAMMA_53:
        return out.glV(in, snk_op - GAMMA_50);
    default:
        ERR.General("", "apply_op_snk()", "Invalid op_id = %d\n", snk_op);
        return out;
    }
}

// compute meson correlation functions point sink, source type is
// determined implicitly by the propagator.
void run_meson_pt(const AllProp &propA,
                  const AllProp &propB,
                  Operator snk_op,
                  Operator src_op,
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

    FILE *fp = Fopen(fn.c_str(), "w");

    for(unsigned k = 0; k < t_size; ++k) {
        if(propA.empty(k, ptype)) continue;
        if(propB.empty(k, ptype)) continue;

        std::vector<Rcomplex> meson(t_size_ap, Rcomplex(0, 0));

#pragma omp parallel
        {
            // threaded results
            std::vector<Rcomplex> tmp(t_size_ap, Rcomplex(0, 0));

#pragma omp for
            for(int i = 0; i < t_scale * lcl_vol; ++i) {
                int x[4];
                compute_coord_ap(x, lcl, i, t_size);
                int t_glb = x[3] + shift;
            
                WilsonMatrix w = propB(i, k, ptype);
                w.hconj();
                WilsonMatrix p[2];

                apply_op_snk(p[0], propA(i, k, ptype), snk_op);
                apply_op_src(p[1], w, src_op);

                tmp[t_glb] += Trace(p[0], p[1]);
            } // sites
#pragma omp critical
            for(int t = 0; t < t_size_ap; ++t) {
                meson[t] += tmp[(t+k)%t_size_ap];
            } // critical, for
        }//omp

        // FIXME
        assert(GJP.Snodes() == 1);
        QMP_sum_double_array((double *)meson.data(), 2 * t_size_ap);

        for(unsigned t = 0; t < t_size_ap; ++t) {
            Fprintf(fp, "%3u %3u %17.10e %17.10e\n", k, t,
                    meson[t].real(), meson[t].imag());
        }
    } //k

    Fclose(fp);
}

// compute meson correlation functions wall sink, source type is
// determined implicitly by the propagator.
void run_meson_wall(const AllProp &propA,
                    const AllProp &propB,
                    Operator snk_op,
                    Operator src_op,
                    const std::string &fn,
                    PROP_TYPE ptype)
{
    const int t_scale = ptype == PROP_PA ? 2 : 1;
    const int t_size = GJP.TnodeSites() * GJP.Tnodes();
    const int t_size_ap = t_scale * t_size;

    std::vector<std::vector<WilsonMatrix> > Awsnk, Bwsnk;
    run_wall_snk(&Awsnk, propA, ptype);
    run_wall_snk(&Bwsnk, propB, ptype);

    FILE *fp = Fopen(fn.c_str(), "w");

    for(unsigned src = 0; src < t_size; ++src) {
        if(Awsnk[src].empty()) continue;
        if(Bwsnk[src].empty()) continue;

        std::vector<Rcomplex> meson(t_size_ap, Rcomplex(0, 0));
        for(unsigned dt = 0; dt < t_size_ap; ++dt) {
            unsigned snk = (src + dt) % t_size_ap;

            WilsonMatrix w = Bwsnk[src][snk];
            w.hconj();

            WilsonMatrix p[2];
            apply_op_snk(p[0], Awsnk[src][snk], snk_op);
            apply_op_src(p[1], w, src_op);
            
            meson[dt] += Trace(p[0], p[1]);
        } // dt

        for(unsigned t = 0; t < t_size_ap; ++t) {
            Fprintf(fp, "%3u %3u %17.10e %17.10e\n", src, t,
                    meson[t].real(), meson[t].imag());
        }
    } // src
    
    Fclose(fp);
}

static inline const WilsonMatrix &apply_op_dis(WilsonMatrix &out,
                                               const WilsonMatrix &in,
                                               Operator op)
{
    WilsonMatrix tmp;
    switch(op) {
    case GAMMA_0:
    case GAMMA_1:
    case GAMMA_2:
    case GAMMA_3:
        return out.glV(in, op - GAMMA_0);
    case GAMMA_5:
        return out.glV(in, -5);
    case ID:
        return out = in;
    case GAMMA_05:
    case GAMMA_15:
    case GAMMA_25:
    case GAMMA_35:
        return out.glA(in, op - GAMMA_05);
    case GAMMA_50:
    case GAMMA_51:
    case GAMMA_52:
    case GAMMA_53:
        tmp.glA(in, op - GAMMA_50);
        return out = -1*tmp;
    default:
        ERR.General("", "apply_op_dis()", "Invalid op_id = %d\n", op);
        return out;
    }
}

// compute meson correlation functions wall sink, source type is
// determined implicitly by the propagator.
//
// Disconnected diagram
void run_meson_disc(const AllProp &propL,
                    const AllProp &propR,
                    Operator l_op,
                    Operator r_op,
                    const std::string &fn,
                    PROP_TYPE ptype)
{
    const int t_scale = ptype == PROP_PA ? 2 : 1;
    const int t_size = GJP.TnodeSites() * GJP.Tnodes();
    const int t_size_ap = t_scale * t_size;

    std::vector<std::vector<WilsonMatrix> > Lwsnk, Rwsnk;
    run_wall_snk(&Lwsnk, propL, ptype);
    run_wall_snk(&Rwsnk, propR, ptype);

    FILE *fp = Fopen(fn.c_str(), "w");

    for(unsigned src = 0; src < t_size_ap; ++src) {
        WilsonMatrix lp[2] = {0, 0};

        if(!Lwsnk[src].empty()) {
            apply_op_dis(lp[0], Lwsnk[src][src], l_op);
        }
        if(!Rwsnk[src].empty()) {
            apply_op_dis(lp[1], Rwsnk[src][src], r_op);
        }

        Rcomplex p0t = lp[0].Trace();
        Rcomplex p1t = lp[1].Trace();

        Fprintf(fp, "%3u %17.10e %17.10e %17.10e %17.10e\n", src,
                p0t.real(), p0t.imag(), p1t.real(), p1t.imag());
    } // src
    
    Fclose(fp);
}    

CPS_END_NAMESPACE
