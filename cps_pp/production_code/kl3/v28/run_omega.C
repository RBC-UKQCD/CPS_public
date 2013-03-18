// -*- mode:c++; c-basic-offset:4 -*-
#include <vector>
#include <string>
#include <alg/qpropw.h>
#include <util/qcdio.h>
#include "prop_container.h"
#include "my_util.h"
#include "run_wsnk.h"

CPS_START_NAMESPACE

// Computes the following,
// out = C gamma_op in gamma_t gamma_op^\dag C^\dag gamma_t
static inline WilsonMatrix
apply_op(const WilsonMatrix &in, Operator op)
{
    WilsonMatrix ret;

    ret.grV(in, 3);
    switch(op) {
    case GAMMA_0:
    case GAMMA_1:
    case GAMMA_2:
    case GAMMA_3:
        ret.gl(op - GAMMA_0);
        ret.gr(op - GAMMA_0);
        break;
    case GAMMA_5:
        ret.gl(-5);
        ret.gr(-5);
        break;
    case ID:
        break;
    case GAMMA_05:
    case GAMMA_15:
    case GAMMA_25:
    case GAMMA_35:
    case GAMMA_50:
    case GAMMA_51:
    case GAMMA_52:
    case GAMMA_53:
    default:
        ERR.General("", "apply_op()", "Not implemented yet for op_id = %d\n", op);
        break;
    }

    ret.gl(1);
    ret.gl(3);
    ret.gr(1);
    return ret;
}

enum ContractionType {Type0, Type1};

// Type0: e_{abc}e_{def} qi_{ad} * Tr^spin [qj_{be} qk^T_{cf}]
// Type1: e_{abc}e_{def} qi_{ad} * qj_{be}^T qk_{cf}
static SpinMatrix
color_contract(const WilsonMatrix q[3], ContractionType ct)
{
    SpinMatrix a(0.0);
    static const int eidx[6][3] = {
        {0, 1, 2},  // +
        {1, 2, 0},  // +
        {2, 0, 1},  // +
        {0, 2, 1},  // -
        {2, 1, 0},  // -
        {1, 0, 2},  // -
    };

    for(int i = 0; i < 6; ++i) {
        for(int j = 0; j < 6; ++j) {
            bool positive = i < 3 == j < 3;
            SpinMatrix x[3] = {
                q[0].ColorComponent(eidx[i][0], eidx[j][0]),
                q[1].ColorComponent(eidx[i][1], eidx[j][1]),
                q[2].ColorComponent(eidx[i][2], eidx[j][2]),
            };

            SpinMatrix v = ct == Type0
                ? x[0] * TraceTranspose(x[1], x[2])
                : x[0] * x[1].transpose() * x[2];
            a = positive ? a + v : a - v;
        }
    }
    
    return a;
}

// compute omega correlation functions point sink, source type is
// determined implicitly by the propagator.
void run_omega_pt(const AllProp &prop,
                  Operator op,
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
        if(prop.empty(k, ptype)) continue;

        std::vector<SpinMatrix> baryon[2] = {
            std::vector<SpinMatrix>(t_size_ap, SpinMatrix(0.0)),
            std::vector<SpinMatrix>(t_size_ap, SpinMatrix(0.0)),
        };

#pragma omp parallel
        {
            // threaded results
            std::vector<SpinMatrix> tmp[2] = {
                std::vector<SpinMatrix>(t_size_ap, SpinMatrix(0.0)),
                std::vector<SpinMatrix>(t_size_ap, SpinMatrix(0.0)),
            };

#pragma omp for
            for(int i = 0; i < t_scale * lcl_vol; ++i) {
                int x[4];
                compute_coord_ap(x, lcl, i, t_size);
                int t_glb = x[3] + shift;
            
                const WilsonMatrix p[3] = {
                    prop(i, k, ptype),
                    apply_op(prop(i, k, ptype), op),
                    prop(i, k, ptype),
                };

                tmp[0][t_glb] += color_contract(p, Type0);
                tmp[1][t_glb] += color_contract(p, Type1);
            } // sites
#pragma omp critical
            for(int t = 0; t < t_size_ap; ++t) {
                if(ptype == PROP_A && t + k >= t_size_ap) {
                    baryon[0][t] -= tmp[0][(t+k)%t_size_ap];
                    baryon[1][t] -= tmp[1][(t+k)%t_size_ap];
                } else {
                    baryon[0][t] += tmp[0][(t+k)%t_size_ap];
                    baryon[1][t] += tmp[1][(t+k)%t_size_ap];
                }
            } // critical, for
        }//omp

        // FIXME
        assert(GJP.Snodes() == 1);
        assert(sizeof(Float) == sizeof(double));

        size_t length = t_size_ap * sizeof(SpinMatrix) / sizeof(double);
        QMP_sum_double_array((double *)baryon[0].data(), length);
        QMP_sum_double_array((double *)baryon[1].data(), length);

        for(unsigned t = 0; t < t_size_ap; ++t) {
            Fprintf(fp, "%3u %3u", k, t);
            for(unsigned type = 0; type < 2; ++type) {
                for(unsigned i = 0; i < 16; ++i) {
                    unsigned r = i / 4;
                    unsigned c = i % 4;
                    
                    const Rcomplex &val = baryon[type][t](r, c);
                    Fprintf(fp, " %17.10e %17.10e",
                            val.real(), val.imag());
                }//element
            }//type
            Fprintf(fp, "\n");
        }//time slices
    } //k

    Fclose(fp);
}

CPS_END_NAMESPACE
