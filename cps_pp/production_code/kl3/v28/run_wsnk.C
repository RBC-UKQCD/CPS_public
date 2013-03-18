// -*- mode:c++; c-basic-offset:4 -*-
#include <vector>
#include <alg/qpropw.h>
#include "prop_container.h"

CPS_START_NAMESPACE

// Computing the phase factor, assuming antiperiodic BC
static Rcomplex cal_phaseA(const int p[3], const int x[3])
{
    double phase = 0;
    for(int i = 0; i < 3; ++i) {
        phase += p[i] * x[i] * 3.1415926535897932384626433832795028842L / GJP.Sites(i);
    }
    return Rcomplex(cos(phase), sin(phase));
}

// compute wall sink propagator, with momentum
void run_wall_snk(std::vector<std::vector<WilsonMatrix> > *wsnk,
                  const AllProp &prop,
                  PROP_TYPE ptype, const int *p = NULL)
{
    const int t_scale = ptype == PROP_PA ? 2 : 1;
    const int t_size = GJP.TnodeSites() * GJP.Tnodes();
    const int t_size_ap = t_scale * t_size;
    const int lcl[4] = {GJP.XnodeSites(), GJP.YnodeSites(),
                        GJP.ZnodeSites(), GJP.TnodeSites(),};
    const int lcl_vol = lcl[0] * lcl[1] * lcl[2] * lcl[3];
    const int shift[4] = {
        GJP.XnodeSites() * GJP.XnodeCoor(),
        GJP.YnodeSites() * GJP.YnodeCoor(),
        GJP.ZnodeSites() * GJP.ZnodeCoor(),
        GJP.TnodeSites() * GJP.TnodeCoor()
    };

    wsnk->assign(t_size_ap, std::vector<WilsonMatrix>());
    for(unsigned tsrc = 0; tsrc < t_size_ap; ++tsrc) {
        if(prop.empty(tsrc, ptype)) continue;
        (*wsnk)[tsrc].assign(t_size_ap, WilsonMatrix(0.));

#pragma omp parallel
        {
            std::vector<WilsonMatrix> w(t_size_ap, 0.);

#pragma omp for
            for(int i = 0; i < t_scale * lcl_vol; ++i) {
                int x[4];
                compute_coord_ap(x, lcl, i, t_size);
                for(unsigned k = 0; k < 4; ++k) {
                    x[k] += shift[k];
                }
                int t_glb = x[3];

                if(p != NULL) {
                    w[t_glb] += cal_phaseA(p, x) * prop(i, tsrc, ptype);
                } else {
                    w[t_glb] += prop(i, tsrc, ptype);
                }
            } // sites

#pragma omp critical
            for(int t = 0; t < t_size_ap; ++t) {
                (*wsnk)[tsrc][t] += w[t];
            } // critial, for
        }//omp

        assert(GJP.Snodes() == 1);
        QMP_sum_double_array((double *)(*wsnk)[tsrc].data(),
                             t_size_ap * sizeof(WilsonMatrix) / sizeof(double));
    }//tsrc
}

CPS_END_NAMESPACE
