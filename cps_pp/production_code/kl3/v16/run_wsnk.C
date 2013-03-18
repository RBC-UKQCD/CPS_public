// -*- mode:c++; c-basic-offset:4 -*-
#include <cstring>

#include <qmp.h>
#include <util/gjp.h>
#include <util/rcomplex.h>
#include <util/qcdio.h>
#include <alg/qpropw.h>
#include <vector>
#include <cassert>
#include <cmath>

#include <omp.h>

#include "my_util.h"
#include "run_wsnk.h"

USING_NAMESPACE_CPS
using namespace std;

// Computing the phase factor, assuming antiperiodic BC
static Rcomplex cal_phaseA(const int p[3], const int x[3])
{
    double phase = 0;
    for(int i = 0; i < 3; ++i) {
        phase += p[i] * x[i] * 3.1415926535897932384626433832795028842L / GJP.Sites(i);
    }
    return Rcomplex(cos(phase), sin(phase));
}

// compute wall sink propagator
void run_wall_snk(vector<vector<WilsonMatrix> > *wsnk,
                  const AllProp &prop,
                  PROP_TYPE ptype, const int *p)
{
    const int t_scale = ptype == PROP_PA ? 2 : 1;
    const int t_size = GJP.TnodeSites() * GJP.Tnodes();
    const int t_size_ap = t_scale * t_size;
    const int lcl[4] = {GJP.XnodeSites(), GJP.YnodeSites(),
                        GJP.ZnodeSites(), GJP.TnodeSites(),};
    const int lcl_vol = lcl[0] * lcl[1] * lcl[2] * lcl[3];
    const int shift = GJP.TnodeSites() * GJP.TnodeCoor();

    wsnk->assign(t_size_ap, vector<WilsonMatrix>());
    for(unsigned tsrc = 0; tsrc < t_size_ap; ++tsrc) {
        if(prop.empty(tsrc, ptype)) continue;
        (*wsnk)[tsrc].assign(t_size_ap, WilsonMatrix(0.));

#pragma omp parallel
        {
            vector<WilsonMatrix> w(t_size_ap, 0.);

#pragma omp for
            for(int i = 0; i < t_scale * lcl_vol; ++i) {
                int x[4];
                compute_coord_ap(x, lcl, i, t_size);
                int t_glb = x[3] + shift;

                if(p != NULL) {
                    w[t_glb] += cal_phaseA(p, x) * prop(tsrc, ptype)(i, ptype);
                } else {
                    w[t_glb] += prop(tsrc, ptype)(i, ptype);
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
