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
#include "run_pion.h"
#include "run_wsnk.h"

USING_NAMESPACE_CPS
using namespace std;

// compute wall sink propagator
void run_wall_snk(vector<vector<WilsonMatrix> > *wsnk,
                  const AllProp &prop)
{
    const int t_size_glb = 2 * GJP.TnodeSites() * GJP.Tnodes();
    const int lcl[4] = {GJP.XnodeSites(), GJP.YnodeSites(),
                        GJP.ZnodeSites(), GJP.TnodeSites(),};
    const int lcl_vol = lcl[0] * lcl[1] * lcl[2] * lcl[3];
    const int shift = GJP.TnodeSites() * GJP.TnodeCoor();

    wsnk->assign(t_size_glb, vector<WilsonMatrix>());
    for(unsigned tsrc = 0; tsrc < t_size_glb; ++tsrc) {
        if(prop.empty(tsrc)) continue;
        (*wsnk)[tsrc].assign(t_size_glb, WilsonMatrix(0.));

#pragma omp parallel
        {
            vector<WilsonMatrix> w(t_size_glb, 0.);

#pragma omp for
            for(int i = 0; i < 2 * lcl_vol; ++i) {
                int x[4];
                compute_coord_ap(x, lcl, i, t_size_glb / 2);

                w[x[3] + shift] += prop[tsrc][i];
            } // sites

#pragma omp critical
            for(int t = 0; t < t_size_glb; ++t) {
                (*wsnk)[tsrc][t] += w[t];
            } // critial, for
        }//omp

        assert(GJP.Snodes() == 1);
        QMP_sum_double_array((double *)(*wsnk)[tsrc].data(),
                             t_size_glb * sizeof(WilsonMatrix) / sizeof(double));
    }//tsrc
}
