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

USING_NAMESPACE_CPS
using namespace std;

// compute Bk correlation functions, source type is
// determined implicitly by the propagator.
void run_bk(const AllProp &lpropA,
            const AllProp &spropA,
            const AllProp &lpropB,
            const AllProp &spropB,
            const char fn_t[], int traj)
{
    const int t_size_glb = 2 * GJP.TnodeSites() * GJP.Tnodes();
    const int lcl[4] = {GJP.XnodeSites(), GJP.YnodeSites(),
                        GJP.ZnodeSites(), GJP.TnodeSites(),};
    const int lcl_vol = lcl[0] * lcl[1] * lcl[2] * lcl[3];
    const int shift = GJP.TnodeSites() * GJP.TnodeCoor();

    char buf[1024];
    sprintf(buf, fn_t, traj);
    FILE *fp = Fopen(buf, "w");

    for(int sep = 0; sep < t_size_glb; ++sep) {
        for(int ka = 0; ka < t_size_glb / 2; ++ka) {
            int kb = (ka + sep) % t_size_glb;
            if(lpropA.empty(ka)) continue;
            if(spropA.empty(ka)) continue;
            if(lpropB.empty(kb)) continue;
            if(spropB.empty(kb)) continue;

            // bk[0] = Tr(loop0) * Tr(loop1)
            // bk[1] = Tr(loop0 * loop1)
            vector<Rcomplex> bk0(t_size_glb, Rcomplex(0, 0));
            vector<Rcomplex> bk1(t_size_glb, Rcomplex(0, 0));

#pragma omp parallel
            {
                // threaded results
                vector<Rcomplex> t0(t_size_glb, Rcomplex(0, 0));
                vector<Rcomplex> t1(t_size_glb, Rcomplex(0, 0));

#pragma omp for
                for(int i = 0; i < 2 * lcl_vol; ++i) {
                    int x[4];
                    compute_coord_ap(x, lcl, i, t_size_glb / 2);
                    int t_glb = x[3] + shift;

                    // K & Kb
                    WilsonMatrix K[2]  = {lpropA[ka][i], spropA[ka][i]};
                    WilsonMatrix Kb[2] = {lpropB[kb][i], spropB[kb][i]};

                    K[1].hconj();
                    Kb[1].hconj();

                    WilsonMatrix lines[2] = {K[0] * K[1], Kb[0] * Kb[1]};

                    // bk[0] = Tr(loop0) * Tr(loop1)
                    // bk[1] = Tr(loop0 * loop1)
                    for(int mu = 0; mu < 4; ++mu) {
                        WilsonMatrix tmp[2];
                        // AA term
                        tmp[0].glA(lines[0], mu);
                        tmp[1].glA(lines[1], mu);

                        t0[t_glb] += tmp[0].Trace() * tmp[1].Trace();
                        t1[t_glb] += Trace(tmp[0], tmp[1]);

                        // VV term
                        tmp[0].glV(lines[0], mu);
                        tmp[1].glV(lines[1], mu);

                        t0[t_glb] += tmp[0].Trace() * tmp[1].Trace();
                        t1[t_glb] += Trace(tmp[0], tmp[1]);
                    }
                } // sites

#pragma omp critical
                {
                    for(int t = 0; t < t_size_glb; ++t) {
                        bk0[t] += t0[(t+ka)%t_size_glb];
                        bk1[t] += t1[(t+ka)%t_size_glb];
                    }
                }
            }//omp

            
            assert(GJP.Snodes() == 1);
            // FIXME
            QMP_sum_double_array((double *)bk0.data(), 2 * t_size_glb);
            QMP_sum_double_array((double *)bk1.data(), 2 * t_size_glb);

            for(int t = 0; t < t_size_glb; ++t) {
                Fprintf(fp, "%3d %3d %3d %17.10e %17.10e %17.10e %17.10e\n",
                        sep, ka, t, 
                        bk0[t].real(), bk0[t].imag(),
                        bk1[t].real(), bk1[t].imag());
            }
        }} // ka, sep

    Fclose(fp);
}
