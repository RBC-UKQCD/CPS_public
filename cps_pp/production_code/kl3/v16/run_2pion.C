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
#include "run_wsnk.h"
#include "run_2pion.h"

USING_NAMESPACE_CPS
using namespace std;

// 2 pion contractions, D and C diagram only (for I=2 pi pi
// scattering).
//
// Assuming cosine source with antiperiodic boundary condition.
//
// FIXME: there may be one or more factors of 2 in the sink side.
void run_2pionDC(const AllProp &uprop,
                 const AllProp &dprop,
                 const char fn_t[], int traj,
                 PROP_TYPE ptype, const int pd[3])
{
    const int t_scale = ptype == PROP_PA ? 2 : 1;
    const int t_size = GJP.TnodeSites() * GJP.Tnodes();
    const int t_size_ap = t_scale * t_size;

    vector<vector<WilsonMatrix> > us, dsp, dsm;

    run_wall_snk(&us, uprop, ptype);

    run_wall_snk(&dsp, dprop, ptype, pd);

    int pdm[3] = {-pd[0], -pd[1], -pd[2]};
    run_wall_snk(&dsm, dprop, ptype, pdm);

    char buf[1024];
    sprintf(buf, fn_t, ptype_str[ptype], traj);
    FILE *fp = Fopen(buf, "w");

    for(unsigned src = 0; src < t_size; ++src) {
        if( us[src].empty()) continue;
        if(dsp[src].empty()) continue;
        if(dsm[src].empty()) continue;

        vector<Rcomplex> ddiag(t_size_ap, Rcomplex(0, 0));
        vector<Rcomplex> cdiag(t_size_ap, Rcomplex(0, 0));

        for(unsigned dt = 0; dt < t_size_ap; ++dt) {
            unsigned snk = (src + dt) % t_size_ap;

            WilsonMatrix w[2];

            w[0] = dsp[src][snk];
            w[0].hconj();
            w[0] *= us[src][snk];

            w[1] = dsm[src][snk];
            w[1].hconj();
            w[1] *= us[src][snk];

            ddiag[dt] = w[0].Trace() * w[1].Trace();
            cdiag[dt] = Trace(w[0], w[1]);
        } // dt

        for(unsigned t = 0; t < t_size_ap; ++t) {
            Fprintf(fp, "%3u %3u %17.10e %17.10e %17.10e %17.10e\n",
                    src, t,
                    real(ddiag[t]), imag(ddiag[t]),
                    real(cdiag[t]), imag(cdiag[t]));
        }
    } // src
    
    Fclose(fp);
}
