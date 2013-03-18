// -*- mode:c++; c-basic-offset:4 -*-
#include <util/rcomplex.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/lattice/fbfm.h>

#include <omp.h>
#include <cmath>

#include "my_util.h"

USING_NAMESPACE_CPS
using namespace std;

// ----------------------------------------------------------------
// twisted_bc: set twisted boundary condition.
//
// The gauge field must be in CANONICAL order.
//
// add: if true then multiply exp(i*angle), otherwise multiply
// exp(-i*angle). This parameter can be used to add/remove the bc.
// ----------------------------------------------------------------
void twisted_bc(Lattice &lat, const double mom[4], bool add)
{
    Matrix *gauge = (Matrix *)(lat.GaugeField());

    const double PI = 3.1415926535897932384626433832795028842;
    double sign = add ? 1 : -1;

    for(int mu = 0; mu < 4; ++mu) {
        if(mom[mu] == 0) continue;

        double t = 2.0 * PI / GJP.Sites(mu) * mom[mu];
        const Rcomplex cval(cos(t), sign * sin(t));

        int low[4] = { 0, 0, 0, 0 };
        int high[4] = { GJP.XnodeSites(), GJP.YnodeSites(),
                        GJP.ZnodeSites(), GJP.TnodeSites() };

        int hl[4] = { high[0] - low[0], high[1] - low[1],
                      high[2] - low[2], high[3] - low[3] };

        const int hl_sites = hl[0] * hl[1] * hl[2] * hl[3];

#pragma omp parallel for
        for(int i = 0; i < hl_sites; ++i) {
            int x[4];
            compute_coord(x, hl, low, i);

            int off = mu + 4 * compute_id(x, high);
            gauge[off] *= cval;
        }
    }

    Fbfm *fbfm = dynamic_cast<Fbfm *>(&lat);
    if(fbfm != NULL) {
        fbfm->ImportGauge();
    }
}
