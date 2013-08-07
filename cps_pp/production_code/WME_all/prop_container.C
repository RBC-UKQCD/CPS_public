// -*- mode:c++; c-basic-offset:4 -*-
#include <util/gjp.h>
#include <util/verbose.h>
#include <alg/qpropw.h>
#include <util/qio_general.h>
#include <omp.h>
#include <cassert>

#include <qio.h>
#include "my_util.h"
#include "prop_container.h"

USING_NAMESPACE_CPS

// store all propagators I have.
//
//! IMPORTANT: This function assumes gauge fixed wall source!!!
void AllProp::store_all(const std::string &fn_stem, double mass, int traj)const
{
    const char *fname = "store_all()";
    if(prec == SINGLE) {
        VRB.Result("AllProp", fname, "We don't store single precision propagators.\n");
        return;
    }

    // We need to reconstruct the source since we want to store the
    // propagators when the source was destroyed long time ago.

    for(int bc = 0; bc < 2; ++bc) {
        std::string bc_string = bc == 0 ? "P" : "A";
        for(size_t t = 0; t < t_size_glb; ++t) {
            std::string fn = fn_stem + bc_string +
                + "_m" + tostring(mass)
                + "_gfwall_" + tostring(t)
                + "." + tostring(traj);

            if(bc == 0) {
                if(pd[t].empty()) continue;
                store(fn, pd[t], t);
            } else {
                if(ad[t].empty()) continue;
                store(fn, ad[t], t);
            }
        }
    }
}

static void set_wall(cps::Rcomplex *src, int t_src)
{
    const int lcl[4] = {GJP.XnodeSites(), GJP.YnodeSites(),
                        GJP.ZnodeSites(), GJP.TnodeSites(),};
    const int lcl_vol = lcl[0] * lcl[1] * lcl[2] * lcl[3];
    const int shift = GJP.TnodeSites() * GJP.TnodeCoor();

    for(int i = 0; i < lcl_vol; ++i) {
        int x[4];
        compute_coord(x, lcl, i);
        int t_glb = x[3] + shift;

        src[i] = t_src == t_glb ? 1.0 : 0.0;
    }
}

void AllProp::store(const std::string &fn,
                    const std::vector<cps::WilsonMatrix> &prop,
                    int t)const
{
    const size_t sites = GJP.VolNodeSites();

    std::string fn_node = fn + "." + tostring(UniqueID());

    std::FILE *fp = fopen(fn_node.c_str(), "w");

    fwrite(prop.data(), sizeof(WilsonMatrix) * sites, 1, fp);

    fclose(fp);
}
