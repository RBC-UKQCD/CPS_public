// -*- mode:c++; c-basic-offset:4 -*-
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <util/lattice.h>
#include <util/lattice/fbfm.h>

#include <alg/array_arg.h>
#include <alg/alg_int.h>
#include <alg/eigcg_arg.h>
#include <alg/qpropw.h>

#include <util/gjp.h>
#include <util/time_cps.h>
#include <util/verbose.h>
#include <util/error.h>

#undef USE_HDCG

#ifdef USE_HDCG
#include <util/lattice/hdcg_controller.h>
#endif

#include <string>
#include <vector>
#include <cassert>

#include "prop_container.h"
#include "eigcg.h"
#include "my_util.h"
#include "run_mres.h"

static const int SAVE_PROP = 0;

static const char *cname = "";

USING_NAMESPACE_CPS
using namespace std;

static void run_mres_za(const QPropW &qp, const QPropWArg &qp_arg,
                        const string &rdir, int traj)
{
    string mres_fn = rdir + "/mres_"
        + tostring(qp_arg.cg.mass) + "."
        + tostring(traj);

    string za_fn = rdir + "/za_"
        + tostring(qp_arg.cg.mass) + "."
        + tostring(traj);

    run_mres(qp, qp_arg.t, mres_fn.c_str());
    run_za(qp, qp_arg.cg.mass, qp_arg.t, za_fn.c_str());
}

// Note: How many times we solve the volume source depends on how many
// low modes we want to solve. Lattice properties also apply.
//
// For 300 low modes, 1 propagator using mixed solver will be good
// (depends on EigCGArg).
//
// On 48^3 2 solves are needed for 600 low modes.
static void collect_lowmodes(Lattice &lat,
                             QPropWArg &qp_arg,
                             CommonArg &com_prop)
{
    double stop_rsd = qp_arg.cg.stop_rsd;
    double true_rsd = qp_arg.cg.true_rsd;

    qp_arg.cg.stop_rsd = 1e-10;
    qp_arg.cg.true_rsd = 1e-10;

    QPropW4DBoxArg vol_arg;
    for(int mu = 0; mu < 4; ++mu) {
        vol_arg.box_start[mu] = 0;
        vol_arg.box_size[mu] = GJP.Sites(mu);
        vol_arg.mom[mu] = 0;
    }

    // 4 solves for 1500 low modes.
    for(unsigned i = 0; i < 4; ++i) {
        QPropW4DBoxSrc qp_box(lat, &qp_arg, &vol_arg, &com_prop);
    }

    qp_arg.cg.stop_rsd = stop_rsd;
    qp_arg.cg.true_rsd = true_rsd;
}

void run_wall_prop(const char *pname,
		   AllProp *prop_e,
                   AllProp *prop,
                   IntArray &eloc, 
                   Lattice &lat,
                   QPropWArg &qp_arg,
                   EigCGArg *eigcg_arg,
                   int traj,
                   bool do_mres)
{
    const char *fname = "run_wall_prop()";

    // Check boundary condition. We need this to ensure that we are
    // doing P + A and P - A, not A + P and A - P (I think it's OK to
    // skip this check, though).
    if(GJP.Tbc() == BND_CND_APRD) {
        ERR.General(cname, fname, "Boundary condition does not match!\n");
    }

    char buf[256];
    CommonArg com_prop;
    sprintf(buf, "../results/%s.%d", qp_arg.ensemble_label, traj);
    com_prop.set_filename(buf);

    // A only
    for(int bc = 1; bc < 2; ++bc) {
        GJP.Tbc(bc == 0 ? BND_CND_PRD : BND_CND_APRD);
        lat.BondCond();
#ifdef USE_HDCG
	HDCGController<Float> *control = HDCGController<Float>::getInstance();
        if (control) control->freeHDCG();
#endif

        EigCG *eig_cg = NULL;
        if(eigcg_arg) {
            eig_cg = new EigCG(eigcg_arg, Fbfm::use_mixed_solver);
            collect_lowmodes(lat, qp_arg, com_prop);

            const string fn = string("../results") + (bc == 0 ? "EP" : "EA")
                + "/eigH_" + (do_mres ? "wall_" : "twist_")
                + tostring(qp_arg.cg.mass) + "."
                + tostring(traj);

            eig_cg->printH(fn);
        }

        // exact propagators
        if(prop_e != NULL) {
            double stop_rsd = qp_arg.cg.stop_rsd;
            double true_rsd = qp_arg.cg.true_rsd;
            
            qp_arg.cg.stop_rsd = 1e-8;
            qp_arg.cg.true_rsd = 1e-8;
            for(unsigned i = 0; i < eloc.v.v_len; ++i) {
                qp_arg.t = eloc.v.v_val[i];
                VRB.Result(cname, fname, "Solving exact propagator at %d\n", qp_arg.t);
		string t_name = string(pname)+"_e_"+tostring(qp_arg.t)+"_"+tostring(bc)+"."+tostring(traj);
		VRB.Result(cname,fname,"qp_wall.SaveQProp(%s,0)\n",t_name.c_str());
		qp_arg.file = &t_name[0];
		qp_arg.save_prop = SAVE_PROP;

                QPropWWallSrc qp_wall(lat, &qp_arg, &com_prop);
                if(do_mres) {
                    run_mres_za(qp_wall, qp_arg,
                                string("../results") + (bc == 0 ? "EP" : "EA"),
                                traj);
                }
                prop_e->add(qp_wall, qp_arg.t, bc == 0);
            }
            qp_arg.cg.stop_rsd = stop_rsd;
            qp_arg.cg.true_rsd = true_rsd;
	    qp_arg.save_prop = 0;
        }

        // inexact propagators
        for(int t = 0; t < GJP.Sites(3); ++t) {
            qp_arg.t = t;
	    string t_name = string(pname)+"_"+tostring(qp_arg.t)+"_"+tostring(bc)+"."+tostring(traj);
	    VRB.Result("",fname,"qp_wall.SaveQProp(%s,0)\n",t_name.c_str());
		qp_arg.file = &t_name[0];
		qp_arg.save_prop = SAVE_PROP;
            QPropWWallSrc qp_wall(lat, &qp_arg, &com_prop);
            if(do_mres) {
                run_mres_za(qp_wall, qp_arg,
                            string("../results") + (bc == 0 ? "P" : "A"),
                            traj);
            }
            prop->add(qp_wall, qp_arg.t, bc == 0);
	    qp_arg.save_prop = 0;
        }

        delete eig_cg;
        lat.BondCond();
    }

    // Note: If I call lat.BondCond() even times, then there is no
    // overall effect.
    GJP.Tbc(BND_CND_PRD);
}

void run_mom_prop(const char *pname,
		  AllProp *prop_e,
                  AllProp *prop,
                  IntArray &eloc,
                  Lattice &lat,
                  QPropWArg &qp_arg,
                  EigCGArg *eigcg_arg,
                  int traj,
                  const int mom[3])
{
    const char *fname = "run_mom_prop()";

    // Ensure that all 4 directions have periodic boundary condition.
    // FIXME: This check is not perfect as we have no way detecting
    // how the actual gauge field data were manipulated.
    for(int mu = 0; mu < 4; ++mu) {
        if(GJP.Bc(mu) == BND_CND_APRD) {
            ERR.General(cname, fname, "Boundary condition does not match!\n");
        }
        if(mu < 3 && mom[mu]) {
            GJP.Bc(mu, BND_CND_APRD);
        }
    }


    char buf[256];
    CommonArg com_prop;
    sprintf(buf, "../results/%s.%d", qp_arg.ensemble_label, traj);
    com_prop.set_filename(buf);

    // A only
    for(int bc = 1; bc < 2; ++bc) {
        GJP.Tbc(bc == 0 ? BND_CND_PRD : BND_CND_APRD);
        lat.BondCond();
#ifdef USE_HDCG
        HDCGController<Float> *control = HDCGController<Float>::getInstance();
        if (control) control->freeHDCG();
#endif

        EigCG *eig_cg = NULL;
        if(eigcg_arg) {
            eig_cg = new EigCG(eigcg_arg, Fbfm::use_mixed_solver);
            collect_lowmodes(lat, qp_arg, com_prop);

            const string fn = string("../results") + (bc == 0 ? "EP" : "EA")
                + "/eigH_mom_" + tostring(qp_arg.cg.mass) + "."
                + tostring(traj);

            eig_cg->printH(fn);
        }

        // exact propagators
        if(prop_e != NULL) {
            double stop_rsd = qp_arg.cg.stop_rsd;
            double true_rsd = qp_arg.cg.true_rsd;

            qp_arg.cg.stop_rsd = 1e-8;
            qp_arg.cg.true_rsd = 1e-8;
            for(unsigned i = 0; i < eloc.v.v_len; ++i) {
                qp_arg.t = eloc.v.v_val[i];
                VRB.Result(cname, fname, "Solving exact propagator at %d\n", qp_arg.t);

	    string t_name = string(pname)+"_"+tostring(qp_arg.t)+"_"+tostring(bc)
		+tostring(mom[0])
		+tostring(mom[1])
		+tostring(mom[2])
		+"."+tostring(traj);
		VRB.Result("",fname,"qp_mom.SaveQProp(%s,0)\n",t_name.c_str());
		qp_arg.file = &t_name[0];
		qp_arg.save_prop = SAVE_PROP;
                QPropWMomCosTwistSrc qp_mom(lat, &qp_arg, mom, &com_prop);
                prop_e->add(qp_mom, qp_arg.t, bc == 0);
            }
            qp_arg.cg.stop_rsd = stop_rsd;
            qp_arg.cg.true_rsd = true_rsd;
	    qp_arg.save_prop = 0;
        }

        // inexact propagators
        for(int t = 0; t < GJP.Sites(3); ++t) {
            qp_arg.t = t;
	    string t_name = string(pname)+"_"+tostring(qp_arg.t)+"_"+tostring(bc)
		+tostring(mom[0])
		+tostring(mom[1])
		+tostring(mom[2])
		+"."+tostring(traj);
	    VRB.Result("",fname,"qp_mom.SaveQProp(%s,0)\n",t_name.c_str());
		qp_arg.file = &t_name[0];
		qp_arg.save_prop = SAVE_PROP;
            QPropWMomCosTwistSrc qp_mom(lat, &qp_arg, mom, &com_prop);
            prop->add(qp_mom, qp_arg.t, bc == 0);
        }
	qp_arg.save_prop = 0;

        delete eig_cg;
        lat.BondCond();
    }

    // Note: If I call lat.BondCond() even times, then there is no
    // overall effect.
    for(int mu = 0; mu < 4; ++mu) {
        GJP.Bc(mu, BND_CND_PRD);
    }
}
