// -*- mode:c++; c-basic-offset:4 -*-
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <util/lattice.h>
#include <util/lattice/fbfm.h>

#include <alg/common_arg.h>
#include <alg/alg_int.h>
#include <alg/eigcg_arg.h>
#include <alg/alg_fix_gauge.h>
#include <alg/qpropw.h>

#include <util/gjp.h>
#include <util/time_cps.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/qcdio.h>
#include <util/rcomplex.h>
#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>

#include <omp.h>

#include <string>
#include <vector>
#include <cassert>

#include "eigcg.h"
#include "run_kl3.h"
#include "run_bk.h"
#include "run_mres.h"
#include "run_meson.h"
#include "twisted_bc.h"
#include "prop_container.h"

const char *cname = "";

USING_NAMESPACE_CPS
using namespace std;

DoArg do_arg;
MeasArg meas_arg;
QPropWArg lqpropw_arg;
QPropWArg sqpropw_arg;
QPropW4DBoxArg wall_arg;
FixGaugeArg fix_gauge_arg;
EigCGArg l_eigcg_arg;

// The following are VML files feeding in light quark twist and
// strange quark twist. They have nothing to do with box sources.
QPropW4DBoxArg l_twist_arg;
QPropW4DBoxArg s_twist_arg;

#define decode_vml(arg_name)  do{                                       \
        if ( ! arg_name.Decode(#arg_name".vml", #arg_name) )            \
            ERR.General(cname, fname, "Bad " #arg_name ".vml.\n");      \
    } while(0)  

void decode_vml_all(void)
{
    const char *fname = "decode_vml_all()";

    decode_vml(do_arg);
    decode_vml(meas_arg);
    decode_vml(lqpropw_arg);
    decode_vml(wall_arg);
    decode_vml(sqpropw_arg);
    decode_vml(fix_gauge_arg);
    decode_vml(l_eigcg_arg);

    decode_vml(l_twist_arg);
    decode_vml(s_twist_arg);
}

void load_checkpoint(int traj);
void init_bfm(int *argc, char **argv[]);
void setup(int argc, char *argv[]);

class FixGauge
{
public:
    FixGauge(Lattice &lat, FixGaugeArg &fix_gauge_arg, int traj) {
        char buf[256];
        sprintf(buf, "../results/fg-bc.%d", traj);
        Fclose(Fopen(buf, "w"));
        com_fg.set_filename(buf);
        
        Float dtime = -dclock();
        fg = new AlgFixGauge(lat, &com_fg, &fix_gauge_arg);
        fg->run();
        dtime += dclock();
        VRB.Result("AlgFixGauge", "run()", "takes %17.10e seconds.\n", dtime);
    }

    ~FixGauge() {
        fg->free();
        delete fg;
    }
private:
    CommonArg com_fg;
    AlgFixGauge *fg;
};

void run_mres_za(const QPropW &qp, const QPropWArg &qp_arg, unsigned tsrc, char bc, int traj)
{
    char buf[1024];
    sprintf(buf, "../results/mres_%g_%c.%d", qp_arg.cg.mass, bc, traj);
    run_mres(qp, tsrc, buf);

    sprintf(buf, "../results/za_%g_%c.%d", qp_arg.cg.mass, bc, traj);
    run_za(qp, qp_arg.cg.mass, tsrc, buf);
}

void run_prop(AllProp *prop, Lattice &lat,
              QPropWArg &qp_arg, QPropW4DBoxArg &box_arg,
              EigCGArg *eigcg_arg, int traj, bool do_mres)
{
    const char *fname = "run_light_prop()";
    const int t_size_glb = GJP.TnodeSites() * GJP.Tnodes();

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

    // P + A
    for(int bc = 0; bc < 2; ++bc) {
        GJP.Tbc(bc == 0 ? BND_CND_PRD : BND_CND_APRD);
        lat.BondCond();

        EigCG *eig_cg = NULL;
        if(eigcg_arg) {
            eig_cg = new EigCG(eigcg_arg, Fbfm::use_mixed_solver);
        }

        for(int t = 0; t < t_size_glb; ++t) {
            box_arg.box_start[3] = t;
            QPropW4DBoxSrc qp_box(lat, &qp_arg, &box_arg, &com_prop);
            if(do_mres) {
                run_mres_za(qp_box, qp_arg, t, bc == 0 ? 'P' : 'A', traj);
            }

            prop->add(qp_box, t, bc == 0);
        }

        delete eig_cg;
        lat.BondCond();
    }

    // Note: If I call lat.BondCond() even times, then there is no
    // overall effect.
    GJP.Tbc(BND_CND_PRD);
}

// stw: twisting angle of the strange quark (connecting the operator
// and the kaon).
//
// ltw: twisting angle of the light quark (connecting the operator and
// the pion).
void run_all(Lattice &lat,
             const double stw[4], // strange quark twists
             const double ltw[4], // light quark twists
             int traj)
{
    const char *fname = "run_all()";

    Float dtime0 = dclock();
    //////////////////////////////////////////////////////////////////////
    // 1. props: strange, strange twisted, light, light twisted
    AllProp sprop, stwst, lprop, ltwst;

    {
        FixGauge fg(lat, fix_gauge_arg, traj);

        // light quark
        run_prop(&lprop, lat, lqpropw_arg, wall_arg, &l_eigcg_arg, traj, true );
        twisted_bc(lat, ltw, true);
        run_prop(&ltwst, lat, lqpropw_arg, wall_arg, &l_eigcg_arg, traj, false);
        twisted_bc(lat, ltw, false);

        // strange quark
        run_prop(&sprop, lat, sqpropw_arg, wall_arg, NULL,         traj, true );
        twisted_bc(lat, stw, true);
        run_prop(&stwst, lat, sqpropw_arg, wall_arg, NULL,         traj, false);
        twisted_bc(lat, stw, false);
    }

    Float dtime1 = dclock();
    //////////////////////////////////////////////////////////////////////
    // 2. meson contractions
    run_meson_pt(lprop, lprop, GAMMA_5, GAMMA_5, "../results/pion-00WP.%d", traj);
    run_meson_pt(lprop, ltwst, GAMMA_5, GAMMA_5, "../results/pion-01WP.%d", traj);
    run_meson_pt(sprop, lprop, GAMMA_5, GAMMA_5, "../results/kaon-00WP.%d", traj);
    run_meson_pt(stwst, lprop, GAMMA_5, GAMMA_5, "../results/kaon-10WP.%d", traj);
    run_meson_pt(sprop, ltwst, GAMMA_5, GAMMA_5, "../results/kaon-01WP.%d", traj);
    
    run_meson_wall(lprop, lprop, GAMMA_5, GAMMA_5, "../results/pion-00WW.%d", traj);
    run_meson_wall(lprop, ltwst, GAMMA_5, GAMMA_5, "../results/pion-01WW.%d", traj);
    run_meson_wall(sprop, lprop, GAMMA_5, GAMMA_5, "../results/kaon-00WW.%d", traj);
    run_meson_wall(stwst, lprop, GAMMA_5, GAMMA_5, "../results/kaon-10WW.%d", traj);
    run_meson_wall(sprop, ltwst, GAMMA_5, GAMMA_5, "../results/kaon-01WW.%d", traj);

    // eta eta' contractions
    //
    // We share the light-light propagator with pion contractions.
    run_meson_wall(sprop, sprop, GAMMA_5, GAMMA_5, "../results/ss-00WW.%d", traj);
    run_meson_dis_wall(lprop, sprop, GAMMA_5, GAMMA_5, "../results/ls-dis-00WW.%d", traj);

    // f_K and f_pi measurements
    run_meson_pt(lprop, lprop, GAMMA_35, GAMMA_5, "../results/fp-00WP.%d", traj);
    run_meson_pt(sprop, lprop, GAMMA_35, GAMMA_5, "../results/fk-00WP.%d", traj);

    // rho meson
    run_meson_pt(lprop, lprop, GAMMA_0, GAMMA_0, "../results/rho-x-00WP.%d", traj);
    run_meson_pt(lprop, lprop, GAMMA_1, GAMMA_1, "../results/rho-y-00WP.%d", traj);
    run_meson_pt(lprop, lprop, GAMMA_2, GAMMA_2, "../results/rho-z-00WP.%d", traj);
    run_meson_wall(lprop, lprop, GAMMA_0, GAMMA_0, "../results/rho-x-00WW.%d", traj);
    run_meson_wall(lprop, lprop, GAMMA_1, GAMMA_1, "../results/rho-y-00WW.%d", traj);
    run_meson_wall(lprop, lprop, GAMMA_2, GAMMA_2, "../results/rho-z-00WW.%d", traj);

    Float dtime2 = dclock();
    //////////////////////////////////////////////////////////////////////
    // 3. Kl3
    run_kl3(sprop, lprop, lprop, "../results/kl3-00-g%s.%d", traj);
    run_kl3(sprop, lprop, ltwst, "../results/kl3-01-g%s.%d", traj);
    run_kl3(stwst, lprop, lprop, "../results/kl3-10-g%s.%d", traj);

    // The following contractions are used for Z_V.
    // Useful for the double ratio method or the UKQCD method.
    run_kl3(lprop, lprop, lprop, "../results/zpa-00-g%s.%d", traj);
    run_kl3(sprop, lprop, sprop, "../results/zka-00-g%s.%d", traj);
    run_kl3(lprop, sprop, lprop, "../results/zkb-00-g%s.%d", traj);

    // Since line 1 and 3 both carry momentum, are their directions
    // consistent?
    // I think they are.
    // run_kl3(ltwst, lprop, ltwst, "../results/zpa-11-g%s.%d", traj);
    // run_kl3(lprop, ltwst, lprop, "../results/zpb-11-g%s.%d", traj);
    // run_kl3(stwst, lprop, stwst, "../results/zka-11-g%s.%d", traj);
    // run_kl3(lprop, stwst, lprop, "../results/zkb-11-g%s.%d", traj);

    Float dtime3 = dclock();
    //////////////////////////////////////////////////////////////////////
    // 4. Bk
    run_bk(lprop, sprop, lprop, sprop, "../results/bk.%d", traj);

    Float dtime4 = dclock();

    VRB.Result(cname, fname, "prop  = %17.10e seconds\n", dtime1 - dtime0);
    VRB.Result(cname, fname, "meson = %17.10e seconds\n", dtime2 - dtime1);
    VRB.Result(cname, fname, "Kl3   = %17.10e seconds\n", dtime3 - dtime2);
    VRB.Result(cname, fname, "Bk    = %17.10e seconds\n", dtime4 - dtime3);
}

int main(int argc,char *argv[])
{
    const char *fname = "main()";
    setup(argc, argv);

    int traj = meas_arg.TrajStart;
    const int m_int = meas_arg.TrajIncrement;
    const int m_limit = meas_arg.TrajLessThanLimit;

    for(int conf = 0; conf < m_limit; ++conf) {
        load_checkpoint(traj);

        GnoneFbfm lat;

        const double ltw[4] = {
            l_twist_arg.mom[0], l_twist_arg.mom[1],
            l_twist_arg.mom[2], l_twist_arg.mom[3],
        };
        const double stw[4] = {
            s_twist_arg.mom[0], s_twist_arg.mom[1],
            s_twist_arg.mom[2], s_twist_arg.mom[3],
        };

        VRB.Result(cname, fname,
                   "l quark twist = %17.10e %17.10e %17.10e %17.10e\n",
                   ltw[0], ltw[1], ltw[2], ltw[3]);
        VRB.Result(cname, fname,
                   "s quark twist = %17.10e %17.10e %17.10e %17.10e\n",
                   stw[0], stw[1], stw[2], stw[3]);

        run_all(lat, stw, ltw, traj);
        traj += m_int;
    }

    VRB.Result(cname, fname, "Program ended normally.\n");
    End();
}

void load_checkpoint(int traj)
{
    const char *fname = "load_checkpoint()";

    char lat_file[256];
    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
    sprintf(lat_file, "%s.%d", meas_arg.GaugeStem, traj);
    QioArg rd_arg(lat_file, 0.001);
    rd_arg.ConcurIONumber = meas_arg.IOconcurrency;
    ReadLatticeParallel rl;
    rl.read(lat,rd_arg);
    if(!rl.good()) ERR.General(cname,fname,"Failed read lattice %s\n",lat_file);
    LatticeFactory::Destroy();
}

void setup(int argc, char *argv[])
{
    const char *fname = "setup()";

    Start(&argc, &argv);

    if(argc < 2) {
        ERR.General(cname, fname, "Must provide VML directory.\n");
    }

    if(chdir(argv[1]) != 0) {
        ERR.General(cname, fname, "Changing directory to %s failed.\n", argv[1]);
    }

    decode_vml_all();

    if(chdir(meas_arg.WorkDirectory) != 0) {
        ERR.General(cname, fname, "Changing directory to %s failed.\n", meas_arg.WorkDirectory);
    }
    VRB.Result(cname, fname, "Reading VML files successfully.\n");

    GJP.Initialize(do_arg);
    LRG.Initialize();

    init_bfm(&argc, &argv);
}

void init_bfm(int *argc, char **argv[])
{
    Chroma::initialize(argc, argv);
    multi1d<int> nrow(Nd);

    for(int i = 0; i< Nd; ++i)
        nrow[i] = GJP.Sites(i);

    Layout::setLattSize(nrow);
    Layout::create();

    // Fbfm::bfm_arg.solver = HtCayleyTanh;
    // Fbfm::bfm_arg.precon_5d = 0;
    Fbfm::bfm_arg.solver = DWF;
    Fbfm::bfm_arg.precon_5d = 1;
    // Fbfm::bfm_arg.solver = HmCayleyTanh;
    // Fbfm::bfm_arg.precon_5d = 0;

    Fbfm::bfm_arg.Ls = GJP.SnodeSites();
    Fbfm::bfm_arg.M5 = GJP.DwfHeight();
    Fbfm::bfm_arg.mass = 0.1;
    Fbfm::bfm_arg.residual = 1e-8;
    Fbfm::bfm_arg.max_iter = 10000;
    Fbfm::bfm_arg.Csw = 0.0;

    Fbfm::bfm_arg.node_latt[0] = QDP::Layout::subgridLattSize()[0];
    Fbfm::bfm_arg.node_latt[1] = QDP::Layout::subgridLattSize()[1];
    Fbfm::bfm_arg.node_latt[2] = QDP::Layout::subgridLattSize()[2];
    Fbfm::bfm_arg.node_latt[3] = QDP::Layout::subgridLattSize()[3];

    multi1d<int> procs = QDP::Layout::logicalSize();

    Fbfm::bfm_arg.local_comm[0] = procs[0] > 1 ? 0 : 1;
    Fbfm::bfm_arg.local_comm[1] = procs[1] > 1 ? 0 : 1;
    Fbfm::bfm_arg.local_comm[2] = procs[2] > 1 ? 0 : 1;
    Fbfm::bfm_arg.local_comm[3] = procs[3] > 1 ? 0 : 1;

    Fbfm::bfm_arg.ncoor[0] = 0;
    Fbfm::bfm_arg.ncoor[1] = 0;
    Fbfm::bfm_arg.ncoor[2] = 0;
    Fbfm::bfm_arg.ncoor[3] = 0;

    // mobius_scale = b + c in Andrew's notation
    bfmarg::mobius_scale = 2.;
    bfmarg::Threads(64);
    bfmarg::Reproduce(0);
    bfmarg::ReproduceChecksum(0);
    bfmarg::ReproduceMasterCheck(0);
    bfmarg::Verbose(0);

    Fbfm::use_mixed_solver = true;
}
