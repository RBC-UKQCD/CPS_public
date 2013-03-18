// -*- mode:c++; c-basic-offset:4 -*-

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#include <util/lattice.h>
#include <util/lattice/fbfm.h>

#include <alg/array_arg.h>
#include <alg/common_arg.h>
#include <alg/eigcg_arg.h>
#include <alg/alg_fix_gauge.h>
#include <alg/qpropw_arg.h>
#include <alg/meas_arg.h>

#include <util/gjp.h>
#include <util/time_cps.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/qcdio.h>
#include <util/ReadLatticePar.h>

#include <string>
#include <vector>
#include <cassert>

#include <qmp.h>

#include "prop_container.h"
#include "eigcg.h"
#include "my_util.h"
#include "run_kl3.h"
#include "run_2pion.h"
#include "run_k2pipi.h"
#include "run_bk.h"
#include "run_mres.h"
#include "run_meson.h"
#include "run_omega.h"
#include "run_prop.h"
#include "twisted_bc.h"

static const char *cname = "";

USING_NAMESPACE_CPS
using namespace std;

DoArg do_arg;
MeasArg meas_arg;
QPropWArg lqpropw_arg;
QPropWArg sqpropw_arg;
FixGaugeArg fix_gauge_arg;
EigCGArg l_eigcg_arg;

// Integer arrays, for setting time locations of exact propagators
IntArray l_ut_loc; // untwisted light
IntArray l_tw_loc; // twisted light 
IntArray s_ut_loc; // untwisted strange
IntArray s_tw_loc; // twisted strange

// d quark momentum for K -> pi pi, for 48^3 this should be {1, 1, 1}.
IntArray d_mom_kpp;

// l and s quark twists
FloatArray l_twist_arg;
FloatArray s_twist_arg;

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
    decode_vml(sqpropw_arg);
    decode_vml(fix_gauge_arg);
    decode_vml(l_eigcg_arg);

    decode_vml(l_ut_loc);
    decode_vml(l_tw_loc);
    decode_vml(s_ut_loc);
    decode_vml(s_tw_loc);

    decode_vml(d_mom_kpp);

    decode_vml(l_twist_arg);
    decode_vml(s_twist_arg);
}

void load_checkpoint(int traj);
void init_bfm(int *argc, char **argv[]);
void setup(int argc, char *argv[]);

void run_contractions(const AllProp &sprop, const AllProp &stwst,
                      const AllProp &lprop, const AllProp &ltwst,
                      const string &rdir,
                      int traj, PROP_TYPE ptype)
{
    const char *fname = "run_contractions()";

    const string trajs = string(".") + tostring(traj);
    
    //////////////////////////////////////////////////////////////////////
    // 2. meson contractions

    // pion and kaon
    run_meson_pt(lprop, lprop, GAMMA_5, GAMMA_5, rdir + "/pion-00WP" + trajs, ptype);
    run_meson_pt(lprop, ltwst, GAMMA_5, GAMMA_5, rdir + "/pion-01WP" + trajs, ptype);
    run_meson_pt(sprop, lprop, GAMMA_5, GAMMA_5, rdir + "/kaon-00WP" + trajs, ptype);
    run_meson_pt(stwst, lprop, GAMMA_5, GAMMA_5, rdir + "/kaon-10WP" + trajs, ptype);
    run_meson_pt(sprop, ltwst, GAMMA_5, GAMMA_5, rdir + "/kaon-01WP" + trajs, ptype);
    
    run_meson_wall(lprop, lprop, GAMMA_5, GAMMA_5, rdir + "/pion-00WW" + trajs, ptype);
    run_meson_wall(lprop, ltwst, GAMMA_5, GAMMA_5, rdir + "/pion-01WW" + trajs, ptype);
    run_meson_wall(sprop, lprop, GAMMA_5, GAMMA_5, rdir + "/kaon-00WW" + trajs, ptype);
    run_meson_wall(stwst, lprop, GAMMA_5, GAMMA_5, rdir + "/kaon-10WW" + trajs, ptype);
    run_meson_wall(sprop, ltwst, GAMMA_5, GAMMA_5, rdir + "/kaon-01WW" + trajs, ptype);

    // scalar meson (sigma) contractions.
    run_meson_wall(lprop, lprop, ID,      ID,      rdir + "/sigma-00WW" + trajs, ptype);
    run_meson_disc(lprop, lprop, ID,      ID,      rdir + "/sigma-dis-00WW" + trajs, ptype);

    // eta eta' contractions
    //
    // We share the light-light propagator with pion contractions.
    run_meson_wall(sprop, sprop, GAMMA_5, GAMMA_5, rdir + "/ss-00WW" + trajs, ptype);
    run_meson_disc(lprop, sprop, GAMMA_5, GAMMA_5, rdir + "/ls-dis-00WW" + trajs, ptype);

    // f_K and f_pi measurements
    run_meson_pt(lprop, lprop, GAMMA_35, GAMMA_5, rdir + "/fp-00WP" + trajs, ptype);
    run_meson_pt(sprop, lprop, GAMMA_35, GAMMA_5, rdir + "/fk-00WP" + trajs, ptype);

    run_meson_pt(lprop, lprop, GAMMA_5, GAMMA_35, rdir + "/fpr-00WP" + trajs, ptype);
    run_meson_pt(sprop, lprop, GAMMA_5, GAMMA_35, rdir + "/fkr-00WP" + trajs, ptype);

    run_meson_pt(lprop, lprop, GAMMA_35, GAMMA_35, rdir + "/ap-00WP" + trajs, ptype);
    run_meson_pt(sprop, lprop, GAMMA_35, GAMMA_35, rdir + "/ak-00WP" + trajs, ptype);

    run_meson_wall(lprop, lprop, GAMMA_35, GAMMA_5, rdir + "/fp-00WW" + trajs, ptype);
    run_meson_wall(sprop, lprop, GAMMA_35, GAMMA_5, rdir + "/fk-00WW" + trajs, ptype);

    // rho meson
    run_meson_pt(lprop, lprop, GAMMA_0, GAMMA_0, rdir + "/rho-x-00WP" + trajs, ptype);
    run_meson_pt(lprop, lprop, GAMMA_1, GAMMA_1, rdir + "/rho-y-00WP" + trajs, ptype);
    run_meson_pt(lprop, lprop, GAMMA_2, GAMMA_2, rdir + "/rho-z-00WP" + trajs, ptype);
    run_meson_wall(lprop, lprop, GAMMA_0, GAMMA_0, rdir + "/rho-x-00WW" + trajs, ptype);
    run_meson_wall(lprop, lprop, GAMMA_1, GAMMA_1, rdir + "/rho-y-00WW" + trajs, ptype);
    run_meson_wall(lprop, lprop, GAMMA_2, GAMMA_2, rdir + "/rho-z-00WW" + trajs, ptype);

    //////////////////////////////////////////////////////////////////////
    // 2. Omega baryon
    run_omega_pt(sprop, GAMMA_0, rdir + "/sss-x-00WP" + trajs, ptype);
    run_omega_pt(sprop, GAMMA_1, rdir + "/sss-y-00WP" + trajs, ptype);
    run_omega_pt(sprop, GAMMA_2, rdir + "/sss-z-00WP" + trajs, ptype);
    run_omega_pt(sprop, GAMMA_3, rdir + "/sss-t-00WP" + trajs, ptype);
    run_omega_pt(sprop, GAMMA_5, rdir + "/sss-5-00WP" + trajs, ptype);

    //////////////////////////////////////////////////////////////////////
    // 3. Kl3
    run_kl3(sprop, lprop, lprop, rdir + "/kl3-00" + trajs, ptype);
    run_kl3(sprop, lprop, ltwst, rdir + "/kl3-01" + trajs, ptype);
    run_kl3(stwst, lprop, lprop, rdir + "/kl3-10" + trajs, ptype);

    // The following contractions are used for Z_V.
    // Useful for the double ratio method or the UKQCD method.
    run_kl3(lprop, lprop, lprop, rdir + "/zpa-00" + trajs, ptype);
    run_kl3(sprop, lprop, sprop, rdir + "/zka-00" + trajs, ptype);
    run_kl3(lprop, sprop, lprop, rdir + "/zkb-00" + trajs, ptype);
    run_kl3(sprop, sprop, sprop, rdir + "/zss-00" + trajs, ptype);

    // Since line 1 and 3 both carry momentum, are their directions
    // consistent?
    // I think they are.
    run_kl3(ltwst, lprop, ltwst, rdir + "/zpa-11" + trajs, ptype);
    run_kl3(lprop, ltwst, lprop, rdir + "/zpb-11" + trajs, ptype);
    run_kl3(stwst, lprop, stwst, rdir + "/zka-11" + trajs, ptype);
    run_kl3(lprop, stwst, lprop, rdir + "/zkb-11" + trajs, ptype);

    //////////////////////////////////////////////////////////////////////
    // 4. Bk
    run_bk(lprop, sprop, lprop, sprop, rdir + "/bk" + trajs, ptype);
}

void run_k2pipi_contractions(const AllProp &sprop,
                             const AllProp &uprop,
                             const AllProp &dprop,
                             const string &rdir,
                             int traj, const int mom[3],
                             PROP_TYPE ptype)
{
    const string trajs = string(".") + tostring(traj);

    // zero momentum pi pi scattering
    int zmom[3] = {0, 0, 0};
    run_2pionDC(uprop, uprop, rdir + "/2pion000" + trajs, ptype, zmom);

    for(unsigned i = 0; i < 8; ++i) {
        int p[3];
        p[0] = (i & 1) ? -mom[0] : mom[0];
        p[1] = (i & 2) ? -mom[1] : mom[1];
        p[2] = (i & 4) ? -mom[2] : mom[2];

        const string fn = rdir + "/2pion"
            + tostring(p[0]) + tostring(p[1]) + tostring(p[2])
            + trajs;

        run_2pionDC(uprop, dprop, fn, ptype, p);
    }

    // K->pipi without momentum.
    run_k2pipi(sprop, uprop, uprop, rdir + "/k2pipi-0" + trajs, ptype);

    // K->pipi with momentum.
    run_k2pipi(sprop, uprop, dprop, rdir + "/k2pipi-1" + trajs, ptype);
}

class FixGauge
{
public:
    FixGauge(cps::Lattice &lat,
             cps::FixGaugeArg &fix_gauge_arg,
             int traj)
    {
        char buf[256];
        sprintf(buf, "../results/fg-bc.%d", traj);
        cps::Fclose(cps::Fopen(buf, "w"));
        com_fg.set_filename(buf);
        
        fg = new cps::AlgFixGauge(lat, &com_fg, &fix_gauge_arg);
        fg->run();
    }

    ~FixGauge() {
        fg->free();
        delete fg;
    }
private:
    cps::CommonArg com_fg;
    cps::AlgFixGauge *fg;
};

// stw: twisting angle of the strange quark (connecting the operator
// and the kaon).
//
// ltw: twisting angle of the light quark (connecting the operator and
// the pion).
void run_all(Lattice &lat,
             const double stw[4], // strange quark twists, for Kl3
             const double ltw[4], //   light quark twists, for Kl3
             const int mom[3],    // momentum for the d quark, used by K2pipi
             int traj)
{
    const char *fname = "run_all()";

    //////////////////////////////////////////////////////////////////////
    // 1. props: strange, strange twisted, light, light twisted

    // exact strange
    AllProp sprop(AllProp::DOUBLE), stwst(AllProp::DOUBLE);
    // exact light
    AllProp lprop_e(AllProp::DOUBLE), ltwst_e(AllProp::DOUBLE);
    // sloppy light, single precision
    AllProp lprop(AllProp::SINGLE), ltwst(AllProp::SINGLE);

    Float dtime0 = dclock();

    FixGauge fg(lat, fix_gauge_arg, traj);

    Float dtime1 = dclock();

    // l untwisted
    run_wall_prop(&lprop_e, &lprop, l_ut_loc, lat, lqpropw_arg, &l_eigcg_arg, traj, true );

    // l twisted
    twisted_bc(lat, ltw, true);
    run_wall_prop(&ltwst_e, &ltwst, l_tw_loc, lat, lqpropw_arg, &l_eigcg_arg, traj, false);
    twisted_bc(lat, ltw, false);

    // s untwisted
    run_wall_prop(NULL, &sprop, s_ut_loc, lat, sqpropw_arg, NULL, traj, true );

    // s twisted
    // twisted_bc(lat, stw, true);
    // run_wall_prop(NULL, &stwst, s_tw_loc, lat, sqpropw_arg, NULL, traj, false);
    // twisted_bc(lat, stw, false);

    Float dtime2 = dclock();

    run_contractions(sprop, stwst, lprop_e, ltwst_e, "../resultsEPA", traj, PROP_PA);
    run_contractions(sprop, stwst, lprop_e, ltwst_e, "../resultsEP",  traj, PROP_P);
    run_contractions(sprop, stwst, lprop_e, ltwst_e, "../resultsEA",  traj, PROP_A);

    run_contractions(sprop, stwst, lprop, ltwst, "../resultsPA", traj, PROP_PA);
    run_contractions(sprop, stwst, lprop, ltwst, "../resultsP",  traj, PROP_P);
    run_contractions(sprop, stwst, lprop, ltwst, "../resultsA",  traj, PROP_A);

    Float dtime3 = dclock();

    ////////////////////////////////////////////////////////////////////////
    // I=2 K to pi pi
    // free unwanted propagators to save some memory.
    ltwst_e.clear();
    ltwst.clear();
    stwst.clear();

    // twisted light for K -> pi pi
    run_mom_prop(&ltwst_e, &ltwst, l_tw_loc, lat, lqpropw_arg, &l_eigcg_arg, traj, mom);

    Float dtime4 = dclock();

    run_k2pipi_contractions(sprop, lprop_e, ltwst_e, "../resultsEPA", traj, mom, PROP_PA);
    run_k2pipi_contractions(sprop, lprop_e, ltwst_e, "../resultsEP",  traj, mom, PROP_P);
    run_k2pipi_contractions(sprop, lprop_e, ltwst_e, "../resultsEA",  traj, mom, PROP_A);

    run_k2pipi_contractions(sprop, lprop, ltwst, "../resultsPA", traj, mom, PROP_PA);
    run_k2pipi_contractions(sprop, lprop, ltwst, "../resultsP",  traj, mom, PROP_P);
    run_k2pipi_contractions(sprop, lprop, ltwst, "../resultsA",  traj, mom, PROP_A);

    Float dtime5 = dclock();

    VRB.Result(cname, fname, "fix gauge    = %17.10e seconds\n", dtime1 - dtime0);
    VRB.Result(cname, fname, "kl3 prop     = %17.10e seconds\n", dtime2 - dtime1);
    VRB.Result(cname, fname, "kl3          = %17.10e seconds\n", dtime3 - dtime2);
    VRB.Result(cname, fname, "k2pipi prop  = %17.10e seconds\n", dtime4 - dtime3);
    VRB.Result(cname, fname, "k2pipi       = %17.10e seconds\n", dtime5 - dtime4);
    VRB.Result(cname, fname, "total        = %17.10e seconds\n", dtime5 - dtime0);

    //////////////////////////////////////////////////////////////////////
    // store propagators
    // lprop_e.store_all("lprop_raw_", lqpropw_arg.cg.mass, traj);
    
    // Float dtime6 = dclock();

    // VRB.Result(cname, fname, "store prop   = %17.10e seconds\n", dtime6 - dtime5);
}

// Shift the locations of where exact propagators are calculated. This
// is done by shifting the following arrays,
//
// l_ut_loc, l_tw_loc
// s_ut_loc, s_tw_loc
void do_shift(int traj)
{
    const int t_size = GJP.TnodeSites() * GJP.Tnodes();

    // We shift the exact solutions by a random amount uniformly
    // distributed between [0,T).
    int shift = drand48() * t_size;

    // Make sure we have the same number on all nodes.
    QMP_broadcast(&shift, sizeof(int));
    assert(shift >= 0 && shift < t_size);

    static int shift_acc = 0;

    shift_acc = (shift_acc + shift) % t_size;
    // printf("traj = %d, Shift on the exact propagators = %d\n", traj, shift_acc);
    VRB.Result(cname, "do_shift()",
               "traj = %d, Shift on the exact propagators = %d\n",
               traj, shift_acc);

    for(unsigned i = 0; i < l_ut_loc.v.v_len; ++i) {
        l_ut_loc.v.v_val[i] = (l_ut_loc.v.v_val[i] + shift) % t_size;
    }
    for(unsigned i = 0; i < l_tw_loc.v.v_len; ++i) {
        l_tw_loc.v.v_val[i] = (l_tw_loc.v.v_val[i] + shift) % t_size;
    }
    for(unsigned i = 0; i < s_ut_loc.v.v_len; ++i) {
        s_ut_loc.v.v_val[i] = (s_ut_loc.v.v_val[i] + shift) % t_size;
    }
    for(unsigned i = 0; i < s_tw_loc.v.v_len; ++i) {
        s_tw_loc.v.v_val[i] = (s_tw_loc.v.v_val[i] + shift) % t_size;
    }
}

int main(int argc,char *argv[])
{
    const char *fname = "main()";
    srand48(time(NULL));
    setup(argc, argv);

    int traj = meas_arg.TrajStart;
    const int m_int = meas_arg.TrajIncrement;
    const int m_limit = meas_arg.TrajLessThanLimit;
    for(int conf = 0; conf < m_limit; ++conf) {
        // shift the exact propagators
        do_shift(traj);
        load_checkpoint(traj);

        GnoneFbfm lat;

        // NOTE: there are 4 twists but 3 momenta. This is just
        // because of how code is written. The t twist is normally
        // zero.
        assert(l_twist_arg.Floats.Floats_len == 4);
        assert(s_twist_arg.Floats.Floats_len == 4);
        assert(d_mom_kpp.v.v_len == 3);

        const double *ltw = l_twist_arg.Floats.Floats_val;
        const double *stw = s_twist_arg.Floats.Floats_val;
        const int *dmom = d_mom_kpp.v.v_val;

        VRB.Result(cname, fname,
                   "l quark twist (kl3) = %17.10e %17.10e %17.10e %17.10e\n",
                   ltw[0], ltw[1], ltw[2], ltw[3]);
        VRB.Result(cname, fname,
                   "s quark twist (kl3) = %17.10e %17.10e %17.10e %17.10e\n",
                   stw[0], stw[1], stw[2], stw[3]);
        VRB.Result(cname, fname,
                   "d quark mom  (k2pp) = %d %d %d\n",
                   dmom[0], dmom[1], dmom[2]);

        run_all(lat, stw, ltw, dmom, traj);
        traj += m_int;
    }

    VRB.Result(cname, fname, "Program ended normally.\n");
    End();
}

void load_checkpoint(int traj)
{
    const char *fname = "load_checkpoint()";

    char lat_file[256];
    GnoneFnone lat;

    sprintf(lat_file, "%s.%d", meas_arg.GaugeStem, traj);
    QioArg rd_arg(lat_file, 0.001);
    rd_arg.ConcurIONumber = meas_arg.IOconcurrency;
    ReadLatticeParallel rl;
    rl.read(lat,rd_arg);
    if(!rl.good()) ERR.General(cname,fname,"Failed read lattice %s\n",lat_file);
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
    bfmarg::Verbose(1);

    Fbfm::use_mixed_solver = true;
}
