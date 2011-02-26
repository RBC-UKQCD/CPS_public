#include<config.h>

#include<alg/cg_arg.h>
#include<alg/hmd_arg.h>
#include<alg/do_arg.h>
#include<alg/w_spect_arg.h>
#include<alg/meas_arg.h>
#include<alg/pot_arg.h>
#include<alg/pbp_arg.h>
#include<alg/array_arg.h>
#include<alg/fix_gauge_arg.h>
#include<alg/s_spect_arg.h>

#include<alg/hmc_arg.h>
#include<alg/int_arg.h>
#include<alg/eig_arg.h>
#include<alg/qpropw_arg.h>
#include<alg/alg_nuc3pt.h>

//--------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

CgArg  cg_arg;
HmdArg hmd_arg;
EvoArg evo_arg;
DoArg  do_arg;
WspectArg  ws_arg;
WspectOutput ws_out;
PotArg pa;
PbpArg pbp;
FixGaugeArg fga;
MeasArg Meas;
MeasTask Task;
HmcArg hmc_arg;
ActionQuotientArg quo_arg;
ActionFermionArg frm_arg;
ActionBosonArg bsn_arg;
ActionRationalArg rat_arg;
ActionRationalSplitArg split_arg;
ActionRationalQuotientArg rat_quo_arg;
ActionGaugeArg gauge_arg;
IntABArg ab_arg;
EigArg eig_arg;
FloatArray float_array;
ParamArray param_array;
StagQuarkSrc stag_quark_src;
StagQuarkArg stag_quark_arg;
StagMesonArg stag_meson_arg;
StagMomMesonArg stag_mom_meson_arg;
StagNucleonArg stag_nucleon_arg;
StagNonLocalArg stag_non_local_arg;
NLStagMesonArg nlstag_meson_arg;
QPropWArg qpropw_arg;
Nuc3ptArg nuc3pt_arg;

int main(int argc, char *argv[])
{ 
  Start(&argc, &argv);

  Meas.TaskList.TaskList_len = 1;
  Meas.TaskList.TaskList_val = &Task;

  Meas.GaugeStem = "";
  Meas.RNGStem = "";

  Task.OutputFilestem = "";
  Task.ArgFilename = "arg_stem";

  ws_out.cg = "";
  ws_out.cg2 = "filename";
  ws_out.pbp = "filename";
  ws_out.mid_point = "filename";
  ws_out.nucleon = "filename";
  ws_out.nucleon_prime = "filename";
  ws_out.delta_x = "filename";
  ws_out.delta_y = "filename";
  ws_out.delta_z = "filename";
  ws_out.delta_t = "filename";

  evo_arg.ensemble_id = "id";
  evo_arg.ensemble_label = "label";
  evo_arg.creator = "creator";
  evo_arg.gauge_file_stem = "ckpoint_lat";
  evo_arg.rng_file_stem = "ckpoint_rng";
  evo_arg.plaquette_stem = "plaq";
  evo_arg.evo_stem = "hmd";
  evo_arg.work_directory = "";

  ws_arg.CgArgFile = "CgFile";

  rat_arg.resize(4);
  rat_arg.resize(0, 10, 16);
  rat_arg.resize(1, 10, 16);
  rat_arg.resize(2, 9, 14);
  rat_arg.resize(3, 6, 9);
  bsn_arg.resize(1);
  frm_arg.resize(1);
  split_arg.resize(2);
  quo_arg.resize(1);
  rat_quo_arg.resize(1);
  rat_quo_arg.resize(0, 9, 15, 6, 9);
  float_array.resize(4);
  param_array.resize(3);
  param_array.params.params_val[0].name="mass";
  param_array.params.params_val[1].name="epsilon";
  param_array.params.params_val[2].name="u0";
  stag_quark_src.type = S_QUARK_POINT;
  for(int i = 0; i < 4; i++)
    {
      stag_quark_src.origin[i] = 0;
      stag_quark_src.end[i] = 0;
    }
  stag_quark_src.dir = HDM_T;
  stag_quark_arg.qid = 0;
  stag_quark_arg.cg.mass = 0.1;
  stag_quark_arg.cg.max_num_iter = 5000;
  stag_quark_arg.cg.stop_rsd = 1e-10;
  stag_quark_arg.cg.true_rsd = 1e-10;
  stag_quark_arg.cg.RitzMatOper = MATPCDAG_MATPC;
  stag_quark_arg.cg.Inverter = CG;
  stag_quark_arg.cg.bicgstab_n = 1;
  stag_quark_arg.src = stag_quark_src;
  stag_quark_arg.sln = LOCAL;
  stag_meson_arg.qid0 = 0;
  stag_meson_arg.qid1 = 1;
  stag_meson_arg.dir = HDM_T;
  stag_meson_arg.meson_buf = 0;
  stag_mom_meson_arg.qid0 = 0;
  stag_mom_meson_arg.qid1 = 1;
  stag_mom_meson_arg.dir = HDM_T;
  stag_mom_meson_arg.no_of_momenta = 4;
  stag_mom_meson_arg.meson_buf = 0;
  stag_nucleon_arg.qid0 = 0;
  stag_nucleon_arg.qid1 = 1;
  stag_nucleon_arg.qid2 = 2;
  stag_nucleon_arg.dir = HDM_T;
  stag_nucleon_arg.nucleon_buf = 0;
  stag_non_local_arg.qid0 = 0;
  stag_non_local_arg.qid1 = 1;
  stag_non_local_arg.qid2 = 2;
  stag_non_local_arg.dir = HDM_T;
  stag_non_local_arg.nlocal_buf = 0;
  for(int i = 0; i < 8; i++)
    nlstag_meson_arg.qid0[i] = i;
  nlstag_meson_arg.dir = HDM_T;
  nlstag_meson_arg.nlocal_buf = 0;
  
  qpropw_arg.Encode("qpropw_arg.vml","qpropw_arg");

  nuc3pt_arg.cname="AlgNuc3pt";
  nuc3pt_arg.ensemble_label="AlgNuc3pt";
  nuc3pt_arg.Encode("nuc3pt_arg.vml","nuc3pt_arg");


  End();

  return(0);
}
