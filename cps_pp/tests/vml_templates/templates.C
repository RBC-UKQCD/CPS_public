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

#include<alg/hmc_arg.h>
#include<alg/int_arg.h>
#include<alg/eig_arg.h>

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


int main(int argc, char *argv[])
{ 

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

  cg_arg.Encode("cg_arg.vml","cg_arg");
  hmd_arg.Encode("hmd_arg.vml","hmd_arg");
  evo_arg.Encode("evo_arg.vml","evo_arg");
  do_arg.Encode("do_arg.vml","do_arg");
  ws_arg.Encode("w_spect_arg.vml","w_spect_arg");
  ws_out.Encode("w_spect_output.vml","w_spect_output");
  Meas.Encode("meas_arg.vml","meas_arg");
  pa.Encode("pot_arg.vml","pot_arg");
  pbp.Encode("pbp_arg.vml","pbp_arg");
  fga.Encode("fix_gauge_arg.vml","fix_gauge_arg");
  hmc_arg.Encode("hmc_arg.vml","hmc_arg");
  quo_arg.Encode("quo_arg.vml","quo_arg");
  frm_arg.Encode("frm_arg.vml","frm_arg");
  bsn_arg.Encode("bsn_arg.vml","bsn_arg");
  rat_arg.Encode("rat_arg.vml","rat_arg");
  rat_quo_arg.Encode("rat_quo_arg.vml","rat_quo_arg");
  split_arg.Encode("split_arg.vml","split_arg");
  gauge_arg.Encode("gauge_arg.vml","gauge_arg");
  ab_arg.Encode("ab_arg.vml","ab_arg");
  eig_arg.Encode("eig_arg.vml","eig_arg");
  float_array.Encode("float_array.vml","float_array");

  vml_markup_type(VML_XML);
  
  cg_arg.Encode("cg_arg.xml","cg_arg");
  hmd_arg.Encode("hmd_arg.xml","hmd_arg");
  evo_arg.Encode("evo_arg.xml","evo_arg");
  do_arg.Encode("do_arg.xml","do_arg");
  ws_arg.Encode("w_spect_arg.xml","w_spect_arg");
  ws_out.Encode("w_spect_output.xml","w_spect_output");
  hmc_arg.Encode("hmc_arg.xml","hmc_arg");
  quo_arg.Encode("quo_arg.xml","quo_arg");
  frm_arg.Encode("frm_arg.xml","frm_arg");
  bsn_arg.Encode("bsn_arg.xml","bsn_arg");
  rat_arg.Encode("rat_arg.xml","rat_arg");
  rat_quo_arg.Encode("rat_quo_arg.xml","rat_quo_arg");
  split_arg.Encode("split_arg.xml","split_arg");
  gauge_arg.Encode("gauge_arg.xml","gauge_arg");
  ab_arg.Encode("ab_arg.xml","ab_arg");
  eig_arg.Encode("eig_arg.xml","eig_arg");

  return(0);
}
