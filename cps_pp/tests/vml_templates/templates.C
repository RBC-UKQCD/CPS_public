#include<config.h>

#include<alg/cg_arg.h>
#include<alg/hmd_arg.h>
#include<alg/do_arg.h>
#include<alg/w_spect_arg.h>
#include<alg/meas_arg.h>
#include<alg/pot_arg.h>
#include<alg/pbp_arg.h>
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
RhmcPolesState rhmc_poles;
MeasArg Meas;
MeasTask Task;
HmcArg hmc_arg;
ActionFermionArg frm_arg;
ActionBosonArg bsn_arg;
ActionRationalArg rat_arg;
ActionRationalSplitArg split_arg;
ActionGaugeArg gauge_arg;
IntABArg ab_arg;
EigArg eig_arg;

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

  hmd_arg.rhmc_poles_file = "filename";

  evo_arg.ensemble_id = "id";
  evo_arg.ensemble_label = "label";
  evo_arg.creator = "creator";
  evo_arg.gauge_file_stem = "ckpoint_lat";
  evo_arg.rng_file_stem = "ckpoint_rng";
  evo_arg.plaquette_stem = "plaq";
  evo_arg.evo_stem = "hmd";
  evo_arg.work_directory = "";

  ws_arg.CgArgFile = "CgFile";

  rat_arg.resize(2);
  rat_arg.resize(0, 10, 16);
  rat_arg.resize(1, 7, 13);
  bsn_arg.resize(0);
  frm_arg.resize(1);
  split_arg.resize(2);

  double masses[3];
  masses[0] =0.1;
  masses[1] =0.2;
  masses[2] =0.3;
  eig_arg.Mass.Mass_len=3;
  eig_arg.Mass.Mass_val=masses;
  

  cg_arg.Encode("cg_arg.vml","cg_arg");
  hmd_arg.Encode("hmd_arg.vml","hmd_arg");
  evo_arg.Encode("evo_arg.vml","evo_arg");
  do_arg.Encode("do_arg.vml","do_arg");
  ws_arg.Encode("w_spect_arg.vml","w_spect_arg");
  ws_out.Encode("w_spect_output.vml","w_spect_output");
  rhmc_poles.Encode("rhmc_poles_state.vml","rhmc_poles_state");
  Meas.Encode("meas_arg.vml","meas_arg");
  pa.Encode("pot_arg.vml","pot_arg");
  pbp.Encode("pbp_arg.vml","pbp_arg");
  fga.Encode("fix_gauge_arg.vml","fix_gauge_arg");
  hmc_arg.Encode("hmc_arg.vml","hmc_arg");
  frm_arg.Encode("frm_arg.vml","frm_arg");
  bsn_arg.Encode("bsn_arg.vml","bsn_arg");
//  rat_arg.Encode("rat_arg.vml","rat_arg");
  split_arg.Encode("split_arg.vml","split_arg");
  gauge_arg.Encode("gauge_arg.vml","gauge_arg");
  ab_arg.Encode("ab_arg.vml","ab_arg");
  eig_arg.Encode("eig_arg.vml","eig_arg");
#if 0

  vml_markup_type(VML_XML);
  
  cg_arg.Encode("cg_arg.xml","cg_arg");
  hmd_arg.Encode("hmd_arg.xml","hmd_arg");
  evo_arg.Encode("evo_arg.xml","evo_arg");
  do_arg.Encode("do_arg.xml","do_arg");
  ws_arg.Encode("w_spect_arg.xml","w_spect_arg");
  ws_out.Encode("w_spect_output.xml","w_spect_output");
  rhmc_poles.Encode("rhmc_poles_state.xml","rhmc_poles_state");
  hmc_arg.Encode("hmc_arg.xml","hmc_arg");
  frm_arg.Encode("frm_arg.xml","frm_arg");
  bsn_arg.Encode("bsn_arg.xml","bsn_arg");
  rat_arg.Encode("rat_arg.xml","rat_arg");
  split_arg.Encode("split_arg.xml","split_arg");
  gauge_arg.Encode("gauge_arg.xml","gauge_arg");
  ab_arg.Encode("ab_arg.xml","ab_arg");
  eig_arg.Encode("eig_arg.xml","eig_arg");
#endif

  return(0);
}
