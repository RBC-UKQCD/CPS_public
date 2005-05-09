#include<config.h>

#include<alg/cg_arg.h>
#include<alg/hmd_arg.h>
#include<alg/do_arg.h>
#include<alg/w_spect_arg.h>

//--------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

CgArg  cg_arg;
HmdArg hmd_arg;
EvoArg evo_arg;
DoArg  do_arg;
WspectArg  ws_arg;
WspectOutput ws_out;
RhmcPolesState rhmc_poles;

int main(int argc, char *argv[])
{ 

  ws_out.cg = "filename";
  ws_out.cg2 = "filename";
  ws_out.pbp = "filename";
  ws_out.mid_point = "filename";
  ws_out.a0_p = "filename";
  ws_out.a1 = "filename";
  ws_out.b1 = "filename";
  ws_out.pion = "filename";
  ws_out.pion_prime = "filename";
  ws_out.rho = "filename";
  ws_out.rho_prime = "filename";
  ws_out.a0 = "filename";
  ws_out.a0_prime = "filename";
  ws_out.a1_x = "filename";
  ws_out.a1_y = "filename";
  ws_out.a1_z = "filename";
  ws_out.b1_x = "filename";
  ws_out.b1_y = "filename";
  ws_out.b1_z = "filename";
  ws_out.rho_x = "filename";
  ws_out.rho_y = "filename";
  ws_out.rho_z = "filename";
  ws_out.rho_x_prime = "filename";
  ws_out.rho_y_prime = "filename";
  ws_out.rho_z_prime = "filename";
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

    cg_arg.Encode("cg_arg.vml","cg_arg");
    hmd_arg.Encode("hmd_arg.vml","hmd_arg");
    evo_arg.Encode("evo_arg.vml","evo_arg");
    do_arg.Encode("do_arg.vml","do_arg");
    ws_arg.Encode("w_spect_arg.vml","w_spect_arg");
    ws_out.Encode("w_spect_output.vml","w_spect_output");
    rhmc_poles.Encode("rhmc_poles_state.vml","rhmc_poles_state");
    
    vml_markup_type(VML_XML);
    
    cg_arg.Encode("cg_arg.xml","cg_arg");
    hmd_arg.Encode("hmd_arg.xml","hmd_arg");
    evo_arg.Encode("evo_arg.xml","evo_arg");
    do_arg.Encode("do_arg.xml","do_arg");
    ws_arg.Encode("w_spect_arg.xml","w_spect_arg");
    ws_out.Encode("w_spect_output.xml","w_spect_output");
    rhmc_poles.Encode("rhmc_poles_state.xml","rhmc_poles_state");

     
 return(0);
}
