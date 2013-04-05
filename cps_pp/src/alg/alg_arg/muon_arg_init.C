// Initialize the muon argument structure
#include <alg/muon_arg.h>
CPS_START_NAMESPACE
MuonArg::MuonArg():
  u1lat(0),
  source_time(0),
  oper_time_start(0),
  oper_time_end(0),
  operator_gamma(0),
  n_source(1),
  source_inc(1),
  loop_mass(0.1),
  line_mass(0.1),
  charge(0),
  GaugeFix(0),
  tINC(1),
  ptINC(1),
  conf(0),
  vp_kind(CONSERVED_LOCAL)
{;}
CPS_END_NAMESPACE
