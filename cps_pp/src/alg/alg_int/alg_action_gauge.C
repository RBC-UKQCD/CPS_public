#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_action_gauge.C
//
// AlgActionGauge is represents the pure gauge contribution to the QCD
// action
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<alg/alg_hmd.h>
#include<util/lattice.h>
#include<util/vector.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/time_cps.h>
#include<alg/alg_int.h>
CPS_START_NAMESPACE

AlgActionGauge::AlgActionGauge(AlgMomentum &mom, ActionGaugeArg &g_arg)
    : AlgAction(mom, g_arg.action_arg), cname("AlgActionGauge")
{
  int_type = INT_GAUGE;
  gauge_arg = &g_arg;
  gluon = g_arg.gluon;

}

AlgActionGauge::~AlgActionGauge() {

}

//!< Heat Bath for the gauge action (i.e., does nothing)
void AlgActionGauge::heatbath() {

}

//!< Calculate gauge contribution to the Hamiltonian
Float AlgActionGauge::energy() {
  Float dtime = -dclock();

  char *fname = "energy()";
  Lattice &lat = 
    LatticeFactory::Create(F_CLASS_NONE, gluon);
  Float h = lat.GhamiltonNode();
  LatticeFactory::Destroy();

  dtime += dclock();
  print_flops(cname, fname, 0, dtime);

  return h;
}

void AlgActionGauge::prepare_fg(Matrix * force, Float dt_ratio)
{
  Float dtime = -dclock();
  const char *fname = "prepare_fg(M*,F)";

  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, gluon);
  Fdt = lat.EvolveMomGforce(force, dt_ratio);
  if (force_measure == FORCE_MEASURE_YES) {
    char label[200];
    sprintf(label, "%s:", force_label);
    Fdt.print(dt_ratio, label);
  }
  LatticeFactory::Destroy();

  dtime += dclock();
  print_flops(cname, fname, 0, dtime);
}

//!< evolve method evolves the momentum due to the gauge force
void AlgActionGauge::evolve(Float dt, int steps) 
{
  Float dtime = -dclock();
  const char *fname = "evolve()";
  //!< Create an appropriate lattice
  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, gluon);

  for (int i=0; i<steps; i++) {
    Fdt = lat.EvolveMomGforce(mom, dt);

    if (force_measure == FORCE_MEASURE_YES) {
      char label[200];
      sprintf(label, "%s:", force_label);
      Fdt.print(dt, label);
    }
  }

  LatticeFactory::Destroy();
  dtime += dclock();
  print_flops(cname, fname, 0, dtime);
}

//!< Dummy methods
void AlgActionGauge::cost(CgStats *cg_stats_global) {

}

void AlgActionGauge::init() {

}

CPS_END_NAMESPACE
