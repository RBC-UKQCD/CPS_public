#include <config.h>
#include <alg/qpropw_arg.h>
CPS_START_NAMESPACE

QPropWArg::QPropWArg()
: 
  file(NULL),
  t(0),
  gauge_fix_src(0),
  gauge_fix_snk(0),
  store_midprop(0),
  save_prop(0),
  do_half_fermion(0)
{}

QPropWPointArg::QPropWPointArg()
: x(0),y(0),z(0)
{}

QPropWBoxArg::QPropWBoxArg()
:box_start(0),box_end(0)
{}

QPropWRandArg::QPropWRandArg()
:rng(NORAND),seed(1111)
{}

#if 1
QPropWExpArg::QPropWExpArg()
:exp_A(1.2),exp_B(0.1),exp_C(8)
{}
#endif

CPS_END_NAMESPACE
