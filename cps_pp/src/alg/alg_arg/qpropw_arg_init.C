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
  save_ls_prop(0),
  do_half_fermion(0),
//  bstart(0),
//  bend(0),
//  gauss_N(30),
//  gauss_W(4.35), 
  SeqSmearSink(POINT)
{};

QPropWPointArg::QPropWPointArg()
: x(0),y(0),z(0)
{};

QPropWBoxArg::QPropWBoxArg()
:box_start(0),box_end(0)
{};

// added by Hantao
QPropW4DBoxArg::QPropW4DBoxArg()
{}

QPropWRandArg::QPropWRandArg()
:rng(NORAND),seed(1111)
{};

#if 1
QPropWExpArg::QPropWExpArg()
:exp_A(1.2),exp_B(0.1),exp_C(8)
{};
#endif

QPropWGaussArg::QPropWGaussArg()
  :gauss_N(30),gauss_W(4.35),
   gauss_link_smear_type(GKLS_NONE),
   gauss_link_smear_N(0),
  gauss_link_smear_coeff(0)
{};

CPS_END_NAMESPACE
