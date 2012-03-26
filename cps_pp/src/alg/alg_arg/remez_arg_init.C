#include <alg/remez_arg.h>
CPS_START_NAMESPACE
RemezArg::RemezArg()
:
approx_type(RATIONAL_APPROX_POWER),
valid_approx(0),
lambda_low(0),
lambda_high(0),
degree(0),
field_type(FERMION),
power_num(0),
power_den(0),
error(0),
norm(0),
precision(0),
norm_inv(0),
delta_m(0)
{
//  RationalApproxType approx_type;

//  FieldType field_type;

for(int i =0; i< MAX_RAT_DEGREE;i++){

  residue[i]=0.;
  pole[i]=0.;

  residue_inv[i]=0.;
  pole_inv[i]=0.;
}


};
CPS_END_NAMESPACE
