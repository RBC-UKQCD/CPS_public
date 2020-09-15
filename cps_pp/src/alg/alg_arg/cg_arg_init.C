/*
 * put your initializations here
 */

#include "alg/mobius_arg.h"
#include "alg/cg_arg.h"
CPS_START_NAMESPACE
CgArg::CgArg()
{
  Inverter=CG;
}

MobiusArg::MobiusArg()
{
   mobius_b_coeff = 1.0;
   mobius_c_coeff = 0.0;
   M5 = 1.0;
};

CPS_END_NAMESPACE
