#define INCLUDED_ALG_INT_FACT_H

#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief  Definitions of the AlgIntFactory class and derived classes.

*/
//------------------------------------------------------------------


#ifndef INCLUDED_REMEZ_ARG_H
#define INCLUDED_REMEZ_ARG_H

CPS_END_NAMESPACE
#include <alg/enum.h>

CPS_START_NAMESPACE

class RemezArg{

 public:
  int valid_approx;

  int degree;

  Float lambda_low;
  Float lambda_high;

  int power_num;
  int power_den;

  FieldType field_type;

  Float error;
  Float norm;
  Float residue[MAX_RAT_DEGREE];
  Float pole[MAX_RAT_DEGREE];

  Float norm_inv;
  Float residue_inv[MAX_RAT_DEGREE];
  Float pole_inv[MAX_RAT_DEGREE];

  long precision;
};

#endif

CPS_END_NAMESPACE
