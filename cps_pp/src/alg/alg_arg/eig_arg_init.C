#include <alg/eig_arg.h>

CPS_START_NAMESPACE

EigArg::EigArg(){
  Mass.Mass_len=0;
  Mass.Mass_val=NULL;
}

EigArg::~EigArg(){
//  if(Mass.Mass_val)
//    delete[] Mass.Mass_val;
}

void EigArg::resize ( u_int n_floats){
  Mass.Mass_len=n_floats;
//  if(Mass.Mass_val)
//    delete[] Mass.Mass_val;
  if(n_floats>0)
    Mass.Mass_val= new Float[n_floats];
  else
    Mass.Mass_val= NULL;
}
CPS_END_NAMESPACE
