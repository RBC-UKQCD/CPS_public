#include <alg/array_arg.h>

CPS_START_NAMESPACE

FloatArray::FloatArray(){
  Floats.Floats_len=0;
  Floats.Floats_val=NULL;
}

#if 0
FloatArray::~FloatArray(){
  if(Floats.Floats_val)
    delete[] Floats.Floats_val;
}
#endif

void FloatArray::resize ( int n_floats){
  Floats.Floats_len=n_floats;
  if(Floats.Floats_val)
    delete[] Floats.Floats_val;
  if(n_floats>0)
    Floats.Floats_val= new Float[n_floats];
  else
    Floats.Floats_val= NULL;
}
CPS_END_NAMESPACE
