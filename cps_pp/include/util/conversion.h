#ifndef _CONVERSION_H
#define _CONVERSION_H
#include <comms/glb.h>

CPS_START_NAMESPACE
#if 1
#define NOINLINE_MACRO __attribute((noinline))
#else
#define NOINLINE_MACRO 
#endif

inline void moveFloattofloat (float *out, Float * in, size_t f_size)
{
  Float  sum=0.;
#pragma omp parallel for reduction(+:sum)
  for (size_t i = 0; i < f_size; i++) {
      out[i]=in[i];
//    flt = (float) in[i];
//    out[i] = flt;
      sum +=out[i]*out[i];
  }
  glb_sum(&sum);
//  VRB.Result("","moveFloattofloat()","norm=%e\n",sum);
};

inline void movefloattoFloat (Float * out, float *in, size_t f_size)
{
//  float flt;
  Float  sum=0.;
#pragma omp parallel for reduction(+:sum)
  for (size_t i = 0; i < f_size; i++) {
      out[i]=in[i];
//    flt = in[i];
//    out[i] = (Float) flt;
      sum +=out[i]*out[i];
  }
  glb_sum(&sum);
//  VRB.Result("","moveFloattofloat()","norm=%e\n",sum);
};
CPS_END_NAMESPACE

#endif
