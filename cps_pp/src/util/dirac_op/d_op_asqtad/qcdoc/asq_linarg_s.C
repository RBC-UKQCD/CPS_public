#include "asq_data_types.h"
#include "asqtad_int.h"

extern "C" 
void asq_vaxpy3_cpp(Float *res,Float *scale, Float *mult, Float *add,int
nvec){
  Float *dest = res;

  for(int i = 0;i<nvec;i++){
    *dest++ = *scale * (*mult++) + (*add++);
    *dest++ = *scale * (*mult++) + (*add++);
    *dest++ = *scale * (*mult++) + (*add++);
    *dest++ = *scale * (*mult++) + (*add++);
    *dest++ = *scale * (*mult++) + (*add++);
    *dest++ = *scale * (*mult++) + (*add++);
  }
}

extern "C" 
void asq_vaxpy3_norm_cpp(Float *res,Float *scale, Float *mult, Float *add,int
nvec, Float *norm){
  Float *dest = res;
  *norm = 0.;
  register Float temp;

  for(int i = 0;i<nvec;i++){
    temp = *scale * (*mult++) + (*add++); *norm += temp*temp; (*dest++)=temp;
    temp = *scale * (*mult++) + (*add++); *norm += temp*temp; (*dest++)=temp;
    temp = *scale * (*mult++) + (*add++); *norm += temp*temp; (*dest++)=temp;
    temp = *scale * (*mult++) + (*add++); *norm += temp*temp; (*dest++)=temp;
    temp = *scale * (*mult++) + (*add++); *norm += temp*temp; (*dest++)=temp;
    temp = *scale * (*mult++) + (*add++); *norm += temp*temp; (*dest++)=temp;
  }
}

extern "C" 
void asq_vaxmy_cpp(Float *scale,Float *mult, Float *sub,
  int nvec){
  register Float temp;

  for(int i = 0;i<nvec;i++){
    *sub = *scale * (*mult++) - (*sub); sub++;
    *sub = *scale * (*mult++) - (*sub); sub++;
    *sub = *scale * (*mult++) - (*sub); sub++;
    *sub = *scale * (*mult++) - (*sub); sub++;
    *sub = *scale * (*mult++) - (*sub); sub++;
    *sub = *scale * (*mult++) - (*sub); sub++;
  }
}

extern "C" 
void asq_vaxmy_vxdot_cpp(Float *scale,Float *mult, Float *sub,
  int nvec, Float *norm){
  Float *dest = sub;
  *norm = 0.;
  register Float temp;

  for(int i = 0;i<nvec;i++){
    temp = *scale * (*mult) - (*sub++); *norm += temp*(*mult++); (*dest++)=temp;
    temp = *scale * (*mult) - (*sub++); *norm += temp*(*mult++); (*dest++)=temp;
    temp = *scale * (*mult) - (*sub++); *norm += temp*(*mult++); (*dest++)=temp;
    temp = *scale * (*mult) - (*sub++); *norm += temp*(*mult++); (*dest++)=temp;
    temp = *scale * (*mult) - (*sub++); *norm += temp*(*mult++); (*dest++)=temp;
    temp = *scale * (*mult) - (*sub++); *norm += temp*(*mult++); (*dest++)=temp;
  }
}

extern "C"
void copy_buffer_cpp(int n, long src, long dest, long ptable){
  int *offset = (int *)ptable;
  Float *d_p = (Float *)dest;
  for(int i = 0;i<n;i++){
    Float *s_p = (Float *)src;
    s_p +=offset[i];
//    printf("s_p=%p d_p=%p offset=%x\n",s_p,d_p,offset[i]);
    *d_p++ = *s_p++;
    *d_p++ = *s_p++;
    *d_p++ = *s_p++;
    *d_p++ = *s_p++;
    *d_p++ = *s_p++;
    *d_p++ = *s_p++;
  }
}

