#ifndef DWF_INTERNAL
#define DWF_INTERNAL
extern "C" {
void vaxpy3_norm (Float *out,Float *scalep,
                 Float *InScale,
                 Float *Add,
                 int len,
                 Float *res);
 
void vaxpy3 (Float *res,Float *scalep,Float *InScale, Float *Add, int len);
}
#endif
