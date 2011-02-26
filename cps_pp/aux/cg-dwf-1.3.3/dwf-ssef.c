#line 4283 "dwf.nw"
#define Nc  3      /* Number of colors */
#define DIM 4      /* number of dimensions */
#define Fd  4      /* Fermion representation dimension */
#line 4322 "dwf.nw"
#define Vs 4
typedef float REAL;
/*typedef REAL vReal __attribute__((mode(V4SF),aligned(16))); */
#include <xmmintrin.h>
typedef __m128 vReal;
#line 4295 "dwf.nw"
typedef struct {
    float re, im;
} scalar_complex;

typedef struct {
   vReal re, im;
} vector_complex;
#line 4331 "dwf.nw"
static inline vReal
vmk_1(double a)                                        /* return (a a ... a) */
{
#if 1
     float b = a;
     vReal v = _mm_load_ss((float *)&b);
     asm("shufps\t$0,%0,%0" : "+x" (v));
#else
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a; r[1] = a; r[2] = a; r[3] = a;
#endif
  return v;
}
#line 4343 "dwf.nw"
static inline vReal
vmk_n1(double a, double b)                             /* return (a ... a b) */
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a; r[1] = a; r[2] = a; r[3] = b;
  return v;
}
#line 4355 "dwf.nw"
static inline vReal
vmk_1n(double a, double b)                             /* return (a b ... b) */
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a; r[1] = b; r[2] = b; r[3] = b;
  return v;
}
#line 4368 "dwf.nw"
static inline vReal
vmk_fn(double a, double b)                  /* return (a a*b ... a*b^(Vs-1)) */
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a; r[1] = a*b; r[2] = a*b*b; r[3] = a*b*b*b;
  return v;
}
#line 4379 "dwf.nw"
static inline vReal
vmk_bn(double a, double b)                  /* return (a^(Vs-1)*b ... a*b b) */
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a*a*a*b; r[1] = a*a*b; r[2] = a*b; r[3] = b;
  return v;
}
#line 4391 "dwf.nw"
static inline double
vsum(vReal v)                                         /* return sum(i, [i]v) */
{
  REAL *r = (REAL *)&v;
  return r[0] + r[1] + r[2] + r[3];
}
#line 4401 "dwf.nw"
static inline vReal
vput_0(vReal a, double b)                     /* return (b [1]a ... [Vs-1]a) */
{
   REAL *v = (REAL *)&a;
   v[0] = b;
   return a;
}
#line 4412 "dwf.nw"
static inline vReal
vput_n(vReal a, double b)                     /* return ([0]a ... [Vs-2]a b) */
{
   REAL *v = (REAL *)&a;
   v[3] = b;
   return a;
}
#line 4424 "dwf.nw"
static inline vReal
shift_up1(vReal a, vReal b)                /* return ([1]a ... [Vs-1]a [0]b) */
{
   vReal z;
   REAL *X = (REAL *)&a;
   REAL *Y = (REAL *)&b;
   REAL *Z = (REAL *)&z;

   Z[0] = X[1]; Z[1] = X[2]; Z[2] = X[3]; Z[3] = Y[0];
   return z;
}
#line 4438 "dwf.nw"
static inline vReal
shift_upN(vReal a, vReal b)             /* return ([Vs-1]a [0]b ... [Vs-2]b) */
{
   vReal z;
   REAL *X = (REAL *)&a;
   REAL *Y = (REAL *)&b;
   REAL *Z = (REAL *)&z;

   Z[0] = X[3]; Z[1] = Y[0]; Z[2] = Y[1]; Z[3] = Y[2];
   return z;
}
#line 4309 "dwf.nw"
#include "dwf-ssef.h"
#define MACHINE "sse float"
#define L3(n) MIT_ssef_##n
#define PAD16(size) (15+(size))
#define ALIGN16(addr) ((void *)(~15 & (15 + (size_t)(addr))))
#line 3008 "dwf.nw"
#define BLOCKOF_YA(j,t,c,ri) BLOCKOF4_YA(j,t,c,ri)
#define BLOCKOF_YB(j,t,c,ri) BLOCKOF4_YB(j,t,c,ri)
#line 4315 "dwf.nw"
#include "dwf.c"
