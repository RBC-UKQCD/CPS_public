#line 4283 "dwf.nw"
#define Nc  3      /* Number of colors */
#define DIM 4      /* number of dimensions */
#define Fd  4      /* Fermion representation dimension */
#line 4469 "dwf.nw"
#define Vs 2
typedef double REAL;
/* typedef double __attribute__((mode(V2DF),aligned(16))) vReal; */
#include <xmmintrin.h>
typedef __m128d vReal;
#line 4295 "dwf.nw"
typedef struct {
    float re, im;
} scalar_complex;

typedef struct {
   vReal re, im;
} vector_complex;
#line 4478 "dwf.nw"
static inline vReal
vmk_1(double a)                                        /* return (a a ... a) */
{
     vReal v;
     REAL *w = (REAL *)&v;
     w[0] = w[1] = a;
     return v;
}
#line 4490 "dwf.nw"
static inline vReal
vmk_n1(double a, double b)                             /* return (a ... a b) */
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a; r[1] = b;
  return v;
}
#line 4502 "dwf.nw"
static inline vReal
vmk_1n(double a, double b)                             /* return (a b ... b) */
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a; r[1] = b;
  return v;
}
#line 4515 "dwf.nw"
static inline vReal
vmk_fn(double a, double b)                  /* return (a a*b ... a*b^(Vs-1)) */
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a; r[1] = a*b;
  return v;
}
#line 4526 "dwf.nw"
static inline vReal
vmk_bn(double a, double b)                  /* return (a^(Vs-1)*b ... a*b b) */
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a*b; r[1] = b;
  return v;
}
#line 4538 "dwf.nw"
static inline double
vsum(vReal v)                                         /* return sum(i, [i]v) */
{
  REAL *r = (REAL *)&v;
  return r[0] + r[1];
}
#line 4548 "dwf.nw"
static inline vReal
vput_0(vReal a, double b)                     /* return (b [1]a ... [Vs-1]a) */
{
   REAL *v = (REAL *)&a;
   v[0] = b;
   return a;
}
#line 4559 "dwf.nw"
static inline vReal
vput_n(vReal a, double b)                     /* return ([0]a ... [Vs-2]a b) */
{
   REAL *v = (REAL *)&a;
   v[1] = b;
   return a;
}
#line 4571 "dwf.nw"
static inline vReal
shift_up1(vReal a, vReal b)                /* return ([1]a ... [Vs-1]a [0]b) */
{
   vReal r;
   REAL *x = (REAL *)&a;
   REAL *y = (REAL *)&b;
   REAL *z = (REAL *)&r;
   z[0] = x[1];
   z[1] = y[0];

   return r;
}
#line 4586 "dwf.nw"
static inline vReal
shift_upN(vReal a, vReal b)             /* return ([Vs-1]a [0]b ... [Vs-2]b) */
{
   vReal r;
   REAL *x = (REAL *)&a;
   REAL *y = (REAL *)&b;
   REAL *z = (REAL *)&r;
   z[0] = x[1];
   z[1] = y[0];

   return r;
}
#line 4456 "dwf.nw"
#include "dwf-ssed.h"
#define MACHINE "sse double"
#define L3(n) MIT_ssed_##n
#define PAD16(size) (15+(size))
#define ALIGN16(addr) ((void *)(~15 & (15 + (size_t)(addr))))
#line 3028 "dwf.nw"
#define BLOCKOF_YA(j,t,c,ri) BLOCKOF2_YA(j,t,c,ri)
#define BLOCKOF_YB(j,t,c,ri) BLOCKOF2_YB(j,t,c,ri)
#line 4462 "dwf.nw"
#include "dwf.c"
