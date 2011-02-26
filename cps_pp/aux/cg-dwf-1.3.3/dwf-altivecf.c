#line 4283 "dwf.nw"
#define Nc  3      /* Number of colors */
#define DIM 4      /* number of dimensions */
#define Fd  4      /* Fermion representation dimension */
#line 4617 "dwf.nw"
#include <altivec.h>
typedef float REAL;
typedef vector float vReal;
#define Vs 4
#line 4295 "dwf.nw"
typedef struct {
    float re, im;
} scalar_complex;

typedef struct {
   vReal re, im;
} vector_complex;
#line 4628 "dwf.nw"
static inline vReal
vmk_1(double a)                                        /* return (a a ... a) */
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a; r[1] = a; r[2] = a; r[3] = a;
  return v;
}
#line 4640 "dwf.nw"
static inline vReal
vmk_n1(double a, double b)                             /* return (a ... a b) */
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a; r[1] = a; r[2] = a; r[3] = b;
  return v;
}
#line 4652 "dwf.nw"
static inline vReal
vmk_1n(double a, double b)                             /* return (a b ... b) */
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a; r[1] = b; r[2] = b; r[3] = b;
  return v;
}
#line 4665 "dwf.nw"
static inline vReal
vmk_fn(double a, double b)                  /* return (a a*b ... a*b^(Vs-1)) */
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a; r[1] = a*b; r[2] = a*b*b; r[3] = a*b*b*b;
  return v;
}
#line 4676 "dwf.nw"
static inline vReal
vmk_bn(double a, double b)                  /* return (a^(Vs-1)*b ... a*b b) */
{
  vReal v;
  REAL *r = (REAL *)&v;
  r[0] = a*a*a*b; r[1] = a*a*b; r[2] = a*b; r[3] = b;
  return v;
}
#line 4688 "dwf.nw"
static inline double
vsum(vReal v)                                         /* return sum(i, [i]v) */
{
  REAL *r = (REAL *)&v;
  return r[0] + r[1] + r[2] + r[3];
}
#line 4698 "dwf.nw"
static inline vReal
vput_0(vReal a, double b)                     /* return (b [1]a ... [Vs-1]a) */
{
   REAL *v = (REAL *)&a;
   v[0] = b;
   return a;
}
#line 4709 "dwf.nw"
static inline vReal
vput_n(vReal a, double b)                     /* return ([0]a ... [Vs-2]a b) */
{
   REAL *v = (REAL *)&a;
   v[3] = b;
   return a;
}
#line 4720 "dwf.nw"
static inline vReal
shift_up1(vReal a, vReal b)                /* return ([1]a ... [Vs-1]a [0]b) */
{
   return vec_sld(a, b, 4);
}
#line 4728 "dwf.nw"
static inline vReal
shift_upN(vReal a, vReal b)             /* return ([Vs-1]a [0]b ... [Vs-2]b) */
{
   return vec_sld(a, b, 12);
}
#line 4605 "dwf.nw"
#include "dwf-altivecf.h"
#define MACHINE "altivec float"
#define L3(n) MIT_altivecf_##n
#define PAD16(size) (15+(size))
#define ALIGN16(addr) ((void *)(~15 & (15 + (size_t)(addr))))
#line 3008 "dwf.nw"
#define BLOCKOF_YA(j,t,c,ri) BLOCKOF4_YA(j,t,c,ri)
#define BLOCKOF_YB(j,t,c,ri) BLOCKOF4_YB(j,t,c,ri)
#line 4611 "dwf.nw"
#include "dwf.c"
