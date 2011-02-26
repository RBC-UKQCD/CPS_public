#line 115 "interface.nw"
#ifndef _DWF_SSEF_H
#define _DWF_SSEF_H
#define L3(n) MIT_ssef_##n
#line 153 "interface.nw"
#include <stdlib.h>
#line 162 "interface.nw"
#if defined (__cplusplus)
extern "C" {
#endif
#line 176 "interface.nw"
typedef struct L3(DWF_Fermion)   L3(DWF_Fermion);
typedef struct L3(DWF_Gauge)     L3(DWF_Gauge);
#line 191 "interface.nw"
typedef double (*L3(DWF_gauge_reader))(const void     *OuterGauge,
                                       void           *env,
                                       const int       lattice_addr[4],
                                       int             dim,
                                       int             a,
                                       int             b,
                                       int             re_im);
#line 215 "interface.nw"
typedef double (*L3(DWF_fermion_reader))(const void     *OuterFermion,
                                         void           *env,
                                         const int       lattice_addr[5],
                                         int             color,
                                         int             dirac,
                                         int             re_im);
#line 229 "interface.nw"
typedef void (*L3(DWF_fermion_writer))(void              *OuterFemrion,
                                       void              *env,
                                       const int          lattice_addr[5],
                                       int                color,
                		       int                dirac,
                                       int                re_im,
                                       double             value);
#line 245 "interface.nw"
int L3(DWF_init)(const int         lattice[5],
                 void           *(*allocator)(size_t size),
                 void            (*deallocator)(void *));
#line 272 "interface.nw"
void L3(DWF_fini)(void);
#line 292 "interface.nw"
L3(DWF_Gauge) *L3(DWF_load_gauge)(const void            *OuterGauge_U,
                                  const void            *OuterGauge_V,
                                  void                  *env,
                                  L3(DWF_gauge_reader)   reader);
#line 303 "interface.nw"
void L3(DWF_delete_gauge)(L3(DWF_Gauge) *);
#line 312 "interface.nw"
L3(DWF_Fermion) *L3(DWF_load_fermion)(const void              *OuterFermion,
                                      void                    *env,
                                      L3(DWF_fermion_reader)   reader);
#line 321 "interface.nw"
L3(DWF_Fermion) *L3(DWF_allocate_fermion)(void);
#line 326 "interface.nw"
void L3(DWF_delete_fermion)(L3(DWF_Fermion) *);
#line 333 "interface.nw"
void L3(DWF_save_fermion)(void                    *OuterFermion,
                          void                    *env,
                          L3(DWF_fermion_writer)   writer,
			  L3(DWF_Fermion)         *CGfermion);
#line 347 "interface.nw"
int L3(DWF_cg_solver)(L3(DWF_Fermion)       *result,
                      double                *out_eps,
                      int                   *out_iter,
                      const L3(DWF_Gauge)   *gauge,
                      double                 M_0,
                      double                 m_f,
                      const L3(DWF_Fermion) *guess,
                      const L3(DWF_Fermion) *rhs,
                      double                 eps,
                      int                    min_iter,
                      int                    max_iter);
#line 373 "interface.nw"
void L3(DWF_Dirac_Operator)(L3(DWF_Fermion)        *chi,
                            const L3(DWF_Gauge)    *gauge,
                            double                  M_0,
                            double                  m_f,
                            const L3(DWF_Fermion)  *psi);
#line 383 "interface.nw"
void L3(DWF_Dirac_Operator_conjugate)(L3(DWF_Fermion)        *chi,
                                      const L3(DWF_Gauge)    *gauge,
                                      double                  M_0,
                                      double                  m_f,
                                      const L3(DWF_Fermion)  *psi);
#line 394 "interface.nw"
void L3(DWF_Add_Fermion)(L3(DWF_Fermion)         *psi,
	                 const L3(DWF_Fermion)   *phi,
			 double                   a,
	                 const L3(DWF_Fermion)   *eta);
#line 403 "interface.nw"
void L3(DWF_Fermion_Dot_Product)(double                  *v_re,
	                         double                  *v_im,
                                 const L3(DWF_Fermion)   *psi,
                                 const L3(DWF_Fermion)   *phi);
#line 167 "interface.nw"
#if defined (__cplusplus)
}
#endif
#line 119 "interface.nw"
#undef L3
#endif
