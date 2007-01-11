#ifndef INCLUDED_PT_QCDOC_H
#define INCLUDED_PT_QCDOC_H
#include "asq_data_types.h"
#include "pt_int.h"

static unsigned long PEC = 0xb0000000;
static unsigned long PLB = 0xb0000000;

#define TESTING
#undef TESTING
#undef CPP
//CPS_START_NAMESPACE
//void dirac_cmv_jcw_agg_cpp( int sites, long chi, long u,long in, long out);

//External function definitions
extern "C"{

  //matrix multiply for checkerboarded fields
  void pt_cmm_cpp(int sites, long u, long in, long out, long gauge_field);
  void pt_cmm_dag_cpp(int sites, long u, long in, long out, long gauge_field);

  //------------------------------------------------------------------------
  //C++ routines
#ifdef CPP
  //Matrix multiplication for full matrix fields
  void cmm_agg_cpp(gauge_agg *chi, matrix *phi, matrix *result, int counter);
  void cmv_agg_cpp( int sites, long u,long in, long out);
  #define partrans_cmm_agg(A,B,C,D) cmm_agg_cpp(A,B,C,D)
  #define partrans_cmv_agg(A,B,C,D) cmv_agg_cpp(A,B,C,D)

  //matrix vector multiply for checkerboarded fields
  void pt_cmv_cpp(int sites, ind_agg *agg, double *gauge_field, double *src, double *dest);
  void pt_cmv_dag_cpp(int sites, ind_agg *agg, double *gauge_field, double *src, double *dest);
  void pt_cmv_pad_cpp(int sites, ind_agg *agg, double *gauge_field, double *src, double *dest);
  void pt_cmv_dag_pad_cpp(int sites, ind_agg *agg, double *gauge_field, double *src, double *dest);
  #define partrans_cmv(A,B,C,D,E) pt_cmv_cpp(A,B,C,D,E)
  #define partrans_cmv_dag(A,B,C,D,E) pt_cmv_dag_cpp(A,B,C,D,E)
  #define partrans_cmv_pad(A,B,C,D,E) pt_cmv_pad_cpp(A,B,C,D,E)
  #define partrans_cmv_dag_pad(A,B,C,D,E) pt_cmv_dag_pad_cpp(A,B,C,D,E)
  //--------------------------------------------------------------------------
  //Assembly Routines
#else
  //Matrix multiplication for full matrix fields
  void pt_cmm_agg(gauge_agg *chi, matrix *phi,matrix *result, int counter);
  //void cmm_agg(gauge_agg *chi, matrix *phi,matrix *result, int counter);
  void pt_asqtad_agg( int sites, long chi, long u,long in, long out);
  void pt_asqtad_agg_s( int sites, long chi, long u,long in, long out);
  #define partrans_cmm_agg(A,B,C,D) pt_cmm_agg(A,B,C,D)
  #define partrans_cmv_agg(A,B,C,D) pt_asqtad_agg(A,0,B,C,D)

  void pt_cmv(int count, ind_agg *ind, double *gauge, double *src, double *dest);
  void pt_cmv_pad(int count, ind_agg *ind, double *gauge, double *src, double *dest);
  void pt_cmv_dag(int count, ind_agg *ind, double *gauge, double *src, double *dest);
  void pt_cmv_dag_pad(int count, ind_agg *ind, double *gauge, double *src, double *dest);
  void pt_cmv_s(int count, ind_agg *ind, float *gauge, float *src, float *dest);
  void pt_cmv_pad_s(int count, ind_agg *ind, float *gauge, float *src, float *dest);
  void pt_cmv_dag_s(int count, ind_agg *ind, float *gauge, float *src, float *dest);
  void pt_cmv_dag_pad_s(int count, ind_agg *ind, float *gauge, float *src, float *dest);
  #define partrans_cmv(A,B,C,D,E) pt_cmv(A,B,C,D,E)
  #define partrans_cmv_dag(A,B,C,D,E) pt_cmv_dag(A,B,C,D,E)
  #define partrans_cmv_pad(A,B,C,D,E) pt_cmv_pad(A,B,C,D,E)
  #define partrans_cmv_dag_pad(A,B,C,D,E) pt_cmv_dag_pad(A,B,C,D,E)
#endif
  //--------------------------------------------------------------------------

  void pt_copy(int count, ind_agg *ind, double *src, double *dest);
  void pt_copy_pad(int count, ind_agg *ind, double *src, double *dest);
  void pt_copy_s(int count, ind_agg *ind, float *src, float *dest);
  void pt_copy_pad_s(int count, ind_agg *ind, float *src, float *dest);

  void pt_copy_buffer(int n, long src, long dest, long ptable);
  // Assembler copying routines
  void copy_matrix(IFloat *res, IFloat *src, int *length, 
		   unsigned long *res_ptr, unsigned long *src_ptr);
  void copy_gauge(IFloat *res, struct gauge_agg *src, int *length,
		  unsigned long *res_ptr);
  // This is perhaps overkill but gives a couple of extra flops
  // cross_look - all input fields are lookup and sum to result
  void cross_look(IFloat *result, Float *fac, const IFloat *chi, const IFloat *phi,  
	     int counter, unsigned long *src, unsigned long *dest);
  // cross_lin - one input field is linear and sum to result
  void cross_lin(IFloat *result, Float *fac, const IFloat *chi, const IFloat *phi,  
	     int counter, unsigned long *dest, unsigned long *dest);
  // cross_over_look - all input fields are lookup and overwrite result
  void cross_over_look(IFloat *result, Float *fac, const IFloat *chi, const IFloat *phi,  
	     int counter, unsigned long *src, unsigned long *dest);
  // cross_over_lin - one input field is linear and overwrite result
  void cross_over_lin(IFloat *result, Float *fac, const IFloat *chi, const IFloat *phi,  
	     int counter, unsigned long *dest, unsigned long *dest);
  //Copies a vectors from v to u
  void copy_vector(IFloat *u, IFloat *v, int *length, unsigned long *dest, unsigned long *src);

  //---------------------------------------------------------------------------
  
  //---------------------------------------------------------------------------

  void m1m2_lookup(matrix *result, matrix *m1, matrix *m2, int length,
		   unsigned long *dest, unsigned long *dest, unsigned long *src);
  void m1m2_lookup_copy(matrix *result2, matrix *result, matrix *m1, matrix *m2, 
			int length, unsigned long *dest2,  
			unsigned long *dest, unsigned long *dest, 
			unsigned long *src);
  void m1m2_lin_copy(matrix *result2, matrix *result, matrix *m1, matrix *m2, 
		     int length, unsigned long *dest2,
		     unsigned long *dest, unsigned long *dest);
  
}
inline  void pt_cmm_agg_print(gauge_agg *chi, matrix *phi,matrix *result, int counter){
   printf("pt_cmm_agg(%p %p %p %d)\n",chi,phi,result,counter);
//    for(int i =0;i<2*counter;i++){
//      printf("%d: %d %d\n",i,chi[i].src,chi[i].dest);
//    }
   printf("pt_cmm_agg(%p %p %p %d) done \n",chi,phi,result,counter);
}

inline  void cross_over_lin_cpp(IFloat *result, Float *fac, const IFloat *chi,
 const IFloat *phi,  int counter, unsigned long *src, unsigned long *dest){
    printf("cross_over_lin(%p %0.4f %p %p %d %p %p)\n",
    result,*fac,chi,phi,counter,src,dest);
    for(int i =0;i<counter;i++){
      printf("%d: %d %d\n",i,src[i],dest[i]);
    }
    cross_over_lin(result,fac,chi,phi,counter,src,dest);
    printf("cross_over_lin(%p %0.4f %p %p %d %p %p) done\n");
}

inline  void cross_over_look_cpp(IFloat *result, Float *fac, const IFloat *chi,
 const IFloat *phi,  int counter, unsigned long *src, unsigned long *dest){
    printf("cross_over_look(%p %0.4f %p %p %d %p %p)\n",
    result,*fac,chi,phi,counter,src,dest);
    for(int i =0;i<counter;i++){
      printf("%d: %d %d\n",i,src[i],dest[i]);
    }
    cross_over_look(result,fac,chi,phi,counter,src,dest);
    printf("cross_over_look(%p %0.4f %p %p %d %p %p) done\n");
}

#ifdef ASQD_SINGLE
#define pt_asqtad_agg(A,B,C,D,E) pt_asqtad_agg_s(A,B,C,D,E)
#define pt_cmv(A,B,C,D,E) pt_cmv_s(A,B,C,D,E)
#define pt_cmv_dag(A,B,C,D,E) pt_cmv_dag_s(A,B,C,D,E)
#define pt_cmv_pad(A,B,C,D,E) pt_cmv_pad_s(A,B,C,D,E)
#define pt_cmv_dag_pad(A,B,C,D,E) pt_cmv_dag_pad_s(A,B,C,D,E)
#define pt_copy_pad(A,B,C,D) pt_copy_pad_s(A,B,C,D)
#define pt_copy(A,B,C,D) pt_copy_s(A,B,C,D)
#endif

#endif
