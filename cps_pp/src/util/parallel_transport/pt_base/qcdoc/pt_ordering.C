/*! \file
  \brief  Definition of parallel transport definitions for QCDOC.
  
  $Id: pt_ordering.C,v 1.1 2006-07-05 18:13:49 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2006-07-05 18:13:49 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qcdoc/pt_ordering.C,v 1.1 2006-07-05 18:13:49 chulwoo Exp $
//  $Id: pt_ordering.C,v 1.1 2006-07-05 18:13:49 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_ordering.C,v $
//  $Revision: 1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qcdoc/pt_ordering.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
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


//dest=src
void PT::cpy (IFloat *dest, IFloat *src){
  for(int i=0;i<18;i++)
    dest[i]=src[i];
}

//dest=src.Dagger()
void PT::dag_cpy (IFloat *dest, IFloat *src){
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++){
      dest[2*(3*i+j)]=src[2*(3*j+i)];
      dest[2*(3*i+j)+1]=-src[2*(3*j+i)+1];
    }
}

//Returns lexical index associated with coordinate x[4]
//where the 0th coordinate runs fastest, 3rd coordinate runs slowest
int PT::lex_xyzt(int *x){
//  printf("lex_xyzt(%d %d %d %d)\n",x[0],x[1],x[2],x[3]);
  int result = x[0] + size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3]));
  return result;
}

//Returns checkerboard index associated with coordinate x[4]
int PT::lex_xyzt_cb_o(int *x){
//  printf("lex_xyzt_cb_o(%d %d %d %d)\n",x[0],x[1],x[2],x[3]);
  int result = x[0] + size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3]));
  if ( (x[0]+x[1]+x[2]+x[3]+evenodd)%2 == 0) result = result/2+vol/2;
  else result = result/2;
  return result;
}

//Returns checkerboard index associated with coordinate x[4]
int PT::lex_xyzt_cb_e(int *x){
//  printf("lex_xyzt_cb_o(%d %d %d %d)\n",x[0],x[1],x[2],x[3]);
  int result = x[0] + size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3]));
  if ( (x[0]+x[1]+x[2]+x[3]+evenodd)%2 == 1) result = result/2+vol/2;
  else result = result/2;
  return result;
}


//---------------------------------------------------------------------------
//Returns index for fields in the STAG storage order on lattice
//sites of a given parity
int PT:: lex_txyz_cb(int *x)
{
//  printf("lex_txyz_cb(%d %d %d %d)\n",x[0],x[1],x[2],x[3]);
  int result = x[3]+size[3]*(x[0]+size[0]*(x[1]+size[1]*x[2]));
  return result/2;
}
//---------------------------------------------------------------------------

//Returns index associated with x[4] for txyz ordering
int PT:: lex_txyz(int *x){
  return  (x[3] + size[3]*(x[0]+size[0]*(x[1]+size[1]*x[2])))/2 ;
}

//Returns first index associated with the surface x[3] = 0, on
//a checkerboarded lattice
int PT:: LexSurface(int *x){
  return  (x[0]+size[0]*(x[1]+size[1]*x[2]))/2 ;
}

//Returns index associated with gauge link in the mu direction 
//and coordinate x
int PT::lex_g_xyzt(int *x, int mu){
  int temp =  lex_xyzt(x);
  return (temp*NDIM + mu);
}

//---------------------------------------------------------------------------
//Returns block ordering for the gauge fields, where all directions
//are stored in one block and sites ordered txyz
int PT::lex_g_txyz(int *x, int mu){
  int temp = mu*vol+x[3] + size[3]*(x[0]+size[0]*(x[1]+size[1]*x[2]));
  return temp;
}

//Returns block ordering for the gauge fields, where all directions
//are stored in one block with sites checkerboarded txyz
int PT::lex_g_txyz_cb(int *x, int mu){
  int result = (x[3]+size[3]*(x[0]+size[0]*(x[1]+size[1]*x[2])))/2;
  return mu*vol + result + ((x[0]+x[1]+x[2]+x[3])%2)*vol/2;
}
//---------------------------------------------------------------------------

//Returns index associated with gauge link in the mu direction and 
//coordinate x for checkerboarded storage
int PT:: lex_g_xyzt_cb_o(int *x, int mu){
  int temp =  lex_xyzt_cb_o(x);
  return (mu*vol+temp);
}
int PT:: lex_g_xyzt_cb_e(int *x, int mu){
  int temp =  lex_xyzt_cb_e(x);
  return (mu*vol+temp);
}
//CPS_END_NAMESPACE
