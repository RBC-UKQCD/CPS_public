#ifdef USE_QMP
#include <stdio.h>
#include <math.h>
#include "asq_data_types.h"
#include "pt_int.h"
#undef USE_SSE
#ifdef USE_SSE
#include "inline_sse.h"
#endif


//----------------------------------------------------------------------------
//C++ implementation of vector copy for ParTransAsqtad.shift_field_vec()
//
//Parameters
//
//u - Pointer to shifted field (destination)
//v - Pointer to source field
//length - number of vectors to copy
//dest - pointer to struct holding destination index
//src - point to struct holding source index

#if 1
void inline PRINT(char *format,...) {}
#else
#define PRINT printf
#endif

extern "C" {
void copy_vector(IFloat *u, IFloat *v, int *length, unsigned long *dest, unsigned long *src)
{
 int total_copies = *(length);
 for(int i = 0; i < total_copies; i++)
   for(int j = 0; j < 6; j++)
     *(u + dest[i]*6+j) = *(v + src[i]*6+j);
}


//-----------------------------------------------------------------------------
//C++ implementation of matrix multiply for checkerboarded field
//
//Parameters
//
//sites - Total number of matrix multiplies
//u - Pointer to structs that hold source, destination indexes for parallel 
//    transport, also holds the index for the appropriate SU(3) gauge link
//a - Matrix field that is to be transported
//tmpfrm - Transported matrix field
//gauge_field - pointer to the gauge_field

#if TARGET != BGL && TARGET != BGP
  void copy_gauge(IFloat *res, struct gauge_agg *src, int *length,
		  unsigned long *res_ptr)
{
  IFloat *fp1, *uu;
  int s,c,d;
  struct gauge_agg *agg = (gauge_agg *)src;
  int sites = *length;
  unsigned long tmpfrm = (unsigned long)res;
  for(s=0;s<sites;s++)
    {
      //Assignment of source,destination,gauge field index
//      fp1 = (IFloat *)(tmpfrm + 3*agg[s].dest);
      fp1 = (IFloat *)(res + res_ptr[s]*18);
//      fp0 = (IFloat *)((long)a + 3*agg[s].src);
      uu = &(agg[s].mat[0]);
    
      //c indexes the row of the matrix
      //d indexes the column of the matrix
      for(d = 0; d<3; d++)
	for(c=0; c<3; c++)
	  {
	    //Re part
	    *(fp1+6*c+2*d) = *(uu+6*c+2*d);
	    //Im part 
	    *(fp1+6*c+1+2*d) = *(uu+6*c+1+2*d);
	  }
    }
}

  // Assembler copying routines
void copy_matrix(IFloat *res, IFloat *src, int *length, 
		   unsigned long *res_ptr, unsigned long *src_ptr)
{
  IFloat *fp0, *fp1;
  int s,c,d;
  struct gauge_agg_cb *agg = (gauge_agg_cb *)src;
  int sites = *length;
  unsigned long tmpfrm = (unsigned long)res;
  for(s=0;s<sites;s++)
    {
      //Assignment of source,destination,gauge field index
//      fp1 = (IFloat *)(tmpfrm + 3*agg[s].dest);
      fp1 = (IFloat *)(res + res_ptr[s]*18);
      fp0 = (IFloat *)(src + src_ptr[s]*18);
    
      //c indexes the row of the matrix
      //d indexes the column of the matrix
      for(d = 0; d<3; d++)
	for(c=0; c<3; c++)
	  {
	    //Re part
	    *(fp1+6*c+2*d) = *(fp0+6*c+2*d);
	    //Im part 
	    *(fp1+6*c+1+2*d) = *(fp0+6*c+1+2*d);
	  }
    }
}
#endif
}

extern "C"
void pt_cmm_cpp(int sites, long u, long a, long tmpfrm, long gauge_field)
{
  IFloat *fp0, *fp1, *uu;
  int s,c,d;
  struct gauge_agg_cb *agg = (gauge_agg_cb *)u;

  for(s=0;s<sites;s++)
    {
      //Assignment of source,destination,gauge field index
      fp1 = (IFloat *)(tmpfrm + 3*agg[s].dest);
      fp0 = (IFloat *)((long)a + 3*agg[s].src);
      uu = (IFloat *)(gauge_field + agg[s].gauge);
    
      //c indexes the row of the matrix
      //d indexes the column of the matrix
      for(d = 0; d<3; d++)
	for(c=0; c<3; c++)
	  {
	    //Re part
	    *(fp1+6*c+2*d) = *(fp0+2*d) * *(uu+6*c) - *(fp0+2*d+1) * *(uu+6*c+1) + 
	      *(fp0+2*d+6) * *(uu+6*c+2) - *(fp0+2*d+7)* *(uu+6*c+3) + 
	      *(fp0+2*d+12) * *(uu+6*c+4) - *(fp0+2*d+13) * *(uu+6*c+5);
	    //Im part 
	    *(fp1+6*c+1+2*d) = *(fp0+2*d) * *(uu+6*c+1) + *(fp0+1+2*d) * *(uu+6*c) + 
	      *(fp0+2*d+6) * *(uu+6*c+3) + *(fp0+2*d+7)* *(uu+6*c+2) + 
	      *(fp0+2*d+12) * *(uu+6*c+5) + *(fp0+2*d+13) * *(uu+6*c+4);
	  }
    }
}

extern "C"
void pt_cmm_dag_cpp(int sites, long u, long a, long tmpfrm, long gauge_field)
{
  IFloat *fp0, *fp1, *uu;
  int s,c,d;
  struct gauge_agg_cb *agg = (gauge_agg_cb *)u;

  for(s=0;s<sites;s++)
    {
      //Assignment of source,destination,gauge field index
      fp1 = (IFloat *)(tmpfrm + 3*agg[s].dest);
      fp0 = (IFloat *)((long)a + 3*agg[s].src);
      uu = (IFloat *)(gauge_field + agg[s].gauge);
    
      //c indexes the row of the matrix
      //d indexes the column of the matrix
      for(d = 0; d<3; d++)
	for(c=0; c<3; c++)
	  {
	    //Re part
	    *(fp1+6*c+2*d) = *(fp0+2*d) * *(uu+2*c) + *(fp0+2*d+1) * *(uu+2*c+1) +
	      *(fp0+2*d+6) * *(uu+2*c+6) + *(fp0+2*d+7) * *(uu+2*c+7) +
	      *(fp0+2*d+12) * *(uu+2*c+12) + *(fp0+2*d+13) * *(uu+2*c+13);
	    //Im part
	    *(fp1+6*c+2*d+1) = *(fp0+2*d+1) * *(uu+2*c) - *(fp0+2*d) * *(uu+2*c+1) +
	      *(fp0+2*d+7) * *(uu+2*c+6) - *(fp0+2*d+6) * *(uu+2*c+7) +
	      *(fp0+2*d+13) * *(uu+2*c+12) - *(fp0+2*d+12) * *(uu+2*c+13);
	  }
    }
}


//----------------------------------------------------------------------------

extern "C"
void cmm_agg_cpp(gauge_agg *u, matrix *a, matrix *tmpfrm, int sites)
{
  IFloat *fp0, *fp1, *uu;
  int s, c,i;
  struct gauge_agg *agg = u;
  const int SITE_LEN = 72;
  const int MATRIX_SIZE = 18;
  int *ch;

  for (s = 0; s< sites; s++)
    {
      //Assignment of source, destination indexes
      fp1 = (IFloat *)((long long)tmpfrm + 3*agg[s].dest );
      fp0 = (IFloat *)((long long)a + 3*agg[s].src);
      uu = &(agg[s].mat[0]);
#if 0
      printf("dest=%d src=%d\n", agg[s].dest, agg[s].src);

      for(i=0;i<18;i++){
      printf("%d %0.4e,",i,uu[i]);
      if(i % 6 == 5 ) printf("\n");
      }
#endif
      for (int d=0; d<3; d++)
      for (c=0; c<3; c++){
	//Re part
	*(fp1+6*c+2*d) = *(fp0+2*d) * *(uu+6*c) - *(fp0+2*d+1) * *(uu+6*c+1) + 
	  *(fp0+2*d+6) * *(uu+6*c+2) - *(fp0+2*d+7)* *(uu+6*c+3) + 
	  *(fp0+2*d+12) * *(uu+6*c+4) - *(fp0+2*d+13) * *(uu+6*c+5);
	//Im part 
	*(fp1+6*c+1+2*d) = *(fp0+2*d) * *(uu+6*c+1) + *(fp0+1+2*d) * *(uu+6*c) + 
	  *(fp0+2*d+6) * *(uu+6*c+3) + *(fp0+2*d+7)* *(uu+6*c+2) + 
	  *(fp0+2*d+12) * *(uu+6*c+5) + *(fp0+2*d+13) * *(uu+6*c+4);
//	      printf("%e %e ",*(fp1+2*c),*(fp1+2*c+1));
      }
//                 printf("\n");
    }
}

extern "C"
void pt_copy_cpp(int sites, ind_agg *u, double *a, double *tmpfrm)
{
  IFloat *fp0, *fp1, *uu, *tmp;
  int s,c,d;
  struct gauge_agg_cb *agg = u;
  PRINT("sites=%d u=%p src=%p dest=%p\n",sites,u,a,tmpfrm);

  for(s=0;s<sites;s++)
    {
      //Assignment of source, destination, gauge field indexes
      fp1 = (IFloat *)((long)tmpfrm + agg[s].dest);
      fp0 = (IFloat *)((long)a + agg[s].src);

      //d indexes the row of the gauge link matrix
      for(d = 0; d<3; d++)
	{
	  tmp = uu+6*d;
	  //Re part
	  *(fp1+2*d) = *(fp0+2*d); 
	  //Im part 
	  *(fp1+1+2*d) = *(fp0+1+2*d); 
	}
    }
}

extern "C"
void pt_copy_buffer_cpp(int n, double * src, double *dest, int *ptable)
{
	for(int j = 0; j < n;j++)
	for(int k = 0; k < 6;k++)
		*(dest+j*6+k) = *(src + *(ptable+j)+ k);
}

extern "C"
void pt_copy_pad_cpp(int sites, ind_agg *u, double *a, double *tmpfrm)
{ pt_copy_cpp(sites,u,a,tmpfrm);}

extern "C"
void pt_cmv_cpp(int sites, ind_agg *u, double *gauge_field, double *a, double *tmpfrm)
{
  int s,c,d;
  struct gauge_agg_cb *agg = u;
  //printf("index=%p src=%p dest=%p gauge=%p\n",u,a,tmpfrm,gauge_field);
  PRINT("sites=%d\n",sites);

  for(s=0;s<sites;s++)
    {
#ifdef USE_SSE
  IFloat *tmp;
  QLA_D3_ColorVector *fp0, *fp1;
  QLA_D3_ColorMatrix *uu ;
//  printf("index=%p src=%p dest=%p gauge=%p\n",u,a,tmpfrm,gauge_field);
      fp1 = (QLA_D3_ColorVector *)((unsigned long)tmpfrm + agg[s].dest);
      fp0 = (QLA_D3_ColorVector *)((unsigned long)a + agg[s].src);
      uu = (QLA_D3_ColorMatrix *)((unsigned long)gauge_field + agg[s].gauge);
      QLA_D3_V_eq_M_times_V(fp1,uu,fp0)
#else
  IFloat *fp0, *fp1, *uu, *tmp;
      //Assignment of source, destination, gauge field indexes
      fp1 = (IFloat *)((long)tmpfrm + agg[s].dest);
      fp0 = (IFloat *)((long)a + agg[s].src);
      uu = (IFloat *)((long)gauge_field + agg[s].gauge);
      if(fabs(*fp0)>1e-6){
        PRINT("a[%d]=%e\n",agg[s].src/48,*fp0);
        PRINT("uu[%d]=%e\n",agg[s].gauge/144,*uu);
      }

      //d indexes the row of the gauge link matrix
      for(d = 0; d<3; d++)
	{
	  tmp = uu+6*d;
	  //Re part
	  *(fp1+2*d) = *(fp0) * *(tmp) - *(fp0+1) * *(tmp+1) + 
	    *(fp0+2) * *(tmp+2) - *(fp0+3)* *(tmp+3) + 
	    *(fp0+4) * *(tmp+4) - *(fp0+5) * *(tmp+5);
	  //Im part 
	  *(fp1+1+2*d) = *(fp0) * *(tmp+1) + *(fp0+1) * *(tmp) + 
	    *(fp0+2) * *(tmp+3) + *(fp0+3)* *(tmp+2) + 
	    *(fp0+4) * *(tmp+5) + *(fp0+5) * *(tmp+4);
	}
      if(fabs(*fp1)>1e-6){
        PRINT("%p %p tmpfrm[%d]=%e\n",tmpfrm,fp1, agg[s].dest/48,*fp1);
      }
#endif

      //printf("agg[s].dest = %d  agg[s].src = %d, agg[s].gauge = %d\n", agg[s].dest, agg[s].src, agg[s].gauge);
#if 0
      for(int d = 0; d < 6; d++)
      {
	printf("*(uu+%d) = %e %e %e\n",d,*(uu+3*d),*(uu+3*d+1),*(uu+3*d+2));
      	printf("*(fp0+%d) = %e\n",d,*(fp0+d));
      	printf("*(fp1+%d) = %e\n",d,*(fp1+d));
      }
#endif

    }
}


extern "C"
void pt_cmv_dag_cpp(int sites, ind_agg *u, double *gauge_field, double *a, double *tmpfrm)
{
  IFloat *fp0, *fp1, *uu, *tmp;
  int s,c,d;
  struct gauge_agg_cb *agg = u;
  PRINT("sites=%d\n",sites);

  for(s=0;s<sites;s++)
    {
      //Assignment of source, destination, gauge field indexes
      fp1 = (IFloat *)((long)tmpfrm + agg[s].dest);
      fp0 = (IFloat *)((long)a + agg[s].src);
      uu = (IFloat *)((long)gauge_field + agg[s].gauge);
      if(fabs(*fp0)>1e-6){
        PRINT("a[%d]=%e\n",agg[s].src/48,*fp0);
        PRINT("uu[%d]=%e\n",agg[s].gauge/144,*uu);
      }

      //d indexes the row of the gauge link matrix
      for(d = 0; d<3; d++)
	{
	  tmp = uu+2*d;
	  //Re part
	  *(fp1+2*d) = *(fp0) * *(tmp) + *(fp0+1) * *(tmp+1) +
	    *(fp0+2) * *(tmp+6) + *(fp0+3) * *(tmp+7) +
	    *(fp0+4) * *(tmp+12) + *(fp0+5) * *(tmp+13);
	  //Im Part
	  *(fp1+2*d+1) = *(fp0+1) * *(tmp) - *(fp0) * *(tmp+1) +
	    *(fp0+3) * *(tmp+6) - *(fp0+2) * *(tmp+7) +
	    *(fp0+5) * *(tmp+12) - *(fp0+4) * *(tmp+13);
	}
      if(fabs(*fp1)>1e-6){
        PRINT("%p %p tmpfrm[%d]=%e\n",tmpfrm,fp1, agg[s].dest/48,*fp1);
      }
    }
}

extern "C"
void pt_cmv_pad_cpp(int sites, ind_agg *u, double *gauge_field, double *a, double *tmpfrm)
{
  int s,c,d;
  struct gauge_agg_cb *agg = u;
  PRINT("sites=%d\n",sites);

  for(s=0;s<sites;s++)
    {
#ifdef USE_SSE
  IFloat *tmp;
  QLA_D3_ColorVector *fp0, *fp1;
  QLA_D3_ColorMatrix *uu ;
//  printf("index=%p src=%p dest=%p gauge=%p\n",u,a,tmpfrm,gauge_field);
      fp1 = (QLA_D3_ColorVector *)((unsigned long)tmpfrm + agg[s].dest);
      fp0 = (QLA_D3_ColorVector *)((unsigned long)a + agg[s].src);
      uu = (QLA_D3_ColorMatrix *)((unsigned long)gauge_field + agg[s].gauge);
      QLA_D3_V_eq_M_times_V(fp1,uu,fp0)
#else
  IFloat *fp0, *fp1, *uu, *tmp;

      //Assignment of source, destination, gauge field indexes
      fp1 = (IFloat *)((long)tmpfrm + agg[s].dest);
      fp0 = (IFloat *)((long)a + agg[s].src);
      uu = (IFloat *)((long)gauge_field + agg[s].gauge);
      if(fabs(*fp0)>1e-6){
        PRINT("a[%d]=%e\n",agg[s].src/48,*fp0);
        PRINT("uu[%d]=%e\n",agg[s].gauge/144,*uu);
      }

      //d indexes the row of the gauge link matrix
      for(d = 0; d<3; d++)
	{
	  tmp = uu+6*d;
	  //Re part
	  *(fp1+2*d) = *(fp0) * *(tmp) - *(fp0+1) * *(tmp+1) + 
	    *(fp0+2) * *(tmp+2) - *(fp0+3)* *(tmp+3) + 
	    *(fp0+4) * *(tmp+4) - *(fp0+5) * *(tmp+5);
	  //Im part 
	  *(fp1+1+2*d) = *(fp0) * *(tmp+1) + *(fp0+1) * *(tmp) + 
	    *(fp0+2) * *(tmp+3) + *(fp0+3)* *(tmp+2) + 
	    *(fp0+4) * *(tmp+5) + *(fp0+5) * *(tmp+4);
	}
      if(fabs(*fp1)>1e-6){
        PRINT("%p %p tmpfrm[%d]=%e\n",tmpfrm,fp1, agg[s].dest/64,*fp1);
      }
#endif
    }
}

extern "C"
void pt_cmv_dag_pad_cpp(int sites, ind_agg *u, double *gauge_field, double *a, double *tmpfrm)
{
  IFloat *fp0, *fp1, *uu, *tmp;
  int s,c,d;
  struct gauge_agg_cb *agg = u;
  PRINT("sites=%d\n",sites);

  for(s=0;s<sites;s++)
    {
      //Assignment of source, destination, gauge field indexes
      fp1 = (IFloat *)((long)tmpfrm + agg[s].dest);
      fp0 = (IFloat *)((long)a + agg[s].src);
      uu = (IFloat *)((long)gauge_field + agg[s].gauge);
      if(fabs(*fp0)>1e-6){
        PRINT("a[%d]=%e\n",agg[s].src/48,*fp0);
        PRINT("uu[%d]=%e\n",agg[s].gauge/144,*uu);
      }

      //printf("agg[s].dest = %d \t agg[s].src = %ld \t agg[s].gauge=%d\n",agg[s].dest/sizeof(IFloat), agg[s].src/sizeof(IFloat),agg[s].gauge/sizeof(IFloat));
      //printf("*(fp0) = %e  *(uu) = %e\n", *(fp0),*(uu));

      //d indexes the row of the gauge link matrix
      for(d = 0; d<3; d++)
	{
	  tmp = uu+2*d;
	  //Re part
	  *(fp1+2*d) = *(fp0) * *(tmp) + *(fp0+1) * *(tmp+1) +
	    *(fp0+2) * *(tmp+6) + *(fp0+3) * *(tmp+7) +
	    *(fp0+4) * *(tmp+12) + *(fp0+5) * *(tmp+13);
	  //Im Part
	  *(fp1+2*d+1) = *(fp0+1) * *(tmp) - *(fp0) * *(tmp+1) +
	    *(fp0+3) * *(tmp+6) - *(fp0+2) * *(tmp+7) +
	    *(fp0+5) * *(tmp+12) - *(fp0+4) * *(tmp+13);
	}
      if(fabs(*fp1)>1e-6){
        PRINT("%p %p tmpfrm[%d]=%e\n",tmpfrm,fp1, agg[s].dest/64,*fp1);
      }
    }
}
//-----------------------------------------------------------------------------

extern "C"
void cmv_agg_cpp( int sites, long u, long a, 
                long tmpfrm)
{
  IFloat *fp0, *fp1, *uu;
  int s, c,i;
  struct gauge_agg *agg = (struct gauge_agg*)u;
  const int SITE_LEN = 72;
  const int MATRIX_SIZE = 18;
  int *ch;

      uu = &(agg[0].mat[0]);
      PRINT("in=%p out=%p ",a,tmpfrm);
      PRINT("src=%p dest=%p uu=%p\n",agg[0].src,agg[0].dest,uu);
      fflush(stdout);
      //printf("uu[%d]= ",0);
      //for(i=0;i<18;i++){
      //printf("%0.4e ",uu[i]);
      //if(i % 6 == 5 ) printf("\n");
      //}
  for (s = 0; s< sites; s++)
	{

      fp1 = (IFloat *)(tmpfrm + agg[s].dest);
      fp0 = (IFloat *)((long)a + agg[s].src);
      uu = &(agg[s].mat[0]);
      if(agg[s].src< 0 || agg[s].src>48*64){
         printf("src=%p\n",agg[s].src);exit(-40);
      }
      if(agg[s].dest< 0 || agg[s].dest>48*64){
         printf("dest=%p\n",agg[s].dest);exit(-40);
      }
//      printf("src=%d dest=%d uu=%p\n",agg[0].src,agg[0].dest,uu);fflush(stdout);
//      printf("fp1=%p fp0=%p\n",fp1,fp0);
//      bzero(fp1,sizeof(IFloat)*6);
#if 0
      printf("uu[%d]= ",s);
      for(i=0;i<18;i++){
      printf("%d %0.4e,",i,uu[i]);
      if(i % 6 == 5 ) printf("\n");
      }
      printf("fp0[%d]= ",s);
      for(i=0;i<6;i++){
      printf("%d %0.4e,",i,fp0[i]);
      if(i % 6 == 5 ) printf("\n");
      }
      printf("fp1[%d]= ",s);
      for(i=0;i<6;i++){
      printf("%d %0.4e,",i,fp1[i]);
      if(i % 6 == 5 ) printf("\n");
      }
      fflush(stdout);
#endif
      for (int d=0; d<3; d++){
	//Re part
	*(fp1+2*d) = *(fp0) * *(uu+6*d) - *(fp0+1) * *(uu+6*d+1) + 
	  *(fp0+2) * *(uu+6*d+2) - *(fp0+3)* *(uu+6*d+3) + 
	  *(fp0+4) * *(uu+6*d+4) - *(fp0+5) * *(uu+6*d+5);
	//Im part 
	*(fp1+1+2*d) = *(fp0) * *(uu+6*d+1) + *(fp0+1) * *(uu+d*6) + 
	  *(fp0+2) * *(uu+6*d+3) + *(fp0+3)* *(uu+6*d+2) + 
	  *(fp0+4) * *(uu+6*d+5) + *(fp0+5) * *(uu+6*d+4);
//      printf("%e %e ",*(fp1+2*c),*(fp1+2*c+1));
      }
 //     printf("\n");
    }

}
#endif
