#if 0
#include <util/gjp.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/vector.h>
#include <sysfunc.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#endif
#include <config.h>
#include <util/vector.h>
#include <stdio.h>
typedef double IFloat;

struct gauge_agg{
  int src;
  int dest;
  IFloat mat[18];
};

//---------------------------------------------------------------
//Checkerboarding
struct gauge_agg_cb{
  int src;
  int dest;
  int gauge;
};
//---------------------------------------------------------------


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

extern "C"
void pt_cmm_cpp(int sites, long u, long a, long tmpfrm, long gauge_field)
{
  IFloat *fp0, *fp1, *uu;
  int s,c,d;
  struct gauge_agg_cb *agg = (struct gauge_agg_cb*)u;

  for(s=0;s<sites;s++)
    {
      //Assignment of source,destination,gauge field index
      fp1 = (IFloat *)(tmpfrm + 3*agg[s].dest);
      fp0 = (IFloat *)((int)a + 3*agg[s].src);
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
  struct gauge_agg_cb *agg = (struct gauge_agg_cb*)u;

  for(s=0;s<sites;s++)
    {
      //Assignment of source,destination,gauge field index
      fp1 = (IFloat *)(tmpfrm + 3*agg[s].dest);
      fp0 = (IFloat *)((int)a + 3*agg[s].src);
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
void cmm_agg_cpp( int sites, long chi, long u, long a, 
                long tmpfrm)
{
  IFloat *fp0, *fp1, *uu;
  int s, c,i;
  struct gauge_agg *agg = (struct gauge_agg*)u;
  const int SITE_LEN = 72;
  const int MATRIX_SIZE = 18;
  int *ch;

  for (s = 0; s< sites; s++)
    {
      //Assignment of source, destination indexes
      fp1 = (IFloat *)(tmpfrm + 3*agg[s].dest );
      fp0 = (IFloat *)((int)a + 3*agg[s].src);
      uu = &(agg[s].mat[0]);
#if 0
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
	//      printf("%e %e ",*(fp1+2*c),*(fp1+2*c+1));
      }
      //           printf("\n");
    }
}

//-----------------------------------------------------------------------------
//C++ implementation of matrix vector multiply for checkerboarded field
//
//Parameters
//
//sites - Total number of matrix multiplies
//u - Pointer to structs that hold source, destination indexes for parallel 
//    transport, also holds the index for the appropriate SU(3) gauge link
//a - Vector field that is to be transported
//tmpfrm - Transported vector field
//gauge_field - pointer to the gauge field

extern "C"
void pt_cmv_cpp(int sites, long u, long a, long tmpfrm, long gauge_field)
{
  IFloat *fp0, *fp1, *uu, *tmp;
  int s,c,d;
  struct gauge_agg_cb *agg = (struct gauge_agg_cb*)u;

  for(s=0;s<sites;s++)
    {
      //Assignment of source, destination, gauge field indexes
      fp1 = (IFloat *)(tmpfrm + agg[s].dest);
      fp0 = (IFloat *)((int)a + agg[s].src);
      uu = (IFloat *)(gauge_field + agg[s].gauge);

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

      #if 0
      printf("agg[s].dest = %d  agg[s].src = %d, agg[s].gauge = %d\n", agg[s].dest, agg[s].src, agg[s].gauge);
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
void pt_cmv_dag_cpp(int sites, long u, long a, long tmpfrm, long gauge_field)
{
  IFloat *fp0, *fp1, *uu, *tmp;
  int s,c,d;
  struct gauge_agg_cb *agg = (struct gauge_agg_cb*)u;

  for(s=0;s<sites;s++)
    {
      //Assignment of source, destination, gauge field indexes
      fp1 = (IFloat *)(tmpfrm + agg[s].dest);
      fp0 = (IFloat *)((int)a + agg[s].src);
      uu = (IFloat *)(gauge_field + agg[s].gauge);

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
    }
}

extern "C"
void pt_cmv_pad_cpp(int sites, long u, long a, long tmpfrm, long gauge_field)
{
  IFloat *fp0, *fp1, *uu, *tmp;
  int s,c,d;
  struct gauge_agg_cb *agg = (struct gauge_agg_cb*)u;

  for(s=0;s<sites;s++)
    {
      //Assignment of source, destination, gauge field indexes
      fp1 = (IFloat *)(tmpfrm + agg[s].dest);
      fp0 = (IFloat *)((int)a + agg[s].src);
      uu = (IFloat *)(gauge_field + agg[s].gauge);

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
    }
}

extern "C"
void pt_cmv_dag_pad_cpp(int sites, long u, long a, long tmpfrm, long gauge_field)
{
  IFloat *fp0, *fp1, *uu, *tmp;
  int s,c,d;
  struct gauge_agg_cb *agg = (struct gauge_agg_cb*)u;
  //printf("Number of sites = %d\n",sites);

  for(s=0;s<sites;s++)
    {
      //Assignment of source, destination, gauge field indexes
      fp1 = (IFloat *)(tmpfrm + agg[s].dest);
      fp0 = (IFloat *)((int)a + agg[s].src);
      uu = (IFloat *)(gauge_field + agg[s].gauge);

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
      //printf("in=%p out=%p ",a,tmpfrm);
      //printf("src=%p dest=%p uu=%p\n",agg[0].src,agg[0].dest,uu);
      //fflush(stdout);
      //printf("uu[%d]= ",0);
      //for(i=0;i<18;i++){
      //printf("%0.4e ",uu[i]);
      //if(i % 6 == 5 ) printf("\n");
      //}
  for (s = 0; s< sites; s++)
	{

      fp1 = (IFloat *)(tmpfrm + agg[s].dest);
      fp0 = (IFloat *)((int)a + agg[s].src);
      uu = &(agg[s].mat[0]);
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
