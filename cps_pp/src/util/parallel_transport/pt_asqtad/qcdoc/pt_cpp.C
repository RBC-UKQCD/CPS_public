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
      fp1 = (IFloat *)(tmpfrm + 3*agg[s].dest );
      fp0 = (IFloat *)((int)a + 3*agg[s].src);
      uu = &(agg[s].mat[0]);
#if 0
      for(i=0;i<18;i++){
      printf("%d %0.4e,",i,uu[i]);
      if(i % 6 == 5 ) printf("\n");
      }
#endif
//      printf("src=%d dest=%d ",agg[s].src,agg[s].dest);
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
 //     printf("\n");
    }
}
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
      printf("in=%p out=%p ",a,tmpfrm);
      printf("src=%p dest=%p uu=%p\n",agg[0].src,agg[0].dest,uu);fflush(stdout);
      printf("uu[%d]= ",0);
      for(i=0;i<18;i++){
      printf("%0.4e ",uu[i]);
      if(i % 6 == 5 ) printf("\n");
      }
  for (s = 0; s< sites; s++)
	{
      fp1 = (IFloat *)(tmpfrm + agg[s].dest );
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
	*(fp1+1+2*d) = *(fp0) * *(uu+6*d+1) + *(fp0+1+2*d) * *(uu+d*6) + 
	  *(fp0+2) * *(uu+6*d+3) + *(fp0+3)* *(uu+6*d+2) + 
	  *(fp0+4) * *(uu+6*d+5) + *(fp0+5) * *(uu+6*d+4);
//      printf("%e %e ",*(fp1+2*c),*(fp1+2*c+1));
      }
 //     printf("\n");
    }
}
