#if 0
#include <util/gjp.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/vector.h>
#include <sysfunc.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#endif
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
