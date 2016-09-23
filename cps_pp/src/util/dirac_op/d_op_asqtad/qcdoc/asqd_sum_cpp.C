/*!\file
    Asqtad Dirac operator code for QCDOC.

    $Id: asqd_sum_cpp.C,v 1.3 2008/02/08 18:35:07 chulwoo Exp $
*/

#include "asq_data_types.h"
#include "asqtad_int.h"
#if 0
#include <util/gjp.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/vector.h>
#include <sysfunc_cps.h>
//#include <stdio.h>


struct gauge_agg{
  int src;
  int dest;
  IFloat mat[18];
};

extern void SCUTransCRC(SCUDirArg *X, SCUDirArg *R);

void dirac_cmv_nl_mod_cpp( int sites, long chi, long u, long a, 
                long tmpfrm)
{
  IFloat *fp0, *fp1, *uu;
  int s, c;
  const int MATRIX_SIZE = 18;
  int *ch;
  ch = (int *)chi;
  for (s = 0; s< sites; s++)
    {
      fp1 = (IFloat *)(tmpfrm + (int)*(ch + s ));
      fp0 = (IFloat *)((int)a +  s*48 );
      uu = (IFloat *)u + MATRIX_SIZE * s;
      for (c=0; c<3; c++){
	//Re part
	*(fp1+2*c) += *fp0 * *(uu+6*c) - *(fp0+1) * *(uu+6*c+1) + 
	  *(fp0+2) * *(uu+6*c+2) - *(fp0+3)* *(uu+6*c+3) + 
	  *(fp0+4) * *(uu+6*c+4) - *(fp0+5) * *(uu+6*c+5);
	//Im part 
	*(fp1+2*c+1) += *fp0 * *(uu+6*c+1) + *(fp0+1) * *(uu+6*c) + 
	  *(fp0+2) * *(uu+6*c+3) + *(fp0+3)* *(uu+6*c+2) + 
	  *(fp0+4) * *(uu+6*c+5) + *(fp0+5) * *(uu+6*c+4);
      }
    }
}

void dirac_cmv_mod_cpp( int sites, long chi, long u, long a, 
                long tmpfrm)
{
  IFloat *fp0, *fp1, *uu;
  int s, c;
  const int MATRIX_SIZE = 18;
  int *ch;
  ch = (int *)chi;
  for (s = 0; s< sites; s++)
    {
      fp1 = (IFloat *)(tmpfrm + (int)*(ch + 2 * s + 1));
      fp0 = (IFloat *)((int)a +  (int)*(ch + 2 * s ));
      uu = (IFloat *)u + MATRIX_SIZE * s;
      for (c=0; c<3; c++){
	//Re part
	*(fp1+2*c) += *fp0 * *(uu+6*c) - *(fp0+1) * *(uu+6*c+1) + 
	  *(fp0+2) * *(uu+6*c+2) - *(fp0+3)* *(uu+6*c+3) + 
	  *(fp0+4) * *(uu+6*c+4) - *(fp0+5) * *(uu+6*c+5);
	//Im part 
	*(fp1+2*c+1) += *fp0 * *(uu+6*c+1) + *(fp0+1) * *(uu+6*c) + 
	  *(fp0+2) * *(uu+6*c+3) + *(fp0+3)* *(uu+6*c+2) + 
	  *(fp0+4) * *(uu+6*c+5) + *(fp0+5) * *(uu+6*c+4);
      }
    }
}

void dirac_cmv_jcw_agg_cpp( int sites, long chi, long u, long a, 
                long tmpfrm)
{
  IFloat *fp0, *fp1, *uu;
  int s, c,i;
  struct gauge_agg *agg = (struct gauge_agg*)u;
//  const int SITE_LEN = 72;
//  const int MATRIX_SIZE = 18;
  int *ch;
  printf("dirac_cmv_jcw_agg_cpp:\n");
  for (s = 0; s< sites; s++)
    {
      fp1 = (IFloat *)(tmpfrm + agg[s].dest );
      fp0 = (IFloat *)((int)a + agg[s].src);
      uu = &(agg[s].mat[0]);
   if( (unsigned long)fp1 >0xc0000000 ||
      (unsigned long)fp0 >0xc0000000 ||
      (unsigned long)uu >0xc0000000){
      printf("fp1=%p fp0=%p uu=%p out of bound!!\n",fp1,fp0,uu);
      exit(4);
   }

	for (c=0; c<3; c++){
	    //Re part
	    *(fp1+2*c) = *fp0 * *(uu+6*c) - *(fp0+1) * *(uu+6*c+1) + 
		*(fp0+2) * *(uu+6*c+2) - *(fp0+3)* *(uu+6*c+3) + 
		*(fp0+4) * *(uu+6*c+4) - *(fp0+5) * *(uu+6*c+5);
	    //Im part 
	    *(fp1+2*c+1) = *fp0 * *(uu+6*c+1) + *(fp0+1) * *(uu+6*c) + 
		*(fp0+2) * *(uu+6*c+3) + *(fp0+3)* *(uu+6*c+2) + 
		*(fp0+4) * *(uu+6*c+5) + *(fp0+5) * *(uu+6*c+4);
	}
#if 0
	printf("dirac_cmv_jcw_agg_cpp:\n");
	if(fp0[0]*fp0[0]>1e-5){
	    printf("src=%p(%d) dest=%p(%d) ",fp0,agg[s].src,fp1,agg[s].dest);
	    fflush(stdout);

	    printf("src=");
	    for(i=0;i<6;i++){
		printf("%d %0.4e,",i,fp0[i]);
		if(i % 6 == 5 ) printf("\n");
	    }

	    printf("uu=");
	    for(i=0;i<18;i++){
		printf("%d %0.4e,",i,uu[i]);
		if(i % 6 == 5 ) printf("\n");
	    }

	    printf("dest=");
	    for (c=0; c<3; c++)
		printf("%e %e ",*(fp1+2*c),*(fp1+2*c+1));
	    printf("\n");
	}
#endif

#if 0
   if(fp0[0]*fp0[0]>1e-5){
      printf("src=%p(%d) dest=%p(%d) uu=%p",fp0,agg[s].src,fp1,agg[s].dest,uu);
      fflush(stdout);
#if 1
      printf("src=");
      for(i=0;i<6;i++){
      printf("%d %0.4e,",i,fp0[i]);
      if(i % 6 == 5 ) printf("\n");
      }
#endif
#if 1
      printf("uu=");
      for(i=0;i<18;i++){
      printf("%d %0.4e,",i,uu[i]);
      if(i % 6 == 5 ) printf("\n");
      }
#endif
      printf("dest=");
      for (c=0; c<3; c++)
      printf("%e %e ",*(fp1+2*c),*(fp1+2*c+1));
      printf("\n");
  }
#endif
    }
  
}


void dirac_cmv_cpp( int sites, long chi, long u, long a, 
                long tmpfrm)
{
  IFloat *fp0, *fp1, *uu;
  int s, c;
  const int MATRIX_SIZE = 18;
  int *ch;
  ch = (int *)chi;

  for (s = 0; s< sites; s++)
    {
      fp1 = (IFloat *)(tmpfrm + (int)*(ch + 2 * s + 1));
      fp0 = (IFloat *)((int)a +  (int)*(ch + 2 * s ));
      uu = (IFloat *)u + MATRIX_SIZE * s;
      for (c=0; c<3; c++){
	//Re part
	*(fp1+2*c) = *fp0 * *(uu+6*c) - *(fp0+1) * *(uu+6*c+1) + 
	  *(fp0+2) * *(uu+6*c+2) - *(fp0+3)* *(uu+6*c+3) + 
	  *(fp0+4) * *(uu+6*c+4) - *(fp0+5) * *(uu+6*c+5);
	//Im part 
	*(fp1+2*c+1) = *fp0 * *(uu+6*c+1) + *(fp0+1) * *(uu+6*c) + 
	  *(fp0+2) * *(uu+6*c+3) + *(fp0+3)* *(uu+6*c+2) + 
	  *(fp0+4) * *(uu+6*c+5) + *(fp0+5) * *(uu+6*c+4);
      }
    }
}

void dirac_cmv_nl( int sites, long chi, long u, long a, 
                long tmpfrm)
{
  IFloat *fp0, *fp1, *uu;
  int s, c;

  const int MATRIX_SIZE = 18;
  int *ch;
  ch = (int *)chi;

  for (s = 0; s< sites; s++)
    {
	  fp1 = (IFloat *)tmpfrm + (int)*(ch + 2 * s + 1);
	  fp0 = (IFloat *)((int)/*(long)*/0 +(int)*(ch + 2 * s ));
  	  uu = (IFloat *)u + MATRIX_SIZE * s;
//  	  printf("s_nl=%d, fp0(src) = %d, fp1(sol) = %d, uu = %d \n", s, fp0, fp1, uu);
	  for (c=0; c<3; c++){
/*
        if( fabs(*(fp0+0))>1e-6 || fabs(*(fp0+1))>1e-6 )
           printf("%e %e\n",*(fp0+0), *(fp0+1));
        if( fabs(*(fp0+2))>1e-6 || fabs(*(fp0+3))>1e-6 )
           printf("%e %e\n",*(fp0+2), *(fp0+3));
        if( fabs(*(fp0+4))>1e-6 || fabs(*(fp0+5))>1e-6 )
           printf("%e %e\n",*(fp0+4), *(fp0+5));
*/
	    //Re part
	    *(fp1+2*c) = *fp0 * *(uu+6*c) - *(fp0+1) * *(uu+6*c+1) + 
	      *(fp0+2) * *(uu+6*c+2) - *(fp0+3)* *(uu+6*c+3) + 
	      *(fp0+4) * *(uu+6*c+4) - *(fp0+5) * *(uu+6*c+5);
	    //Im part 
	    *(fp1+2*c+1) = *fp0 * *(uu+6*c+1) + *(fp0+1) * *(uu+6*c) + 
	      *(fp0+2) * *(uu+6*c+3) + *(fp0+3)* *(uu+6*c+2) + 
	      *(fp0+4) * *(uu+6*c+5) + *(fp0+5) * *(uu+6*c+4);
	  }
    }
}
#endif

void asqd_sum_acc_cpp(int s, long chi, long tmpfrm, long b)
{
  int i,j;
  IFloat *fp0, *fp1;
  int *ch;
  ch = (int *)chi;

  for (i = 0; i < s; i++){
    fp0 = (IFloat *)((int)b + (int)*(ch +17 * i));
    for (j = 0; j < 16; j++){
      fp1 = (IFloat *)(tmpfrm + (int)*(ch + 17 * i + j + 1));
//        printf(" *( chi + %d ) = %x \n", 9*i + j+1, 
//    	       *( ch + 9 * i+ j+1));
      *fp0 += *fp1;
      *(fp0+1) += *(fp1+1);
      *(fp0+2) += *(fp1+2);
      *(fp0+3) += *(fp1+3);
      *(fp0+4) += *(fp1+4);
      *(fp0+5) += *(fp1+5);
    }
  }
}

#if 0

void dirac_sum2_64_cpp(int s, long chi, long tmpfrm, long b)
{
  int i,j;
  IFloat *fp0, *fp1;
  int *ch;
  ch = (int *)chi;
  double sum;


  for (i = 0; i < s; i++){
    sum = 0.0;
    fp0 = (IFloat *)((int)b + 48*i);
    fp1 = (IFloat *)(tmpfrm + i*64*16);
    *fp0 = *fp1;
    *(fp0+1) = *(fp1+1);
    *(fp0+2) = *(fp1+2);
    *(fp0+3) = *(fp1+3);
    *(fp0+4) = *(fp1+4);
    *(fp0+5) = *(fp1+5);

    for (j = 1; j < 16; j++){
      fp1 = (IFloat *)(tmpfrm + 64*(i*16+j) );
      *fp0 += *fp1;
      *(fp0+1) += *(fp1+1);
      *(fp0+2) += *(fp1+2);
      *(fp0+3) += *(fp1+3);
      *(fp0+4) += *(fp1+4);
      *(fp0+5) += *(fp1+5);
    }
  }
}

void dirac_sum2_cpp(int s, long chi, long tmpfrm, long b)
{
  int i,j;
  IFloat *fp0, *fp1;
  int *ch;
  ch = (int *)chi;


  for (i = 0; i < s; i++){
    fp0 = (IFloat *)((int)b + 48*i);
    fp1 = (IFloat *)(tmpfrm + i*48*16);
    *fp0 = *fp1;
    *(fp0+1) = *(fp1+1);
    *(fp0+2) = *(fp1+2);
    *(fp0+3) = *(fp1+3);
    *(fp0+4) = *(fp1+4);
    *(fp0+5) = *(fp1+5);

    for (j = 1; j < 16; j++){
      fp1 = (IFloat *)(tmpfrm + 48*(i*16+j) );
      *fp0 += *fp1;
      *(fp0+1) += *(fp1+1);
      *(fp0+2) += *(fp1+2);
      *(fp0+3) += *(fp1+3);
      *(fp0+4) += *(fp1+4);
      *(fp0+5) += *(fp1+5);
    }
  }
}

void dirac_sum_cpp(int s, long chi, long tmpfrm, long b)
{
  int i,j;
  IFloat *fp0, *fp1;
  int *ch;
  ch = (int *)chi;


  for (i = 0; i < s; i++){
    fp0 = (IFloat *)((int)b + (int)*(ch +17 * i));
    fp1 = (IFloat *)(tmpfrm + (int)*(ch + 17 * i + 1));
    *fp0 = *fp1;
    *(fp0+1) = *(fp1+1);
    *(fp0+2) = *(fp1+2);
    *(fp0+3) = *(fp1+3);
    *(fp0+4) = *(fp1+4);
    *(fp0+5) = *(fp1+5);

    for (j = 1; j < 16; j++){
      fp1 = (IFloat *)(tmpfrm + (int)*(ch + 17 * i + j + 1));
      *fp0 += *fp1;
      *(fp0+1) += *(fp1+1);
      *(fp0+2) += *(fp1+2);
      *(fp0+3) += *(fp1+3);
      *(fp0+4) += *(fp1+4);
      *(fp0+5) += *(fp1+5);
    }
  }
}

void dirac_cmm_jcw_agg_cpp( int sites, long chi, long u, long a, 
                long tmpfrm)
{
  IFloat *fp0, *fp1, *uu;
  int s, c;
  struct gauge_agg *agg = (struct gauge_agg*)u;
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

#if 0
void dirac_SCU( SCUDirArg ** Xarg, SCUDirArg ** Rarg, IFloat * a,
  int a_odd, int * Xoffset, int ** ToffsetP, int ** ToffsetM,
  int * countP, int * countM )
{
  int i,j;
  IFloat *ar;
  IFloat *ax;


  for (i=0; i<countP[a_odd]; i++)
    {
      ar = (IFloat *)(*(Rarg+4))->Addr() + (6 * i);
      ax = a + (int)*(ToffsetP[a_odd]+i);
      for (j=0; j< 6; j++)
	{
	  *ar++ = *ax++;
	}
    } 

  for (i=0; i<countM[a_odd]; i++)
    {
      ar = (IFloat *)(*(Rarg))->Addr() + 6 * i;
      ax = a + (int)*(ToffsetM[a_odd]+i);
      for (j=0; j< 6; j++)
	{
	  *ar++ = *ax++;
	}
    }   

  for (i=1; i<4; i++)
    {
      (*(Xarg+i))->Addr(a+Xoffset[i]);
//    printf("Transfering in %d direction \n", i);
      SCUTransCRC(*(Xarg+i), *(Rarg+i+4));
      (*(Xarg+i+4))->Addr(a+Xoffset[i+4]);
//    printf("Transfering in %d direction \n", i+4);
      SCUTransCRC(*(Xarg+i+4), *(Rarg+i));

    }
}
#endif

CPS_END_NAMESPACE
#endif

