/*!----------------------------------------------------------------------
  $Id: dirac_serial_cpp.C,v 1.5 2004-09-21 18:11:16 chulwoo Exp $

----------------------------------------------------------------------*/


#include <config.h>
#include <stdio.h>
#include <util/data_types.h>

CPS_START_NAMESPACE
struct gauge_agg{
  int src;
  int dest;
  IFloat mat[18];
};


void dirac_cmv_nl_mod_cpp( int sites, long chi, long u, long a, 
                long tmpfrm)
{
  IFloat *fp0, *fp1, *uu;
  int s, c;
  const int SITE_LEN = 72;
  const int MATRIX_SIZE = 18;
  int *ch;
  ch = (int *)chi;
  for (s = 0; s< sites; s++)
    {
      fp1 = (IFloat *)(tmpfrm + (int)*(ch + s ));
      fp0 = (IFloat *)(a +  s*6*sizeof(IFloat) );      
      
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
  const int SITE_LEN = 72;
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
  const int SITE_LEN = 72;
  const int MATRIX_SIZE = 18;
  int *ch;

  for (s = 0; s< sites; s++)
    {
      fp1 = (IFloat *)(tmpfrm + agg[s].dest );
      fp0 = (IFloat *)((int)a + agg[s].src);
      uu = &(agg[s].mat[0]);
#if 0
      for(i=0;i<18;i++){
      printf("%d %0.4e,",i,uu[i]);
      if(i % 6 == 5 ) printf("\n");
      }
#endif
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


void dirac_cmv_cpp( int sites, long chi, long u, long a, 
                long tmpfrm)
{
  IFloat *fp0, *fp1, *uu;
  int s, c;
  const int SITE_LEN = 72;
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
  const int SITE_LEN = 72;
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

void dirac_sum_acc_cpp(int s, long chi, long tmpfrm, long b)
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

void dirac_sum2_64_cpp(int s, long chi, long tmpfrm, long b)
{
  int i,j;
  IFloat *fp0, *fp1;


  for (i = 0; i < s; i++){
//     fp0 = (IFloat *)((int)b + 48*i);
//     fp1 = (IFloat *)(tmpfrm + i*64*16);
    fp0 = (IFloat *)((int)b + 6*sizeof(IFloat)*i);
    fp1 = (IFloat *)(tmpfrm + i*8*sizeof(IFloat)*16);
    *fp0 = *fp1;
    *(fp0+1) = *(fp1+1);
    *(fp0+2) = *(fp1+2);
    *(fp0+3) = *(fp1+3);
    *(fp0+4) = *(fp1+4);
    *(fp0+5) = *(fp1+5);

    for (j = 1; j < 16; j++){
//       fp1 = (IFloat *)(tmpfrm + 64*(i*16+j) );
      fp1 = (IFloat *)(tmpfrm + 8*sizeof(IFloat)*(i*16+j) );	
      *fp0 += *fp1;
      *(fp0+1) += *(fp1+1);
      *(fp0+2) += *(fp1+2);
      *(fp0+3) += *(fp1+3);
      *(fp0+4) += *(fp1+4);
      *(fp0+5) += *(fp1+5);
    }
  }
}

void dirac_sum2_acc_cpp(int s, long chi, long tmpfrm, long b)
{
  int i,j;
  IFloat *fp0, *fp1;


  for (i = 0; i < s; i++){
    fp0 = (IFloat *)((int)b + 6*sizeof(IFloat)*i);
    fp1 = (IFloat *)(tmpfrm + i*6*sizeof(IFloat)*16);

    
    *fp0 = *fp1;
    *(fp0+1) += *(fp1+1);
    *(fp0+2) += *(fp1+2);
    *(fp0+3) += *(fp1+3);
    *(fp0+4) += *(fp1+4);
    *(fp0+5) += *(fp1+5);

    for (j = 1; j < 16; j++){
      fp1 = (IFloat *)(tmpfrm + 6*sizeof(IFloat)*(i*16+j) );	
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
 

  for (i = 0; i < s; i++){
    fp0 = (IFloat *)(b + 6*sizeof(IFloat)*i);
    fp1 = (IFloat *)(tmpfrm + i*6*sizeof(IFloat)*16);
    *fp0 = *fp1;
    *(fp0+1) = *(fp1+1);
    *(fp0+2) = *(fp1+2);
    *(fp0+3) = *(fp1+3);
    *(fp0+4) = *(fp1+4);
    *(fp0+5) = *(fp1+5);

    for (j = 1; j < 16; j++){
      fp1 = (IFloat *)(tmpfrm + 6*sizeof(IFloat)*(i*16+j) );	
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

void copy_buffer_cpp(int n, long src, long dest, long ptable){
  int i,j;
  IFloat *fp1,*fp0;
  int * chi = (int *)ptable;
  for(i=0;i<n;i++){
    fp1 = (IFloat *) (src + chi[i]*sizeof(IFloat));
    fp0 = (IFloat *) (dest + i*6*sizeof(IFloat));
    for(j=0;j<6;j++)
      *(fp0+j) = *(fp1+j);
  }
}

CPS_END_NAMESPACE
