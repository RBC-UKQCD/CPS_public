//--------------------------------------------------------------------
//  Comparison testing utilities for the testing framework.
//
//
//
//  $Id: testing_framework_cmp.C,v 1.4 2004-08-17 03:33:16 chulwoo Exp $
//--------------------------------------------------------------------


#include <util/qcdio.h>
#include <stdlib.h>	// exit()
#include <math.h>
#include <config.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <comms/glb.h>
#include <alg/w_ginfo.h>

CPS_START_NAMESPACE

static inline int almost_zero(Float x)
{
  int ans ; 

  if( fabs((double) x)  < 0.0000000000001 )
    {
      ans = 1 ;
    }
    else
      { 
	ans = 0 ; 
      }


  return ans ; 
}


void compare_array_relative(Float* pion_corr_A, 
			    Float* pion_corr_B, 
			    Float tol, 
			    int time_size)
{
  printf("Start_regression_test\n") ; 

  int pass_test = 1 ; 
  for(int t= 0 ; t < time_size ; ++t)
    {
      Float pA = pion_corr_A[t] ; 
      Float pB = pion_corr_B[t] ; 

      Float rel ;
      if( ! almost_zero(pA)  )
	{
	 rel = fabs( (pA-pB)/pA   ) ; 
	}
      else
	{
	 rel = fabs(pB) ; 
	}
      if( rel > tol) 
	pass_test = 0 ; 

      printf("|(%g - %g) / %g | = %g\n",pA,pB,pA,rel ) ; 
    }

  // write some text for the test
    printf("test_type       = relative_float_array\n") ; 
    printf("test_toleration = %g\n",tol) ; 

  if (pass_test) 
    printf("test_result = PASS\n") ; 
  else 
    printf("test_result = FAIL\n") ; 


  printf("End_regression_test\n") ; 
}


void compare_float_relative(Float pion_corr_A, 
			    Float pion_corr_B, 
			    Float tol)
{
  printf("Start_regression_test\n") ; 
  int pass_test = 1 ; 
  
  Float pA = pion_corr_A ; 
  Float pB = pion_corr_B ; 

 Float rel ;
 if( ! almost_zero(pA)  )
   {
     rel = fabs( (pA-pB)/pA   ) ; 
   }
 else
   {
     rel = fabs(pB) ; 
   }
 if( rel > tol) 
   pass_test = 0 ; 
 
 printf("|(%g - %g) / %g | = %g\n",pA,pB,pA,rel ) ; 


  // write some text for the test
    printf("test_type       = relative_float\n") ; 
    printf("test_toleration = %g\n",tol) ; 

  if (pass_test) 
    printf("test_result = PASS\n") ; 
  else 
    printf("test_result = FAIL\n") ; 


  printf("End_regression_test\n") ; 
}


CPS_END_NAMESPACE
