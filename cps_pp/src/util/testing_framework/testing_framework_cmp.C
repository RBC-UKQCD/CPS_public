//--------------------------------------------------------------------
/*!\file
  \brief Comparison testing utilities for the testing framework.

  $Id: testing_framework_cmp.C,v 1.6 2004-09-02 16:57:12 zs Exp $
*/
//--------------------------------------------------------------------

#include <config.h>
#include <util/data_types.h>
#include <math.h>
#include <stdio.h>
#include <util/verbose.h>
#include <util/error.h>

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

/*!
  Computes the relative difference between each number in the first array
  and the corresponding number in the second array. The differences are
  printed to \c stdout. If all differences are less than a specified
  tolerance, "test_result = PASS" is printed to \c stdout; otherwise
  "test_result = FAIL" is printed.
  
  \param pion_corr_A An floating-point array
  \param pion_corr_B Another floating-point array
  \param tol The tolerance
  \param time_size The length of the arrays
*/

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

/*!
  Computes the relative difference between two numbers
  The difference is printed to \c stdout.
  If the difference is less than a specified
  tolerance, "test_result = PASS" is printed to \c stdout; otherwise
  "test_result = FAIL" is printed.
  
  \param pion_corr_A An floating-point array
  \param pion_corr_B Another floating-point array
  \param tol The tolerance
*/

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
