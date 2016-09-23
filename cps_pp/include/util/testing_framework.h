/*!\file
  \brief Function prototypes for the testing framework

  $Id: testing_framework.h,v 1.4 2004/09/02 16:58:02 zs Exp $
*/

#ifndef TESTING_FRAMEWORK_INC
#define TESTING_FRAMEWORK_INC

#include<config.h>
CPS_START_NAMESPACE


/*! \defgroup testing Testing framework utilities
  @{ */

// ------ light-weight hadron spectrum -------

//! Simple pseudoscalar correlator.

void staggered_local_pion(Lattice &lat, Float mass, 
			  IFloat* pion_corr, int time_size);

//
//  ----- comparison functions -----
//

//! Comparison of two arrays

void compare_array_relative(Float* pion_corr_A, 
			    Float* pion_corr_B, 
			    Float tol, 
			    int time_size) ;

//! Comparison of two floating-point numbers

void compare_float_relative(Float pion_corr_A,
                            Float pion_corr_B,
                            Float tol) ;

/*! @} */

CPS_END_NAMESPACE

#endif

