#ifndef TESTING_FRAMEWORK_INC
#define TESTING_FRAMEWORK_INC

#include<config.h>
CPS_START_NAMESPACE

//
// Function prototypes for the testing framework
//


//
//  ----- light weight hadron spectrum -----
//

void staggered_local_pion(Lattice &lat, Float mass, 
			  IFloat* pion_corr, int time_size);

//
//  ----- comparison functions -----
//

void compare_array_relative(Float* pion_corr_A, 
			    Float* pion_corr_B, 
			    Float tol, 
			    int time_size) ;

void compare_float_relative(Float pion_corr_A,
                            Float pion_corr_B,
                            Float tol) ;



CPS_END_NAMESPACE

#endif

