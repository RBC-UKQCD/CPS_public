#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief Utility routines for the p4 fermion Dirac operator
*/

/****************************************************************************/
/*                                                                          */
/* p4.h                                                                     */
/*                                                                          */
/* C header file for the p4 improved fermion library.                       */
/*                                                                          */
/****************************************************************************/

#ifndef INCLUDED_P4_LIB_H
#define INCLUDED_P4_LIB_H

//------------------------------------------------------------------
// The staggered p4 dslash operator
//------------------------------------------------------------------
extern "C" void p4_dirac(Vector *f_out, Vector *f_in, int cb, int dag);

//------------------------------------------------------------------
// Initialize all global variables and address tables needed by
// staggered dirac() function.
//------------------------------------------------------------------
extern "C" void p4_dirac_init(const void *gauge_field_addr);
extern "C" void p4_dirac_init_g ();

//------------------------------------------------------------------
// Destroy all address tables needed by staggered dirac() function
//------------------------------------------------------------------
extern "C" void p4_destroy_dirac_buf();
extern "C" void p4_destroy_dirac_buf_g();

//------------------------------------------------------------------
// The staggered p4 dMdmu operator
//------------------------------------------------------------------
extern "C" void p4_dMdmu(Vector *f_out, Vector *f_in, int cb, int dag, int order);


#endif

CPS_END_NAMESPACE
