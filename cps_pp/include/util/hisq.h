#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief Utility routines for the hisq fermion Dirac operator
*/

/****************************************************************************/
/*                                                                          */
/* hisq.h                                                                     */
/*                                                                          */
/* C header file for the hisq improved fermion library.                       */
/*                                                                          */
/****************************************************************************/

#ifndef INCLUDED_HISQ_LIB_H
#define INCLUDED_HISQ_LIB_H

//------------------------------------------------------------------
// The staggered hisq dslash operator
//------------------------------------------------------------------
extern "C" void hisq_dirac(Vector *f_out, Vector *f_in, int cb, int dag);

//------------------------------------------------------------------
// Initialize all global variables and address tables needed by
// staggered dirac() function.
//------------------------------------------------------------------
extern "C" void hisq_dirac_init(const void *gauge_field_addr);
extern "C" void hisq_dirac_init_g ();

//------------------------------------------------------------------
// Destroy all address tables needed by staggered dirac() function
//------------------------------------------------------------------
extern "C" void hisq_destroy_dirac_buf();
extern "C" void hisq_destroy_dirac_buf_g();

//------------------------------------------------------------------
// The staggered hisq dMdmu operator
//------------------------------------------------------------------
extern "C" void hisq_dMdmu(Vector *f_out, Vector *f_in, int cb, int dag, int order);


#endif

CPS_END_NAMESPACE
