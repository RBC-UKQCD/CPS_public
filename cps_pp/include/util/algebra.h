#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief  Auxilary algebra.
 
  $Id: algebra.h,v 1.2 2004/09/21 18:07:15 chulwoo Exp $
*/
//------------------------------------------------------------------
 
 
#ifndef INCLUDED_ALGEBRA_H
#define INCLUDED_ALGEBRA_H

void m_multiply2r( Float *AB, Float *B );
void m_multiply2l( Float *AB, Float *A );
void m_multiply3( Float *AB, Float *A, Float *B );
void m_add( Float *AplusB, Float *A, Float *B );
void m_equal( Float *A, Float *B );
void m_identity( Float *A );
void m_conjugate( Float *A );
void m_invert( Float *matrix );
void m_rand( Float* eps, Float squeeze );
Float m_determinantR( Float* mtx );
Float m_determinantI( Float* mtx );
Float absR( Float x );
void m_subtract( Float* AminusB, Float* A, Float* B );
void m_zero( Float* x );
void m_print( Float *A );
static const int MATRIXSIZE = 18;

#endif
CPS_END_NAMESPACE
