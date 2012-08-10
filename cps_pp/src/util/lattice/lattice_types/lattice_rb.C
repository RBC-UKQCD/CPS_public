#include<config.h>

#ifdef USE_BFM

CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of GnoneFasqtad class.

  $Id: lattice_rb.C,v 1.3 2012-08-10 14:05:33 chulwoo Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/lattice/fbfm.h>
#include <util/verbose.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// Iwasaki gauge action + BFM fermion action
//------------------------------------------------------------------
GimprRectFbfm::GimprRectFbfm():cname("GimprRectFbfm")
{
}

GimprRectFbfm::~GimprRectFbfm()
{
}

CPS_END_NAMESPACE

#endif
