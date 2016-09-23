#include<config.h>

#ifdef USE_BFM

CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of GnoneFasqtad class.

*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/lattice/fbfm.h>
#include <util/verbose.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// No gauge action + BFM fermion action
//------------------------------------------------------------------
GnoneFbfm::GnoneFbfm():cname("GnoneFbfm")
{
}

GnoneFbfm::~GnoneFbfm()
{
}

CPS_END_NAMESPACE

#endif
