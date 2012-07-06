#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of GnoneFasqtad class.

  $Id: lattice_nb.C,v 1.2 2012-07-06 20:22:08 chulwoo Exp $
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
