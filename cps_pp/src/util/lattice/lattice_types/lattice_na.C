#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of GnoneFasqtad class.

  $Id: lattice_na.C,v 1.5 2004-08-18 11:58:05 zs Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// No gauge action + asqtad staggered fermion action
//------------------------------------------------------------------
GnoneFasqtad::GnoneFasqtad()
	: Gnone(),Fasqtad(),Fsmear(3)
{
  cname = "GnoneFasqtad";
  char *fname = "GnoneFasqtad()";
  VRB.Func(cname,fname);
}

GnoneFasqtad::~GnoneFasqtad()
{
  char *fname = "~GnoneFasqtad()";
  VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
