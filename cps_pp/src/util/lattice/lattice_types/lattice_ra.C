#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of GimprRectFasqtad class.

  $Id: lattice_ra.C,v 1.5 2004/08/18 11:58:05 zs Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Wilson gauge action + asqtadstaggered fermion action
//------------------------------------------------------------------
GimprRectFasqtad::GimprRectFasqtad()
	: Fsmear(3)
{
  cname = "GimprRectFasqtad";
  char *fname = "GimprRectFasqtad()";
  VRB.Func(cname,fname);
}

GimprRectFasqtad::~GimprRectFasqtad()
{
  char *fname = "~GimprRectFasqtad()";
  VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
