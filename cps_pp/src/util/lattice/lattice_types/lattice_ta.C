#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of GtadpoleRectFasqtad class.

  $Id: lattice_ta.C,v 1.2 2005-05-12 20:47:01 chulwoo Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Wilson gauge action + asqtad staggered fermion action
//------------------------------------------------------------------
GtadpoleRectFasqtad::GtadpoleRectFasqtad()
	: Fsmear(3)
{
  cname = "GtadpoleRectFasqtad";
  char *fname = "GtadpoleRectFasqtad()";
  VRB.Func(cname,fname);
}

GtadpoleRectFasqtad::~GtadpoleRectFasqtad()
{
  char *fname = "~GtadpoleRectFasqtad()";
  VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
