#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of the GpowerRectFasqtad class.

  $Id: lattice_va.C,v 1.5 2004-08-18 11:58:06 zs Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// PowerRect gauge action + asqtad staggered fermion action
//------------------------------------------------------------------
GpowerRectFasqtad::GpowerRectFasqtad()
	: Fsmear(3)
{
  cname = "GpowerRectFasqtad";
  char *fname = "GpowerRectFasqtad()";
  VRB.Func(cname,fname);
}

GpowerRectFasqtad::~GpowerRectFasqtad()
{
  char *fname = "~GpowerRectFasqtad()";
  VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
