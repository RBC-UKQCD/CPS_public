#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of the GpowerRectFasqtad class.

  $Id: lattice_va.C,v 1.3 2003-10-30 05:39:12 cwj Exp $
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
