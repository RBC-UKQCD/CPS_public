#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of GimprRectFasqtad class.

  $Id: lattice_ra.C,v 1.3 2003-10-30 05:39:12 cwj Exp $
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
