#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of GimprRectFstagAsqtad class.

  $Id: lattice_ra.C,v 1.2 2003-10-23 13:38:59 zs Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Wilson gauge action + asqtadstaggered fermion action
//------------------------------------------------------------------
GimprRectFstagAsqtad::GimprRectFstagAsqtad()
{
  cname = "GimprRectFstagAsqtad";
  char *fname = "GimprRectFstagAsqtad()";
  VRB.Func(cname,fname);
}

GimprRectFstagAsqtad::~GimprRectFstagAsqtad()
{
  char *fname = "~GimprRectFstagAsqtad()";
  VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
