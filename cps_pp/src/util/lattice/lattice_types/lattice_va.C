#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of the GpowerRectFstagAsqtad class.

  $Id: lattice_va.C,v 1.2 2003-10-23 13:38:59 zs Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// PowerRect gauge action + asqtad staggered fermion action
//------------------------------------------------------------------
GpowerRectFstagAsqtad::GpowerRectFstagAsqtad()
{
  cname = "GpowerRectFstagAsqtad";
  char *fname = "GpowerRectFstagAsqtad()";
  VRB.Func(cname,fname);
}

GpowerRectFstagAsqtad::~GpowerRectFstagAsqtad()
{
  char *fname = "~GpowerRectFstagAsqtad()";
  VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
