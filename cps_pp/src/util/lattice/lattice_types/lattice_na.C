#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of GnoneFstagAsqtad class.

  $Id: lattice_na.C,v 1.2 2003-10-23 13:38:59 zs Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// No gauge action + asqtad staggered fermion action
//------------------------------------------------------------------
GnoneFstagAsqtad::GnoneFstagAsqtad()
{
  cname = "GnoneFstagAsqtad";
  char *fname = "GnoneFstagAsqtad()";
  VRB.Func(cname,fname);
}

GnoneFstagAsqtad::~GnoneFstagAsqtad()
{
  char *fname = "~GnoneFstagAsqtad()";
  VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
