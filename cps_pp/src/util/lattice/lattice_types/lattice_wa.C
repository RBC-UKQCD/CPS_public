#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of the GwilsonFstagAsqtad class.

  $Id: lattice_wa.C,v 1.2 2003-10-23 13:38:59 zs Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Wilson gauge action + asqtad staggered fermion action
//------------------------------------------------------------------
GwilsonFstagAsqtad::GwilsonFstagAsqtad()
{
  cname = "GwilsonFstagAsqtad";
  char *fname = "GwilsonFstagAsqtad()";
  VRB.Func(cname,fname);
}

GwilsonFstagAsqtad::~GwilsonFstagAsqtad()
{
  char *fname = "~GwilsonFstagAsqtad()";
  VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
