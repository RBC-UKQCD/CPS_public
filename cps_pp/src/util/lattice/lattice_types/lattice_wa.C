#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of the GwilsonFasqtad class.

  $Id: lattice_wa.C,v 1.3 2003-10-30 05:39:12 cwj Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Wilson gauge action + asqtad staggered fermion action
//------------------------------------------------------------------
GwilsonFasqtad::GwilsonFasqtad()
{
  cname = "GwilsonFasqtad";
  char *fname = "GwilsonFasqtad()";
  VRB.Func(cname,fname);
}

GwilsonFasqtad::~GwilsonFasqtad()
{
  char *fname = "~GwilsonFasqtad()";
  VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
