#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of the GwilsonFasqtad class.

  $Id: lattice_wa.C,v 1.4 2004-08-09 07:47:24 chulwoo Exp $
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
	: Fsmear(3)
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
