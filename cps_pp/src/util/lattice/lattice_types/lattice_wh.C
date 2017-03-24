#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of the GwilsonFhisq class.

*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Wilson gauge action + hisq staggered fermion action
//------------------------------------------------------------------
GwilsonFhisq::GwilsonFhisq()
	: Fsmear(1)
{
  cname = "GwilsonFhisq";
  char *fname = "GwilsonFhisq()";
  VRB.Func(cname,fname);
}

GwilsonFhisq::~GwilsonFhisq()
{
  char *fname = "~GwilsonFhisq:()";
  VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
