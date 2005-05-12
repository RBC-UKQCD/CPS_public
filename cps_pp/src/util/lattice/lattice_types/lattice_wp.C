#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of the GwilsonFp4 class.

*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Wilson gauge action + p4 staggered fermion action
//------------------------------------------------------------------
GwilsonFp4::GwilsonFp4()
	: Fsmear(1)
{
  cname = "GwilsonFp4";
  char *fname = "GwilsonFp4()";
  VRB.Func(cname,fname);
}

GwilsonFp4::~GwilsonFp4()
{
  char *fname = "~GwilsonFp4()";
  VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
