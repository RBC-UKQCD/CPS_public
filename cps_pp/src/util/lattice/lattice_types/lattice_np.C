#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of GnoneFp4 class.

  $Id: lattice_np.C,v 1.2 2004/12/08 20:52:14 chulwoo Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// No gauge action + P4 improved staggered fermion action
//------------------------------------------------------------------
GnoneFp4::GnoneFp4()
	: Gnone(),Fp4(),Fsmear(1)
{
  cname = "GnoneFp4";
  char *fname = "GnoneFp4()";
  VRB.Func(cname,fname);
}

GnoneFp4::~GnoneFp4()
{
  char *fname = "~GnoneFp4()";
  VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
