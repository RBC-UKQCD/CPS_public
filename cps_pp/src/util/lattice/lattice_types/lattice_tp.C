#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of GtadpoleRectFp4 class.

  $Id: lattice_tp.C,v 1.2 2005-05-12 20:47:01 chulwoo Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Wilson gauge action + p4staggered fermion action
//------------------------------------------------------------------
GtadpoleRectFp4::GtadpoleRectFp4()
	: Fsmear(1)
{
  cname = "GtadpoleRectFp4";
  char *fname = "GtadpoleRectFp4()";
  VRB.Func(cname,fname);
}

GtadpoleRectFp4::~GtadpoleRectFp4()
{
  char *fname = "~GtadpoleRectFp4()";
  VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
