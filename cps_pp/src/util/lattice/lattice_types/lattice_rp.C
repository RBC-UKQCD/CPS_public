#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of GimprRectFp4 class.

  $Id: lattice_rp.C,v 1.2 2005-05-12 20:47:01 chulwoo Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Wilson gauge action + p4staggered fermion action
//------------------------------------------------------------------
GimprRectFp4::GimprRectFp4()
	: Fsmear(1)
{
  cname = "GimprRectFp4";
  char *fname = "GimprRectFp4()";
  VRB.Func(cname,fname);
}

GimprRectFp4::~GimprRectFp4()
{
  char *fname = "~GimprRectFp4()";
  VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
