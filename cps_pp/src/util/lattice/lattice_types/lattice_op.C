#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of GimprOLSymFp4 class.

*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// One Loop Symanzik improved gauge action with P4 staggered fermion
//------------------------------------------------------------------
GimprOLSymFp4::GimprOLSymFp4()
	: Fsmear(1)
{
  cname = "GimprOLSymFp4"; 
  char *fname = "GimprOLSymFp4()";
  VRB.Func(cname,fname);
}

GimprOLSymFp4::~GimprOLSymFp4()
{
    char *fname = "~GimprOLSymFp4()";
    VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
