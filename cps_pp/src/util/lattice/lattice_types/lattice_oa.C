#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of GimprOLSymFasqtad class.

  $Id: lattice_oa.C,v 1.4 2004-08-18 11:58:05 zs Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// One Loop Symanzik improved gauge action with Asqtad staggered fermion
//------------------------------------------------------------------
GimprOLSymFasqtad::GimprOLSymFasqtad()
	: Fsmear(3)
{
  cname = "GimprOLSymFasqtad"; 
  char *fname = "GimprOLSymFasqtad()";
  VRB.Func(cname,fname);
}

GimprOLSymFasqtad::~GimprOLSymFasqtad()
{
    char *fname = "~GimprOLSymFasqtad()";
    VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
