#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of GpowerPlaqFasqtad class.

  $Id: lattice_pa.C,v 1.4 2004-08-09 07:47:24 chulwoo Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// PowerPlaq gauge action + asqtad staggered fermion action
//------------------------------------------------------------------
GpowerPlaqFasqtad::GpowerPlaqFasqtad()
	: Fsmear(3)
{
  cname = "GpowerPlaqFasqtad"; 
  char *fname = "GpowerPlaqFasqtad()";
  VRB.Func(cname,fname);
}

GpowerPlaqFasqtad::~GpowerPlaqFasqtad()
{
    char *fname = "~GpowerPlaqFasqtad()";
    VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
