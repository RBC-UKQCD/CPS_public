#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of GpowerPlaqFasqtad class.

  $Id: lattice_pa.C,v 1.3 2003-10-30 05:39:12 cwj Exp $
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
