#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Implementation of GpowerPlaqFstagAsqtad class.

  $Id: lattice_pa.C,v 1.2 2003-10-23 13:38:59 zs Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// PowerPlaq gauge action + asqtad staggered fermion action
//------------------------------------------------------------------
GpowerPlaqFstagAsqtad::GpowerPlaqFstagAsqtad()
{
  cname = "GpowerPlaqFstagAsqtad"; 
  char *fname = "GpowerPlaqFstagAsqtad()";
  VRB.Func(cname,fname);
}

GpowerPlaqFstagAsqtad::~GpowerPlaqFstagAsqtad()
{
    char *fname = "~GpowerPlaqFstagAsqtad()";
    VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
