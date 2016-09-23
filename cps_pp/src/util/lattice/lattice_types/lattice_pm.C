#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of GpowerPlaqFdwf class.

  $Id: lattice_pm.C,v 1.2 2011/03/21 21:04:50 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2011/03/21 21:04:50 $
//  $Header: /space/cvs/cps/cps++/src/util/lattice/lattice_types/lattice_pm.C,v 1.2 2011/03/21 21:04:50 chulwoo Exp $
//  $Id: lattice_pm.C,v 1.2 2011/03/21 21:04:50 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $Revision: 1.2 $
//  $Source: /space/cvs/cps/cps++/src/util/lattice/lattice_types/lattice_pm.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// lattice_pd.C
//
// This class has double inheritance. The virtual
// base class is Lattice with two types of derived classes. One type
// is relevant to the gauge part and has a name that begins with 
// "G". The other type is relevant to the fermion part and has
// a name that begins with "F". The class below inherits from one 
// GpowerPlaq and Fdwf.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// PowerPlaq gauge action -- domain wall fermion action
//------------------------------------------------------------------
GpowerPlaqFmdwf::GpowerPlaqFmdwf()
{
  cname = "GpowerPlaqFmdwf";
  char *fname = "GpowerPlaqFmdwf()";
  VRB.Func(cname,fname);

  //???
}

GpowerPlaqFmdwf::~GpowerPlaqFmdwf()
{
  char *fname = "~GpowerPlaqFmdwf()";
  VRB.Func(cname,fname);

  //???
}



CPS_END_NAMESPACE
