#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of the GwilsonFwilson class.

  $Id: lattice_ww.C,v 1.5 2006/06/11 05:35:06 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2006/06/11 05:35:06 $
//  $Header: /space/cvs/cps/cps++/src/util/lattice/lattice_types/lattice_ww.C,v 1.5 2006/06/11 05:35:06 chulwoo Exp $
//  $Id: lattice_ww.C,v 1.5 2006/06/11 05:35:06 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $Revision: 1.5 $
//  $Source: /space/cvs/cps/cps++/src/util/lattice/lattice_types/lattice_ww.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// lattice_ww.C
//
// This class has double inheritance. The virtual
// base class is Lattice with two types of derived classes. One type
// is relevant to the gauge part and has a name that begins with 
// "G". The other type is relevant to the fermion part and has
// a name that begins with "F". The class below inherits from one 
// Gwilson and Fwilson.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Wilson gauge action -- wilson fermion action
//------------------------------------------------------------------
GwilsonFwilson::GwilsonFwilson()
:Gwilson(),Fwilson()
{
  cname = "GwilsonFwilson";
  char *fname = "GwilsonFwilson()";
  VRB.Func(cname,fname);

  //???
}

GwilsonFwilson::~GwilsonFwilson()
{
  char *fname = "~GwilsonFwilson()";
  VRB.Func(cname,fname);

  //???
}



CPS_END_NAMESPACE
