#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of the GwilsonFstag class.

  $Id: lattice_ws.C,v 1.3 2004-06-04 21:14:14 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:14 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_types/lattice_ws.C,v 1.3 2004-06-04 21:14:14 chulwoo Exp $
//  $Id: lattice_ws.C,v 1.3 2004-06-04 21:14:14 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_types/lattice_ws.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// lattice_ws.C
//
// This class has double inheritance. The virtual
// base class is Lattice with two types of derived classes. One type
// is relevant to the gauge part and has a name that begins with 
// "G". The other type is relevant to the fermion part and has
// a name that begins with "F". The class below inherits from one 
// Gwilson and Fstag.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Wilson gauge action -- staggered fermion action
//------------------------------------------------------------------
GwilsonFstag::GwilsonFstag()
{
  cname = "GwilsonFstag";
  char *fname = "GwilsonFstag()";
  VRB.Func(cname,fname);

  //???
}

GwilsonFstag::~GwilsonFstag()
{
  char *fname = "~GwilsonFstag()";
  VRB.Func(cname,fname);

  //???
}



CPS_END_NAMESPACE
