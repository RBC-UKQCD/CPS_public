#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of GnoneFstag class.

  $Id: lattice_ns.C,v 1.3 2004-06-04 21:14:13 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:13 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_types/lattice_ns.C,v 1.3 2004-06-04 21:14:13 chulwoo Exp $
//  $Id: lattice_ns.C,v 1.3 2004-06-04 21:14:13 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_types/lattice_ns.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// lattice_ns.C
//
// This class has double inheritance. The virtual
// base class is Lattice with two types of derived classes. One type
// is relevant to the gauge part and has a name that begins with 
// "G". The other type is relevant to the fermion part and has
// a name that begins with "F". The class below inherits from one 
// Gnone and Fstag.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// No gauge action -- staggered fermion action
//------------------------------------------------------------------
GnoneFstag::GnoneFstag()
{
  cname = "GnoneFstag";
  char *fname = "GnoneFstag()";
  VRB.Func(cname,fname);

  //???
}

GnoneFstag::~GnoneFstag()
{
  char *fname = "~GnoneFstag()";
  VRB.Func(cname,fname);

  //???
}



CPS_END_NAMESPACE
