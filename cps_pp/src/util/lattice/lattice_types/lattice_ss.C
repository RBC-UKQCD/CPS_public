#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of GimprOLSymFstag class methods.

  $Id: lattice_ss.C,v 1.2 2005-03-07 00:33:41 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2005-03-07 00:33:41 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_types/lattice_ss.C,v 1.2 2005-03-07 00:33:41 chulwoo Exp $
//  $Id: lattice_ss.C,v 1.2 2005-03-07 00:33:41 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_types/lattice_ss.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// lattice_ss.C
//
// This class has double inheritance. The virtual
// base class is Lattice with two types of derived classes. One type
// is relevant to the gauge part and has a name that begins with 
// "G". The other type is relevant to the fermion part and has
// a name that begins with "F". The class below inherits from one 
// GimprOLSym and Fstag.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Wilson gauge action -- staggered fermion action
//------------------------------------------------------------------
GimprOLSymFstag::GimprOLSymFstag()
{
  cname = "GimprOLSymFstag";
  char *fname = "GimprOLSymFstag()";
  VRB.Func(cname,fname);

  //???
}

GimprOLSymFstag::~GimprOLSymFstag()
{
  char *fname = "~GimprOLSymFstag()";
  VRB.Func(cname,fname);

  //???
}



CPS_END_NAMESPACE
