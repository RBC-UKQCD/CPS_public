#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of GimprOLSymFnone

  class methods.

  $Id: lattice_sn.C,v 1.3 2004-06-04 21:14:13 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:13 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_impr_OLSym/lattice_sn.C,v 1.3 2004-06-04 21:14:13 chulwoo Exp $
//  $Id: lattice_sn.C,v 1.3 2004-06-04 21:14:13 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_impr_OLSym/lattice_sn.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// lattice_sn.C
//
// This class has double inheritance. The virtual
// base class is Lattice with two types of derived classes. One type
// is relevant to the gauge part and has a name that begins with 
// "G". The other type is relevant to the fermion part and has
// a name that begins with "F". The class below inherits from one 
// GimprOLSym and Fnone.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Wilson gauge action -- no fermion action
//------------------------------------------------------------------
GimprOLSymFnone::GimprOLSymFnone()
{
  cname = "GimprOLSymFnone";
  char *fname = "GimprOLSymFnone()";
  VRB.Func(cname,fname);

  //???
}

GimprOLSymFnone::~GimprOLSymFnone()
{
  char *fname = "~GimprOLSymFnone()";
  VRB.Func(cname,fname);

  //???
}



CPS_END_NAMESPACE
