#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of GimprRectFnone class.

  $Id: lattice_rn.C,v 1.4 2004-08-18 11:58:05 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:58:05 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_types/lattice_rn.C,v 1.4 2004-08-18 11:58:05 zs Exp $
//  $Id: lattice_rn.C,v 1.4 2004-08-18 11:58:05 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_types/lattice_rn.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// lattice_wn.C
//
// This class has double inheritance. The virtual
// base class is Lattice with two types of derived classes. One type
// is relevant to the gauge part and has a name that begins with 
// "G". The other type is relevant to the fermion part and has
// a name that begins with "F". The class below inherits from one 
// GimprRect and Fnone.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Wilson gauge action -- no fermion action
//------------------------------------------------------------------
GimprRectFnone::GimprRectFnone()
{
  cname = "GimprRectFnone";
  char *fname = "GimprRectFnone()";
  VRB.Func(cname,fname);

  //???
}

GimprRectFnone::~GimprRectFnone()
{
  char *fname = "~GimprRectFnone()";
  VRB.Func(cname,fname);

  //???
}



CPS_END_NAMESPACE
