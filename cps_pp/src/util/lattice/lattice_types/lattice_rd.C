#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of GimprRectFdwf class.

  $Id: lattice_rd.C,v 1.4 2004-08-18 11:58:05 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:58:05 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_types/lattice_rd.C,v 1.4 2004-08-18 11:58:05 zs Exp $
//  $Id: lattice_rd.C,v 1.4 2004-08-18 11:58:05 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_types/lattice_rd.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// lattice_wd.C
//
// This class has double inheritance. The virtual
// base class is Lattice with two types of derived classes. One type
// is relevant to the gauge part and has a name that begins with 
// "G". The other type is relevant to the fermion part and has
// a name that begins with "F". The class below inherits from one 
// GimprRect and Fdwf.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Wilson gauge action -- domain wall fermion action
//------------------------------------------------------------------
GimprRectFdwf::GimprRectFdwf()
{
  cname = "GimprRectFdwf";
  char *fname = "GimprRectFdwf()";
  VRB.Func(cname,fname);

  //???
}

GimprRectFdwf::~GimprRectFdwf()
{
  char *fname = "~GimprRectFdwf()";
  VRB.Func(cname,fname);

  //???
}



CPS_END_NAMESPACE
