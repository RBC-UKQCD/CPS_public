#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of GnoneFnone class.

  $Id: lattice_nn.C,v 1.5 2012-08-10 14:05:33 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-08-10 14:05:33 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_types/lattice_nn.C,v 1.5 2012-08-10 14:05:33 chulwoo Exp $
//  $Id: lattice_nn.C,v 1.5 2012-08-10 14:05:33 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_types/lattice_nn.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// lattice_nn.C
//
// This class has double inheritance. The virtual
// base class is Lattice with two types of derived classes. One type
// is relevant to the gauge part and has a name that begins with 
// "G". The other type is relevant to the fermion part and has
// a name that begins with "F". The class below inherits from one 
// Gnone and Fnone.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// No gauge action -- no fermion action
//------------------------------------------------------------------
GnoneFnone::GnoneFnone()
{
  cname = "GnoneFnone";
  char *fname = "GnoneFnone()";
  VRB.Func(cname,fname);

  //???
}

GnoneFnone::~GnoneFnone()
{
  char *fname = "~GnoneFnone()";
  VRB.Func(cname,fname);

  //???
}



CPS_END_NAMESPACE
