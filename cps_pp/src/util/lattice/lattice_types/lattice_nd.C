#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of GnoneFdwf class.

  $Id: lattice_nd.C,v 1.3 2004-06-04 21:14:13 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:13 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_types/lattice_nd.C,v 1.3 2004-06-04 21:14:13 chulwoo Exp $
//  $Id: lattice_nd.C,v 1.3 2004-06-04 21:14:13 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_types/lattice_nd.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// lattice_nd.C
//
// This class has double inheritance. The virtual
// base class is Lattice with two types of derived classes. One type
// is relevant to the gauge part and has a name that begins with 
// "G". The other type is relevant to the fermion part and has
// a name that begins with "F". The class below inherits from one 
// Gnone and Fdwf.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// No gauge action -- domain wall fermion action
//------------------------------------------------------------------
GnoneFdwf::GnoneFdwf()
{
  cname = "GnoneFdwf";
  char *fname = "GnoneFdwf()";
  VRB.Func(cname,fname);

  //???
}

GnoneFdwf::~GnoneFdwf()
{
  char *fname = "~GnoneFdwf()";
  VRB.Func(cname,fname);

  //???
}



CPS_END_NAMESPACE
