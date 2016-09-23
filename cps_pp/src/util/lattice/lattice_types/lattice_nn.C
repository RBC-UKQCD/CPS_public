#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of GnoneFnone class.

*/
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
