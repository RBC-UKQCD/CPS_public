#include<config.h>
CPS_START_NAMESPACE

//-------------------------------------------------------------------
//
// lattice_tn.C
//
// This class has double inheritance. The virtual
// base class is Lattice with two types of derived classes. One type
// is relevant to the gauge part and has a name that begins with 
// "G". The other type is relevant to the fermion part and has
// a name that begins with "F". The class below inherits from one 
// GtadpoleRect and Fnone.
//
//------------------------------------------------------------------


CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Tadpole Rect gauge action -- no fermion action
//------------------------------------------------------------------
GtadpoleRectFnone::GtadpoleRectFnone()
{
  cname = "GtadpoleRectFnone";
  char *fname = "GtadpoleRectFnone()";
  VRB.Func(cname,fname);
}

GtadpoleRectFnone::~GtadpoleRectFnone()
{
  char *fname = "~GtadpoleRectFnone()";
  VRB.Func(cname,fname);
}



CPS_END_NAMESPACE
