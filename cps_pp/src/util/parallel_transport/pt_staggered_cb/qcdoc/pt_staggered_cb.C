#include <config.h>
#include <stdio.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
//
// pt_staggered_cb.C
//
// ParTransStaggered is derived from the ParTransStagTypes class.
// ParTransStaggered implements a parallel transporter for staggered
// actions where the checkerboarded storage scheme is used.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/pt.h>
#include <util/gjp.h>
CPS_START_NAMESPACE



//------------------------------------------------------------------

/*!
  \param latt The lattice object containing the gauge field on which the
  parallel tranport operates.
  \post The gauge field storage order is converted to STAG order if it is not
  already.
*/

//All the necessary initializations are performed in ParTransStagTypes
//Conversion of storage order is also done by ParTransStagtypes
static StrOrdType old_str_ord2;
ParTransStaggered_cb::ParTransStaggered_cb(Lattice & latt) :
			 ParTransStagTypes(latt)
{
  cname = "ParTransStaggered_cb";
  char *fname = "ParTransStaggered_cb(L&,V*,V*,CgArg*,CnvFrmType)";
  VRB.Func(cname,fname);
  old_str_ord2 = lat.StrOrd();
  if (lat.StrOrd() != STAG_BLOCK){
    lat.Convert(STAG_BLOCK);
  }
}

/*!
  \post The gauge field storage order is converted back to CANONICAL order
  if that was how it was when this object was created.
 */

//------------------------------------------------------------------
//All necessary memory de-allocation is performed by the
//destructor for ParTransStagTypes

ParTransStaggered_cb::~ParTransStaggered_cb() {
  char *fname = "~ParTransStaggered_cb()";
  VRB.Func(cname,fname);
  if ( old_str_ord2 !=STAG_BLOCK){
    lat.Convert(old_str_ord2);
  }
}

CPS_END_NAMESPACE
