#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of ParTransStagTypes class constructor and destructor.

*/

CPS_END_NAMESPACE
#include <util/pt.h>
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
/*!
  \param latt The lattice on which this operation is defined
 */
//------------------------------------------------------------------

static StrOrdType old_str_ord;
ParTransWilsonTypes::ParTransWilsonTypes(Lattice & latt) :
                                   ParTrans(latt)
{
  cname = "ParTransStagTypes";
  char *fname = "ParTransStagTypes(Lattice&)";
  VRB.Func(cname,fname);
  if (lat.StrOrd() != WILSON && lat.StrOrd() != CANONICAL
      && lat.StrOrd() != DWF_4D_EOPREC  && lat.StrOrd() != DWF_4D_EOPREC_EE ){
    old_str_ord = lat.StrOrd();
    lat.Convert(CANONICAL);
  }
  pt_init(lat);
  pt_init_g();
}


//------------------------------------------------------------------
ParTransWilsonTypes::~ParTransWilsonTypes() {
  char *fname = "~ParTransStagTypes()";
  VRB.Func(cname,fname);
  pt_delete_g();
  pt_delete();
  if (lat.StrOrd() != WILSON && lat.StrOrd() != CANONICAL 
      && lat.StrOrd() != DWF_4D_EOPREC  && lat.StrOrd() != DWF_4D_EOPREC_EE ){
    lat.Convert(old_str_ord);
  }
}

CPS_END_NAMESPACE
