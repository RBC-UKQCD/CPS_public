#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of ParTransStagTypes class constructor and destructor.

*/
//--------------------------------------------------------------------
//
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/pt.h>
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/dirac_op.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
/*!
  \param latt The lattice on which this operation is defined
 */
//------------------------------------------------------------------

static StrOrdType old_str_ord;
ParTransGauge::ParTransGauge(Lattice & latt) :
                                   ParTrans(latt)
{
  cname = "ParTransGauge";
  char *fname = "ParTransGauge(Lattice&)";
  VRB.Func(cname,fname);
  if (lat.StrOrd() != WILSON && lat.StrOrd() != CANONICAL &&
      lat.StrOrd() != DWF_4D_EOPREC  && lat.StrOrd() != DWF_4D_EOPREC_EE ){
    old_str_ord = lat.StrOrd();
    lat.Convert(CANONICAL);
  }
  if(GJP.Gparity() && DiracOp::scope_lock ==0)
    BondCond(lat, gauge_field); //remove fermionic boundary conditions on gauge links (they are applied in the ParTrans constructor)

  pt_init(lat);
  pt_init_g();
}


//------------------------------------------------------------------
ParTransGauge::~ParTransGauge() {
  char *fname = "~ParTransGauge()";
  VRB.Func(cname,fname);
  pt_delete_g();
  pt_delete();
  if (lat.StrOrd() != WILSON && lat.StrOrd() != CANONICAL &&
      lat.StrOrd() != DWF_4D_EOPREC && lat.StrOrd() != DWF_4D_EOPREC_EE  ){

  if(GJP.Gparity() && DiracOp::scope_lock ==0)
    BondCond(lat, gauge_field); //remove fermionic boundary conditions on gauge links (they are applied in the ParTrans constructor)

    lat.Convert(old_str_ord);
  }
}

CPS_END_NAMESPACE
