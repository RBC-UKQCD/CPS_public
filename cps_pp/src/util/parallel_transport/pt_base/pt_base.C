#include <config.h>
#include <util/qcdio.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of ParTrans class methods.
  
  $Id: pt_base.C,v 1.10 2004/08/18 11:58:06 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: ckelly $
//  $Date: 2004/08/18 11:58:06 $
//  $Header: /space/cvs/cps/cps++/src/util/parallel_transport/pt_base/pt_base.C,v 1.10 2004/08/18 11:58:06 zs Exp $
//  $Id: pt_base.C,v 1.10 2004/08/18 11:58:06 zs Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: pt_base.C,v $
//  $Revision: 1.10.470.1 $
//  $Source: /space/cvs/cps/cps++/src/util/parallel_transport/pt_base/pt_base.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


CPS_END_NAMESPACE
#include <util/pt.h>
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/dirac_op.h>
#include <util/error.h>
#include <util/gjp.h>
#include <comms/cbuf.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// Initialize static variables
//------------------------------------------------------------------
int ParTrans::scope_lock = 0;
double ParTrans::PTflops = 0;
int ParTrans:: bc[4] = {0,0,0,0};	// boundary condition on this node
int ParTrans::node_sites[5] = {0,0,0,0,0};

void ParTrans::BondCond(Lattice& lat, Matrix *u_base)
{
  int uconj_offset = 4*GJP.VolNodeSites();
  for(int u = 0; u < 4; ++u) {
    if(bc[u]) {
      int x[4];
      for(x[3] = 0; x[3] < node_sites[3]; ++x[3]) {
	for(x[2] = 0; x[2] < node_sites[2]; ++x[2]) {
	  for(x[1] = 0; x[1] < node_sites[1]; ++x[1]) {
	    for(x[0] = 0; x[0] < node_sites[0]; ++x[0]) {
	      if(x[u]==node_sites[u]-1) {
		int site_off = lat.GsiteOffset(x);
	        Matrix *m = u_base+site_off+u;
	        *m *= -1;

		if(GJP.Gparity()){
		  //for example, APBC in time direction
		  Matrix *m = u_base+uconj_offset+site_off+u;
		  *m *= -1;
		}	

	      }
	    }
          }
	}
      }
    }

    if(GJP.Bc(u)==BND_CND_GPARITY && GJP.NodeCoor(u) == GJP.Nodes(u)-1){
      int x[4];
      //also put minus signs on outwards facing links of U* field
      for(x[0] = 0; x[0] < node_sites[0]; ++x[0]) {
        for(x[1] = 0; x[1] < node_sites[1]; ++x[1]) {
          for(x[2] = 0; x[2] < node_sites[2]; ++x[2]) {
            for(x[3] = 0; x[3] < node_sites[3]; ++x[3]) {
	      if(x[u]==node_sites[u]-1) {
		int site_off = lat.GsiteOffset(x);
	        Matrix *m = u_base+uconj_offset+site_off+u;
	        *m *= -1;
	      }
	    }
          }
	}
      }
    }


  }
}


//------------------------------------------------------------------
/*!
  Only one instance of this class is allowed to be in existence at
  any time.
  \param latt The lattice containing the gauge field on which this operation
  is defined
 */
//------------------------------------------------------------------
ParTrans::ParTrans(Lattice & latt) :
    lat(latt)
{
  cname = "ParTrans";
  char *fname = "ParTrans(Lattice&)";
  VRB.Clock(cname,fname,"Just entered\n");
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Check and set the scope_lock
  //----------------------------------------------------------------
  if(scope_lock != 0){
    ERR.General(cname,fname,
		"Only one ParTrans object is allowed to be on scope\n");
  }
  scope_lock = 1;

  //----------------------------------------------------------------
  // Set the gauge field pointer
  //----------------------------------------------------------------
  gauge_field = lat.GaugeField();

  //----------------------------------------------------------------
  // Initialize boundary condition on this node
  //----------------------------------------------------------------
  bc[0] = 0;
  bc[1] = 0;
  bc[2] = 0;
  bc[3] = 0;
  if(GJP.Xbc() == BND_CND_APRD)
  	bc[0] = GJP.XnodeCoor() == (GJP.Xnodes()-1) ? 1 : 0;
  if(GJP.Ybc() == BND_CND_APRD)
  	bc[1] = GJP.YnodeCoor() == (GJP.Ynodes()-1) ? 1 : 0;
  if(GJP.Zbc() == BND_CND_APRD)
  	bc[2] = GJP.ZnodeCoor() == (GJP.Znodes()-1) ? 1 : 0;
  if(GJP.Tbc() == BND_CND_APRD)
  	bc[3] = GJP.TnodeCoor() == (GJP.Tnodes()-1) ? 1 : 0;

  node_sites[0] = GJP.XnodeSites();
  node_sites[1] = GJP.YnodeSites();
  node_sites[2] = GJP.ZnodeSites();
  node_sites[3] = GJP.TnodeSites();
  node_sites[4] = GJP.SnodeSites();

  //----------------------------------------------------------------
  // turn on the boundary condition, if not inside a dirac operator
  //----------------------------------------------------------------
  if(DiracOp::scope_lock || lat.BcApplied() ) bc_already_applied=1;
  else bc_already_applied=0;
  if (!bc_already_applied) lat.BondCond();

//  lat.Convert(STAG);


  // Added in by Ping for anisotropic lattices
  //------------------------------------------------------------------
  // No data manipulations if the scaling factor is 1.0
  // [built into lat.MltFloat].
  {
    Float factor = GJP.XiVXi() / GJP.XiV();
    lat.MltFloat(factor, GJP.XiDir());
  }
}


//------------------------------------------------------------------

ParTrans::~ParTrans() {
  char *fname = "~ParTrans()";

  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // turn off the boundary condition
  //----------------------------------------------------------------
//  if(DiracOp::scope_lock ==0)
  if (!bc_already_applied)
  lat.BondCond();

  VRB.Clock(cname,fname,"Exiting\n");


  //----------------------------------------------------------------
  // Release the scope_lock
  //----------------------------------------------------------------
  scope_lock = 0;



  // Added in by Ping for anisotropic lattices
  //------------------------------------------------------------------
  // No data manipulations if the scaling factor is 1.0
  // [built into lat.MltFloat].
  {
    Float factor = GJP.XiV() / GJP.XiVXi();
    lat.MltFloat(factor, GJP.XiDir());
  }
}

CPS_END_NAMESPACE
