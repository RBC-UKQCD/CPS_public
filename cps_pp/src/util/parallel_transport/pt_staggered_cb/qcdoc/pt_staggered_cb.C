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
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/vector.h>
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/pt.h>
#include <util/error.h>
#include <util/stag.h>
#include <comms/cbuf.h>
#include <comms/glb.h>
#include <math.h>
CPS_START_NAMESPACE



//------------------------------------------------------------------

ParTransStaggered_cb::ParTransStaggered_cb(Lattice & latt) :
			 ParTransStagTypes(latt)
{
  cname = "ParTransStaggered_cb";
  char *fname = "ParTransStaggered_cb(L&,V*,V*,CgArg*,CnvFrmType)";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Do the necessary conversions
  //----------------------------------------------------------------
//    lat.Convert(STAG);

  //----------------------------------------------------------------
  // Set the node checkerboard size of the fermion field
  //----------------------------------------------------------------
  f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / 2;

//  pt_init(lat.StrOrd(),lat.GaugeField());
//  pt_init_g();

#if 0
  //----------------------------------------------------------------
  // Allocate memory for the temporary fermion vector frm_tmp.
  //----------------------------------------------------------------
  frm_tmp = (Vector *) smalloc(f_size_cb * sizeof(Float));
  if(frm_tmp == 0)
    ERR.Pointer(cname,fname, "frm_tmp");
  VRB.Smalloc(cname,fname, "frm_tmp", 
	      frm_tmp, f_size_cb * sizeof(Float));
#endif
  //printf("%s:%s end\n",cname,fname);
}


//------------------------------------------------------------------
ParTransStaggered_cb::~ParTransStaggered_cb() {
  char *fname = "~ParTransStaggered_cb()";
  VRB.Func(cname,fname);

//    lat.Convert(CANONICAL);

  //----------------------------------------------------------------
  // Free memory
  //----------------------------------------------------------------
//  pt_delete_g();
//  pt_delete();
}

CPS_END_NAMESPACE
