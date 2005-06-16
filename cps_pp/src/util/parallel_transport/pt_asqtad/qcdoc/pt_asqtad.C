#include <config.h>
#include <stdio.h>
CPS_START_NAMESPACE
/*! \file
  \brief    Definition of ParTransAsqtad class methods for QCDOC.

  $Id: pt_asqtad.C,v 1.10 2005-06-16 14:36:50 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2005-06-16 14:36:50 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_asqtad/qcdoc/pt_asqtad.C,v 1.10 2005-06-16 14:36:50 chulwoo Exp $
//  $Id: pt_asqtad.C,v 1.10 2005-06-16 14:36:50 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_asqtad.C,v $
//  $Revision: 1.10 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_asqtad/qcdoc/pt_asqtad.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// pt_asqtad.C
//
// ParTransAsqtad is derived from the ParTransStagTypes class.
// ParTransAsqtad is the front end for a library that contains
// all Dirac operators associated with Staggered fermions.
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

ParTransAsqtad::ParTransAsqtad(Lattice & latt) :
			 ParTransStagTypes(latt)
{
  cname = "ParTransAsqtad";
  char *fname = "ParTransAsqtad(L&,V*,V*,CgArg*,CnvFrmType)";
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
//  printf("%s:%s end\n",cname,fname);
}


//------------------------------------------------------------------
ParTransAsqtad::~ParTransAsqtad() {
  char *fname = "~ParTransAsqtad()";
  VRB.Func(cname,fname);

//    lat.Convert(CANONICAL);

  //----------------------------------------------------------------
  // Free memory
  //----------------------------------------------------------------
//  pt_delete_g();
//  pt_delete();
}

CPS_END_NAMESPACE
