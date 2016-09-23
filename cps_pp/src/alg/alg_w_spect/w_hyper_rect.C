#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004/08/18 11:57:40 $
//  $Header: /space/cvs/cps/cps++/src/alg/alg_w_spect/w_hyper_rect.C,v 1.6 2004/08/18 11:57:40 zs Exp $
//  $Id: w_hyper_rect.C,v 1.6 2004/08/18 11:57:40 zs Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: w_hyper_rect.C,v $
//  $Revision: 1.6 $
//  $Source: /space/cvs/cps/cps++/src/alg/alg_w_spect/w_hyper_rect.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <alg/w_all.h>
#include <util/error.h>
CPS_START_NAMESPACE



char * WspectHyperRectangle::d_class_name = "WspectHyperRectangle";


//---------------------------------------------------------------------------
// WspectHyperRectangle::CTOR
//---------------------------------------------------------------------------
WspectHyperRectangle::WspectHyperRectangle(int dir, int glb_coord)
  : d_dir(dir)
{
  // Check for invalid arguments
  //-------------------------------------------------------------------------
  if (dir < 0       || dir >= LORENTZs   ||           
      glb_coord < 0 || glb_coord >= glb_sites[d_dir]) {
    ERR.General(d_class_name, 	ctor_str, 
		"%s %d %d\n", 	out_range_str, dir, glb_coord);
  }

  for (int l = 0; l < LORENTZs; ++l) {
    d_lcl_min[l] = 0;
    d_lcl_max[l] = lcl_sites[l] - 1;
    d_glb_min[l] = lcl2glb_offset[l];
    d_glb_max[l] = d_lcl_max[l] + lcl2glb_offset[l];    
  }

  d_is_on_node = (glb_coord >= d_glb_min[d_dir] && 
		  glb_coord <= d_glb_max[d_dir]);

  d_glb_min[d_dir] = d_glb_max[d_dir] = glb_coord;
  d_lcl_min[d_dir] = d_lcl_max[d_dir] = glb_coord - lcl2glb_offset[d_dir];
}




CPS_END_NAMESPACE
