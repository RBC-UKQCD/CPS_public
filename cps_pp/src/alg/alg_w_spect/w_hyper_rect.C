#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:00 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_w_spect/w_hyper_rect.C,v 1.5 2004-06-04 21:14:00 chulwoo Exp $
//  $Id: w_hyper_rect.C,v 1.5 2004-06-04 21:14:00 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: w_hyper_rect.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_w_spect/w_hyper_rect.C,v $
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
