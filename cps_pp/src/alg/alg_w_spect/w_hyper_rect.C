#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-01-13 20:39:02 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_w_spect/w_hyper_rect.C,v 1.3 2004-01-13 20:39:02 chulwoo Exp $
//  $Id: w_hyper_rect.C,v 1.3 2004-01-13 20:39:02 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2.10.1  2003/11/06 00:26:06  cwj
//  *** empty log message ***
//
//  Revision 1.1.1.1  2003/11/04 05:04:59  chulwoo
//
//  starting again
//
//
//  Revision 1.2  2003/07/24 16:53:54  zs
//  Addition of documentation via doxygen:
//  doxygen-parsable comment blocks added to many source files;
//  New target in makefile and consequent alterations to configure.in;
//  New directories and files under the doc directory.
//
//  Revision 1.2  2001/06/19 18:11:34  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:00  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: w_hyper_rect.C,v $
//  $Revision: 1.3 $
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
