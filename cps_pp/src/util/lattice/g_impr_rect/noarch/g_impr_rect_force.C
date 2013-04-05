#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of GimprRect class.

  $Id: g_impr_rect_force.C,v 1.4 2013-04-05 17:51:14 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:51:14 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_impr_rect/noarch/g_impr_rect_force.C,v 1.4 2013-04-05 17:51:14 chulwoo Exp $
//  $Id: g_impr_rect_force.C,v 1.4 2013-04-05 17:51:14 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_impr_rect/noarch/g_impr_rect_force.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdlib.h>
#include <math.h>
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/gw_hb.h>
#include <util/time_cps.h>
//#include <comms/nga_reg.h>
#include <comms/glb.h>
#include <comms/cbuf.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------------------
enum { MATRIX_SIZE = 18 };


//------------------------------------------------------------------------------
// static variables used only inside this file
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//  CRAM temp buffer
//------------------------------------------------------------------------------
static Matrix mt0;
static Matrix mt1;
static Matrix mt2;
//static Matrix mt3;
static Matrix *mp0 = &mt0;		// ihdot
static Matrix *mp1 = &mt1;
static Matrix *mp2 = &mt2;

#define PROFILE
//------------------------------------------------------------------------------
// EvolveMomGforce(Matrix *mom, Float dt):
// It evolves the canonical momentum mom by dt
// using the pure gauge force.
//------------------------------------------------------------------------------


ForceArg GimprRect::EvolveMomGforce(Matrix *mom, Float dt){
  char *fname = "EvolveMomGforce(M*,F)";
  VRB.Func(cname,fname);

  Float L1=0.0;
  Float L2=0.0;
  Float Linf=0.0;

#ifdef PROFILE
  Float time = -dclock();
  ForceFlops = 0;
#endif
  
  setCbufCntrlReg(4, CBUF_MODE4);

  int x[4];
  
  for(x[0] = 0; x[0] < GJP.XnodeSites(); ++x[0])
  for(x[1] = 0; x[1] < GJP.YnodeSites(); ++x[1])
  for(x[2] = 0; x[2] < GJP.ZnodeSites(); ++x[2])
  for(x[3] = 0; x[3] < GJP.TnodeSites(); ++x[3]) {

    int uoff = GsiteOffset(x);

    for (int mu = 0; mu < 4; ++mu) {
      GforceSite(*mp0, x, mu);

      IFloat *ihp = (IFloat *)(mom+uoff+mu);
      IFloat *dotp = (IFloat *)mp0;
      fTimesV1PlusV2(ihp, dt, dotp, ihp, 18);
      Float norm = ((Matrix*)dotp)->norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);
   }
  }
  ForceFlops +=GJP.VolNodeSites()*4*18*2;
#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,ForceFlops,time);
#endif

  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();

  VRB.FuncEnd(cname,fname);
  return ForceArg(dt*L1, dt*sqrt(L2), dt*Linf);
}

CPS_END_NAMESPACE
