#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:13:59 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_rect/alg_rect.C,v 1.5 2004-06-04 21:13:59 chulwoo Exp $
//  $Id: alg_rect.C,v 1.5 2004-06-04 21:13:59 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_rect.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_rect/alg_rect.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// alg_rect.C
//
// AlgRect is derived from Alg and it measures the average
// value of the rectangle. 
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdio.h>
#include <alg/alg_rect.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
#include <comms/glb.h>
CPS_START_NAMESPACE

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))

//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
AlgRect::AlgRect(Lattice& latt, 
             CommonArg *c_arg,
             NoArg *arg) : 
             Alg(latt, c_arg) 
{
  cname = "AlgRect";
  char *fname = "AlgRect(L&,CommonArg*,NoArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_rect_arg = arg;

  // Calculate normalization factor
  //----------------------------------------------------------------
  int total_sites = GJP.VolNodeSites() * GJP.Xnodes() *
                    GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
  norm_fac = 1.0 / ( 12.0 * Float(total_sites) );

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgRect::~AlgRect() {
  char *fname = "~AlgRect()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
//
//------------------------------------------------------------------
void AlgRect::run()
{
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  char *fname = "run()";
  VRB.Func(cname, fname);

  //----------------------------------------------------------------
  // Check if anisotropy is present and exit 
  //----------------------------------------------------------------
  if(GJP.XiBare() != 1 ) {
    ERR.General(cname,fname,
    "XiBare=%g : Not implemented for anisotropy\n", GJP.XiBare());
  }

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  Float r_sum    =  0.0 ;
  Float r_sq_sum =  0.0 ;
  Float r_min    =  3.0 ;
  Float r_max    = -3.0 ;
  int x[4];
 
  for (x[0] = 0; x[0] < GJP.XnodeSites(); ++x[0]) 
  for (x[1] = 0; x[1] < GJP.YnodeSites(); ++x[1]) 
  for (x[2] = 0; x[2] < GJP.ZnodeSites(); ++x[2]) 
  for (x[3] = 0; x[3] < GJP.TnodeSites(); ++x[3]) 
  for (int mu = 0; mu < 4; ++mu)
  for (int nu = 0; nu < 4; ++nu)
  if ( mu != nu )  {

     Float tmp_flt = lat.ReTrRect(x, mu, nu) ;

     r_sum    += tmp_flt ;
     r_sq_sum += tmp_flt * tmp_flt ;
     r_min     = min(r_min, tmp_flt) ;
     r_max     = max(r_max, tmp_flt) ;
  }

  glb_sum(&r_sum) ;
  glb_sum(&r_sq_sum) ;
  glb_min(&r_min) ;
  glb_max(&r_max) ;

  Float one_third = 1.0 / 3.0 ;

  r_sum    *= one_third ;
  r_sq_sum *= one_third * one_third ;
  r_min    *= one_third ;
  r_max    *= one_third ;

  Float r_var = norm_fac * (1.0 + norm_fac) \
                * (r_sq_sum - norm_fac*r_sum*r_sum) ;

  r_sum *= norm_fac ;

  x[0]=0 ; x[1]=0 ; x[2]=0 ; x[3]=0 ;

  // Print out results
  //----------------------------------------------------------------
  if(common_arg->results != 0){
    FILE *fp;
    if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    fprintf(fp, "%e %e %e %e %e\n", IFloat(r_sum), IFloat(r_var),
      IFloat(r_max), IFloat(r_min),
      IFloat(one_third*lat.ReTrRect(x, 0, 1)) ) ;
    fclose(fp);
  }

}

CPS_END_NAMESPACE
