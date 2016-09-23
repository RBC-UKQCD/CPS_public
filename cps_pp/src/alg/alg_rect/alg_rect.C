#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2005/06/23 18:01:01 $
//  $Header: /space/cvs/cps/cps++/src/alg/alg_rect/alg_rect.C,v 1.8 2005/06/23 18:01:01 chulwoo Exp $
//  $Id: alg_rect.C,v 1.8 2005/06/23 18:01:01 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: alg_rect.C,v $
//  $Revision: 1.8 $
//  $Source: /space/cvs/cps/cps++/src/alg/alg_rect/alg_rect.C,v $
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
#include <util/qcdio.h>
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
//! Performs the computation.
/*!
  The real trace of the rectangle is averaged over the lattice and normalised
  to unity.
  \post The following data are written to the file specified in the
  common_arg structure:
  -# mean rectangle
  -# variance of mean rectangle
  -# maximum rectangle
  -# minimum rectangle
  -# mean temporal 1x2 rectangle (the short axis is the time direction)
  -# mean temporal 2x1 rectangle (the long axis is the time direction)
  -# mean spatial rectangle
*/
//------------------------------------------------------------------
void AlgRect::run()
{
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  char *fname = "run()";
  VRB.Func(cname, fname);

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  // Modified by Schmidt for anisotropic lattices
  //----------------------------------------------------------------
  if (GJP.XiBare() != 1.0 || GJP.XiV() != 1 || GJP.XiVXi() != 1   ) {
    Float ave_time1  = lat.AveReTrRectXi1();
    Float ave_time2  = lat.AveReTrRectXi2();
    Float ave_space = lat.AveReTrRectNoXi();
    Float ave_both  = ((ave_time1+ave_time2)/2 + ave_space) / 2;  
 
    if (common_arg->results != 0){
      FILE *fp;
      if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
        ERR.FileA(cname,fname, (char *)common_arg->results);
      }
      Fprintf(fp, "%0.16e %0.16e %0.16e %0.16e \n", 
              (IFloat)ave_space, (IFloat)ave_time1, (IFloat)ave_time2, (IFloat)ave_both);  
      Fclose(fp);  
    }
    return;
  }


  Float r_temporal1, r_temporal2, r_spacial;
  r_temporal1 = r_temporal2 = r_spacial = 0.0;
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

     if(nu==DIR_T) {
       r_temporal1+=tmp_flt;
     } else {
       if (mu==DIR_T){
	 r_temporal2+=tmp_flt;
       } else {
	 r_spacial+=tmp_flt;
       }
     }
  }

  glb_sum(&r_sum) ;
  glb_sum(&r_sq_sum) ;
  glb_min(&r_min) ;
  glb_max(&r_max) ;
  glb_sum(&r_temporal1) ;
  glb_sum(&r_temporal2) ;
  glb_sum(&r_spacial) ;

  Float one_third = 1.0 / 3.0 ;

  r_sum    *= one_third ;
  r_sq_sum *= one_third * one_third ;
  r_min    *= one_third ;
  r_max    *= one_third ;

  Float r_var = norm_fac * (1.0 + norm_fac) \
                * (r_sq_sum - norm_fac*r_sum*r_sum) ;

  r_sum *= norm_fac ;
  r_temporal1 *= 4.0*norm_fac*one_third;
  r_temporal2 *= 4.0*norm_fac*one_third;
  r_spacial   *= 2.0*norm_fac*one_third;

  // Print out results
  //----------------------------------------------------------------
  if(common_arg->results != 0){
    FILE *fp;
    if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    Fprintf(fp, "%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e \n", IFloat(r_sum), IFloat(r_var),
      IFloat(r_max), IFloat(r_min),
      IFloat(r_temporal1), IFloat(r_temporal2), IFloat(r_spacial) ) ;
    Fclose(fp);
  }

}

CPS_END_NAMESPACE
