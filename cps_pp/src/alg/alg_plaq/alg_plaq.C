#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief Definitions of the AlgPlaq class methods.
  
  $Id: alg_plaq.C,v 1.10 2005-05-12 20:14:10 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2005-05-12 20:14:10 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_plaq/alg_plaq.C,v 1.10 2005-05-12 20:14:10 chulwoo Exp $
//  $Id: alg_plaq.C,v 1.10 2005-05-12 20:14:10 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_plaq.C,v $
//  $Revision: 1.10 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_plaq/alg_plaq.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// alg_plaq.C
//
// AlgPlaq is derived from Alg and it measures the average
// value of the plaquette. 
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <alg/alg_plaq.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <comms/glb.h>
#include <util/qcdio.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
/*!
  \param latt The lattice on which to compute the plaquette.
  \param c_arg The common argument structure for all algorithms.
  \param arg A dummy argument structure.
 */
//------------------------------------------------------------------
AlgPlaq::AlgPlaq(Lattice& latt, 
	     CommonArg *c_arg,
	     NoArg *arg) : 
	     Alg(latt, c_arg) 
{

    cname = "AlgPlaq";
  char *fname = "AlgPlaq(L&,CommonArg*,NoArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_plaq_arg = arg;


  // Calculate normalization factor
  //----------------------------------------------------------------
  int total_sites = GJP.VolNodeSites() * GJP.Xnodes() *
                    GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
  norm_fac = 1.0 / ( 6.0 * Float(total_sites) );

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgPlaq::~AlgPlaq() {
  char *fname = "~AlgPlaq()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
//! Performs the computation.
/*!
  The real trace of the plaquette is averaged over the lattice and normalised
  to unity.
  \post The following data are written to the file specified in the
  common_arg structure:
  -# mean plaquette
  -# variance of mean plaquette
  -# maximum plaquette
  -# minimum plaquette
  -# mean temporal plaquette
  -# mean spatial plaquette
*/
//------------------------------------------------------------------
void AlgPlaq::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);

#if TARGET==cpsMPI
  using MPISCU::fprintf;
#endif
  
  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  // Modified by Ping Chen for anisotropic lattices
  //----------------------------------------------------------------
  if (GJP.XiBare() != 1.0 || GJP.XiV() != 1 || GJP.XiVXi() != 1   ) {
    Float ave_time  = lat.AveReTrPlaqXi();
    Float ave_space = lat.AveReTrPlaqNoXi();
    Float ave_both  = (ave_time + ave_space) / 2;  
 
    if (common_arg->results != 0){
      FILE *fp;
      if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
        ERR.FileA(cname,fname, (char *)common_arg->results);
      }
      Fprintf(fp, "%f \t %f \t %f \n", 
              (IFloat)ave_space, (IFloat)ave_time, (IFloat)ave_both);  
      Fclose(fp);  
    }
    return;
  }

  Float p_temporal, p_spatial;
  p_temporal = p_spatial = 0.0;

  Float p_sum    =  0.0 ;
  Float p_sq_sum =  0.0 ;
  Float p_min    =  3.0 ;
  Float p_max    = -3.0 ;


  int x[4];
  for (x[0] = 0; x[0] < GJP.XnodeSites(); ++x[0]) {
    for (x[1] = 0; x[1] < GJP.YnodeSites(); ++x[1]) {
      for (x[2] = 0; x[2] < GJP.ZnodeSites(); ++x[2]) {
  	for (x[3] = 0; x[3] < GJP.TnodeSites(); ++x[3]) {
  	  for (int mu = 0; mu < 3; ++mu) {
    	    for (int nu = mu+1; nu < 4; ++nu) {
      		Float tmp_flt = lat.ReTrPlaq(x, mu, nu) ;
      		p_sum    += tmp_flt ;
      		p_sq_sum += tmp_flt * tmp_flt ;
      		p_min    =  p_min<tmp_flt ? p_min : tmp_flt;
      		p_max    =  p_max>tmp_flt ? p_max : tmp_flt;
		if(nu==DIR_T||mu==DIR_T) p_temporal += tmp_flt;
		else p_spatial += tmp_flt;
	    }
	  }
	}
      }
    }
  }

  glb_sum(&p_sum) ;
  glb_sum(&p_sq_sum) ;
  glb_min(&p_min) ;
  glb_max(&p_max) ;
  glb_sum(&p_spatial);
  glb_sum(&p_temporal);
  

  Float one_third = 1.0 / 3.0 ;

  p_sum    *= one_third ;
  p_sq_sum *= one_third * one_third ;
  p_min    *= one_third ;
  p_max    *= one_third ;

  Float p_var = norm_fac * (1.0 + norm_fac) \
                * (p_sq_sum - norm_fac*p_sum*p_sum) ;
  p_sum *= norm_fac ;
  p_spatial *= 2.0*norm_fac*one_third;
  p_temporal *= 2.0*norm_fac*one_third;


  // Print out results
  //----------------------------------------------------------------
  if(common_arg->results != 0){
    FILE *fp;

    if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) 
	ERR.FileA(cname,fname, (char *)common_arg->results);
    
    Fprintf(fp, "%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e\n",
	     p_sum, p_var,
	     p_max, p_min,
	     p_temporal, p_spatial);
    Fclose(fp);
  }

}






CPS_END_NAMESPACE
