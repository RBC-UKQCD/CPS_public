/*!\file
  \brief Implementation of AlgEqState class methods.

  $Id: alg_eq_state.C,v 1.8 2004/09/02 17:00:12 zs Exp $
*/

#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004/09/02 17:00:12 $
//  $Header: /space/cvs/cps/cps++/src/alg/alg_eq_state/alg_eq_state.C,v 1.8 2004/09/02 17:00:12 zs Exp $
//  $Id: alg_eq_state.C,v 1.8 2004/09/02 17:00:12 zs Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: alg_eq_state.C,v $
//  $Revision: 1.8 $
//  $Source: /space/cvs/cps/cps++/src/alg/alg_eq_state/alg_eq_state.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

//------------------------------------------------------------------
//
// alg_eq_state.C
//
// AlgEqState is derived from Alg and it measures the sum, normalized
// by volume, of the plaquettes on particular hyperplane(s).
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <alg/alg_eq_state.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
#include <comms/glb.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
/*!
  \param latt The lattice object containing the gauge filed on which to
  compute the plaquette.
  \param arg   Container for parameters specific to this algorithm.
  \param c_arg Container for generic algorithm parameters.
 */
AlgEqState::AlgEqState(Lattice& latt, 
	     CommonArg *c_arg,
	     EqStateArg *arg) : 
	     Alg(latt, c_arg) 
{
  cname = "AlgEqState";
  char *fname = "AlgEqState(L&,CommonArg*,EqStateArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_eq_state_arg = arg;


  // Calculate normalization factor
  //----------------------------------------------------------------
  int total_sites = GJP.VolNodeSites() * GJP.Xnodes() *
                    GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
  norm_fac = 1.0 / (18.0 * Float(total_sites));
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgEqState::~AlgEqState() {
  char *fname = "~AlgEqState()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
/*!
  \post If a file is specified in CommonArg, the plaquette values in the
  perpendicular and then parallel to the specified direction are appended
  to this file.
 */
//------------------------------------------------------------------
void AlgEqState::run()
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
  if (GJP.XiBare() != 1.0 || GJP.XiV() != 1 || GJP.XiVXi() != 1 ) {

    if (GJP.XiDir() != alg_eq_state_arg->dir)
      ERR.General(cname,fname,
	 "EOS direction is different with anisotropic direction\n");

    Float ave_time  = lat.AveReTrPlaqXi();
    Float ave_space = lat.AveReTrPlaqNoXi();

    if (common_arg->results != 0){
      FILE *fp;
      if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
        ERR.FileA(cname,fname, (char *)common_arg->results);
      }
      Fprintf(fp, "%e\t%e\n", (IFloat)ave_space, (IFloat)ave_time);  
      Fclose(fp);  
    }
    return;
  }

  // Calculate the sum of plaquettes
  //----------------------------------------------------------------
  Float plaq_perpen = 0.0;
  Float plaq_parall = 0.0;
  int x[4];
    
  for(x[0] = 0; x[0] < GJP.XnodeSites(); ++x[0]) 
      for(x[1] = 0; x[1] < GJP.YnodeSites(); ++x[1]) 
	  for(x[2] = 0; x[2] < GJP.ZnodeSites(); ++x[2]) 
	      for(x[3] = 0; x[3] < GJP.TnodeSites(); ++x[3]) {
  
		  for (int mu = 0; mu<3; ++mu) for(int nu = mu+1; nu<4; ++nu){
		      Float re_tr_plaq = lat.ReTrPlaq(x,mu,nu);
		      if (mu == alg_eq_state_arg->dir || nu == alg_eq_state_arg->dir) 
			  plaq_parall += re_tr_plaq;
		      else 
			  plaq_perpen += re_tr_plaq;
			  
		  }
	  
	      }
	  
      
  

  glb_sum(&plaq_perpen);
  glb_sum(&plaq_parall);

  plaq_perpen *= norm_fac;
  plaq_parall *= norm_fac;


  // Print out results
  //----------------------------------------------------------------
  if(common_arg->results != 0){
    FILE *fp;
    if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    Fprintf(fp, "%e\t%e\n", IFloat(plaq_perpen), IFloat(plaq_parall));
    Fclose(fp);
  }
}






CPS_END_NAMESPACE
