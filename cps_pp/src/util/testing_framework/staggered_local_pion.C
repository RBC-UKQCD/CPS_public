/*!\file
  \brief Pion correlation function for the testing framework.

  $Id: staggered_local_pion.C,v 1.6 2004-12-11 20:58:06 chulwoo Exp $
*/
//--------------------------------------------------------------------
// Compute the local pion correlator for staggered fermions.
// This is a light weight code for the testing framework.
//
// To write the code to compute the pion correlator, I used
// the example in void AlgThreePt::spectrum.
//
//
//  Starting point
//  $Date: 2004-12-11 20:58:06 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/testing_framework/staggered_local_pion.C,v 1.6 2004-12-11 20:58:06 chulwoo Exp $
//
//-------------------------------------------------------------------

#include <config.h>
#include <util/qcdio.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <comms/glb.h>
#include <alg/w_ginfo.h>


#include<assert.h>

CPS_START_NAMESPACE

/*!
  The correlation function is computed using a point source at the origin.
  \param lat The lattice object containing the gauge configuration.
  It is assumed that there are no spin degrees of freedom, so this should
  be a class derived from FstagTypes.
  \param mass Valence quark mass parameter
  \param pion_corr The timesliced pseudoscalar correlation funtion
  \param The number of timeslices over which to compute the correlator.
  \post The correlator is written to the array and also printed to \c stdout.
*/
void staggered_local_pion(Lattice &lat, Float mass, 
			  IFloat* pion_corr, int time_size)
{
  printf("Start of propagator inversion\n") ; 
  printf("Physical mass = %f\n",mass) ; 

  int cg_iter ;
  Float true_res;
  IFloat *src ;
  IFloat *sln ; 

  // the cg data structure is the only one needed

  WspectArg w_spect_arg;
  CgArg cg;
  cg.mass = mass ;
  cg.stop_rsd = 1.0E-12;
  cg.max_num_iter = 5000;

  WspectGinfo *con_ ; 
  // constructor is protected, so need hack

  assert(con_->COLORs == 3 ); 
  const int COLORs = con_->COLORs ;

  //   alg/alg_w_spect/w_quark.C
  int v_size = GJP.VolNodeSites() * (con_->COLORs *con_->COMPLEXs) ; 
  int local_size = GJP.VolNodeSites() ; 

  // 
  //  memory for colour vector over the lattice
  //
  src = (IFloat*)smalloc(sizeof(IFloat) * v_size  )  ;
  sln = (IFloat*)smalloc(sizeof(IFloat) * v_size  )  ;

  //  int time_size=GJP.Tnodes()*GJP.TnodeSites();
  //IFloat* pion_corr=(IFloat*)smalloc(time_size*sizeof(IFloat));

  int t ;
  for(t=0; t<time_size; t++)
    pion_corr[t] = 0.0 ; 

  int shift_t=GJP.TnodeCoor()*GJP.TnodeSites();
  int vol= (GJP.VolNodeSites()/GJP.TnodeSites());



  // need spin as well
  for(int color=0 ; color < 3 ;++color)
    {
      printf("Starting inversion for color  = %d\n",color) ; 


      // set initial guess for the inverter
      for(int is = 0 ; is < v_size ; ++ is)
	sln[is] = 0.0 ; 

      // create the local source (ignore || stuff first)
      for(int is = 0 ; is < v_size ; ++ is)
	src[is] = 0.0 ; 

      // set (0,0,0,0) = 1
      //
      // fermion_field[t][z][y][x][colour][comp]
      //

      // more work not parallel
      *(src + 0 + 2*(color)) = 1.0 ; 

      cg_iter = lat.FmatInv((Vector *)sln, (Vector *)src, 
			    &cg, &true_res, CNV_FRM_YES);

      printf("CG iters %d\n",cg_iter); 
      printf("True residual = %e\n",true_res); 

      //
      // compute the pion
      //

      for(int sink_color=0 ; sink_color < 3 ;++sink_color)
	  {
	      for(int is=0; is< GJP.VolNodeSites(); is++)
	      {
		t=is/vol;
		t+=shift_t;

		//
		// fermion_field[t][z][y][x][colour][comp]
		//

		int pt = 2*(sink_color + COLORs*is );
		IFloat 
		  tmp_re = *(sln + pt) ; 

		IFloat 
		  tmp_im = *(sln + 1 + pt ) ; 

		pion_corr[t] += 
		  tmp_re * tmp_re  
		  + tmp_im * tmp_im ; 

	      } // end loop over lattice volume on this node


	      } // end loop over sink color


    } // end loop over source color

  // sum the correlators over the nodes
  for(t=0; t<time_size; t++)
    glb_sum (&pion_corr[t]) ; 

  // write the pion correlator to the screem
  for(t=0; t<time_size; t++)
    printf("pion[%d] = %f\n",t,pion_corr[t] ); 

  // free memory 
  sfree(src) ;
  sfree(sln) ;


  printf("END of propagator inversion \n") ; 

}


CPS_END_NAMESPACE
