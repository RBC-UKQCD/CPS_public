#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:45 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_plaq/alg_plaq.C,v 1.1.1.1 2003-06-22 13:34:45 mcneile Exp $
//  $Id: alg_plaq.C,v 1.1.1.1 2003-06-22 13:34:45 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.7  2002/03/11 22:25:40  anj
//  This should now be the correct, fully merged code from our two versions. Anj
//
//  Revision 1.4.2.1  2002/03/08 16:34:57  anj
//  Checking in the Columbia code branch on tag Columbia4_1_1_test-branch, to be
//  merged with the UKQCD head branch shortly.  Anj
//
//  Revision 1.4  2001/08/16 10:49:40  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:11:28  anj
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
//  Revision 1.2  2001/05/25 06:15:59  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: alg_plaq.C,v $
//  $Revision: 1.1.1.1 $
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
#include <stdlib.h>	// exit()
#include <stdio.h>
#include<alg/alg_plaq.h>
#include<alg/common_arg.h>
#include<alg/no_arg.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<comms/glb.h>
CPS_START_NAMESPACE

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))

//------------------------------------------------------------------
// Constructor 
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
//
//------------------------------------------------------------------
void AlgPlaq::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);

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
      if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
        ERR.FileA(cname,fname, (char *)common_arg->results);
      }
      fprintf(fp, "%f \t %f \t %f \n", 
              (IFloat)ave_space, (IFloat)ave_time, (IFloat)ave_both);  
      fclose(fp);  
    }
    return;
  }
  

  Float p_sum    =  0.0 ;
  Float p_sq_sum =  0.0 ;
  Float p_min    =  3.0 ;
  Float p_max    = -3.0 ;
  int x[4];

  int mu; 
  for (x[0] = 0; x[0] < GJP.XnodeSites(); ++x[0]) {
    for (x[1] = 0; x[1] < GJP.YnodeSites(); ++x[1]) {
      for (x[2] = 0; x[2] < GJP.ZnodeSites(); ++x[2]) {
  	for (x[3] = 0; x[3] < GJP.TnodeSites(); ++x[3]) {
  	  for (mu = 0; mu < 3; ++mu) {
    	    for (int nu = mu+1; nu < 4; ++nu) {
      		Float tmp_flt = lat.ReTrPlaq(x, mu, nu) ;
      		p_sum    += tmp_flt ;
      		p_sq_sum += tmp_flt * tmp_flt ;
      		p_min    =  min(p_min, tmp_flt) ;
      		p_max    =  max(p_max, tmp_flt) ;
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

  Float one_third = 1.0 / 3.0 ;

  p_sum    *= one_third ;
  p_sq_sum *= one_third * one_third ;
  p_min    *= one_third ;
  p_max    *= one_third ;

  Float p_var = norm_fac * (1.0 + norm_fac) \
                * (p_sq_sum - norm_fac*p_sum*p_sum) ;
  p_sum *= norm_fac ;

  x[0]=0 ; x[1]=0 ; x[2]=0 ; x[3]=0 ;

  // Print out results
  //----------------------------------------------------------------
  if(common_arg->results != 0){
    FILE *fp;
    if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    fprintf(fp, "%e %e %e %e %e\n", IFloat(p_sum),
                                 IFloat(p_var),
                                 IFloat(p_max),
                                 IFloat(p_min),
				 IFloat(one_third*lat.ReTrPlaq(x, 0, 1)) ) ;
    fclose(fp);
  }

}





CPS_END_NAMESPACE
