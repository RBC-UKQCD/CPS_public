#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_eq_state/alg_eq_state.C,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Id: alg_eq_state.C,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:49:38  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:11:21  anj
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
//  $RCSfile: alg_eq_state.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_eq_state/alg_eq_state.C,v $
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
// Reduntant #include <stdlib.h>	// exit()
#include <stdio.h>
#include <alg/alg_eq_state.h>
#include <alg/common_arg.h>
#include <alg/eq_state_arg.h>
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
//
//------------------------------------------------------------------
void AlgEqState::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);

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
      if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
        ERR.FileA(cname,fname, (char *)common_arg->results);
      }
      fprintf(fp, "%e\t%e\n", (IFloat)ave_space, (IFloat)ave_time);  
      fclose(fp);  
    }
    return;
  }

  // Calculate the sum of plaquettes
  //----------------------------------------------------------------
  Float plaq_perpen = 0.0;
  Float plaq_parall = 0.0;
  int x[4];
    
  for(x[0] = 0; x[0] < GJP.XnodeSites(); ++x[0]) {
    for(x[1] = 0; x[1] < GJP.YnodeSites(); ++x[1]) {
      for(x[2] = 0; x[2] < GJP.ZnodeSites(); ++x[2]) {
        for(x[3] = 0; x[3] < GJP.TnodeSites(); ++x[3]) {
  
          for (int mu = 0; mu < 3; ++mu) {
            for(int nu = mu+1; nu < 4; ++nu) {
	      Float re_tr_plaq = lat.ReTrPlaq(x,mu,nu);
	      if (mu == alg_eq_state_arg->dir || nu == alg_eq_state_arg->dir) {
		plaq_parall += re_tr_plaq;
	      } else {
		plaq_perpen += re_tr_plaq;
	      }
            }
          }
        }
      }
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
    if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    fprintf(fp, "%e\t%e\n", IFloat(plaq_perpen), IFloat(plaq_parall));
    fclose(fp);
  }
}






CPS_END_NAMESPACE
