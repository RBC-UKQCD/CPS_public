#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-08-17 03:33:11 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_wline/alg_wline.C,v 1.6 2004-08-17 03:33:11 chulwoo Exp $
//  $Id: alg_wline.C,v 1.6 2004-08-17 03:33:11 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_wline.C,v $
//  $Revision: 1.6 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_wline/alg_wline.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// alg_wline.C
//
// AlgWline is derived from Alg and it measures the average
// value of the Wilson line for each direction.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <alg/alg_wline.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
#include <comms/scu.h>
#include <comms/glb.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
AlgWline::AlgWline(Lattice& latt, 
	     CommonArg *c_arg,
	     NoArg *arg) : 
	     Alg(latt, c_arg) 
{
  cname = "AlgWline";
  char *fname = "AlgWline(L&,CommonArg*,NoArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_wline_arg = arg;

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgWline::~AlgWline() {
  char *fname = "~AlgWline()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
//
//------------------------------------------------------------------
void AlgWline::run()
{
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif

  char *fname = "run()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();


  // Compute the Wilson Line for each direction
  //----------------------------------------------------------------

  Matrix accum_link, in_link, out_link ;
  Matrix *gauge=lat.GaugeField() ;

  int num_nodes[4]
    = { GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes() } ;

  int node_sites[4]
    = { GJP.XnodeSites(), GJP.YnodeSites(), GJP.ZnodeSites(), GJP.TnodeSites() } ;

  Float norm[4] ;

  norm[0] = 1.0 / (Float)(   node_sites[1] * node_sites[2] * node_sites[3]
                           *  num_nodes[1] *  num_nodes[2] *  num_nodes[3] ) ;
  norm[1] = 1.0 / (Float)(   node_sites[0] * node_sites[2] * node_sites[3]
                           *  num_nodes[0] *  num_nodes[2] *  num_nodes[3] ) ;
  norm[2] = 1.0 / (Float)(   node_sites[0] * node_sites[1] * node_sites[3]
                           *  num_nodes[0] *  num_nodes[1] *  num_nodes[3] ) ;
  norm[3] = 1.0 / (Float)(   node_sites[0] * node_sites[1] * node_sites[2]
                           *  num_nodes[0] *  num_nodes[1] *  num_nodes[2] ) ;

  // Added for anisotropic lattices
  Float norm_factor = 1/GJP.XiBare();
  // End modification

  for (int mu=0; mu<4; mu++) {
    VRB.Debug(cname, fname, "Begin Direction = %i\n", mu) ;

    Complex wline[4] = { 0.0, 0.0, 0.0, 0.0 } ;
    int x[4] ;

    for (x[(mu+1)%4]=0; x[(mu+1)%4] < node_sites[(mu+1)%4]; x[(mu+1)%4]++) 
    for (x[(mu+2)%4]=0; x[(mu+2)%4] < node_sites[(mu+2)%4]; x[(mu+2)%4]++) 
    for (x[(mu+3)%4]=0; x[(mu+3)%4] < node_sites[(mu+3)%4]; x[(mu+3)%4]++) {

      // for these three coords fixed, calc local wline

      accum_link.UnitMatrix() ;

      for ( x[mu]=0; x[mu] < node_sites[mu]; x[mu]++ ) {
        out_link.DotMEqual(accum_link, gauge[mu+lat.GsiteOffset(x)]) ;

	// Added for anisotropic lattices
	if (GJP.XiBare() != 1.0 && mu == GJP.XiDir())
	  out_link *= norm_factor;
	// End modification 

        accum_link = out_link ;
      }

      // now pass local segments between nodes and accum
      // TRICK: out_link already contains local segment

      for (int node_cntr=1; node_cntr < num_nodes[mu]; node_cntr++) {
        getMinusData((IFloat *)&in_link, (IFloat *)&out_link,
          sizeof(Matrix)/sizeof(IFloat), mu) ;
        out_link.DotMEqual(in_link, accum_link) ;
        accum_link = out_link ;
        out_link = in_link ;
      }

      // compute the characters for the 3, 6, 8, and 10 reps.

      wline[0] += accum_link.Char3() ;
      wline[1] += accum_link.Char6() ;
      wline[2] += accum_link.Char8() ;
      wline[3] += accum_link.Char10() ;

    } // end for loop over coords != mu

    slice_sum((Float *)wline, 4*sizeof(Complex)/sizeof(IFloat), mu) ;

    wline[0] *= norm[mu] ;
    wline[1] *= norm[mu] ;
    wline[2] *= norm[mu] ;
    wline[3] *= norm[mu] ;

    // Print out results
    //----------------------------------------------------------------

    VRB.Debug(cname, fname, "%e %e %e %e %e %e %e %e\n",
        wline[0].real(), wline[0].imag(),
        wline[1].real(), wline[1].imag(),
        wline[2].real(), wline[2].imag(),
        wline[3].real(), wline[3].imag() );
 
    if(common_arg->results != 0){
      FILE *fp;
      if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
        ERR.FileA(cname,fname, (char *)common_arg->results);
      }
      Fprintf(fp, "%e %e %e %e %e %e %e %e\n",
        wline[0].real(), wline[0].imag(),
        wline[1].real(), wline[1].imag(),
        wline[2].real(), wline[2].imag(),
        wline[3].real(), wline[3].imag() );
      Fclose(fp);
    }

    VRB.Debug(cname, fname, "End Direction = %i\n", mu) ;
  } // end for mu
}

CPS_END_NAMESPACE
