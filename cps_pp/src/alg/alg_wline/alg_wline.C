#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief Implementation of AlgWline class methods.

  $Id: alg_wline.C,v 1.12 2007-10-30 20:40:34 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2007-10-30 20:40:34 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_wline/alg_wline.C,v 1.12 2007-10-30 20:40:34 chulwoo Exp $
//  $Id: alg_wline.C,v 1.12 2007-10-30 20:40:34 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_wline.C,v $
//  $Revision: 1.12 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_wline/alg_wline.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

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
#if TARGET == BGL
//#include <sys/bgl/bgl_sys_all.h>
#endif
CPS_START_NAMESPACE

//------------------------------------------------------------------
/*!
  \param latt The Lattice object containg the gauge field on which to compute the %Wilson lines.
  \param c_arg Container for generic parameters. .
  \param arg Empty parameter container.
*/
//------------------------------------------------------------------
AlgWline::AlgWline(Lattice& latt, 
	     CommonArg *c_arg,
	     NoArg *arg) : 
	     Alg(latt, c_arg) 
{
  cname = "AlgWline";
  const char *fname = "AlgWline";
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
  const char *fname = "~AlgWline";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
/*!
  The characters of the result for the 3, 6, 8, and 10 representations
  are calculated.
  If an output file is specified in the CommonArg argument, then the real
  and imaginary parts of these are written to it all on one line per direction.
*/  
//------------------------------------------------------------------
void AlgWline::run()
{
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif

  const char *fname = "run";
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


  for (int mu=0; mu<4; mu++) {
    VRB.Debug(cname, fname, "Begin Direction = %i\n", mu) ;

    Complex wline[4] CPS_FLOAT_ALIGN;
    wline[0] = 0.0;
    wline[1] = 0.0;
    wline[2] = 0.0;
    wline[3] = 0.0;
//    for(int i =0;i<4;i++)
//    printf("Node %d: wline[%d]= %e %e\n",UniqueID(),i,wline[i].real(),wline[i].imag());
     

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

//    for(int i =0;i<4;i++)
//    printf("Node %d: wline[%d]= %e %e\n",UniqueID(),i,wline[i].real(),wline[i].imag());
//    slice_sum((Float *)wline, 4*sizeof(Complex)/sizeof(IFloat), mu) ;
    slice_sum((Float *)wline, 8, mu) ;
//    for(int i =0;i<4;i++)
//    printf("Node %d: after slice_sum : wline[%d]= %e %e\n",UniqueID(),i,wline[i].real(),wline[i].imag());

    wline[0] *= norm[mu] ;
    wline[1] *= norm[mu] ;
    wline[2] *= norm[mu] ;
    wline[3] *= norm[mu] ;

    // Print out results
    //----------------------------------------------------------------

    VRB.Debug(cname, fname, "%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e\n",
        wline[0].real(), wline[0].imag(),
        wline[1].real(), wline[1].imag(),
        wline[2].real(), wline[2].imag(),
        wline[3].real(), wline[3].imag() );
 
    if(common_arg->filename != 0){
      FILE *fp;
       if( (fp = Fopen(common_arg->filename, "a")) == NULL ) {
         ERR.FileA(cname,fname, (char *)common_arg->filename);
       }
       Fprintf(fp, "%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e\n",
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
