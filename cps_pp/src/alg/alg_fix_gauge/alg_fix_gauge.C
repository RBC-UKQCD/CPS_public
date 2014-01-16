#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief Implementation of AlgFixGauge class methods.

  $Id: alg_fix_gauge.C,v 1.11 2007-06-25 15:49:20 chulwoo Exp $
*/
//------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2007-06-25 15:49:20 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_fix_gauge/alg_fix_gauge.C,v 1.11 2007-06-25 15:49:20 chulwoo Exp $
//  $Id: alg_fix_gauge.C,v 1.11 2007-06-25 15:49:20 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_fix_gauge.C,v $
//  $Revision: 1.11 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_fix_gauge/alg_fix_gauge.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


CPS_END_NAMESPACE
//#include <stdlib.h>	// exit()
//#include <util/qcdio.h>
#include <alg/alg_fix_gauge.h>
// #include <alg/common_arg.h>
// #include <alg/fix_gauge_arg.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
CPS_START_NAMESPACE



//------------------------------------------------------------------
/*!
  \param latt The lattice object containing the gauge field to be fixed.
  \param c_arg Generic algorithm parameters
  \param arg Gauge fixing parameters.
*/
//------------------------------------------------------------------
AlgFixGauge::AlgFixGauge(Lattice& latt, 
			 CommonArg *c_arg,
			 FixGaugeArg *arg) :
			 Alg(latt, c_arg) 
{
  cname = "AlgFixGauge";
  const char *fname = "AlgFixGauge";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_fix_gauge_arg = arg;

}


//------------------------------------------------------------------
/*!
  \note The destructor does not free the memory allocated.
  Use AlgFixGauge::free for this or Lattice::FixGaugeFree
*/
//------------------------------------------------------------------
AlgFixGauge::~AlgFixGauge() {
  const char *fname = "~AlgFixGauge";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// Allocates memory and constructs the gauge fixing matrices
//------------------------------------------------------------------
void AlgFixGauge::run()
{
  const char *fname = "run";
  VRB.Func(cname,fname);

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  // Set up arguments
  //----------------------------------------------------------------
  FixGaugeType fix = alg_fix_gauge_arg->fix_gauge_kind;
  int start = alg_fix_gauge_arg->hyperplane_start;
  int step = alg_fix_gauge_arg->hyperplane_step;
  int num = alg_fix_gauge_arg->hyperplane_num;
  int lattice_dir_size=0;
  int *h_planes = 0;

  // If coulomb gauge
  //----------------------------------------------------------------
  if( (fix == FIX_GAUGE_COULOMB_X) ||
      (fix == FIX_GAUGE_COULOMB_Y) ||
      (fix == FIX_GAUGE_COULOMB_Z) ||
      (fix == FIX_GAUGE_COULOMB_T)  ){ 

      switch(fix){
      case FIX_GAUGE_COULOMB_X:
	  lattice_dir_size = GJP.XnodeSites() * GJP.Xnodes();
	  break;
      case FIX_GAUGE_COULOMB_Y:
	  lattice_dir_size = GJP.YnodeSites() * GJP.Ynodes();
	  break;
      case FIX_GAUGE_COULOMB_Z:
	  lattice_dir_size = GJP.ZnodeSites() * GJP.Znodes();
	  break;
      case FIX_GAUGE_COULOMB_T:
	  lattice_dir_size = GJP.TnodeSites() * GJP.Tnodes();
	  break;
      case FIX_GAUGE_NONE:
      case FIX_GAUGE_LANDAU:
	  break;
      }

    if(start+step*(num-1) >= lattice_dir_size)
	ERR.General(cname, fname,
		    "The coordinate of the last hyperplane (%d+%d*%d) is greater than the global lattice size.", 
		    start, step, num-1, lattice_dir_size);

    h_planes = (int *) smalloc(num * sizeof(int));
    if(h_planes == 0) ERR.Pointer(cname,fname, "h_planes");
    VRB.Smalloc(cname,fname, "h_planes", h_planes, num * sizeof(int));
    
    for(int i=0; i<num; i++) h_planes[i] = start + step * i;

  }

  // Allocate gauge fixing matrices and set them to 1
  //----------------------------------------------------------------
  lat.FixGaugeAllocate(fix, num, h_planes);


  // Calculate the gauge fixing matrices
  //----------------------------------------------------------------
// added to make it possible to explicitly bypass gauge fixin:w
  if ( alg_fix_gauge_arg->max_iter_num > 0 ) 
  lat.FixGauge(alg_fix_gauge_arg->stop_cond, 
	       alg_fix_gauge_arg->max_iter_num);


  VRB.Sfree(cname,fname, "h_planes",h_planes);
  sfree(h_planes);

}

//------------------------------------------------------------------
// Free the memory of the gauge fixing matrices.
//------------------------------------------------------------------
void AlgFixGauge::free()
{
  char *fname = "free()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  // Free gauge fixing matrices
  //----------------------------------------------------------------
  lat.FixGaugeFree();
}


CPS_END_NAMESPACE
