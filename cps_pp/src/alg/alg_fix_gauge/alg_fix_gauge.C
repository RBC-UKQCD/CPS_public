#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-08-17 03:33:08 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_fix_gauge/alg_fix_gauge.C,v 1.7 2004-08-17 03:33:08 chulwoo Exp $
//  $Id: alg_fix_gauge.C,v 1.7 2004-08-17 03:33:08 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_fix_gauge.C,v $
//  $Revision: 1.7 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_fix_gauge/alg_fix_gauge.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// AlgFixGauge is derived from Alg and is relevant to the 
// gauge fixing algorithms. The type of glue and fermion is
// determined by the argument to the constructor.
//
//------------------------------------------------------------------


CPS_END_NAMESPACE
#include <stdlib.h>	// exit()
#include <util/qcdio.h>
#include <alg/alg_fix_gauge.h>
#include <alg/common_arg.h>
#include <alg/fix_gauge_arg.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
CPS_START_NAMESPACE



//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
AlgFixGauge::AlgFixGauge(Lattice& latt, 
			 CommonArg *c_arg,
			 FixGaugeArg *arg) :
			 Alg(latt, c_arg) 
{
  cname = "AlgFixGauge";
  char *fname = "AlgFixGauge(L&,CommonArg*,FixGaugeArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_fix_gauge_arg = arg;


  //???
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgFixGauge::~AlgFixGauge() {
  char *fname = "~AlgFixGauge()";
  VRB.Func(cname,fname);


  //???
}


//------------------------------------------------------------------
// Allocates memory and constructs the gauge fixing matrices
//------------------------------------------------------------------
void AlgFixGauge::run()
{
  char *fname = "run()";
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
  int lattice_dir_size;
  int *h_planes = 0;

  // If coulomb gauge
  //----------------------------------------------------------------
  if( (fix == FIX_GAUGE_COULOMB_X) ||
      (fix == FIX_GAUGE_COULOMB_Y) ||
      (fix == FIX_GAUGE_COULOMB_Z) ||
      (fix == FIX_GAUGE_COULOMB_T)  ){ 
    if(fix == FIX_GAUGE_COULOMB_X) 
      lattice_dir_size = GJP.XnodeSites() * GJP.Xnodes();
    else if(fix == FIX_GAUGE_COULOMB_Y) 
      lattice_dir_size = GJP.YnodeSites() * GJP.Ynodes();
    else if(fix == FIX_GAUGE_COULOMB_Z) 
      lattice_dir_size = GJP.ZnodeSites() * GJP.Znodes();
    else 
      lattice_dir_size = GJP.TnodeSites() * GJP.Tnodes();

    if(start+step*(num-1) >= lattice_dir_size){
      ERR.General(cname, fname, 
		  "Wrong fix_gauge_arg, start = %d, step = %d, num = %d\n",
		  start, step, num);
    }

    h_planes = (int *) smalloc(num * sizeof(int));
    if(h_planes == 0)
      ERR.Pointer(cname,fname, "h_planes");
    VRB.Smalloc(cname,fname, "h_planes", h_planes, num * sizeof(int));
    
    for(int i=0; i<num; i++){
      h_planes[i] = start + step * i;
    }

  }

  // Allocate gauge fixing matrices
  //----------------------------------------------------------------
  lat.FixGaugeAllocate(fix, num, h_planes);


  // Calculate the gauge fixing matrices
  //----------------------------------------------------------------
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
