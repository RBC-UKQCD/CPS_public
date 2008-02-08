#include <config.h>

#include <stdio.h>
#include <stdlib.h>	// exit()

#if TARGET==QCDOC
#include <sysfunc_cps.h>
#endif

#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <alg/alg_ghb.h>
#include <alg/alg_plaq.h>
#include <alg/do_arg.h>
#include <alg/common_arg.h>
#include <alg/ghb_arg.h>
#include <alg/no_arg.h>
#include <alg/fix_gauge_arg.h>
#include <alg/alg_fix_gauge.h>
#include <alg/fourierprop_arg.h>
#include <alg/qpropw.h>
#include <alg/alg_fourier_prop.h>
#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>

USING_NAMESPACE_CPS

int main(int argc,char *argv[])
{

  char *cname = argv[0] ;
  char *fname = "main()" ;

  if(argc<2) {
    ERR.General(cname, fname, "Please provide the lattice name to load\n");
  }

  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------

  DoArg do_arg;

  //=========================
  // number of masses to run
  //=========================

  int num_mass;

  //=========================
  // array of the masses
  //=========================

  Float *mass;

  //do_arg.start_seed_kind = START_SEED_INPUT;

  int i = 0;

  //-----------------------------------------------------------------
  // Command Line Input
  //-----------------------------------------------------------------
  
  //======================
  // Machine orientation
  //======================

  do_arg.x_nodes = SizeX();
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();
  do_arg.s_nodes = SizeS();

  do_arg.x_node_sites = 16/SizeX();
  do_arg.y_node_sites = 16/SizeY();
  do_arg.z_node_sites = 16/SizeZ();
  do_arg.t_node_sites = 32/SizeT();
  do_arg.s_node_sites = 16/SizeS();

  do_arg.x_sites = 16;
  do_arg.y_sites = 16;
  do_arg.z_sites = 16;
  do_arg.t_sites = 32;
  do_arg.s_sites = 16;

  //======================
  // Boundary Conditions
  //======================

  do_arg.x_bc = BND_CND_PRD ;
  do_arg.y_bc = BND_CND_PRD ;
  do_arg.z_bc = BND_CND_PRD ;
  do_arg.t_bc = BND_CND_PRD ;

  //==========================
  // type of lattice to load
  //==========================

  do_arg.start_conf_kind = START_CONF_LOAD;

  //==============================
  // Position of Lattice to Load 
  //==============================

  do_arg.start_conf_load_addr = 
    (u_long)smalloc(do_arg.x_node_sites * do_arg.y_node_sites *
		    do_arg.z_node_sites * do_arg.t_node_sites * 4 * sizeof(Matrix));

  do_arg.start_conf_filename = argv[1];

  //======================
  // seed value and type
  //======================

  do_arg.start_seed_kind = START_SEED_FIXED ;
  do_arg.start_seed_value = 3323;

  //======================
  // physics parameters
  //======================

  //  do_arg.colors = 3;  //not needed any more
 
  do_arg.beta = 6.0;

  do_arg.dwf_height = 1.7;
  
  //=========================
  // Gauge fixing arguments
  //=========================

  FixGaugeArg fix_arg;
  fix_arg.fix_gauge_kind   = FIX_GAUGE_LANDAU;
  fix_arg.hyperplane_start = 0;
  fix_arg.hyperplane_step  = 0;
  fix_arg.hyperplane_num   = 0;

  fix_arg.stop_cond        = 1e-8;

  fix_arg.max_iter_num = 0;  // disable GF to compare with QCDSP
  //  fix_arg.max_iter_num     = 20000;  // normally should use this value

  //========================
  // Momenta to calculate
  //========================

  int max_px;
  int max_py;
  int max_pz;
  int max_pt;

  max_px = 2;
  max_py = 2;
  max_pz = 2;
  max_pt = 4;
 
  //==============================
  // fourier transform parameters
  //==============================

  FourierPropArg fourierprop_arg;
  
  fourierprop_arg.cg.max_num_iter = 5001;
  fourierprop_arg.cg.stop_rsd     = 1.0E-8;

  //====================
  // Hard-wired options
  //====================

  fourierprop_arg.results = "gffprop_cd.dat";
  fourierprop_arg.x_src           = 0;
  fourierprop_arg.y_src           = 0;
  fourierprop_arg.z_src           = 0;
  fourierprop_arg.t_src           = 0;

  
  //================
  // masses to run
  //================
  
  num_mass = 5;

  //===================================
  // allocate the array for the masses
  //===================================
  // delete at the end of the program
  //===================================

  mass = new Float[num_mass];

  //=================
  // read in masses 
  //=================

  mass[0] = 0.01;
  mass[1] = 0.02;
  mass[2] = 0.03;
  mass[3] = 0.04;
  mass[4] = 0.05;

  //==========================
  // Initialize the GJP class
  //==========================

  GJP.Initialize(do_arg);

  //===================
  // Set verbose level
  //===================

  VRB.Level(VERBOSE_RESULT_LEVEL);
  
  //================================
  // Initialize argument structures
  //================================

  CommonArg common_arg;
  CommonArg common_arg_plaq;

  NoArg plaq_arg;

  GwilsonFdwf lat;

  {

    // Read lattice
    ReadLatticeParallel rd;
    rd.read(lat, do_arg.start_conf_filename);

    //    printf("check loading\n");

    //    WriteLatticeParallel wt;
    //    wt.write(lat, "check_lat.dat");

    printf("Plaq\n");

    common_arg_plaq.results = (char*)"plaq.dat";
    AlgPlaq plaq(lat,&common_arg_plaq,&plaq_arg);
    plaq.run();
  }
  
  //===============
  // fix the gauge
  //===============
    
  {

    printf("Fix gauge\n");

    AlgFixGauge fix_gauge(lat,&common_arg,&fix_arg);
    fix_gauge.run();
  } 
  
  //===========================
  // Contruct array of momenta
  //===========================
  
  printf("Construct array\n");

  fourierprop_arg.plist.size(  ( 2 * max_px + 1 )
                               * ( 2 * max_py + 1 )
                               * ( 2 * max_pz + 1 )
                               * ( 2 * max_pt + 1 ) );

  int four_sum(0);
  
  int p1,p2,p3,p4;

  for ( p1=-max_px; p1<=max_px; ++p1)
    {
      for ( p2=-max_py; p2 <= max_py ; ++p2 )
	{
	  for ( p3=-max_pz ; p3<=max_pz ; ++p3)
	    {
	      for ( p4=-max_pt; p4<=max_pt ; ++p4)
		{
		  fourierprop_arg.plist[four_sum]= FourMom(p1,
                                                           p2,
                                                           p3,
                                                           p4);
		  four_sum++;
		  if (four_sum > fourierprop_arg.plist.size() ){
		    printf("plist.size() exceeded\n");
		    exit(13);
		  }
		}
	    }
	}
    }
  
  //====================
  // fourier transform
  //====================
  {
    printf("DFT\n");

    common_arg.results = (char*)"info.dat";
    
    for( int i=0; i< num_mass; i++)
      {
        
        fourierprop_arg.cg.mass = mass[i];
        AlgFourierProp fprop(lat,&common_arg,&fourierprop_arg);
        fprop.run();
        
      } //mass
  }
  //=====================================================
  // deallocate the mass array and the space for lattice
  //=====================================================

  delete[] mass;
  sfree((void*)do_arg.start_conf_load_addr);


  // useless return code
  return 0;
}








