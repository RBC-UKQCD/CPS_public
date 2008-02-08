#include<config.h>
/*----------------------------------------------------------*/
/*!  A simple text program for the qcdio library.

  This requires a suitable UKQCD format dataset, which it 
  will then load.

  A.N.Jackson: ajackson@epcc.ed.ac.uk                      
  -----------------------------------------------------------
  CVS keywords
 
  $Author: chulwoo $
  $Date: 2008-02-08 18:35:08 $
  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/qcdio/main.C,v 1.8 2008-02-08 18:35:08 chulwoo Exp $
  $Id: main.C,v 1.8 2008-02-08 18:35:08 chulwoo Exp $
  $Name: not supported by cvs2svn $
  $Locker:  $
  $RCSfile: main.C,v $
  $Revision: 1.8 $
  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/qcdio/main.C,v $
  $State: Exp $  */
/*----------------------------------------------------------*/


#include <util/qcdio.h>
#include <stdlib.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/do_arg.h>
#include<alg/common_arg.h>
#include<util/qcdio.h>
#ifdef PARALLEL
#include<comms/sysfunc_cps.h>
#endif

CPS_START_NAMESPACE
void qcdio_create_dummy_data( void ) {
  char* fprefix = "D52C202K3500U007000";
  char* fpar = "D52C202K3500U007000.par", fname[200], number[200];
  int latX = 16, latY = 16, latZ = 16, latT = 32;
  int t, z, y, x, dir, mu, nu, comp, bytes_out;
  double num;
  FILE* pfp;
  FILE* bfp;

// Make the .par parameter file:
  pfp = Fopen( fpar, "w" );
  Fprintf(pfp,"beta 0.51999998E+01\n");
  Fprintf(pfp,"latt[X] %i\n",latX);
  Fprintf(pfp,"latt[Y] %i\n",latY);
  Fprintf(pfp,"latt[Z] %i\n",latZ);
  Fprintf(pfp,"latt[T] %i\n",latT);
// The following items are members of a full .par file, but are not
// currently used in the CPS:
/*
start_type                    0
start_sweep                7000
rng_seed                      0
bc[X]                      periodic
bc[Y]                      periodic
bc[Z]                      periodic
bc[T]                      antiperiodic
calc_par                      0
kappa                      0.13500001E+00
clover                     0.20171001E+01
beta_shift                 0.00000000E+00
kappa_shift                0.00000000E+00
clover_shift               0.00000000E+00
mixing_parameter           0.10000000E+01
timestep                   0.78125000E-02
gauge_updates                 2
sweep_length                 32
num_sweeps                  905
save_int                    100
checkpoint_int                1
next_checkpoint_num           1
checkpoint_cycle_len          2
time_limit_in_mins          700
load_momenta               yes
solver_type                bicgstab
restart_solver_type        cg
guess_strategy_x           zero
guess_strategy_y           zero
min_iter                     10
max_iter                   1000
restarts                      2
print_freq                    0
omega                      0.11000000E+01
target_residue             0.10000000E-06
relax_md_residue_by        0.10000000E+01
precondition               yes
validate_plaquettes        yes
validate_tplaquettes       yes
validate_tcsum             yes
input_gauge_size              4
byte_swap_input_gauge      no
output_gauge_size             4
byte_swap_output_gauge     no
output_gauge_path          /ework2/c1/NF2/BETA52/CLOVER202/V16X32/KAPPA3500/gauge/
output_fe_path             /ework2/c1/NF2/BETA52/CLOVER202/V16X32/KAPPA3500/fe/
output_mom_path            /ework2/c1/NF2/BETA52/CLOVER202/V16X32/KAPPA3500/momentum/
input_gauge_path           /ework2/c1/NF2/BETA52/CLOVER202/V16X32/KAPPA3500/gauge/
input_fe_path              /ework2/c1/NF2/BETA52/CLOVER202/V16X32/KAPPA3500/fe/
input_mom_path             /ework2/c1/NF2/BETA52/CLOVER202/V16X32/KAPPA3500/momentum/
checkpoint_path            /ework2/c1/NF2/BETA52/CLOVER202/V16X32/KAPPA3500/checkpoint/
input_stem                 D52C202K3500U007000
output_prefix              D52C202K3500
*/
  Fclose(pfp);

// Now write out each timeslice file in turn, filling every file with simple
// unique numbers for every position:
  num = 0.0;
  for( t = 0; t < latT; t++ ) {
    // Open the file:
    sprintf(fname,"%sT%02i",fprefix,t);
    printf("Writing %s...\n",fname);
    bfp = Fopen( fname, "wb" );
    for( z = 0; z < latZ; z++ ) {
    for( y = 0; y < latY; y++ ) {
    for( x = 0; x < latX; x++ ) {
    for( dir = 0; dir < 4; dir++ ) {
    for( mu = 0; mu < 3; mu++ ) {
    for( nu = 0; nu < 2; nu++ ) {
    for( comp = 0; comp < 2; comp++ ) {
      sprintf(number,"%i%i%i%i.%i%i%i%i",x,y,z,t,dir,mu,nu,comp);
      sscanf(number,"%lf",&num);
      bytes_out = fwrite((char*)&num, 1, sizeof(double), bfp );
      if( bytes_out != sizeof(double) ) {
        printf("Failed while writing to file %s!\n",fname);
	exit(1);
      }
    }}}}}}}

    Fclose(bfp);
  }
}
CPS_END_NAMESPACE

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
USING_NAMESPACE_CPS

int main( int argc, char** argv ) {
  DoArg do_arg;



  do_arg.x_node_sites = 8;
  do_arg.y_node_sites = 8;
  do_arg.z_node_sites = 8;
  do_arg.t_node_sites = 16;
  do_arg.s_node_sites = 0;

  do_arg.x_nodes = 2;
  do_arg.y_nodes = 2;
  do_arg.z_nodes = 2;
  do_arg.t_nodes = 2;
  do_arg.s_nodes = 1;
#else
  do_arg.x_node_sites = 16;
  do_arg.y_node_sites = 16;
  do_arg.z_node_sites = 16;
  do_arg.t_node_sites = 32;
  do_arg.s_node_sites = 0;
  
  do_arg.x_nodes = 1;
  do_arg.y_nodes = 1;
  do_arg.z_nodes = 1;
  do_arg.t_nodes = 1;
  do_arg.s_nodes = 1;
#endif

  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_DISORD;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.beta = 6.0;
  do_arg.dwf_height = 0.9;


  GJP.Initialize(do_arg);
  
//--
  
  GwilsonFwilson lat;

  // Create the dummy data files:
  if( UniqueID() == 0 ) {
    qcdio_create_dummy_data();
  }
  
  // Avoid normalizing the SU(3) matricies as this will break this test:
  qcdio_set_normalize(0); 
  
  // Load the array, of doubles, no byte-swap and no transpose:
  printf("Loading the data...\n");
  qload_gauge("D52C202K3500U007000", lat, sizeof(double), 0, 0);

  // Check the data has been loaded up correctly...
  printf("Checking the data...\n");
  int site[4]; 
  const Matrix *link;
  Complex num;
  double trunum, latnum;
  char number[200];
 
  for( site[3] = 0; site[3] < GJP.TnodeSites(); site[3]++ ) {

   for( site[0] = 0; site[0] < GJP.XnodeSites(); site[0]++ ) {
   for( site[1] = 0; site[1] < GJP.YnodeSites(); site[1]++ ) {
   for( site[2] = 0; site[2] < GJP.ZnodeSites(); site[2]++ ) {

    for( int dirn = 0; dirn < 4; dirn ++ ) { 
    link = lat.GetBufferedLink(site,dirn);
    for( int i = 0; i < 3; i++ ) {
      for( int j = 0; j < 2; j++ ) {
       num = (*link)(i,j);
       for( int comp = 0; comp < 2; comp++ ) {
        // Create the composite unique ID for this point in the lattice:
	sprintf(number,"%i%i%i%i.%i%i%i%i",
	site[0] + GJP.XnodeSites()*GJP.XnodeCoor(),
	site[1] + GJP.YnodeSites()*GJP.YnodeCoor(),
	site[2] + GJP.ZnodeSites()*GJP.ZnodeCoor(),
	site[3] + GJP.TnodeSites()*GJP.TnodeCoor(),
	dirn,i,j,comp);
        sscanf(number,"%lf",&trunum);
	// Choose the lattice number, either the complex or real part:
        if( comp == 0 ) {
	  latnum = num.real();
	} else {
	  latnum = num.imag();
	}
        // Check the numbers agree, and complain if not:
	if( latnum != trunum ) { 
          qprintf_allid(" Error! T[%i,%i,%i,%i,%i,%i,%i] %f,%f\n", 
	  site[0] + GJP.XnodeSites()*GJP.XnodeCoor(),
	  site[1] + GJP.YnodeSites()*GJP.YnodeCoor(),
	  site[2] + GJP.ZnodeSites()*GJP.ZnodeCoor(),
	  site[3] + GJP.TnodeSites()*GJP.TnodeCoor(),
	  dirn,
	  i,
	  j,
	  num.real(), num.imag() );
	  exit(1);
	}
       }
      }
    }

   }
 
  }}}

  }

  // Save the configuration out again:
  printf("Saving the data...\n");
  qsave_gauge("out.D52C202K3500U007000", lat, sizeof(double), 0, 0);

#if INCLUDE_MPI_SCU == 1
  SCUCommsFinalize(); /* Not present on QCDSP */
#endif
  
  printf("Test complete.\n");

  return(0);
}

