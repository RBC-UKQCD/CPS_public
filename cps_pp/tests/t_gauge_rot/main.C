//
// This is an example program that shows how
// to perform a random gauge trasformation on the 
// gauge fields. A random gauge transformation is 
// done on a unit gauge configuration. The plaquette
// is computed and the link trace.
//
//
// The gauge transformation code was orginally written
// by Chris Dawson (using existing code from the BNL group).
//

#include <config.h>
#include <stdio.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/random.h>
#include <alg/alg_plaq.h>
#include <alg/do_arg.h>
#include <alg/common_arg.h>
//#include <alg/alg_smear.h>
//#include <alg/alg_tcharge.h>
#include <util/ReadLattice.h>
//#include <util/command_line.h>
#include <alg/alg_rnd_gauge.h>


CPS_START_NAMESPACE
GlobalJobParameter GJP;
LatRanGen          LRG;
Verbose            VRB;
Error              ERR;
CPS_END_NAMESPACE

USING_NAMESPACE_CPS


void link_trace(Lattice& lattice)
{
  Float linktrace(0);
  int is;
  Matrix *m =  lattice.GaugeField(); 
  for(is=0;is< GJP.VolSites()*4; is++){
    linktrace += m->ReTr();
    m++;
  }
  linktrace /= (GJP.VolSites()*12.0);
  printf(" link trace %e\n",linktrace );
}

void setup_do_arg(DoArg& do_arg) ; 

const int nx = 4 ;
const int ny = 4 ;
const int nz = 4 ;
const int nt = 4 ;


int main(int argc,char *argv[])
{
  printf("Demonstrate random gauge transform code\n") ; 
  //  CommandLine::is(argc,argv);
  
  //  ReadLattice rl(CommandLine::arg());
  ReadLattice rl("u_TEST_4223.101") ; 


  
  setup_do_arg(rl.do_arg) ; 

  GJP.Initialize(rl.do_arg);
  
  GwilsonFnone    lattice;
  CommonArg    common_arg;
  NoArg           no_arg ;
  common_arg.results = (void*)"out.dat";
  

  AlgPlaq plaq(lattice,&common_arg,&no_arg);

  plaq.run();
  link_trace(lattice);

  AlgRandomGauge rnd( lattice, &common_arg );

  rnd.set_theta(1);
  rnd.run ();

  AlgRotateGauge rot( lattice, &common_arg );
  printf("Performed random gauge transform\n") ; 

  rot.run ();

  plaq.run();
  link_trace(lattice);

  rnd.free();


  return 0;
}



void setup_do_arg(DoArg& do_arg)
{


#ifdef PARALLEL
#if 1
  do_arg.x_node_sites = nx/SizeX();
  do_arg.y_node_sites = ny/SizeY();
  do_arg.z_node_sites = nz/SizeZ();
  do_arg.t_node_sites = nt/SizeT();
#else
  do_arg.x_node_sites = 2;
  do_arg.y_node_sites = 2;
  do_arg.z_node_sites = 2;
  do_arg.t_node_sites = 2;
#endif
  do_arg.s_node_sites = 0;
  do_arg.x_nodes = SizeX();
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();
  do_arg.s_nodes = 1;
#else
  do_arg.x_node_sites = nx ;
  do_arg.y_node_sites = ny ;
  do_arg.z_node_sites = nz ;
  do_arg.t_node_sites = nt ;
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
  do_arg.t_bc = BND_CND_PRD;
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.beta = 5.6;

  printf("node_sites[x,y,z,t] = [%d,%d,%d,%d]  \n",
	 do_arg.x_node_sites  , do_arg.y_node_sites , 
	 do_arg.z_node_sites  , do_arg.t_node_sites );

  printf("_nodes[x,y,z,t] = [%d,%d,%d,%d]\n",do_arg.x_nodes,
	 do_arg.y_nodes,
	 do_arg.z_nodes,
	 do_arg.t_nodes) ; 

}


