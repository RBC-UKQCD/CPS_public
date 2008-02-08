///  $Id: main.C,v 1.13 2008-02-08 18:35:09 chulwoo Exp $
///  Demonstrate the random gauge transformation code.
///


/*!------------------------------

 This is an example program that shows how
 to perform a random gauge trasformation on the 
 gauge fields. A random gauge transformation is 
 done on a unit gauge configuration. The plaquette
 is computed and the link trace.

 The gauge transformation code was orginally written
 by Chris Dawson (using existing code from the BNL group).

*/

#include <config.h>
#include <util/qcdio.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/random.h>
#include <alg/alg_plaq.h>
#include <alg/do_arg.h>
#include <util/ReadLattice.h>
#include <alg/alg_rnd_gauge.h>
#include<util/testing_framework.h>
#include<util/dump_xml.h>
#include <comms/sysfunc_cps.h>

USING_NAMESPACE_CPS


static Float link_trace(Lattice& lattice)
{
#if TARGET==cpsMPI
  using MPISCU::printf;
  using MPISCU::fprintf;
#endif
    
  Float linktrace(0);
  int is;
  Matrix *m =  lattice.GaugeField(); 
  for(is=0;is< GJP.VolSites()*4; is++){
    linktrace += m->ReTr();
    m++;
  }
  linktrace /= (GJP.VolSites()*12.0);
  printf(" link trace %e\n",linktrace );


  return linktrace ; 
}

static void setup_do_arg(DoArg& do_arg) ; 

static const int nx = 4 ;
static const int ny = 4 ;
static const int nz = 4 ;
static const int nt = 4 ;


int main(int argc,char *argv[])
{

    DoArg do_arg;
    setup_do_arg(do_arg);

#if TARGET==cpsMPI
  MPISCU::set_pe_grid(do_arg.x_nodes, do_arg.y_nodes, do_arg.z_nodes, do_arg.t_nodes);
  using MPISCU::printf;
  using MPISCU::fprintf;
#endif

  printf("Demonstrate random gauge transform code\n") ;
  printf("node_sites[x,y,z,t] = [%d,%d,%d,%d]  \n",
	 do_arg.x_node_sites  , do_arg.y_node_sites , 
	 do_arg.z_node_sites  , do_arg.t_node_sites );
  printf("_nodes[x,y,z,t] = [%d,%d,%d,%d]\n",do_arg.x_nodes,
	 do_arg.y_nodes,
	 do_arg.z_nodes,
	 do_arg.t_nodes) ; 
  
  GJP.Initialize(do_arg);
  
  GwilsonFnone    lattice;
  CommonArg    common_arg;
  NoArg           no_arg ;
  common_arg.results = (void*)"out.dat";
  

  AlgPlaq plaq(lattice,&common_arg,&no_arg);

  plaq.run();
  Float tmp = link_trace(lattice); 

  AlgRandomGauge rnd( lattice, &common_arg );

  rnd.set_theta(1);
  rnd.run ();

  AlgRotateGauge rot( lattice, &common_arg );
  printf("Performed random gauge transform\n") ; 

  rot.run ();

  plaq.run();
  Float new_trace = link_trace(lattice);

  // write out as XML
  dump_xml output("candidate_t_gauge_rot.xml", "gauge") ;

  output.write(new_trace,"link_trace"); 

  output.close(); 

  //------------------------------ 
  // testing framework - compare the test against the old value
  const Float tol = 0.0001 ; 
  const Float old_trace = 0.6241053;
  compare_float_relative(new_trace,old_trace,tol) ; 
  // ----- end of testing framework -------------

  rnd.free();


  return 0;
}



void setup_do_arg(DoArg& do_arg)
{

  do_arg.x_node_sites = nx/SizeX();
  do_arg.y_node_sites = ny/SizeY();
  do_arg.z_node_sites = nz/SizeZ();
  do_arg.t_node_sites = nt/SizeT();
  do_arg.s_node_sites = 0;
  do_arg.x_nodes = SizeX();
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_PRD;
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.beta = 5.6;

}


