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
#include <alg/alg_smear.h>
#include <alg/alg_tcharge.h>
#include <singlenode/ReadLattice.h>
#include <util/command_line.h>
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


int main(int argc,char *argv[])
{
  CommandLine::is(argc,argv);
  
  ReadLattice rl(CommandLine::arg());

  rl.do_arg.beta=6.0;
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

  rot.run ();

  plaq.run();
  link_trace(lattice);

  rnd.free();


  return 0;
}
