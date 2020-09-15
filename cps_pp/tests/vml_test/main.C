#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <util/qcdio.h>
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
#include <comms/scu.h>
#include <comms/glb.h>

#include <util/lattice.h>
#include <util/time_cps.h>
#include <alg/do_arg.h>
#include <alg/no_arg.h>
#include <alg/common_arg.h>
#include <alg/hmd_arg.h>
#include <alg/alg_plaq.h>
#include <alg/alg_rnd_gauge.h>
#include <alg/threept_arg.h>
#include <alg/threept_prop_arg.h>
#include <alg/alg_threept.h>
#include <util/smalloc.h>
#if(0==1)
 #include <ReadLattice.h>
 #include <WriteLattice.h>
#endif

#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>

#include <util/command_line.h>
#include <sstream>

#include<unistd.h>
#include<config.h>

#include <alg/qpropw.h>
#include <alg/qpropw_arg.h>

#include <alg/alg_fix_gauge.h>
#include <alg/fix_gauge_arg.h>

#include <util/data_shift.h>

#include <alg/prop_attribute_arg.h>
#include <alg/gparity_contract_arg.h>
#include <alg/propmanager.h>
#include <alg/alg_gparitycontract.h>
#include <alg/prop_dft.h>

using namespace std;
USING_NAMESPACE_CPS



int main(int argc,char *argv[])
{
  Start(&argc,&argv);
  CommandLine::is(argc,argv);

  if(argc!=4){ printf("Usage: ./executable <type> <filename> <classname>\n"); exit(0); }
  
#define TEST(TYPE)\
  {\
    TYPE arg; \
    if(!arg.Decode(argv[2],argv[3])){ \
      ERR.General("","main()","Failed to decode %s\n",argv[3]); exit(-1); \
    }else{ printf("Pass\n"); exit(0); } \
  }
  string type(argv[1]);
  if(type == string("GparityContractArg") ){
    TEST(GparityContractArg);
  }else{
    printf("Unknown type, please add\n"); exit(-1);
  }
  return 0;
}




