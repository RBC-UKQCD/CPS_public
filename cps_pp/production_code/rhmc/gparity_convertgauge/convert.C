//Read a standard no-Gparity lattice and save it as a G-parity lattice. This can be used for example to convert a configuration
//from an existing thermalized ensemble into a form usable by the G-parity code as a start point for an evolution.

#include<config.h>

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#include<util/lattice.h>
#include<util/lattice/fbfm.h>
#include<util/random.h>
#include<util/time_cps.h>

#include<alg/alg_hmc.h>
#include<alg/bfm_arg.h>
#include<alg/common_arg.h>
#include<alg/hmc_arg.h>
#include<alg/hmd_arg.h>

#include<alg/alg_int.h>
#include<alg/int_arg.h>
#include<alg/alg_wline.h>
#include<alg/no_arg.h>
#include<alg/do_arg.h>
#include<alg/alg_plaq.h>
#include<alg/alg_pbp.h>
#include<alg/alg_tcharge.h>
#include<alg/alg_smear.h>
#include<alg/ape_smear_arg.h>
#include<alg/pbp_arg.h>
#include<alg/alg_remez.h>
#include<alg/alg_w_spect.h>
#include<alg/array_arg.h>
#include<alg/alg_fix_gauge.h>

#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/qcdio.h>
#include<util/WriteLatticePar.h>
#include<util/ReadLatticePar.h>
#include<util/qioarg.h>

#include <util/lat_cont.h>

#include <chroma.h>
#include <omp.h>
#include <pthread.h>
//-------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

const char *cname = "";

DoArg do_arg_in;
DoArg do_arg_out;

#define decode_vml(arg_name)  do{                                       \
        if ( ! arg_name.Decode(#arg_name".vml", #arg_name) )            \
            ERR.General(cname, fname, "Bad " #arg_name ".vml.\n");      \
    } while(0)

void decode_vml_all(void)
{
  const char *fname = "decode_vml_all()";

  decode_vml(do_arg_in);
  decode_vml(do_arg_out);
}



void setup(int argc, char *argv[])
{
  const char *fname = "setup()";

  Start(&argc, &argv);

  if(argc < 2) {
    ERR.General(cname, fname, "Must provide VML directory.\n");
  }

  if(chdir(argv[1]) != 0) {
    ERR.General(cname, fname, "Changing directory to %s failed.\n", argv[1]);
  }

  decode_vml_all();

  VRB.Result(cname, fname, "Reading VML files successfully.\n");

  GJP.Initialize(do_arg_in);
  LRG.Initialize();

#define CHECKEQ(WHAT) if( do_arg_in.WHAT != do_arg_out.WHAT ) ERR.General(cname,fname,"Nodes and node sites must be equal")

  CHECKEQ(x_node_sites);
  CHECKEQ(y_node_sites);
  CHECKEQ(z_node_sites);
  CHECKEQ(t_node_sites);
  CHECKEQ(s_node_sites);
  CHECKEQ(x_nodes);
  CHECKEQ(y_nodes);
  CHECKEQ(z_nodes);
  CHECKEQ(t_nodes);
  CHECKEQ(s_nodes);
}

int main(int argc, char *argv[])
{
    const char *fname = "main()";

    setup(argc, argv);
    if(!UniqueID()){
      printf("Original boundary conditions:\n");
      for(int i=0;i<4;i++) printf("%s\n",BndCndType_map[(int)GJP.Bc(i)].name);
      printf("\n");
    }

    long lat_sz_nogp = GJP.VolNodeSites()*18*4;
    Float* gauge_nogp = (Float*)pmalloc(lat_sz_nogp * sizeof(Float));

    {
      Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE); //will load the lattice
      memcpy((void*)gauge_nogp, (void*)lat.GaugeField(), lat_sz_nogp *sizeof(Float));
      lat.FreeGauge(); //free memory and reset
      LatticeFactory::Destroy();
    }

    do_arg_out.start_conf_kind = START_CONF_ORD;
    GJP.Initialize(do_arg_out);

    if(!UniqueID()){
      printf("New boundary conditions:\n");
      for(int i=0;i<4;i++) printf("%s\n",BndCndType_map[(int)GJP.Bc(i)].name);
      printf("\n");
    }
    
    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
    memcpy((void*)lat.GaugeField(),(void*)gauge_nogp,lat_sz_nogp *sizeof(Float));

    lat.CopyConjMatrixField(lat.GaugeField(),4);
    
    {
      QioArg wt_arg(do_arg_out.start_conf_filename,0.001);
      
      wt_arg.ConcurIONumber=1;
      WriteLatticeParallel wl;
      wl.setHeader("converted","converted",0);
      wl.write(lat,wt_arg);
      
      if(!wl.good())
        ERR.General(cname,fname,"Failed write lattice %s",do_arg_out.start_conf_filename);
    }      
    LatticeFactory::Destroy();

    End();
    VRB.Result(cname, fname, "Program ended normally.\n");
    
    return 0;
}
