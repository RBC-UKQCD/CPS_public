#include<config.h>

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#include<util/lattice.h>
#include<util/random.h>
#include<util/time_cps.h>
#include<util/command_line.h>

#include<alg/alg_hmc.h>
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

//#include <chroma.h>
#include <omp.h>
#include <pthread.h>
//-------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

const char *cname = "";

HmcArg hmc_arg;

ActionGaugeArg gauge_arg;

// define all integrators; vml file arguments will not change
IntABArg ab1_arg;
//IntABArg ab2_arg;
//IntABArg ab3_arg;

EvoArg evo_arg;
DoArg do_arg;
PbpArg pbp_arg;
NoArg no_arg;
ApeSmearArg ape_arg;

void checkpoint(int traj);

#define decode_vml(arg_name)  do{                                       \
        if ( ! arg_name.Decode(#arg_name".vml", #arg_name) )            \
            ERR.General(cname, fname, "Bad " #arg_name ".vml.\n");      \
    } while(0)

void decode_vml_all(void)
{
    const char *fname = "decode_vml_all()";

    decode_vml(do_arg);
    decode_vml(hmc_arg);
    decode_vml(evo_arg);
    decode_vml(gauge_arg);
    // decode_vml(quo_tm_arg);
//    decode_vml(quo_arg);
//    decode_vml(rat_quo_arg);
    decode_vml(ab1_arg);
//    decode_vml(ab2_arg);
    // decode_vml(ab3_arg);
//    decode_vml(pbp_arg);
      decode_vml(ape_arg);
}

void truncate_it(CommonArg *common_arg, const char stem[], int traj);
void measure_plaq(CommonArg &common_arg);
void measure_wline(CommonArg &common_arg);
void measure_tc(CommonArg &common_arg, int n_smear, int cycle);
//void measure_pbp(CommonArg &common_arg, int traj);
void run_hmc(CommonArg &common_arg, int traj, AlgIntAB &int_ab);


void setup(int argc, char *argv[])
{
    const char *fname = "setup()";

    Start(&argc, &argv);
    CommandLine::is(argc,argv);

    if(argc < 2) {
        ERR.General(cname, fname, "Must provide VML directory.\n");
    }

    if(chdir(CommandLine::arg()) != 0) {
        ERR.General(cname, fname, "Changing directory to %s failed.\n", argv[1]);
    }

    decode_vml_all();

    if(chdir(evo_arg.work_directory) != 0) {
        ERR.General(cname, fname, "Changing directory to %s failed.\n", evo_arg.work_directory);
    }
    VRB.Result(cname, fname, "Reading VML files successfully.\n");

    GJP.Initialize(do_arg);
    LRG.setSerial();
    LRG.Initialize();


    if ( do_arg.start_conf_kind != START_CONF_FILE ) {
        checkpoint(evo_arg.traj_start);
    }
}

int main(int argc, char *argv[])
{
    const char *fname = "main()";
    int if_noisy=1;

    setup(argc, argv);

if (if_noisy)
{   
    GwilsonFnone lat;
    lat.delta_beta = CommandLine::arg_as_Float();
    lat.deltaS_offset = CommandLine::arg_as_Float();
    int block[]={1,1,1,1};
    lat.SetSigmaBlock(block);
    VRB.Result("","main()","delta_beta=%g deltaS_offset=%g\n",lat.delta_beta ,  lat.deltaS_offset);
}

    //////////////////////////////////////////////////////////////////////
    // creating actions and integrators
    AlgMomentum mom;
    AlgActionGauge gauge(mom, gauge_arg);

    IntABArg sum_arg;
    sum_arg.A_steps = 1;
    sum_arg.B_steps = 1;
    sum_arg.level = TOP_LEVEL_INTEGRATOR;


    //!< Construct numerical integrators
    AlgIntAB &ab1 = AlgIntAB::Create(mom, gauge, ab1_arg);
    //////////////////////////////////////////////////////////////////////

    int traj = evo_arg.traj_start;
    for(int conf = 0; conf< evo_arg.gauge_configurations; ++conf) {
        CommonArg common_arg_plaq;
        CommonArg common_arg_pbp;
        CommonArg common_arg_hmc;
        CommonArg common_arg_top;
        CommonArg common_arg_wline;

        truncate_it(&common_arg_plaq , evo_arg.plaquette_stem      , traj);
        truncate_it(&common_arg_wline, "../results/alg_wline/wline", traj);
        truncate_it(&common_arg_pbp  , evo_arg.pbp_stem            , traj);
        truncate_it(&common_arg_hmc  , evo_arg.evo_stem            , traj);
//	    measure_tc(common_arg_top ,30,2);
if (if_noisy)
{
	    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, gauge_arg.gluon);
	    lat.SigmaHeatbath();
	    LatticeFactory::Destroy();
}

        // Inner trajectory loop
        for(int i = 0; i < evo_arg.gauge_unload_period; ++i, ++traj) {

            measure_plaq(common_arg_plaq);
//            measure_wline(common_arg_wline);
//            measure_pbp(common_arg_pbp, traj);

            run_hmc(common_arg_hmc, traj, ab1);
        }//End of inter-cfg sweep
        truncate_it(&common_arg_top  , "../results/alg_top/top"    , traj);
	measure_tc(common_arg_top ,30,2);

        checkpoint(traj);

    } //End config loop

    // AlgIntAB::Destroy(ab3);
//    AlgIntAB::Destroy(ab2);
    AlgIntAB::Destroy(ab1);

    End();

    VRB.Result(cname, fname, "Program ended normally.\n");
    
    return 0;
}

#undef encode_vml
#define encode_vml(arg_name, traj) do{                                  \
        char vml_file[256];                                             \
        sprintf(vml_file, #arg_name".%d", traj);                        \
        if( !arg_name.Encode(vml_file, #arg_name) ){                    \
            ERR.General(cname, fname, #arg_name " encoding failed.\n"); \
        }                                                               \
    }while(0)

void checkpoint(int traj)
{
    const char *fname="checkpoint()";

    char lat_file[256];
    char rng_file[256];

    Float time = -dclock();

    // Save this config to disk
    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);

    sprintf(lat_file,"%s.%d",evo_arg.gauge_file_stem,traj);
    QioArg wt_arg(lat_file,0.001);

    wt_arg.ConcurIONumber=evo_arg.io_concurrency;
    WriteLatticeParallel wl;
    wl.setHeader(evo_arg.ensemble_id,evo_arg.ensemble_label,traj);
    wl.write(lat,wt_arg);

    if(!wl.good())
        ERR.General(cname,fname,"Failed write lattice %s",lat_file);

    LatticeFactory::Destroy();

    // Save the RNG's
    sprintf(rng_file,"%s.%d",evo_arg.rng_file_stem,traj);
    if ( !LRG.Write(rng_file) )
        ERR.General(cname,fname,"Failed write RNG file %s",rng_file);

    // Update the parameter files for restart
    do_arg.start_seed_filename = rng_file;
    do_arg.start_seed_kind = START_SEED_FILE;
    do_arg.start_conf_filename = lat_file;
    do_arg.start_conf_kind = START_CONF_FILE;
    evo_arg.traj_start     = traj;

    encode_vml(hmc_arg, traj);
    encode_vml(gauge_arg, traj);
    // encode_vml(quo_tm_arg, traj);
//    encode_vml(quo_arg, traj);
//    encode_vml(rat_quo_arg, traj);
    encode_vml(ab1_arg, traj);
//    encode_vml(ab2_arg, traj);
    // encode_vml(ab3_arg, traj);
    encode_vml(pbp_arg, traj);
    encode_vml(do_arg, traj);
    encode_vml(evo_arg, traj);

    time += dclock();
    print_flops("","checkpoint()",0,time);
}

void truncate_it(CommonArg *common_arg, const char stem[], int traj)
{
    char fnbuf[1024];
    sprintf(fnbuf, "%s.%d", stem, traj);
    FILE *truncate_it = Fopen(fnbuf, "w");
    Fclose(truncate_it);
    common_arg->set_filename(fnbuf);
}

void measure_plaq(CommonArg &common_arg)
{
    const char *fname = "measure_plaq()";

    Float dtime = -dclock();

    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_WILSON);
    AlgPlaq plaq(lat, &common_arg, &no_arg);
    plaq.run();
    LatticeFactory::Destroy();

    dtime += dclock();
    print_flops("AlgPlaq", "run()", 0, dtime);	
}

void measure_wline(CommonArg &common_arg)
{
    const char *fname = "measure_wline()";

    Float dtime = -dclock();

    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_WILSON);
    AlgWline wline(lat, &common_arg, &no_arg);
    wline.run();
    LatticeFactory::Destroy();

    dtime += dclock();
    print_flops("AlgWline", "run()", 0, dtime);	
}

void measure_tc(CommonArg &common_arg, int n_smear, int cycle)
{
    const char *fname = "measure_tc()";

    Float dtime = -dclock();

    // calculate the topological charge. Need to copy the lattice since
    // we need to smear it first. Use Chulwoo's "lattice container"
    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_WILSON);
    LatticeContainer lat_cont;
    // copy lattice
    lat_cont.Get(lat);

    //----- mess up lattice -------------------------
    AlgApeSmear ape(lat, &common_arg, &ape_arg);
    AlgTcharge  tcharge(lat, &common_arg);
    for (int i = 0; i < cycle; ++i) {
        VRB.Result(cname,fname,"%i\n",i);
        VRB.Result(cname,fname,"   running tcharge\n"); tcharge.run();
	for(int j=0; j<n_smear; j++){
        VRB.Result(cname,fname,"   running ape\n"); ape.run();
        }
    }
    tcharge.run();
    // restore the lattice
    lat_cont.Set(lat);
    LatticeFactory::Destroy();

    dtime += dclock();
    print_flops("AlgTcharge", "run()", 0, dtime);
}


void run_hmc(CommonArg &common_arg, int traj, AlgIntAB &int_ab)
{
    const char *fname = "run_hmc()";

    Float dtime = -dclock();

    if ( (evo_arg.reproduce_interval > 0) &&
         (traj % evo_arg.reproduce_interval) == 0 ) {
        VRB.Result(cname,fname,"Running traj %d with reproduction\n",traj);
        hmc_arg.reproduce = REPRODUCE_YES;
    } else {
        VRB.Result(cname,fname,"Running traj %d without reproduction\n",traj);
        hmc_arg.reproduce = REPRODUCE_NO;	
    }
    
    AlgHmc hmc(int_ab, common_arg, hmc_arg);
    hmc.run();

    dtime += dclock();
    print_flops("AlgHmc", "run()", 0, dtime);
}
