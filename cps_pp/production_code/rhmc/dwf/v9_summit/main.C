#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<config.h>
#include<util/gjp.h>
#include<util/timer.h>
#include<util/lattice.h>

#include<alg/hmc_arg.h>
#include<alg/hmd_arg.h>
#include<alg/no_arg.h>
#include<alg/alg_hmc.h>
#include<alg/alg_pbp.h>
#include<alg/alg_plaq.h>


USING_NAMESPACE_CPS
using namespace std;

const char *cname = "";

HmcArg hmc_arg;

ActionGaugeArg gauge_arg;
ActionQuotientArg quo_arg;
//ActionRationalQuotientArg rat_quo_arg;
ActionEOFAArg             eofa_arg_s;

// define all integrators; vml file arguments will not change
IntABArg ab1_arg;
IntABArg ab2_arg;
//IntABArg ab3_arg;

EvoArg evo_arg;
DoArg do_arg;
DoArgExt doext_arg;
PbpArg pbp_arg;
NoArg no_arg;
//ApeSmearArg ape_arg;

void checkpoint(int traj);

#undef encode_vml
#define encode_vml(arg_name, traj) do{                                  \
        char vml_file[256];                                             \
        sprintf(vml_file, #arg_name".%d", traj);                        \
        if( !arg_name.Encode(vml_file, #arg_name) ){                    \
            ERR.General(cname, fname, #arg_name " encoding failed.\n"); \
        }                                                               \
    }while(0)

#undef decode_vml
#define decode_vml(arg_name)  do{                                       \
        if ( ! arg_name.Decode(#arg_name".vml", #arg_name) )            \
            ERR.General(cname, fname, "Bad " #arg_name ".vml.\n");      \
    } while(0)

void decode_vml_all(void)
{
    const char *fname = "decode_vml_all()";

    decode_vml(do_arg);
    encode_vml(do_arg,0);
    decode_vml(doext_arg);
    encode_vml(doext_arg,0);
    decode_vml(hmc_arg);
    encode_vml(hmc_arg,0);
    decode_vml(evo_arg);
    encode_vml(evo_arg,0);
    decode_vml(gauge_arg);
    encode_vml(gauge_arg,0);
    decode_vml(quo_arg);
    encode_vml(quo_arg,0);
    decode_vml(eofa_arg_s);
    encode_vml(eofa_arg_s,0);
    decode_vml(ab1_arg);
    encode_vml(ab1_arg,0);
    decode_vml(ab2_arg);
    encode_vml(ab2_arg,0);
    decode_vml(pbp_arg);
    encode_vml(pbp_arg,0);
}

void truncate_it(CommonArg *common_arg, const char stem[], int traj);
void measure_plaq(CommonArg &common_arg);
void measure_pbp(CommonArg &common_arg, int traj);
void run_hmc(CommonArg &common_arg, int traj, AlgIntAB &int_ab);

void init(int *argc, char **argv[]){
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
#ifdef USE_QUDA
   if ( !QudaParam.Decode("quda_arg.vml","QudaParam") )
        { printf("Bum quda_arg\n"); exit(-1);}
#endif

    if(chdir(evo_arg.work_directory) != 0) {
        ERR.General(cname, fname, "Changing directory to %s failed.\n", evo_arg.work_directory);
    }
    VRB.Result(cname, fname, "Reading VML files successfully.\n");

    GJP.Initialize(do_arg);
    GJP.InitializeExt(doext_arg);
    GJP.ZMobius_PC_Type(ZMOB_PC_SYM1 );
    LRG.setSerial();
    LRG.Initialize();
 int threads;
 char * nthr_str = getenv("OMP_NUM_THREADS");
 if(nthr_str) sscanf(nthr_str,"%d",&threads);
 if(!UniqueID()) printf("nthreads=%d\n",threads);
 GJP.SetNthreads(threads);
 printf("GJP.Nthreads()=%d\n",GJP.Nthreads());

    init(&argc, &argv);

     if ( do_arg.start_conf_kind != START_CONF_FILE ) {
         checkpoint(evo_arg.traj_start);
     }
}

int main(int argc, char *argv[])
{
    const char *fname = "main()";

    setup(argc, argv);

    //////////////////////////////////////////////////////////////////////
    // creating actions and integrators
    AlgMomentum mom;
    AlgActionGauge gauge(mom, gauge_arg);
    AlgActionQuotient quotient(mom, quo_arg);
    AlgActionEOFA     s_quark(mom, eofa_arg_s,false);

    IntABArg sum_arg;
    sum_arg.A_steps = 1;
    sum_arg.B_steps = 1;
    sum_arg.level = EMBEDDED_INTEGRATOR;
    AlgIntSum sum(quotient, s_quark, sum_arg);

    //!< Construct numerical integrators
    AlgIntAB &ab1 = AlgIntAB::Create(mom, gauge, ab1_arg);
    AlgIntAB &ab2 = AlgIntAB::Create(ab1, sum,   ab2_arg);
    //////////////////////////////////////////////////////////////////////

    int traj = evo_arg.traj_start;
    for(int conf = 0; conf< evo_arg.gauge_configurations; ++conf) {
        CommonArg common_arg_plaq;
        CommonArg common_arg_pbp;
        CommonArg common_arg_hmc;
        // CommonArg common_arg_top;
        CommonArg common_arg_wline;

        truncate_it(&common_arg_plaq , evo_arg.plaquette_stem      , traj);
        truncate_it(&common_arg_pbp  , evo_arg.pbp_stem            , traj);
        truncate_it(&common_arg_hmc  , evo_arg.evo_stem            , traj);

        // Inner trajectory loop
        for(int i = 0; i < evo_arg.gauge_unload_period; ++i, ++traj) {

            measure_plaq(common_arg_plaq);
            measure_pbp(common_arg_pbp, traj);

            run_hmc(common_arg_hmc, traj, ab2);
            Timer::display("End of HMC");
        }//End of inter-cfg sweep

        checkpoint(traj);

    } //End config loop

    // AlgIntAB::Destroy(ab3);


    AlgIntAB::Destroy(ab2);
    AlgIntAB::Destroy(ab1);
    VRB.Result(cname, fname, "Program ended normally.\n");
    
    End();
    return 0;
}

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
    encode_vml(quo_arg, traj);
    encode_vml(eofa_arg_s, traj);
    encode_vml(ab1_arg, traj);
    encode_vml(ab2_arg, traj);
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

void measure_pbp(CommonArg &common_arg, int traj)
{
    const char *fname = "measure_pbp()";

    // fix pbp_arg
    pbp_arg.src_u_s = 0;
    pbp_arg.src_l_s = GJP.SnodeSites() * GJP.Snodes() - 1;
    pbp_arg.snk_u_s = GJP.SnodeSites() * GJP.Snodes() - 1;
    pbp_arg.snk_l_s = 0;

    const int g_int = evo_arg.gauge_unload_period;
    if (traj % g_int == 0 && evo_arg.measure_pbp) {
        Float dtime = -dclock();

        LRGState rng_state;
        rng_state.GetStates();

        Lattice &lat = LatticeFactory::Create(F_CLASS_MOBIUS, G_CLASS_NONE);
        AlgPbp pbp(lat, &common_arg, &pbp_arg);

        for(int pbp_counter = 0; pbp_counter < evo_arg.measure_pbp; pbp_counter++) {
            pbp.run();
        }
        LatticeFactory::Destroy();
        rng_state.SetStates();

        dtime += dclock();
        print_flops("AlgPbp", "run()", 0, dtime);	
    }
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
