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


HmcArg hmc_arg;

ActionGaugeArg gauge_arg;
ActionQuotientArg quo_arg;

//do the twisted mass term with two RHMCs
ActionRationalQuotientArg rat_quo_tm_1_arg;
ActionRationalQuotientArg rat_quo_tm_2_arg;

ActionRationalQuotientArg rat_quo_arg;

// define all integrators; vml file arguments will not change
IntABArg ab1_arg;
IntABArg ab2_arg;
IntABArg ab3_arg;

//BfmArg b_arg;
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
    char *fname = "decode_vml_all()";

    //decode_vml(b_arg);
    decode_vml(do_arg);
    decode_vml(hmc_arg);
    decode_vml(evo_arg);
    decode_vml(gauge_arg);
    decode_vml(rat_quo_tm_1_arg);
    decode_vml(rat_quo_tm_2_arg);
    decode_vml(quo_arg);
    decode_vml(rat_quo_arg);
    decode_vml(ab1_arg);
    decode_vml(ab2_arg);
    decode_vml(ab3_arg);
    decode_vml(pbp_arg);
    // decode_vml(ape_arg);
}

void truncate_it(CommonArg *common_arg, const char stem[], int traj);
void measure_plaq(CommonArg &common_arg);
void measure_wline(CommonArg &common_arg);
void measure_tc(CommonArg &common_arg, int cycle);
void measure_pbp(CommonArg &common_arg, int traj);
void run_hmc(CommonArg &common_arg, int traj, AlgIntAB &int_ab);
void config_repro(const char* repro_config, const Float &tolerance = 1e-09);

void init_bfm(int *argc, char **argv[])
{
    Chroma::initialize(argc, argv);
    multi1d<int> nrow(Nd);

    for(int i = 0; i< Nd; ++i)
        nrow[i] = GJP.Sites(i);

    Layout::setLattSize(nrow);
    Layout::create();

    // Fbfm::bfm_args[0].solver = HtCayleyTanh;
    // Fbfm::bfm_args[0].precon_5d = 0;
    // Fbfm::bfm_args[0].solver = DWF;
    // Fbfm::bfm_args[0].precon_5d = 1;
    Fbfm::bfm_args[0].solver = HmCayleyTanh;
    Fbfm::bfm_args[0].precon_5d = 0;

    printf("Fbfm::bfm_args[0].solver = %d\n",Fbfm::bfm_args[0].solver);

    //VRB.Result( "cps", "init_bfm", "-> getenv(\"use_mixed_solver\")\n") ;

    // mixed-precision CG *based on environment variable*, *true by default*
    char* use_mixed_solver_env = getenv( "use_mixed_solver" );
    Fbfm::use_mixed_solver = true;
    if ( use_mixed_solver_env && strcmp( use_mixed_solver_env, "false" ) == 0 )
      Fbfm::use_mixed_solver = false;
    VRB.Result( "cps", "init_bfm", "Fbfm::use_mixed_solver: %d\n", Fbfm::use_mixed_solver );

    Fbfm::bfm_args[0].Ls = GJP.SnodeSites();
    Fbfm::bfm_args[0].M5 = GJP.DwfHeight();
    Fbfm::bfm_args[0].mass = 0.1;
    Fbfm::bfm_args[0].twistedmass = 0.0;
    Fbfm::bfm_args[0].residual = 1e-8;
    Fbfm::bfm_args[0].max_iter = 10000;
    Fbfm::bfm_args[0].Csw = 0.0;
    Fbfm::bfm_args[0].time_report_iter = -16;

    Fbfm::bfm_args[0].node_latt[0] = QDP::Layout::subgridLattSize()[0];
    Fbfm::bfm_args[0].node_latt[1] = QDP::Layout::subgridLattSize()[1];
    Fbfm::bfm_args[0].node_latt[2] = QDP::Layout::subgridLattSize()[2];
    Fbfm::bfm_args[0].node_latt[3] = QDP::Layout::subgridLattSize()[3];

    multi1d<int> procs = QDP::Layout::logicalSize();

    Fbfm::bfm_args[0].local_comm[0] = procs[0] > 1 ? 0 : 1;
    Fbfm::bfm_args[0].local_comm[1] = procs[1] > 1 ? 0 : 1;
    Fbfm::bfm_args[0].local_comm[2] = procs[2] > 1 ? 0 : 1;
    Fbfm::bfm_args[0].local_comm[3] = procs[3] > 1 ? 0 : 1;

    Fbfm::bfm_args[0].ncoor[0] = 0;
    Fbfm::bfm_args[0].ncoor[1] = 0;
    Fbfm::bfm_args[0].ncoor[2] = 0;
    Fbfm::bfm_args[0].ncoor[3] = 0;

    Fbfm::current_arg_idx = 0;

    // mobius_scale = b + c in Andrew's notation
    bfmarg::mobius_scale = 2.0; //b_arg.mobius_scale;
#if TARGET == BGQ
    bfmarg::Threads(32);
#else
    bfmarg::Threads(1);
#endif 
    omp_set_num_threads(bfmarg::threads);

    bfmarg::Reproduce(0);
    bfmarg::ReproduceChecksum(0);
    bfmarg::ReproduceMasterCheck(0);
    bfmarg::Verbose(1);
}



//By setting quo_tm_arg.bi_arg.fermion = F_CLASS_BFM_TYPE2 we enable using Fbfm to do the Twisted Mass fermions also
void init_bfm_wilsontm(){
  if(!UniqueID()) printf("Doing DSDR term with Fbfm also\n");

  Fbfm::bfm_args[1].solver = WilsonTM;
  Fbfm::bfm_args[1].precon_5d = 0;

  Fbfm::bfm_args[1].Ls = 1;
  Fbfm::bfm_args[1].M5 = 0;
  Fbfm::bfm_args[1].mass = 0.1; //These are overwriten with the correct values within the evolution code
  Fbfm::bfm_args[1].twistedmass = 0.1;
  Fbfm::bfm_args[1].residual = 1e-8;
  Fbfm::bfm_args[1].max_iter = 10000;
  Fbfm::bfm_args[1].Csw = 0.0;
  Fbfm::bfm_args[1].time_report_iter = -16;

  Fbfm::bfm_args[1].node_latt[0] = QDP::Layout::subgridLattSize()[0];
  Fbfm::bfm_args[1].node_latt[1] = QDP::Layout::subgridLattSize()[1];
  Fbfm::bfm_args[1].node_latt[2] = QDP::Layout::subgridLattSize()[2];
  Fbfm::bfm_args[1].node_latt[3] = QDP::Layout::subgridLattSize()[3];

  multi1d<int> procs = QDP::Layout::logicalSize();

  Fbfm::bfm_args[1].local_comm[0] = procs[0] > 1 ? 0 : 1;
  Fbfm::bfm_args[1].local_comm[1] = procs[1] > 1 ? 0 : 1;
  Fbfm::bfm_args[1].local_comm[2] = procs[2] > 1 ? 0 : 1;
  Fbfm::bfm_args[1].local_comm[3] = procs[3] > 1 ? 0 : 1;

  Fbfm::bfm_args[1].ncoor[0] = 0;
  Fbfm::bfm_args[1].ncoor[1] = 0;
  Fbfm::bfm_args[1].ncoor[2] = 0;
  Fbfm::bfm_args[1].ncoor[3] = 0;
}






void setup(int *argc, char*** argv){
  const char *fname = "setup()";

  Start(argc, argv);

  if(*argc < 2) {
    ERR.General(cname, fname, "Must provide VML directory.\n");
  }

  if(chdir( (*argv)[1]) != 0) {
    ERR.General(cname, fname, "Changing directory to %s failed.\n", argv[1]);
  }

  decode_vml_all();

  if(chdir(evo_arg.work_directory) != 0) {
    ERR.General(cname, fname, "Changing directory to %s failed.\n", evo_arg.work_directory);
  }
  VRB.Result(cname, fname, "Reading VML files successfully.\n");

#if TARGET == BGQ
  LRG.setSerial();
#endif

  GJP.Initialize(do_arg);
  LRG.Initialize();

  VRB.Result(cname, fname, "VRB.Level(%d)\n", do_arg.verbose_level);
  VRB.Level(do_arg.verbose_level);

  init_bfm(argc, argv);
  VRB.Result(cname, fname, "Mobius scale (2c+1) = %.10f\n", bfmarg::mobius_scale);

  //Use Fbfm to do the Twisted Mass Wilson part also
  if(rat_quo_tm_1_arg.bi_arg.fermion == F_CLASS_BFM_TYPE2) ERR.General(cname,fname,"BFM wilsonTM disabled in this executable");
  else if(rat_quo_tm_1_arg.bi_arg.fermion != F_CLASS_WILSON_TM) ERR.General(cname,fname,"quo_tm_arg fermion type can only be either F_CLASS_WILSON_TM or F_CLASS_BFM_TYPE2");
}


int main(int argc, char *argv[])
{
    const char *fname = "main()";

    setup(&argc, &argv);

    //CK: Option to run through 'gauge_unload_period' trajectories then compare the resulting gauge field
    //    with some other loaded field. Job will continue after comparison returns true
    char* repro_config;
    bool do_config_repro = false;
    Float repro_tolerance = 1e-09;

    int i=1;
    while(i<argc){
      char* cmd = argv[i];  
      if( strncmp(cmd,"-repro_config",20) == 0){
	if(i==argc-1){ printf("-repro_config requires an argument\n"); exit(-1); }
	repro_config = argv[i+1];
	if(!UniqueID()) printf("Running reproduction test against config %s at end of first unload period\n",repro_config);
	do_config_repro = true;
	i+=2;
      }else if( strncmp(cmd,"-repro_tolerance",20) == 0){
	if(i==argc-1){ printf("-repro_tolerance requires an argument\n"); exit(-1); }
	repro_tolerance = atof(argv[i+1]);
	if(repro_tolerance == 0.0) ERR.General(cname,fname,"Either you set repro_tolerance to 0.0 or the str->double conversion failed, either way that's an invalid input!");
	if(!UniqueID()) printf("Set reproduction tolerance to %e\n",repro_tolerance);
	i+=2;
      }else i++;
    }

    //VRB.ElapsedTime("CPS", fname);

    // OpenMP test
    VRB.Result( "CPS", "main", "omp_get_num_threads[1] -> %d\n", omp_get_num_threads() );
    
    #pragma omp parallel
    {
      if ( UniqueID() == 0 && omp_get_thread_num() == 0 ) {
        VRB.Result( "CPS", "main", "omp_get_num_threads[2] -> %d\n", omp_get_num_threads() );
      }
    }

    VRB.Result( "CPS", "main", "omp_get_num_threads[3] -> %d\n", omp_get_num_threads() );

    //////////////////////////////////////////////////////////////////////
    // creating actions and integrators
    AlgMomentum mom;
    AlgActionGauge gauge(mom, gauge_arg);
    AlgActionRationalQuotient rat_quo_tm_1(mom, rat_quo_tm_1_arg);
    AlgActionRationalQuotient rat_quo_tm_2(mom, rat_quo_tm_2_arg);
    AlgActionQuotient quotient(mom, quo_arg);
    AlgActionRationalQuotient rat_quo(mom, rat_quo_arg);

    IntABArg sum_arg;
    sum_arg.A_steps = 1;
    sum_arg.B_steps = 1;
    sum_arg.level = EMBEDDED_INTEGRATOR;
    AlgIntSum sum(quotient, rat_quo, sum_arg);

    AlgIntSum sum_tm(rat_quo_tm_1,rat_quo_tm_2, sum_arg);


    //!< Construct numerical integrators
    AlgIntAB &ab1 = AlgIntAB::Create(mom, gauge,       ab1_arg);
    AlgIntAB &ab2 = AlgIntAB::Create(ab1, sum_tm, ab2_arg);
    AlgIntAB &ab3 = AlgIntAB::Create(ab2, sum,         ab3_arg);
    //////////////////////////////////////////////////////////////////////

    int traj = evo_arg.traj_start;
    for(int conf = 0; conf< evo_arg.gauge_configurations; ++conf) {
        CommonArg common_arg_plaq;
        CommonArg common_arg_pbp;
        CommonArg common_arg_hmc;
        // CommonArg common_arg_top;
        // CommonArg common_arg_wline;

        truncate_it(&common_arg_plaq , evo_arg.plaquette_stem      , traj);
        // truncate_it(&common_arg_top  , "../results/alg_top/top"    , traj);
        // truncate_it(&common_arg_wline, "../results/alg_wline/wline", traj);
        truncate_it(&common_arg_pbp  , evo_arg.pbp_stem            , traj);
        truncate_it(&common_arg_hmc  , evo_arg.evo_stem            , traj);

        // Inner trajectory loop
        for(int i = 0; i < evo_arg.gauge_unload_period; ++i, ++traj) {

            measure_plaq(common_arg_plaq);
            // measure_wline(common_arg_wline);
            // measure_tc(common_arg_top, 20);

	    printf("prepbp Fbfm::bfm_args[0].solver = %d\n",Fbfm::bfm_args[0].solver);
	    printf("prepbp cur idx %d, Fbfm::bfm_args[cur_idx].solver = %d\n",Fbfm::current_arg_idx, Fbfm::bfm_args[Fbfm::current_arg_idx].solver);

            measure_pbp(common_arg_pbp, traj);

            //VRB.ElapsedTime("CPS", "main[2]");
            VRB.Result( "CPS", "main", "omp_get_num_threads[4] -> %d", omp_get_num_threads() );
            run_hmc(common_arg_hmc, traj, ab3);
        }//End of inter-cfg sweep
	
	if(do_config_repro){
	  config_repro(repro_config, repro_tolerance);
	  do_config_repro = false;
	}	 
        checkpoint(traj);

    } //End config loop

    AlgIntAB::Destroy(ab3);
    AlgIntAB::Destroy(ab2);
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


void config_repro(const char* repro_config, const Float &tolerance){
  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  LatticeContainer cfg_store;
  cfg_store.Get(lat);
  
  {
    ReadLatticeParallel readLat;
    if(!UniqueID()) printf("Reading: configuration %s for reproduction test\n",repro_config);
    readLat.read(lat,repro_config);
  }
  
  Float* compare_to = (Float*)lat.GaugeField();
  Float* generated = (Float*)cfg_store.GaugeField();
  
  Float fail(0);

  long gauge_size = 18*4*GJP.VolNodeSites();

  for(int i=0;i<gauge_size;i++){
    int rem = i;
    int midx = rem % 18; rem/=18;
    int mu = rem % 4; rem/=4;
    int x[4];
    for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

    if(fabs(generated[i] - compare_to[i]) > tolerance){
      printf("Reproduce comparison fail: midx %d mu %d (%d %d %d %d). Generated: %.10f Loaded: %.10f\n",midx,mu,x[0],x[1],x[2],x[3],generated[i],compare_to[i]);
      fail=1;
    }
  }
  glb_sum(&fail);

  if(fail > 0 && !UniqueID()){ printf("Failed reproduction test\n"); exit(-1); }
  else if(!UniqueID()) printf("Passed reproduction test\n");

  cfg_store.Set(lat);
  LatticeFactory::Destroy();
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
#if TARGET == BGQ
  wl.setSerial();
#endif
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

    //encode_vml(b_arg, traj);
    encode_vml(hmc_arg, traj);
    encode_vml(gauge_arg, traj);
    encode_vml(rat_quo_tm_1_arg, traj);
    encode_vml(rat_quo_tm_2_arg, traj);
    encode_vml(quo_arg, traj);
    encode_vml(rat_quo_arg, traj);
    encode_vml(ab1_arg, traj);
    encode_vml(ab2_arg, traj);
    encode_vml(ab3_arg, traj);
    encode_vml(pbp_arg, traj);
    // encode_vml(ape_arg, traj);
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

void measure_tc(CommonArg &common_arg, int cycle)
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
        VRB.Result(cname,fname,"   running ape\n"); ape.run();
        VRB.Result(cname,fname,"   running ape\n"); ape.run();
        VRB.Result(cname,fname,"   running ape\n"); ape.run();
    }
    tcharge.run();
    // restore the lattice
    lat_cont.Set(lat);
    LatticeFactory::Destroy();

    dtime += dclock();
    print_flops("AlgTcharge", "run()", 0, dtime);
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

        Lattice &lat = LatticeFactory::Create(F_CLASS_BFM, G_CLASS_NONE);
        VRB.Result( "cps", "measure_pbp", "LatticeFactory::Create(F_CLASS_BFM, G_CLASS_NONE)" );

	printf("prepbp post lat create cur idx %d, Fbfm::bfm_args[cur_idx].solver = %d\n",Fbfm::current_arg_idx, Fbfm::bfm_args[Fbfm::current_arg_idx].solver);
	fflush(stdout);

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
