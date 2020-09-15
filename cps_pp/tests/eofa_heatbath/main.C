#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
#include <comms/scu.h>
#include <comms/glb.h>

#include <util/lattice.h>
#include <util/lattice/fbfm.h>
#include <util/random.h>
#include <util/time_cps.h>

#include <alg/alg_hmc.h>
#include <alg/alg_plaq.h>
#include <alg/alg_rect.h>
#include <alg/common_arg.h>
#include <alg/hmc_arg.h>
#include <alg/hmd_arg.h>

#include <alg/alg_int.h>
#include <alg/int_arg.h>

#include <alg/no_arg.h>
#include <alg/do_arg.h>
#include <alg/alg_plaq.h>
#include <alg/alg_pbp.h>
#include <alg/pbp_arg.h>
#include <alg/alg_remez.h>
#include <alg/array_arg.h>

#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/qcdio.h>
#include <util/WriteLatticePar.h>
#include <util/ReadLatticePar.h>
#include <util/qioarg.h>
#include <util/command_line.h>

#include <omp.h>
#include <pthread.h>

#ifdef HAVE_BFM
#include <chroma.h>
#endif
//--------------------------------------------------------------

#define HAVE_DOEXT

USING_NAMESPACE_CPS
using namespace std;

const char* cname = "";

DoArg do_arg;
#ifdef HAVE_DOEXT
DoArgExt doext_arg;
#endif
EvoArg evo_arg;
ActionEOFAArg eofa_arg;

#ifdef USE_BFM
void init_bfm(int *argc, char **argv[], int nthrds)
{
	cps_qdp_init(argc, argv);
	Chroma::initialize(argc, argv);
	multi1d<int> nrow(Nd);

	for(int i=0; i<Nd; ++i){ nrow[i] = GJP.Sites(i); }

	Layout::setLattSize(nrow);
	Layout::create();

	Fbfm::current_arg_idx = 0;
	bfmarg &bfm_arg = Fbfm::bfm_args[0];

	// mixed-precision CG *based on environment variable*, *true by default*
	char* use_mixed_solver_env = getenv("use_mixed_solver");
	Fbfm::use_mixed_solver = true;
	if(use_mixed_solver_env && (strcmp(use_mixed_solver_env,"false") == 0)){
		Fbfm::use_mixed_solver = false;
	}
	VRB.Result("CPS", "init_bfm", "Fbfm::use_mixed_solver: %d\n", Fbfm::use_mixed_solver);

	// bfm_arg.solver = HtCayleyEOFA; // DWF
	bfm_arg.solver = HmCayleyEOFA; // Mobius
	// Fbfm::use_eofa = true;
	Fbfm::use_eofa_4d_precond = true;
	Fbfm::use_chebyshev_precond_Minv = false;

	bfm_arg.precon_5d = 0;
	bfm_arg.Ls = GJP.SnodeSites();
	bfm_arg.M5 = GJP.DwfHeight();
	bfm_arg.mass = 0.1;
	bfm_arg.residual = 1e-8;
	bfm_arg.max_iter = 10000;
	bfm_arg.Csw = 0.0;

	bfm_arg.node_latt[0] = QDP::Layout::subgridLattSize()[0];
	bfm_arg.node_latt[1] = QDP::Layout::subgridLattSize()[1];
	bfm_arg.node_latt[2] = QDP::Layout::subgridLattSize()[2];
	bfm_arg.node_latt[3] = QDP::Layout::subgridLattSize()[3];

	multi1d<int> procs = QDP::Layout::logicalSize();
	multi1d<int> ncoor = QDP::Layout::nodeCoord();

	for(int i=0; i<4; ++i){ bfm_arg.ncoor[i] = ncoor[i]; }

	bfm_arg.local_comm[0] = procs[0] > 1 ? 0 : 1;
	bfm_arg.local_comm[1] = procs[1] > 1 ? 0 : 1;
	bfm_arg.local_comm[2] = procs[2] > 1 ? 0 : 1;
	bfm_arg.local_comm[3] = procs[3] > 1 ? 0 : 1;

	if(GJP.Gparity()){
		bfm_arg.gparity = 1;
		if(!UniqueID()){ printf("G-parity directions: "); }
		for(int d=0; d<3; d++){
			if(GJP.Bc(d) == BND_CND_GPARITY){
				bfm_arg.gparity_dir[d] = 1;
				if(!UniqueID()){ printf("%d ", d); }
			}
		}
		for(int d=0; d<4; d++){ bfm_arg.nodes[d] = procs[d]; }
		if(!UniqueID()){ printf("\n"); }
	}

	// mobius_scale = b + c in Andrew's notation
#ifdef USE_NEW_BFM_GPARITY
	bfm_arg.mobius_scale = 32.0/12.0;
#else
	bfmarg::mobius_scale = 32.0/12.0;
#endif

	bfmarg::Threads(nthrds);
	bfmarg::Reproduce(0);
	bfmarg::ReproduceChecksum(0);
	bfmarg::ReproduceMasterCheck(0);
	bfmarg::Verbose(0);
}
#endif

inline int Chdir(const char* dir){
	if(chdir(dir)!=0){
		printf("cannot change directory to %s\n", dir);
		exit(-1);
	}
	return 0;
}

int main(int argc, char *argv[])
{ 
	char *fname = "main()";
	Float dtime;
  
	CommonArg common_arg_plaq;
	CommonArg common_arg_pbp;
	CommonArg common_arg_hmc;

	Start(&argc, &argv);
	CommandLine::is(argc, argv);

	if(!do_arg.Decode(CommandLine::arg(), "do_arg")){ 
		do_arg.Encode("bum_arg", "bum_arg");
		printf("Bum do_arg\n"); 
		exit(-1);
	}
  #ifdef HAVE_DOEXT
  if(!doext_arg.Decode(CommandLine::arg(), "doext_arg")){ printf("Bum doext_arg\n"); }
  #endif
	if(!evo_arg.Decode(CommandLine::arg(), "evo_arg")){ printf("Bum evo_arg\n"); exit(-1); }
	if(!eofa_arg.Decode(CommandLine::arg(), "eofa_arg")){ printf("Bum eofa_arg\n"); exit(-1); }

  #ifdef USE_QUDA
  if(!QudaParam.Decode(CommandLine::arg(), "QudaParam")){ printf("Bum quda_arg\n"); }
  #endif

	Chdir(evo_arg.work_directory);

	GJP.Initialize(do_arg);
  #ifdef HAVE_DOEXT
  GJP.InitializeExt(doext_arg);
  #endif
  GJP.ZMobius_PC_Type(ZMOB_PC_ORIG);
  LRG.setSerial();
	LRG.Initialize();

  #ifdef USE_BFM
	init_bfm(&argc, &argv, nthreads);
  #else
  {
    int threads;
    const char* nthr_str = getenv("OMP_NUM_THREADS");
    if(nthr_str) {
      sscanf(nthr_str, "%d", &threads);
      GJP.SetNthreads(threads);
    }
    if(!UniqueID()){ printf("GJP.Nthreads() = %d\n", GJP.Nthreads()); }
  }
  #endif

	// Outer config loop
	int traj = evo_arg.traj_start;
	int g_int = evo_arg.gauge_unload_period;
  
	//!< Create fictitous Hamiltonian (mom + action)
	AlgMomentum mom;
	int n_step = eofa_arg.num_mass.num_mass_len;

	char lat_file[256];
	char gauge_file[256];
	char plaq_file[256];
	char vml_file[256];
	Float *rw_fac, *norm;
	rw_fac = new Float[n_step];
	norm = new Float[n_step];
	Float rw;
  
	sprintf(gauge_file, "gauge.H");
	bool if_bound = 0;
	if(eofa_arg.eigen.eigen_measure == EIGEN_MEASURE_YES){ if_bound = 1; }
  
	for(int conf=0; conf<evo_arg.gauge_configurations; conf++){

		sprintf(lat_file, "%s.%d", evo_arg.gauge_file_stem, traj);
		sprintf(plaq_file, "plaq.%d", traj);
		sprintf(vml_file, "eofa_arg.%d", traj);
    if(!eofa_arg.Encode(vml_file, "eofa_arg")){
			printf("bad eofa_arg encode\n");
			exit(-1);
		}

		{
			GwilsonFmobius lat;
			if(do_arg.start_conf_kind == START_CONF_FILE){ ReadLatticeParallel(lat, lat_file); }
			NoArg no_arg;
			CommonArg common_arg;
			common_arg.set_filename(plaq_file);
			AlgPlaq plaq(lat, &common_arg, &no_arg);
			plaq.run();
		}

    AlgActionEOFA s_quark(mom, eofa_arg);
    s_quark.heatbath();
    s_quark.set_skip_force(true);
    s_quark.evolve(0.0,1); // fake evolve step since integrator caches energy
    s_quark.energy(); // recompute energy after heatbath

		traj += g_int;
	} //End config loop

	End();
	
	return(0);
}
