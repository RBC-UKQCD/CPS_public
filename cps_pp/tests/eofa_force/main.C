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
#include <alg/alg_remez.h>
#include <alg/common_arg.h>
#include <alg/hmc_arg.h>
#include <alg/hmd_arg.h>

#include <alg/no_arg.h>
#include <alg/do_arg.h>
#include <alg/alg_int.h>
#include <alg/int_arg.h>
#include <alg/pbp_arg.h>
#include <alg/array_arg.h>

#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/qcdio.h>
#include <util/ReadLatticePar.h>
#include <util/qioarg.h>
#include <util/command_line.h>

#include <omp.h>
#include <pthread.h>

#ifdef HAVE_BFM
#include <chroma.h>
#endif

#define HAVE_DOEXT

USING_NAMESPACE_CPS

static const char* cname = "";

// global argument structures
DoArg do_arg;
#ifdef HAVE_DOEXT
DoArgExt doext_arg;
#endif
EvoArg evo_arg;

void setup(int argc, char **argv);
void load_checkpoint(int traj);

#ifdef USE_BFM
void init_bfm(int* argc, char ***argv, int nthrds);
#endif

inline int Chdir(const char* dir)
{
	const char *fname = "Chdir(char*)";
	if(chdir(dir)!=0){
		ERR.General(cname, fname, "Cannot change directory to %s\n", dir);
	}
	return 0;
}

// link_perturbation_magnitude should be a small number, say 1e-5
// The smaller the perturbation the better the agreement you can
// expect between the force and the numerical derivative of the action.
void test_force(AlgActionEOFA& integrator, Float link_perturbation_size)
{
	const char* fname = "test_force";

	VRB.Result(cname, fname, "Comparing force to numerical derivative of action.\n");
	VRB.Result(cname, fname, "link_perturbation_size = %e\n", link_perturbation_size);

	// integrator.init();
	integrator.heatbath();

	Float ham_initial = integrator.energy();
	glb_sum(&ham_initial);

	// compute force on each link
	int num_links = 4*GJP.VolNodeSites();
	Matrix *force = (Matrix*) smalloc(num_links*sizeof(Matrix), "force", fname, cname);
	for(int i=0; i<num_links; i++){ force[i].ZeroMatrix(); }
	integrator.prepare_fg(force, 1.0);

	// Make a small random perturbation to every link.
	// For every link we generate a random 3x3 traceless antihermitian
	// matrix with small entries. The magnitude of the entries is
	// determined by link_perturbation_magnitude. Then each link
	// gets multiplied by the exponential of the perturbation matrix
	Matrix *link_perturbation = (Matrix*) smalloc(num_links*sizeof(Matrix), "link_perturbation", fname, cname);
	for(int i=0; i<num_links; i++){
		Matrix *lp = &link_perturbation[i];
		for(int j=0; j<18; j++){
			((Float*)lp)[j] = link_perturbation_size * 2 * (-0.5 + drand48());
		}
		lp->TrLessAntiHermMatrix();
	}
  
	{
		GnoneFnone lat;
		// U_{x,u} = exp(link_perturbation) * U_{x,u}
		lat.EvolveGfield(link_perturbation, 1.0);
	}

	// The force is supposed to be the derivative of the energy with respect
	// to the links, so we can use it to predict the change in the action
	// under the given perturbation of the links.
	// predicted change = -Tr [link_perturbation * T^a d^a_x,u f]
	Float predicted_change = 0.0;
	for(int i=0; i<num_links; i++){
		Matrix tmp;
		tmp.DotMEqual(link_perturbation[i], force[i]);
		predicted_change += tmp.Tr().real();
	}
	glb_sum(&predicted_change);

	// Need to do a fake evolve step because some integrators cache
	// the initial energy. To make them actually recompute the energy
	// we need to insert a fake evolve step.
	VRB.Result(cname, fname, "Doing dummy evolve by zero time units...\n");
  integrator.set_skip_force(true);
	integrator.evolve(0.0, 1);
	VRB.Result(cname, fname, "Finished dummy evolve.\n");
  
  // compute the actual hamiltonian of the perturbed gauge field
	Float ham_final = integrator.energy();
	glb_sum(&ham_final);
  
	Float actual_change = ham_final - ham_initial;

	// The two different ways of computing the change in the energy should agree
	// The number of digits of agreement should increase if you make link_perturbation_size smaller.
	VRB.Result(cname, fname, "Initial hamiltonian                        = %0.16e\n", ham_initial);
	VRB.Result(cname, fname, "Final hamiltonian                          = %0.16e\n", ham_final);
	VRB.Result(cname, fname, "Change in hamiltonian                      = %0.16e\n", actual_change);
	VRB.Result(cname, fname, "Change in hamiltonian predicted from force = %0.16e\n", predicted_change);
	VRB.Result(cname, fname, "Agreement between these numbers should improve as link_perturbation_size goes to zero.\n");

	sfree(force, "force", fname, cname);
	sfree(link_perturbation, "link_perturbation", fname, cname);
}

void construct_integrator_and_test_force(Float link_perturbation_size) 
{
	const char* fname = "construct_integrator";

	ActionEOFAArg eofa_arg;
	if(!eofa_arg.Decode(CommandLine::arg(), "eofa_arg")){ 
		ERR.General(cname, fname, "Bum eofa_arg.\n"); 
	}
  
	AlgMomentum mom;
	AlgActionEOFA s_quark(mom, eofa_arg);

	test_force(s_quark, link_perturbation_size);
}

int main(int argc, char **argv)
{
	const char* fname = "main()";

	setup(argc, argv);

	int traj = evo_arg.traj_start;
	int g_int = evo_arg.gauge_unload_period;

	char lat_file[256];
	char gauge_file[256];
	char plaq_file[256];

	sprintf(gauge_file, "gauge.H");
	for(int conf=0; conf<evo_arg.gauge_configurations; conf++){
		sprintf(lat_file, "%s.%d", evo_arg.gauge_file_stem, traj);
		sprintf(plaq_file, "plaq.%d", traj);
	}

	{
		GwilsonFdwf lat;
		if(do_arg.start_conf_kind == START_CONF_FILE){ ReadLatticeParallel(lat, lat_file); }
		NoArg no_arg;
		CommonArg common_arg;
		common_arg.set_filename(plaq_file);
		AlgPlaq plaq(lat, &common_arg, &no_arg);
		plaq.run();
	}

	Float link_perturbation_size = CommandLine::arg_as_Float();
  VRB.Result(cname, fname, "Got: link_perturbation_size = %g\n", link_perturbation_size);
	construct_integrator_and_test_force(link_perturbation_size);

	VRB.Result(cname, fname, "Program exiting normally\n");

	End();

	return 0;
}

void setup(int argc, char **argv)
{
	const char* fname = "setup()";

	Start(&argc, &argv);
	CommandLine::is(argc,argv);

	if(!do_arg.Decode(CommandLine::arg(), "do_arg")){ ERR.General(cname, fname, "Bum do_arg\n"); }
  #ifdef HAVE_DOEXT
  if(!doext_arg.Decode(CommandLine::arg(), "doext_arg")){ ERR.General(cname, fname, "Bum doext_arg\n"); }
  #endif
	if(!evo_arg.Decode(CommandLine::arg(),"evo_arg")){ ERR.General(cname, fname, "Bum evo_arg\n"); }
  #ifdef USE_QUDA
  if(!QudaParam.Decode(CommandLine::arg(), "QudaParam")){ ERR.General(cname, fname, "Bum quda_arg\n"); }
  #endif

	Chdir(evo_arg.work_directory);

	GJP.Initialize(do_arg);
  #ifdef HAVE_DOEXT
  GJP.InitializeExt(doext_arg);
  #endif
  GJP.ZMobius_PC_Type(ZMOB_PC_ORIG);
  LRG.setSerial();
	LRG.Initialize();

  int nthreads;
  const char* nthr_str = getenv("OMP_NUM_THREADS");
  if(nthr_str) {
    sscanf(nthr_str, "%d", &nthreads);
	  GJP.SetNthreads(nthreads);
  } else {
    nthreads = 1;
  }
  VRB.Result(cname, fname, "GJP.Nthreads() = %d\n", GJP.Nthreads());
  
  #ifdef USE_BFM
  #ifdef USE_NEW_BFM_GPARITY
	Fbfm::bfm_args[0].threads = nthreads;
  #endif
	init_bfm(&argc, &argv, nthreads);
  #endif
}

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
	
	// mixed-precision CG *based on environmental variable*, *true by default*
	char* use_mixed_solver_env = getenv("use_mixed_solver");
	Fbfm::use_mixed_solver = true;
	if(use_mixed_solver_env && (strcmp(use_mixed_solver_env,"false")==0)){
		Fbfm::use_mixed_solver = false;
	}
	VRB.Result("CPS", "init_bfm", "Fbfm::use_mixed_solver: %d\n", Fbfm::use_mixed_solver);

	// Set for EOFA
	bfm_arg.solver = HtCayleyEOFA; // DWF
	// bfm_arg.solver = HmCayleyEOFA; // Mobius
	// Fbfm::use_eofa = true;
	Fbfm::use_eofa_4d_precond = true;
	Fbfm::use_chebyshev_precond_Minv = false;

	bfm_arg.precon_5d = 0;
	bfm_arg.Ls = GJP.SnodeSites();
	bfm_arg.M5 = GJP.DwfHeight();
	bfm_arg.mass = 0.1;
	bfm_arg.residual = 1.0e-08;
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

	// mobius scale = b + c in Andrew's notation
#ifdef USE_NEW_BFM_GPARITY
	bfm_arg.mobius_scale = 1.0;
#else
	bfmarg::mobius_scale = 1.0;
#endif

	bfmarg::Threads(nthrds);
	bfmarg::Reproduce(0);
	bfmarg::ReproduceChecksum(0);
	bfmarg::ReproduceMasterCheck(0);
	bfmarg::Verbose(0);
}
#endif
