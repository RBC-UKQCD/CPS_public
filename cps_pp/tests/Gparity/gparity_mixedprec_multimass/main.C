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
#include <util/dirac_op.h>
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

#include <alg/lanc_arg.h>
#include <alg/prop_attribute_arg.h>
#include <alg/gparity_contract_arg.h>
#include <alg/propmanager.h>
#include <alg/alg_gparitycontract.h>

#include <util/gparity_singletodouble.h>

//some piece of **** defines these elsewhere, so the bfm header gets screwed up
#undef ND
#undef SPINOR_SIZE
#undef HALF_SPINOR_SIZE
#undef GAUGE_SIZE
#undef Nmu
#undef Ncb
#undef NMinusPlus
#undef Minus
#undef Plus
#undef DaggerYes
#undef DaggerNo
#undef SingleToDouble
#undef DoubleToSingle
#undef Odd
#undef Even


#include <util/lattice/bfm_evo.h>
#include <util/lattice/bfm_eigcg.h>
#include <util/lattice/fbfm.h>
#include <util/wilson.h>
#include <util/verbose.h>
#include <util/gjp.h>
#include <util/error.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <util/enum_func.h>
#include <util/sproj_tr.h>
#include <util/time_cps.h>
#include <util/lattice/fforce_wilson_type.h>

#include<sstream>

#include<omp.h>

#ifdef HAVE_BFM
#include <chroma.h>
#endif

using namespace std;
USING_NAMESPACE_CPS

void setup_double_latt(Lattice &double_latt, Matrix* orig_gfield, bool gparity_X, bool gparity_Y){
  //orig latt ( U_0 U_1 ) ( U_2 U_3 ) ( U_4 U_5 ) ( U_6 U_7 )
  //double tatt ( U_0 U_1 U_2 U_3 ) ( U_4 U_5 U_6 U_7 ) ( U_0* U_1* U_2* U_3* ) ( U_4* U_5* U_6* U_7* )

  Matrix *dbl_gfield = double_latt.GaugeField();

  if(!UniqueID()){ printf("Setting up 1f lattice.\n"); fflush(stdout); }
  SingleToDoubleLattice lattdoubler(gparity_X,gparity_Y,orig_gfield,double_latt);
  lattdoubler.Run();
  if(!UniqueID()){ printf("Finished setting up 1f lattice\n"); fflush(stdout); }
}
void setup_double_rng(bool gparity_X, bool gparity_Y){
  //orig 4D rng 2 stacked 4D volumes
  //orig ([R_0 R_1][R'_0 R'_1])([R_2 R_3][R'_2 R'_3])([R_4 R_5][R'_4 R'_5])([R_6 R_7][R'_6 R'_7])
  //double (R_0 R_1 R_2 R_3)(R_4 R_5 R_6 R_7)(R'_0 R'_1 R'_2 R'_3)(R'_4 R'_5 R'_6 R'_7)
  
  //orig 5D rng 2 stacked 4D volumes per ls/2 slice (ls/2 as only one RNG per 2^4 block)

  SingleToDouble4dRNG fourDsetup(gparity_X,gparity_Y);
  SingleToDouble5dRNG fiveDsetup(gparity_X,gparity_Y);
  
  LRG.Reinitialize(); //reset the LRG and prepare for doubled lattice form
  
  if(!UniqueID()){ printf("Setting up 1f 4D RNG\n"); fflush(stdout); }
  fourDsetup.Run();      
  if(!UniqueID()){ printf("Setting up 1f 5D RNG\n"); fflush(stdout); }
  fiveDsetup.Run();    
}
void setup_double_matrixfield(Matrix* double_mat, Matrix* orig_mat, int nmat_per_site, bool gparity_X, bool gparity_Y){
  if(!UniqueID()){ printf("Setting up 1f matrix field.\n"); fflush(stdout); }
  SingleToDoubleMatrixField doubler(gparity_X,gparity_Y,nmat_per_site,orig_mat,double_mat);
  doubler.Run();
  if(!UniqueID()){ printf("Finished setting up 1f matrixfield\n"); fflush(stdout); }
}
void setup_double_5d_vector(Vector *double_vect, Vector* orig_vect, bool gparity_X, bool gparity_Y){
  if(!UniqueID()){ printf("Setting up 1f vector field.\n"); fflush(stdout); }
  SingleToDouble5dVectorField doubler(gparity_X, gparity_Y, orig_vect, double_vect, CANONICAL);
  doubler.Run();
  if(!UniqueID()){ printf("Finished setting up 1f vector field\n"); fflush(stdout); }
}
  
void GaugeTransformU(Matrix *gtrans, Lattice &lat);

void convert_ferm_cpsord_sord(Float *cps, Float* &sord, bfm_evo<Float> &bfm){
  Fermion_t handle[2] = { bfm.allocFermion(), bfm.allocFermion() };
  bfm.cps_impexFermion(cps,handle,1);
  
  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();
  if(GJP.Gparity()) f_size*=2;
  sord = (Float *)pmalloc(sizeof(Float) * f_size);
  bfm.cps_impexFermion_s(sord,handle,0);

  bfm.freeFermion(handle[0]);
  bfm.freeFermion(handle[1]);
}
void convert_ferm_sord_cpsord(Float *sord, Float* &cps, bfm_evo<Float> &bfm){
  Fermion_t handle[2] = { bfm.allocFermion(), bfm.allocFermion() };
  bfm.cps_impexFermion_s(sord,handle,1);
  
  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();
  if(GJP.Gparity()) f_size*=2;
  cps = (Float *)pmalloc(sizeof(Float) * f_size);
  bfm.cps_impexFermion(cps,handle,0);

  bfm.freeFermion(handle[0]);
  bfm.freeFermion(handle[1]);
}


void setup_bfmargs(bfmarg &dwfa, const BfmSolver &solver){
  printf("Setting up bfmargs\n");

   int nthreads = 1; 
#if TARGET == BGQ
   nthreads = 64;
#endif
   omp_set_num_threads(nthreads);

  dwfa.node_latt[0]  = GJP.XnodeSites();
  dwfa.node_latt[1]  = GJP.YnodeSites();
  dwfa.node_latt[2]  = GJP.ZnodeSites();
  dwfa.node_latt[3]  = GJP.TnodeSites();
  
  multi1d<int> ncoor(4);
  multi1d<int> procs(4);
  for(int i=0;i<4;i++){ ncoor[i] = GJP.NodeCoor(i); procs[i] = GJP.Nodes(i); }

  if(GJP.Gparity()){
    dwfa.gparity = 1;
    printf("G-parity directions: ");
    for(int d=0;d<3;d++)
      if(GJP.Bc(d) == BND_CND_GPARITY){ dwfa.gparity_dir[d] = 1; printf("%d ",d); }
      else dwfa.gparity_dir[d] = 0;
    for(int d=0;d<4;d++){
      dwfa.nodes[d] = procs[d];
      dwfa.ncoor[d] = ncoor[d];
    }
    printf("\n");
  }

  dwfa.verbose=1;
  dwfa.reproduce=0;
  bfmarg::Threads(nthreads);
  bfmarg::Reproduce(0);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(1);

  for(int mu=0;mu<4;mu++){
    if ( procs[mu]>1 ) {
      dwfa.local_comm[mu] = 0;
      printf("Non-local comms in direction %d\n",mu);
    } else { 
      dwfa.local_comm[mu] = 1;
      printf("Local comms in direction %d\n",mu);
    }
  }

  dwfa.precon_5d = 1;
  if(solver == HmCayleyTanh){
    dwfa.precon_5d = 0; //mobius uses 4d preconditioning
    dwfa.mobius_scale = 2.0; //b = 0.5(scale+1) c=0.5(scale-1), hence this corresponds to b=1.5 and c=0.5, the params used for the 48^3
  }
  dwfa.Ls   = GJP.SnodeSites();
  dwfa.solver = solver;
  dwfa.M5   = toDouble(GJP.DwfHeight());
  dwfa.mass = toDouble(0.001);
  dwfa.Csw  = 0.0;
  dwfa.max_iter = 5000;
  dwfa.residual = 1e-08;
  printf("Finished setting up bfmargs\n");
}

Float* rand_5d_canonical_fermion(Lattice &lat){
  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();
  if(GJP.Gparity()) f_size*=2;
  Float *v1 = (Float *)pmalloc(sizeof(Float) * f_size);
  printf("Making random gaussian 5d vector\n");
  lat.RandGaussVector((Vector*)v1, 0.5, 2, CANONICAL, FIVE_D);
  printf("Finished making random gaussian vector\n");
  return v1;
}

#include <bfm_vmx.h>

static int no_gparity_test(Lattice* lattice, const BfmSolver &solver, const double &min_fp_resid){
  //1) Test the single shift inverter by comparing to the multi-shift inverter with a single shift
  
  //template<typename Float>
  //int threaded_CGNE_MdagM_plus_shift(Fermion_t psi, Fermion_t src, Float shift, bfm_evo<Float> &bfm)
  
  bfmarg dwfa;
  setup_bfmargs(dwfa,solver);
  
  bfm_evo<double> bfm_d;
  bfm_d.init(dwfa);
  bfm_d.verbose = 1;
    
  lattice->BondCond();
  Float* gauge = (Float*) lattice->GaugeField();
  bfm_d.cps_importGauge(gauge);  


  LatRanGen LRGbak(LRG);
  Float* v1 = rand_5d_canonical_fermion(*lattice);
  printf("Restoring RNG\n"); fflush(stdout);
  LRG = LRGbak;

  printf("Allocating fermions\n"); fflush(stdout);
  Fermion_t src[2] = {bfm_d.allocFermion(), bfm_d.allocFermion()}; //odd/even

  printf("Impexing random vector to bfm src vectors\n"); fflush(stdout);
  bfm_d.cps_impexFermion(v1,src,1);
  
  Fermion_t src_copy = bfm_d.allocFermion();
  bfm_d.copy(src_copy,src[0]);

  printf("Allocating sol_1\n");
  Fermion_t sol_1 = bfm_d.allocFermion();
  printf("Setting it to zero\n");
#pragma omp parallel
  {
    bfm_d.set_zero(sol_1); //start from same zero source as the multi-shift
  }   
 
  double shift = 0.5;
  bfm_d.time_report_iter = 1;

  printf("Starting mixed_cg::threaded_CGNE_MdagM_plus_shift\n"); fflush(stdout);
#pragma omp parallel
  {
    mixed_cg::threaded_CGNE_MdagM_plus_shift<double>(sol_1, src[0], shift, bfm_d);
  }  
  printf("Finished mixed_cg::threaded_CGNE_MdagM_plus_shift\n"); fflush(stdout);

  Fermion_t sol_2[1] = { bfm_d.allocFermion() };
  
  double alpha = 1.0;
  double residual = bfm_d.residual;

#pragma omp parallel
  {
    bfm_d.CGNE_prec_MdagM_multi_shift(sol_2,src[0],&shift,&alpha,1,&residual,0);
  }  

  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();
  long f_size_cb = f_size/2;
  Float* sol_1_cps = (Float *)pmalloc(sizeof(Float) * f_size_cb);
  Float* sol_2_cps = (Float *)pmalloc(sizeof(Float) * f_size_cb);
  bfm_d.cps_impexcbFermion(sol_1_cps,sol_1,0,1);
  bfm_d.cps_impexcbFermion(sol_2_cps,sol_2[0],0,1);

  bfm_d.freeFermion(sol_1);
  bfm_d.freeFermion(sol_2[0]);
  //delete[] sol_2;

  bool fail(false);
  
  for(int i=0;i<f_size_cb;i++){
    // int rem = i;
    // int midx = rem % 24; rem/=24;

    // int x[5];
    // for(int j=0;j<5;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

    if( fabs(sol_1_cps[i] - sol_2_cps[i]) > 1e-07 ){ 
      // printf("MdagM+shift invert test fail (%d %d %d %d %d, %d): %f %f\n",x[0],x[1],x[2],x[3],x[4],
      // 	     midx, sol_1_cps[i],sol_2_cps[i]); 
      printf("MdagM+shift invert test fail %d: %f %f\n",i, sol_1_cps[i],sol_2_cps[i]); 
      fail = true; 
    }//else printf("MdagM+shift invert test pass %d: %f %f\n",i, sol_1_cps[i],sol_2_cps[i]); 
  }
  if(fail){
    printf("Failed MdagM+shift invert test\n"); exit(-1);
  }else printf("Passed MdagM+shift invert test\n");

  pfree(sol_1_cps);
  pfree(sol_2_cps);


  //2) Test and compare the mixed precision multi-shift inverter

  bfm_d.comm_end();

  bfm_evo<float> bfm_f;
  bfm_f.init(dwfa);
  bfm_f.verbose = 1;

  bfm_f.cps_importGauge(gauge);  
  bfm_f.comm_end();
  bfm_d.comm_init();
    
  Fermion_t sol_mm_mixed[4] = { bfm_d.allocFermion(), bfm_d.allocFermion(), bfm_d.allocFermion(), bfm_d.allocFermion() };
  Fermion_t sol_mm_std[4] = { bfm_d.allocFermion(), bfm_d.allocFermion(), bfm_d.allocFermion(), bfm_d.allocFermion() };
  
  double mass[] = {0.002,0.004,0.006,0.008};
  double alpha_ms[] = {1.0,1.0,1.0,1.0};
  double resid[] = {1e-08, 1e-08, 1e-08, 1e-08};
  
  int max_cycle = 1000;
  
  struct timeval start_mixed,stop_mixed,start_std,stop_std;
  gettimeofday(&start_mixed,NULL);

  #pragma omp parallel
  {
    if(!UniqueID()) printf("Doing mixed precision multi-shift\n");
    mixed_cg::threaded_cg_mixed_MdagM_multi_shift(sol_mm_mixed,src[0],mass,alpha_ms,bfm_d,bfm_f,4,resid,min_fp_resid,0,max_cycle);
  }
  gettimeofday(&stop_mixed,NULL);
  
  gettimeofday(&start_std,NULL);
  #pragma omp parallel
  {
    if(!UniqueID()) printf("Doing double precision multi-shift\n");
    bfm_d.CGNE_prec_MdagM_multi_shift(sol_mm_std,src[0],mass,alpha_ms,4,resid,0);
  }
  gettimeofday(&stop_std,NULL);

  for(int i=0;i<shift;i++){
    struct timeval diff_mixed;
    timersub(&stop_mixed,&start_mixed,&diff_mixed);
    struct timeval diff_std;
    timersub(&stop_std,&start_std,&diff_std);
     
    if(!UniqueID()) printf("Mixed-precision time %d.%6.6d s, double precision time %d.%6.6d s\n",diff_mixed.tv_sec,diff_mixed.tv_usec,diff_std.tv_sec,diff_std.tv_usec);
  }  

  Float* sol_mm_mixed_cps = (Float *)pmalloc(sizeof(Float) * f_size_cb);
  Float* sol_mm_std_cps = (Float *)pmalloc(sizeof(Float) * f_size_cb);
  for(int shift=0;shift<4;shift++){
    bfm_d.cps_impexcbFermion(sol_mm_mixed_cps,sol_mm_mixed[shift],0,1);
    bfm_d.cps_impexcbFermion(sol_mm_std_cps,sol_mm_std[shift],0,1);

    for(int i=0;i<f_size_cb;i++){
      if( fabs(sol_mm_mixed_cps[i] - sol_mm_std_cps[i]) > 1e-07 ){ 
	printf("Mixed-prec multi-mass test fail. Shift %d, offset %d: %.10e %.10e, diff %.10e\n",shift, i,sol_mm_mixed_cps[i],sol_mm_std_cps[i], sol_mm_mixed_cps[i]-sol_mm_std_cps[i]); 
	fail = true; 
      }
    }
    if(fail){
      printf("Mixed-prec multi-mass test fail. Shift %d\n",shift); exit(-1);
    }else printf("Passed Mixed-prec multi-mass test. Shift %d.\n",shift);
  }


  {
    //Test the version with guesses by first doing a single precision multi-shift and then using the solutions as guesses
    int nshift = 4;
    for(int i=0;i<nshift;i++) resid[i] = 1e-03;

    gettimeofday(&start_std,NULL);
#pragma omp parallel
    {
      if(!UniqueID()) printf("Doing single precision multi-shift for multi-shift with guesses test\n");
      Fermion_t src_f = bfm_f.threadedAllocFermion();
      mixed_cg::threaded_convFermion(src_f, src[0], bfm_f, bfm_d);
      mixed_cg::switch_comm(bfm_f, bfm_d);
      
      Fermion_t sol_f[nshift];
      for(int i=0;i<nshift;i++) sol_f[i] = bfm_f.threadedAllocFermion();
     
      bfm_f.CGNE_prec_MdagM_multi_shift(sol_f,src_f,mass,alpha_ms,nshift,resid,0);
      
      for(int i=0;i<nshift;i++) mixed_cg::threaded_convFermion(sol_mm_std[i], sol_f[i], bfm_d, bfm_f);
      mixed_cg::switch_comm(bfm_d, bfm_f);
    }

    //Do the double precision version with guesses from single prec solve
    for(int i=0;i<nshift;i++) resid[i] = 1e-08;

#pragma omp parallel
    {
      if(!UniqueID()) printf("Doing double precision multi-shift with guesses\n");
      mixed_cg::CGNE_prec_MdagM_multi_shift_with_guesses(sol_mm_mixed, src[0], sol_mm_std, mass, alpha_ms, nshift, resid, 0, bfm_d);
    }
    gettimeofday(&stop_std,NULL);

    struct timeval diff_std;
    timersub(&stop_std,&start_std,&diff_std);
     
    if(!UniqueID()) printf("Mixed-precision restarted solve time %d.%6.6d s\n",diff_std.tv_sec,diff_std.tv_usec);
  }

  {
    //Test multi-source method. Do a single precision multi-shift solve and use the residuals as sources for double prec solve 
    int nshift = 4;
    for(int i=0;i<nshift;i++) resid[i] = 1e-03;

    gettimeofday(&start_std,NULL);
    Fermion_t multi_src[nshift];

#pragma omp parallel
    {
      if(!UniqueID()) printf("Doing single precision multi-shift for multiple sources test\n");
      Fermion_t src_f = bfm_f.threadedAllocFermion();
      mixed_cg::threaded_convFermion(src_f, src[0], bfm_f, bfm_d);
      mixed_cg::switch_comm(bfm_f, bfm_d);
      
      Fermion_t sol_f[nshift];
      for(int i=0;i<nshift;i++) sol_f[i] = bfm_f.threadedAllocFermion();
     
      bfm_f.CGNE_prec_MdagM_multi_shift(sol_f,src_f,mass,alpha_ms,nshift,resid,0);
      
      for(int i=0;i<nshift;i++) mixed_cg::threaded_convFermion(sol_mm_std[i], sol_f[i], bfm_d, bfm_f);
      mixed_cg::switch_comm(bfm_d, bfm_f);

      Fermion_t tmp = bfm_d.threadedAllocFermion();
      Fermion_t tmp2 = bfm_d.threadedAllocFermion();
      Fermion_t tmp3 = bfm_d.threadedAllocFermion();
      Fermion_t tmp4 = bfm_d.threadedAllocFermion();
      
      for(int i=0;i<nshift;i++){
	multi_src[i] = bfm_d.threadedAllocFermion();
	bfm_d.copy(multi_src[i],src[0]);

	mixed_cg::threaded_convFermion(tmp, sol_f[i], bfm_d, bfm_f);
	mixed_cg::MdagMplusShift(tmp,tmp2,mass[i], tmp3,tmp4, bfm_d);
	bfm_d.axpy(multi_src[i],  tmp2, multi_src[i], -1.0); //residual
      }
      bfm_d.freeFermion(tmp); bfm_d.freeFermion(tmp2); bfm_d.freeFermion(tmp3); bfm_d.freeFermion(tmp4);      
    }

    //Do the double precision version with guesses from single prec solve
    for(int i=0;i<nshift;i++) resid[i] = 1e-08;

#pragma omp parallel
    {
      mixed_cg::CGNE_prec_MdagM_multi_shift_multi_src<double>(sol_mm_mixed, multi_src, mass, alpha_ms, nshift, resid, 0, bfm_d);
    }
    gettimeofday(&stop_std,NULL);

    struct timeval diff_std;
    timersub(&stop_std,&start_std,&diff_std);
     
    if(!UniqueID()) printf("Mixed-precision restarted multi-src solve time %d.%6.6d s\n",diff_std.tv_sec,diff_std.tv_usec);
    for(int i=0;i<nshift;i++) bfm_d.freeFermion(multi_src[i]);
  }

  {
    //Test restarted mixed prec multi-shift
    int nshift = 4;
    double fresid[nshift];
    for(int i=0;i<nshift;i++){
      resid[i] = 1e-08;
      fresid[i] = 1e-04;
    }

    gettimeofday(&start_std,NULL);

#pragma omp parallel
    {
      if(!UniqueID()) printf("Doing restarted multi-mass test\n");
      mixed_cg::threaded_cg_mixed_multi_shift_MdagM(sol_mm_mixed,src[0],mass,alpha_ms,nshift,resid,fresid,0,bfm_d,bfm_f,100);
    }
  }





  //Free mem
  for(int shift=0;shift<4;shift++){
    bfm_d.freeFermion(sol_mm_mixed[shift]);
    bfm_d.freeFermion(sol_mm_std[shift]);
  }
  pfree(sol_mm_mixed_cps);
  pfree(sol_mm_std_cps);

  lattice->BondCond();


  // //Try using Fbfm and old DWF code to get the same MD force vectors

  // if(solver != DWF) return 0;

  // bool fail;

  // bfmarg dwfa;
  // setup_bfmargs(dwfa,solver);
  // int mom_size = 18*4*GJP.VolNodeSites();
  // long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();

  // delete lattice; //temporarily delete the lattice object to remove the scope lock. The memory allocated to the gauge field remains
  // Fbfm::bfm_arg = dwfa;
  // GnoneFbfm *fbfm = new GnoneFbfm;
  

  // //Generate 2 random CANONICAL ordered fermions
  // LatRanGen LRGbak(LRG);
  // Float *v1 = rand_5d_canonical_fermion(*fbfm);
  // Float *v2 = rand_5d_canonical_fermion(*fbfm);
  // LRG = LRGbak;

  // //Convert to bfm format
  // Fermion_t in[2] = {fbfm->bd.allocFermion(), fbfm->bd.allocFermion()};
  // Fermion_t in2[2] = {fbfm->bd.allocFermion(), fbfm->bd.allocFermion()};
  // fbfm->bd.cps_impexFermion(v1,in,1);
  // fbfm->bd.cps_impexFermion(v2,in2,1);


  // //expects phi to be in cps checkerboard format. We can use in[0] and in2[0] as above, and export from bfm to cps format
  // size_t f_size_cb = f_size / 2;
  // Vector* phi1 = (Vector *)pmalloc(sizeof(Float) * f_size_cb);
  // Vector* phi2 = (Vector *)pmalloc(sizeof(Float) * f_size_cb);

  // fbfm->bd.cps_impexcbFermion((Float *)phi1, in[0], 0, 1); //export
  // fbfm->bd.cps_impexcbFermion((Float *)phi2, in2[0], 0, 1);

  // Float *chi = (Float*)phi1;
	
  // Float *v1_fbfm;
  // Float *v2_fbfm;
  // Float *rho = (Float *)pmalloc(sizeof(Float)*f_size_cb);
  // {
  //   Float *v1_fbfm_tmp = (Float *)pmalloc(sizeof(Float) * f_size);
  //   Float *v2_fbfm_tmp = (Float *)pmalloc(sizeof(Float) * f_size);
    
  //   fbfm->MatPc( (Vector*)rho, (Vector*)chi, 0.5, DAG_NO); // rho [ODD] = Mprec chi [ODD]
  //   fbfm->CalcHmdForceVecsBilinear(v1_fbfm_tmp, v2_fbfm_tmp, (Vector*)rho, (Vector*)chi, 0.5); 
    
  //   /*
  //    * From BFM:
  //    * DWF : Mee = Moo = (5-M5)
  //    *       Meo = -1/2 Ddwf_eo (5d hopping term)
  //    *       Moe = -1/2 Ddwf_oe (5d hopping term)
  //    * Mprec = Moo-MoeMee^{-1}Meo
  //    */
  //   //CalcHmdForceVecsBilinear calculates the following:
  //   // v2e      =  Bee * 1/(5-M5) * Meo phi2
  //   // v2o = Boo phi2
  //   // v1e  =  1/(5-M5) Meo^dag phi1
  //   // v1o = 1oo phi1
  //   //For DWF, Boo = Bee = 1 and phi1 = rho,  phi2 = chi in the above naming convention
  //   //hence the net result is
  //   //v1 = ( rho, 1/(5-M5) Meo^dag rho )     v2 = ( chi, 1/(5-M5) * Meo chi )
  //   //i.e.
  //   //v1 = ( rho, -1/[2(5-M5)] Deo^dag rho )     v2 = ( chi, -1/[2(5-M5)] * Deo chi )
  //   //with rho = Mprec chi
  //   //v1 = ( Mprec chi, -1/[2(5-M5)] Deo^dag Mprec chi )     v2 = ( chi, -1/[2(5-M5)] * Deo chi )
  //   //  ((Above vectors are in  (odd,even) format))
    

  //   //convert to CPS canonical ordering
  //   convert_ferm_sord_cpsord(v1_fbfm_tmp, v1_fbfm, fbfm->bd);
  //   convert_ferm_sord_cpsord(v2_fbfm_tmp, v2_fbfm, fbfm->bd);
  //   pfree(v1_fbfm_tmp);
  //   pfree(v2_fbfm_tmp);
  // }


  // Float* mom_test5 = (Float *)pmalloc( sizeof(Float) * mom_size);
  // for(int i=0;i<mom_size;i++) mom_test5[i] = 0.0;
  // fbfm->EvolveMomFforce( (Matrix*)mom_test5, (Vector*)chi, 0.5, 1.0);
  
  // delete fbfm;
  // lattice = new GwilsonFdwf; //put the lattice back where it was
 

  // {
  //   CgArg cg_arg ;
  //   cg_arg.mass = 0.5;
    
  //   Float *v1_dopdwf = (Float *)pmalloc(sizeof(Float)*f_size);
  //   Float *v2_dopdwf = (Float *)pmalloc(sizeof(Float)*f_size);
    
  //   /*  DiracOp(Lattice& latt, Vector *f_field_out, Vector *f_field_in, CgArg *arg,  CnvFrmType convert) */
  //   DiracOpDwf* dwf = new DiracOpDwf(*lattice, (Vector*)v2_dopdwf, (Vector*)v1_dopdwf, &cg_arg, CNV_FRM_YES) ;
    
  //   {
  //     //Test preconditioned matrix
  //     Float *rho_dopdwf = (Float *)pmalloc(sizeof(Float)*f_size_cb);
  //     dwf->MatPc( (Vector*)rho_dopdwf, (Vector*)chi);
      
  //     //Note: CPS Mprec = 1 - 1/[2(5-M5)]^2 Doe Deo
  //     //which differs from BFM Mprec = (5-M5) - 1/[4(5-M5)] Deo Doe
  //     //by a normalization factor Mprec [CPS] = 1/(5-M5) Mprec [BFM]

  //     Float norm = 5-GJP.DwfHeight();
  //     fail = false;
  //     for(int i=0;i<f_size_cb;i++){
  // 	if( fabs( norm*rho_dopdwf[i] - rho[i]) > 1e-07 ){ printf("|Fail Mprec test Fbfm vs old code %d: %f %f\n",i,norm*rho_dopdwf[i], rho[i]); fail = true; }
  //     }
  //     if(fail){ printf("Failed Mprec test Fbfm vs old code test\n"); exit(-1); }
  //     else printf("Passed Mprec test Fbfm vs old code test\n");
  //   }

  //   //Net result is:
  //   // f_in = ( -kappa^2 Mprec chi, -kappa^2 D_eo^dag Mprec chi )    f_out = ( chi, D_eo chi )
  //   // where kappa = 1/[2(5-M5)]
  //   // where chi = psi1
  //   // These are converted into CANONICAL format when DiracOpDWF is destroyed
  //   dwf->CalcHmdForceVecs((Vector*)chi) ;
  //   delete dwf;

  //   /* Comparing the CPS and BFM HMD force vectors:
  //    * BFM :  v1 = ( Mprec chi, -1/[2(5-M5)] Deo^dag Mprec chi )     v2 = ( chi, -1/[2(5-M5)] * Deo chi )
  //    * CPS :  f_in = ( -1/[2(5-M5)]^2 Mprec chi, -1/[2(5-M5)]^2 D_eo Mprec chi )    f_out = ( chi, D_eo chi )
  //    * and using Mprec [CPS] = 1/(5-M5) Mprec [BFM]
  //    * CPS :  f_in = ( -1/[2(5-M5)]^2 1/(5-M5) Mprec[BFM] chi, -1/[2(5-M5)]^2 1/(5-M5) D_eo^dag Mprec [BFM] chi )    f_out = ( chi, D_eo chi )

  //    * We expect normalization differences   
  //    * v1 [ODD,BFM] = -[2(5-M5)]^2 (5-M5) f_in [ODD,CPS]
  //    * v1 [EVEN,BFM] = 2(5-M5)^2
  //    * v2 [ODD,BFM] = f_out [ODD,CPS]
  //    * v2 [EVEN,BFM] = -1/[2(5-M5)] f_out [EVEN,CPS]

  //    * BFM and CPS use 5D preconditioning, i.e. cb = (x+y+z+t+s)&0x1
  //    */

  //   fail = false;
  //   for(int i=0;i<f_size;i++){
  //     int rem = i;
  //     int midx = rem % 24; rem/=24;

  //     int x[5];
  //     for(int j=0;j<5;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

  //     int cb = ( x[0]+x[1]+x[2]+x[3]+x[4] )&0x1;
	    	    
  //     Float _5mM5 = 5-GJP.DwfHeight();
  //     Float norm_v1; //multiply CPS vector components by normalization factors
  //     Float norm_v2;

  //     if(cb == 1){ //odd part
  // 	norm_v1 = -4 * _5mM5 * _5mM5 * _5mM5;
  // 	norm_v2 = 1.0;
  //     }else{
  // 	norm_v1 = 2 * _5mM5 * _5mM5;
  // 	norm_v2 = -1.0/( 2 * _5mM5 );
  //     }

  //     if( fabs(norm_v1 * v1_dopdwf[i] - v1_fbfm[i]) > 1e-08 ){ 
  // 	printf("|Fail MD force vec test Fbfm vs old code v1 (%d %d %d %d %d, %d) [cb %d]: %f %f\n",x[0],x[1],x[2],x[3],x[4],
  // 	       midx, cb, norm_v1 * v1_dopdwf[i],v1_fbfm[i]); 
  // 	fail = true; 
  //     }
  //     if( fabs(norm_v2 * v2_dopdwf[i] - v2_fbfm[i]) > 1e-08 ){ 
  // 	printf("|Fail MD force vec test Fbfm vs old code v2 (%d %d %d %d %d, %d) [cb %d]: %f %f\n",x[0],x[1],x[2],x[3],x[4],
  // 	       midx, cb, norm_v2 * v2_dopdwf[i],v2_fbfm[i]); 
  // 	fail = true; 
  //     }
  //   }
  //   if(fail){ printf("Failed MD force vec Fbfm vs old code test\n"); exit(-1); }
  //   else printf("Passed MD force vec Fbfm vs old code test\n");


  //   //OK, try to calculate the same thing again but using the old G-parity DWF evolution code I used for the 16^3 lattice

  //   Float* mom_test6 = (Float *)pmalloc( sizeof(Float) * mom_size);
  //   for(int i=0;i<mom_size;i++) mom_test6[i] = 0.0;
	
  //   lattice->EvolveMomFforce( (Matrix*)mom_test6, (Vector*)chi, 0.5, 1.0); 
  //   Float _5mM5 = 5-GJP.DwfHeight();
  //   Float mom_norm = _5mM5 * _5mM5; //Hantao knows about this normalization difference; he did not attempt to reconcile the differences between the bfm and CPS fermion normalizations
    
  //   fail = false;
  //   for(int i=0;i<mom_size;i++){
  //     int rem = i;
  //     int midx = rem % 18; rem/=18;
  //     int mu = rem % 4; rem/=4;

  //     int x[4];
  //     for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

  //     if(fabs(mom_norm*mom_test6[i] - mom_test5[i])>1e-08){
  // 	printf("FforceWilsonType vs. old DWF code test fail midx %d mu %d (%d %d %d %d): %f %f\n",midx,mu,x[0],x[1],x[2],x[3],mom_norm*mom_test6[i],mom_test5[i]);
  // 	fail=true;
  //     }
  //   }
  //   if(fail){ printf("Failed FforceWilsonType vs. old DWF code test\n"); exit(-1); }
  //   else printf("Passed FforceWilsonType vs. old DWF code test\n");

  //   pfree(mom_test5);
  //   pfree(mom_test6);

  //   pfree(v1_fbfm);
  //   pfree(v2_fbfm);
  //   pfree(v1_dopdwf);
  //   pfree(v2_dopdwf);
  // }

  return 0;
}



#if 1
int main(int argc,char *argv[])
{
  Start(&argc,&argv); //initialises QMP

#ifdef HAVE_BFM
  Chroma::initialize(&argc,&argv);
#endif

  CommandLine::is(argc,argv);

  bool gparity_X(false);
  bool gparity_Y(false);

  int arg0 = CommandLine::arg_as_int(0);
  printf("Arg0 is %d\n",arg0);
  if(arg0==0){
    gparity_X=true;
    printf("Doing G-parity HMC test in X direction\n");
  }else if(arg0==1){
    printf("Doing G-parity HMC test in X and Y directions\n");
    gparity_X = true;
    gparity_Y = true;
  }else{
    printf("Doing No G-parity test\n");
  }

  bool dbl_latt_storemode(false);
  bool save_config(false);
  bool load_config(false);
  bool load_lrg(false);
  bool save_lrg(false);
  char *load_config_file;
  char *save_config_file;
  char *save_lrg_file;
  char *load_lrg_file;
  bool gauge_fix(false);
  bool verbose(false);
  bool skip_gparity_inversion(false);
  bool unit_gauge(false);

  BfmSolver solver = DWF;

  int size[] = {2,2,2,2,2};

  double min_fp_resid = 1e-05;

  int i=2;
  while(i<argc){
    char* cmd = argv[i];  
    if( strncmp(cmd,"-save_config",15) == 0){
      if(i==argc-1){ printf("-save_config requires an argument\n"); exit(-1); }
      save_config=true;
      save_config_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-load_config",15) == 0){
      if(i==argc-1){ printf("-save_config requires an argument\n"); exit(-1); }
      load_config=true;
      load_config_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-latt",10) == 0){
      if(i>argc-6){
	printf("Did not specify enough arguments for 'latt' (require 5 dimensions)\n"); exit(-1);
      }
      size[0] = CommandLine::arg_as_int(i); //CommandLine ignores zeroth input arg (i.e. executable name)
      size[1] = CommandLine::arg_as_int(i+1);
      size[2] = CommandLine::arg_as_int(i+2);
      size[3] = CommandLine::arg_as_int(i+3);
      size[4] = CommandLine::arg_as_int(i+4);
      i+=6;
    }else if( strncmp(cmd,"-save_double_latt",20) == 0){
      dbl_latt_storemode = true;
      i++;
    }else if( strncmp(cmd,"-load_lrg",15) == 0){
      if(i==argc-1){ printf("-load_lrg requires an argument\n"); exit(-1); }
      load_lrg=true;
      load_lrg_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-min_fp_resid",15) == 0){
      std::stringstream ss;
      ss << argv[i+1];
      ss >> min_fp_resid;
      if(UniqueID()) printf("Set minimum floating point residual for mixed-prec multi-mass shift to %e\n",min_fp_resid);
      i+=2;
    }else if( strncmp(cmd,"-save_lrg",15) == 0){
      if(i==argc-1){ printf("-save_lrg requires an argument\n"); exit(-1); }
      save_lrg=true;
      save_lrg_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-gauge_fix",15) == 0){
      gauge_fix=true;
      i++;   
    }else if( strncmp(cmd,"-verbose",15) == 0){
      verbose=true;
      i++;
    }else if( strncmp(cmd,"-skip_gparity_inversion",30) == 0){
      skip_gparity_inversion=true;
      i++;
    }else if( strncmp(cmd,"-unit_gauge",15) == 0){
      unit_gauge=true;
      i++;
    }else if( strncmp(cmd,"-mobius",15) == 0){
      solver= HmCayleyTanh;
      i++;
    }else{
      if(UniqueID()==0) printf("Unrecognised argument: %s\n",cmd);
      exit(-1);
    }
  }
  

  printf("Lattice size is %d %d %d %d\n",size[0],size[1],size[2],size[3],size[4]);

  DoArg do_arg;
  do_arg.x_sites = size[0];
  do_arg.y_sites = size[1];
  do_arg.z_sites = size[2];
  do_arg.t_sites = size[3];
  do_arg.s_sites = size[4];
  do_arg.x_node_sites = 0;
  do_arg.y_node_sites = 0;
  do_arg.z_node_sites = 0;
  do_arg.t_node_sites = 0;
  do_arg.s_node_sites = 0;
  do_arg.x_nodes = 0;
  do_arg.y_nodes = 0;
  do_arg.z_nodes = 0;
  do_arg.t_nodes = 0;
  do_arg.s_nodes = 0;
  do_arg.updates = 0;
  do_arg.measurements = 0;
  do_arg.measurefreq = 0;
  do_arg.cg_reprod_freq = 10;
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_conf_load_addr = 0x0;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.start_seed_filename = "../rngs/ckpoint_rng.0";
  do_arg.start_conf_filename = "../configurations/ckpoint_lat.0";
  do_arg.start_conf_alloc_flag = 6;
  do_arg.wfm_alloc_flag = 2;
  do_arg.wfm_send_alloc_flag = 2;
  do_arg.start_seed_value = 83209;
  do_arg.beta =   2.25;
  do_arg.c_1 =   -3.3100000000000002e-01;
  do_arg.u0 =   1.0000000000000000e+00;
  do_arg.dwf_height =   1.8000000000000000e+00;
  do_arg.dwf_a5_inv =   1.0000000000000000e+00;
  do_arg.power_plaq_cutoff =   0.0000000000000000e+00;
  do_arg.power_plaq_exponent = 0;
  do_arg.power_rect_cutoff =   0.0000000000000000e+00;
  do_arg.power_rect_exponent = 0;
  do_arg.verbose_level = -1202; //VERBOSE_DEBUG_LEVEL; //-1202;
  do_arg.checksum_level = 0;
  do_arg.exec_task_list = 0;
  do_arg.xi_bare =   1.0000000000000000e+00;
  do_arg.xi_dir = 3;
  do_arg.xi_v =   1.0000000000000000e+00;
  do_arg.xi_v_xi =   1.0000000000000000e+00;
  do_arg.clover_coeff =   0.0000000000000000e+00;
  do_arg.clover_coeff_xi =   0.0000000000000000e+00;
  do_arg.xi_gfix =   1.0000000000000000e+00;
  do_arg.gfix_chkb = 1;
  do_arg.asqtad_KS =   0.0000000000000000e+00;
  do_arg.asqtad_naik =   0.0000000000000000e+00;
  do_arg.asqtad_3staple =   0.0000000000000000e+00;
  do_arg.asqtad_5staple =   0.0000000000000000e+00;
  do_arg.asqtad_7staple =   0.0000000000000000e+00;
  do_arg.asqtad_lepage =   0.0000000000000000e+00;
  do_arg.p4_KS =   0.0000000000000000e+00;
  do_arg.p4_knight =   0.0000000000000000e+00;
  do_arg.p4_3staple =   0.0000000000000000e+00;
  do_arg.p4_5staple =   0.0000000000000000e+00;
  do_arg.p4_7staple =   0.0000000000000000e+00;
  do_arg.p4_lepage =   0.0000000000000000e+00;

  if(verbose) do_arg.verbose_level = VERBOSE_DEBUG_LEVEL;

  if(gparity_X) do_arg.x_bc = BND_CND_GPARITY;
  if(gparity_Y) do_arg.y_bc = BND_CND_GPARITY;

  GJP.Initialize(do_arg);

  SerialIO::dbl_latt_storemode = dbl_latt_storemode;
  
  LRG.Initialize(); //usually initialised when lattice generated, but I pre-init here so I can load the state from file

  if(load_lrg){
    if(UniqueID()==0) printf("Loading RNG state from %s\n",load_lrg_file);
    LRG.Read(load_lrg_file,32);
  }
  if(save_lrg){
    if(UniqueID()==0) printf("Writing RNG state to %s\n",save_lrg_file);
    LRG.Write(save_lrg_file,32);
  }
  
  GwilsonFdwf* lattice = new GwilsonFdwf;
					       
  if(!load_config){
    printf("Creating gauge field\n");
    if(!unit_gauge) lattice->SetGfieldDisOrd();
    else lattice->SetGfieldOrd();
  }else{
    ReadLatticeParallel readLat;
    if(UniqueID()==0) printf("Reading: %s (NERSC-format)\n",load_config_file);
    readLat.read(*lattice,load_config_file);
    if(UniqueID()==0) printf("Config read.\n");
  }

  if(save_config){
    if(UniqueID()==0) printf("Saving config to %s\n",save_config_file);

    QioArg wt_arg(save_config_file,0.001);
    
    wt_arg.ConcurIONumber=32;
    WriteLatticeParallel wl;
    wl.setHeader("disord_id","disord_label",0);
    wl.write(*lattice,wt_arg);
    
    if(!wl.good()) ERR.General("main","()","Failed write lattice %s",save_config_file);

    if(UniqueID()==0) printf("Config written.\n");
  }

  if(gauge_fix){
    lattice->FixGaugeAllocate(FIX_GAUGE_COULOMB_T);
    lattice->FixGauge(1e-06,2000);
    if(!UniqueID()){ printf("Gauge fixing finished\n"); fflush(stdout); }
  }
  cps_qdp_init(&argc,&argv);

  if(!gparity_X && !gparity_Y) return no_gparity_test(lattice, solver, min_fp_resid);
  
  return 0;


#if 0
  {
    //check mStarDotMTrans
    Float *a = new Float[18];
    Float *b = new Float[18];
    Float *c = new Float[18];
    Float *d = new Float[18];
    Float *e = new Float[18];

    for(int i=0;i<18;i++){ a[i] = LRG.Grand(); b[i] = LRG.Grand(); }
    
    // mTransDotMStarEqual
    bfm_evo_aux::mStarDotMTransEqual(c,a,b);
    Matrix ma; ma.Conj(a);
    Matrix mb; mb.Trans(b);
    Matrix mc; mc.DotMEqual(ma,mb);
    Float *mcp = (Float*)&mc[0];

    bool fail(false);
    for(int i=0;i<18;i++){
      if(mcp[i] != c[i]){ printf("Mprod check fail %d: %f %f\n",i,mcp[i],c[i]); fail =true; }
    }
    if(fail){ printf("Mprod check fail\n"); exit(-1); }
    else{ printf("Mprod check pass\n"); }

  }



  Float *v1;
  Float *v2;
  Float *mom_test1;

  Float *v1_test2;
  Float *v2_test2;
  
  Float *mom_test3;
  {
    bfmarg dwfa;
    setup_bfmargs(dwfa,solver);

    bfm_evo<double> bfm;
    bfm.init(dwfa);
    LatRanGen LRGbak(LRG);
    v1 = rand_5d_canonical_fermion(*lattice);
    v2 = rand_5d_canonical_fermion(*lattice);
    LRG = LRGbak;

    //v1 and v2 are in CPS canonical ordering, we need them in 's-ordering'
    Float* v1_sord;
    Float* v2_sord;
    convert_ferm_cpsord_sord(v1, v1_sord, bfm);
    convert_ferm_cpsord_sord(v2, v2_sord, bfm);

    int mom_size = 18*4*2*GJP.VolNodeSites();
    mom_test1 = (Float *)pmalloc( sizeof(Float) * mom_size);
    for(int i=0;i<mom_size;i++) mom_test1[i] = 0.0;

    int nthreads = 1; 
#if TARGET == BGQ
    nthreads = 64;
#endif
    omp_set_num_threads(nthreads);

    lattice->BondCond();
    Float* gauge = (Float*) lattice->GaugeField();
    bfm.cps_importGauge(gauge);  

    //first use random vectors to compute force from internal sites
#pragma omp parallel
    {
      int me = omp_get_thread_num();
      for(int mu=0;mu<4;mu++){
	bfm.fforce_internal(mom_test1, gauge, v1_sord, v2_sord, 0.1234, mu, me, nthreads); 
      }
    }

    //second test, try calculating MD force vectors
    //use v1 and v2 from test1 as the 'phi' vectors
    Fermion_t in[2] = {bfm.allocFermion(), bfm.allocFermion()};
    Fermion_t in2[2] = {bfm.allocFermion(), bfm.allocFermion()};
    bfm.cps_impexFermion(v1,in,1);
    bfm.cps_impexFermion(v2,in2,1);

    Fermion_t v1_tmp[2] = {bfm.allocFermion(), bfm.allocFermion()};
    Fermion_t v2_tmp[2] = {bfm.allocFermion(), bfm.allocFermion()};

    //zero the output vecs
    long f_size = (long)2*24 * GJP.VolNodeSites() * GJP.SnodeSites();
    Float *zeroferm = (Float *)pmalloc(sizeof(Float) * f_size);
    for(int i=0;i<f_size;i++) zeroferm[i] = 0.0;
    bfm.cps_impexFermion(zeroferm,v1_tmp,1);
    bfm.cps_impexFermion(zeroferm,v2_tmp,1);

    bfm.calcMDForceVecs(v1_tmp, v2_tmp, in[0], in2[0]); //just use left(=odd?) components

    v1_test2 = (Float *)pmalloc(sizeof(Float) * f_size);
    v2_test2 = (Float *)pmalloc(sizeof(Float) * f_size);
    for(int i=0;i<f_size;i++){v1_test2[i] = 0.0; v2_test2[i] = 0.0; }

    bfm.cps_impexFermion(v1_test2,v1_tmp,0);
    bfm.cps_impexFermion(v2_test2,v2_tmp,0);


    //OK, lets now try to calculate the full force. Again use v1 and v2 as the 'phi' vectors
    mom_test3 = (Float *)pmalloc( sizeof(Float) * mom_size);
    for(int i=0;i<mom_size;i++) mom_test3[i] = 0.0;

    bfm.compute_force(mom_test3, gauge, in[0], in2[0], 1.0);

    bfm.freeFermion(v1_tmp[0]);
    bfm.freeFermion(v1_tmp[1]);
    bfm.freeFermion(v2_tmp[0]);
    bfm.freeFermion(v2_tmp[1]);
    lattice->BondCond(); //undo lattice bcs


    //OK, Try to calculate the same conjugate momentum using the FforceWilsonType routines (developed later but based on the bfm_evo code I believe)
    //Much of the code is the same as that above
    {
      Float* mom_test4 = (Float *)pmalloc( sizeof(Float) * mom_size);
      for(int i=0;i<mom_size;i++) mom_test4[i] = 0.0;

      delete lattice; //temporarily delete the lattice object to remove the scope lock. The memory allocated to the gauge field remains
      Fbfm::bfm_arg = dwfa;
      GnoneFbfm *fbfm = new GnoneFbfm;
      //expects phi to be in cps checkerboard format. We can use in[0] and in2[0] as above, and export from bfm to cps format
      size_t f_size_cb = f_size / 2;
      Vector* phi1 = (Vector *)pmalloc(sizeof(Float) * f_size_cb);
      Vector* phi2 = (Vector *)pmalloc(sizeof(Float) * f_size_cb);

      fbfm->bd.cps_impexcbFermion((Float *)phi1, in[0], 0, 1); //export
      fbfm->bd.cps_impexcbFermion((Float *)phi2, in2[0], 0, 1);

      fbfm->EvolveMomFforceBase((Matrix*)mom_test4, phi1, phi2, 0.5, 1.0);

      //compare mom
      bool fail(false);
      for(int i=0;i<mom_size;i++){
	int rem = i;
	int midx = rem % 18; rem/=18;
	int mu = rem % 4; rem/=4;
	int flav = rem / GJP.VolNodeSites(); rem = rem % GJP.VolNodeSites();

	int x[4];
	for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

	if(fabs(mom_test4[i] - mom_test3[i])>1e-08){
	  printf("FforceWilsonType test fail flav %d midx %d mu %d (%d %d %d %d): %f %f\n",flav,midx,mu,x[0],x[1],x[2],x[3],mom_test4[i],mom_test3[i]);
	  fail=true;
	}
      }
      if(fail){ printf("Failed FforceWilsonType test\n"); exit(-1); }
      else printf("Passed FforceWilsonType test\n");
      pfree(mom_test4);

      if(solver == DWF){
	//test conversion between different formats
	//use v1 = CPS CANONICAL FORMAT
	Fermion_t v1_cb[2] = {fbfm->bd.allocFermion(), fbfm->bd.allocFermion()};
	fbfm->bd.cps_impexFermion(v1, v1_cb,1);

	//try extracting the odd and even checkerboarded parts and check them
	Float *v1_odd = (Float *)pmalloc(sizeof(Float)*f_size_cb);
	Float *v1_even = (Float *)pmalloc(sizeof(Float)*f_size_cb);
	fbfm->bd.cps_impexcbFermion(v1_odd, v1_cb[1],0,1);
	fbfm->bd.cps_impexcbFermion(v1_even, v1_cb[0],0,0);

	//check
	fail = false;
	for(int i=0;i<f_size;i++){
	  //G-parity CANONICAL ordering,  |   s=0   |   s=1   |
	  //                              | f0 | f1 | f0 | f1 |

	  int rem = i;
	  int midx = rem % 24; rem/=24;

	  int x[5];
	  x[4] = rem / (2*GJP.VolNodeSites()) + GJP.SnodeSites()*GJP.SnodeCoor();  rem = rem % (2*GJP.VolNodeSites());	  
	  int flav = rem / GJP.VolNodeSites(); rem = rem % GJP.VolNodeSites();

	  for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

	  int cb = ( x[0]+x[1]+x[2]+x[3]+x[4] )&0x1;
	  Float *vcb = v1_even;
	  if(cb==1) vcb = v1_odd;

	    //For G-parity the WILSON layout is 5d preconditioned
	    //|       s=0       |        s=1        | ......... |        s = 0      |.....
	    //| odd f0 | odd f1 | even f0 | even f1 | ......... | even f0 | even f1 |.....
	    //where the blocks on the lowest line have their *4d* parity indicated. (5d parity) = [(4d parity) + s] % 2
	    //hence the first half of the full WILSON vector had 5d parity odd, and the second half 5d parity even
	  
	  int halffourvol = GJP.VolNodeSites()/2;
	  int s_off = 2*halffourvol;
	  int cboff = midx + 24* ( (x[0]+ GJP.XnodeSites()*(x[1]+GJP.YnodeSites()*(x[2]+GJP.ZnodeSites()*x[3])))/2 + flav * halffourvol + x[4]*s_off );
	  
	  if( fabs(vcb[cboff] - v1[i]) > 1e-08 ){ 
	    if(midx == 0) printf("|Fail vec convert test flav %d (%d %d %d %d %d, %d) [cb %d]: %f %f\n",flav,x[0],x[1],x[2],x[3],x[4],
				 midx,cb, vcb[cboff],v1[i]); 
	    fail = true; 
	  }//else if(midx == 0) printf("Pass vec convert test flav %d (%d %d %d %d %d, %d) [cb %d]: %f %f\n",flav,x[0],x[1],x[2],x[3],x[4],
	  //		     midx,cb, vcb[cboff],v1[i]); 
	}
	if(fail){ printf("Failed vec convert test\n"); exit(-1); }
	else printf("Passed vec convert test\n");

	fbfm->bd.freeFermion(v1_cb[0]);
	fbfm->bd.freeFermion(v1_cb[1]);
	pfree(v1_odd);
	pfree(v1_even);
      }




      if(solver == DWF){
	//OK, run the complete EvolveMomFforce code that includes generating applying the M to phi1 to generate phi2
	Float *chi = (Float*)phi1;
	
	Float *v1_fbfm;
	Float *v2_fbfm;
	Float *rho = (Float *)pmalloc(sizeof(Float)*f_size_cb);
	{
	  Float *v1_fbfm_tmp = (Float *)pmalloc(sizeof(Float) * f_size);
	  Float *v2_fbfm_tmp = (Float *)pmalloc(sizeof(Float) * f_size);
	  
	  fbfm->MatPc( (Vector*)rho, (Vector*)chi, 0.5, DAG_NO); // rho = Moo chi
	  fbfm->CalcHmdForceVecsBilinear(v1_fbfm_tmp, v2_fbfm_tmp, (Vector*)rho, (Vector*)chi, 0.5); 

	  //convert to CPS canonical ordering
	  convert_ferm_sord_cpsord(v1_fbfm_tmp, v1_fbfm, fbfm->bd);
	  convert_ferm_sord_cpsord(v2_fbfm_tmp, v2_fbfm, fbfm->bd);
	  pfree(v1_fbfm_tmp);
	  pfree(v2_fbfm_tmp);
	}


	Float* mom_test5 = (Float *)pmalloc( sizeof(Float) * mom_size);
	for(int i=0;i<mom_size;i++) mom_test5[i] = 0.0;
	fbfm->EvolveMomFforce( (Matrix*)mom_test5, (Vector*)chi, 0.5, 1.0);

	delete fbfm;
	lattice = new GwilsonFdwf; //put the lattice back where it was
      

	{
	  CgArg cg_arg ;
	  cg_arg.mass = 0.5;

	  Float *v1_dopdwf = (Float *)pmalloc(sizeof(Float)*f_size);
	  Float *v2_dopdwf = (Float *)pmalloc(sizeof(Float)*f_size);
	  
	  /*  DiracOp(Lattice& latt, Vector *f_field_out, Vector *f_field_in, CgArg *arg,  CnvFrmType convert) */
	  DiracOpDwf* dwf = new DiracOpDwf(*lattice, (Vector*)v2_dopdwf, (Vector*)v1_dopdwf, &cg_arg, CNV_FRM_YES) ;
	  
	  {
	    //Test preconditioned matrix
	    Float *rho_dopdwf = (Float *)pmalloc(sizeof(Float)*f_size_cb);
	    dwf->MatPc( (Vector*)rho_dopdwf, (Vector*)chi);
	    
	    //Normalization factor: Mprec chi [CPS] = 1/(5-M5) Mprec chi [BFM]  (cf. no_gparity_test at top)
	    Float norm = 5-GJP.DwfHeight(); //multiply cps vector by norm

	    fail = false;
	    for(int i=0;i<f_size_cb;i++){
	      if( fabs( norm*rho_dopdwf[i] - rho[i]) > 1e-08 ){ printf("|Fail Mprec test Fbfm vs old code %d: %f %f\n",i,norm*rho_dopdwf[i], rho[i]); fail = true; }
	    }
	    if(fail){ printf("Failed Mprec test Fbfm vs old code test\n"); exit(-1); }
	    else printf("Passed Mprec test Fbfm vs old code test\n");
	  }

	  dwf->CalcHmdForceVecs((Vector*)chi) ;
	  delete dwf;


	  /* Comparing the CPS and BFM HMD force vectors:  (cf. no_gparity_test at top)
	   * BFM :  v1 = ( Mprec chi, -1/[2(5-M5)] Deo^dag Mprec chi )     v2 = ( chi, -1/[2(5-M5)] * Deo chi )
	   * CPS :  f_in = ( -1/[2(5-M5)]^2 Mprec chi, -1/[2(5-M5)]^2 D_eo Mprec chi )    f_out = ( chi, D_eo chi )
	   * and using Mprec [CPS] = 1/(5-M5) Mprec [BFM]
	   * CPS :  f_in = ( -1/[2(5-M5)]^2 1/(5-M5) Mprec[BFM] chi, -1/[2(5-M5)]^2 1/(5-M5) D_eo^dag Mprec [BFM] chi )    f_out = ( chi, D_eo chi )

	   * We expect normalization differences   
	   * v1 [ODD,BFM] = -[2(5-M5)]^2 (5-M5) f_in [ODD,CPS]
	   * v1 [EVEN,BFM] = 2(5-M5)^2
	   * v2 [ODD,BFM] = f_out [ODD,CPS]
	   * v2 [EVEN,BFM] = -1/[2(5-M5)] f_out [EVEN,CPS]

	   * BFM and CPS use 5D preconditioning, i.e. cb = (x+y+z+t+s)&0x1
	   */


	  fail = false;
	  for(int i=0;i<f_size;i++){
	    int rem = i;
	    int midx = rem % 24; rem/=24;

	    int x[5];
	    x[4] = rem / (2*GJP.VolNodeSites()) + GJP.SnodeSites()*GJP.SnodeCoor();  rem = rem % (2*GJP.VolNodeSites());	  
	    int flav = rem / GJP.VolNodeSites(); rem = rem % GJP.VolNodeSites();

	    for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

	    int cb = ( x[0]+x[1]+x[2]+x[3]+x[4] )&0x1;
	    	    
	    Float _5mM5 = 5-GJP.DwfHeight();
	    Float norm_v1 = 1.0;
	    Float norm_v2 = 1.0;

	    if(cb == 1){ //odd part
	      norm_v1 = -4 * _5mM5 * _5mM5 * _5mM5;
	      norm_v2 = 1.0;
	    }else{
	      norm_v1 = 2 * _5mM5 * _5mM5;
	      norm_v2 = -1.0/( 2 * _5mM5 );
	    }

	    if( fabs(norm_v1 * v1_dopdwf[i] - v1_fbfm[i]) > 1e-08 ){ 
	      if(midx == 0) printf("|Fail MD force vec test Fbfm vs old code v1 flav %d (%d %d %d %d %d, %d) [cb %d]: %f %f\n",flav,x[0],x[1],x[2],x[3],x[4],
		     midx,cb, norm_v1 * v1_dopdwf[i],v1_fbfm[i]); 
	      fail = true; 
	    }

	    if( fabs(norm_v2 * v2_dopdwf[i] - v2_fbfm[i]) > 1e-08 ){ 
	      if(midx == 0) printf("|Fail MD force vec test Fbfm vs old code v2 flav %d (%d %d %d %d %d, %d) [cb %d]: %f %f\n",flav,x[0],x[1],x[2],x[3],x[4],
		     midx,cb, norm_v2 * v2_dopdwf[i],v2_fbfm[i]); 
	      fail = true; 
	    }
	  }
	  if(fail){ printf("Failed MD force vec Fbfm vs old code test\n"); exit(-1); }
	  else printf("Passed MD force vec Fbfm vs old code test\n");

	  pfree(v1_fbfm);
	  pfree(v2_fbfm);
	  pfree(v1_dopdwf);
	  pfree(v2_dopdwf);
	}



	//OK, try to calculate the same thing again but using the old G-parity DWF evolution code I used for the 16^3 lattice

	Float* mom_test6 = (Float *)pmalloc( sizeof(Float) * mom_size);
	for(int i=0;i<mom_size;i++) mom_test6[i] = 0.0;
	
	lattice->EvolveMomFforce( (Matrix*)mom_test6, phi1, 0.5, 1.0); 
	Float _5mM5 = 5-GJP.DwfHeight();
	Float norm = _5mM5 * _5mM5; //Hantao knows about this normalization difference; he did not attempt to reconcile the differences between the bfm and CPS fermion normalizations

	/* Explanation:
	   The pseudofermion fields are generated by multiplying dslash to random vectors that follow the normal distribution. So the pseudofermion has these conventions built in. 
	*/

	fail = false;
	for(int i=0;i<mom_size;i++){
	  int rem = i;
	  int midx = rem % 18; rem/=18;
	  int mu = rem % 4; rem/=4;
	  int flav = rem / GJP.VolNodeSites(); rem = rem % GJP.VolNodeSites();

	  int x[4];
	  for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

	  if(fabs(norm*mom_test6[i] - mom_test5[i])>1e-08){
	    printf("FforceWilsonType vs. old DWF code test fail flav %d midx %d mu %d (%d %d %d %d): %f %f\n",flav,midx,mu,x[0],x[1],x[2],x[3],norm*mom_test6[i],mom_test5[i]);
	    fail=true;
	  }
	}
	if(fail){ printf("Failed FforceWilsonType vs. old DWF code test\n"); exit(-1); }
	else printf("Passed FforceWilsonType vs. old DWF code test\n");

	pfree(mom_test5);
	pfree(mom_test6);
      }else{
	delete fbfm;
	lattice = new GwilsonFdwf;
      }

      pfree(phi1);
      pfree(phi2);
    }

    bfm.freeFermion(in[0]);
    bfm.freeFermion(in[1]);
    bfm.freeFermion(in2[0]);
    bfm.freeFermion(in2[1]);


  }

  if(UniqueID()==0) printf("Starting double lattice section\n");
  
  int array_size = 2*lattice->GsiteSize() * GJP.VolNodeSites() * sizeof(Float);
  Matrix *orig_lattice = (Matrix *) pmalloc(array_size);
  memcpy((void*)orig_lattice, (void*)lattice->GaugeField(), array_size);

  lattice->FreeGauge(); //free memory and reset
  delete lattice; //lattice objects are singleton (scope_lock)
  
  //setup 1f model. Upon calling GJP.Initialize the lattice size will be doubled in the appropriate directions
  //and the boundary condition set to APRD
  if(gparity_X) do_arg.gparity_1f_X = 1;
  if(gparity_Y) do_arg.gparity_1f_Y = 1;

  GJP.Initialize(do_arg);

  if(GJP.Gparity()){ printf("Que?\n"); exit(-1); }
  if(UniqueID()==0) printf("Doubled lattice : %d %d %d %d\n", GJP.XnodeSites()*GJP.Xnodes(),GJP.YnodeSites()*GJP.Ynodes(),
			   GJP.ZnodeSites()*GJP.Znodes(),GJP.TnodeSites()*GJP.Tnodes());
  
#ifdef HAVE_BFM
  {
    QDP::multi1d<int> nrow(Nd);  
    for(int i = 0;i<Nd;i++) nrow[i] = GJP.Sites(i);
    //  multi1d<LatticeFermion> test(Nd);  
    //  nrow=size;
    QDP::Layout::setLattSize(nrow);
    QDP::Layout::create();
  }
#endif
  //cps_qdp_init(argc,argv);

  GwilsonFdwf *doubled_lattice = new GwilsonFdwf;
  setup_double_latt(*doubled_lattice,orig_lattice,gparity_X,gparity_Y);
  setup_double_rng(gparity_X,gparity_Y);
 
  if(gauge_fix){
    doubled_lattice->FixGaugeAllocate(FIX_GAUGE_COULOMB_T);
    doubled_lattice->FixGauge(1e-06,2000);
    if(!UniqueID()){ printf("Gauge fixing finished\n"); fflush(stdout); }
  }
 
  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();

  //convert the random CANONICAL vectors from the 2f setup to 1f form
  Float *v1_dbl = (Float *)pmalloc(sizeof(Float) * f_size);
  setup_double_5d_vector((Vector*)v1_dbl,(Vector*)v1, gparity_X,gparity_Y);
  pfree(v1);

  Float *v2_dbl = (Float *)pmalloc(sizeof(Float) * f_size);
  setup_double_5d_vector((Vector*)v2_dbl,(Vector*)v2, gparity_X,gparity_Y);
  pfree(v2);

  Float *v1_test2_dbl = (Float *)pmalloc(sizeof(Float) * f_size);
  setup_double_5d_vector((Vector*)v1_test2_dbl,(Vector*)v1_test2, gparity_X,gparity_Y);
  pfree(v1_test2);

  Float *v2_test2_dbl = (Float *)pmalloc(sizeof(Float) * f_size);
  setup_double_5d_vector((Vector*)v2_test2_dbl,(Vector*)v2_test2, gparity_X,gparity_Y);
  pfree(v2_test2);


  Float *mom_test1_dbl = (Float *) pmalloc(GJP.VolNodeSites()*18*4*sizeof(Float));
  setup_double_matrixfield((Matrix*)mom_test1_dbl, (Matrix*)mom_test1, 4, gparity_X, gparity_Y);
  pfree(mom_test1);

  Float *mom_test3_dbl = (Float *) pmalloc(GJP.VolNodeSites()*18*4*sizeof(Float));
  setup_double_matrixfield((Matrix*)mom_test3_dbl, (Matrix*)mom_test3, 4, gparity_X, gparity_Y);
  pfree(mom_test3);

  {
    long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();
    Fermion_t in[2];
    Fermion_t in2[2];

    Float* mom_test3_dbllatt;

    bfmarg dwfa;
    {
      setup_bfmargs(dwfa,solver);

      bfm_evo<double> bfm;
      bfm.init(dwfa);

      //v1 and v2 are in CPS canonical ordering, we need them in 's-ordering'
      Float* v1_sord;
      Float* v2_sord;
      convert_ferm_cpsord_sord(v1_dbl, v1_sord, bfm);
      convert_ferm_cpsord_sord(v2_dbl, v2_sord, bfm);

      int mom_size = 18*4*GJP.VolNodeSites();
      Float* mom_test1_dbllat = (Float *)pmalloc( sizeof(Float) * mom_size);
      for(int i=0;i<mom_size;i++) mom_test1_dbllat[i] = 0.0;

      doubled_lattice->BondCond();
      Float* gauge = (Float*) doubled_lattice->GaugeField();
      bfm.cps_importGauge(gauge);

      int nthreads = 1; 
#if TARGET == BGQ
      nthreads = 64;
#endif
      omp_set_num_threads(nthreads);

      //first use random vectors to compute force from internal sites
#pragma omp parallel
      {
	int me = omp_get_thread_num();
	for(int mu=0;mu<4;mu++){
	  bfm.fforce_internal(mom_test1_dbllat, gauge, v1_sord, v2_sord, 0.1234, mu, me, nthreads); 
	}
      }

      bool fail(false);
    
      if(!gparity_Y){ //skip this test for quad lattice
	//compare mom
	for(int i=0;i<mom_size;i++){
	  int rem = i;
	  int midx = rem % 18; rem/=18;
	  int mu = rem % 4; rem/=4;
	  int x[4];
	  for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

	  //printf("XnodeSites is %d, XnodeSites()/2 is %d, this site x coord is %d\n",GJP.XnodeSites(), GJP.XnodeSites()/2 , x[0]);

	  if(GJP.Xnodes() == 1 && (x[0] == GJP.XnodeSites()/2-1 || x[0] == GJP.XnodeSites()-1)) continue; //these were surface sites before we doubled the lattice

	  if(fabs(mom_test1_dbllat[i] - mom_test1_dbl[i])>1e-08){
	    printf("Mom test 1 fail midx %d mu %d (%d %d %d %d): %f %f\n",midx,mu,x[0],x[1],x[2],x[3],mom_test1_dbllat[i],mom_test1_dbl[i]);
	    fail=true;
	  }//else  printf("Mom test 1 pass midx %d mu %d (%d %d %d %d): %f %f\n",midx,mu,x[0],x[1],x[2],x[3],mom_test1_dbllat[i],mom_test1_dbl[i]);
	}
	if(fail){ printf("Failed mom test 1\n"); exit(-1); }
	else printf("Passed mom test 1\n");
      }



      //second test, try calculating MD force vectors
      //use v1 and v2 from test1 as the 'phi' vectors
      in[0] = bfm.allocFermion();       in[1] = bfm.allocFermion();
      in2[0] = bfm.allocFermion();       in2[1] = bfm.allocFermion();

      bfm.cps_impexFermion(v1_dbl,in,1);
      bfm.cps_impexFermion(v2_dbl,in2,1);

      Fermion_t v1_tmp[2] = {bfm.allocFermion(), bfm.allocFermion()};
      Fermion_t v2_tmp[2] = {bfm.allocFermion(), bfm.allocFermion()};
      //zero the output vecs
      
      Float *zeroferm = (Float *)pmalloc(sizeof(Float) * f_size);
      for(int i=0;i<f_size;i++) zeroferm[i] = 0.0;
      bfm.cps_impexFermion(zeroferm,v1_tmp,1);
      bfm.cps_impexFermion(zeroferm,v2_tmp,1);

      bfm.calcMDForceVecs(v1_tmp, v2_tmp, in[0], in2[0]);

      Float* v1_test2_dbllatt = (Float *)pmalloc(sizeof(Float) * f_size);
      Float* v2_test2_dbllatt = (Float *)pmalloc(sizeof(Float) * f_size);
      for(int i=0;i<f_size;i++){v1_test2_dbllatt[i] = 0.0; v2_test2_dbllatt[i] = 0.0; }

      bfm.cps_impexFermion(v1_test2_dbllatt,v1_tmp,0);
      bfm.cps_impexFermion(v2_test2_dbllatt,v2_tmp,0);

      bfm.freeFermion(v1_tmp[0]);
      bfm.freeFermion(v1_tmp[1]);
      bfm.freeFermion(v2_tmp[0]);
      bfm.freeFermion(v2_tmp[1]);
    
      fail = false;
      for(int i=0;i<f_size;i++){
	int rem = i;
	int midx = rem % 24; rem/=24;
	int x[5];
	for(int j=0;j<5;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }
      
	if( fabs(v1_test2_dbllatt[i] - v1_test2_dbl[i]) > 1e-08 ){ printf("|Fail MD force vec test v1 (%d %d %d %d %d, %d): %f %f\n",x[0],x[1],x[2],x[3],x[4],midx,v1_test2_dbllatt[i],v1_test2_dbl[i]); fail = true; }
	//else printf("Pass MD force vec test v1 (%d %d %d %d %d, %d): %f %f\n",x[0],x[1],x[2],x[3],x[4],midx,v1_test2_dbllatt[i],v1_test2_dbl[i]);

	if( fabs(v2_test2_dbllatt[i] - v2_test2_dbl[i]) > 1e-08 ){ printf("|Fail MD force vec test v2 (%d %d %d %d %d, %d): %f %f\n",x[0],x[1],x[2],x[3],x[4],midx,v2_test2_dbllatt[i],v2_test2_dbl[i]); fail = true; }
	//else printf("Pass MD force vec test v2 (%d %d %d %d %d, %d): %f %f\n",x[0],x[1],x[2],x[3],x[4],midx,v2_test2_dbllatt[i],v2_test2_dbl[i]);

      }
      if(fail){ printf("Failed MD force vec test\n"); exit(-1); }
      else printf("Passed MD force vec test\n");


      mom_test3_dbllatt = (Float *)pmalloc( sizeof(Float) * mom_size);
      for(int i=0;i<mom_size;i++) mom_test3_dbllatt[i] = 0.0;

      bfm.compute_force(mom_test3_dbllatt, gauge, in[0], in2[0], 1.0);


      for(int i=0;i<mom_size;i++){
	int rem = i;
	int midx = rem % 18; rem/=18;
	int mu = rem % 4; rem/=4;
	int x[4];
	for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

	if(fabs(mom_test3_dbllatt[i] - mom_test3_dbl[i])>1e-08){
	  printf("Mom test 3 fail midx %d mu %d (%d %d %d %d): %f %f\n",midx,mu,x[0],x[1],x[2],x[3],mom_test3_dbllatt[i],mom_test3_dbl[i]);
	  fail=true;
	}
      }
      if(fail){ printf("Failed mom test 3\n"); exit(-1); }
      else printf("Passed mom test 3\n");


      // bfm.freeFermion(in[0]);
      // bfm.freeFermion(in[1]);
      // bfm.freeFermion(in2[0]);
      // bfm.freeFermion(in2[1]);
      doubled_lattice->BondCond();
    }



    //OK, Try to calculate the same conjugate momentum using the FforceWilsonType routines (developed later but based on the bfm_evo code I believe)
    //Much of the code is the same as that above
    {
      int mom_size = 18*4*GJP.VolNodeSites();
      Float* mom_test4_dbllatt = (Float *)pmalloc( sizeof(Float) * mom_size);
      for(int i=0;i<mom_size;i++) mom_test4_dbllatt[i] = 0.0;

      delete doubled_lattice; //temporarily delete the lattice object to remove the scope lock. The memory allocated to the gauge field remains
      Fbfm::bfm_arg = dwfa;
      GnoneFbfm *fbfm = new GnoneFbfm;
      //expects phi to be in cps checkerboard format. We can use in[0] and in2[0] as above, and export from bfm to cps format
      size_t f_size_cb = f_size / 2;
      Vector* phi1 = (Vector *)pmalloc(sizeof(Float) * f_size_cb);
      Vector* phi2 = (Vector *)pmalloc(sizeof(Float) * f_size_cb);

      fbfm->bd.cps_impexcbFermion((Float *)phi1, in[0], 0, 1); //export
      fbfm->bd.cps_impexcbFermion((Float *)phi2, in2[0], 0, 1);

      fbfm->EvolveMomFforceBase((Matrix*)mom_test4_dbllatt, phi1, phi2, 0.5, 1.0);

      //compare mom
      bool fail(false);
      for(int i=0;i<mom_size;i++){
	int rem = i;
	int midx = rem % 18; rem/=18;
	int mu = rem % 4; rem/=4;

	int x[4];
	for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

	if(fabs(mom_test4_dbllatt[i] - mom_test3_dbllatt[i])>1e-08){
	  printf("FforceWilsonType double latt test fail midx %d mu %d (%d %d %d %d): %f %f\n",midx,mu,x[0],x[1],x[2],x[3],mom_test4_dbllatt[i],mom_test3_dbllatt[i]);
	  fail=true;
	}
      }
      if(fail){ printf("Failed FforceWilsonType double latt test\n"); exit(-1); }
      else printf("Passed FforceWilsonType double latt test\n");
      pfree(mom_test4_dbllatt);
    }
  }


#ifdef HAVE_BFM
  Chroma::finalize();
#endif

  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  
  return 0;
  #endif
}




void GaugeTransformU(Matrix *gtrans, Lattice &lat){
  Matrix recv_buf;
  Matrix tmp;
  //apply the gauge transformation to U
  int nflav = 1;
  if(GJP.Gparity()) nflav = 2;

  for(int flav=0;flav<nflav;flav++){
    for(int t=0;t<GJP.TnodeSites();t++){
      for(int z=0;z<GJP.ZnodeSites();z++){
	for(int y=0;y<GJP.YnodeSites();y++){
	  for(int x=0;x<GJP.XnodeSites();x++){
	    int pos[4] = {x,y,z,t};
	    int v_x_off = x + GJP.XnodeSites()*(y+GJP.YnodeSites()*(z+GJP.ZnodeSites()*t)) + flav*GJP.VolNodeSites();
	    Matrix &v_x = *(gtrans + v_x_off);

	    for(int mu=0;mu<4;mu++){
	      int u_x_off = lat.GsiteOffset(pos) + mu + flav*4*GJP.VolNodeSites();
	      Matrix &u_x = *(lat.GaugeField() + u_x_off);

	      //get V_x+mu
	      int posp[4] = {x,y,z,t};
	      posp[mu] = (posp[mu]+1)%GJP.NodeSites(mu);

	      Matrix *v_xpmu_ptr = gtrans + posp[0] + GJP.XnodeSites()*(posp[1]+GJP.YnodeSites()*(posp[2]+GJP.ZnodeSites()*posp[3])) + flav*GJP.VolNodeSites();
	      if(pos[mu] == GJP.NodeSites(mu)-1){
		//if node is on the left wall, send the opposite flavour 
		if(GJP.Bc(mu) == BND_CND_GPARITY && GJP.NodeCoor(mu) == 0){
		  if(flav == 1)
		    v_xpmu_ptr-= GJP.VolNodeSites();
		  else
		    v_xpmu_ptr+= GJP.VolNodeSites();		  
		}

		//doesnt need to be fast!
		getPlusData((double *)&recv_buf, (double *)v_xpmu_ptr, 18, mu);
		v_xpmu_ptr = &recv_buf; 
	      }

	      //dagger/transpose it
	      Matrix vdag_xpmu;
	      vdag_xpmu.Dagger(*v_xpmu_ptr);

	      //gauge transform link
	      tmp.DotMEqual(v_x,u_x);
	      u_x.DotMEqual(tmp,vdag_xpmu);
	    }
	  }
	}
      }
    }

  }

}

#endif
