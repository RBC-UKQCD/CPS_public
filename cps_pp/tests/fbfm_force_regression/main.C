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

#include<omp.h>

#ifdef HAVE_BFM
#include <chroma.h>
#endif

using namespace std;
USING_NAMESPACE_CPS

static void convert_ferm_cpsord_sord(Float *cps, Float* &sord, bfm_evo<Float> &bfm){
  Fermion_t handle[2] = { bfm.allocFermion(), bfm.allocFermion() };
  bfm.cps_impexFermion(cps,handle,1);
  
  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();
  sord = (Float *)pmalloc(sizeof(Float) * f_size);
  bfm.cps_impexFermion_s(sord,handle,0);

  bfm.freeFermion(handle[0]);
  bfm.freeFermion(handle[1]);
}

static void convert_ferm_sord_cpsord(Float *sord, Float* &cps, bfm_evo<Float> &bfm){
  Fermion_t handle[2] = { bfm.allocFermion(), bfm.allocFermion() };
  bfm.cps_impexFermion_s(sord,handle,1);
  
  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();
  cps = (Float *)pmalloc(sizeof(Float) * f_size);
  bfm.cps_impexFermion(cps,handle,0);

  bfm.freeFermion(handle[0]);
  bfm.freeFermion(handle[1]);
}

void setup_bfmargs(bfmarg &dwfa){
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
  dwfa.Ls   = GJP.SnodeSites();
  dwfa.solver = DWF;
  dwfa.M5   = toDouble(GJP.DwfHeight());
  dwfa.mass = toDouble(0.5);
  dwfa.Csw  = 0.0;
  dwfa.max_iter = 5000;
  dwfa.residual = 1e-08;
  printf("Finished setting up bfmargs\n");
}

Float* rand_5d_canonical_fermion(Lattice &lat){
  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();
  Float *v1 = (Float *)pmalloc(sizeof(Float) * f_size);
  printf("Making random gaussian 5d vector\n");
  lat.RandGaussVector((Vector*)v1, 0.5, 2, CANONICAL, FIVE_D);
  printf("Finished making random gaussian vector\n");
  return v1;
}

static void hmd_force_vec_test(GwilsonFdwf* &lattice){
  //Try using Fbfm and old DWF code to get the same MD force vectors
  bool fail;

  bfmarg dwfa;
  setup_bfmargs(dwfa);
  long mom_size = (long)18*4*GJP.VolNodeSites();
  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();

  delete lattice; //temporarily delete the lattice object to remove the scope lock. The memory allocated to the gauge field remains
  Fbfm::bfm_arg = dwfa;
  GnoneFbfm *fbfm = new GnoneFbfm;
  
  //Generate a random CANONICAL ordered fermion
  LatRanGen LRGbak(LRG);
  Float *v1 = rand_5d_canonical_fermion(*fbfm);
  LRG = LRGbak;

  //Convert to bfm format
  Fermion_t in[2] = {fbfm->bd.allocFermion(), fbfm->bd.allocFermion()};
  fbfm->bd.cps_impexFermion(v1,in,1);

  //expects phi to be in cps checkerboard format. We can extract this easily by exporting from the bfm 'in' object to cps format
  size_t f_size_cb = f_size / 2;
  
  Vector* chi = (Vector *)pmalloc(sizeof(Float) * f_size_cb); //ODD vector
  fbfm->bd.cps_impexcbFermion((Float *)chi, in[1], 0, 1); 

  Float *v1_fbfm;
  Float *v2_fbfm;

  {
    Float *v1_fbfm_tmp = (Float *)pmalloc(sizeof(Float) * f_size);
    Float *v2_fbfm_tmp = (Float *)pmalloc(sizeof(Float) * f_size);

    Float *rho = (Float *)pmalloc(sizeof(Float)*f_size_cb);    
    fbfm->MatPc( (Vector*)rho, (Vector*)chi, 0.5, DAG_NO); // rho [ODD] = Mprec chi [ODD]
    fbfm->CalcHmdForceVecsBilinear(v1_fbfm_tmp, v2_fbfm_tmp, (Vector*)rho, (Vector*)chi, 0.5); 
    
    /*
     * From BFM:
     * DWF : Mee = Moo = (5-M5)
     *       Meo = -1/2 Ddwf_eo (5d hopping term)
     *       Moe = -1/2 Ddwf_oe (5d hopping term)
     * Mprec = Moo-MoeMee^{-1}Meo
     */
    //CalcHmdForceVecsBilinear calculates the following:
    // v2e      =  Bee * 1/(5-M5) * Meo phi2
    // v2o = Boo phi2
    // v1e  =  1/(5-M5) Meo^dag phi1
    // v1o = 1oo phi1
    //For DWF, Boo = Bee = 1 and phi1 = rho,  phi2 = chi in the above naming convention
    //hence the net result is
    //v1 = ( rho, 1/(5-M5) Meo^dag rho )     v2 = ( chi, 1/(5-M5) * Meo chi )
    //i.e.
    //v1 = ( rho, -1/[2(5-M5)] Deo^dag rho )     v2 = ( chi, -1/[2(5-M5)] * Deo chi )
    //with rho = Mprec chi
    //v1 = ( Mprec chi, -1/[2(5-M5)] Deo^dag Mprec chi )     v2 = ( chi, -1/[2(5-M5)] * Deo chi )
    //  ((Above vectors are in  (odd,even) format))
    
    //convert to CPS canonical ordering
    convert_ferm_sord_cpsord(v1_fbfm_tmp, v1_fbfm, fbfm->bd);
    convert_ferm_sord_cpsord(v2_fbfm_tmp, v2_fbfm, fbfm->bd);
    pfree(v1_fbfm_tmp);
    pfree(v2_fbfm_tmp);
  }

  
  delete fbfm;
  lattice = new GwilsonFdwf; //put the lattice back where it was
 

  {
    CgArg cg_arg ;
    cg_arg.mass = 0.5;
    
    Float *v1_dopdwf = (Float *)pmalloc(sizeof(Float)*f_size);
    Float *v2_dopdwf = (Float *)pmalloc(sizeof(Float)*f_size);
    
    /*  DiracOp(Lattice& latt, Vector *f_field_out, Vector *f_field_in, CgArg *arg,  CnvFrmType convert) */
    DiracOpDwf* dwf = new DiracOpDwf(*lattice, (Vector*)v2_dopdwf, (Vector*)v1_dopdwf, &cg_arg, CNV_FRM_YES) ;
    
    //Net result is:
    // f_in = ( -kappa^2 Mprec chi, -kappa^2 D_eo^dag Mprec chi )    f_out = ( chi, D_eo chi )
    // where kappa = 1/[2(5-M5)]
    // where chi = psi1
    // These are converted into CANONICAL format when DiracOpDWF is destroyed
    dwf->CalcHmdForceVecs((Vector*)chi) ;
    delete dwf;

    /* Comparing the CPS and BFM HMD force vectors:
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
      for(int j=0;j<5;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

      int cb = ( x[0]+x[1]+x[2]+x[3]+x[4] )&0x1;
	    	    
      Float _5mM5 = 5-GJP.DwfHeight();
      Float norm_v1; //multiply CPS vector components by normalization factors
      Float norm_v2;

      if(cb == 1){ //odd part
	norm_v1 = -4 * _5mM5 * _5mM5 * _5mM5;
	norm_v2 = 1.0;
      }else{
	norm_v1 = 2 * _5mM5 * _5mM5;
	norm_v2 = -1.0/( 2 * _5mM5 );
      }

      if( fabs(norm_v1 * v1_dopdwf[i] - v1_fbfm[i]) > 1e-08 ){ 
	printf("|Fail MD force vec test Fbfm vs old code v1 (%d %d %d %d %d, %d) [cb %d]: %f %f\n",x[0],x[1],x[2],x[3],x[4],
	       midx, cb, norm_v1 * v1_dopdwf[i],v1_fbfm[i]); 
	fail = true; 
      }
      if( fabs(norm_v2 * v2_dopdwf[i] - v2_fbfm[i]) > 1e-08 ){ 
	printf("|Fail MD force vec test Fbfm vs old code v2 (%d %d %d %d %d, %d) [cb %d]: %f %f\n",x[0],x[1],x[2],x[3],x[4],
	       midx, cb, norm_v2 * v2_dopdwf[i],v2_fbfm[i]); 
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
}


static void force_test(GwilsonFdwf* &lattice){
  //Try using Fbfm and old DWF code to get the same MD force
  bool fail;

  bfmarg dwfa;
  setup_bfmargs(dwfa);
  long mom_size = (long)18*4*GJP.VolNodeSites();
  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();

  delete lattice; //temporarily delete the lattice object to remove the scope lock. The memory allocated to the gauge field remains
  Fbfm::bfm_arg = dwfa;
  GnoneFbfm *fbfm = new GnoneFbfm;
  
  //Generate a random CANONICAL ordered fermion
  LatRanGen LRGbak(LRG);
  Float *v1 = rand_5d_canonical_fermion(*fbfm);
  LRG = LRGbak;

  //Convert to bfm format
  Fermion_t in[2] = {fbfm->bd.allocFermion(), fbfm->bd.allocFermion()};
  fbfm->bd.cps_impexFermion(v1,in,1);

  //expects phi to be in cps checkerboard format. We can extract this easily by exporting from the bfm 'in' object to cps format
  size_t f_size_cb = f_size / 2;
  
  Vector* chi = (Vector *)pmalloc(sizeof(Float) * f_size_cb); //ODD vector
  fbfm->bd.cps_impexcbFermion((Float *)chi, in[1], 0, 1); 

  //Calculate force in bfm
  Float* mom_test_bfm = (Float *)pmalloc( sizeof(Float) * mom_size);
  for(int i=0;i<mom_size;i++) mom_test_bfm[i] = 0.0;
  fbfm->EvolveMomFforce( (Matrix*)mom_test_bfm, (Vector*)chi, 0.5, 1.0);

  delete fbfm;
  lattice = new GwilsonFdwf; //put the lattice back where it was
 

  Float* mom_test_cps = (Float *)pmalloc( sizeof(Float) * mom_size);
  for(int i=0;i<mom_size;i++) mom_test_cps[i] = 0.0;
	
  lattice->EvolveMomFforce( (Matrix*)mom_test_cps, (Vector*)chi, 0.5, 1.0); 
  
  //NORMALIZATION WEIRDNESS: I asked Hantao and he says this is expected; he did not attempt to reconcile the differences between the BFM and CPS normalizations
  Float _5mM5 = 5-GJP.DwfHeight();
  Float mom_norm = _5mM5 * _5mM5;
    
  fail = false;
  for(int i=0;i<mom_size;i++){
    int rem = i;
    int midx = rem % 18; rem/=18;
    int mu = rem % 4; rem/=4;

    int x[4];
    for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

    if(fabs(mom_norm*mom_test_cps[i] - mom_test_bfm[i])>1e-08){
      printf("FforceWilsonType vs. old DWF code test fail midx %d mu %d (%d %d %d %d): %f %f\n",midx,mu,x[0],x[1],x[2],x[3],mom_test_bfm[i],mom_norm*mom_test_cps[i]);
      fail=true;
    }
  }
  if(fail){ printf("Failed FforceWilsonType vs. old DWF code test\n"); exit(-1); }
  else printf("Passed FforceWilsonType vs. old DWF code test\n");

  pfree(mom_test_bfm);
  pfree(mom_test_cps);
  pfree(chi);
  pfree(v1);
}





int main(int argc,char *argv[])
{
  Start(&argc,&argv); //initialises QMP

#ifdef HAVE_BFM
  Chroma::initialize(&argc,&argv);
#endif

  CommandLine::is(argc,argv);

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
  bool unit_gauge(false);

  int size[] = {2,2,2,2,2};

  int i=1;
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
    }else if( strncmp(cmd,"-load_lrg",15) == 0){
      if(i==argc-1){ printf("-load_lrg requires an argument\n"); exit(-1); }
      load_lrg=true;
      load_lrg_file = argv[i+1];
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
    }else if( strncmp(cmd,"-unit_gauge",15) == 0){
      unit_gauge=true;
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


  GJP.Initialize(do_arg);
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

  hmd_force_vec_test(lattice);
  force_test(lattice);

  delete lattice;
#ifdef HAVE_BFM
  Chroma::finalize();
#endif

  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  
  return 0;
}

