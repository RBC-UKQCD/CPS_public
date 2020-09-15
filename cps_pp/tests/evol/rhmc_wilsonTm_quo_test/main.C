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
#include<util/dirac_op.h>

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

// define all integrators; vml file arguments will not change
IntABArg ab1_arg;
IntABArg ab2_arg;

EvoArg evo_arg;
DoArg do_arg;
PbpArg pbp_arg;
NoArg no_arg;

void checkpoint(int traj);

#define decode_vml(arg_name)\
  printf("Decoding %s\n",#arg_name);\
  do{                                       \
        if ( ! arg_name.Decode(#arg_name".vml", #arg_name) )            \
            ERR.General(cname, fname, "Bad " #arg_name ".vml.\n");      \
    } while(0)

void decode_vml_all(void)
{
    char *fname = "decode_vml_all()";

    decode_vml(do_arg);
    decode_vml(hmc_arg);
    decode_vml(evo_arg);
    decode_vml(gauge_arg);
    decode_vml(quo_arg);
    decode_vml(ab1_arg);
    decode_vml(ab2_arg);
    decode_vml(pbp_arg);
}

void truncate_it(CommonArg *common_arg, const char stem[], int traj);
void measure_plaq(CommonArg &common_arg);
void measure_pbp(CommonArg &common_arg, int traj);
void run_hmc(CommonArg &common_arg, int traj, AlgIntAB &int_ab);


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

    if(chdir(evo_arg.work_directory) != 0) {
        ERR.General(cname, fname, "Changing directory to %s failed.\n", evo_arg.work_directory);
    }
    VRB.Result(cname, fname, "Reading VML files successfully.\n");

    GJP.Initialize(do_arg);
    //LRG.Initialize();

    VRB.Result(cname, fname, "VRB.Level(%d)\n", do_arg.verbose_level);
    VRB.Level(do_arg.verbose_level);

    cps_qdp_init(&argc,&argv);
}

Float* rand_4d_canonical_fermion(Lattice &lat){
  long f_size = (long)24 * GJP.VolNodeSites();
  Float *v1 = (Float *)pmalloc(sizeof(Float) * f_size);
  printf("Making random gaussian 4d vector\n");
  lat.RandGaussVector((Vector*)v1, 0.5, 2, CANONICAL, FOUR_D);
  printf("Finished making random gaussian vector\n");
  return v1;
}


 // CgArg cg_arg;
 //  cg_arg.mass = -1.8;
 //  cg_arg.epsilon = 0.5;
 //  cg_arg.max_num_iter = 10000;
 //  cg_arg.stop_rsd = 1e-08;
 //  cg_arg.true_rsd = 1e-08;
 //  cg_arg.Inverter = CG;

int InvCgShift_CPS(Vector *out, 
		   Vector *in, 
		   Float src_norm_sq, 
		   Float *true_res,
		   Float *shift,
		   DiracOp &dop,
		   Lattice& lat, CgArg &cg_arg){
  int itr;                       // Current number of CG iterations
  int max_itr;                       // Max number of CG iterations
  Float stp_cnd;                   // Stop if residual^2 <= stp_cnd
  Float res_norm_sq_prv;          // The previous step |residual|^2
  Float res_norm_sq_cur;           // The current step |residual|^2 
  Float a;
  Float b;
  Float d;
  int i, ic, icb;
  const char *fname = "InvCgShift(V*,V*,F,F*) [Duplicate of d_op_base/noarch version]";
  IFloat *temp;


// Print out input parameters
//------------------------------------------------------------------
  VRB.Input(cname,fname,
	    "stop_rsd = %e\n",IFloat(cg_arg.stop_rsd));
  VRB.Input(cname,fname,
	    "max_num_iter = %d\n",cg_arg.max_num_iter);
  VRB.Input(cname,fname,
	    "mass = %e\n",IFloat(cg_arg.mass));
  VRB.Input(cname,fname,
	    "src_norm_sq = %e\n",IFloat(src_norm_sq));


//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------

// Set the source vector pointer
//------------------------------------------------------------------
  Vector *src = in;

// Set the solution vector pointer
//------------------------------------------------------------------
  Vector *sol = out;

// Set the node checkerboard size of the fermion field
//------------------------------------------------------------------

  size_t f_size_cb;

  if(lat.Fclass() == F_CLASS_CLOVER) {
    f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  } else {
    f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / (lat.FchkbEvl()+1);
  }
    
  if(GJP.Gparity()) f_size_cb*=2;

  Vector *res = (Vector *) smalloc(f_size_cb * sizeof(Float));
  Vector *dir = (Vector *) smalloc(f_size_cb * sizeof(Float));
  Vector *mmp = (Vector *) smalloc(f_size_cb * sizeof(Float));

// If src_norm_sq is not provided calculate it
//------------------------------------------------------------------
  if(src_norm_sq == 0){
    src_norm_sq = src->NormSqNode(f_size_cb); //CK: in G-parity situation we want the norm^2 of the whole 2-flavour double-wrapped source
    glb_sum(&src_norm_sq);
  }
  VRB.Flow(cname,fname,"src_norm_sq=%e\n",src_norm_sq);

// Calculate stopping condition
//------------------------------------------------------------------
  stp_cnd = src_norm_sq * cg_arg.stop_rsd * cg_arg.stop_rsd;
  VRB.Flow(cname,fname, 
	   "stp_cnd =%e\n", IFloat(stp_cnd));

// Make IFloat pointers out of Vector pointers
//------------------------------------------------------------------
  IFloat *f_sol = (IFloat *) sol; 
  IFloat *f_dir = (IFloat *) dir; 
  IFloat *f_res = (IFloat *) res; 
  IFloat *f_mmp = (IFloat *) mmp; 

//------------------------------------------------------------------
// Initial step:
// res = src - MatPcDagMatPc * sol
// dir = res
// if( |res|^2 <= stp_cnd ){ 
//   n_count = 0
//   free memory
//   return
// }
//------------------------------------------------------------------
  Float *in_f =  (Float *) sol;
  // Mmp = MatPcDagMatPc * sol
  dop.MatPcDagMatPc(mmp, sol);
  if (shift){
    mmp -> FTimesV1PlusV2(*shift,sol,mmp, f_size_cb);
  }
  //print_vec( mmp, "mmp");

  // res = src
  res->CopyVec(src, f_size_cb);
  //print_vec( res, "res");

  // res -= mmp
  res->VecMinusEquVec(mmp, f_size_cb);
  //print_vec( res, "res");

  // dir = res
  dir->CopyVec(res, f_size_cb);  
  //print_vec( dir, "dir");

  // res_norm_sq_cur = res * res
  res_norm_sq_cur = res->NormSqNode(f_size_cb);
  //printf("res_norm_sq_cur=%e\n",res_norm_sq_cur);
  glb_sum(&res_norm_sq_cur);

  // if( |res|^2 <= stp_cnd ) we are done
  VRB.Flow(cname,fname,
  	   "|res[0]|^2 = %e\n", IFloat(res_norm_sq_cur));
  itr = 0;
  max_itr = 9999;
  if(res_norm_sq_cur <= stp_cnd) max_itr = 0;
  //printf("max_itr=%d\n",max_itr);


//------------------------------------------------------------------
// Loop over CG iterations
//------------------------------------------------------------------

  for(i=0; i < max_itr; i++){
    itr = itr + 1;
    res_norm_sq_prv = res_norm_sq_cur;

    // mmp = MatPcDagMatPc * dir
    // d = <dir, MatPcDagMatPc*dir>

    dop.MatPcDagMatPc(mmp, dir, &d);
    if (shift){
      mmp -> FTimesV1PlusV2(*shift,dir,mmp, f_size_cb);
      Float dir_sq = dir -> NormSqNode(f_size_cb);
      glb_sum(&dir_sq);
      d += (*shift) * dir_sq;
    }
    //printf("d=%e\n",d);
    //print_vec( mmp, "mmp");
  
    glb_sum(&d);
    VRB.Flow(cname,fname, "d = %e\n", IFloat(d));

    // If d = 0 we are done
    if(d == 0.0) {
      VRB.Warn(cname,fname,"d(%e) = 0.0!!\n",d);
      //	exit(5);
      break;
      //??? or should we give a warning or error? Yes we should, really.
    }

    a = res_norm_sq_prv / d;
    VRB.Flow(cname,fname, "a = %e\n", IFloat(a));

    // Set circular buffer
    //    setCbufCntrlReg(4, CBUF_MODE4);

    // sol = a * dir + sol;
    sol->FTimesV1PlusV2(a, dir, sol, f_size_cb);
    //print_vec( sol, "sol");

    // res = - a * (MatPcDagMatPc * dir) + res;
    res->FTimesV1PlusV2(-a, mmp, res, f_size_cb);
    //print_vec( res, "res");

    // res_norm_sq_cur = res * res
    res_norm_sq_cur = res->NormSqNode(f_size_cb);
    glb_sum(&res_norm_sq_cur);

    // if( |res|^2 <= stp_cnd ) we are done
    VRB.Flow(cname,fname,
	     "|res[%d]|^2 = %e\n", itr, IFloat(res_norm_sq_cur));
    if(res_norm_sq_cur <= stp_cnd) break;

    b = res_norm_sq_cur / res_norm_sq_prv;
    VRB.Flow(cname,fname, "b = %e\n", IFloat(b));

    // dir = b * dir + res;
    dir->FTimesV1PlusV2(b, dir, res, f_size_cb);
    //print_vec( dir, "dir");
  }

  // It has not reached stp_cnd: Issue a warning
  if(itr == cg_arg.max_num_iter - 1){
    VRB.Warn(cname,fname,
	      "CG reached max iterations = %d. |res|^2 = %e\n",
	     itr+1, IFloat(res_norm_sq_cur) );
  }

//------------------------------------------------------------------
// Done. Finish up and return
//------------------------------------------------------------------
  // Calculate and set true residual: 
  // true_res = |src - MatPcDagMatPc * sol| / |src|
  dop.MatPcDagMatPc(mmp, sol);
  if (shift){
    mmp -> FTimesV1PlusV2(*shift,sol,mmp, f_size_cb);
  }
  res->CopyVec(src, f_size_cb);
  res->VecMinusEquVec(mmp, f_size_cb);
  res_norm_sq_cur = res->NormSqNode(f_size_cb);
  glb_sum(&res_norm_sq_cur);
  Float tmp = res_norm_sq_cur / src_norm_sq;
  tmp = sqrt(tmp);
  if(true_res != 0){
    *true_res = tmp;
  }
  VRB.Result(cname,fname,
	     "True |res| / |src| = %e, iter = %d\n", IFloat(tmp), itr+1);

  // Free memory
  VRB.Sfree(cname,fname, "mmp", mmp);
  sfree(mmp);
  VRB.Sfree(cname,fname, "dir", dir);
  sfree(dir);
  VRB.Debug("b ============\n");
  VRB.Sfree(cname,fname, "res", res);
  sfree(res);

  VRB.Debug("a ============\n");

  // Return number of iterations
  return itr+1;
}

int InvCg_CPS(Vector *out, 
	      Vector *in, 
	      Float src_norm_sq, 
	      Float *true_res,
	      DiracOp &dop,
	      Lattice& lat, CgArg &cg_arg){
  return InvCgShift_CPS(out,in,src_norm_sq,true_res,NULL,dop,lat,cg_arg);
}














int InvCGtest(){
  //Test CPS vs. BFM invCG method on twisted mass Dirac operator
  GnoneFwilsonTm* lattice = new GnoneFwilsonTm;

  CgArg cg_arg;
  cg_arg.mass = -1.8;
  cg_arg.epsilon = 0.5;
  cg_arg.max_num_iter = 10000;
  cg_arg.stop_rsd = 1e-08;
  cg_arg.true_rsd = 1e-08;
  cg_arg.Inverter = CG;
  
  DiracOpWilsonTm dop(*lattice, (Vector*)0, (Vector*)0, &cg_arg, CNV_FRM_NO);

  //technically we need a WILSON ordered fermion, but as it is random it doesn't matter
  Float* in = rand_4d_canonical_fermion(*lattice);
  
  size_t f_size_cb =  GJP.VolNodeSites() * 24/2;
  Float *out_1 = (Float*)pmalloc(f_size_cb*sizeof(Float));
  Float *out_2 = (Float*)pmalloc(f_size_cb*sizeof(Float));
  for(int i=0;i<f_size_cb;i++){
    out_1[i] = 0.0; out_2[i] = 0.0;
  }

  Float true_rsd;

  dop.InvCg((Vector*)out_1, (Vector*)in, 0.0, &true_rsd); //bfm version
  
  InvCg_CPS((Vector*)out_2, (Vector*)in, 0.0, &true_rsd, dop, *lattice, cg_arg); 
  
  bool fail(false);
  for(int i=0;i<f_size_cb;i++){
    if( fabs(out_1[i]-out_2[i])>1e-08 ){
      printf("InvCGtest fail %d: %f %f, ratio %f\n",i,out_1[i],out_2[i],out_1[i]/out_2[i]);
      fail=true;
    }
  }
  if(fail){
    printf("Failed InvCg test\n"); exit(-1);
  }else printf("Passed InvCg test\n");
    
  pfree(in);
  pfree(out_1);
  pfree(out_2);
  delete lattice;

  return 0;
}





int main(int argc, char *argv[])
{
    const char *fname = "main()";

    setup(argc, argv);

    //VRB.ElapsedTime("CPS", fname);

    // OpenMP test
    VRB.Result( "CPS", "main", "omp_get_num_threads[1] -> %d", omp_get_num_threads() );
    
    #pragma omp parallel
    {
      if ( UniqueID() == 0 && omp_get_thread_num() == 0 ) {
        VRB.Result( "CPS", "main", "omp_get_num_threads[2] -> %d", omp_get_num_threads() );
      }
    }

    VRB.Result( "CPS", "main", "omp_get_num_threads[3] -> %d", omp_get_num_threads() );

    return InvCGtest();

    
    //////////////////////////////////////////////////////////////////////
    // creating actions and integrators
    AlgMomentum mom;
    AlgActionGauge gauge(mom, gauge_arg);
    AlgActionQuotient quo(mom, quo_arg);

    //!< Construct numerical integrators
    AlgIntAB &ab1 = AlgIntAB::Create(mom, gauge,       ab1_arg);
    AlgIntAB &ab2 = AlgIntAB::Create(ab1, quo,         ab2_arg);
    //////////////////////////////////////////////////////////////////////

    int traj = evo_arg.traj_start;
    for(int conf = 0; conf< evo_arg.gauge_configurations; ++conf) {
        CommonArg common_arg_plaq;
        CommonArg common_arg_pbp;
        CommonArg common_arg_hmc;

        truncate_it(&common_arg_plaq , evo_arg.plaquette_stem      , traj);
        truncate_it(&common_arg_pbp  , evo_arg.pbp_stem            , traj);
        truncate_it(&common_arg_hmc  , evo_arg.evo_stem            , traj);

        // Inner trajectory loop
        for(int i = 0; i < evo_arg.gauge_unload_period; ++i, ++traj) {

            measure_plaq(common_arg_plaq);
            measure_pbp(common_arg_pbp, traj);

            //VRB.ElapsedTime("CPS", "main[2]");
            VRB.Result( "CPS", "main", "omp_get_num_threads[4] -> %d", omp_get_num_threads() );
            run_hmc(common_arg_hmc, traj, ab2);
        }//End of inter-cfg sweep

        checkpoint(traj);

    } //End config loop

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
  return;

    const char *fname = "measure_pbp()";

    // fix pbp_arg
    pbp_arg.src_u_s = 0;
    pbp_arg.src_l_s = 0;
    pbp_arg.snk_u_s = 0;
    pbp_arg.snk_l_s = 0;

    const int g_int = evo_arg.gauge_unload_period;
    if (traj % g_int == 0 && evo_arg.measure_pbp) {
        Float dtime = -dclock();

        LRGState rng_state;
        rng_state.GetStates();

        Lattice &lat = LatticeFactory::Create(F_CLASS_WILSON_TM, G_CLASS_NONE);
        VRB.Result( "cps", "measure_pbp", "LatticeFactory::Create(F_CLASS_WILSON_TM, G_CLASS_NONE)" );

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
