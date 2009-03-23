//
//  gauge + fermion only
//
#include<config.h>

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#include<util/lattice.h>
#include<util/random.h>
#include<util/time_cps.h>

#include<alg/alg_hmc.h>
#include<alg/common_arg.h>
#include<alg/hmc_arg.h>
#include<alg/hmd_arg.h>

#include<alg/alg_int.h>
#include<alg/int_arg.h>

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

#undef USE_SCU_CHECKSUMS
#ifdef USE_SCU_CHECKSUMS
#include <qcdocos/scu_checksum.h>
#endif
//-------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

const int MAX_FILENAME = 256;

HmcArg hmc_arg;
HmcArg hmc_arg_pass;

ActionGaugeArg gauge_arg;
ActionBosonArg bsn_arg;
ActionFermionArg frm_arg;
ActionQuotientArg quo_arg;
ActionQuotientArg quo_tm_arg;
ActionRationalQuotientArg rat_quo_arg;

// define all integrators; vml file arguments will not change
IntABArg ab1_arg;
IntABArg ab2_arg;
IntABArg ab3_arg;
IntABArg ab4_arg;
IntABArg sum_arg;

EvoArg evo_arg;
DoArg do_arg;
PbpArg pbp_arg;
NoArg no_arg;
FloatArray w_spect_mass_list;

void checkpoint(int traj);

int main(int argc, char *argv[]){

  Start(&argc,&argv);

  char plaq_file[256];
  char top_file[256];
  char pbp_file[256];
  char hmc_file[256];

  char *cname=argv[0];
  char *fname="main()";
  Float dtime;

  CommonArg common_arg_plaq;
  CommonArg common_arg_pbp;
  CommonArg common_arg_hmc;
  CommonArg common_arg_top;

  WspectArg w_spect_arg;
  CgArg w_spect_cg_arg;
  WspectOutput w_spect_output;
  FixGaugeArg w_spect_fg_arg;

  if ( argc!=19 ) { 
	if(!UniqueID()){
    printf("Args:\t do_arg.vml hmc_arg.vml evo_arg.vml \n");
    printf("\t gauge_arg.vml bsn_arg.vml frm_arg.vml \n");
    printf("\t quo_tm_arg.vml quo_arg.vml rat_quo_arg.vml \n");
    printf("\t ab1_arg.vml ab2_arg.vml ab3_arg.vml ab4_arg.vml \n");
    printf("\t w_spect_arg.vml w_spect_cg_arg.vml w_spect_mass_list.vml \n");
    printf("\t pbp_arg.vml current_dir \n");
	}
    exit(-1);
  }

  chdir (argv[18]);
// these from Michael (?) and I added one "/"
//gauge_Unload
//gauge_unload



  if ( !do_arg.Decode(argv[1],"do_arg") ) { 
    do_arg.Encode("bum_arg","bum_arg"); printf("Bum do_arg\n"); exit(-1); }
  do_arg.Encode("do_arg.dat", "do_arg");
  if ( !hmc_arg.Decode(argv[2],"hmc_arg")){printf("Bum hmc_arg\n"); exit(-1);}
  hmc_arg.Encode("hmc_arg.dat","hmc_arg");
  if ( !evo_arg.Decode(argv[3],"evo_arg")){printf("Bum evo_arg\n"); exit(-1);}
  evo_arg.Encode("evo_arg.dat","evo_arg");

  if ( !gauge_arg.Decode(argv[4],"gauge_arg")){printf("Bum gauge_arg\n"); exit(-1);}
  gauge_arg.Encode("gauge_arg.dat","gauge_arg");
  if ( !bsn_arg.Decode(argv[5],"bsn_arg")){printf("Bum bsn_arg\n"); exit(-1);}
  bsn_arg.Encode("bsn_arg.dat","bsn_arg");
  if ( !frm_arg.Decode(argv[6],"frm_arg")){printf("Bum frm_arg\n"); exit(-1);}
  frm_arg.Encode("frm_arg.dat","frm_arg");

  if ( !quo_tm_arg.Decode(argv[7],"quo_tm_arg")){printf("Bum quo_tm_arg\n"); exit(-1);}
  quo_tm_arg.Encode("quo_tm_arg.dat","quo_tm_arg");
  if ( !quo_arg.Decode(argv[8],"quo_arg")){printf("Bum quo_arg\n"); exit(-1);}
  quo_arg.Encode("quo_arg.dat","quo_arg");
  if ( !rat_quo_arg.Decode(argv[9],"rat_quo_arg")){printf("Bum rat_quo_arg\n"); exit(-1);}
  rat_quo_arg.Encode("rat_quo_arg.dat","rat_quo_arg");

  if ( !ab1_arg.Decode(argv[10],"ab1_arg")){printf("Bum ab1_arg\n"); exit(-1);}
  ab1_arg.Encode("ab1_arg.dat","ab1_arg");
  if ( !ab2_arg.Decode(argv[11],"ab2_arg")){printf("Bum ab2_arg\n"); exit(-1);}
  ab2_arg.Encode("ab2_arg.dat","ab2_arg");
  if ( !ab3_arg.Decode(argv[12],"ab3_arg")){printf("Bum ab3_arg\n"); exit(-1);}
  ab3_arg.Encode("ab3_arg.dat","ab3_arg");
  if ( !ab4_arg.Decode(argv[13],"ab4_arg")){printf("Bum ab4_arg\n"); exit(-1);}
  ab4_arg.Encode("ab4_arg.dat","ab4_arg");

  if ( !pbp_arg.Decode(argv[14],"pbp_arg")){printf("Bum pbp_arg\n"); exit(-1);}
  pbp_arg.Encode("pbp_arg.dat","pbp_arg");
  if ( !w_spect_arg.Decode(argv[15],"w_spect_arg") )
    { printf("Bum w_spect_arg\n"); exit(-1);}   
  if ( !w_spect_cg_arg.Decode(argv[16],"w_spect_cg_arg") )
    { printf("Bum w_spect_cg_arg\n"); exit(-1);}
  if ( !w_spect_mass_list.Decode(argv[17],"w_spect_mass_list") )
    { printf("Bum w_spect_mass_list\n"); exit(-1);}

//  ApeSmearArg ape_arg("ape_arg.vml");
  ApeSmearArg ape_arg;
  if ( !ape_arg.Decode("ape_arg.vml","ape_arg") )
    { printf("Bum ape_arg\n"); exit(-1);}
  ape_arg.Encode("ape_arg.dat","ape_arg");

  chdir(evo_arg.work_directory);

  // do_arg.verbose_level=VERBOSE_RESULT_LEVEL;
  GJP.Initialize(do_arg);
  //VRB.Level(VERBOSE_FUNC_LEVEL);
  VRB.Level(3);
  LRG.Initialize();

  // Outer config loop
  int traj = evo_arg.traj_start;
  int w_int = evo_arg.measure_w_spect_interval;
  int g_int = evo_arg.gauge_unload_period;

  if( w_int>0 && g_int>0 && w_int%g_int !=0){
    ERR.General("", "main()", "w spec meas int(%d) not divisible by Gauge storing int(%d) \n",w_int,g_int);
  }

  if ( do_arg.start_conf_kind != START_CONF_FILE ) {
    checkpoint(traj);
  }

  hmc_arg_pass = hmc_arg;
  sum_arg.A_steps = 1;
  sum_arg.B_steps = 1;

  VRB.Result(cname,fname,"\n\n\n ******************* start create ACTIONS & INTEGRATORS\n\n\n");
  //!< Create fictitous Hamiltonian (mom + action)

  // must precede all AlgAction allocations
   AlgMomentum mom;
   AlgActionGauge gauge(mom, gauge_arg);
   AlgActionQuotient quotient(mom, quo_arg);
   AlgActionQuotient quotient_tm(mom, quo_tm_arg);
   AlgActionFermion Wfermion(mom, frm_arg);
   AlgActionRationalQuotient rat_quo(mom, rat_quo_arg);

   //!< Construct numerical integrators
   //  embedded integrators:  ab1_arg, ab2_arg, ab3_arg
   //  top-level integrator:  ab4_arg

    AlgIntAB &ab1 = AlgIntAB::Create(mom, gauge,    ab1_arg);
    AlgIntAB &ab2 = AlgIntAB::Create(ab1, quotient_tm,  ab2_arg);
    AlgIntAB &ab3 = AlgIntAB::Create(ab2, rat_quo, ab3_arg);
    AlgIntAB &ab4 = AlgIntAB::Create(ab3, quotient, ab4_arg);


  VRB.Result(cname,fname,"\n\n\n ******************* end create ACTIONS & INTEGRATORS\n\n\n");  

  for(int conf=0; conf< evo_arg.gauge_configurations; conf ++ ) {

    sprintf(plaq_file,"%s.%d",evo_arg.plaquette_stem,traj);
    FILE * truncate_it = Fopen(plaq_file,"w");
    Fclose(truncate_it);
    common_arg_plaq.set_filename(plaq_file);

    sprintf(top_file,"%s.%d","../results/alg_top/top",traj);
    truncate_it = Fopen(top_file,"w");
    Fclose(truncate_it);
    common_arg_top.set_filename(top_file);

    sprintf(pbp_file,"%s.%d",evo_arg.pbp_stem,traj);
	//    printf("pbp_file=%s\n",pbp_file);
    truncate_it = Fopen(pbp_file,"w");
    Fclose(truncate_it);
    common_arg_pbp.set_filename(pbp_file);
    pbp_arg.src_u_s = 0;
    pbp_arg.src_l_s = GJP.SnodeSites() * GJP.Snodes() - 1;
    pbp_arg.snk_u_s = GJP.SnodeSites() * GJP.Snodes() - 1;
    pbp_arg.snk_l_s = 0;

    sprintf(hmc_file,"%s.%d",evo_arg.evo_stem,traj);
    common_arg_hmc.set_filename(hmc_file);

    LRGState rng_state;
    // Inner trajectory loop
    for(int i=0;i<g_int;i++,traj++ ) {
    		
	// Wilson gauge action used for plaquette measurement
	VRB.Result(cname,fname,"\n\n\n ******************* start PLAQUETTE\n\n\n");
	Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_WILSON);
	AlgPlaq plaq(lat, &common_arg_plaq, &no_arg);
  	plaq.run();
	LatticeFactory::Destroy();
	VRB.Result(cname,fname,"\n\n\n ******************* end PLAQUETTE\n\n\n");

	// calculate the topological charge. Need to copy the lattice since
	// we need to smear it first. Use Chulwoo's "lattice container"
	{
	  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_WILSON);
	  LatticeContainer lat_cont;
	  // copy lattice
	  lat_cont.Get(lat);
	  //----- mess up lattice -------------------------
	  //AlgPlaq plaq(lat, &common_arg_top, &no_arg);
	  //plaq.run();
	  AlgApeSmear ape(lat,&common_arg_top,&ape_arg);
	  AlgTcharge  tcharge(lat, &common_arg_top);
	  for (int i(0);i<20;i++)
	    {
	      VRB.Result(cname,fname,"%i\n",i);
	      VRB.Result(cname,fname,"   running tcharge\n"); tcharge.run(); 
	      VRB.Result(cname,fname,"   running ape\n"); ape.run();
	      VRB.Result(cname,fname,"   running ape\n"); ape.run();
	      VRB.Result(cname,fname,"   running ape\n"); ape.run();
	    }
	  tcharge.run();
	  // restore the lattice
	  lat_cont.Set(lat);
	  //plaq.run();
	  LatticeFactory::Destroy();
	}
	if ((traj%g_int == 0) && evo_arg.measure_pbp){
	    dtime = -dclock();
	    VRB.Result(cname,fname,"Running pbp\n");
	    VRB.Result(cname,fname,"\n\n\n ******************* start PBP\n\n\n");
	    rng_state.GetStates();
	    Lattice &lat = LatticeFactory::Create(F_CLASS_DWF, G_CLASS_WILSON);
	    AlgPbp pbp(lat,&common_arg_pbp,&pbp_arg);
	    for(int pbp_counter = 0; pbp_counter < evo_arg.measure_pbp; pbp_counter++) {
	    	pbp.run();
	    }
	    LatticeFactory::Destroy();
	    rng_state.SetStates();
	    VRB.Result(cname,fname,"\n\n\n ******************* end PBP\n\n\n");
	    dtime += dclock();
	    print_flops("AlgPbp","",0,dtime);	
	}
	
	if ( (evo_arg.reproduce_interval > 0) &&
	     (traj % evo_arg.reproduce_interval) == 0 ) { 
		VRB.Result(cname,fname,"Running traj %d with reproduction\n",traj);
		hmc_arg_pass.reproduce = REPRODUCE_YES;
	} else { 
		VRB.Result(cname,fname,"Running traj %d without reproduction\n",traj);
		hmc_arg_pass.reproduce = REPRODUCE_NO;	
	}

	//!< Run hybrid Monte Carlo
#if 0
        VRB.Result(cname,fname,"\n\n\n ******************* start HMC\n\n\n");
        if ( (traj > 19) && (hmc_arg_pass.metropolis == METROPOLIS_NO)) { 
           VRB.Result(cname,fname,"Switching to ACCECPT/REJECT at traj %d \n",traj);
           hmc_arg_pass.metropolis = METROPOLIS_YES;
        }
#endif
	VRB.Result(cname,fname,"\n\n\n ******************* start HMC\n\n\n");
       	AlgHmc hmc(ab4, common_arg_hmc, hmc_arg_pass);
	Float time = -dclock();
	hmc.run();
	VRB.Result(cname,fname,"\n\n\n ******************* end HMC\n\n\n");
	time += dclock();
	print_flops("AlgHmc","run()",0,time);
	    
    }//End of inter-cfg sweep

    //if (evo_arg.measure_w_spect_interval){
    if ((w_int>0) && (traj%w_int==0)){

      dtime = -dclock();
      VRB.Result(cname,fname,"Measuring hadron masses\n");
      char filenames[32][MAX_FILENAME];
      char suffix[MAX_FILENAME];
      CommonArg common_arg;
      CommonArg common_arg_fg;
      char *w_dir = evo_arg.w_spect_directory;
  
      Lattice &lat = LatticeFactory::Create(F_CLASS_DWF, G_CLASS_NONE);
      rng_state.GetStates();

      // some hard coding of VML-lized input parameters.
      w_spect_arg.extended_mesons_on=0;
      w_spect_output.fold = BARYON_RAW;
      w_spect_fg_arg.stop_cond = 1e-8;
      w_spect_fg_arg.max_iter_num = 20000;
      w_spect_fg_arg.hyperplane_start = 0;
      w_spect_fg_arg.hyperplane_step = 1;
      int dir = w_spect_arg.prop_dir;
      w_spect_fg_arg.hyperplane_num = GJP.NodeSites(dir) * GJP.Nodes(dir);

      switch (dir) {
        case 0  : w_spect_fg_arg.fix_gauge_kind = FIX_GAUGE_COULOMB_X ; break ;
        case 1  : w_spect_fg_arg.fix_gauge_kind = FIX_GAUGE_COULOMB_Y ; break ;
        case 2  : w_spect_fg_arg.fix_gauge_kind = FIX_GAUGE_COULOMB_Z ; break ;
        case 3  : w_spect_fg_arg.fix_gauge_kind = FIX_GAUGE_COULOMB_T ; break ;
        default : w_spect_fg_arg.fix_gauge_kind = FIX_GAUGE_NONE      ;
      }
  
  
      int n_mf = w_spect_mass_list.Floats.Floats_len;
      Float *mf = w_spect_mass_list.Floats.Floats_val;
  
      for (int mf_ind=0;mf_ind <n_mf; mf_ind++){
        
        sprintf(suffix,"_mf%0.3f",mf[mf_ind]);
        sprintf(filenames[0],"%s/scalar%s.dat.%d" ,w_dir,suffix,traj);
        sprintf(filenames[1],"%s/vector_x%s.dat.%d" ,w_dir,suffix,traj);
        sprintf(filenames[2],"%s/vector_y%s.dat.%d" ,w_dir,suffix,traj);
        sprintf(filenames[3],"%s/tensor_xy%s.dat.%d",w_dir,suffix,traj);
        sprintf(filenames[4], "%s/vector_z%s.dat.%d",w_dir,suffix,traj);
        sprintf(filenames[5], "%s/tensor_xz%s.dat.%d",w_dir,suffix,traj);
        sprintf(filenames[6], "%s/tensor_yz%s.dat.%d",w_dir,suffix,traj);
        sprintf(filenames[7], "%s/axial_vector_t%s.dat.%d",w_dir,suffix,traj);
        sprintf(filenames[8], "%s/vector_t%s.dat.%d",w_dir,suffix,traj);
        sprintf(filenames[9], "%s/tensor_xt%s.dat.%d",w_dir,suffix,traj);
        sprintf(filenames[10], "%s/tensor_yt%s.dat.%d",w_dir,suffix,traj);
        sprintf(filenames[11], "%s/axial_vector_z%s.dat.%d" ,w_dir,suffix,traj);
        sprintf(filenames[12], "%s/tensor_zt%s.dat.%d" ,w_dir,suffix,traj);
        sprintf(filenames[13], "%s/axial_vector_y%s.dat.%d" ,w_dir,suffix,traj);
        sprintf(filenames[14], "%s/axial_vector_x%s.dat.%d" ,w_dir,suffix,traj);
        sprintf(filenames[15], "%s/pseudoscalar%s.dat.%d",w_dir,suffix,traj);
  
        w_spect_output.meson_name00            = &(filenames[0][0]) ;
        w_spect_output.meson_name01            = &(filenames[1][0]) ;
        w_spect_output.meson_name02            = &(filenames[2][0]) ;
        w_spect_output.meson_name03            = &(filenames[3][0]) ;
        w_spect_output.meson_name04            = &(filenames[4][0]) ;
        w_spect_output.meson_name05            = &(filenames[5][0]) ;
        w_spect_output.meson_name06            = &(filenames[6][0]) ;
        w_spect_output.meson_name07            = &(filenames[7][0]) ;
        w_spect_output.meson_name08            = &(filenames[8][0]) ;
        w_spect_output.meson_name09            = &(filenames[9][0]) ;
        w_spect_output.meson_name10            = &(filenames[10][0]) ;
        w_spect_output.meson_name11            = &(filenames[11][0]) ;
        w_spect_output.meson_name12            = &(filenames[12][0]) ;
        w_spect_output.meson_name13            = &(filenames[13][0]) ;
        w_spect_output.meson_name14            = &(filenames[14][0]) ;
        w_spect_output.meson_name15            = &(filenames[15][0]) ;
    
        int ind_end=16;
        sprintf(filenames[ind_end++],"%s/cg%s.dat.%d" ,w_dir,suffix,traj);
	//        sprintf(filenames[ind_end++],"%s/cg2%s.dat.%d" ,w_dir,suffix,traj);
	//        sprintf(filenames[ind_end++],"%s/pbp%s.dat.%d" ,w_dir,suffix,traj);
        sprintf(filenames[ind_end++],"%s/mid_point%s.dat.%d",w_dir,suffix,traj);
        sprintf(filenames[ind_end++], "%s/a0%s.dat.%d",w_dir,suffix,traj);
	//        sprintf(filenames[ind_end++], "%s/a1%s.dat.%d",w_dir,suffix,traj);
	//        sprintf(filenames[ind_end++], "%s/b1%s.dat.%d",w_dir,suffix,traj);
        sprintf(filenames[ind_end++], "%s/nucleon%s.dat.%d",w_dir,suffix,traj);
        sprintf(filenames[ind_end++], "%s/nucleon_prime%s.dat.%d",w_dir,suffix,traj);
        sprintf(filenames[ind_end++], "%s/delta_x%s.dat.%d",w_dir,suffix,traj);
        sprintf(filenames[ind_end++], "%s/delta_y%s.dat.%d",w_dir,suffix,traj);
        sprintf(filenames[ind_end++], "%s/delta_z%s.dat.%d",w_dir,suffix,traj);
        sprintf(filenames[ind_end++], "%s/delta_t%s.dat.%d",w_dir,suffix,traj);
	//        sprintf(filenames[ind_end++], "%s/rho%s.dat.%d",w_dir,suffix,traj);
	//        sprintf(filenames[ind_end++], "%s/rho_prime%s.dat.%d",w_dir,suffix,traj);
	
        ind_end=16;
        //-------------------------------------------------
        //cg,nucleon and other stuff
        //-------------------------------------------------
        w_spect_output.cg            = &(filenames[ind_end++][0]) ;
	//        w_spect_output.cg2           = &(filenames[ind_end++][0]) ;
	//        w_spect_output.pbp           = &(filenames[ind_end++][0]) ;
        w_spect_output.mid_point     = &(filenames[ind_end++][0]) ;
        w_spect_output.a0_p          = &(filenames[ind_end++][0]) ;
	//        w_spect_output.a1            = &(filenames[ind_end++][0]) ;
	//        w_spect_output.b1            = &(filenames[ind_end++][0]) ;
        w_spect_output.nucleon       = &(filenames[ind_end++][0]) ;
        w_spect_output.nucleon_prime = &(filenames[ind_end++][0]) ;
        w_spect_output.delta_x       = &(filenames[ind_end++][0]) ;
        w_spect_output.delta_y       = &(filenames[ind_end++][0]) ;
        w_spect_output.delta_z       = &(filenames[ind_end++][0]) ;
        w_spect_output.delta_t       = &(filenames[ind_end++][0]) ;
	//        w_spect_output.rho           = &(filenames[ind_end++][0]) ;
	//        w_spect_output.rho_prime     = &(filenames[ind_end++][0]) ;
	
        for(int ind=0;ind<ind_end;ind++){
          filenames[ind][MAX_FILENAME-1]='\0';
          truncate_it = Fopen(filenames[ind],"w");
          Fclose(truncate_it);
        }
	
        common_arg.results = (void *) &w_spect_output ;
        w_spect_cg_arg.mass = mf[mf_ind] ;
	
        if ( ( w_spect_arg.source_kind == BOX_W ) &&	  \
	     ( w_spect_arg.src_box_e[w_spect_arg.prop_dir] !=	\
	       w_spect_arg.src_box_b[w_spect_arg.prop_dir] ) )
          ERR.General(cname, fname, "Invalid box source.\n") ;
	
        if ( (w_spect_arg.source_kind == BOX_W)		\
	     || (w_spect_arg.source_kind == WALL_W) ) {
          AlgFixGauge fg(lat, &common_arg_fg, &w_spect_fg_arg) ;
          fg.run() ;
        }
        AlgWspect ws(lat,&common_arg,&w_spect_arg,&w_spect_cg_arg);
	
        ws.SetCounter(traj, 0) ;
        
        ws.run();
        
        //------------------------------------------------
        // free gauge fixing matrices, if needed.
        //------------------------------------------------
        
        if ( (w_spect_arg.source_kind == BOX_W)		\
	     || (w_spect_arg.source_kind == WALL_W) ) {
          AlgFixGauge fg(lat, &common_arg_fg, &w_spect_fg_arg) ;
          fg.free() ;
        }
      }
      dtime += dclock();
      print_flops("AlgWspect","",0,dtime);
      rng_state.SetStates();
      LatticeFactory::Destroy();
    }// end w spect meas

    checkpoint(traj);
    
  } //End config loop
  
  AlgIntAB::Destroy(ab4);
  AlgIntAB::Destroy(ab3);
  AlgIntAB::Destroy(ab2);
  AlgIntAB::Destroy(ab1);
  
  End();

  return(0);
}

void checkpoint(int traj)
{
  char *cname="cps::";
  char *fname="checkpoint()";

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
    ERR.General(cname,fname,"Failed RNG file %s",rng_file);
  
  // Update the parameter files for restart
  do_arg.start_seed_filename = rng_file;
  do_arg.start_seed_kind = START_SEED_FILE;
  do_arg.start_conf_filename = lat_file;
  do_arg.start_conf_kind = START_CONF_FILE;
  evo_arg.traj_start     = traj;

  char vml_file[256];
  sprintf(vml_file,"do_arg.%d",traj);
  if ( !do_arg.Encode(vml_file,"do_arg") ){
	  printf("bad do_arg encode\n");
	  exit(-1);
  }

  sprintf(vml_file,"hmc_arg.%d",traj);
  if ( !hmc_arg.Encode(vml_file,"hmc_arg") ){
	  printf("bad hmc_arg encode\n");
	  exit(-1);
  }

  sprintf(vml_file,"evo_arg.%d",traj); 
  if ( !evo_arg.Encode(vml_file,"evo_arg")){
	  printf("bad evo_arg encode\n");
	  exit(-1);
  }

  sprintf(vml_file,"quo_arg.%d",traj);
  if ( !quo_arg.Encode(vml_file,"quo_arg") ){
	  printf("bad quo_arg encode\n");
	  exit(-1);
  }

  sprintf(vml_file,"quo_tm_arg.%d",traj);
  if ( !quo_tm_arg.Encode(vml_file,"quo_tm_arg") ){
	  printf("bad quo_tm_arg encode\n");
	  exit(-1);
  }

  sprintf(vml_file,"rat_quo_arg.%d",traj);
  if ( !rat_quo_arg.Encode(vml_file,"rat_quo_arg") ){
	  printf("bad rat_quo_arg encode\n");
	  exit(-1);
  }

  sprintf(vml_file,"gauge_arg.%d",traj);
  if ( !gauge_arg.Encode(vml_file,"gauge_arg") ){
	  printf("bad gauge_arg encode\n");
	  exit(-1);
  }

  sprintf(vml_file,"ab1_arg.%d",traj);
  if ( !ab1_arg.Encode(vml_file,"ab1_arg") ){
	  printf("bad ab1_arg encode\n");
	  exit(-1);
  }

  sprintf(vml_file,"ab2_arg.%d",traj);
  if ( !ab2_arg.Encode(vml_file,"ab2_arg") ){
	  printf("bad ab2_arg encode\n");
	  exit(-1);
  }

  sprintf(vml_file,"ab3_arg.%d",traj);
  if ( !ab3_arg.Encode(vml_file,"ab3_arg") ){
	  printf("bad ab3_arg encode\n");
	  exit(-1);
  }

  sprintf(vml_file,"ab4_arg.%d",traj);
  if ( !ab4_arg.Encode(vml_file,"ab4_arg") ){
	  printf("bad ab4_arg encode\n");
	  exit(-1);
  }

  time += dclock();
  print_flops("","checkpoint()",0,time);

}
