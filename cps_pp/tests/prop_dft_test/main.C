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
#if(0==1)
#include <ReadLattice.h>
#include <WriteLattice.h>
#endif

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

#include <alg/prop_attribute_arg.h>
#include <alg/gparity_contract_arg.h>
#include <alg/propmanager.h>
#include <alg/alg_gparitycontract.h>
#include <alg/prop_dft.h>
#include <util/gparity_singletodouble.h>

#ifdef HAVE_BFM
#include <chroma.h>
#endif

using namespace std;
USING_NAMESPACE_CPS

#define SETUP_ARRAY(OBJ,ARRAYNAME,TYPE,SIZE)		\
  OBJ . ARRAYNAME . ARRAYNAME##_len = SIZE;		\
  OBJ . ARRAYNAME . ARRAYNAME##_val = new TYPE [SIZE]

#define ELEM(OBJ,ARRAYNAME,IDX) OBJ . ARRAYNAME . ARRAYNAME##_val[IDX]


void GaugeTransformU(Matrix *gtrans, Lattice &lat);

bool test_equals(const Rcomplex &a, const Rcomplex &b, const double &eps){
  if( fabs(a.real()-b.real()) > eps || fabs(a.imag()-b.imag()) > eps ) return false;
  return true;
}
void print(const Rcomplex &w){
  printf("(%.4f %.4f) ",w.real(),w.imag());
}

void print(const WilsonMatrix &w){
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      Rcomplex c = w(i,0,j,0);
      printf("(%.4f %.4f) ",c.real(),c.imag());
    }
    printf("\n");
  }
  printf("\n");
}

bool test_equals(const WilsonMatrix &a, const WilsonMatrix &b, const double &eps){
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      for(int aa=0;aa<3;aa++){
	for(int bb=0;bb<3;bb++){
	  Rcomplex ca = a(i,aa,j,bb);
	  Rcomplex cb = b(i,aa,j,bb);
	  if( fabs(ca.real()-cb.real()) > eps || fabs(ca.imag()-cb.imag()) > eps ) return false;
	}
      }
    }
  }
  return true;
}


void print(const Matrix &w){
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      Rcomplex c = w(i,j);
      printf("(%.4f %.4f) ",c.real(),c.imag());
    }
    printf("\n");
  }
  printf("\n");
}

bool test_equals(const Matrix &a, const Matrix &b, const double &eps){
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      const Rcomplex &ca = a(i,j);
      const Rcomplex &cb = b(i,j);
      if( fabs(ca.real()-cb.real()) > eps || fabs(ca.imag()-cb.imag()) > eps ) return false;
    }
  }
  return true;
}


int main(int argc,char *argv[])
{
  Start(&argc,&argv); //initialises QMP

#ifdef HAVE_BFM
  Chroma::initialize(&argc,&argv);
#endif

  CommandLine::is(argc,argv);

  bool save_config(false);
  bool load_config(false);
  bool load_lrg(false);
  bool save_lrg(false);
  char *load_config_file;
  char *save_config_file;
  char *save_lrg_file;
  char *load_lrg_file;
  bool verbose(false);
  bool unit_gauge(false);

  int size[] = {2,2,2,2,2};

  int gfix_type = 0; //0 = Coulomb, 1 Landau

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
    }else if( strncmp(cmd,"-gauge_fix_landau",15) == 0){
      gfix_type = 1;
      i++;   
    }else if( strncmp(cmd,"-gauge_fix_coulomb",15) == 0){
      gfix_type = 0;
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

  if(gfix_type == 0) lattice->FixGaugeAllocate(FIX_GAUGE_COULOMB_T);
  else lattice->FixGaugeAllocate(FIX_GAUGE_LANDAU);

  lattice->FixGauge(1e-06,2000);
  if(!UniqueID()){ printf("Gauge fixing finished\n"); fflush(stdout); }

  //Get gfixmat satisfaction
  {
    IFloat* satisfaction = lattice->GaugeFixCondSatisfaction(lattice->FixGaugePtr(), lattice->FixGaugeKind(), 1e-06);
    if(gfix_type == 0) for(int t=0;t<GJP.TnodeSites();t++) printf("t=%d GFix condition satisfaction %f\n",t,satisfaction[t]);
    else printf("GFix condition satisfaction %f\n",satisfaction[0]);
  }

  //Check new source gauge fixing code
  {
    int c(0),s(0),t(0), flav(0);
    FermionVectorTp f1;
    FermionVectorTp f2;
    f1.ZeroSource();    
    f2.ZeroSource();

    if(gfix_type == 0){
      //Coulomb
      f1.SetWallSource(c, s, t, flav);
      f1.GaugeFixVector(*lattice,s);

      f2.SetWallSource(c, s, t, flav);
      f2.GFWallSource(*lattice, s, 3, t,flav);
    }else{
      //Landau
      f1.SetVolSource(c, s, flav); //use a zero momentum volume source
      f1.GaugeFixVector(*lattice,s);
      
      int p[4] = {0,0,0,0};
      f2.SetLandauGaugeMomentaSource(*lattice, c,s, p, flav);
    }
    bool fail(false);
    for(int f=0; f< 1; f++){
      for(int site = 0; site < GJP.VolNodeSites(); site++){
	for(int spn = 0; spn < 4; spn++){
	  for(int clr = 0; clr < 3 ; clr++){
	    int off = 2*(clr + 3*(spn + 4*(site + GJP.VolNodeSites()*f)));
	    Float *f1s = f1.data()+off;
	    Float *f2s = f2.data()+off;
	    if(fabs(f1s[0]-f2s[0])>1e-12 || fabs(f1s[1]-f2s[1])>1e-12){
	      printf("Err %d %d %d %d: (%f,%f) (%f,%f)\n",clr,spn,site,f,f1s[0],f1s[1],f2s[0],f2s[1]);
	      fail = true;
	    }
	  }
	}
      }
    }
    if(fail){
      printf("Source gauge fixing code check failed\n");
      exit(-1);
    }else printf("Source gauge fixing code check passed\n");
  }

  std::vector<WilsonMatrix> FTprop;
  std::vector<WilsonMatrix> FTprop_nogfix;

  {
    PropManager::clear();
    
    JobPropagatorArgs prop_args2;
    SETUP_ARRAY(prop_args2,props,PropagatorArg,1);
  
    char* names[1] = {"prop"};
    BndCndType bndcnd[1] = {BND_CND_APRD};
    Float masses[1] = {0.1};

    for(int i=0;i<1;i++){
      PropagatorArg &parg = prop_args2.props.props_val[i];
    
      parg.generics.tag = names[i];
      parg.generics.mass = masses[i];
      parg.generics.bc[0] = GJP.Xbc();
      parg.generics.bc[1] = GJP.Ybc();
      //parg.generics.bc[2] = BND_CND_TWISTED; //GJP.Zbc();
      parg.generics.bc[2] = GJP.Zbc();
      parg.generics.bc[3] = bndcnd[i];

      SETUP_ARRAY(parg,attributes,AttributeContainer,3);
    
      //ELEM(parg,attributes,0).type = VOLUME_SOURCE_ATTR;
      ELEM(parg,attributes,0).type = WALL_SOURCE_ATTR;
      WallSourceAttrArg &wallsrcattr = ELEM(parg,attributes,0).AttributeContainer_u.wall_source_attr;
      wallsrcattr.t = 0;

      ELEM(parg,attributes,1).type = CG_ATTR;
      CGAttrArg &cgattr = ELEM(parg,attributes,1).AttributeContainer_u.cg_attr;
      cgattr.max_num_iter = 5000;
      cgattr.stop_rsd = 1e-08;
      cgattr.true_rsd = 1e-08;

      // ELEM(parg,attributes,2).type = MOMENTUM_ATTR;
      // MomentumAttrArg &momarg = ELEM(parg,attributes,2).AttributeContainer_u.momentum_attr;
      // for(int ii=0;ii<3;ii++)
      // 	if(ii==0) momarg.p[ii] = 1;
      // 	else momarg.p[ii] = 0;

      ELEM(parg,attributes,2).type = GAUGE_FIX_ATTR;
      GaugeFixAttrArg &gfarg = ELEM(parg,attributes,2).AttributeContainer_u.gauge_fix_attr;
      gfarg.gauge_fix_src = 1;
      gfarg.gauge_fix_snk = 0;

      // ELEM(parg,attributes,4).type = TWISTED_BC_ATTR;
      // TwistedBcAttrArg &tbcarg = ELEM(parg,attributes,4).AttributeContainer_u.twisted_bc_attr;
      // tbcarg.theta[0] = 0.0; tbcarg.theta[1] = 0.0; tbcarg.theta[2] = 0.3; //units of pi
    }
    if(UniqueID()==0) printf("prop_args contains %d propagators\n", prop_args2.props.props_len);

    PropManager::setup(prop_args2);   
    PropManager::calcProps(*lattice);
    printf("Props inverted\n"); fflush(stdout);
    
    int desired_mom_x[3] = {0,0,0};
    // for(int i=0;i<3;i++) 
    //   if(i==1) desired_mom_x[i] = 2;
    //   else desired_mom_x[i] = 0;
    
    const static Float Pi_const = 3.1415926535897932384626433832795;
    std::vector<Float> momvec(3); 
    for(int i=0;i<3;i++){
      momvec[i] = desired_mom_x[i]* 2*Pi_const/(GJP.Nodes(i)*GJP.NodeSites(i));
    }
    {
      FourierProp<WilsonMatrix> ft;
      ft.add_momentum(momvec);
      FTprop = ft.getFTProp(*lattice,momvec,"prop");
      //write to disk
      ft.write("prop","ftprop.dat",*lattice);
    }
    {
      FourierProp<WilsonMatrix> ft;
      ft.gaugeFixSink(false);
      ft.add_momentum(momvec);
      FTprop_nogfix = ft.getFTProp(*lattice,momvec,"prop");
    }

    //calculate a few bilinears and quadrilinears and test them against dumb versions

    PropagatorBilinear<WilsonMatrix> bilinear;
    bilinear.add_momentum(momvec);

    ContractedBilinear<WilsonMatrix> conbil;
    conbil.add_momentum(momvec);

    {
      //try all combinations
      ContractedBilinear<WilsonMatrix> conbil2;
      conbil2.add_momentum(momvec);
      conbil2.write("prop", PropDFT::Dagger,
		   "prop", PropDFT::None,
		    "allconbil.dat",*lattice);
    }



    ContractedWallSinkBilinear<WilsonMatrix> conwallsinkbil;
    conwallsinkbil.add_momentum(momvec);


    std::vector<WilsonMatrix> wsprop1, wsprop2;
    {
      std::vector<Float> halfmomvec(momvec); for(int ii=0;ii<3;ii++) halfmomvec[ii]/=2.0;
      std::vector<Float> mhalfmomvec(halfmomvec); for(int ii=0;ii<3;ii++) mhalfmomvec[ii] *= -1;
      FourierProp<WilsonMatrix> ft;
      ft.gaugeFixSink(true);
      ft.add_momentum(halfmomvec);
      ft.add_momentum(mhalfmomvec);
      wsprop1 = ft.getFTProp(*lattice,mhalfmomvec,"prop");
      wsprop2 = ft.getFTProp(*lattice,halfmomvec,"prop");
    }


    int ntest = 4;
    int gs[4][2] = { {0,0}, {15,15}, {0,1}, {2,7} };

    for(int test=0;test<ntest;test++){
      //std::vector<std::pair<WilsonMatrix,WilsonMatrix> > quad_result;
      std::vector<WilsonMatrix> bilinear_result;
      std::vector<Rcomplex> conbil_result;
      std::vector<Rcomplex> conwallsinkbil_result;

      //----------
      bilinear_result = bilinear.getBilinear(*lattice,momvec,
					     "prop", PropDFT::Dagger,
					     "prop", PropDFT::None,
					     gs[test][0],0);

      std::ostringstream bilfname; bilfname << "bilinear_test" << test << ".dat";
      
      bilinear.write("prop", PropDFT::Dagger,
		     "prop", PropDFT::None,
		     gs[test][0],0,
		     bilfname.str().c_str(),*lattice);

      //----------
      conbil_result = conbil.getBilinear(*lattice,momvec,
					 "prop", PropDFT::Dagger,
					 "prop", PropDFT::None,
					 gs[test][0],gs[test][1]);

      std::ostringstream conbilfname; conbilfname << "contracted_bilinear_test"<<test << ".dat";

      conbil.write("prop", PropDFT::Dagger,
		   "prop", PropDFT::None,
		   gs[test][0],0,gs[test][1],0,
		   conbilfname.str().c_str(),*lattice);

      //--------
      conwallsinkbil_result = conwallsinkbil.getBilinear(*lattice,momvec,
							 "prop", PropDFT::Dagger,
							 "prop", PropDFT::None,
							 gs[test][0],0,gs[test][1],0);

      std::ostringstream conwallsinkbilfname; conwallsinkbilfname << "contracted_wallsink_bilinear_test"<<test << ".dat";

      conwallsinkbil.write("prop", PropDFT::Dagger,
			   "prop", PropDFT::None,
			   gs[test][0],0,gs[test][1],0,
			   conwallsinkbilfname.str().c_str(),*lattice);



      int global_T = GJP.Tnodes()*GJP.TnodeSites();
      //std::vector<std::pair<WilsonMatrix,WilsonMatrix> > dumb_quad_result(global_T);
      std::vector<WilsonMatrix> dumb_bilinear_result(global_T);
      std::vector<Rcomplex> dumb_conbil_result(global_T);
      std::vector<Rcomplex> dumb_conwallsinkbil_result(global_T);
      std::vector<WilsonMatrix> dumb_FTprop(global_T);

      for(int t=0;t<global_T;t++){ 
	//dumb_quad_result[t].first = 0.0; dumb_quad_result[t].second = 0.0; 
	dumb_bilinear_result[t] = 0.0;
	dumb_conbil_result[t] = 0.0;
	dumb_conwallsinkbil_result[t] = 0.0;
	dumb_FTprop[t] = 0.0;
      }
      

      //wallsink test
      {
	for(int t=0;t<global_T;t++){ 
	  WilsonMatrix mat = wsprop1[t];
	  mat.hconj();
	  AlgGparityContract::qdp_gl(mat,gs[test][0]);      
	  AlgGparityContract::qdp_gr(mat,gs[test][1]);      

	  WilsonMatrix mat2 = wsprop2[t];

	  dumb_conwallsinkbil_result[t] = Trace(mat,mat2);
	}

      }

      PropagatorContainer &prop = PropManager::getProp("prop");

      for(int x=0;x<GJP.VolNodeSites();x++){
	int pos[4];
	int rem = x; 
	pos[0] = rem % GJP.XnodeSites(); rem/=GJP.XnodeSites();
	pos[1] = rem % GJP.YnodeSites(); rem/=GJP.YnodeSites();
	pos[2] = rem % GJP.ZnodeSites(); rem/=GJP.ZnodeSites();
	pos[3] = rem;

	int global_t = pos[3] + GJP.TnodeCoor()*GJP.TnodeSites();
	int glob_sp[3] = { pos[0] + GJP.XnodeCoor()*GJP.XnodeSites(), pos[1] + GJP.YnodeCoor()*GJP.YnodeSites(), pos[2] + GJP.ZnodeCoor()*GJP.ZnodeSites() };

	Float pdotx = 0.0;
	for(int i=0;i<3;i++) pdotx += momvec[i]*glob_sp[i];
	Rcomplex phase(cos(pdotx),sin(pdotx));

	WilsonMatrix bils[2];
	WilsonMatrix tmp = prop.getProp(*lattice).SiteMatrix(x);
	WilsonMatrix tmpdag = tmp; tmpdag.hconj();

	bils[0] = tmpdag;
	AlgGparityContract::qdp_gr(bils[0],gs[test][0]);      
	bils[0] *= tmp;
	bils[0] *= phase;
	
	bils[1] = tmpdag;
	AlgGparityContract::qdp_gr(bils[1],gs[test][1]);      
	bils[1] *= tmp;	


	WilsonMatrix bilbil = tmpdag;
	AlgGparityContract::qdp_gr(bilbil,gs[test][0]);      
	bilbil *= tmp;
	bilbil *= phase;

	dumb_bilinear_result[global_t] += bilbil;

	bilbil = tmpdag;
	AlgGparityContract::qdp_gl(bilbil,gs[test][0]);      
	AlgGparityContract::qdp_gr(bilbil,gs[test][1]);      

	dumb_conbil_result[global_t] += Trace(bilbil,tmp) * phase;

	WilsonMatrix ft = tmp * phase;
	_FourierProp_helper<WilsonMatrix>::mult_gauge_fix_mat(ft,x,*lattice);
	
	dumb_FTprop[global_t] += ft;
      }

      for(int t=0;t<global_T;t++){
	printf("test %d t=%d, dumb conbil: %f %f\n",test,t,dumb_conbil_result[t].real(),dumb_conbil_result[t].imag());
      }


      bool fail = false;
      for(int t=0;t<global_T;t++){
	_FourierProp_helper<WilsonMatrix>::lattice_sum(dumb_bilinear_result[t]);
	_FourierProp_helper<WilsonMatrix>::lattice_sum(dumb_FTprop[t]);
	slice_sum( (Float*)&dumb_conbil_result[t], 2, 99); 
	slice_sum( (Float*)&dumb_conwallsinkbil_result[t], 2, 99); 
	
	if(!test_equals(bilinear_result[t],dumb_bilinear_result[t],1e-08)){
	  printf("t=%d test %d fail bilinear repro\n",t,test);
	  print(bilinear_result[t]);
	  print(dumb_bilinear_result[t]);
	  fail = true;
	}else printf("t=%d test %d pass bilinear repro\n",t,test);

	if(!test_equals(conbil_result[t],dumb_conbil_result[t],1e-08)){
	  printf("t=%d test %d fail conbil repro\n",t,test);
	  fail=true;
	}else printf("t=%d test %d pass conbil repro\n",t,test);

	if(!test_equals(conwallsinkbil_result[t],dumb_conwallsinkbil_result[t],1e-08)){
	  printf("t=%d test %d fail conwallsinkbil repro: library (%f,%f) dumb (%f,%f) \n",t,test,conwallsinkbil_result[t].real(), conwallsinkbil_result[t].imag(), 
		 dumb_conwallsinkbil_result[t].real(), dumb_conwallsinkbil_result[t].imag());
	  fail=true;
	}else printf("t=%d test %d pass conwallsinkbil repro\n",t,test);

	if(!test_equals(FTprop[t],dumb_FTprop[t],1e-08)){
	  printf("t=%d FTprop fail\n",t,test);
	  fail=true;
	}else printf("t=%d FTprop pass\n",t,test);

      }
      if(fail){
	printf("Dumb test fail\n"); exit(-1);
      }else printf("Dumb test pass\n");

      
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





















// void print(const WilsonMatrix &w){
//   for(int i=0;i<4;i++){
//     for(int j=0;j<4;j++){
//       Complex c = w(i,0,j,0);
//       printf("(%.2f %.2f) ",c.real(),c.imag());
//     }
//     printf("\n");
//   }
//   printf("\n");
// }

// bool test_equals(const WilsonMatrix &a, const WilsonMatrix &b, const double &eps){
//   for(int i=0;i<4;i++){
//     for(int j=0;j<4;j++){
//       for(int aa=0;aa<3;aa++){
// 	for(int bb=0;bb<3;bb++){
// 	  Complex ca = a(i,aa,j,bb);
// 	  Complex cb = b(i,aa,j,bb);
// 	  if( fabs(ca.real()-cb.real()) > eps || fabs(ca.imag()-cb.imag()) > eps ) return false;
// 	}
//       }
//     }
//   }
//   return true;
// }





//   if(gparity_X && !gparity_Y){
//     if(GJP.Xnodes()>1){
//       

//       SingleToDoubleLattice lattdoubler(gparity_X,gparity_Y,orig_gfield,double_latt);
//       lattdoubler.Run();

//       if(!UniqueID()){ printf("Finished setting up doubled lattice\n"); fflush(stdout); }
//     }else{
//       //only one node in X-direction
//       //copy data from orig_latt stored on this node
//       int pos[4];
//       for(pos[0]=0;pos[0]<GJP.XnodeSites();pos[0]++){
// 	for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]++){
// 	  for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]++){
// 	    for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]++){
// 	      for(int mu=0;mu<4;mu++){
// 		int dbl_site_idx = 4*(pos[0]+GJP.XnodeSites() *(pos[1]+GJP.YnodeSites()*(pos[2]+GJP.ZnodeSites()*pos[3])))+mu;
// 		if(pos[0] < GJP.XnodeSites()/2){
// 		  //site is stored on-node in orig_gfield
// 		  int orig_site_idx = 4*(pos[0]+GJP.XnodeSites()/2 *(pos[1]+GJP.YnodeSites()*(pos[2]+GJP.ZnodeSites()*pos[3])))+mu;
// 		  dbl_gfield[dbl_site_idx] = orig_gfield[orig_site_idx];
// 		}else{
// 		  //site is stored on-node in orig_gfield but needs complex conjugating
// 		  int orig_site_idx = 4*(pos[0]-GJP.XnodeSites()/2 +GJP.XnodeSites()/2 *(pos[1]+GJP.YnodeSites()*(pos[2]+GJP.ZnodeSites()*pos[3])))+mu;
// 		  dbl_gfield[dbl_site_idx].Conj(orig_gfield[orig_site_idx]);
// 		}
// 	      }
// 	    }
// 	  }
// 	}
//       }
// #if 0
//       printf("Converted single to double lattice. Scan across X-direction of original:\n");
//       for(int mu=0;mu<4;mu++){
// 	for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]++){
// 	  for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]++){
// 	    for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]++){
// 	      for(pos[0]=0;pos[0]<GJP.XnodeSites()/2;pos[0]++){
// 		int orig_site_idx = 4*(pos[0]+GJP.XnodeSites()/2 *(pos[1]+GJP.YnodeSites()*(pos[2]+GJP.ZnodeSites()*pos[3])))+mu;
// 		printf("U(%d,%d,%d,%d)_%d = (%f,%f) ",pos[0],pos[1],pos[2],pos[3],mu,*((IFloat*)(orig_gfield+orig_site_idx)),*((IFloat*)(orig_gfield+orig_site_idx)+1));
// 	      }
// 	      printf("\n");
// 	    }
// 	  }
// 	}
//       }
//       printf("\nDoubled lattice:\n");
//       for(int mu=0;mu<4;mu++){
// 	for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]++){
// 	  for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]++){
// 	    for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]++){
// 	      for(pos[0]=0;pos[0]<GJP.XnodeSites();pos[0]++){
// 		int dbl_site_idx = 4*(pos[0]+GJP.XnodeSites() *(pos[1]+GJP.YnodeSites()*(pos[2]+GJP.ZnodeSites()*pos[3])))+mu;
// 		printf("U(%d,%d,%d,%d)_%d = (%f,%f) ",pos[0],pos[1],pos[2],pos[3],mu,*((IFloat*)(dbl_gfield+dbl_site_idx)),*((IFloat*)(dbl_gfield+dbl_site_idx)+1));
// 	      }
// 	      printf("\n");
// 	    }
// 	  }
// 	}
// 	printf("\n");
//       }
// #endif
//     }


//   }else if(gparity_X && gparity_Y){
//     if(!UniqueID()){ printf("Setting up quad lattice. sizeof(Float) %d sizeof(IFloat) %d\n",sizeof(Float), sizeof(IFloat)); fflush(stdout); }

//       SingleToDoubleLattice lattdoubler(gparity_X,gparity_Y,orig_gfield,double_latt);
//       lattdoubler.Run();

//       if(!UniqueID()){ printf("Finished setting up quad lattice\n"); fflush(stdout); }
//   }

// }












// if(gparity_X && !gparity_Y){
//   if(!UniqueID()) printf("Setting up RNG from original stacked version\n");


// //orig 4D rng 2 stacked 4D volumes
// //orig   ([R_0 R_1][R'_0 R'_1]) ([R_2 R_3][R'_2 R'_3]) ([R_4 R_5][R'_4 R'_5]) ([R_6 R_7][R'_6 R'_7]) ([R_8 R_9][R'_8 R'_9]) ([R_10 R_11][R'_10 R'_11]) ([R_12 R_13][R'_12 R'_13]) ([R_14 R_15][R'_14 R'_15])
// //double (R_0 R_1 R_2 R_3)      (R_4 R_5 R_6 R_7)      (R_8 R_9 R_10 R_11)    (R_12 R_13 R_13 R_15)  (R'_0 R'_1 R'_2 R'_3)  (R'_4 R'_5 R'_6 R'_7)      (R'_8 R'_9 R'_10 R'_11)    (R'_12 R'_13 R'_14 R'_15)

//   if(GJP.Xnodes()>1){
//     SingleToDouble4dRNG fourDsetup;
//     SingleToDouble5dRNG fiveDsetup;

//     LRG.Reinitialize(); //reset the LRG and prepare for doubled lattice form
      
//     if(!UniqueID()){ printf("Setting up 4D RNG\n"); fflush(stdout); }
//     fourDsetup.Run(gparity_X,gparity_Y);      
//     if(!UniqueID()){ printf("Setting up 5D RNG\n"); fflush(stdout); }
//     fiveDsetup.Run(gparity_X,gparity_Y);    
//   }else{    
//     int n_rgen_4d = GJP.VolNodeSites()/16; //applies both to original and doubled latt
//     int n_rgen = n_rgen_4d;
//     if (GJP.SnodeSites()>=2)
// 	n_rgen = GJP.VolNodeSites()*GJP.SnodeSites() / 32;

//     int stk_index_4d_off = n_rgen_4d/2; //offset for R' on 4D orig latt
//     int blocks_per_s_layer = n_rgen /( GJP.SnodeSites() / 2 ); //also same for original and doubled latt
//     int stk_index_5d_off = blocks_per_s_layer/2; //offset for R' on 5D orig latt

//     //copy the originals
//     UGrandomGenerator *ugran_4d_orig = new UGrandomGenerator[n_rgen_4d];
//     for(int i=0;i<n_rgen_4d;i++) ugran_4d_orig[i] = LRG.UGrandGen4D(i);

//     UGrandomGenerator *ugran_orig = new UGrandomGenerator[n_rgen];
//     for(int i=0;i<n_rgen;i++) ugran_orig[i] = LRG.UGrandGen(i);

//     LRG.Reinitialize(); //reset the LRG and prepare for doubled lattice form

  
//     int pos[5];

//     for(pos[4]=0;pos[4]<GJP.SnodeSites();pos[4]+=2){
// 	for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]+=2){
// 	  for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]+=2){
// 	    for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]+=2){
// 	      for(pos[0]=0;pos[0]<GJP.XnodeSites();pos[0]+=2){
// 		//do the 4D RNG
// 		if(pos[4]==0){
// 		  if(pos[0]>=GJP.XnodeSites()/2){
// 		    int orig_idx = (pos[0]-GJP.XnodeSites()/2)/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + stk_index_4d_off;
// 		    int new_idx = pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
// 		    LRG.UGrandGen4D(new_idx) = ugran_4d_orig[orig_idx];
// 		  }else{
// 		    int orig_idx = pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
// 		    int new_idx = pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
// 		    LRG.UGrandGen4D(new_idx) = ugran_4d_orig[orig_idx];
// 		  }
// 		}
// 		//do the 5D RNG
// 		if(pos[0]>=GJP.XnodeSites()/2){
// 		  int orig_idx = pos[4]/2*blocks_per_s_layer + (pos[0]-GJP.XnodeSites()/2)/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + stk_index_5d_off;
// 		  int new_idx = pos[4]/2*blocks_per_s_layer + pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
// 		  LRG.UGrandGen(new_idx) = ugran_orig[orig_idx];
// 		}else{
// 		  int orig_idx = pos[4]/2*blocks_per_s_layer + pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
// 		  int new_idx = pos[4]/2*blocks_per_s_layer + pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
// 		  LRG.UGrandGen(new_idx) = ugran_orig[orig_idx];
// 		}

// 	      }
// 	    }
// 	  }
// 	}
//     }
//     delete[] ugran_4d_orig;
//     delete[] ugran_orig;
//   }//single node

// }//gpx and gpy
// else if(gparity_X && gparity_Y){
//   SingleToDouble4dRNG fourDsetup;
//   SingleToDouble5dRNG fiveDsetup;
    
//   LRG.Reinitialize(); //reset the LRG and prepare for doubled lattice form
    
//   if(!UniqueID()){ printf("Setting up 4D RNG\n"); fflush(stdout); }
//   fourDsetup.Run(gparity_X,gparity_Y);      
//   if(!UniqueID()){ printf("Setting up 5D RNG\n"); fflush(stdout); }
//   fiveDsetup.Run(gparity_X,gparity_Y);  
// }















// #ifdef DOUBLE_RNG_TEST
//   if(gparity_X && !gparity_Y){
//     if(!UniqueID()) printf("Testing RNG\n");

//     //generate rands again and compare to originals
//     if(GJP.Xnodes()==1){
//       int pos[5];

//       for(pos[4]=0;pos[4]<GJP.SnodeSites();pos[4]+=2){
// 	for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]+=2){
// 	  for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]+=2){
// 	    for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]+=2){
// 	      for(pos[0]=0;pos[0]<GJP.XnodeSites();pos[0]+=2){
// 		LRG.AssignGenerator(pos);
// 		if(pos[4]==0){
// 		  int origflav = 0;
// 		  int origidx = pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]+GJP.ZnodeSites()/2*pos[3]));
// 		  if(pos[0]>=GJP.XnodeSites()/2){
// 		    origflav = 1;
// 		    origidx = (pos[0]-GJP.XnodeSites()/2)/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2+GJP.ZnodeSites()/2*pos[3]/2));
// 		  }
// 		  IFloat &origval = orig_4d_rands[origflav][origidx];
// 		  IFloat newval = LRG.Urand(FOUR_D);//do the 4D RNG
// 		  if(origval != newval){ 
// 		    printf("4D RNG disparity: (%d %d %d %d): orig %f new %f\n",pos[0],pos[1],pos[2],pos[3],origval,newval); exit(-1);
// 		  }
// 		}
// 		int origflav = 0;
// 		int origidx = pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2+GJP.ZnodeSites()/2*(pos[3]/2+GJP.SnodeSites()/2*pos[4]/2)));
// 		if(pos[0]>=GJP.XnodeSites()/2){
// 		  origflav = 1;
// 		  origidx = (pos[0]-GJP.XnodeSites()/2)/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2+GJP.ZnodeSites()/2*(pos[3]/2+GJP.SnodeSites()/2*pos[4]/2)));
// 		}
// 		IFloat &origval = orig_5d_rands[origflav][origidx];
// 		IFloat newval = LRG.Urand(FIVE_D);//do the 5D RNG
// 		if(origval != newval){ 
// 		  printf("5D RNG disparity: (%d %d %d %d): orig %f new %f\n",pos[0],pos[1],pos[2],pos[3],origval,newval); exit(-1);
// 		}
// 	      }
// 	    }
// 	  }
// 	}
//       }

//     }else{
//       bool printnode = false;
//       if(GJP.YnodeCoor()==0 && GJP.ZnodeCoor()==0 && GJP.TnodeCoor()==0) printnode=true;

//       { //4D RNG check
// 	int n_rgen_4d = GJP.VolNodeSites()/16;// both flavours (remember volume doubled now)      
// 	int buf_size = n_rgen_4d * sizeof(Float);
// 	Float *recv_buf = (Float *) pmalloc(buf_size);
// 	Float *send_buf = (Float *) pmalloc(buf_size);
    
// 	for(int f=0;f<2;f++){
// 	  int foff = n_rgen_4d/2;
// 	  for(int site =0; site < n_rgen_4d/2; site++){
// 	    send_buf[site+f*foff] = orig_4d_rands[f][site];
// 	    //if(printnode) printf("Node %d: f %d site %d rand %f\n",GJP.XnodeCoor(),f,site,orig_4d_rands[f][site]);
// 	  }
// 	}

// 	// for(int i=0;i<n_rgen_4d;i++){
// 	// 	if(printnode) printf("Node %d: send_buf site %d = %f\n",GJP.XnodeCoor(),i,send_buf[i]);
// 	// }
      
// 	int data_nodecoor_hf1;  //what xnode coor is this nodes data for the first halves currently stored on? (second half is always on the next node)
// 	int data_nodecoor_hf2;
// 	int data_flav = 0; 
// 	int x_origin = GJP.XnodeCoor()*GJP.XnodeSites(); //x position of start of first half
// 	if(GJP.XnodeCoor()>=GJP.Xnodes()/2){  
// 	  x_origin = (GJP.XnodeCoor()-GJP.Xnodes()/2)*GJP.XnodeSites(); 
// 	  data_flav = 1;
// 	}
// 	data_nodecoor_hf1 = (x_origin/(GJP.XnodeSites()/2) ) % GJP.Xnodes();
// 	data_nodecoor_hf2 = (data_nodecoor_hf1+1) % GJP.Xnodes();

// 	Float nodes_unhappy = 1.0;
// 	Float *cur_data_buf = send_buf;
// 	Float *send_buf_p = send_buf;
// 	Float *recv_buf_p = recv_buf;
// 	int xnode_coor_of_buf_data = GJP.XnodeCoor(); 
// 	int got_hf1 = 0;
// 	int got_hf2 = 0;

// 	int pos[5]; pos[4] = 0;
// 	while(nodes_unhappy != 0.0){
// 	  if(xnode_coor_of_buf_data == data_nodecoor_hf1 || xnode_coor_of_buf_data == data_nodecoor_hf2 ){
// 	    if(xnode_coor_of_buf_data == data_nodecoor_hf1 && printnode) printf("Node %d: Buffer contains data from node %d, testing 1st half (flav %d)\n",GJP.XnodeCoor(),xnode_coor_of_buf_data,data_flav);
// 	    else if(xnode_coor_of_buf_data == data_nodecoor_hf2 && printnode) printf("Node %d: Buffer contains data from node %d, testing 2nd half (flav %d)\n",GJP.XnodeCoor(),xnode_coor_of_buf_data,data_flav);

// 	    for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]+=2){
// 	      for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]+=2){
// 		for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]+=2){
// 		  for(pos[0]=0;pos[0]<GJP.XnodeSites();pos[0]+=2){
// 		    LRG.AssignGenerator(pos);

// 		    if(pos[0]>=GJP.XnodeSites()/2 && xnode_coor_of_buf_data == data_nodecoor_hf2){
// 		      int orig_idx = (pos[0]-GJP.XnodeSites()/2)/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + data_flav * n_rgen_4d/2;

// 		      IFloat &origval = cur_data_buf[orig_idx];
// 		      IFloat newval = LRG.Urand(FOUR_D);//do the 4D RNG
// 		      if(origval != newval){ 
// 			printf("Node %d: 4D RNG disparity: (%d %d %d %d): orig %f new %f (idx %d)\n",GJP.XnodeCoor(),pos[0],pos[1],pos[2],pos[3],origval,newval,orig_idx); exit(-1);
// 		      }
// 		      got_hf2 = 1;
// 		    }else if(pos[0]<GJP.XnodeSites()/2 && xnode_coor_of_buf_data == data_nodecoor_hf1){
// 		      int orig_idx = pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + data_flav * n_rgen_4d/2;

// 		      IFloat &origval = cur_data_buf[orig_idx];
// 		      IFloat newval = LRG.Urand(FOUR_D);//do the 4D RNG
// 		      if(origval != newval){ 
// 			printf("Node %d: 4D RNG disparity: (%d %d %d %d): orig %f new %f (idx %d)\n",GJP.XnodeCoor(),pos[0],pos[1],pos[2],pos[3],origval,newval,orig_idx); exit(-1);
// 		      }
// 		      got_hf1 = 1;
// 		    }
// 		  }
// 		}
// 	      }
// 	    }
// 	  }
// 	  if(got_hf1 && got_hf2) nodes_unhappy = 0.0;
// 	  else nodes_unhappy = 1.0;
	
// 	  glb_sum(&nodes_unhappy);
// 	  if(!UniqueID()) printf("nodes_unhappy = %f\n",nodes_unhappy);

// 	  if(nodes_unhappy!=0.0){
// 	    if(!UniqueID()) printf("Passing data left\n");
// 	    cur_data_buf = recv_buf_p;
// 	    getPlusData((IFloat *)recv_buf_p, (IFloat *)send_buf_p, buf_size/sizeof(Float), 0);
// 	    xnode_coor_of_buf_data = (xnode_coor_of_buf_data+1) % GJP.Xnodes();
// 	    //swap buffers over for next send
// 	    Float *tmp = recv_buf_p;
// 	    recv_buf_p = send_buf_p;
// 	    send_buf_p = tmp;

// 	    // for(int i=0;i<n_rgen_4d;i++){
// 	    //   if(printnode) printf("Node %d: post-pass recv_buf site %d = %f\n",GJP.XnodeCoor(),i,cur_data_buf[i]);
// 	    // }
// 	  }


// 	}
// 	pfree(recv_buf);
// 	pfree(send_buf);
//       }//end of 4d RNG check

//       { //5D RNG check
// 	int n_rgen = GJP.SnodeSites()/2*GJP.VolNodeSites()/16;// both flavours (remember volume doubled now)      
// 	int buf_size = n_rgen * sizeof(Float);
// 	Float *recv_buf = (Float *) pmalloc(buf_size);
// 	Float *send_buf = (Float *) pmalloc(buf_size);
	
// 	int pos[5];
// 	for(pos[4]=0;pos[4]<GJP.SnodeSites();pos[4]+=2){
// 	  for(int f=0;f<2;f++){
// 	    for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]+=2){
// 	      for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]+=2){
// 		for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]+=2){
// 		  for(pos[0]=0;pos[0]<GJP.XnodeSites()/2;pos[0]+=2){
// 		    int osite = pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2+GJP.ZnodeSites()/2*(pos[3]/2+GJP.SnodeSites()/2*pos[4]/2)));
// 		    int bsite = pos[4]/2*GJP.VolNodeSites()/16 +  pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2+GJP.ZnodeSites()/2*pos[3]/2)) + f*GJP.VolNodeSites()/32;
// 		    send_buf[bsite] = orig_5d_rands[f][osite];
// 		  }
// 		}
// 	      }
// 	    }
// 	  }
// 	}
      
// 	int data_nodecoor_hf1;  //what xnode coor is this nodes data for the first halves currently stored on? (second half is always on the next node)
// 	int data_nodecoor_hf2;
// 	int data_flav = 0; 
// 	int x_origin = GJP.XnodeCoor()*GJP.XnodeSites(); //x position of start of first half
// 	if(GJP.XnodeCoor()>=GJP.Xnodes()/2){  
// 	  x_origin = (GJP.XnodeCoor()-GJP.Xnodes()/2)*GJP.XnodeSites(); 
// 	  data_flav = 1;
// 	}
// 	data_nodecoor_hf1 = (x_origin/(GJP.XnodeSites()/2) ) % GJP.Xnodes();
// 	data_nodecoor_hf2 = (data_nodecoor_hf1+1) % GJP.Xnodes();

// 	Float nodes_unhappy = 1.0;
// 	Float *cur_data_buf = send_buf;
// 	Float *send_buf_p = send_buf;
// 	Float *recv_buf_p = recv_buf;
// 	int xnode_coor_of_buf_data = GJP.XnodeCoor(); 
// 	int got_hf1 = 0;
// 	int got_hf2 = 0;
	
// 	while(nodes_unhappy != 0.0){
// 	  if(xnode_coor_of_buf_data == data_nodecoor_hf1 || xnode_coor_of_buf_data == data_nodecoor_hf2 ){
// 	    if(xnode_coor_of_buf_data == data_nodecoor_hf1 && printnode) printf("Node %d: Buffer contains data from node %d, testing 1st half (flav %d)\n",GJP.XnodeCoor(),xnode_coor_of_buf_data,data_flav);
// 	    else if(xnode_coor_of_buf_data == data_nodecoor_hf2 && printnode) printf("Node %d: Buffer contains data from node %d, testing 2nd half (flav %d)\n",GJP.XnodeCoor(),xnode_coor_of_buf_data,data_flav);

// 	    for(pos[4]=0;pos[4]<GJP.SnodeSites();pos[4]+=2){
// 	      for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]+=2){
// 		for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]+=2){
// 		  for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]+=2){
// 		    for(pos[0]=0;pos[0]<GJP.XnodeSites();pos[0]+=2){
// 		      LRG.AssignGenerator(pos);

// 		      if(pos[0]>=GJP.XnodeSites()/2 && xnode_coor_of_buf_data == data_nodecoor_hf2){
// 			int orig_idx = pos[4]/2 * GJP.VolNodeSites()/16 + (pos[0]-GJP.XnodeSites()/2)/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + data_flav * GJP.VolNodeSites()/32;

// 			IFloat &origval = cur_data_buf[orig_idx];
// 			IFloat newval = LRG.Urand(FIVE_D);
// 			if(origval != newval){ 
// 			  printf("Node %d: 5D RNG disparity: (%d %d %d %d %d): orig %f new %f (idx %d)\n",GJP.XnodeCoor(),pos[0],pos[1],pos[2],pos[3],pos[4],origval,newval,orig_idx); exit(-1);
// 			}
// 			got_hf2 = 1;
// 		      }else if(pos[0]<GJP.XnodeSites()/2 && xnode_coor_of_buf_data == data_nodecoor_hf1){
// 			int orig_idx = pos[4]/2 * GJP.VolNodeSites()/16 + pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + data_flav * GJP.VolNodeSites()/32;

// 			IFloat &origval = cur_data_buf[orig_idx];
// 			IFloat newval = LRG.Urand(FIVE_D);
// 			if(origval != newval){ 
// 			  printf("Node %d: 5D RNG disparity: (%d %d %d %d %d): orig %f new %f (idx %d)\n",GJP.XnodeCoor(),pos[0],pos[1],pos[2],pos[3],pos[4],origval,newval,orig_idx); exit(-1);
// 			}
// 			got_hf1 = 1;
// 		      }
// 		    }
// 		  }
// 		}
// 	      }
// 	    }
// 	  }
// 	  if(got_hf1 && got_hf2) nodes_unhappy = 0.0;
// 	  else nodes_unhappy = 1.0;
	
// 	  glb_sum(&nodes_unhappy);
// 	  if(!UniqueID()) printf("nodes_unhappy = %f\n",nodes_unhappy);

// 	  if(nodes_unhappy!=0.0){
// 	    if(!UniqueID()) printf("Passing data left\n");
// 	    cur_data_buf = recv_buf_p;
// 	    getPlusData((IFloat *)recv_buf_p, (IFloat *)send_buf_p, buf_size/sizeof(Float), 0);
// 	    xnode_coor_of_buf_data = (xnode_coor_of_buf_data+1) % GJP.Xnodes();
// 	    //swap buffers over for next send
// 	    Float *tmp = recv_buf_p;
// 	    recv_buf_p = send_buf_p;
// 	    send_buf_p = tmp;

// 	    // for(int i=0;i<n_rgen_4d;i++){
// 	    //   if(printnode) printf("Node %d: post-pass recv_buf site %d = %f\n",GJP.XnodeCoor(),i,cur_data_buf[i]);
// 	    // }
// 	  }


// 	}
// 	pfree(recv_buf);
// 	pfree(send_buf);
//       }//end of 5d RNG check
//     }
//     printf("Passed RNG test\n");
//   }

// #endif






// #define DOUBLE_RNG_TEST
// #ifdef DOUBLE_RNG_TEST
//   IFloat orig_4d_rands[2][GJP.VolNodeSites()/16]; //one per RNG
//   IFloat orig_5d_rands[2][GJP.SnodeSites()/2*GJP.VolNodeSites()/16];
//   if(gparity_X && !gparity_Y){
//     //generate fields of 4D and 5D random numbers and reset the generator. After lattice doubling compare the new and old random fields to ensure they are the same
//     LatRanGen LRGbak(LRG);
//     int pos[5];

//     for(int flav=0;flav<2;flav++){
//       for(pos[4]=0;pos[4]<GJP.SnodeSites();pos[4]+=2){
// 	for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]+=2){
// 	  for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]+=2){
// 	    for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]+=2){
// 	      for(pos[0]=0;pos[0]<GJP.XnodeSites();pos[0]+=2){
// 		LRG.AssignGenerator(pos,flav);
// 		if(pos[4]==0){
// 		  int idx = pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2+GJP.ZnodeSites()/2*pos[3]/2));
// 		  orig_4d_rands[flav][idx] = LRG.Urand(FOUR_D);//do the 4D RNG
// 		}
// 		int idx = pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2+GJP.ZnodeSites()/2*(pos[3]/2+GJP.SnodeSites()/2*pos[4]/2)));
// 		orig_5d_rands[flav][idx] = LRG.Urand(FIVE_D);//do the 5D RNG
// 	      }
// 	    }
// 	  }
// 	}
//       }
//     }
//     LRG = LRGbak;
//   }
// #endif
