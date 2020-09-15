#ifndef _TESTS_GPARITY_AMA
#define _TESTS_GPARITY_AMA

CPS_START_NAMESPACE

bool exit_on_error = true;

void setupDoArg(DoArg &do_arg, int size[5], int ngp, bool verbose = true){
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

  BndCndType* bc[3] = { &do_arg.x_bc, &do_arg.y_bc, &do_arg.z_bc };
  for(int i=0;i<ngp;i++){ 
    *(bc[i]) = BND_CND_GPARITY;
  }
}


void setup_bfmargs(BfmArg &dwfa, const BfmSolverType solver, const int nthreads, const double mobius_scale = 2.){
  dwfa.verbose=1;
  dwfa.reproduce=0;
  dwfa.precon_5d = 1;
  if(solver == BFM_HmCayleyTanh){
    dwfa.precon_5d = 0; //mobius uses 4d preconditioning
    dwfa.mobius_scale = mobius_scale;
  }else if(solver != BFM_DWF){
    ERR.General("","setup_bfmargs","Unknown solver\n");
  }

  dwfa.threads = nthreads;
  dwfa.Ls = GJP.SnodeSites();
  dwfa.solver = solver;
  dwfa.M5 = toDouble(GJP.DwfHeight());
  dwfa.mass = toDouble(0.04);
  dwfa.Csw = 0.0;
  dwfa.max_iter = 200000;
  dwfa.residual = 1e-08;
}

void setupFixGaugeArg(FixGaugeArg &r){
  r.fix_gauge_kind = FIX_GAUGE_COULOMB_T;
  r.hyperplane_start = 0;
  r.hyperplane_step = 1;
  r.hyperplane_num = GJP.TnodeSites()*GJP.Tnodes();
  r.stop_cond = 1e-06;
  r.max_iter_num = 60000;
}

void compareProps(QPropW &A, QPropW &B, const std::string label, Float tol = 1e-12){
  Float fail = 0;
  for(int f=0;f<GJP.Gparity()+1;f++){
    for(int x=0;x<GJP.VolNodeSites();x++){
      WilsonMatrix & Am = A.SiteMatrix(x,f);
      WilsonMatrix & Bm = B.SiteMatrix(x,f);

      Float* Amp = (Float*)&Am;
      Float* Bmp = (Float*)&Bm;

      for(int i=0;i<24;i++){
	if(fabs(Amp[i]- Bmp[i]) > tol){
	  printf("Node %d fail f=%d x=%d i=%d : %g vs %g\n",UniqueID(),f,x,i,Amp[i],Bmp[i]);
	  fail = 1.;
	}
      }
    }
  }
  glb_sum(&fail);
  if(fail!=0.0){
    if(!UniqueID()) printf("Prop comparison '%s' failed\n",label.c_str());
    if(exit_on_error) exit(-1);
  }else{
    if(!UniqueID()) printf("Prop comparison '%s' passed\n",label.c_str());
  }
}



PropagatorContainer & computePropagatorOld(const std::string &tag, const double mass, const double stop_prec, const int t, const int flav, const int p[3], const BndCndType time_bc, Lattice &latt, const std::string tag_otherflav = "", BFM_Krylov::Lanczos_5d<double> *deflate = NULL){ 
  if(deflate != NULL) ERR.General("","computePropagatorOld","Deflation not yet implemented\n");

  PropagatorArg parg;
  parg.generics.type = QPROPW_TYPE;
  parg.generics.tag = strdup(tag.c_str());
  parg.generics.mass = mass;
  for(int i=0;i<3;i++)
    parg.generics.bc[i] = GJP.Bc(i);
  parg.generics.bc[3] = time_bc;
  
  int len = 6 + (tag_otherflav != "" ? 1 : 0);
  parg.attributes.attributes_len = len;
  parg.attributes.attributes_val = new AttributeContainer[len];

  parg.attributes.attributes_val[0].type = WALL_SOURCE_ATTR;
  WallSourceAttrArg &wa = parg.attributes.attributes_val[0].AttributeContainer_u.wall_source_attr;
  wa.t = t;

  parg.attributes.attributes_val[1].type = GPARITY_FLAVOR_ATTR;
  GparityFlavorAttrArg &fa = parg.attributes.attributes_val[1].AttributeContainer_u.gparity_flavor_attr;
  fa.flavor = flav;

  parg.attributes.attributes_val[2].type = CG_ATTR;
  CGAttrArg &cga = parg.attributes.attributes_val[2].AttributeContainer_u.cg_attr;
  cga.max_num_iter = 10000;
  cga.stop_rsd = stop_prec;
  cga.true_rsd = stop_prec;

  parg.attributes.attributes_val[3].type = GAUGE_FIX_ATTR;
  GaugeFixAttrArg &gfa = parg.attributes.attributes_val[3].AttributeContainer_u.gauge_fix_attr;
  gfa.gauge_fix_src = 1;
  gfa.gauge_fix_snk = 0;

  parg.attributes.attributes_val[4].type = MOMENTUM_ATTR;
  MomentumAttrArg &ma = parg.attributes.attributes_val[4].AttributeContainer_u.momentum_attr;
  memcpy(ma.p,p,3*sizeof(int));

  parg.attributes.attributes_val[5].type = STORE_MIDPROP_ATTR;
  
  //STORE_MIDPROP_ATTR
  //store_midprop_attr
  //StoreMidpropAttrArg

  if(tag_otherflav != ""){
    parg.attributes.attributes_val[6].type = GPARITY_OTHER_FLAV_PROP_ATTR;
    GparityOtherFlavPropAttrArg &gofa = parg.attributes.attributes_val[6].AttributeContainer_u.gparity_other_flav_prop_attr;
    gofa.tag = strdup(tag_otherflav.c_str());
  }
  PropagatorContainer::fbfm_assume_bc_applied = false; //for AMA code the lattice by default does not have the BC applied.
  return PropManager::addProp(parg);
}

PropagatorContainer & computeCombinedPropagatorOld(const std::string &tag, const double mass, const PropCombination comb, const std::string &tag_A, const std::string &tag_B,  const std::string tag_otherflav = ""){
  PropagatorArg parg;
  parg.generics.type = QPROPW_TYPE;
  parg.generics.tag = strdup(tag.c_str());
  parg.generics.mass = mass;
  for(int i=0;i<4;i++)
    parg.generics.bc[i] = GJP.Bc(i);

  int len = 1 + (tag_otherflav != "" ? 1 : 0);
  parg.attributes.attributes_len = len;
  parg.attributes.attributes_val = new AttributeContainer[len];

  parg.attributes.attributes_val[0].type = PROP_COMBINATION_ATTR;
  PropCombinationAttrArg &pc = parg.attributes.attributes_val[0].AttributeContainer_u.prop_combination_attr;
  pc.prop_A = strdup(tag_A.c_str());
  pc.prop_B = strdup(tag_B.c_str());
  pc.combination = comb;

  //Other necessary attributes will be copied over if not specified here

  if(tag_otherflav != ""){
    parg.attributes.attributes_val[1].type = GPARITY_OTHER_FLAV_PROP_ATTR;
    GparityOtherFlavPropAttrArg &gofa = parg.attributes.attributes_val[1].AttributeContainer_u.gparity_other_flav_prop_attr;
    gofa.tag = strdup(tag_otherflav.c_str());
  }
  return PropManager::addProp(parg);
}

inline double reldiff(const double a, const double b){
  double reldiff = fabs( 2.*(a-b)/(a+b) );
  if(a==0.0 && b==0.0) reldiff =  0.;
  return reldiff;
}


inline double largest_rel_diff_reim(const Complex &ca, const Complex &cb){
  double reldiffre = reldiff(ca.real(),cb.real());
  double reldiffim = reldiff(ca.imag(),cb.imag());
  return reldiffre > reldiffim ? reldiffre : reldiffim;
}


bool test_equals(const SpinColorFlavorMatrix &a, const SpinColorFlavorMatrix &b, const double &eps){
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      for(int aa=0;aa<3;aa++){
	for(int bb=0;bb<3;bb++){
	  for(int f0=0;f0<2;f0++){
	    for(int f1=0;f1<2;f1++){
	    
	      Complex ca = a(i,aa,f0,j,bb,f1);
	      Complex cb = b(i,aa,f0,j,bb,f1);
	      if(largest_rel_diff_reim(ca,cb) > eps){		
		printf("FAIL %d %d %d %d %d %d (%g %g) (%g %g), reldiff (%g, %g)\n",i,aa,f0,j,bb,f1,ca.real(),ca.imag(),cb.real(),cb.imag(), reldiff(ca.real(),cb.real()),reldiff(ca.imag(),cb.imag()) );
		return false;
	      }
	    }
	  }
	}
      }
    }
  }
  return true;
}

void print2(const SpinColorFlavorMatrix &a,const SpinColorFlavorMatrix &b){
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      for(int aa=0;aa<3;aa++){
	for(int bb=0;bb<3;bb++){
	  for(int f0=0;f0<2;f0++){
	    for(int f1=0;f1<2;f1++){
	      Complex ca = a(i,aa,f0,j,bb,f1);
	      Complex cb = b(i,aa,f0,j,bb,f1);
	      printf("%d %d %d %d %d %d (%.4f %.4f) (%.4f %.4f)\n",i,aa,f0,j,bb,f1,ca.real(),ca.imag(),cb.real(),cb.imag());
	    }
	  }
	}
      }
    }
  }
}

bool test_equals(const std::complex<double> &a, const std::complex<double> &b, const double eps){
  return largest_rel_diff_reim(a,b) < eps;
}
bool test_equals(const double a, const double b, const double eps){
  return reldiff(a,b) < eps;
}

bool test_equals(const Matrix &a, const Matrix &b, const double &eps){
  for(int aa=0;aa<3;aa++){
    for(int bb=0;bb<3;bb++){
      Complex ca = a(aa,bb);
      Complex cb = b(aa,bb);
      if(largest_rel_diff_reim(ca,cb) > eps){		
	printf("FAIL %d %d (%g %g) (%g %g), reldiff (%g, %g)\n",aa,bb,ca.real(),ca.imag(),cb.real(),cb.imag(), reldiff(ca.real(),cb.real()),reldiff(ca.imag(),cb.imag()) );
	return false;
      }
    }
  }
  return true;
}


QPropWPointSrc* computePointPropagator(const double mass, const double stop_prec, const int t, const int flav, const int p[3], const BndCndType time_bc, const bool store_midprop, 
				Lattice &latt,  BFM_Krylov::Lanczos_5d<double> *deflate = NULL){ 
  if(!UniqueID()) printf("Computing point propagator\n");
  LatticeContainer lat_bak; lat_bak.Get(latt);

  multi1d<float> *eval_conv = NULL;

  if(deflate != NULL){
    if(latt.Fclass() != F_CLASS_BFM && latt.Fclass() != F_CLASS_BFM_TYPE2)
      ERR.General("","computePropagator","Deflation only implemented for Fbfm\n");
    if(Fbfm::use_mixed_solver){
      //Have to convert evals to single prec
      eval_conv = new multi1d<float>(deflate->bl.size());
      for(int i=0;i<eval_conv->size();i++) eval_conv->operator[](i) = deflate->bl[i];
      dynamic_cast<Fbfm&>(latt).set_deflation<float>(&deflate->bq,eval_conv,0); //last argument is really obscure - it's the number of eigenvectors subtracted from the solution to produce a high-mode inverse - we want zero here
    }else dynamic_cast<Fbfm&>(latt).set_deflation(&deflate->bq,&deflate->bl,0);
  }

  CommonArg c_arg;
  
  CgArg cg;
  cg.mass = mass;
  cg.max_num_iter = 10000;
  cg.stop_rsd = stop_prec;
  cg.true_rsd = stop_prec;
  cg.RitzMatOper = NONE;
  cg.Inverter = CG;
  cg.bicgstab_n = 0;

  QPropWArg qpropw_arg;
  qpropw_arg.cg = cg;
  qpropw_arg.x = 0;
  qpropw_arg.y = 0;
  qpropw_arg.z = 0;
  qpropw_arg.t = t;
  qpropw_arg.flavor = flav; 
  qpropw_arg.ensemble_label = "ens";
  qpropw_arg.ensemble_id = "ens_id";
  qpropw_arg.StartSrcSpin = 0;
  qpropw_arg.EndSrcSpin = 4;
  qpropw_arg.StartSrcColor = 0;
  qpropw_arg.EndSrcColor = 3;
  qpropw_arg.gauge_fix_src = 1;
  qpropw_arg.gauge_fix_snk = 0;
  qpropw_arg.store_midprop = store_midprop ? 1 : 0; //for mres

  //Switching boundary conditions is poorly implemented in Fbfm (and probably FGrid)
  //For traditional lattice types the BC is applied by modififying the gauge field when the Dirac operator is created and reverting when destroyed. This only ever happens internally - no global instance of the Dirac operator exists
  //On the other hand, Fbfm does all its inversion internally and doesn't instantiate a CPS Dirac operator. We therefore have to manually force Fbfm to change its internal gauge field by applying BondCond
  bool is_wrapper_type = ( latt.Fclass() == F_CLASS_BFM || latt.Fclass() == F_CLASS_BFM_TYPE2 ); //I hate this!

  BndCndType init_tbc = GJP.Tbc();
  BndCndType target_tbc = time_bc;

  GJP.Bc(3,target_tbc);
  if(is_wrapper_type) latt.BondCond();  //Apply new BC to internal gauge fields

  QPropWPointSrc* ret = new QPropWPointSrc(latt,&qpropw_arg,&c_arg);

  //Restore the BCs
  if(is_wrapper_type) latt.BondCond();  //unapply existing BC
  GJP.Bc(3,init_tbc);

  if(deflate != NULL) dynamic_cast<Fbfm&>(latt).unset_deflation();
  if(eval_conv !=NULL) delete eval_conv;

  //Check to make sure the lattice has not been somehow changed by this function
  Float fail = 0.;
  for(int i=0;i<(GJP.Gparity()+1)*4*GJP.VolNodeSites();i++)
    if(!test_equals(latt.GaugeField()[i], lat_bak.GaugeField()[i],1e-12))
      fail = 1.0;
  glb_sum_five(&fail);
  if(fail != 0.)
    ERR.General("","computePointPropagator","Lattice checksum fail\n");

  return ret;
}

QPropW* computePropagator(const SourceType src_type, const double mass, const double stop_prec, const int t, const int flav, const int p[3], const BndCndType time_bc, const bool store_midprop, 
			  Lattice &latt,  BFM_Krylov::Lanczos_5d<double> *deflate = NULL){ 
  if(src_type == WALL) return (QPropW*)computePropagator(mass,stop_prec,t,flav,p,time_bc,store_midprop,latt,deflate);
  else if(src_type == POINT) return (QPropW*)computePointPropagator(mass,stop_prec,t,flav,p,time_bc,store_midprop,latt,deflate);
  else{
    ERR.General("","computePropagator","Unknown type\n");
  }
}

void test_lattice_doubler(Lattice &lattice, const DoArg &do_arg, const bool fix_gauge_unity){
  if(!UniqueID()) printf("Testing lattice doubler\n");

  SourceType src_type = POINT; //WALL;

  //Generate P+A propagator
  int p[3] = {0,0,0}; for(int i=0;i<3;i++) if(GJP.Bc(i) == BND_CND_GPARITY) p[i] = 1;
  int tsrc = 0;

  QPropW* prop_f0_P = computePropagator(src_type,0.04,1e-10,tsrc,0, p, BND_CND_PRD, false, lattice);
  QPropW* prop_f0_A = computePropagator(src_type,0.04,1e-10,tsrc,0, p, BND_CND_APRD, false, lattice);
  QPropW* prop_f1_P = NULL;
  QPropW* prop_f1_A = NULL;
  if(GJP.Gparity()){
    prop_f1_P = computePropagator(src_type,0.04,1e-10,tsrc,1, p, BND_CND_PRD, false, lattice);
    prop_f1_A = computePropagator(src_type,0.04,1e-10,tsrc,1, p, BND_CND_APRD, false, lattice);
  }
  PropWrapper prop_P(prop_f0_P,prop_f1_P);
  PropWrapper prop_A(prop_f0_A,prop_f1_A);

  PropWrapper prop_F = PropWrapper::combinePA(prop_P,prop_A,CombinationF);
  PropWrapper prop_B = PropWrapper::combinePA(prop_P,prop_A,CombinationB);

  // //Compute tr(G^dag G) at each time for comparison
  // fVector<Rcomplex> Fcor(2*Lt); //Put F and B together using F(t+Lt) = B(t)
  // for(int t=0;t<GJP.TnodeSites();t++){
  //   int t_glb_F = t + GJP.TnodeCoor()*GJP.TnodeSites();
  //   int t_glb_B = t_glb_F + Lt;
  //   for(int x=0;x<GJP.VolNodeSites()/GJP.TnodeSites();x++){
  int Lt = GJP.TnodeSites()*GJP.Tnodes(); //pre-doubling
  fVector<Rcomplex> Fnorms(2*Lt); //Put F and B together using F(t+Lt) = B(t)
  fVector<Rcomplex> Fnorms_f00(2*Lt);
  {
    SpinColorFlavorMatrix tmp_scf;
    WilsonMatrix tmp_sc;
    int vol3d = GJP.VolNodeSites()/GJP.TnodeSites();

    QPropW* prop_F_f0 = prop_F.getPtr(0);
    QPropW* prop_B_f0 = prop_B.getPtr(0);

    for(int t=0;t<GJP.TnodeSites();t++){
      int t_glb_F = t + GJP.TnodeCoor()*GJP.TnodeSites();
      int t_glb_B = t_glb_F + Lt;
      for(int x=0;x<vol3d;x++){
	if(GJP.Gparity()){
	  prop_F.siteMatrix(tmp_scf,x + vol3d*t);
	  Fnorms(t_glb_F) += tmp_scf.norm();
	  
	  prop_B.siteMatrix(tmp_scf,x + vol3d*t);
	  Fnorms(t_glb_B) += tmp_scf.norm();
	}else{
	  prop_F.siteMatrix(tmp_sc,x + vol3d*t);
	  Fnorms(t_glb_F) += tmp_sc.norm();
	  
	  prop_B.siteMatrix(tmp_sc,x + vol3d*t);
	  Fnorms(t_glb_B) += tmp_sc.norm();
	}

	tmp_sc = prop_F_f0->SiteMatrix(x + vol3d*t,0);
	Fnorms_f00(t_glb_F) += tmp_sc.norm();

	tmp_sc = prop_B_f0->SiteMatrix(x + vol3d*t,0);
	Fnorms_f00(t_glb_B) += tmp_sc.norm();
      }
    }
    Fnorms.nodeSum();
    Fnorms_f00.nodeSum();
  }


  //Double the lattice T extent
  lattice.FixGaugeFree();
  LatticeTimeDoubler doubler;
  doubler.doubleLattice(lattice,do_arg);
    
  if(!UniqueID()){ printf("Lattice double complete, testing\n"); fflush(stdout); }

  //Test by shifting the doubled lattice by Lt/2 and making sure the fields match
  int lat_size = 18*4*GJP.VolNodeSites()*(GJP.Gparity()+1);
  Lt = GJP.TnodeSites()*GJP.Tnodes(); //after doubling
  LatticeContainer bak;
  bak.Get(lattice);
    
  Float fail = 0.0;

  if(GJP.Tnodes() > 1){
    if(!UniqueID()) printf("Multi Tnode test\n");
    Float* buf1 = (Float*)malloc(lat_size * sizeof(Float));
    Float* buf2 = (Float*)malloc(lat_size * sizeof(Float));
    memcpy((void*)buf1, (void*)lattice.GaugeField(), lat_size*sizeof(Float));

    Float *send = buf2;
    Float *recv = buf1;    
    for(int i=0;i<GJP.Tnodes()/2;i++){
      Float* tmp = send;
      send = recv;
      recv = tmp;
      getPlusData(recv, send, lat_size, 3);
    }

    Float* orig_p = (Float*)bak.GaugeField();
    for(int i=0;i<lat_size;i++){
      if(orig_p[i] != recv[i]){
	printf("Node %d index %d expected %g got %g\n",UniqueID(),i,orig_p[i],recv[i]); fflush(stdout); fail = 1.0;;
      }
    }
    glb_sum_five(&fail);
    free(buf1);
    free(buf2);
  }else{
    if(!UniqueID()) printf("Single Tnode test\n");
    Tshift4D((Float*)lattice.GaugeField(), 18*4, GJP.TnodeSites()/2);
    Float* orig_p = (Float*)bak.GaugeField();
    Float* lat_p = (Float*)lattice.GaugeField();
    for(int i=0;i<lat_size;i++){
      if(orig_p[i] != lat_p[i]){
	printf("Node %d index %d expected %g got %g\n",UniqueID(),i,orig_p[i],lat_p[i]); fflush(stdout); fail = 1.0;;
      }
    }
    glb_sum_five(&fail);
  }
  if(fail != 0.0){
    if(!UniqueID()){ printf("Lattice double test fail\n"); exit(-1); }
  }else{
    if(!UniqueID()){ printf("Lattice double test pass\n"); }
  }
  //Regenerate gauge fixing matrix
  FixGaugeArg fix_gauge_arg; setupFixGaugeArg(fix_gauge_arg);
  CommonArg common_arg;
  AlgFixGauge fix_gauge(lattice,&common_arg,&fix_gauge_arg);
  if(fix_gauge_unity){
    gaugeFixUnity(lattice,fix_gauge_arg);
  }else{
    fix_gauge.run();
  }

  //Examine P propagators and compare to P+A
  QPropW* prop_f0_P_dbl = computePropagator(src_type,0.04,1e-10,tsrc,0, p, BND_CND_PRD, false, lattice);
  QPropW* prop_f1_P_dbl = GJP.Gparity() ? computePropagator(src_type,0.04,1e-10,tsrc,1, p, BND_CND_PRD, false, lattice) : NULL;
  PropWrapper prop_P_dbl(prop_f0_P_dbl,prop_f1_P_dbl);

  //Simple start, just compare matrix norms over the lattice
  fVector<Rcomplex> Pnorms(Lt);
  fVector<Rcomplex> Pnorms_f00(Lt);
  {
    SpinColorFlavorMatrix tmp_scf;
    WilsonMatrix tmp_sc;
    int vol3d = GJP.VolNodeSites()/GJP.TnodeSites();
    
    QPropW* prop_P_f0 = prop_P_dbl.getPtr(0);

    for(int t=0;t<GJP.TnodeSites();t++){
      int t_glb = t + GJP.TnodeCoor()*GJP.TnodeSites();

      for(int x=0;x<vol3d;x++){
	if(GJP.Gparity()){
	  prop_P_dbl.siteMatrix(tmp_scf,x + vol3d*t);
	  Pnorms(t_glb) += tmp_scf.norm();
	}else{
	  prop_P_dbl.siteMatrix(tmp_sc,x + vol3d*t);
	  Pnorms(t_glb) += tmp_sc.norm();
	}	
	tmp_sc = prop_P_f0->SiteMatrix(x + vol3d*t,0);
	Pnorms_f00(t_glb) += tmp_sc.norm();
      }
    }
    Pnorms.nodeSum();
    Pnorms_f00.nodeSum();
  }

  if(!UniqueID()){
    printf("Prop norms\n");
    for(int t=0;t<Lt;t++){      
      Float f = Fnorms(t).real();
      Float p = Pnorms(t).real();
      Float reldiff = 2.*(f-p)/(f+p);
      if(f==0. && p==0.) reldiff = 0.;
      printf("%d %g %g  diff %g\n", f, p, reldiff);
    }

    printf("Prop_00 norms\n");
    for(int t=0;t<Lt;t++){      
      Float f = Fnorms_f00(t).real();
      Float p = Pnorms_f00(t).real();
      Float reldiff = 2.*(f-p)/(f+p);
      if(f==0. && p==0.) reldiff = 0.;
      printf("%d %g %g  diff %g\n", f, p, reldiff);
    }

  }


  // //Compare tr( G^dag G )
  // fVector<Rcomplex> Fcor(Lt);
  // fVector<Rcomplex> Pcor(Lt);
  // for(int t=0;t<GJP.TnodeSites();t++){
  //   int t_glb = t + GJP.TnodeCoor()*GJP.TnodeSites();
  //   PropWrapper &FBuse = t_glb < Lt/2

  //   for(int x=0;x<GJP.VolNodeSites()/GJP.TnodeSites();x++){
      



  exit(0);
}


int run_tests(int argc,char *argv[])
{
  Start(&argc, &argv);
  int ngp;
  { std::stringstream ss; ss << argv[1]; ss >> ngp; }

  if(!UniqueID()) printf("Doing G-parity in %d directions\n",ngp);

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
  int nthreads = 1;

  BfmSolverType solver = BFM_DWF;//DWF; HmCayleyTanh
  double mobius_scale = 2.;

  bool lanc_args_setup = false;
  Fbfm::use_mixed_solver = false;
  double tol = 1e-8;

  bool test_t_doubler = false;
  bool fix_gauge_unity = false;

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
      for(int d=0;d<5;d++)
	size[d] = toInt(argv[i+1+d]);
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
    }else if( strncmp(cmd,"-verbose",15) == 0){
      verbose=true;
      i++;
    }else if( strncmp(cmd,"-test_t_doubler",15) == 0){
      test_t_doubler = true;
      if(!UniqueID()) printf("Testing t doubler\n");
      i++;
    }else if( strncmp(cmd,"-fix_gauge_unity",15) == 0){
      fix_gauge_unity = true;
      if(!UniqueID()) printf("Setting gauge fixing matrices to unity\n");
      i++;
    }else if( strncmp(cmd,"-disable_exit_on_error",15) == 0){
      exit_on_error = false;
      i++;
    }else if( strncmp(cmd,"-nthread",15) == 0){
      nthreads = toInt(argv[i+1]);
      i+=2;
    }else if( strncmp(cmd,"-unit_gauge",15) == 0){
      unit_gauge=true;
      i++;
    }else if( strncmp(cmd,"-mobius",15) == 0){
      solver = BFM_HmCayleyTanh;
      i++;
    }else if( strncmp(cmd,"-tolerance",15) == 0){
      std::stringstream ss; ss << argv[i+1];
      ss >> tol;
      if(!UniqueID()) printf("Set tolerance to %g\n",tol);
      i+=2;
    }else if( strncmp(cmd,"-mobius_scale",15) == 0){
      std::stringstream ss; ss << argv[i+1];
      ss >> mobius_scale;
      if(!UniqueID()) printf("Set Mobius scale to %g\n",mobius_scale);
      i+=2;
    /* }else if( strncmp(cmd,"-load_lanc_arg",20) == 0){ */
    /*   if(!lanc_arg.Decode("lanc_arg.vml","lanc_arg")){VRB.Result("","","Can't open lanc_arg.vml!\n");exit(1);} */
    /*   if(!lanc_arg_s.Decode("lanc_arg_s.vml","lanc_arg")){VRB.Result("","","Can't open lanc_arg_s.vml!\n");exit(1);} */
    /*   lanc_args_setup = true; */
    /*   i++; */
    }else if( std::string(argv[i]) == "-use_mixed_solver" ){
      Fbfm::use_mixed_solver = true;
      i++;
    }else{
      if(UniqueID()==0) printf("Unrecognised argument: %s\n",cmd);
      exit(-1);
    }
  }

  printf("Lattice size is %d %d %d %d\n",size[0],size[1],size[2],size[3],size[4]);

  CommonArg common_arg;
  DoArg do_arg;  setupDoArg(do_arg,size,ngp,verbose);
  if(test_t_doubler) do_arg.t_bc = BND_CND_PRD;

  GJP.Initialize(do_arg);

#if TARGET == BGQ
  LRG.setSerial();
#endif
  LRG.Initialize(); //usually initialised when lattice generated, but I pre-init here so I can load the state from file

  //cps_qdp_init(&argc,&argv);
  //Chroma::initialize(&argc,&argv);

  BfmArg bfm_arg;
  setup_bfmargs(bfm_arg,solver,nthreads,mobius_scale);
  init_fbfm(&argc,&argv,bfm_arg);

  if(UniqueID()==0)
    if(Fbfm::use_mixed_solver) printf("Using Fbfm mixed precision solver\n");
    else printf("Using Fbfm double precision solver\n");

  GnoneFbfm lattice;
  //GwilsonFdwf lattice;

  if(load_lrg){
    if(UniqueID()==0) printf("Loading RNG state from %s\n",load_lrg_file);
    LRG.Read(load_lrg_file,32);
  }
  if(save_lrg){
    if(UniqueID()==0) printf("Writing RNG state to %s\n",save_lrg_file);
    LRG.Write(save_lrg_file,32);
  }					       
  if(!load_config){
    printf("Creating gauge field\n");
    if(!unit_gauge) lattice.SetGfieldDisOrd();
    else lattice.SetGfieldOrd();
  }else{
    ReadLatticeParallel readLat;
    if(UniqueID()==0) printf("Reading: %s (NERSC-format)\n",load_config_file);
    readLat.read(lattice,load_config_file);
    if(UniqueID()==0) printf("Config read.\n");
  }
  if(save_config){
    if(UniqueID()==0) printf("Saving config to %s\n",save_config_file);

    QioArg wt_arg(save_config_file,0.001);
    
    wt_arg.ConcurIONumber=32;
    WriteLatticeParallel wl;
    wl.setHeader("disord_id","disord_label",0);
    wl.write(lattice,wt_arg);
    
    if(!wl.good()) ERR.General("main","()","Failed write lattice %s",save_config_file);

    if(UniqueID()==0) printf("Config written.\n");
  }

  FixGaugeArg fix_gauge_arg; setupFixGaugeArg(fix_gauge_arg);

  AlgFixGauge fix_gauge(lattice,&common_arg,&fix_gauge_arg);
  if(fix_gauge_unity){
    gaugeFixUnity(lattice,fix_gauge_arg);
  }else{
    fix_gauge.run();
  }


  if(test_t_doubler){
    test_lattice_doubler(lattice,do_arg,fix_gauge_unity);
  }

  //Here on is Gparity only
  assert(ngp > 0);

  //Generate props
  double prec = 1e-6;
  ThreeMomentum p;
  for(int i=0;i<3;i++) if(GJP.Bc(i) == BND_CND_GPARITY) p(i) = 1;
  ThreeMomentum mp(-p);

  int tsrc = 0;
  PropagatorContainer &prop_f0_pplus = computePropagatorOld("prop_f0_pplus",0.01,prec,tsrc,0,p.ptr(),BND_CND_APRD,lattice, "prop_f1_pplus");
  PropagatorContainer &prop_f1_pplus = computePropagatorOld("prop_f1_pplus",0.01,prec,tsrc,1,p.ptr(),BND_CND_APRD,lattice, "prop_f0_pplus");

  PropagatorContainer &prop_f0_pminus = computePropagatorOld("prop_f0_pminus",0.01,prec,tsrc,0,mp.ptr(),BND_CND_APRD,lattice, "prop_f1_pminus");
  PropagatorContainer &prop_f1_pminus = computePropagatorOld("prop_f1_pminus",0.01,prec,tsrc,1,mp.ptr(),BND_CND_APRD,lattice, "prop_f0_pminus");

  // lattice.BondCond();
  // PropagatorContainer::fbfm_assume_bc_applied = true;

  PropManager::calcProps(lattice);

  //Check complex conj relation basic
  for(int i=0;i<GJP.VolNodeSites();i++){
    SpinColorFlavorMatrix pm_calc, pm_flip;
    //pw_prop_pplus.siteMatrix(pm_calc,i); //p- formed from p- f0 and p- f1

    pm_calc.generate(prop_f0_pplus.convert<QPropWcontainer>(), prop_f1_pplus.convert<QPropWcontainer>(), lattice,i);
    pm_flip.generate_from_cconj_pair(prop_f0_pplus.convert<QPropWcontainer>(), prop_f0_pminus.convert<QPropWcontainer>(), lattice,i);
    
    bool eq = test_equals(pm_calc,pm_flip,tol);
    double fail = eq ? 0. : 1.;
    glb_sum(&fail);
    if(fail != 0.0){
      if(!UniqueID()){
	print2(pm_calc,pm_flip);	
      }
      if(!UniqueID()) printf("cconj reln test plus failed\n");
      if(exit_on_error) exit(-1);
    }
  }
  if(!UniqueID()) printf("cconj reln test plus passed\n");

  //Check we can reproduce the propagator using the new code
  QPropWMomSrc* prop_f0_pplus_test = computePropagator(0.01, prec, tsrc, 0, p, BND_CND_APRD, false,lattice);
  
  compareProps(QPropWcontainer::verify_convert(prop_f0_pplus,"","").getProp(lattice), *prop_f0_pplus_test, "f0 p+ test", 1e-12);

  PropagatorContainer &prop_f0_pplus_P = computePropagatorOld("prop_f0_pplus_P",0.01,prec,tsrc,0,p.ptr(),BND_CND_PRD,lattice, "prop_f1_pplus_P");
  QPropWMomSrc* prop_f0_pplus_P_test = computePropagator(0.01, prec, tsrc, 0, p, BND_CND_PRD, false, lattice);

  compareProps(QPropWcontainer::verify_convert(prop_f0_pplus_P,"","").getProp(lattice), *prop_f0_pplus_P_test, "f0 p+ P test", 1e-12);

  //Made propwrappers
  PropWrapper pw_prop_pplus(&QPropWcontainer::verify_convert(prop_f0_pplus,"","").getProp(lattice),
			    &QPropWcontainer::verify_convert(prop_f1_pplus,"","").getProp(lattice));
  
  PropWrapper pw_prop_pminus(&QPropWcontainer::verify_convert(prop_f0_pminus,"","").getProp(lattice),
			     &QPropWcontainer::verify_convert(prop_f1_pminus,"","").getProp(lattice));

  PropSiteMatrixStandard pw_prop_pplus_getter(pw_prop_pplus, BND_CND_APRD, tsrc);
  PropSiteMatrixStandard pw_prop_pminus_getter(pw_prop_pminus, BND_CND_APRD, tsrc);

  //Check that the propagator Tbc combinations are computed properly
  {
    PropagatorContainer &prop_f1_pplus_P = computePropagatorOld("prop_f1_pplus_P",0.01,prec,tsrc,1,p.ptr(),BND_CND_PRD,lattice, "prop_f0_pplus_P");
    
    PropagatorContainer & prop_f0_pplus_F = computeCombinedPropagatorOld("prop_f0_pplus_F",0.01,A_PLUS_B,"prop_f0_pplus_P","prop_f0_pplus",  "prop_f1_pplus_F");
    PropagatorContainer & prop_f1_pplus_F = computeCombinedPropagatorOld("prop_f1_pplus_F",0.01,A_PLUS_B,"prop_f1_pplus_P","prop_f1_pplus",  "prop_f0_pplus_F");

    PropagatorContainer & prop_f0_pplus_B = computeCombinedPropagatorOld("prop_f0_pplus_F",0.01,A_MINUS_B,"prop_f0_pplus_P","prop_f0_pplus",  "prop_f1_pplus_B");
    PropagatorContainer & prop_f1_pplus_B = computeCombinedPropagatorOld("prop_f1_pplus_F",0.01,A_MINUS_B,"prop_f1_pplus_P","prop_f1_pplus",  "prop_f0_pplus_B");

  
    PropWrapper pw_prop_pplus_P(&QPropWcontainer::verify_convert(prop_f0_pplus_P,"","").getProp(lattice),
				&QPropWcontainer::verify_convert(prop_f1_pplus_P,"","").getProp(lattice));
    PropWrapper pw_prop_pplus_F = PropWrapper::combinePA(pw_prop_pplus_P,pw_prop_pplus,CombinationF);
    PropWrapper pw_prop_pplus_B = PropWrapper::combinePA(pw_prop_pplus_P,pw_prop_pplus,CombinationB);

    compareProps(QPropWcontainer::verify_convert(prop_f0_pplus_F,"","").getProp(lattice), *pw_prop_pplus_F.getPtr(0), "f0 F test", 1e-12);
    compareProps(QPropWcontainer::verify_convert(prop_f1_pplus_F,"","").getProp(lattice), *pw_prop_pplus_F.getPtr(1), "f1 F test", 1e-12);
    compareProps(QPropWcontainer::verify_convert(prop_f0_pplus_B,"","").getProp(lattice), *pw_prop_pplus_B.getPtr(0), "f0 B test", 1e-12);
    compareProps(QPropWcontainer::verify_convert(prop_f1_pplus_B,"","").getProp(lattice), *pw_prop_pplus_B.getPtr(1), "f1 B test", 1e-12);
    if(!UniqueID()) printf("Passed prop comb test\n");
  }
  

  //Test pion 2pt against old code
  int L[4] = {GJP.XnodeSites()*GJP.Xnodes(), GJP.YnodeSites()*GJP.Ynodes(), GJP.ZnodeSites()*GJP.Znodes(), GJP.TnodeSites()*GJP.Tnodes()};


  ContractedBilinear<SpinColorFlavorMatrix> conbil;
  std::vector<Float> p_pi_plus(3,0); //sink phase exp(i p.x)
  for(int i=0;i<ngp;i++) p_pi_plus[i] = M_PI/double(L[i]);

  std::vector<Float> p_zero(3,0);

  conbil.add_momentum(p_pi_plus);
  conbil.add_momentum(p_zero);
  conbil.calculateBilinears(lattice, "prop_f0_pminus", OpDagger, "prop_f0_pplus", OpNone);
  conbil.calculateBilinears(lattice, "prop_f0_pplus", OpDagger, "prop_f0_pplus", OpNone);

  Pion2PtSinkOp snk_op_new[] = { AX, AY, AZ, AT, P };
  int snk_spn_old[] = {1,2,4,8,0};

  //Test axial, pseudo sink pion LW 2pt funtions
  for(int oo=0;oo<5;oo++){
    //sigma3(1+sigma2) = sigma3 -i sigma1
    //Note, ordering of operators in ContractedBilinear is source, sink
    const std::vector<Rcomplex> pps3 =  conbil.getBilinear(lattice,p_pi_plus,"prop_f0_pminus", OpDagger, "prop_f0_pplus", OpNone,
  							    0, 3, snk_spn_old[oo], 3);
    const std::vector<Rcomplex> pps1 =  conbil.getBilinear(lattice,p_pi_plus,"prop_f0_pminus", OpDagger, "prop_f0_pplus", OpNone,
  							    0, 1, snk_spn_old[oo], 3);
    std::vector<Rcomplex> ppconbil(pps3);
    for(int i=0;i<pps3.size();i++) ppconbil[i] = 0.25*(pps3[i] + Complex(0,-1)*pps1[i]);
    
    ThreeMomentum p_psibar = p;
    ThreeMomentum p_psi = p;

    fMatrix<Rcomplex> ppnew(1,L[3]);
    pionTwoPointLWGparity(ppnew,0,snk_op_new[oo],p_psibar,p,pw_prop_pminus_getter,pw_prop_pplus_getter);
    
    double fail = 0;
    if(!UniqueID()) printf("Operator %d\n",oo);
    for(int t=0;t<L[3];t++){
      if(!test_equals(ppconbil[t],ppnew(0,t),tol)){
  	if(!UniqueID()) printf("FAIL t=%d, old=(%g,%g) new=(%g,%g), reldiff=(%g,%g)\n",t,
  			       std::real(ppconbil[t]),std::imag(ppconbil[t]),
  			       std::real(ppnew(0,t)),std::imag(ppnew(0,t)),
  			       reldiff(ppconbil[t].real(),ppnew(0,t).real()), reldiff(ppconbil[t].imag(),ppnew(0,t).imag()) );
  	fail = 1.;
      }else if(!UniqueID()) printf("PASS t=%d, old=(%g,%g) new=(%g,%g) reldiff=(%g %g)\n",t,
				   std::real(ppconbil[t]),std::imag(ppconbil[t]),std::real(ppnew(0,t)),std::imag(ppnew(0,t)),reldiff(ppconbil[t].real(),ppnew(0,t).real()), reldiff(ppconbil[t].imag(),ppnew(0,t).imag()));
    }
    glb_sum(&fail);
    if(fail){
      if(!UniqueID()) printf("Failed comparison of operator %d\n",oo);
      if(exit_on_error) exit(-1);
    }
  }

  //Test pion WW 2pt function
  {
    //In this test we consider the following contraction
    //The quarks each have the same sink momenta as they do at the source (up to the necessary - sign in the phase)
    // \sum_{x1,x2,y1,y2} 
    //<
    //  [[ exp(-i[-p2].x1)\bar\psi(x1,t) A exp(-i[-p1].x2)\psi(x2,t) ]]
    //  *
    //  [[ exp(-i p1.y1)\bar\psi(y1,0) B exp(-i p2.y2)\psi(y2,0) ]]
    //>
    // = 
    //Tr( 
    //   [\sum_{x1,y2} exp(-i[-p2].x1) exp(-i p2.y2) G(y2,0;x1,t)] A
    //   *
    //   [\sum_{x2,y1} exp(-i[-p1].x2) exp(-i p1.y1) G(x2,t;y1,0)] B
    //  )
    // = 
    //Tr( 
    //   g5 [\sum_{x1,y2} exp(-ip2.x1) exp(-i[-p2].y2) G(x1,t;y2,0)]^\dagger g5 A
    //   *
    //   [\sum_{x2,y1} exp(-i[-p1].x2) exp(-i p1.y1) G(x2,t;y1,0)] B
    //  )

    //We use p1=p2=p

    std::vector<Float> p_phys_units(3);
    p.latticeUnits(&p_phys_units[0]);

    std::vector<Float> mp_phys_units(3);
    mp.latticeUnits(&mp_phys_units[0]);

    ContractedWallSinkBilinearSpecMomentum<SpinColorFlavorMatrix> conwsbil;

    //Old code applies phases exp(+p.x) at sink whereas new code applies exp(-ip.x) everywhere
    //thus we need to swap out momentum conventions
    //Old code sink momenta are the phases applied outside the g5-hermiticity parentheses
    //Tr( 
    //   \sum_x1 exp(-i[-p2].x1) g5 [\sum_y2 exp(-i[-p2].y2) G(x1,t;y2,0)]^\dagger g5 A
    //   *
    //   \sum_x2 exp(-i[-p1].x2) [\sum_{x2,y1}  exp(-i p1.y1) G(x2,t;y1,0)] B
    //  )

    std::pair< std::vector<Float>,std::vector<Float> > mompair( p_phys_units, p_phys_units );
    conwsbil.add_momentum(mompair);
    conwsbil.calculateBilinears(lattice, "prop_f0_pminus", OpDagger, "prop_f0_pplus", OpNone);
    
    //Compute gauge-fixed wall sink propagators with new code
    WallSinkPropSiteMatrixStandard<SpinColorFlavorMatrix> ws_prop_dag_getter(pw_prop_pminus, BND_CND_APRD, tsrc, &p_phys_units[0], lattice); //daggered prop is given phase +p2=p inside g5-herm parentheses
    WallSinkProp<SpinColorFlavorMatrix> &ws_prop_dag = ws_prop_dag_getter.getWallSinkProp();

    WallSinkPropSiteMatrixStandard<SpinColorFlavorMatrix> ws_prop_undag_getter(pw_prop_pplus, BND_CND_APRD, tsrc, &mp_phys_units[0], lattice);
    WallSinkProp<SpinColorFlavorMatrix> &ws_prop_undag = ws_prop_undag_getter.getWallSinkProp();
    
    fMatrix<Rcomplex> ppnew(1,L[3]);

    ThreeMomentum p_psi_src = p;
    ThreeMomentum p_psi_snk = mp;

    pionTwoPointPPWWGparity(ppnew, 0, p_psi_snk, p_psi_src, ws_prop_dag_getter, ws_prop_undag_getter);

    const std::vector<Rcomplex> pps3s3 =  conwsbil.getBilinear(lattice,mompair,"prop_f0_pminus", OpDagger, "prop_f0_pplus", OpNone,
  							    0, 3, 0, 3);
    const std::vector<Rcomplex> pps1s3 =  conwsbil.getBilinear(lattice,mompair,"prop_f0_pminus", OpDagger, "prop_f0_pplus", OpNone,
  							    0, 1, 0, 3);
    const std::vector<Rcomplex> pps3s1 =  conwsbil.getBilinear(lattice,mompair,"prop_f0_pminus", OpDagger, "prop_f0_pplus", OpNone,
  							    0, 3, 0, 1);
    const std::vector<Rcomplex> pps1s1 =  conwsbil.getBilinear(lattice,mompair,"prop_f0_pminus", OpDagger, "prop_f0_pplus", OpNone,
  							    0, 1, 0, 1);

    //s3(1+s2) = s3 - is1  at source
    //s3(1-s2) = s3 + is1  at sink
    
    //[s3-is1)]o[s3+is1] = [s3 o s3] + [s1 o s1] -i [s1 o s3] +i [s3 o s1]
    std::vector<Rcomplex> ppconbil(L[3]);
    for(int i=0;i<L[3];i++) ppconbil[i] = 0.5*0.5*0.5*( pps3s3[i] + pps1s1[i] + Complex(0,-1)*pps1s3[i] + Complex(0,1)*pps3s1[i] );
        
    double fail = 0;
    if(!UniqueID()) printf("PPWW\n");
    for(int t=0;t<L[3];t++){
      if(!test_equals(ppconbil[t],ppnew(0,t),tol)){
  	if(!UniqueID()) printf("FAIL t=%d, old=(%g,%g) new=(%g,%g), reldiff=(%g,%g)\n",t,
  			       std::real(ppconbil[t]),std::imag(ppconbil[t]),
  			       std::real(ppnew(0,t)),std::imag(ppnew(0,t)),
  			       reldiff(ppconbil[t].real(),ppnew(0,t).real()), reldiff(ppconbil[t].imag(),ppnew(0,t).imag()));
  	fail = 1.;
      }else if(!UniqueID()) printf("PASS t=%d, old=(%g,%g) new=(%g,%g)\n",t,std::real(ppconbil[t]),std::imag(ppconbil[t]),std::real(ppnew(0,t)),std::imag(ppnew(0,t)) );
    }
    glb_sum(&fail);
    if(fail){
      if(!UniqueID()) printf("Failed comparison of PPWW\n");
      if(exit_on_error) exit(-1);
    }
  }



  //Test pseudoscalar flavor singlet against old code. It's stationary. We use momentum assignment    \bar\psi(p) \gamma^5 \psi(-p)  . The psi will be daggered as part of g5-hermiticity op, swapping its momentum.
  //Need projector (1-\sigma_2) because of negative \psi momentum
  {
    const std::vector<Rcomplex> pp1 =  conbil.getBilinear(lattice,p_zero,"prop_f0_pplus", OpDagger, "prop_f0_pplus", OpNone,
							  0, 0, 0, 0);
    const std::vector<Rcomplex> pps2 =  conbil.getBilinear(lattice,p_zero,"prop_f0_pplus", OpDagger, "prop_f0_pplus", OpNone,
     							    0, 2, 0, 0);
    std::vector<Rcomplex> ppconbil(L[3]);
    for(int i=0;i<L[3];i++) ppconbil[i] = 0.25*(pp1[i] - pps2[i]);

    ThreeMomentum p_psibar = p;
    ThreeMomentum p_psi = mp;

    fMatrix<Rcomplex> ppnew(1,L[3]);
    lightFlavorSingletLWGparity(ppnew,0,p_psibar,p_psi,pw_prop_pplus_getter,pw_prop_pplus_getter);
    
    double fail = 0;
    if(!UniqueID()) printf("Pseudoscalar flavor singlet\n");
    for(int t=0;t<L[3];t++){
      if(ppconbil[t].imag() > tol || ppnew(0,t).imag() > tol){
	if(!UniqueID()) printf("FAIL t=%d, old=(%g,%g) new=(%g,%g), imag parts should be zero within tol\n",t,
			       std::real(ppconbil[t]),std::imag(ppconbil[t]),
			       std::real(ppnew(0,t)),std::imag(ppnew(0,t)));
	fail = 1.;
      }else if(!test_equals(ppconbil[t].real(),ppnew(0,t).real(),tol)){ 
	if(!UniqueID()) printf("FAIL t=%d, old=(%g,%g) new=(%g,%g), reldiff=(%g,%g)\n",t,
			       std::real(ppconbil[t]),std::imag(ppconbil[t]),
			       std::real(ppnew(0,t)),std::imag(ppnew(0,t)),
			       reldiff(ppconbil[t].real(),ppnew(0,t).real()), reldiff(ppconbil[t].imag(),ppnew(0,t).imag()));
	fail = 1.;
      }else if(!UniqueID()) printf("PASS t=%d, old=(%g,%g) new=(%g,%g)\n",t,std::real(ppconbil[t]),std::imag(ppconbil[t]),std::real(ppnew(0,t)),std::imag(ppnew(0,t)) );
    }
    glb_sum(&fail);
    if(fail){
      if(!UniqueID()) printf("Failed comparison of pseudoscalar flavor singlet\n");
      if(exit_on_error) exit(-1);
    }
  }





  //Test 
  for(int i=0;i<GJP.VolNodeSites();i++){
    SpinColorFlavorMatrix pm_calc, pm_flip;
    pw_prop_pminus.siteMatrix(pm_calc,i); //p- formed from p- f0 and p- f1
    pm_flip.generate_from_cconj_pair(QPropWcontainer::verify_convert(prop_f0_pminus,"",""), QPropWcontainer::verify_convert(prop_f0_pplus,"",""), lattice,i);

    double fail = 0.0;

    for(int scf=0;scf<24*24;scf++){
      int rem = scf;
      int s1 = rem % 4; rem /= 4;
      int c1 = rem % 3; rem /= 3;
      int f1 = rem % 2; rem /= 2;
      int s2 = rem % 4; rem /= 4;
      int c2 = rem % 3; rem /= 3;
      int f2 = rem % 2; rem /= 2;
      
      const Complex &v_calc = pm_calc(s1,c1,f1,s2,c2,f2);
      const Complex &v_flip = pm_flip(s1,c1,f1,s2,c2,f2);
      if(fabs(std::real(v_calc - v_flip))>1e-12 ||
	 fabs(std::imag(v_calc - v_flip))>1e-12){
	if(!UniqueID()) printf("cconj reln test fail %d (%d %d %d %d %d %d): (%g, %g) vs (%g, %g)\n",i,s1,c1,f2,s2,c2,f2,std::real(v_calc),std::imag(v_calc),std::real(v_flip),std::imag(v_flip));
	fail = 1.;
      }
    }
    glb_sum(&fail);
    if(fail != 0.0){
      if(!UniqueID()) printf("cconj reln test failed\n");
      if(exit_on_error) exit(-1);
    }
  }
  if(!UniqueID()) printf("Passed cconj reln test 2 (-p)\n");

  PropWrapper pw_prop_pminus_fromflip(&QPropWcontainer::verify_convert(prop_f0_pplus,"","").getProp(lattice),
				      &QPropWcontainer::verify_convert(prop_f1_pplus,"","").getProp(lattice), true);

  //Test the mom flip
  for(int i=0;i<GJP.VolNodeSites();i++){
    SpinColorFlavorMatrix pm_calc, pm_flip;
    pw_prop_pminus.siteMatrix(pm_calc,i);
    pw_prop_pminus_fromflip.siteMatrix(pm_flip,i);

    double fail = 0.0;

    for(int scf=0;scf<24*24;scf++){
      int rem = scf;
      int s1 = rem % 4; rem /= 4;
      int c1 = rem % 3; rem /= 3;
      int f1 = rem % 2; rem /= 2;
      int s2 = rem % 4; rem /= 4;
      int c2 = rem % 3; rem /= 3;
      int f2 = rem % 2; rem /= 2;
      
      const Complex &v_calc = pm_calc(s1,c1,f1,s2,c2,f2);
      const Complex &v_flip = pm_flip(s1,c1,f1,s2,c2,f2);
      if(fabs(std::real(v_calc - v_flip))>1e-12 ||
	 fabs(std::imag(v_calc - v_flip))>1e-12){
	if(!UniqueID()) printf("Flip test fail %d (%d %d %d %d %d %d): (%g, %g) vs (%g, %g)\n",i,s1,c1,f2,s2,c2,f2,std::real(v_calc),std::imag(v_calc),std::real(v_flip),std::imag(v_flip));
	fail = 1.;
      }
    }
    glb_sum(&fail);
    if(fail != 0.0){
      if(!UniqueID()) printf("Flip test failed\n");
      if(exit_on_error) exit(-1);
    }
  }
  if(!UniqueID()) printf("Passed mom flip test\n");


  int Lt = L[3];

  //Kaon stuff
  //Assign momentum -p to the strange quark. This is the one that is daggered so we need +p for the propagator
  PropagatorContainer &propH_f0_pplus = computePropagatorOld("propH_f0_pplus",0.04,prec,tsrc,0,p.ptr(),BND_CND_APRD,lattice, "propH_f1_pplus");
  PropagatorContainer &propH_f1_pplus = computePropagatorOld("propH_f1_pplus",0.04,prec,tsrc,1,p.ptr(),BND_CND_APRD,lattice, "propH_f0_pplus");

  PropWrapper pw_propH_pplus(&QPropWcontainer::verify_convert(propH_f0_pplus,"","").getProp(lattice),
			     &QPropWcontainer::verify_convert(propH_f1_pplus,"","").getProp(lattice));
  PropSiteMatrixStandard pw_propH_pplus_getter(pw_propH_pplus, BND_CND_APRD, tsrc);

  //Test BK against old code
  //Generate props for sink kaon
  int tsnk = Lt-1;
  PropagatorContainer &propH_f0_pplus_tsnk = computePropagatorOld("propH_f0_pplus_tsnk",0.04,prec,tsnk,0,p.ptr(),BND_CND_APRD,lattice, "propH_f1_pplus_tsnk");
  PropagatorContainer &propH_f1_pplus_tsnk = computePropagatorOld("propH_f1_pplus_tsnk",0.04,prec,tsnk,1,p.ptr(),BND_CND_APRD,lattice, "propH_f0_pplus_tsnk");

  PropWrapper pw_propH_pplus_tsnk(&QPropWcontainer::verify_convert(propH_f0_pplus_tsnk,"","").getProp(lattice),
				  &QPropWcontainer::verify_convert(propH_f1_pplus_tsnk,"","").getProp(lattice));
  PropSiteMatrixStandard pw_propH_pplus_tsnk_getter(pw_propH_pplus_tsnk, BND_CND_APRD, tsnk);

  PropagatorContainer &prop_f0_pplus_tsnk = computePropagatorOld("prop_f0_pplus_tsnk",0.01,prec,tsnk,0,p.ptr(),BND_CND_APRD,lattice, "prop_f1_pplus_tsnk");
  PropagatorContainer &prop_f1_pplus_tsnk = computePropagatorOld("prop_f1_pplus_tsnk",0.01,prec,tsnk,1,p.ptr(),BND_CND_APRD,lattice, "prop_f0_pplus_tsnk");

  PropWrapper pw_prop_pplus_tsnk(&QPropWcontainer::verify_convert(prop_f0_pplus_tsnk,"","").getProp(lattice),
				 &QPropWcontainer::verify_convert(prop_f1_pplus_tsnk,"","").getProp(lattice));
  PropSiteMatrixStandard pw_prop_pplus_tsnk_getter(pw_prop_pplus_tsnk, BND_CND_APRD, tsnk);

  {
    fMatrix<Rcomplex> bk_new(1,Lt);
    bool do_flavor_project = false;

    ThreeMomentum p_psi_h_t0 = mp; //- the momentum of the daggered prop
    ThreeMomentum p_psi_h_t1 = mp; //- the momentum of the daggered prop

    GparityBK(bk_new, tsrc, tsnk,
	      pw_propH_pplus_getter, pw_prop_pplus_getter, p_psi_h_t0,
	      pw_propH_pplus_tsnk_getter, pw_prop_pplus_tsnk_getter, p_psi_h_t1,
	      do_flavor_project);
    
    ContractionTypeOVVpAA old_args;
    old_args.prop_H_t0 = strdup("propH_f0_pplus");
    old_args.prop_L_t0 = strdup("prop_f0_pplus");
    old_args.prop_H_t1 = strdup("propH_f0_pplus_tsnk");
    old_args.prop_L_t1 = strdup("prop_f0_pplus_tsnk");
    AlgGparityContract gpcon(lattice, common_arg);

    CorrelationFunction bk_old("BK",CorrelationFunction::THREADED);
    gpcon.contract_OVVpAA_gparity(bk_old, old_args);

    std::vector<Rcomplex> bk_old2(Lt); //sum contractions and fix normalization (old contractions had coeff of 1/2 when it should have been 2)
    for(int t=0;t<Lt;t++){
      bk_old2[t] = 4.*bk_old(0,t);
      for(int i=1;i<4;i++) bk_old2[t] += 4.*bk_old(i,t);
    }

    double fail = 0;
    if(!UniqueID()) printf("BK\n");
    for(int t=0;t<L[3];t++){
      if(!test_equals(bk_old2[t],bk_new(0,t),tol)){
	if(!UniqueID()) printf("FAIL t=%d, old=(%g,%g) new=(%g,%g), reldiff=(%g,%g)\n",t,
			       std::real(bk_old2[t]),std::imag(bk_old2[t]),
			       std::real(bk_new(0,t)),std::imag(bk_new(0,t)),
			       reldiff(bk_old2[t].real(),bk_new(0,t).real()), reldiff(bk_old2[t].imag(),bk_new(0,t).imag()) );
	fail = 1.;
      }else if(!UniqueID()) printf("PASS t=%d, old=(%g,%g) new=(%g,%g)\n",t,std::real(bk_old2[t]),std::imag(bk_old2[t]),std::real(bk_new(0,t)),std::imag(bk_new(0,t)) );
    }
    glb_sum(&fail);
    if(fail){
      if(!UniqueID()) printf("Failed comparison of BK\n");
      if(exit_on_error) exit(-1);
    }
  }

  //Check J5 and J5q contractions
  {
    ContractionTypeMres old_args;
    old_args.prop = strdup("prop_f0_pplus");
    
    CorrelationFunction pion_old("J5",1,CorrelationFunction::THREADED), j5_q_old("J5q",1,CorrelationFunction::THREADED);
    AlgGparityContract gpcon(lattice, common_arg);
    gpcon.measure_mres_gparity(old_args, pion_old, j5_q_old);
    
    std::vector<Rcomplex> pion_old2(Lt), j5_q_old2(Lt);
    for(int t=0;t<Lt;t++){
      pion_old2[t] = 0.5*pion_old(0,t); //new normalization
      j5_q_old2[t] = 0.5*j5_q_old(0,t);
    }

    fMatrix<Rcomplex> pion_new(1,Lt), j5_q_new(1,Lt);

    ThreeMomentum p_psibar = p;
    ThreeMomentum p_psi = p;

    J5Gparity(pion_new,tsrc,p_psibar,p_psi,pw_prop_pminus_getter,pw_prop_pplus_getter,SPLANE_BOUNDARY,false); //disable flavor project to compare with old code
    J5Gparity(j5_q_new,tsrc,p_psibar,p_psi,pw_prop_pminus_getter,pw_prop_pplus_getter,SPLANE_MIDPOINT,false);

    double fail = 0;
    if(!UniqueID()) printf("J5\n");
    for(int t=0;t<L[3];t++){
      if(!test_equals(pion_old2[t],pion_new(0,t),tol)){
	if(!UniqueID()) printf("FAIL t=%d, old=(%g,%g) new=(%g,%g), reldiff=(%g,%g)\n",t,
			       std::real(pion_old2[t]),std::imag(pion_old2[t]),
			       std::real(pion_new(0,t)),std::imag(pion_new(0,t)),
			       reldiff(pion_old2[t].real(),pion_new(0,t).real()), reldiff(pion_old2[t].imag(),pion_new(0,t).imag()));
	fail = 1.;
      }else if(!UniqueID()) printf("PASS t=%d, old=(%g,%g) new=(%g,%g)\n",t,std::real(pion_old2[t]),std::imag(pion_old2[t]),std::real(pion_new(0,t)),std::imag(pion_new(0,t)) );
    }
    glb_sum(&fail);
    if(fail){
      if(!UniqueID()) printf("Failed comparison of J5\n");
      if(exit_on_error) exit(-1);
    }

    if(!UniqueID()) printf("J5q\n");
    for(int t=0;t<L[3];t++){
      if(!test_equals(j5_q_old2[t],j5_q_new(0,t),tol)){
	if(!UniqueID()) printf("FAIL t=%d, old=(%g,%g) new=(%g,%g), reldiff=(%g,%g)\n",t,
			       std::real(j5_q_old2[t]),std::imag(j5_q_old2[t]),
			       std::real(j5_q_new(0,t)),std::imag(j5_q_new(0,t)),
			       reldiff(j5_q_old2[t].real(),j5_q_new(0,t).real()), reldiff(j5_q_old2[t].imag(),j5_q_new(0,t).imag()) );
	fail = 1.;
      }else if(!UniqueID()) printf("PASS t=%d, old=(%g,%g) new=(%g,%g)\n",t,std::real(j5_q_old2[t]),std::imag(j5_q_old2[t]),std::real(j5_q_new(0,t)),std::imag(j5_q_new(0,t)) );
    }
    glb_sum(&fail);
    if(fail){
      if(!UniqueID()) printf("Failed comparison of J5q\n");
      if(exit_on_error) exit(-1);
    }
  }


  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  return 0;
}





CPS_END_NAMESPACE

#endif
