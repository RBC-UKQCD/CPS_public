#ifndef _BENCHMARK_MESONFIELD_H
#define _BENCHMARK_MESONFIELD_H

CPS_START_NAMESPACE

inline int toInt(const char* a){
  std::stringstream ss; ss << a; int o; ss >> o;
  return o;
}

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



void randomMatrix(SpinColorFlavorMatrix &A, CPSspinColorFlavorMatrix<cps::ComplexD> &B){
  for(int s1=0;s1<4;s1++)
    for(int s2=0;s2<4;s2++)
      for(int c1=0;c1<3;c1++)
	for(int c2=0;c2<3;c2++)
	  for(int f1=0;f1<2;f1++)
	    for(int f2=0;f2<2;f2++){
	      cps::ComplexD tmp;
	      _testRandom<cps::ComplexD>::rand(&tmp,1, 3.0, -3.0);
	      A(s1,c1,f1,s2,c2,f2) = tmp;
	      B(s1,s2)(c1,c2)(f1,f2) = tmp;
	    }
}
void randomMatrix(SpinColorFlavorMatrix &A){
  for(int s1=0;s1<4;s1++)
    for(int s2=0;s2<4;s2++)
      for(int c1=0;c1<3;c1++)
	for(int c2=0;c2<3;c2++)
	  for(int f1=0;f1<2;f1++)
	    for(int f2=0;f2<2;f2++){
	      cps::ComplexD tmp;
	      _testRandom<cps::ComplexD>::rand(&tmp,1, 3.0, -3.0);
	      A(s1,c1,f1,s2,c2,f2) = tmp;
	    }
}
void randomMatrix(CPSspinColorFlavorMatrix<cps::ComplexD> &B){
  for(int s1=0;s1<4;s1++)
    for(int s2=0;s2<4;s2++)
      for(int c1=0;c1<3;c1++)
	for(int c2=0;c2<3;c2++)
	  for(int f1=0;f1<2;f1++)
	    for(int f2=0;f2<2;f2++){
	      cps::ComplexD tmp;
	      _testRandom<cps::ComplexD>::rand(&tmp,1, 3.0, -3.0);
	      B(s1,s2)(c1,c2)(f1,f2) = tmp;
	    }
}


void benchmarkTrace(const int ntests, const double tol){
  typedef CPSsquareMatrix<CPSsquareMatrix<CPSsquareMatrix<cps::ComplexD,2>,3>,4> SCFmat;
 typedef CPSsquareMatrix<cps::ComplexD,3> Cmat;
    
  //Test they give the same answer
  {
    SpinColorFlavorMatrix old_mat;
    CPSspinColorFlavorMatrix<cps::ComplexD> new_mat;
    randomMatrix(old_mat,new_mat);

    cps::Complex cp =  old_mat.Trace();
    cps::Complex gd = new_mat.Trace();
      
    bool fail = false;
    double rdiff = fabs(gd.real()-cp.real());
    double idiff = fabs(gd.imag()-cp.imag());
    if(rdiff > tol|| idiff > tol){
      printf("Fail: Trace Grid (%g,%g) CPS (%g,%g) Diff (%g,%g)\n",gd.real(),gd.imag(), cp.real(),cp.imag(), cp.real()-gd.real(), cp.imag()-gd.imag());
      fail = true;
    }
    if(fail) ERR.General("","","Trace test failed\n");
    else printf("Trace pass\n");
  }
  {
    //Benchmark it    
    SpinColorFlavorMatrix old_mat;
    CPSspinColorFlavorMatrix<cps::ComplexD> new_mat;
    
    //SpinFlavorTrace of SpinColorFlavorMatrix
    Float total_time_old = 0.;
    cps::ComplexD tmp_old;
    for(int iter=0;iter<ntests;iter++){
      randomMatrix(old_mat);      
      total_time_old -= dclock();
      tmp_old = old_mat.Trace();
      total_time_old += dclock();
    }
    Float total_time_new = 0.;
    cps::ComplexD tmp_new;
    for(int iter=0;iter<ntests;iter++){
      randomMatrix(new_mat);
      total_time_new -= dclock();
      tmp_new = new_mat.Trace();
      total_time_new += dclock();
    }
    printf("Trace: Avg time new code %d iters: %g secs\n",ntests,total_time_new/ntests);
    printf("Trace: Avg time old code %d iters: %g secs\n",ntests,total_time_old/ntests);
  }
}



void benchmarkSpinFlavorTrace(const int ntests, const double tol){
  typedef CPSsquareMatrix<CPSsquareMatrix<CPSsquareMatrix<cps::ComplexD,2>,3>,4> SCFmat;
  typedef CPSsquareMatrix<cps::ComplexD,3> Cmat;
    
  //Test they give the same answer
  {
    SpinColorFlavorMatrix old_mat;
    CPSspinColorFlavorMatrix<cps::ComplexD> new_mat;
    randomMatrix(old_mat,new_mat);

    Matrix tmp_mat_old =  old_mat.SpinFlavorTrace();
    Cmat tmp_mat_new = new_mat.SpinFlavorTrace();//TraceTwoIndices<0,2>();
      
    bool fail = false;
    for(int c1=0;c1<3;c1++)
      for(int c2=0;c2<3;c2++){
	cps::ComplexD gd = tmp_mat_old(c1,c2);
	cps::ComplexD cp = tmp_mat_new(c1,c2);
	  
	double rdiff = fabs(gd.real()-cp.real());
	double idiff = fabs(gd.imag()-cp.imag());
	if(rdiff > tol|| idiff > tol){
	  printf("Fail: SFtrace Grid (%g,%g) CPS (%g,%g) Diff (%g,%g)\n",gd.real(),gd.imag(), cp.real(),cp.imag(), cp.real()-gd.real(), cp.imag()-gd.imag());
	  fail = true;
	}
      }
    if(fail) ERR.General("","","SFtrace test failed\n");
    else printf("SFtrace pass\n");
  }
  {
    //Benchmark it    
    SpinColorFlavorMatrix old_mat;
    CPSspinColorFlavorMatrix<cps::ComplexD> new_mat;
    
    //SpinFlavorTrace of SpinColorFlavorMatrix
    Float total_time_old = 0.;
    Matrix tmp_mat_old;
    for(int iter=0;iter<ntests;iter++){
      randomMatrix(old_mat);      
      total_time_old -= dclock();
      tmp_mat_old = old_mat.SpinFlavorTrace();
      total_time_old += dclock();
    }
    Float total_time_new = 0.;
    Cmat tmp_mat_new;
    for(int iter=0;iter<ntests;iter++){
      randomMatrix(new_mat);
      total_time_new -= dclock();
      //tmp_mat_new.zero();
      //_PartialDoubleTraceImpl<Cmat,CPSspinColorFlavorMatrix<cps::ComplexD>,0,2>::doit(tmp_mat_new,new_mats[iter]);
      tmp_mat_new = new_mat.TraceTwoIndices<0,2>();
      total_time_new += dclock();
    }
    printf("SFtrace: Avg time new code %d iters: %g secs\n",ntests,total_time_new/ntests);
    printf("SFtrace: Avg time old code %d iters: %g secs\n",ntests,total_time_old/ntests);
  }
}


void benchmarkTraceProd(const int ntests, const double tol){
  typedef CPSsquareMatrix<CPSsquareMatrix<CPSsquareMatrix<cps::ComplexD,2>,3>,4> SCFmat;
  typedef CPSsquareMatrix<cps::ComplexD,3> Cmat;
    
  //Test they give the same answer
  {
    SpinColorFlavorMatrix old_mat1, old_mat2;
    CPSspinColorFlavorMatrix<cps::ComplexD> new_mat1, new_mat2;
    randomMatrix(old_mat1,new_mat1);
    randomMatrix(old_mat2,new_mat2);

    cps::ComplexD tr_old = Trace(old_mat1,old_mat2);
    cps::ComplexD tr_new = Trace(new_mat1,new_mat2);
      
    bool fail = false;
    cps::ComplexD &gd = tr_new;
    cps::ComplexD &cp = tr_old;
	  
    double rdiff = fabs(gd.real()-cp.real());
    double idiff = fabs(gd.imag()-cp.imag());
    if(rdiff > tol|| idiff > tol){
      printf("Fail: Prodtrace Grid (%g,%g) CPS (%g,%g) Diff (%g,%g)\n",gd.real(),gd.imag(), cp.real(),cp.imag(), cp.real()-gd.real(), cp.imag()-gd.imag());
      fail = true;
    }
    if(fail) ERR.General("","","Prodtrace test failed\n");
    else printf("Prodtrace pass\n");
  }
  //Benchmark it    
  {
    SpinColorFlavorMatrix old_mat1, old_mat2;
    CPSspinColorFlavorMatrix<cps::ComplexD> new_mat1, new_mat2;
    
    Float total_time_old = 0.;
    cps::ComplexD tr_old;
    for(int iter=0;iter<ntests;iter++){
      randomMatrix(old_mat1);
      randomMatrix(old_mat2);
      total_time_old -= dclock();
      tr_old = Trace(old_mat1,old_mat2);
      total_time_old += dclock();
    }
    Float total_time_new = 0.;
    cps::ComplexD tr_new;
    for(int iter=0;iter<ntests;iter++){
      randomMatrix(new_mat1);
      randomMatrix(new_mat2);
      total_time_new -= dclock();
      tr_new = Trace(new_mat1,new_mat2);
      total_time_new += dclock();
    }
    printf("Prodtrace: Avg time new code %d iters: %g secs\n",ntests,total_time_new/ntests);
    printf("Prodtrace: Avg time old code %d iters: %g secs\n",ntests,total_time_old/ntests);
  }

}


void benchmarkColorTranspose(const int ntests, const double tol){
  typedef CPSsquareMatrix<CPSsquareMatrix<CPSsquareMatrix<cps::ComplexD,2>,3>,4> SCFmat;
  typedef CPSsquareMatrix<cps::ComplexD,3> Cmat;
    
  //Test they give the same answer
  {
    SpinColorFlavorMatrix old_mat;
    CPSspinColorFlavorMatrix<cps::ComplexD> new_mat;
    randomMatrix(old_mat,new_mat);

    old_mat.transpose_color();
    new_mat.TransposeOnIndex<1>();
      
    bool fail = false;
    for(int s1=0;s1<4;s1++)
      for(int s2=0;s2<4;s2++)	
	for(int c1=0;c1<3;c1++)
	  for(int c2=0;c2<3;c2++)
	    for(int f1=0;f1<2;f1++)
	      for(int f2=0;f2<2;f2++){
		cps::ComplexD gd = new_mat(s1,s2)(c1,c2)(f1,f2);
		cps::ComplexD cp = old_mat(s1,c1,f1,s2,c2,f2);
	  
		double rdiff = fabs(gd.real()-cp.real());
		double idiff = fabs(gd.imag()-cp.imag());
		if(rdiff > tol|| idiff > tol){
		  printf("Fail: colortranspose Grid (%g,%g) CPS (%g,%g) Diff (%g,%g)\n",gd.real(),gd.imag(), cp.real(),cp.imag(), cp.real()-gd.real(), cp.imag()-gd.imag());
		  fail = true;
		}
	      }
    if(fail) ERR.General("","","colortranspose test failed\n");
    else printf("colortranspose pass\n");
  }
  {
    //Benchmark it    
    SpinColorFlavorMatrix old_mat;
    CPSspinColorFlavorMatrix<cps::ComplexD> new_mat;
    
    Float total_time_old = 0.;

    for(int iter=0;iter<ntests;iter++){
      randomMatrix(old_mat);      
      total_time_old -= dclock();
      old_mat.transpose_color();
      total_time_old += dclock();
    }
    Float total_time_new = 0.;

    for(int iter=0;iter<ntests;iter++){
      randomMatrix(new_mat);
      total_time_new -= dclock();
      new_mat.TransposeOnIndex<1>();
      total_time_new += dclock();
    }
    printf("colorTranspose: Avg time new code %d iters: %g secs\n",ntests,total_time_new/ntests);
    printf("colorTranspose: Avg time old code %d iters: %g secs\n",ntests,total_time_old/ntests);
  }
}

void multGammaLeftOld(SpinColorFlavorMatrix &M, const int whichGamma, const int i, const int mu){
  assert(whichGamma == 1 || whichGamma==2);
  static int g1[8] = {0,1,0,1,2,3,2,3};
  static int g2[8] = {1,0,3,2,1,0,3,2};

  int gg = whichGamma == 1 ? g1[i] : g2[i];
  switch(gg){
  case 0:
    M.pl(F0).gl(mu);
    break;
  case 1:
    M.pl(F0).glAx(mu);
    break;
  case 2:
    M.pl(F1).gl(mu); M *= -1.0;
    break;
  case 3:
    M.pl(F1).glAx(mu); M *= -1.0;
    break;
  default:
    ERR.General("ComputeKtoPiPiGparityBase","multGammaLeft","Invalid idx\n");
    break;
  }
}

void benchmarkmultGammaLeft(const int ntests, const double tol){
  typedef CPSsquareMatrix<CPSsquareMatrix<CPSsquareMatrix<cps::ComplexD,2>,3>,4> SCFmat;
  
  //Test they give the same answer
  for(int gamma=1;gamma<=2;gamma++)
    for(int i=0;i<8;i++)
      for(int mu=0;mu<4;mu++){
	SpinColorFlavorMatrix old_mat;
	CPSspinColorFlavorMatrix<cps::ComplexD> new_mat;
	randomMatrix(old_mat,new_mat);

	multGammaLeftOld(old_mat,gamma,i,mu);
	ComputeKtoPiPiGparityBase::multGammaLeft(new_mat,gamma,i,mu);
      
	bool fail = false;
	for(int s1=0;s1<4;s1++)
	  for(int s2=0;s2<4;s2++)	
	    for(int c1=0;c1<3;c1++)
	      for(int c2=0;c2<3;c2++)
		for(int f1=0;f1<2;f1++)
		  for(int f2=0;f2<2;f2++){
		    cps::ComplexD gd = new_mat(s1,s2)(c1,c2)(f1,f2);
		    cps::ComplexD cp = old_mat(s1,c1,f1,s2,c2,f2);
		    
		    double rdiff = fabs(gd.real()-cp.real());
		    double idiff = fabs(gd.imag()-cp.imag());
		    if(rdiff > tol|| idiff > tol){
		      printf("Fail: multGammaLeft(%d,%d,%d) Grid (%g,%g) CPS (%g,%g) Diff (%g,%g)\n",gamma,i,mu,gd.real(),gd.imag(), cp.real(),cp.imag(), cp.real()-gd.real(), cp.imag()-gd.imag());
		      fail = true;
		    }
		  }
	if(fail) ERR.General("","","multGammaLeft test %d %d %d failed\n",gamma,i,mu);
	//else printf("multGammaLeft %d %d %d pass\n",gamma,i,mu);
      }

  //Benchmark
  for(int gamma=1;gamma<=2;gamma++)
    for(int i=0;i<8;i++)
      for(int mu=0;mu<4;mu++){
	SpinColorFlavorMatrix old_mat;
	CPSspinColorFlavorMatrix<cps::ComplexD> new_mat;

	Float total_time_old = 0.;
	
	for(int iter=0;iter<ntests;iter++){
	  randomMatrix(old_mat);      
	  total_time_old -= dclock();
	  multGammaLeftOld(old_mat,gamma,i,mu);
	  total_time_old += dclock();
	}
	Float total_time_new = 0.;
	
	for(int iter=0;iter<ntests;iter++){
	  randomMatrix(new_mat);
	  total_time_new -= dclock();
	  ComputeKtoPiPiGparityBase::multGammaLeft(new_mat,gamma,i,mu);
	  total_time_new += dclock();
	}
	printf("multGammaLeft %d %d %d: Avg time new code %d iters: %g secs\n",gamma,i,mu,ntests,total_time_new/ntests);
	printf("multGammaLeft %d %d %d: Avg time old code %d iters: %g secs\n\n",gamma,i,mu,ntests,total_time_old/ntests);
      }  
}

template<typename mf_Complex, typename grid_Complex>
bool compare(const CPSspinColorFlavorMatrix<mf_Complex> &orig, const CPSspinColorFlavorMatrix<grid_Complex> &grid, const double tol){
  bool fail = false;
  
  mf_Complex gd;
  for(int sl=0;sl<4;sl++)
    for(int cl=0;cl<3;cl++)
      for(int fl=0;fl<2;fl++)
	for(int sr=0;sr<4;sr++)
	  for(int cr=0;cr<3;cr++)
	    for(int fr=0;fr<2;fr++){
	      gd = Reduce( grid(sl,sr)(cl,cr)(fl,fr) );
	      const mf_Complex &cp = orig(sl,sr)(cl,cr)(fl,fr);
	      
	      double rdiff = fabs(gd.real()-cp.real());
	      double idiff = fabs(gd.imag()-cp.imag());
	      if(rdiff > tol|| idiff > tol){
		printf("Fail: Grid (%g,%g) CPS (%g,%g) Diff (%g,%g)\n",gd.real(),gd.imag(), cp.real(),cp.imag(), cp.real()-gd.real(), cp.imag()-gd.imag());
		fail = true;
	      }
	    }
  return !fail;
}

template<typename mf_Complex>
bool compare(const CPSspinColorFlavorMatrix<mf_Complex> &orig, const CPSspinColorFlavorMatrix<mf_Complex> &newimpl, const double tol){
  bool fail = false;
  
  for(int sl=0;sl<4;sl++)
    for(int cl=0;cl<3;cl++)
      for(int fl=0;fl<2;fl++)
	for(int sr=0;sr<4;sr++)
	  for(int cr=0;cr<3;cr++)
	    for(int fr=0;fr<2;fr++){
	      const mf_Complex &gd = newimpl(sl,sr)(cl,cr)(fl,fr);
	      const mf_Complex &cp = orig(sl,sr)(cl,cr)(fl,fr);
	      
	      double rdiff = fabs(gd.real()-cp.real());
	      double idiff = fabs(gd.imag()-cp.imag());
	      if(rdiff > tol|| idiff > tol){
		printf("Fail: Newimpl (%g,%g) Orig (%g,%g) Diff (%g,%g)\n",gd.real(),gd.imag(), cp.real(),cp.imag(), cp.real()-gd.real(), cp.imag()-gd.imag());
		fail = true;
	      }
	    }
  return !fail;
}


template< typename mf_Complex, int SiteSize, typename FlavorPolicy, typename AllocPolicy>
void printRow(const CPSfield<mf_Complex,SiteSize,FourDpolicy,FlavorPolicy,AllocPolicy> &field, const int dir, const std::string &comment,
	       typename my_enable_if< _equal<typename ComplexClassify<mf_Complex>::type, complex_double_or_float_mark>::value, const int>::type = 0
	       ){
  int L = GJP.Nodes(dir)*GJP.NodeSites(dir);
  std::vector<mf_Complex> buf(L,0.);

  int other_dirs[3]; int aa=0;
  for(int i=0;i<4;i++)
    if(i!=dir) other_dirs[aa++] = i;

  
  if(GJP.NodeCoor(other_dirs[0]) == 0 && GJP.NodeCoor(other_dirs[1]) == 0 && GJP.NodeCoor(other_dirs[2]) == 0){
    for(int x=GJP.NodeCoor(dir)*GJP.NodeSites(dir); x < (GJP.NodeCoor(dir)+1)*GJP.NodeSites(dir); x++){
      int lcoor[4] = {0,0,0,0};
      lcoor[dir] = x - GJP.NodeCoor(dir)*GJP.NodeSites(dir);
      
      mf_Complex const* site_ptr = field.site_ptr(lcoor);
      buf[x] = *site_ptr;
    }
  }
  globalSumComplex(buf.data(),L);

  
  if(!UniqueID()){
    printf("%s: (",comment.c_str()); fflush(stdout);
    for(int x=0;x<L;x++){
      if(x % GJP.NodeSites(dir) == 0 && x!=0)
	printf(")(");
      
      printf("[%f,%f] ",buf[x].real(),buf[x].imag());
    }
    printf(")\n"); fflush(stdout);
  }
}

#ifdef USE_GRID
template< typename mf_Complex, int SiteSize, typename FlavorPolicy, typename AllocPolicy>
void printRow(const CPSfield<mf_Complex,SiteSize,FourDSIMDPolicy,FlavorPolicy,AllocPolicy> &field, const int dir, const std::string &comment,
	       typename my_enable_if< _equal<typename ComplexClassify<mf_Complex>::type, grid_vector_complex_mark>::value, const int>::type = 0
	       ){
  typedef typename mf_Complex::scalar_type ScalarComplex;
  NullObject null_obj;
  CPSfield<ScalarComplex,SiteSize,FourDpolicy,FlavorPolicy,StandardAllocPolicy> tmp(null_obj);
  tmp.importField(field);
  printRow(tmp,dir,comment);
}
#endif


void testCyclicPermute(){
  NullObject null_obj;
  {//4D
    CPSfield<cps::ComplexD,1,FourDpolicy,FixedFlavorPolicy<1>,StandardAllocPolicy> from(null_obj);
    CPSfield<cps::ComplexD,1,FourDpolicy,FixedFlavorPolicy<1>,StandardAllocPolicy> tmp1(null_obj);
    CPSfield<cps::ComplexD,1,FourDpolicy,FixedFlavorPolicy<1>,StandardAllocPolicy> tmp2(null_obj);

    from.testRandom();

    for(int dir=0;dir<4;dir++){
      for(int pm=-1;pm<=1;pm+=2){
	if(!UniqueID()) printf("Testing 4D permute in direction %c%d\n",pm == 1 ? '+' : '-',dir);
	//permute in incr until we cycle all the way around
	tmp1 = from;
	CPSfield<cps::ComplexD,1,FourDpolicy,FixedFlavorPolicy<1>,StandardAllocPolicy> *send = &tmp1;
	CPSfield<cps::ComplexD,1,FourDpolicy,FixedFlavorPolicy<1>,StandardAllocPolicy> *recv = &tmp2;

	int shifted = 0;
	printRow(from,dir,"Initial line      ");

	int total = GJP.Nodes(dir)*GJP.NodeSites(dir);
	int incr = GJP.NodeSites(dir)/2;
	int perm = 0;
	while(shifted < total){
	  cyclicPermute(*recv,*send,dir,pm,incr);
	  shifted += incr;
	  std::ostringstream comment; comment << "After perm " << perm++ << " by incr " << incr;
	  printRow(*recv,dir,comment.str());
	        
	  if(shifted < total)
	    std::swap(send,recv);
	}
	printRow(*recv,dir,"Final line      ");
      
	int coor[4];
	for(coor[0]=0;coor[0]<GJP.XnodeSites();coor[0]++){
	  for(coor[1]=0;coor[1]<GJP.YnodeSites();coor[1]++){
	    for(coor[2]=0;coor[2]<GJP.ZnodeSites();coor[2]++){
	      for(coor[3]=0;coor[3]<GJP.TnodeSites();coor[3]++){
		cps::ComplexD const* orig = from.site_ptr(coor);
		cps::ComplexD const* permd = recv->site_ptr(coor);
		if(orig->real() != permd->real() || orig->imag() != permd->imag()){
		  printf("Error node coor (%d,%d,%d,%d) (%d,%d,%d,%d) : (%g,%g) vs (%g,%g) diff (%g,%g)\n",GJP.XnodeCoor(),GJP.YnodeCoor(),GJP.ZnodeCoor(),GJP.TnodeCoor(),coor[0],coor[1],coor[2],coor[3],orig->real(),orig->imag(),permd->real(),permd->imag(), orig->real()-permd->real(),orig->imag()-permd->imag());
		}
	      }
	    }
	  }
	}
      }
    }
  }//End 4D

  {//3D
    CPSfield<cps::ComplexD,1,SpatialPolicy,FixedFlavorPolicy<1>,StandardAllocPolicy> from(null_obj);
    CPSfield<cps::ComplexD,1,SpatialPolicy,FixedFlavorPolicy<1>,StandardAllocPolicy> tmp1(null_obj);
    CPSfield<cps::ComplexD,1,SpatialPolicy,FixedFlavorPolicy<1>,StandardAllocPolicy> tmp2(null_obj);

    from.testRandom();

    for(int dir=0;dir<3;dir++){
      for(int pm=-1;pm<=1;pm+=2){
	if(!UniqueID()) printf("Testing 3D permute in direction %c%d\n",pm == 1 ? '+' : '-',dir);
	//permute in incr until we cycle all the way around
	tmp1 = from;
	CPSfield<cps::ComplexD,1,SpatialPolicy,FixedFlavorPolicy<1>,StandardAllocPolicy> *send = &tmp1;
	CPSfield<cps::ComplexD,1,SpatialPolicy,FixedFlavorPolicy<1>,StandardAllocPolicy> *recv = &tmp2;

	int shifted = 0;
	int total = GJP.Nodes(dir)*GJP.NodeSites(dir);
	int incr = GJP.NodeSites(dir)/2;
	int perm = 0;
	while(shifted < total){
	  cyclicPermute(*recv,*send,dir,pm,incr);
	  shifted += incr;
	  if(shifted < total)
	    std::swap(send,recv);
	}
      
	int coor[3];
	for(coor[0]=0;coor[0]<GJP.XnodeSites();coor[0]++){
	  for(coor[1]=0;coor[1]<GJP.YnodeSites();coor[1]++){
	    for(coor[2]=0;coor[2]<GJP.ZnodeSites();coor[2]++){
	      cps::ComplexD const* orig = from.site_ptr(coor);
	      cps::ComplexD const* permd = recv->site_ptr(coor);
	      if(orig->real() != permd->real() || orig->imag() != permd->imag()){
		printf("Error node coor (%d,%d,%d,%d) (%d,%d,%d) : (%g,%g) vs (%g,%g) diff (%g,%g)\n",GJP.XnodeCoor(),GJP.YnodeCoor(),GJP.ZnodeCoor(),GJP.TnodeCoor(),coor[0],coor[1],coor[2],orig->real(),orig->imag(),permd->real(),permd->imag(), orig->real()-permd->real(),orig->imag()-permd->imag());
	      }
	    }
	  }
	}	
      }
    }
  }//End 3D

#ifdef USE_GRID

  {//4D
    typedef FourDSIMDPolicy::ParamType simd_params;
    simd_params sp;
    FourDSIMDPolicy::SIMDdefaultLayout(sp, Grid::vComplexD::Nsimd() );
  
    CPSfield<Grid::vComplexD,1,FourDSIMDPolicy,FixedFlavorPolicy<1>,Aligned128AllocPolicy> from_grid(sp);
    CPSfield<Grid::vComplexD,1,FourDSIMDPolicy,FixedFlavorPolicy<1>,Aligned128AllocPolicy> tmp1_grid(sp);
    CPSfield<Grid::vComplexD,1,FourDSIMDPolicy,FixedFlavorPolicy<1>,Aligned128AllocPolicy> tmp2_grid(sp);

    CPSfield<cps::ComplexD,1,FourDpolicy,FixedFlavorPolicy<1>,StandardAllocPolicy> from(null_obj);
    CPSfield<cps::ComplexD,1,FourDpolicy,FixedFlavorPolicy<1>,StandardAllocPolicy> tmp1(null_obj);
    from.testRandom();
    from_grid.importField(from);

    for(int dir=0;dir<4;dir++){
      for(int pm=-1;pm<=1;pm+=2){
	if(!UniqueID()) printf("Testing 4D permute in direction %c%d with SIMD layout\n",pm == 1 ? '+' : '-',dir);
	//permute in incr until we cycle all the way around
	tmp1_grid = from_grid;
	CPSfield<Grid::vComplexD,1,FourDSIMDPolicy,FixedFlavorPolicy<1>,Aligned128AllocPolicy> *send = &tmp1_grid;
	CPSfield<Grid::vComplexD,1,FourDSIMDPolicy,FixedFlavorPolicy<1>,Aligned128AllocPolicy> *recv = &tmp2_grid;

	int shifted = 0;
	printRow(from_grid,dir,"Initial line      ");

	int total = GJP.Nodes(dir)*GJP.NodeSites(dir);
	int incr = GJP.NodeSites(dir)/2;
	int perm = 0;
	while(shifted < total){
	  cyclicPermute(*recv,*send,dir,pm,incr);
	  shifted += incr;
	  std::ostringstream comment; comment << "After perm " << perm++ << " by incr " << incr;
	  printRow(*recv,dir,comment.str());
	        
	  if(shifted < total)
	    std::swap(send,recv);
	}
	printRow(*recv,dir,"Final line      ");

	tmp1.importField(*recv);
      
	int coor[4];
	for(coor[0]=0;coor[0]<GJP.XnodeSites();coor[0]++){
	  for(coor[1]=0;coor[1]<GJP.YnodeSites();coor[1]++){
	    for(coor[2]=0;coor[2]<GJP.ZnodeSites();coor[2]++){
	      for(coor[3]=0;coor[3]<GJP.TnodeSites();coor[3]++){
		cps::ComplexD const* orig = from.site_ptr(coor);
		cps::ComplexD const* permd = tmp1.site_ptr(coor);
		if(orig->real() != permd->real() || orig->imag() != permd->imag()){
		  printf("Error node coor (%d,%d,%d,%d) (%d,%d,%d,%d) : (%g,%g) vs (%g,%g) diff (%g,%g)\n",GJP.XnodeCoor(),GJP.YnodeCoor(),GJP.ZnodeCoor(),GJP.TnodeCoor(),coor[0],coor[1],coor[2],coor[3],orig->real(),orig->imag(),permd->real(),permd->imag(), orig->real()-permd->real(),orig->imag()-permd->imag());
		}
	      }
	    }
	  }
	}
      }
    }
  }

  {//3D
    typedef ThreeDSIMDPolicy::ParamType simd_params;
    simd_params sp;
    ThreeDSIMDPolicy::SIMDdefaultLayout(sp, Grid::vComplexD::Nsimd() );
  
    CPSfield<Grid::vComplexD,1,ThreeDSIMDPolicy,FixedFlavorPolicy<1>,Aligned128AllocPolicy> from_grid(sp);
    CPSfield<Grid::vComplexD,1,ThreeDSIMDPolicy,FixedFlavorPolicy<1>,Aligned128AllocPolicy> tmp1_grid(sp);
    CPSfield<Grid::vComplexD,1,ThreeDSIMDPolicy,FixedFlavorPolicy<1>,Aligned128AllocPolicy> tmp2_grid(sp);

    CPSfield<cps::ComplexD,1,SpatialPolicy,FixedFlavorPolicy<1>,StandardAllocPolicy> from(null_obj);
    CPSfield<cps::ComplexD,1,SpatialPolicy,FixedFlavorPolicy<1>,StandardAllocPolicy> tmp1(null_obj);
    from.testRandom();
    from_grid.importField(from);

    for(int dir=0;dir<3;dir++){
      for(int pm=-1;pm<=1;pm+=2){
	if(!UniqueID()) printf("Testing 3D permute in direction %c%d with SIMD layout\n",pm == 1 ? '+' : '-',dir);
	//permute in incr until we cycle all the way around
	tmp1_grid = from_grid;
	CPSfield<Grid::vComplexD,1,ThreeDSIMDPolicy,FixedFlavorPolicy<1>,Aligned128AllocPolicy> *send = &tmp1_grid;
	CPSfield<Grid::vComplexD,1,ThreeDSIMDPolicy,FixedFlavorPolicy<1>,Aligned128AllocPolicy> *recv = &tmp2_grid;

	int shifted = 0;
	int total = GJP.Nodes(dir)*GJP.NodeSites(dir);
	int incr = GJP.NodeSites(dir)/2;
	int perm = 0;
	while(shifted < total){
	  cyclicPermute(*recv,*send,dir,pm,incr);
	  shifted += incr;
	  if(shifted < total)
	    std::swap(send,recv);
	}
	tmp1.importField(*recv);
      
	int coor[3];
	for(coor[0]=0;coor[0]<GJP.XnodeSites();coor[0]++){
	  for(coor[1]=0;coor[1]<GJP.YnodeSites();coor[1]++){
	    for(coor[2]=0;coor[2]<GJP.ZnodeSites();coor[2]++){
	      cps::ComplexD const* orig = from.site_ptr(coor);
	      cps::ComplexD const* permd = tmp1.site_ptr(coor);
	      if(orig->real() != permd->real() || orig->imag() != permd->imag()){
		printf("Error node coor (%d,%d,%d,%d) (%d,%d,%d) : (%g,%g) vs (%g,%g) diff (%g,%g)\n",GJP.XnodeCoor(),GJP.YnodeCoor(),GJP.ZnodeCoor(),GJP.TnodeCoor(),coor[0],coor[1],coor[2],orig->real(),orig->imag(),permd->real(),permd->imag(), orig->real()-permd->real(),orig->imag()-permd->imag());
	      }
	    }
	  }
	}
      }

    }
  }
#endif

  if(!UniqueID()){ printf("Passed permute test\n"); fflush(stdout); }
} 




template<typename A2Apolicies>
void demonstrateFFTreln(const A2AArg &a2a_args){
  //Demonstrate relation between FFTW fields
  A2AvectorW<A2Apolicies> W(a2a_args);
  A2AvectorV<A2Apolicies> V(a2a_args);
  W.testRandom();
  V.testRandom();

  int p1[3] = {1,1,1};
  int p5[3] = {5,1,1};

  twist<typename A2Apolicies::FermionFieldType> twist_p1(p1);
  twist<typename A2Apolicies::FermionFieldType> twist_p5(p5);
    
  A2AvectorVfftw<A2Apolicies> Vfftw_p1(a2a_args);
  Vfftw_p1.fft(V,&twist_p1);

  A2AvectorVfftw<A2Apolicies> Vfftw_p5(a2a_args);
  Vfftw_p5.fft(V,&twist_p5);

  //f5(n) = f1(n+1)
  for(int i=0;i<Vfftw_p1.getNmodes();i++)
    cyclicPermute(Vfftw_p1.getMode(i), Vfftw_p1.getMode(i), 0, -1, 1);
    
  printRow(Vfftw_p1.getMode(0),0, "T_-1 V(p1) T_-1");
  printRow(Vfftw_p5.getMode(0),0, "V(p5)          ");

  for(int i=0;i<Vfftw_p1.getNmodes();i++)
    assert( Vfftw_p1.getMode(i).equals( Vfftw_p5.getMode(i), 1e-7, true ) );

  A2AvectorWfftw<A2Apolicies> Wfftw_p1(a2a_args);
  Wfftw_p1.fft(W,&twist_p1);

  A2AvectorWfftw<A2Apolicies> Wfftw_p5(a2a_args);
  Wfftw_p5.fft(W,&twist_p5);

  for(int i=0;i<Wfftw_p1.getNmodes();i++)
    cyclicPermute(Wfftw_p1.getMode(i), Wfftw_p1.getMode(i), 0, -1, 1);

  printRow(Wfftw_p1.getMode(0),0, "T_-1 W(p1) T_-1");
  printRow(Wfftw_p5.getMode(0),0, "W(p5)          ");

  for(int i=0;i<Wfftw_p1.getNmodes();i++)
    assert( Wfftw_p1.getMode(i).equals( Wfftw_p5.getMode(i), 1e-7, true ) );

  if(!UniqueID()) printf("Passed FFT relation test\n");
}

template<typename ParamType, typename mf_Complex>
struct defaultFieldParams{
  static void get(ParamType &into){}
};

template<int N, typename mf_Complex>
struct defaultFieldParams< SIMDdims<N>, mf_Complex >{
  static void get(SIMDdims<N> &into){
    SIMDpolicyBase<N>::SIMDdefaultLayout(into, mf_Complex::Nsimd(), 2);
  }
};

template<typename A2Apolicies>
void testA2AvectorFFTrelnGparity(const A2AArg &a2a_args,Lattice &lat){
  assert(GJP.Gparity());

  if(lat.FixGaugeKind() == FIX_GAUGE_NONE){
    FixGaugeArg fix_gauge_arg;
    fix_gauge_arg.fix_gauge_kind = FIX_GAUGE_COULOMB_T;
    fix_gauge_arg.hyperplane_start = 0;
    fix_gauge_arg.hyperplane_step = 1;
    fix_gauge_arg.hyperplane_num = GJP.Tnodes()*GJP.TnodeSites();
    fix_gauge_arg.stop_cond = 1e-08;
    fix_gauge_arg.max_iter_num = 10000;

    CommonArg common_arg;
  
    AlgFixGauge fix_gauge(lat,&common_arg,&fix_gauge_arg);
    fix_gauge.run();
  }
  
  //Demonstrate relation between FFTW fields
  typedef typename A2Apolicies::ComplexType mf_Complex;
  typedef typename A2AvectorWfftw<A2Apolicies>::FieldInputParamType FieldInputParamType;
  FieldInputParamType fp; defaultFieldParams<FieldInputParamType, mf_Complex>::get(fp);
  
  A2AvectorW<A2Apolicies> W(a2a_args,fp);
  W.testRandom();

  int p_p1[3];
  GparityBaseMomentum(p_p1,+1);

  int p_m1[3];
  GparityBaseMomentum(p_m1,-1);

  //Perform base FFTs
  //twist<typename A2Apolicies::FermionFieldType> twist_p1(p_p1);
  //twist<typename A2Apolicies::FermionFieldType> twist_m1(p_m1);

  gaugeFixAndTwist<typename A2Apolicies::FermionFieldType> twist_p1(p_p1,lat);
  gaugeFixAndTwist<typename A2Apolicies::FermionFieldType> twist_m1(p_m1,lat);
  
  A2AvectorWfftw<A2Apolicies> Wfftw_p1(a2a_args,fp);
  Wfftw_p1.fft(W,&twist_p1);

  A2AvectorWfftw<A2Apolicies> Wfftw_m1(a2a_args,fp);
  Wfftw_m1.fft(W,&twist_m1);


  int p[3];  
  A2AvectorWfftw<A2Apolicies> result(a2a_args,fp);
  A2AvectorWfftw<A2Apolicies> compare(a2a_args,fp);
  
  //Get twist for first excited momentum in p1 set
  {
    memcpy(p,p_p1,3*sizeof(int));
    p[0] = 5;
    //twist<typename A2Apolicies::FermionFieldType> twist_p(p);
    gaugeFixAndTwist<typename A2Apolicies::FermionFieldType> twist_p(p,lat);    
    compare.fft(W,&twist_p);

    result.getTwistedFFT(p, &Wfftw_p1, &Wfftw_m1);

    if(!UniqueID()) printf("Testing p=(%d,%d,%d). Should require permute of 1 in x direction\n",p[0],p[1],p[2]);

    printRow(result.getMode(0),0,  "Result ");
    printRow(compare.getMode(0),0, "Compare");
    
    for(int i=0;i<compare.getNmodes();i++)
      assert( compare.getMode(i).equals( result.getMode(i), 1e-8, true ) );
  }
  //Get twist for first negative excited momentum in p1 set
  if(GJP.Bc(1) == BND_CND_GPARITY){
    memcpy(p,p_p1,3*sizeof(int));
    p[1] = -3;
    //twist<typename A2Apolicies::FermionFieldType> twist_p(p);
    gaugeFixAndTwist<typename A2Apolicies::FermionFieldType> twist_p(p,lat);    
    compare.fft(W,&twist_p);

    result.getTwistedFFT(p, &Wfftw_p1, &Wfftw_m1);

    if(!UniqueID()) printf("Testing p=(%d,%d,%d). Should require permute of 1 in y direction\n",p[0],p[1],p[2]);

    printRow(result.getMode(0),0,  "Result ");
    printRow(compare.getMode(0),0, "Compare");
    
    for(int i=0;i<compare.getNmodes();i++)
      assert( compare.getMode(i).equals( result.getMode(i), 1e-8, true ) );
  }
  //Try two directions
  if(GJP.Bc(1) == BND_CND_GPARITY){
    memcpy(p,p_p1,3*sizeof(int));
    p[0] = -3;
    p[1] = -3;
    //twist<typename A2Apolicies::FermionFieldType> twist_p(p);
    gaugeFixAndTwist<typename A2Apolicies::FermionFieldType> twist_p(p,lat);    
    compare.fft(W,&twist_p);

    result.getTwistedFFT(p, &Wfftw_p1, &Wfftw_m1);

    if(!UniqueID()) printf("Testing p=(%d,%d,%d). Should require permute of 1 in y direction\n",p[0],p[1],p[2]);

    printRow(result.getMode(0),0,  "Result ");
    printRow(compare.getMode(0),0, "Compare");
    
    for(int i=0;i<compare.getNmodes();i++)
      assert( compare.getMode(i).equals( result.getMode(i), 1e-7, true ) );
  }
  //Try 3 directions
  if(GJP.Bc(1) == BND_CND_GPARITY && GJP.Bc(2) == BND_CND_GPARITY){
    memcpy(p,p_p1,3*sizeof(int));
    p[0] = -3;
    p[1] = -3;
    p[2] = -3;
    //twist<typename A2Apolicies::FermionFieldType> twist_p(p);
    gaugeFixAndTwist<typename A2Apolicies::FermionFieldType> twist_p(p,lat);    
    compare.fft(W,&twist_p);

    result.getTwistedFFT(p, &Wfftw_p1, &Wfftw_m1);

    if(!UniqueID()) printf("Testing p=(%d,%d,%d). Should require permute of 1 in y direction\n",p[0],p[1],p[2]);

    printRow(result.getMode(0),0,  "Result ");
    printRow(compare.getMode(0),0, "Compare");
    
    for(int i=0;i<compare.getNmodes();i++)
      assert( compare.getMode(i).equals( result.getMode(i), 1e-7, true ) );
  }
  //Get twist for first excited momentum in m1 set
  {
    memcpy(p,p_m1,3*sizeof(int));
    p[0] = 3;
    //twist<typename A2Apolicies::FermionFieldType> twist_p(p);
    gaugeFixAndTwist<typename A2Apolicies::FermionFieldType> twist_p(p,lat);   
    compare.fft(W,&twist_p);

    result.getTwistedFFT(p, &Wfftw_p1, &Wfftw_m1);

    if(!UniqueID()) printf("Testing p=(%d,%d,%d). Should require permute of 1 in x direction\n",p[0],p[1],p[2]);

    printRow(result.getMode(0),0,  "Result ");
    printRow(compare.getMode(0),0, "Compare");
    
    for(int i=0;i<compare.getNmodes();i++)
      assert( compare.getMode(i).equals( result.getMode(i), 1e-8, true ) );
  }
  
}



template<typename T>
struct _printit{
  static void printit(const T d[], const int n){
    for(int i=0;i<n;i++){
      std::cout << d[i] << " ";
    }
    std::cout << std::endl;
  }
};

template<typename T>
struct _printit<std::complex<T> >{
  static void printit(const std::complex<T> d[], const int n){
    for(int i=0;i<n;i++){
      std::cout << '[' << d[i].real() << ',' << d[i].imag() << "] ";
    }
    std::cout << std::endl;
  }
};

  
template<typename T>
void printit(const T d[], const int n){
  _printit<T>::printit(d,n);
}



template<typename T>
void printvType(const T& v){
  typedef typename T::scalar_type S;
  int Nsimd = T::Nsimd();
  S to[Nsimd];
  vstore(v,to);
  printit(to,Nsimd);
}


template<typename T>
struct _rand{
  inline static T rand(){
    return LRG.Urand();
  }
};



template<typename T>
struct _rand<std::complex<T> >{
  inline static std::complex<T> rand(){
    return std::complex<T>(LRG.Urand(),LRG.Urand());
  }
};

template<typename T>
T randomvType(){
  T out;
  typedef typename T::scalar_type S;
  int Nsimd = T::Nsimd();
  S r[Nsimd];
  for(int i=0;i<Nsimd;i++) r[i] = _rand<S>::rand();
  vset(out,r);
  return out;
}

CPS_END_NAMESPACE
#include<alg/a2a/mesonfield_computemany.h>
CPS_START_NAMESPACE

template<typename A2Apolicies>
void testMultiSource(const A2AArg &a2a_args,Lattice &lat){
  assert(GJP.Gparity());
  
  if(lat.FixGaugeKind() == FIX_GAUGE_NONE){
    FixGaugeArg fix_gauge_arg;
    fix_gauge_arg.fix_gauge_kind = FIX_GAUGE_COULOMB_T;
    fix_gauge_arg.hyperplane_start = 0;
    fix_gauge_arg.hyperplane_step = 1;
    fix_gauge_arg.hyperplane_num = GJP.Tnodes()*GJP.TnodeSites();
    fix_gauge_arg.stop_cond = 1e-08;
    fix_gauge_arg.max_iter_num = 10000;

    CommonArg common_arg;
  
    AlgFixGauge fix_gauge(lat,&common_arg,&fix_gauge_arg);
    fix_gauge.run();
  }
  
  //Demonstrate relation between FFTW fields
  typedef typename A2Apolicies::ComplexType mf_Complex;
  typedef typename A2AvectorWfftw<A2Apolicies>::FieldInputParamType FieldInputParamType;
  FieldInputParamType fp; defaultFieldParams<FieldInputParamType, mf_Complex>::get(fp);
  
  A2AvectorW<A2Apolicies> W(a2a_args,fp);
  A2AvectorV<A2Apolicies> V(a2a_args,fp);
  W.testRandom();
  V.testRandom();

  int p[3];
  GparityBaseMomentum(p,+1);
  ThreeMomentum pp(p);

  GparityBaseMomentum(p,-1);
  ThreeMomentum pm(p);

  ThreeMomentum pp3 = pp * 3;
  ThreeMomentum pm3 = pm * 3;

  
  typedef typename A2AflavorProjectedExpSource<typename A2Apolicies::SourcePolicies>::FieldParamType SrcFieldParamType;
  typedef typename A2AflavorProjectedExpSource<typename A2Apolicies::SourcePolicies>::ComplexType SrcComplexType;
  SrcFieldParamType sfp; defaultFieldParams<SrcFieldParamType, SrcComplexType>::get(sfp);

  typedef A2AflavorProjectedExpSource<typename A2Apolicies::SourcePolicies> ExpSrcType;
  typedef A2AflavorProjectedHydrogenSource<typename A2Apolicies::SourcePolicies> HydSrcType;
  
  ExpSrcType _1s_src(2.0, pp.ptr(), sfp);
  HydSrcType _2s_src(2,0,0, 2.0, pp.ptr(), sfp);

  typedef SCFspinflavorInnerProduct<15,mf_Complex,ExpSrcType,true,false> ExpInnerType;
  typedef SCFspinflavorInnerProduct<15,mf_Complex,HydSrcType,true,false> HydInnerType;
  
  ExpInnerType _1s_inner(sigma3, _1s_src);
  HydInnerType _2s_inner(sigma3, _2s_src);

  A2AvectorWfftw<A2Apolicies> Wfftw_pp(a2a_args,fp);
  Wfftw_pp.gaugeFixTwistFFT(W,pp.ptr(),lat);

  A2AvectorVfftw<A2Apolicies> Vfftw_pp(a2a_args,fp);
  Vfftw_pp.gaugeFixTwistFFT(V,pp.ptr(),lat);
  
  std::vector< A2AmesonField<A2Apolicies,A2AvectorWfftw,A2AvectorVfftw> > mf_std_1s_pp_pp;
  A2AmesonField<A2Apolicies,A2AvectorWfftw,A2AvectorVfftw>::compute(mf_std_1s_pp_pp, Wfftw_pp, _1s_inner, Vfftw_pp);

  typedef GparityFlavorProjectedBasicSourceStorage<A2Apolicies, ExpInnerType> ExpStorageType;
  
  ExpStorageType exp_store_1s_pp_pp(_1s_inner,_1s_src);
  exp_store_1s_pp_pp.addCompute(0,0,pp,pp);

  std::vector< A2AvectorW<A2Apolicies> const*> Wspecies(1, &W);
  std::vector< A2AvectorV<A2Apolicies> const*> Vspecies(1, &V);

  std::cout << "Start 1s ExpStorage compute\n";
  ComputeMesonFields<A2Apolicies,ExpStorageType>::compute(exp_store_1s_pp_pp,Wspecies,Vspecies,lat);

  int Lt = GJP.Tnodes()*GJP.TnodeSites();
  for(int t=0;t<Lt;t++){
    if(!UniqueID()) printf("Comparing test 1 t=%d\n",t);
    assert( exp_store_1s_pp_pp[0][t].equals(mf_std_1s_pp_pp[t],1e-10,true) );
  }
  if(!UniqueID()) printf("Passed equivalence test 1\n");

  typedef Elem<ExpSrcType,Elem<HydSrcType,ListEnd> > SrcList;
  typedef A2AmultiSource<SrcList> MultiSrcType;
  typedef SCFspinflavorInnerProduct<15,mf_Complex,MultiSrcType,true,false> ExpHydMultiInnerType;

  MultiSrcType exp_hyd_multi_src;
  exp_hyd_multi_src.template getSource<0>().setup(2.0,pp.ptr(),sfp);
  exp_hyd_multi_src.template getSource<1>().setup(2,0,0, 2.0, pp.ptr(), sfp);
  
  ExpHydMultiInnerType exp_hyd_multi_inner(sigma3,exp_hyd_multi_src);

  typedef GparityFlavorProjectedBasicSourceStorage<A2Apolicies, HydInnerType> HydStorageType;
  HydStorageType exp_store_2s_pp_pp(_2s_inner,_2s_src);
  exp_store_2s_pp_pp.addCompute(0,0,pp,pp);
  exp_store_2s_pp_pp.addCompute(0,0,pm,pp);
  exp_store_2s_pp_pp.addCompute(0,0,pp3,pp);

  
  ComputeMesonFields<A2Apolicies,HydStorageType>::compute(exp_store_2s_pp_pp,Wspecies,Vspecies,lat);

  
  typedef GparityFlavorProjectedMultiSourceStorage<A2Apolicies, ExpHydMultiInnerType> ExpHydMultiStorageType;
  ExpHydMultiStorageType exp_store_1s_2s_pp_pp(exp_hyd_multi_inner, exp_hyd_multi_src);
  exp_store_1s_2s_pp_pp.addCompute(0,0,pp,pp);

  std::cout << "Start 1s/2s ExpHydMultiStorage compute\n";
  ComputeMesonFields<A2Apolicies,ExpHydMultiStorageType>::compute(exp_store_1s_2s_pp_pp,Wspecies,Vspecies,lat);
  
  for(int t=0;t<Lt;t++){
    if(!UniqueID()) printf("Comparing test 2 t=%d\n",t);
    assert( exp_store_1s_2s_pp_pp(0,0)[t].equals(mf_std_1s_pp_pp[t],1e-10,true) );
  }
  if(!UniqueID()) printf("Passed equivalence test 2\n");
  for(int t=0;t<Lt;t++){
    if(!UniqueID()) printf("Comparing test 3 t=%d\n",t);
    assert( exp_store_1s_2s_pp_pp(1,0)[t].equals(exp_store_2s_pp_pp[0][t],1e-10,true) );
  }
  if(!UniqueID()) printf("Passed equivalence test 3\n");

  
}


template<typename A2Apolicies>
void testMfFFTreln(const A2AArg &a2a_args,Lattice &lat){
  assert(GJP.Gparity());

  if(lat.FixGaugeKind() == FIX_GAUGE_NONE){
    FixGaugeArg fix_gauge_arg;
    fix_gauge_arg.fix_gauge_kind = FIX_GAUGE_COULOMB_T;
    fix_gauge_arg.hyperplane_start = 0;
    fix_gauge_arg.hyperplane_step = 1;
    fix_gauge_arg.hyperplane_num = GJP.Tnodes()*GJP.TnodeSites();
    fix_gauge_arg.stop_cond = 1e-08;
    fix_gauge_arg.max_iter_num = 10000;

    CommonArg common_arg;
  
    AlgFixGauge fix_gauge(lat,&common_arg,&fix_gauge_arg);
    fix_gauge.run();
  }
  typedef typename A2Apolicies::ComplexType mf_Complex;
  typedef typename A2Apolicies::SourcePolicies SourcePolicies;
  
  typedef typename A2AvectorWfftw<A2Apolicies>::FieldInputParamType FieldInputParamType;
  FieldInputParamType fp; defaultFieldParams<FieldInputParamType, mf_Complex>::get(fp);

  typedef typename A2Apolicies::SourcePolicies::DimensionPolicy::ParamType SrcInputParamType;
  SrcInputParamType sp; defaultFieldParams<SrcInputParamType, mf_Complex>::get(sp);

  A2AvectorW<A2Apolicies> W(a2a_args,fp);
  A2AvectorV<A2Apolicies> V(a2a_args,fp);
  W.testRandom();
  V.testRandom();

  std::vector< A2AvectorW<A2Apolicies> const* > Wspecies(1,&W);
  std::vector< A2AvectorV<A2Apolicies> const* > Vspecies(1,&V);
  
  
  int pp[3]; GparityBaseMomentum(pp,+1); //(1,1,1)
  int pm[3]; GparityBaseMomentum(pm,-1); //(-1,-1,-1)

  //M_ij^{4a+k,4b+l} =  \sum_{n=0}^{L-1} \Omega^{\dagger,4a+k}_i(n) \Gamma \gamma(n) N^{4b+l}_j(n)     (1)
  //                    \sum_{n=0}^{L-1} \Omega^{\dagger,k}_i(n-a-b) \Gamma \gamma(n-b) N^l_j(n)         (2)
  
  //\Omega^{\dagger,k}_i(n) = [ \sum_{x=0}^{L-1} e^{-2\pi i nx/L} e^{- (-k) \pi ix/2L} W_i(x) ]^\dagger
  //N^l_j(n) = \sum_{x=0}^{L-1} e^{-2\pi ix/L} e^{-l \pi ix/2L} V_i(x)

  //Use a state with total momentum 0; k=1 l=-1 a=-1 b=1  so total momentum  -3 + 3  = 0

  int a = -1;
  int b = 1;
  int k = 1;
  int l = -1;

  assert(a+b == 0); //don't want to permute W V right now
  
  //For (1) 
  int p1w[3] = { -(4*a+k), pm[1],pm[1] };  //fix other momenta to first allowed
  int p1v[3] = { 4*b+l, pm[1],pm[1] };
  
  //For (2)
  int p2w[3] = {-k, pm[1],pm[1]};
  int p2v[3] = {l, pm[1],pm[1]};
  
  typedef A2AflavorProjectedExpSource<SourcePolicies> SrcType;
  typedef SCFspinflavorInnerProduct<0,mf_Complex,SrcType,true,false> InnerType; //unit matrix spin structure
  typedef GparityFlavorProjectedBasicSourceStorage<A2Apolicies, InnerType> StorageType;

  SrcType src1(2., pp, sp);
  SrcType src2(2., pp, sp);
  cyclicPermute( src2.getSource(), src2.getSource(), 0, 1, b);

  InnerType inner1(sigma0,src1);
  InnerType inner2(sigma0,src2);
  StorageType mf_store1(inner1,src1);
  StorageType mf_store2(inner2,src2);

  mf_store1.addCompute(0,0, ThreeMomentum(p1w), ThreeMomentum(p1v) );
  mf_store1.addCompute(0,0, ThreeMomentum(p2w), ThreeMomentum(p2v) );
  
  mf_store2.addCompute(0,0, ThreeMomentum(p2w), ThreeMomentum(p2v) );

  ComputeMesonFields<A2Apolicies,StorageType>::compute(mf_store1,Wspecies,Vspecies,lat);
  ComputeMesonFields<A2Apolicies,StorageType>::compute(mf_store2,Wspecies,Vspecies,lat);

  printf("Testing mf relation\n"); fflush(stdout);
  assert( mf_store1[0][0].equals( mf_store2[0][0], 1e-6, true) );
  printf("MF Relation proven\n");

  // StorageType mf_store3(inner1);
  // mf_store3.addCompute(0,0, ThreeMomentum(p1w), ThreeMomentum(p1v), true );
#if 1
  
  typedef GparitySourceShiftInnerProduct<mf_Complex,SrcType,flavorMatrixSpinColorContract<0,mf_Complex,true,false> > ShiftInnerType;
  typedef GparityFlavorProjectedShiftSourceStorage<A2Apolicies, ShiftInnerType> ShiftStorageType;
  
  SrcType src3(2., pp, sp);
  ShiftInnerType shift_inner(sigma0,src3);
  ShiftStorageType mf_shift_store(shift_inner,src3);
  mf_shift_store.addCompute(0,0, ThreeMomentum(p1w), ThreeMomentum(p1v) );
  mf_shift_store.addCompute(0,0, ThreeMomentum(p2w), ThreeMomentum(p2v) );
  int nc = mf_shift_store.nCompute();
  printf("Number of optimized computations: %d\n",nc);

  ComputeMesonFields<A2Apolicies,ShiftStorageType>::compute(mf_shift_store,Wspecies,Vspecies,lat);

  assert( mf_shift_store[0][0].equals( mf_store1[0][0], 1e-6, true) );
  assert( mf_shift_store[1][0].equals( mf_store1[1][0], 1e-6, true) );
  printf("Passed test of shift storage for single source type\n");

  typedef Elem<SrcType, Elem<SrcType,ListEnd > > SrcList;
  typedef A2AmultiSource<SrcList> MultiSrcType;
  typedef GparitySourceShiftInnerProduct<mf_Complex,MultiSrcType,flavorMatrixSpinColorContract<0,mf_Complex,true,false> > ShiftMultiSrcInnerType;
  typedef GparityFlavorProjectedShiftSourceStorage<A2Apolicies, ShiftMultiSrcInnerType> ShiftMultiSrcStorageType;

  MultiSrcType multisrc;
  multisrc.template getSource<0>().setup(3.,pp, sp);
  multisrc.template getSource<1>().setup(2.,pp, sp);
  ShiftMultiSrcInnerType shift_inner_multisrc(sigma0,multisrc);
  ShiftMultiSrcStorageType mf_shift_multisrc_store(shift_inner_multisrc, multisrc);
  mf_shift_multisrc_store.addCompute(0,0, ThreeMomentum(p1w), ThreeMomentum(p1v) );
  mf_shift_multisrc_store.addCompute(0,0, ThreeMomentum(p2w), ThreeMomentum(p2v) );
  
  ComputeMesonFields<A2Apolicies,ShiftMultiSrcStorageType>::compute(mf_shift_multisrc_store,Wspecies,Vspecies,lat);

  assert( mf_shift_multisrc_store(1,0)[0].equals( mf_store1[0][0], 1e-6, true) );
  assert( mf_shift_multisrc_store(1,1)[0].equals( mf_store1[1][0], 1e-6, true) );
  
  
#endif
}
  
//  static void ComputeKtoPiPiGparityBase::multGammaLeft(CPSspinColorFlavorMatrix<ComplexType> &M, const int whichGamma, const int i, const int mu){





template<typename A2Apolicies>
void testFFTopt(){
  typedef typename A2Apolicies::FermionFieldType::FieldSiteType mf_Complex;
  typedef typename A2Apolicies::FermionFieldType::FieldDimensionPolicy DimensionPolicy;
  typedef typename A2Apolicies::FermionFieldType::FieldAllocPolicy AllocPolicy;
  
  typedef CPSfield<mf_Complex,12,DimensionPolicy,OneFlavorPolicy, AllocPolicy> FieldType;
  typedef typename FieldType::InputParamType FieldInputParamType;
  FieldInputParamType fp; defaultFieldParams<FieldInputParamType, mf_Complex>::get(fp);

  bool do_dirs[4] = {1,1,0,0};
  
  FieldType in(fp);
  in.testRandom();

  FieldType out1(fp);
  fft(out1,in,do_dirs);

  FieldType out2(fp);
  fft_opt(out2,in,do_dirs);

  assert( out1.equals(out2, 1e-8, true ) );
  printf("Passed FFT test\n");

  //Test inverse
  FieldType inv(fp);
  fft_opt(inv,out2,do_dirs,true);

  assert( inv.equals(in, 1e-8, true ) );
  printf("Passed FFT inverse test\n");  
}

template<typename A2Apolicies>
void testA2AFFTinv(const A2AArg &a2a_args,Lattice &lat){
  assert(GJP.Gparity());

  if(lat.FixGaugeKind() == FIX_GAUGE_NONE){
    FixGaugeArg fix_gauge_arg;
    fix_gauge_arg.fix_gauge_kind = FIX_GAUGE_COULOMB_T;
    fix_gauge_arg.hyperplane_start = 0;
    fix_gauge_arg.hyperplane_step = 1;
    fix_gauge_arg.hyperplane_num = GJP.Tnodes()*GJP.TnodeSites();
    fix_gauge_arg.stop_cond = 1e-08;
    fix_gauge_arg.max_iter_num = 10000;

    CommonArg common_arg;
  
    AlgFixGauge fix_gauge(lat,&common_arg,&fix_gauge_arg);
    fix_gauge.run();
  }
  typedef typename A2Apolicies::ComplexType mf_Complex;
  typedef typename A2Apolicies::SourcePolicies SourcePolicies;
  
  typedef typename A2AvectorWfftw<A2Apolicies>::FieldInputParamType FieldInputParamType;
  FieldInputParamType fp; defaultFieldParams<FieldInputParamType, mf_Complex>::get(fp);

  typedef typename A2Apolicies::SourcePolicies::DimensionPolicy::ParamType SrcInputParamType;
  SrcInputParamType sp; defaultFieldParams<SrcInputParamType, mf_Complex>::get(sp);

  A2AvectorW<A2Apolicies> W(a2a_args,fp);
  A2AvectorV<A2Apolicies> V(a2a_args,fp);
  W.testRandom();
  V.testRandom();

  int pp[3]; GparityBaseMomentum(pp,+1); //(1,1,1)
  int pm[3]; GparityBaseMomentum(pm,-1); //(-1,-1,-1)
  
  A2AvectorVfftw<A2Apolicies> Vfft(a2a_args,fp);
  Vfft.fft(V);

  A2AvectorV<A2Apolicies> Vrec(a2a_args,fp);
  Vfft.inversefft(Vrec);

  for(int i=0;i<V.getNmodes();i++){
    assert( Vrec.getMode(i).equals( V.getMode(i), 1e-08, true) ); 
  }
  if(!UniqueID()) printf("Passed V fft/inverse test\n");

  A2AvectorWfftw<A2Apolicies> Wfft(a2a_args,fp);
  Wfft.fft(W);

  A2AvectorW<A2Apolicies> Wrec(a2a_args,fp);
  Wfft.inversefft(Wrec);

  for(int i=0;i<W.getNl();i++){
    assert( Wrec.getWl(i).equals( W.getWl(i), 1e-08, true) ); 
  }
  if(!UniqueID()) printf("Passed Wl fft/inverse test\n"); 

  for(int i=0;i<W.getNhits();i++){
    assert( Wrec.getWh(i).equals( W.getWh(i), 1e-08, true) ); 
  }
  if(!UniqueID()) printf("Passed Wh fft/inverse test\n"); 
}

#ifdef USE_GRID
template<typename vComplexType, bool conj_left, bool conj_right>
class GridVectorizedSpinColorContractBasic{
public:
  inline static vComplexType g5(const vComplexType *const l, const vComplexType *const r){
    const static int sc_size =12;
    const static int half_sc = 6;

    vComplexType v3; zeroit(v3);

    for(int i = half_sc; i < sc_size; i++){ 
      v3 -= MconjGrid<vComplexType,conj_left,conj_right>::doit(l+i,r+i);
    }
    for(int i = 0; i < half_sc; i ++){ 
      v3 += MconjGrid<vComplexType,conj_left,conj_right>::doit(l+i,r+i);
    }
    return v3;
  }
};

template<typename T>
typename my_enable_if< is_complex_double_or_float<typename T::scalar_type>::value, bool>::type 
vTypeEquals(const T& a, const T &b, const double tolerance = 1e-12, bool verbose = false){
  typedef typename T::scalar_type S;
  int Nsimd = T::Nsimd();
  S ato[Nsimd];
  vstore(a,ato);
  S bto[Nsimd];
  vstore(b,bto);
  
  bool eq = true;
  for(int i=0;i<Nsimd;i++)
    if( fabs(ato[i].real() - bto[i].real()) > tolerance || fabs(ato[i].imag() - bto[i].imag()) > tolerance ){
      if(verbose && !UniqueID()){	
	double rdiff = fabs(ato[i].real() - bto[i].real());
	double idiff = fabs(ato[i].imag() - bto[i].imag());
	printf("Mismatch index %d: (%g,%g) vs (%g,%g) with diffs (%g,%g)\n",i,ato[i].real(),bto[i].real(),ato[i].imag(),bto[i].imag(),rdiff,idiff);
      }
      eq = false; break;
    }

  if(!eq && verbose && !UniqueID()){
    printf("NOT EQUAL:\n");
    printit(ato,Nsimd);
    printit(bto,Nsimd);
  }    
  return eq;
}

#endif

template<typename mf_Complex>
void testGridg5Contract(){
#ifdef USE_GRID
  Grid::Vector<mf_Complex> vec1(12);
  Grid::Vector<mf_Complex> vec2(12);
  for(int i=0;i<12;i++){
    vec1[i] = randomvType<mf_Complex>();
    vec2[i] = randomvType<mf_Complex>();
  }

  mf_Complex a = GridVectorizedSpinColorContractBasic<mf_Complex,true,false>::g5(vec1.data(),vec2.data());
  mf_Complex b = GridVectorizedSpinColorContract<mf_Complex,true,false>::g5(vec1.data(),vec2.data());
  assert(vTypeEquals(a,b,1e-6,true) == true);
  if(!UniqueID()){ printf("Passed g5 contract repro\n"); fflush(stdout); }
#endif
}

template<typename A2Apolicies>
void testVVdag(Lattice &lat){
  typedef typename A2Apolicies::FermionFieldType::FieldSiteType mf_Complex;
  typedef typename A2Apolicies::FermionFieldType::FieldDimensionPolicy DimensionPolicy;
  typedef typename A2Apolicies::FermionFieldType::FieldAllocPolicy AllocPolicy;
  typedef CPSfermion4D<mf_Complex,DimensionPolicy,OneFlavorPolicy, AllocPolicy> FieldType;
  
  typedef typename FieldType::InputParamType FieldInputParamType;
  FieldInputParamType fp; defaultFieldParams<FieldInputParamType, mf_Complex>::get(fp);
  
  if(lat.FixGaugeKind() == FIX_GAUGE_NONE){
    FixGaugeArg fix_gauge_arg;
    fix_gauge_arg.fix_gauge_kind = FIX_GAUGE_COULOMB_T;
    fix_gauge_arg.hyperplane_start = 0;
    fix_gauge_arg.hyperplane_step = 1;
    fix_gauge_arg.hyperplane_num = GJP.Tnodes()*GJP.TnodeSites();
    fix_gauge_arg.stop_cond = 1e-08;
    fix_gauge_arg.max_iter_num = 10000;

    CommonArg common_arg;
  
    AlgFixGauge fix_gauge(lat,&common_arg,&fix_gauge_arg);
    fix_gauge.run();
  }

  FieldType a(fp);
  a.testRandom();
  
  FieldType Va(a);
  Va.gaugeFix(lat,true,false); //parallel, no dagger

  printRow(a,0,"a");
  printRow(Va,0,"Va");
  
  FieldType VdagVa(Va);
  VdagVa.gaugeFix(lat,true,true);

  printRow(VdagVa,0,"VdagVa");

  assert( VdagVa.equals(a, 1e-8, true) );

  FieldType diff = VdagVa - a;
  printRow(diff,0,"diff");

  double n2 = diff.norm2();
  printf("Norm diff = %g\n",n2);

  FieldType zro(fp); zro.zero();

  assert( diff.equals(zro,1e-12,true));
}


template<typename ManualAllocA2Apolicies>
void testDestructiveFFT(const A2AArg &a2a_args,Lattice &lat){
  assert(GJP.Gparity());

  if(lat.FixGaugeKind() == FIX_GAUGE_NONE){
    FixGaugeArg fix_gauge_arg;
    fix_gauge_arg.fix_gauge_kind = FIX_GAUGE_COULOMB_T;
    fix_gauge_arg.hyperplane_start = 0;
    fix_gauge_arg.hyperplane_step = 1;
    fix_gauge_arg.hyperplane_num = GJP.Tnodes()*GJP.TnodeSites();
    fix_gauge_arg.stop_cond = 1e-08;
    fix_gauge_arg.max_iter_num = 10000;

    CommonArg common_arg;
  
    AlgFixGauge fix_gauge(lat,&common_arg,&fix_gauge_arg);
    fix_gauge.run();
  }
  typedef typename ManualAllocA2Apolicies::FermionFieldType FermionFieldType;
  typedef typename ManualAllocA2Apolicies::SourcePolicies SourcePolicies;
  typedef typename ManualAllocA2Apolicies::ComplexType mf_Complex;
  
  typedef typename A2AvectorWfftw<ManualAllocA2Apolicies>::FieldInputParamType FieldInputParamType;
  FieldInputParamType fp; defaultFieldParams<FieldInputParamType, mf_Complex>::get(fp);

  typedef typename ManualAllocA2Apolicies::SourcePolicies::DimensionPolicy::ParamType SrcInputParamType;
  SrcInputParamType sp; defaultFieldParams<SrcInputParamType, mf_Complex>::get(sp);

  A2AvectorW<ManualAllocA2Apolicies> W(a2a_args,fp);
  A2AvectorV<ManualAllocA2Apolicies> V(a2a_args,fp);
  
  for(int i=0;i<V.getNmodes();i++) assert( &V.getMode(i) == NULL);
  V.allocModes();
  for(int i=0;i<V.getNmodes();i++) assert( &V.getMode(i) != NULL);
  
  V.testRandom();

  W.allocModes();
  for(int i=0;i<W.getNl();i++) assert( &W.getWl(i) != NULL);
  for(int i=0;i<W.getNhits();i++) assert( &W.getWh(i) != NULL);
  W.testRandom();

  
  A2AvectorV<ManualAllocA2Apolicies> Vcopy = V;
  A2AvectorW<ManualAllocA2Apolicies> Wcopy = W;
  
  int pp[3]; GparityBaseMomentum(pp,+1); //(1,1,1)
  int pm[3]; GparityBaseMomentum(pm,-1); //(-1,-1,-1)

  gaugeFixAndTwist<FermionFieldType> fft_op(pp,lat);  
  reverseGaugeFixAndTwist<FermionFieldType> invfft_op(pp,lat);
  
  A2AvectorVfftw<ManualAllocA2Apolicies> Vfft(a2a_args,fp); //no allocation yet performed
  Vfft.destructivefft(V, &fft_op);

  for(int i=0;i<V.getNmodes();i++) assert( &V.getMode(i) == NULL);
  for(int i=0;i<Vfft.getNmodes();i++) assert( &Vfft.getMode(i) != NULL);

  
  A2AvectorV<ManualAllocA2Apolicies> Vrec(a2a_args,fp);
  Vfft.destructiveInversefft(Vrec, &invfft_op);

  for(int i=0;i<Vrec.getNmodes();i++) assert( &Vrec.getMode(i) != NULL);
  for(int i=0;i<Vfft.getNmodes();i++) assert( &Vfft.getMode(i) == NULL); 

  for(int i=0;i<Vrec.getNmodes();i++) assert( Vrec.getMode(i).equals( Vcopy.getMode(i), 1e-08, true) );

  
  printf("Passed V destructive fft/inverse test\n");
   
  A2AvectorWfftw<ManualAllocA2Apolicies> Wfft(a2a_args,fp);
  Wfft.destructiveGaugeFixTwistFFT(W,pp,lat);

  for(int i=0;i<W.getNl();i++) assert( &W.getWl(i) == NULL);
  for(int i=0;i<W.getNhits();i++) assert( &W.getWh(i) == NULL);
  
  for(int i=0;i<Wfft.getNmodes();i++) assert( &Wfft.getMode(i) != NULL);
  
  A2AvectorW<ManualAllocA2Apolicies> Wrec(a2a_args,fp);
  Wfft.destructiveUnapplyGaugeFixTwistFFT(Wrec, pp,lat);
  
  for(int i=0;i<Wfft.getNmodes();i++) assert( &Wfft.getMode(i) == NULL);

  for(int i=0;i<Wrec.getNl();i++) assert( &Wrec.getWl(i) != NULL);
  for(int i=0;i<Wrec.getNhits();i++) assert( &Wrec.getWh(i) != NULL);
  
  for(int i=0;i<Wrec.getNl();i++){
    assert( Wrec.getWl(i).equals( Wcopy.getWl(i), 1e-08, true) ); 
  }
  if(!UniqueID()) printf("Passed Wl destructive fft/inverse test\n"); 

  for(int i=0;i<Wrec.getNhits();i++){
    assert( Wrec.getWh(i).equals( Wcopy.getWh(i), 1e-08, true) ); 
  }
  if(!UniqueID()) printf("Passed Wh destructive fft/inverse test\n"); 
  
  
}


void testA2AallocFree(const A2AArg &a2a_args,Lattice &lat){
  typedef A2ApoliciesSIMDdoubleManualAlloc A2Apolicies;
  //typedef A2ApoliciesDoubleManualAlloc A2Apolicies;

  typedef typename A2Apolicies::FermionFieldType FermionFieldType;
  typedef typename A2Apolicies::SourcePolicies SourcePolicies;
  typedef typename A2Apolicies::ComplexType mf_Complex;
  
  typedef typename A2AvectorWfftw<A2Apolicies>::FieldInputParamType FieldInputParamType;
  FieldInputParamType fp; defaultFieldParams<FieldInputParamType, mf_Complex>::get(fp);

  A2AvectorVfftw<A2Apolicies> Vfft(a2a_args,fp);
  double size =  A2AvectorVfftw<A2Apolicies>::Mbyte_size(a2a_args,fp);
  
  for(int i=0;i<100;i++){
    if(!UniqueID()) printf("Pre-init\n");
    printMem(); fflush(stdout);

    if(!UniqueID()) printf("Expected size %f MB\n",size);
    
    if(!UniqueID()) printf("Post-init\n");
    printMem(); fflush(stdout);

    Vfft.allocModes();

    for(int i=0;i<Vfft.getNmodes();i++){
      assert(&Vfft.getMode(i) != NULL);
      Vfft.getMode(i).zero();
    }
    if(!UniqueID()) printf("Post-alloc\n");
    printMem(); fflush(stdout);

    Vfft.freeModes();

    for(int i=0;i<Vfft.getNmodes();i++)
      assert(&Vfft.getMode(i) == NULL);
    
    if(!UniqueID()) printf("Post-free\n");
    printMem(); fflush(stdout);
  }

  exit(0);
}


template<typename GridA2Apolicies>
void benchmarkMFcontractKernel(const int ntests, const int nthreads){
#ifdef USE_GRID
  // GridVectorizedSpinColorContract benchmark
  typedef typename GridA2Apolicies::ComplexType GVtype;
  typedef typename GridA2Apolicies::ScalarComplexType GCtype;
  const int nsimd = GVtype::Nsimd();      

  FourDSIMDPolicy::ParamType simd_dims;
  FourDSIMDPolicy::SIMDdefaultLayout(simd_dims,nsimd,2);
  
  NullObject n;
  CPSfield<GCtype,12,FourDpolicy,OneFlavorPolicy> a(n); a.testRandom();
  CPSfield<GCtype,12,FourDpolicy,OneFlavorPolicy> b(n); b.testRandom();
  CPSfield<GVtype,12,FourDSIMDPolicy,OneFlavorPolicy,Aligned128AllocPolicy> aa(simd_dims); aa.importField(a);
  CPSfield<GVtype,12,FourDSIMDPolicy,OneFlavorPolicy,Aligned128AllocPolicy> bb(simd_dims); bb.importField(b);
  CPSfield<GVtype,1,FourDSIMDPolicy,OneFlavorPolicy,Aligned128AllocPolicy> cc(simd_dims);

  int ntests_scaled = ntests * 1000;
  printf("Max threads %d\n",omp_get_max_threads());
#ifdef TIMERS_OFF
  printf("Timers are OFF\n"); fflush(stdout);
#else
  printf("Timers are ON\n"); fflush(stdout);
#endif
  __itt_resume();

  for(int oloop=0; oloop < 100; oloop++){
    double t0 = Grid::usecond();

#pragma omp parallel //avoid thread creation overheads
    {
      int me = omp_get_thread_num();
      int work, off;
      thread_work(work, off, aa.nfsites(), me, omp_get_num_threads());
	
      GVtype *abase = aa.fsite_ptr(off);
      GVtype *bbase = bb.fsite_ptr(off);
      GVtype *cbase = cc.fsite_ptr(off);

      for(int test=0;test<ntests_scaled;test++){
	GVtype *ai = abase;
	GVtype *bi = bbase;
	GVtype *ci = cbase;
	__SSC_MARK(0x1);
	for(int i=0;i<work;i++){
	  *ci = GridVectorizedSpinColorContract<GVtype,true,false>::g5(ai,bi);
	  ai += 12;
	  bi += 12;
	  ci += 1;
	}
	__SSC_MARK(0x2);
      }
    }

    double t1 = Grid::usecond();
    double dt = t1 - t0;
      
    int FLOPs = 12*6*nsimd //12 vectorized conj(a)*b
      + 12*2*nsimd; //12 vectorized += or -=
    double total_FLOPs = double(FLOPs) * double(aa.nfsites()) * double(ntests_scaled);
      
    double flops = total_FLOPs/dt; //dt in us   dt/(1e-6 s) in Mflops
    std::cout << "GridVectorizedSpinColorContract( conj(a)*b ): New code " << ntests_scaled << " tests over " << nthreads << " threads: Time " << dt << " usecs  flops " << flops/1e3 << " Gflops\n";

  }
  __itt_detach();
#endif
}


template<typename ScalarA2Apolicies, typename GridA2Apolicies>
void testVVgridOrig(const A2AArg &a2a_args, const int ntests, const int nthreads, const double tol){
  std::cout << "Starting vv test/timing\n";

  const int nsimd = GridA2Apolicies::ComplexType::Nsimd();      

  FourDSIMDPolicy::ParamType simd_dims;
  FourDSIMDPolicy::SIMDdefaultLayout(simd_dims,nsimd,2);
      
  A2AvectorWfftw<ScalarA2Apolicies> W(a2a_args);
  A2AvectorVfftw<ScalarA2Apolicies> V(a2a_args);
    
  A2AvectorWfftw<GridA2Apolicies> Wgrid(a2a_args, simd_dims);
  A2AvectorVfftw<GridA2Apolicies> Vgrid(a2a_args, simd_dims);

  W.testRandom();
  V.testRandom();
  Wgrid.importFields(W);
  Vgrid.importFields(V);
  
  typedef typename GridA2Apolicies::ComplexType grid_Complex;
  typedef typename ScalarA2Apolicies::ComplexType mf_Complex;
      
  Float total_time = 0.;
  Float total_time_orig = 0.;
  CPSspinColorFlavorMatrix<mf_Complex> orig_sum[nthreads];
  CPSspinColorFlavorMatrix<grid_Complex> grid_sum[nthreads];

  CPSspinColorFlavorMatrix<mf_Complex> orig_tmp[nthreads];
  CPSspinColorFlavorMatrix<grid_Complex> grid_tmp[nthreads];

  int orig_3vol = GJP.VolNodeSites()/GJP.TnodeSites();
  int grid_3vol = Vgrid.getMode(0).nodeSites(0) * Vgrid.getMode(0).nodeSites(1) *Vgrid.getMode(0).nodeSites(2);
      
  for(int iter=0;iter<ntests;iter++){
    for(int i=0;i<nthreads;i++){
      orig_sum[i].zero(); grid_sum[i].zero();
    }
	
    for(int top = 0; top < GJP.TnodeSites(); top++){
      //std::cout << "top " << top << std::endl;
      //std::cout << "Starting orig\n";
      total_time_orig -= dclock();	  
#pragma omp parallel for
      for(int xop=0;xop<orig_3vol;xop++){
	int me = omp_get_thread_num();
	mult(orig_tmp[me], V, W, xop, top, false, true);
	orig_sum[me] += orig_tmp[me];
      }
      total_time_orig += dclock();
      //std::cout << "Starting Grid\n";
      total_time -= dclock();
#pragma omp parallel for
      for(int xop=0;xop<grid_3vol;xop++){
	int me = omp_get_thread_num();
	mult(grid_tmp[me], Vgrid, Wgrid, xop, top, false, true);
	grid_sum[me] += grid_tmp[me];
      }
      total_time += dclock();	  
    }
    for(int i=1;i<nthreads;i++){
      orig_sum[0] += orig_sum[i];
      grid_sum[0] += grid_sum[i];
    }

	
    bool fail = false;
	
    typename GridA2Apolicies::ScalarComplexType gd;
    for(int sl=0;sl<4;sl++)
      for(int cl=0;cl<3;cl++)
	for(int fl=0;fl<2;fl++)
	  for(int sr=0;sr<4;sr++)
	    for(int cr=0;cr<3;cr++)
	      for(int fr=0;fr<2;fr++){
		gd = Reduce( grid_sum[0](sl,sr)(cl,cr)(fl,fr) );
		const mf_Complex &cp = orig_sum[0](sl,sr)(cl,cr)(fl,fr);

		double rdiff = fabs(gd.real()-cp.real());
		double idiff = fabs(gd.imag()-cp.imag());
		if(rdiff > tol|| idiff > tol){
		  printf("Fail: Iter %d Grid (%g,%g) CPS (%g,%g) Diff (%g,%g)\n",iter, gd.real(),gd.imag(), cp.real(),cp.imag(), cp.real()-gd.real(), cp.imag()-gd.imag());
		  fail = true;
		}
	      }

    if(fail) ERR.General("","","Standard vs Grid implementation test failed\n");
  }

  printf("vv: Avg time new code %d iters: %g secs\n",ntests,total_time/ntests);
  printf("vv: Avg time old code %d iters: %g secs\n",ntests,total_time_orig/ntests);
}

template<typename ScalarA2Apolicies, typename GridA2Apolicies>
void testvMvGridOrig(const A2AArg &a2a_args, const int ntests, const int nthreads, const double tol){
  //#define CPS_VMV
  //#define GRID_VMV
  //#define CPS_SPLIT_VMV
#define GRID_SPLIT_VMV
  //#define CPS_SPLIT_VMV_XALL
  //#define GRID_SPLIT_VMV_XALL

  std::cout << "Starting vMv benchmark\n";

  const int nsimd = GridA2Apolicies::ComplexType::Nsimd();      

  FourDSIMDPolicy::ParamType simd_dims;
  FourDSIMDPolicy::SIMDdefaultLayout(simd_dims,nsimd,2);
      
  A2AvectorWfftw<ScalarA2Apolicies> W(a2a_args);
  A2AvectorVfftw<ScalarA2Apolicies> V(a2a_args);
    
  A2AvectorWfftw<GridA2Apolicies> Wgrid(a2a_args, simd_dims);
  A2AvectorVfftw<GridA2Apolicies> Vgrid(a2a_args, simd_dims);

  W.testRandom();
  V.testRandom();
  Wgrid.importFields(W);
  Vgrid.importFields(V);
  
  A2AmesonField<ScalarA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> mf;
  A2AmesonField<GridA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> mf_grid;
  mf.setup(W,V,0,0);
  mf_grid.setup(Wgrid,Vgrid,0,0);     
  mf.testRandom();
  for(int i=0;i<mf.getNrows();i++)
    for(int j=0;j<mf.getNcols();j++)
      mf_grid(i,j) = mf(i,j); //both are scalar complex
  
  typedef typename GridA2Apolicies::ComplexType grid_Complex;
  typedef typename ScalarA2Apolicies::ComplexType mf_Complex;

      
  Float total_time = 0.;
  Float total_time_orig = 0.;
  Float total_time_split_orig = 0.;
  Float total_time_split_grid = 0.;
  Float total_time_split_orig_xall = 0.;
  Float total_time_split_grid_xall = 0.;
       
  CPSspinColorFlavorMatrix<mf_Complex> orig_sum[nthreads];
  CPSspinColorFlavorMatrix<grid_Complex> grid_sum[nthreads];

  CPSspinColorFlavorMatrix<mf_Complex> orig_tmp[nthreads];
  CPSspinColorFlavorMatrix<grid_Complex> grid_tmp[nthreads];

  CPSspinColorFlavorMatrix<mf_Complex> orig_sum_split[nthreads];
  CPSspinColorFlavorMatrix<grid_Complex> grid_sum_split[nthreads];

  CPSspinColorFlavorMatrix<mf_Complex> orig_sum_split_xall[nthreads];
  CPSspinColorFlavorMatrix<grid_Complex> grid_sum_split_xall[nthreads];

      
  int orig_3vol = GJP.VolNodeSites()/GJP.TnodeSites();
  int grid_3vol = Vgrid.getMode(0).nodeSites(0) * Vgrid.getMode(0).nodeSites(1) *Vgrid.getMode(0).nodeSites(2);

  mult_vMv_split<ScalarA2Apolicies, A2AvectorVfftw, A2AvectorWfftw, A2AvectorVfftw, A2AvectorWfftw> vmv_split_orig;
  mult_vMv_split<GridA2Apolicies, A2AvectorVfftw, A2AvectorWfftw, A2AvectorVfftw, A2AvectorWfftw> vmv_split_grid;

  std::vector<CPSspinColorFlavorMatrix<mf_Complex>> orig_split_xall_tmp(orig_3vol);
  Grid::Vector<CPSspinColorFlavorMatrix<grid_Complex> > grid_split_xall_tmp(grid_3vol);
      
  for(int iter=0;iter<ntests;iter++){
    for(int i=0;i<nthreads;i++){
      orig_sum[i].zero(); grid_sum[i].zero();
      orig_sum_split[i].zero(); grid_sum_split[i].zero();
      orig_sum_split_xall[i].zero(); grid_sum_split_xall[i].zero();
    }
	
    for(int top = 0; top < GJP.TnodeSites(); top++){
#ifdef CPS_VMV
      //ORIG VMV
      total_time_orig -= dclock();	  
#pragma omp parallel for
      for(int xop=0;xop<orig_3vol;xop++){
	int me = omp_get_thread_num();
	mult(orig_tmp[me], V, mf, W, xop, top, false, true);
	orig_sum[me] += orig_tmp[me];
      }
      total_time_orig += dclock();
#endif
#ifdef GRID_VMV
      //GRID VMV
      total_time -= dclock();
#pragma omp parallel for
      for(int xop=0;xop<grid_3vol;xop++){
	int me = omp_get_thread_num();
	mult(grid_tmp[me], Vgrid, mf_grid, Wgrid, xop, top, false, true);
	grid_sum[me] += grid_tmp[me];
      }
      total_time += dclock();
#endif

#ifdef CPS_SPLIT_VMV
      //SPLIT VMV
      total_time_split_orig -= dclock();	  
      vmv_split_orig.setup(V, mf, W, top);

#pragma omp parallel for
      for(int xop=0;xop<orig_3vol;xop++){
	int me = omp_get_thread_num();
	vmv_split_orig.contract(orig_tmp[me], xop, false, true);
	orig_sum_split[me] += orig_tmp[me];
      }
      total_time_split_orig += dclock();
#endif

#ifdef GRID_SPLIT_VMV
      //SPLIT VMV GRID
      total_time_split_grid -= dclock();	  
      vmv_split_grid.setup(Vgrid, mf_grid, Wgrid, top);

#pragma omp parallel for
      for(int xop=0;xop<grid_3vol;xop++){
	int me = omp_get_thread_num();
	vmv_split_grid.contract(grid_tmp[me], xop, false, true);
	grid_sum_split[me] += grid_tmp[me];
      }
      total_time_split_grid += dclock();
#endif

#ifdef CPS_SPLIT_VMV_XALL	  	 
      //SPLIT VMV THAT DOES IT FOR ALL SITES
      total_time_split_orig_xall -= dclock();	  
      vmv_split_orig.setup(V, mf, W, top);
      vmv_split_orig.contract(orig_split_xall_tmp, false, true);
#pragma omp parallel for
      for(int xop=0;xop<orig_3vol;xop++){
	int me = omp_get_thread_num();
	orig_sum_split_xall[me] += orig_split_xall_tmp[xop];
      }
      total_time_split_orig_xall += dclock();
#endif

#ifdef GRID_SPLIT_VMV_XALL
      //SPLIT VMV GRID THAT DOES IT FOR ALL SITES
      total_time_split_grid_xall -= dclock();	  
      vmv_split_grid.setup(Vgrid, mf_grid, Wgrid, top);
      vmv_split_grid.contract(grid_split_xall_tmp, false, true);
#pragma omp parallel for
      for(int xop=0;xop<grid_3vol;xop++){
	int me = omp_get_thread_num();	    
	grid_sum_split_xall[me] += grid_split_xall_tmp[xop];
      }
      total_time_split_grid_xall += dclock();
#endif	  
    }//end top loop
    for(int i=1;i<nthreads;i++){
      orig_sum[0] += orig_sum[i];
      grid_sum[0] += grid_sum[i];
      orig_sum_split[0] += orig_sum_split[i];
      grid_sum_split[0] += grid_sum_split[i];
      orig_sum_split_xall[0] += orig_sum_split_xall[i];
      grid_sum_split_xall[0] += grid_sum_split_xall[i];  
    }
#ifdef CPS_VMV
    if(iter == 0){
#  ifdef GRID_VMV
      if(!compare(orig_sum[0],grid_sum[0],tol)) ERR.General("","","Standard vs Grid implementation test failed\n");
#  endif
#  ifdef CPS_SPLIT_VMV
      if(!compare(orig_sum[0],orig_sum_split[0],tol)) ERR.General("","","Standard vs Split implementation test failed\n");
#  endif
#  ifdef GRID_SPLIT_VMV
      if(!compare(orig_sum[0],grid_sum_split[0],tol)) ERR.General("","","Standard vs Grid Split implementation test failed\n");
#  endif
#  ifdef CPS_SPLIT_VMV_XALL
      if(!compare(orig_sum[0],orig_sum_split_xall[0],tol)) ERR.General("","","Standard vs Split xall implementation test failed\n");
#  endif
#  ifdef GRID_SPLIT_VMV_XALL
      if(!compare(orig_sum[0],grid_sum_split_xall[0],tol)) ERR.General("","","Standard vs Grid split xall implementation test failed\n");
#  endif
    }
#endif
  }
#ifdef CPS_VMV
  printf("vMv: Avg time old code %d iters: %g secs\n",ntests,total_time_orig/ntests);
#endif
#ifdef GRID_VMV
  printf("vMv: Avg time new code %d iters: %g secs\n",ntests,total_time/ntests);
#endif
#ifdef CPS_SPLIT_VMV
  printf("vMv: Avg time old code split %d iters: %g secs\n",ntests,total_time_split_orig/ntests);
#endif
#ifdef GRID_SPLIT_VMV
  printf("vMv: Avg time new code split %d iters: %g secs\n",ntests,total_time_split_grid/ntests);
#endif
#ifdef CPS_SPLIT_VMV_XALL
  printf("vMv: Avg time old code split xall %d iters: %g secs\n",ntests,total_time_split_orig_xall/ntests);
#endif
#ifdef CPS_SPLIT_VMV_XALL
  printf("vMv: Avg time new code split xall %d iters: %g secs\n",ntests,total_time_split_grid_xall/ntests);
#endif
}

template<typename ScalarA2Apolicies>
void testMesonFieldReadWrite(const A2AArg &a2a_args){
  A2AvectorWfftw<ScalarA2Apolicies> W(a2a_args);
  A2AvectorVfftw<ScalarA2Apolicies> V(a2a_args);
  
  A2AmesonField<ScalarA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> mf;
  mf.testRandom();
  
  {
    mf.write("mesonfield.dat",FP_IEEE64BIG);
    A2AmesonField<ScalarA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> mfr;
    mfr.read("mesonfield.dat");
    assert( mfr.equals(mf,1e-18,true));
    if(!UniqueID()) printf("Passed mf single IO test\n");
  }
  {
    A2AmesonField<ScalarA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> mfa;
    A2AmesonField<ScalarA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> mfb;
    A2AmesonField<ScalarA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> mfc;
    mfa.setup(W,V,0,0);
    mfb.setup(W,V,1,1);
    mfc.setup(W,V,2,2);		
      
    mfa.testRandom();
    mfb.testRandom();
    mfc.testRandom();

    std::ofstream *fp = !UniqueID() ? new std::ofstream("mesonfield_many.dat") : NULL;

    mfa.write(fp,FP_IEEE64BIG);
    mfb.write(fp,FP_IEEE64LITTLE);
    mfc.write(fp,FP_IEEE64BIG);

    if(!UniqueID()) fp->close();

    A2AmesonField<ScalarA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> mfra;
    A2AmesonField<ScalarA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> mfrb;
    A2AmesonField<ScalarA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> mfrc;

    std::ifstream *ifp = !UniqueID() ? new std::ifstream("mesonfield_many.dat") : NULL;

    mfra.read(ifp);
    mfrb.read(ifp);
    mfrc.read(ifp);

    if(!UniqueID()) ifp->close();

    assert( mfra.equals(mfa,1e-18,true) );
    assert( mfrb.equals(mfb,1e-18,true) );
    assert( mfrc.equals(mfc,1e-18,true) );
    if(!UniqueID()) printf("Passed mf multi IO test\n");
  }
  {
    std::vector< A2AmesonField<ScalarA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> > mfv(3);
    for(int i=0;i<3;i++){
      mfv[i].setup(W,V,i,i);
      mfv[i].testRandom();
    }
    A2AmesonField<ScalarA2Apolicies,A2AvectorWfftw,A2AvectorVfftw>::write("mesonfield_vec.dat", mfv, FP_IEEE64LITTLE);
	
    std::vector< A2AmesonField<ScalarA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> > mfrv;
    A2AmesonField<ScalarA2Apolicies,A2AvectorWfftw,A2AvectorVfftw>::read("mesonfield_vec.dat", mfrv);

    for(int i=0;i<3;i++)
      assert( mfrv[i].equals(mfv[i], 1e-18, true) );
    if(!UniqueID()) printf("Passed mf vector IO test\n");
  }
}	


template<typename ScalarA2Apolicies, typename GridA2Apolicies>
void testMFcontract(const A2AArg &a2a_args, const int nthreads, const double tol){
  std::cout << "Starting MF contraction test\n";

  const int nsimd = GridA2Apolicies::ComplexType::Nsimd();      

  FourDSIMDPolicy::ParamType simd_dims;
  FourDSIMDPolicy::SIMDdefaultLayout(simd_dims,nsimd,2);
      
  A2AvectorWfftw<ScalarA2Apolicies> W(a2a_args);
  A2AvectorVfftw<ScalarA2Apolicies> V(a2a_args);
    
  A2AvectorWfftw<GridA2Apolicies> Wgrid(a2a_args, simd_dims);
  A2AvectorVfftw<GridA2Apolicies> Vgrid(a2a_args, simd_dims);
  
  std::vector<A2AmesonField<ScalarA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> > mf;
  std::vector<A2AmesonField<GridA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> > mf_grid;
  
  typedef typename GridA2Apolicies::ComplexType grid_Complex;
  typedef typename ScalarA2Apolicies::ComplexType mf_Complex;

  typedef typename GridA2Apolicies::ScalarComplexType Ctype;
  typedef typename Ctype::value_type Ftype;
  
  ThreeDSIMDPolicy::ParamType simd_dims_3d;
  ThreeDSIMDPolicy::SIMDdefaultLayout(simd_dims_3d,nsimd);

  typedef typename GridA2Apolicies::SourcePolicies GridSrcPolicy;    
  int p[3] = {1,1,1};
  A2AflavorProjectedExpSource<GridSrcPolicy> src_grid(2.0,p,simd_dims_3d);
  typedef SCFspinflavorInnerProduct<15,typename GridA2Apolicies::ComplexType,A2AflavorProjectedExpSource<GridSrcPolicy> > GridInnerProduct;
  GridInnerProduct mf_struct_grid(sigma3,src_grid);
  
  A2AflavorProjectedExpSource<> src(2.0,p);
  typedef SCFspinflavorInnerProduct<15,typename ScalarA2Apolicies::ComplexType,A2AflavorProjectedExpSource<> > StdInnerProduct;
  StdInnerProduct mf_struct(sigma3,src);

  W.testRandom();
  V.testRandom();
  Wgrid.importFields(W);
  Vgrid.importFields(V);
      
  //A2AmesonField<GridA2Apolicies,A2AvectorWfftw,A2AvectorVfftw>::compute(mf_grid,Wgrid,mf_struct_grid,Vgrid,0);
  {
    typedef typename std::vector<A2AmesonField<GridA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> >::allocator_type Allocator;
    typedef SingleSrcVectorPoliciesSIMD<GridA2Apolicies,A2AvectorWfftw,A2AvectorVfftw,Allocator,GridInnerProduct> VectorPolicies;
    mfComputeGeneral<GridA2Apolicies,A2AvectorWfftw,A2AvectorVfftw, GridInnerProduct, VectorPolicies> cg;
    cg.compute(mf_grid,Wgrid,mf_struct_grid,Vgrid, true);
  }

  
  A2AmesonField<ScalarA2Apolicies,A2AvectorWfftw,A2AvectorVfftw>::compute(mf,W,mf_struct,V);

  bool fail = false;
  for(int t=0;t<mf.size();t++){
    for(int i=0;i<mf[t].size();i++){
      const Ctype& gd = mf_grid[t].ptr()[i];
      const Ctype& cp = mf[t].ptr()[i];
      Ftype rdiff = fabs(gd.real()-cp.real());
      Ftype idiff = fabs(gd.imag()-cp.imag());
      if(rdiff > tol|| idiff > tol){
	printf("Fail: t %d Grid (%g,%g) CPS (%g,%g) Diff (%g,%g)\n",t, gd.real(),gd.imag(), cp.real(),cp.imag(), cp.real()-gd.real(), cp.imag()-gd.imag());
	fail = true;
      }
    }
  }
  if(fail) ERR.General("","","Standard vs Grid implementation test failed\n");
  else if(!UniqueID()){ printf("Passed MF contraction test\n"); fflush(stdout); }
}



template<typename ScalarA2Apolicies, typename GridA2Apolicies>
void benchmarkMFcontract(const A2AArg &a2a_args, const int ntests, const int nthreads){
  typedef typename GridA2Apolicies::SourcePolicies GridSrcPolicy;
  typedef typename ScalarA2Apolicies::ScalarComplexType Ctype;
  typedef typename Ctype::value_type Ftype;

  const int nsimd = GridA2Apolicies::ComplexType::Nsimd();      

  FourDSIMDPolicy::ParamType simd_dims;
  FourDSIMDPolicy::SIMDdefaultLayout(simd_dims,nsimd,2);

  A2AvectorWfftw<ScalarA2Apolicies> W(a2a_args);
  A2AvectorVfftw<ScalarA2Apolicies> V(a2a_args);
  
  A2AvectorWfftw<GridA2Apolicies> Wgrid(a2a_args, simd_dims);
  A2AvectorVfftw<GridA2Apolicies> Vgrid(a2a_args, simd_dims);
  A2AmesonField<GridA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> mf_grid;

  ThreeDSIMDPolicy::ParamType simd_dims_3d;
  ThreeDSIMDPolicy::SIMDdefaultLayout(simd_dims_3d,nsimd);

  printf("Nsimd = %d, SIMD dimensions:\n", nsimd);
  for(int i=0;i<4;i++)
    printf("%d ", simd_dims[i]);
  printf("\n");
  
  int p[3] = {1,1,1};
  A2AflavorProjectedExpSource<GridSrcPolicy> src_grid(2.0,p,simd_dims_3d);
  typedef SCFspinflavorInnerProduct<15,typename GridA2Apolicies::ComplexType,A2AflavorProjectedExpSource<GridSrcPolicy> > GridInnerProduct;
  GridInnerProduct mf_struct_grid(sigma3,src_grid);
      
  std::cout << "Starting all-time mesonfield contract benchmark\n";
  if(!UniqueID()) printf("Using outer blocking bi %d bj %d bp %d\n",BlockedMesonFieldArgs::bi,BlockedMesonFieldArgs::bj,BlockedMesonFieldArgs::bp);
  if(!UniqueID()) printf("Using inner blocking bi %d bj %d bp %d\n",BlockedMesonFieldArgs::bii,BlockedMesonFieldArgs::bjj,BlockedMesonFieldArgs::bpp);

  Float total_time = 0.;
  std::vector<A2AmesonField<GridA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> > mf_grid_t;

  W.testRandom();
  V.testRandom();
  Wgrid.importFields(W);
  Vgrid.importFields(V);
      
  CALLGRIND_START_INSTRUMENTATION ;
  CALLGRIND_TOGGLE_COLLECT ;
  
  typedef typename std::vector<A2AmesonField<GridA2Apolicies,A2AvectorWfftw,A2AvectorVfftw> >::allocator_type Allocator;
  typedef SingleSrcVectorPoliciesSIMD<GridA2Apolicies,A2AvectorWfftw,A2AvectorVfftw,Allocator,GridInnerProduct> VectorPolicies;
  mfComputeGeneral<GridA2Apolicies,A2AvectorWfftw,A2AvectorVfftw, GridInnerProduct, VectorPolicies> cg;
  
  for(int iter=0;iter<ntests;iter++){
    total_time -= dclock();

    cg.compute(mf_grid_t,Wgrid,mf_struct_grid,Vgrid, true);

    
    //A2AmesonField<GridA2Apolicies,A2AvectorWfftw,A2AvectorVfftw>::compute(mf_grid_t,Wgrid,mf_struct_grid,Vgrid);
    total_time += dclock();
  }

  CALLGRIND_TOGGLE_COLLECT ;
  CALLGRIND_STOP_INSTRUMENTATION ;

  int g5_FLOPs = 12*6*nsimd + 12*2*nsimd;//4 flav * 12 vectorized conj(a)*b  + 12 vectorized += or -=         
  int siteFmat_FLOPs = 3*nsimd;  //1 vectorized z.im*-1, 1 vectorized -1*z                                                                                                                             
  int s3_FLOPs = 4*nsimd; //2 vectorized -1*z                                                                                                                                                          
  int TransLeftTrace_FLOPs = nsimd*4*6 + nsimd*3*2; //4 vcmul + 3vcadd                                                                                                                                 
  int reduce_FLOPs = (nsimd - 1)*2; //nsimd-1 cadd                                                                                                                                                    

  double FLOPs_per_site = 0.;
  for(int t=GJP.TnodeCoor()*GJP.TnodeSites(); t<(GJP.TnodeCoor()+1)*GJP.TnodeSites(); t++){
    const int nl_l = mf_grid_t[t].getRowParams().getNl();
    const int nl_r = mf_grid_t[t].getColParams().getNl();

    int t_lcl = t-GJP.TnodeCoor()*GJP.TnodeSites();

    for(int i = 0; i < mf_grid_t[t].getNrows(); i++){
      modeIndexSet i_high_unmapped; if(i>=nl_l) mf_grid_t[t].getRowParams().indexUnmap(i-nl_l,i_high_unmapped);
      SCFvectorPtr<typename GridA2Apolicies::FermionFieldType::FieldSiteType> lscf = Wgrid.getFlavorDilutedVect(i,i_high_unmapped,0,t_lcl); //dilute flavor in-place if it hasn't been already \
                                                                                                                                                                                                           
      for(int j = 0; j < mf_grid_t[t].getNcols(); j++) {
	modeIndexSet j_high_unmapped; if(j>=nl_r) mf_grid_t[t].getColParams().indexUnmap(j-nl_r,j_high_unmapped);
	SCFvectorPtr<typename GridA2Apolicies::FermionFieldType::FieldSiteType> rscf = Vgrid.getFlavorDilutedVect(j,j_high_unmapped,0,t_lcl);

	for(int a=0;a<2;a++)
	  for(int b=0;b<2;b++)
	    if(!lscf.isZero(a) && !rscf.isZero(b))
	      FLOPs_per_site += g5_FLOPs;
	FLOPs_per_site += siteFmat_FLOPs + s3_FLOPs + TransLeftTrace_FLOPs + reduce_FLOPs;
      }
    }
  }
  const typename GridA2Apolicies::FermionFieldType &mode0 = Wgrid.getMode(0);
  const int size_3d = mode0.nodeSites(0)*mode0.nodeSites(1)*mode0.nodeSites(2);
  double total_FLOPs = double(FLOPs_per_site) * double(size_3d) * double(ntests);

  printf("MF contract all t: Avg time new code %d iters: %g secs. Avg flops %g Gflops\n",ntests,total_time/ntests, total_FLOPs/total_time/1e9);      
}




CPS_END_NAMESPACE

#endif
