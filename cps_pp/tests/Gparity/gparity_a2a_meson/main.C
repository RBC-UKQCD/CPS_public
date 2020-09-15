//Generate meson contractions with A2A propagators for comparison with Daiqian
//NOTE: You will need to link against libfftw3 and libfftw3_threads

#include <alg/a2a/lanc_arg.h>
#include <alg/eigen/Krylov_5d.h>



#include<chroma.h>
//bfm headers
#include<actions/ferm/invert/syssolver_linop_cg_array.h>
#include<bfm.h>
#include<bfm_qdp.h>
#include<bfm_cg.h>
#include<bfm_mprec.h>


//cps headers
#include<alg/fermion_vector.h>
#include<alg/do_arg.h>
#include<alg/meas_arg.h>
#include<util/qioarg.h>
#include<util/ReadLatticePar.h>

//c++ classes
#include<sys/stat.h>
#include<util/qcdio.h>
//#include<fftw3.h>

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
#include <iostream>

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

#include <util/gparity_singletodouble.h>

#include <alg/a2a/alg_a2a.h>
#include <alg/a2a/MesonField.h>
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

#include<set>
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

#include <alg/prop_dft.h>

#include <string>
#include <sstream>
#include<omp.h>

#ifdef HAVE_BFM
#include <chroma.h>
#endif

using namespace std;
USING_NAMESPACE_CPS

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


void setup_bfmargs(bfmarg &dwfa, const BfmSolver &solver = DWF){
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
  dwfa.mass = toDouble(0.5);
  dwfa.Csw  = 0.0;
  dwfa.max_iter = 6000;
  dwfa.residual = 1e-10;
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

void lanczos_arg(LancArg &into, const Float &mass, const bool &precon){
  into.mass = mass;
  into.stop_rsd = 1e-06;
  into.qr_rsd = 1e-14; ///convergence of intermediate QR solves, defaults to 1e-14
  into.EigenOper = DDAGD;
  into.precon = precon; //also try this with true
  into.N_get = 4;///Want K converged vectors
  into.N_use = 11;///Dimension M Krylov space
  into.N_true_get = 4;//Actually number of eigen vectors you will get
  into.ch_ord = 3;///Order of Chebyshev polynomial
  into.ch_alpha = 5.5;///Spectral radius
  into.ch_beta = 0.5;///Spectral offset (ie. find eigenvalues of magnitude less than this)
  into.ch_sh = false;///Shifting or not
  into.ch_mu = 0;///Shift the peak
  into.lock = false;///Use locking transofrmation or not
  into.maxits =10000;///maxiterations
  into.fname = "Lanczos";
}


void lanczos_container_arg(LanczosContainerArg &into, const char* tag, const Float &mass,const bool &lanc_precon, const BfmSolver &solver){
  into.tag = strdup(tag);
  into.cg_max_iter = 10000;
  into.cg_residual = 1e-08;
  into.cg_precon_5d = 1;
  if(solver == HmCayleyTanh){
    into.cg_precon_5d = 0; //mobius uses 4d preconditioning
    into.mobius_scale = 2.0;
  }
  if(solver == DWF) into.solver = BFM_DWF;
  else if(solver == HmCayleyTanh) into.solver = BFM_HmCayleyTanh;
  else{
    printf("lanczos_container_arg(): unknown solver\n");
    exit(-1);
  }
  into.tbc = BND_CND_APRD;
  lanczos_arg(into.lanc_arg,mass,lanc_precon);
}



void a2a_arg(A2AArg &into, const char* lanc_tag, const int &flavor_dilution, const RandomType &rand_type, const int &nl, const bool &do_gauge_fix){
  into.lanczos_tag = strdup(lanc_tag);
  into.nl = nl;
  into.nhits = 1;
  into.rand_type = rand_type;
  into.src_width = 1;
  into.dilute_flavor = flavor_dilution;
  into.do_gauge_fix = do_gauge_fix;
}

void setup_gfix_args(FixGaugeArg &r){
  r.fix_gauge_kind = FIX_GAUGE_COULOMB_T;
  r.hyperplane_start = 0;
  r.hyperplane_step = 1;
  r.hyperplane_num = GJP.TnodeSites();
  r.stop_cond = 1e-06;
  r.max_iter_num = 6000;
}

#define SETUP_ARRAY(OBJ,ARRAYNAME,TYPE,SIZE)	\
  OBJ . ARRAYNAME . ARRAYNAME##_len = SIZE; \
  OBJ . ARRAYNAME . ARRAYNAME##_val = new TYPE [SIZE]

#define ELEM(OBJ,ARRAYNAME,IDX) OBJ . ARRAYNAME . ARRAYNAME##_val[IDX]

void setup_a2a_prop_args(PropagatorArg &p, const char* tag, const Float &mass, const char* lanc_tag, const int &flavor_dilution, const int &nl, const bool &do_gauge_fix){
  p.generics.tag = strdup(tag);
  p.generics.type = A2A_PROP_TYPE;
  p.generics.mass = mass;
  for(int i=0;i<4;i++) p.generics.bc[i] = GJP.Bc(i);
  
  SETUP_ARRAY(p,attributes,AttributeContainer,1);
  AttributeContainer &attr = ELEM(p,attributes,0);
  attr.type = A2A_ATTR;
  A2AAttrArg & a2a_args = attr.AttributeContainer_u.a2a_attr;
  a2a_arg(a2a_args,lanc_tag,flavor_dilution,UONE,nl,do_gauge_fix); //UONE random numbers, nl low modes 
}
    
static void global_coord(const int &site, int *into_vec){
  int rem = site;
  for(int i=0;i<4;i++){
    into_vec[i] = rem % GJP.NodeSites(i) + GJP.NodeCoor(i)*GJP.NodeSites(i);
    rem /= GJP.NodeSites(i);
  }
}

inline static bool ratio_diff_match(const Float &a, const Float &b, const Float &tolerance){
  Float rat = a/b;
  if(rat < 0.0) return false; //opposite signs
  else return fabs(rat-1.0) <= tolerance;
}

static void set_p(Float into[3], char sgn){ //units of pi/L
  for(int i=0;i<3;i++){
    if(GJP.Bc(i) == BND_CND_GPARITY)
      into[i] = (sgn == '+' ? 0.5 : -0.5);
    else into[i] = 0.0;
  }  
}

static void setup_ktopipi_args(ContractionTypeKtoPiPi &into, const char *prop_L, const char* prop_H, char p_pi1_sign, char p_pi2_sign){
  into.prop_L = strdup(prop_L);
  into.prop_H = strdup(prop_H);
  
  set_p(into.p_qpi1.p1,p_pi1_sign);
  set_p(into.p_qpi1.p2,p_pi1_sign);

  set_p(into.p_qpi2.p1,p_pi2_sign);
  set_p(into.p_qpi2.p2,p_pi2_sign);

  set_p(into.p_qK,'+');

  into.gparity_use_transconv_props = 1;
  into.pion_source.type = BOX_3D_SMEARING;
  into.pion_source.A2ASmearing_u.box_3d_smearing.side_length = GJP.Xnodes()*GJP.XnodeSites();
  
  into.kaon_source.type = BOX_3D_SMEARING;
  into.kaon_source.A2ASmearing_u.box_3d_smearing.side_length = GJP.Xnodes()*GJP.XnodeSites();
  
  into.t_sep_pi_k = 2;
  into.t_sep_pion = 1;
  
  into.file = "type1.dat";
}

static void set_p(Float into[3], std::vector<int> q){ //input integers in units of pi/2L
  into[0] = 0.5*q[0];
  into[1] = 0.5*q[1];
  into[2] = 0.5*q[2];
}


static void setup_ktopipi_multimom_args(ContractionTypeKtoPiPi &into, const char *prop_L, const char* prop_H, std::vector<int> p_pi1q1, std::vector<int> p_pi1q2,  std::vector<int> p_pi2q1, std::vector<int> p_pi2q2){
  into.prop_L = strdup(prop_L);
  into.prop_H = strdup(prop_H);
  
  set_p(into.p_qpi1.p1, p_pi1q1);
  set_p(into.p_qpi1.p2, p_pi1q2);

  set_p(into.p_qpi2.p1, p_pi2q1);
  set_p(into.p_qpi2.p2, p_pi2q2);

  set_p(into.p_qK,'+');

  into.gparity_use_transconv_props = 1;
  into.pion_source.type = BOX_3D_SMEARING;
  into.pion_source.A2ASmearing_u.box_3d_smearing.side_length = GJP.Xnodes()*GJP.XnodeSites();
  
  into.kaon_source.type = BOX_3D_SMEARING;
  into.kaon_source.A2ASmearing_u.box_3d_smearing.side_length = GJP.Xnodes()*GJP.XnodeSites();
  
  into.t_sep_pi_k = 2;
  into.t_sep_pion = 1;
  
  into.file = "type1.dat";
}





class Testing{
public:
  static void test_prod1(Gparity_KtoPiPi &gpcon, const ContractionTypeKtoPiPi &args, Lattice &lat){
    if(!gpcon.setup_called) gpcon.setup(args,lat);

    int x_op_loc = 0;
    int tpi2 = 1;
    //Form SpinColorFlavorMatrix prod1 = vL_i(\vec xop, top) [[\sum_{\vec xpi2} wL_i^dag(\vec xpi2, tpi2) S2 vL_j(\vec xpi2, tpi2)]] wL_j^dag(\vec xop,top)
    SpinColorFlavorMatrix prod1;
    MesonField2::contract_vleft_wright(prod1, *gpcon.prop_L, x_op_loc, *gpcon.prop_L, x_op_loc, *gpcon.wdagL_S2_vL_pi2, tpi2);

    printf("Prod1  (0,0,0,0,0,0): %.12le %.12le    (1,1,1,1,1,1): %.12le %.12le\n", 
	   prod1(0,0,0,0,0,0).real(), prod1(0,0,0,0,0,0).imag(),
	   prod1(1,1,1,1,1,1).real(), prod1(1,1,1,1,1,1).imag()
	   );

  }

  static void test_conLLLH(Gparity_KtoPiPi &gpcon, const ContractionTypeKtoPiPi &args, Lattice &lat){
    if(!gpcon.setup_called) gpcon.setup(args,lat);

    int t_K = 0;
    int t_pi1 = 1;

    RangeSpecificT sett2tK(t_K);
    MesonField2 con_LLLH;
    MesonField2::combine_mf_wv_ww(con_LLLH, *gpcon.wdagL_S2_vL_pi1, *gpcon.wdagL_wH, sett2tK);

    
    FILE *fp = fopen("con_LLLH.ck","w");
    for(int i=0;i<con_LLLH.get_size(MesonField2::Left);i++){
      for(int j=0;j<con_LLLH.get_size(MesonField2::Right);j++){
	std::complex<double> con = *con_LLLH.mf_val(i,j,t_pi1);
	fprintf(fp,"%d %d %.12le %.12le\n",i,j,4.0*con.real(),4.0*con.imag());
      }
    }
    fclose(fp);

    std::complex<double> con00 = *con_LLLH.mf_val(0,0,t_pi1);
    std::complex<double> con35 = *con_LLLH.mf_val(3,5,t_pi1);

    printf("conLLLH  (i,j) = (0,0) : %.12le %.12le   and  (i,j) = (3,5) : %.12le %.12le\n", con00.real(),con00.imag(), con35.real(), con35.imag() );


    int x_op_loc = 0;

    SpinColorFlavorMatrix prod2; 
    Testing::contract_vleft_vright_testcopy(prod2, *gpcon.prop_L, x_op_loc, *gpcon.prop_H, x_op_loc, con_LLLH, t_pi1, t_K);

    printf("Prod2 pre-g5 (0,0,0,0,0,0): %.12le %.12le    (1,1,1,1,1,1): %.12le %.12le\n", 
	   prod2(0,0,0,0,0,0).real(), prod2(0,0,0,0,0,0).imag(),
	   prod2(1,1,1,1,1,1).real(), prod2(1,1,1,1,1,1).imag()
	   );

    prod2.gr(-5);

    printf("Prod2  (0,0,0,0,0,0): %.12le %.12le    (1,1,1,1,1,1): %.12le %.12le\n", 
	   prod2(0,0,0,0,0,0).real(), prod2(0,0,0,0,0,0).imag(),
	   prod2(1,1,1,1,1,1).real(), prod2(1,1,1,1,1,1).imag()
	   );


  }


  static void contract_vleft_vright_testcopy(SpinColorFlavorMatrix &result, 
						 A2APropbfm &prop_left,  const int &x,
						 A2APropbfm &prop_right, const int &z,
						 MesonField2 & mf,  const int &t, const int &t_vright){

    result *= 0.0;
    for(int j=0;j<mf.get_size(MesonField2::Right);j++){
      for(int i=0;i<mf.get_size(MesonField2::Left);i++){
	int I = prop_left.idx_v(i,t);
	int J = prop_right.idx_v(j,t_vright);

	printf("i_size = %d,  j_size = %d,  i=%d, j=%d -> I=%d, J=%d\n",
	       mf.get_size(MesonField2::Left),mf.get_size(MesonField2::Right),
	       i,j,I,J);
	       

	std::complex<double>* M = mf.mf_val(i,j,t);
      
	printf("M(i=%d,j=%d,t=%d) = %.12le , %.12le\n",i,j,t, M->real(), M->imag());
	
	for(int fl = 0; fl<2; fl++){
	  std::complex<double>* left_v =  prop_left.get_v(I, x,fl);

	  printf("left_v(I=%d,x=%d,fl=%d) = %.12le , %.12le\n",I,x,fl, left_v->real(), left_v->imag());

	  for(int fr = 0; fr<2; fr++){
	    std::complex<double>* right_v =  prop_right.get_v(J, z,fr);

	    printf("right_v(J=%d,z=%d,fr=%d) = %.12le , %.12le\n",J,z,fr, right_v->real(), right_v->imag());

	    for(int cr=0;cr<3;cr++){
	      for(int sr=0;sr<4;sr++){
		int sc_off_R = cr+3*sr;

		for(int cl=0;cl<3;cl++){
		  for(int sl=0;sl<4;sl++){
		    int sc_off_L = cl+3*sl;
		  
		    result(sl,cl,fl,sr,cr,fr) += left_v[sc_off_L] * (*M) * std::conj(right_v[sc_off_R]);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }



  static void type1_propHg5conj_test(Gparity_KtoPiPi &gpcon, const ContractionTypeKtoPiPi &args, Lattice &lat){
    if(!gpcon.setup_called) gpcon.setup(args,lat);
    
    int t_sep_pi_k = 2;
    int t_sep_pion = 1;
    
    int tkspec = 0;
    
    static const int &c_start = 1;
    static const int &n_con = 6;

    int t_size = GJP.Tnodes()*GJP.TnodeSites();

    std::vector< MesonField2 > con_LLLH(t_size);
    for(int t_K = 0 ; t_K < t_size; t_K++){
      if(t_K != tkspec) continue;
      RangeSpecificT sett2tK(t_K);
      MesonField2::combine_mf_wv_ww(con_LLLH[t_K], *gpcon.wdagL_S2_vL_pi1, *gpcon.wdagL_wH, sett2tK);
    }

    int x_op_loc = 0;

    int t_op_loc = x_op_loc / (GJP.VolNodeSites()/GJP.TnodeSites());
    int t_op_glob = t_op_loc + GJP.TnodeCoor()*GJP.TnodeSites();

    //Average over all tK
    int tK = tkspec;

    //We need to average over choices  1) tpi1 = tK + t_sep_pi_k, tpi2 = tpi1 + t_sep_pion  and  2) tpi2 = tK + t_sep_pi_k, tpi1 = tpi2 + t_sep_pion 
    int tpi1_c[2] = { (tK + t_sep_pi_k) % t_size              , (tK + t_sep_pi_k + t_sep_pion) % t_size };
    int tpi2_c[2] = { (tK + t_sep_pi_k + t_sep_pion) % t_size , (tK + t_sep_pi_k) % t_size              };

    int c = 0;
    int mumax = 3;

    printf("t_K = %d,  t_pi1 = %d,  t_pi2 = %d\n",tK,tpi1_c[c],tpi2_c[c]);

    //Form SpinColorFlavorMatrix prod1 = vL_i(\vec xop, top ; tpi2) [\sum_{\vec xpi2} wL_i^dag(\vec xpi2, tpi2) S2 vL_j(\vec xpi2, tpi2; top)] wL_j^dag(\vec xop,top)
    SpinColorFlavorMatrix prod1;
    MesonField2::contract_vleft_wright(prod1, *gpcon.prop_L, x_op_loc, *gpcon.prop_L, x_op_loc, *gpcon.wdagL_S2_vL_pi2, tpi2_c[c]);

    printf("Prod1  (0,0,0,0,0,0): %.12le %.12le    (1,1,1,1,1,1): %.12le %.12le\n", 
	   prod1(0,0,0,0,0,0).real(), prod1(0,0,0,0,0,0).imag(),
	   prod1(1,1,1,1,1,1).real(), prod1(1,1,1,1,1,1).imag()
	   );

    //Form SpinColorFlavorMatrix prod2 = vL_i(\vec xop, top ; tpi1) con_LLLH_{ij}(tpi1)[tK] vH_j^dag(\vec xop, top; t_K) \gamma^5
    SpinColorFlavorMatrix prod2; 
    MesonField2::contract_vleft_vright(prod2, *gpcon.prop_L, x_op_loc, *gpcon.prop_H, x_op_loc, con_LLLH[tK], tpi1_c[c], tK);
    prod2.gr(-5);

    printf("Prod2  (0,0,0,0,0,0): %.12le %.12le    (1,1,1,1,1,1): %.12le %.12le\n", 
	   prod2(0,0,0,0,0,0).real(), prod2(0,0,0,0,0,0).imag(),
	   prod2(1,1,1,1,1,1).real(), prod2(1,1,1,1,1,1).imag()
	   );

    int g1idx = 0; //M(0,V)
    int g2idx = 1; //M(0,A)

#define USE_DAIQIANS_NEW_DEFINITIONS  

    //Note Daiqian renamed some of his C functions such that the type-II and type-I forms look more similar
#ifdef USE_DAIQIANS_NEW_DEFINITIONS  
    int dd[6] = {1,3,4,6,2,5};
#else
    int dd[6] = {1,2,3,4,5,6};
#endif	 
    
    std::complex<double> C[7];

    for(int mu=0;mu<=mumax;mu++){
      C[dd[0]] += 0.5*Trace( gpcon.Gamma[g1idx][mu], prod1 ) * Trace( gpcon.Gamma[g2idx][mu], prod2 );
      C[dd[1]] += 0.5*( SpinFlavorTrace( gpcon.Gamma[g1idx][mu], prod1 ) * SpinFlavorTrace( gpcon.Gamma[g2idx][mu], prod2 ) ).Trace();
      C[dd[2]] += 0.5*Trace(gpcon.Gamma[g1idx][mu]*prod1 , gpcon.Gamma[g2idx][mu]*prod2);
      C[dd[3]] += 0.5*Trace( ColorTrace( gpcon.Gamma[g1idx][mu], prod1 ), ColorTrace( gpcon.Gamma[g2idx][mu], prod2 ) );
      C[dd[4]] += 0.5*( SpinFlavorTrace( gpcon.Gamma[g1idx][mu], prod1 ) * Transpose(SpinFlavorTrace( gpcon.Gamma[g2idx][mu], prod2 )) ).Trace();
      C[dd[5]] += 0.5*Trace(gpcon.Gamma[g1idx][mu]*prod1 , ColorTranspose(gpcon.Gamma[g2idx][mu]*prod2) );
    }

    printf("C1 %.12le %.12le\n",C[1].real(),C[1].imag());
    printf("C2 %.12le %.12le\n",C[2].real(),C[2].imag());
    printf("C3 %.12le %.12le\n",C[3].real(),C[3].imag());
    printf("C4 %.12le %.12le\n",C[4].real(),C[4].imag());
    printf("C5 %.12le %.12le\n",C[5].real(),C[5].imag());
    printf("C6 %.12le %.12le\n",C[6].real(),C[6].imag());

  }


  static void type4_test(Gparity_KtoPiPi &gpcon, const ContractionTypeKtoPiPi &args, Lattice &lat){
    if(!gpcon.setup_called) gpcon.setup(args,lat);
    
    int t_sep_pi_k = 2;
    int t_sep_pion = 1;
    
    int tkspec = 0;
    std::vector<CorrelationFunction> into;

    static const int &c_start = 23;
    static const int &n_con = 10;
    gpcon.setup_resultvec(n_con,c_start,into);

    int t_size = GJP.Tnodes()*GJP.TnodeSites();

    //A separate meson 'blob'
    //blob = \sum_{\vec x_pi1,\vec x_pi2} -0.5 * Tr( \prop^L(x_pi1,x_pi2) S_2 \prop^L(x_pi2,x_pi1) S_2 )
    //     = \sum_{\vec x_pi1,\vec x_pi2} -0.5 * Tr( [[ wL^dag(x_pi2) S_2 vL(x_pi2) ]] [[ wL^dag(x_pi1) S_2 vL(x_pi1) ]] )
    //     = \sum_{\vec x_pi1,\vec x_pi2} -0.5 * Tr( [[ wL^dag(x_pi1) S_2 vL(x_pi1) ]] [[ wL^dag(x_pi2) S_2 vL(x_pi2) ]] )
    //we need both  t_pi2 = (t_pi1 + sep) % T (index 0 of array)   and t_pi2 = (t_pi1 - sep + T) % T   (index 1 of array)  

    RangeT1plusDelta fix_tpi2_plus_sep(t_sep_pion);
    RangeT1plusDelta fix_tpi2_minus_sep(-t_sep_pion);
  
    CorrelationFunction blob_tpi1_plus("blob",1,CorrelationFunction::THREADED);
    MesonField2::contract_specify_t2range( *gpcon.wdagL_S2_vL_pi1, *gpcon.wdagL_S2_vL_pi2, 0, fix_tpi2_plus_sep, blob_tpi1_plus);
  
    CorrelationFunction blob_tpi1_minus("blob",1,CorrelationFunction::THREADED);
    MesonField2::contract_specify_t2range( *gpcon.wdagL_S2_vL_pi1, *gpcon.wdagL_S2_vL_pi2, 0, fix_tpi2_minus_sep, blob_tpi1_minus);

    //Sum over t_pi1 and average over pion permutations
    std::complex<double> blob(0.0);
    for(int t_pi1=0;t_pi1<t_size;t_pi1++){
      blob += blob_tpi1_plus(0,t_pi1);
      blob += blob_tpi1_minus(0,t_pi1);
    }
    blob *= -0.5*0.5;  //one factor of 1/2 from blob definition, one from average over permutations

    if(!UniqueID()){
      printf("Pion blob %.12le %.12le\n", blob.real(),blob.imag());
    }

#if 0
    //Each contraction of this type is made up of different trace combinations of two objects (below for simplicity we ignore the fact that the two vectors in the meson fields are allowed to vary in position relative to each other):
    //1) \prop^L(x_op,x_K) \gamma^5 \prop^H(x_K,x_op)
    //   we use g5-hermiticity on the strange prop
    //  \prop^L(x_op,x_K)  [ \prop^H(x_op,x_K) ]^dag \gamma^5
    //= vL(x_op) [[ wL^dag(x_K) wH(x_K) ]] vH^dag(x_op) \gamma_5
  
    //2) \prop^L(x_op,x_op)   OR   \prop^H(x_op,x_op)
  
    //Loop over xop
    int n_threads = bfmarg::threads;
    omp_set_num_threads(n_threads);
#pragma omp parallel for 
    for(int x_op_loc = 0; x_op_loc < GJP.VolNodeSites(); x_op_loc++){
      int me = omp_get_thread_num();

      int t_op_loc = x_op_loc / (GJP.VolNodeSites()/GJP.TnodeSites());
      int t_op_glob = t_op_loc + GJP.TnodeCoor()*GJP.TnodeSites();

      //prod2 = \prop^L(x_op,x_op)   OR   \prop^H(x_op,x_op)   =    v(x_op) w^dag(x_op)
      SpinColorFlavorMatrix prod2_L, prod2_H;
      MesonField2::contract_vw(prod2_L,*prop_L,x_op_loc,*prop_L,x_op_loc);
      MesonField2::contract_vw(prod2_H,*prop_H,x_op_loc,*prop_H,x_op_loc);

      //Average over all tK
      for(int tK = 0; tK < t_size; tK++){
	if(tK_vals != NULL && !int_in_vec(tK,*tK_vals) ) continue;

	//prod1 = vL(x_op) [[ wL^dag(x_K) wH(x_K) ]] vH^dag(x_op) \gamma_5
	SpinColorFlavorMatrix prod1;
	MesonField2::contract_vleft_vright(prod1, *prop_L, x_op_loc, *prop_H, x_op_loc, *wdagL_wH , tK, tK);
	prod1.gr(-5);

	for(int g1idx = 0; g1idx < 4; ++g1idx){
	  for(int g2idx = 0; g2idx < 4; ++g2idx){
	    std::complex<double>* corrs[n_con];
	    for(int cc=0;cc<n_con;cc++) corrs[cc] = & into[result_map(cc, g1idx, g2idx)](me,0,t_op_glob);
	  
	    //I is the index of the contraction in the range [1,32]
#define C(I) *corrs[I-c_start] 
	  
	    //Sum over mu for each \Gamma_1 and \Gamma_2	      
	    for(int mu = 0; mu < 4; ++mu){ //remember to include 0.5 from average over the swapping of the pion source timeslices
	      C(23) += blob * Trace( Gamma[g1idx][mu], prod1 ) * Trace( Gamma[g2idx][mu], prod2_L );
	      C(24) += blob * ( SpinFlavorTrace( Gamma[g1idx][mu], prod1 ) * Transpose(SpinFlavorTrace( Gamma[g2idx][mu], prod2_L )) ).Trace();
	      C(25) += blob * ( SpinFlavorTrace( Gamma[g1idx][mu], prod1 ) * SpinFlavorTrace( Gamma[g2idx][mu], prod2_L ) ).Trace();  
	      C(26) += blob * Trace(Gamma[g1idx][mu]*prod1 , Gamma[g2idx][mu]*prod2_L);
	      C(27) += blob * Trace(Gamma[g1idx][mu]*prod1 , ColorTranspose(Gamma[g2idx][mu]*prod2_L) );
	      C(28) += blob * Trace( ColorTrace( Gamma[g1idx][mu], prod1 ), ColorTrace( Gamma[g2idx][mu], prod2_L ) );
	    
	      C(29) += blob * Trace( Gamma[g1idx][mu], prod1 ) * Trace( Gamma[g2idx][mu], prod2_H );
	      C(30) += blob * ( SpinFlavorTrace( Gamma[g1idx][mu], prod1 ) * SpinFlavorTrace( Gamma[g2idx][mu], prod2_H ) ).Trace();
	      C(31) += blob * Trace(Gamma[g1idx][mu]*prod1 , Gamma[g2idx][mu]*prod2_H);
	      C(32) += blob * Trace( ColorTrace( Gamma[g1idx][mu], prod1 ), ColorTrace( Gamma[g2idx][mu], prod2_H ) );
	    }
	  
#undef C
	  }
	}
      }
    }//end of site loop	
    
    for(int i=0;i<into.size();i++)
      into[i].sumLattice();
  
#endif
  }

};

std::vector<int> three_vec(int i, int j, int k){
  std::vector<int> out(3);
  out[0] = i; out[1] = j; out[2] = k;
  return out;
}
std::vector<int> minus_vec(const std::vector<int> &in){
  std::vector<int> out(in);
  for(int i=0;i<out.size();i++) out[i] = -out[i];
  return out;
}
std::vector<int> permute(const std::vector<int> &in, int idx1, int idx2){
  //pairwise swap
  std::vector<int> out(in);
  out[idx1] = in[idx2];
  out[idx2] = in[idx1];
  return out;
}

void compute_type(std::vector<CorrelationFunction> &into, const int &type, const int &tk, const int &t_sep_pi_k, const int &tsep_pion, Gparity_KtoPiPi &gpcon){
  std::vector<int> tkv(1,tk);
	
  if(type == 1){
    gpcon.type1(t_sep_pi_k, tsep_pion, into, &tkv);
  }else if(type == 2){
    gpcon.type2(t_sep_pi_k, tsep_pion, into, &tkv);
  }else if(type == 3){
    gpcon.type3(t_sep_pi_k, tsep_pion, into, &tkv);    
  }else if(type == 4){
    gpcon.type4(t_sep_pi_k, tsep_pion, into, &tkv);    
  }else if(type == 5){
    gpcon.psvertex_type3(t_sep_pi_k, tsep_pion, into, &tkv); 
  }else if(type == 6){
    gpcon.psvertex_type4(t_sep_pi_k, tsep_pion, into, &tkv); 
  }else{
    printf("compute_type(..) : Invalid contraction type\n");
    exit(-1);
  }
}

void daiqian_order_result(std::vector< std::vector< std::vector<cps::Complex> > > &ordered, std::vector<CorrelationFunction> &orig_output, const int &type, const int &tk){
  //Daiqian prints out
  //<t_K> <(t_op-t_K)%T> <C1(0,V,0,A)> <C1(0,A,0,V)> <C1(0,V,1,A)> <C1(0,A,1,V)> <C1(1,V,0,A)> <C1(1,A,0,V)> <C1(1,V,1,A)> <C1(1,A,1,V)> <C2.............
  //where each C contains both the real and imaginary parts
  //And prints 0 when  t_op - t_K outside the range  0 < t_op - t_K < delta    where delta is the pi-K separation (2 here)
  //const static std::string GammaNames[4] = { "M_{0,V}","M_{0,A}", "M_{1,V}", "M_{1,A}" };
  const int _0V(0);
  const int _0A(1);
  const int _1V(2);
  const int _1A(3);
  const int delta(2);
  const int Lt( GJP.Tnodes()*GJP.TnodeSites() );

  int ncon =  orig_output.size()/4/4 ;
  if(type > 4) ncon = 1;

  ordered.resize(Lt); //[tt][c][gammacomb]

  //We want the second time index to be tt= (top-tk+Lt) % Lt in order
  for(int top=0;top<tsize;++top){
    int tt = (top-tk+Lt )%Lt;
    ordered[tt].resize(ncon);

    for(int c=0;c<ncon;c++){
      int size = 8;
      if(type > 4) size = 2;

      std::vector<cps::Complex> &v = ordered[tt][c];
      v.resize( size,cps::Complex(0.0,0.0) );
	    
      if(type <= 4){
	int m[8] = { Gparity_KtoPiPi::result_map(c,_0V,_0A), Gparity_KtoPiPi::result_map(c,_0A,_0V),
		     Gparity_KtoPiPi::result_map(c,_0V,_1A), Gparity_KtoPiPi::result_map(c,_0A,_1V),
		     Gparity_KtoPiPi::result_map(c,_1V,_0A), Gparity_KtoPiPi::result_map(c,_1A,_0V),
		     Gparity_KtoPiPi::result_map(c,_1V,_1A), Gparity_KtoPiPi::result_map(c,_1A,_1V) };
	      
	v[0] = orig_output[m[0]](0,top);
	v[1] = orig_output[m[1]](0,top);
	      
	v[2] = orig_output[m[2]](0,top);
	v[3] = orig_output[m[3]](0,top);
	      
	v[4] = orig_output[m[4]](0,top);
	v[5] = orig_output[m[5]](0,top);
	      
	v[6] = orig_output[m[6]](0,top);
	v[7] = orig_output[m[7]](0,top);
      }else{
	v[0] = orig_output[0](0,top);
	v[1] = orig_output[1](0,top);
      }

    }
  }
}


using namespace Chroma;
using namespace cps;

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

int cout_time(char *);
void ReadGaugeField(const MeasArg &meas_arg);
void bfm_init(bfm_evo<double> &dwf,double mq);

#define Printf if(!UniqueID()) printf

int main (int argc,char **argv )
{
  Start(&argc, &argv);

#ifdef HAVE_BFM
  Chroma::initialize(&argc,&argv);
#endif

  CommandLine::is(argc,argv);

  bool gparity = false;
  bool gparity_dirs[3] = {false,false,false};

  int arg0 = CommandLine::arg_as_int(0);
  if(arg0<4){
    Printf("Number of G-parity directions: %d\n",arg0);
    for(int i=0;i<arg0;i++)
      gparity_dirs[i] = true;
    gparity=true;
  }else{
    Printf("Expect 0,1,2 or 3 G-parity directions\n");
    exit(-1);
  }

  bool save_config(false);
  bool load_config(false);
  bool load_lrg(false);
  bool save_lrg(false);
  char *load_config_file = NULL;
  char *save_config_file = NULL;
  char *save_lrg_file;
  char *load_lrg_file;
  bool verbose(false);
  bool unit_gauge(false);

  int size[] = {2,2,2,2,2};

  int dilute_flavor = 1;

  int i=2;
  while(i<argc){
    char* cmd = argv[i];  
    if( strncmp(cmd,"-save_config",15) == 0){
      if(i==argc-1){ Printf("-save_config requires an argument\n"); exit(-1); }
      save_config=true;
      save_config_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-load_config",15) == 0){
      if(i==argc-1){ Printf("-save_config requires an argument\n"); exit(-1); }
      load_config=true;
      load_config_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-latt",10) == 0){
      if(i>argc-6){
	Printf("Did not specify enough arguments for 'latt' (require 5 dimensions)\n"); exit(-1);
      }
      size[0] = CommandLine::arg_as_int(i); //CommandLine ignores zeroth input arg (i.e. executable name)
      size[1] = CommandLine::arg_as_int(i+1);
      size[2] = CommandLine::arg_as_int(i+2);
      size[3] = CommandLine::arg_as_int(i+3);
      size[4] = CommandLine::arg_as_int(i+4);
      i+=6;
    }else if( strncmp(cmd,"-load_lrg",15) == 0){
      if(i==argc-1){ Printf("-load_lrg requires an argument\n"); exit(-1); }
      load_lrg=true;
      load_lrg_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-save_lrg",15) == 0){
      if(i==argc-1){ Printf("-save_lrg requires an argument\n"); exit(-1); }
      save_lrg=true;
      save_lrg_file = argv[i+1];
      i+=2;  
    }else if( strncmp(cmd,"-verbose",15) == 0){
      verbose=true;
      i++;
    }else if( strncmp(cmd,"-dont_dilute_flavor",15) == 0){
      Printf("Running without flavor dilution\n");
      dilute_flavor = 0;
      i++;
    }else if( strncmp(cmd,"-unit_gauge",15) == 0){
      unit_gauge=true;
      i++;
    }else{
      Printf("Unrecognised argument: %s\n",cmd);
      exit(-1);
    }
  }
  
  Printf("Lattice size is %d %d %d %d\n",size[0],size[1],size[2],size[3],size[4]);

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

  BndCndType* bc[3] = { &do_arg.x_bc, &do_arg.y_bc, &do_arg.z_bc };
  for(int i=0;i<3;i++)
    if(gparity_dirs[i]) *bc[i] = BND_CND_GPARITY;

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
    LatRanGen LRGbak = LRG; //don't let the gauge generation modify the RNG
    printf("Creating gauge field\n");
    if(!unit_gauge) lattice->SetGfieldDisOrd();
    else lattice->SetGfieldOrd();
    LRG = LRGbak;
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
  lattice->FixGaugeAllocate(FIX_GAUGE_COULOMB_T);
  lattice->FixGauge(1e-06,2000);

  cps_qdp_init(&argc,&argv);
  
  // LRG.AssignGenerator(0,0);
  // printf("INIT RAND %f\n",LRG.Urand(FOUR_D));


#if TARGET == BGQ
  omp_set_num_threads(64);
#else
  omp_set_num_threads(1);
#endif

  //Light quark
  Float mass = 0.1;

  #define DO_GAUGE_FIX

#ifdef DO_GAUGE_FIX
  bool do_gauge_fix = true;
#else
  bool do_gauge_fix = false;
#endif

  LanczosContainerArg lanc_arg;
  lanczos_container_arg(lanc_arg, "lanc", mass, true, HmCayleyTanh);

  PropManager::addLanczos(lanc_arg);

  PropagatorArg prop_arg;
  setup_a2a_prop_args(prop_arg,"a2aprop", mass,"lanc",dilute_flavor,4,do_gauge_fix); //4 low modes

  PropManager::addProp(prop_arg);

  //Strange quark (don't need another Lanczos instance; just use 0 low modes)
  Float mass_s = 0.3;
  PropagatorArg prop_arg_s;
  setup_a2a_prop_args(prop_arg_s,"a2aprop_s", mass_s,"lanc",dilute_flavor,0,do_gauge_fix); //0 low modes

  PropManager::addProp(prop_arg_s);

  //Calculate the A2A props
  PropManager::calcProps(*lattice);

  {
    LatRanGen LRGbak(LRG);
    LRG.AssignGenerator(0,0);
    printf("RAND POST PROP GEN %f\n",LRG.Urand(FOUR_D));
    LRG = LRGbak;
  }

  //Do the first test with the pion 2-point function
  A2APropbfm &a2a_prop = PropManager::getProp("a2aprop").convert<A2ApropContainer>().getProp(*lattice);
  
  //Print the stochastic field for the high modes
  int whsize = 2 * GJP.VolNodeSites();
  Float* wh = (Float*)a2a_prop.get_wh();
  
  {
    FILE *fp = fopen("wh0.ck","w");
    for(int i=0;i<whsize;i++){
      fprintf(fp,"%.16f\n",wh[i]);
    }
    fclose(fp);
  }
  //Do the v fields
  for(int vidx=0;vidx< a2a_prop.v_modes();vidx++){
    std::ostringstream fn;
    fn << "v" << vidx << ".ck";
    
    Float* v = (Float*)a2a_prop.get_v(vidx);
    int vsize = 2 * 4*3*2 * GJP.VolNodeSites();
    
    FILE *fp=fopen(fn.str().c_str(),"w");

    for(int i=0;i<vsize;i++){
      fprintf(fp,"%.16f\n",v[i]);
      //if(vidx==3) v[i] = -v[i];
    }
    fclose(fp);
  }
  //Do the w low-mode fields
  for(int widx=0;widx< a2a_prop.get_nl();widx++){
    std::ostringstream fn;
    fn << "w" << widx << ".ck";
    
    Float* w = (Float*)a2a_prop.get_wl(widx);
    int wsize = 2 * 4*3*2 * GJP.VolNodeSites();
    
    FILE *fp=fopen(fn.str().c_str(),"w");

    for(int i=0;i<wsize;i++){
      fprintf(fp,"%.16f\n",w[i]);
      //if(widx==3) w[i] = -w[i];
    }
    fclose(fp);
  }

  #define USE_TRANSCONV_FIELDS
  #define IMPOSE_MOMENTUM

  if(0){
    //Do contraction using Daiqian's translationally covariant field source
    //\sum_{x,x',y,y'} e^{-ipx}e^{ipy}e^{-ipx'}e^{ipy'} 0.5 tr{  G_{x,y} (\sigma_3+i\sigma_1)g5 G_{y',x'} g5(\sigma_3-i\sigma_1) } 
    //\sum_{x,x',y,y'} e^{-ipx}e^{ipy}e^{-ipx'}e^{ipy'} 0.5 tr{  v(x)w^dag(y) (\sigma_3+i\sigma_1)g5 v(y') w^dag(x') g5(\sigma_3-i\sigma_1) } 
    //\sum_{x,x',y,y'} e^{-ipx}e^{ipy}e^{-ipx'}e^{ipy'} 0.5 tr{  w^dag(y) (\sigma_3+i\sigma_1)g5 v(y')    w^dag(x') g5(\sigma_3-i\sigma_1)v(x) } 

    //x is the source coordinate
    //Code does   w^dag_i(tsrc) .. v_j(tsrc;tsnk)  w^dag_j(tsnk) .. v_i(tsnk;tsrc)
    //which is opposite order to above. Change ordering; 
    //\sum_{x,x',y,y'} e^{-ipx}e^{ipy}e^{-ipx'}e^{ipy'} 0.5 tr{ [ w^dag(x') g5(\sigma_3-i\sigma_1)v(x) ] [ w^dag(y) (\sigma_3+i\sigma_1)g5 v(y') ] } 

    SpinColorFlavorMatrix s3g5(gamma5,sigma3);
    SpinColorFlavorMatrix s1g5(gamma5,sigma1);
    static std::complex<Float> _i(0.0,1.0);

    SpinColorFlavorMatrix lmat = s3g5*0.5 - s1g5*_i*0.5;
    SpinColorFlavorMatrix rmat = s3g5 + s1g5*_i;

    MFqdpMatrix mf_mat_left(MFstructure::W, MFstructure::V, true, false, lmat);
    MFqdpMatrix mf_mat_right(MFstructure::W, MFstructure::V, true, false, rmat);
  
    MFBasicSource src(MFBasicSource::BoxSource, (double)(GJP.XnodeSites()*GJP.Xnodes()) );
  
    Float p[3] = {0.0,0.0,0.0};
    for(int i=0;i<3;++i) 
      if(GJP.Bc(i) == BND_CND_GPARITY)
	p[i] = 0.5; //units of pi/L
  
    const Float mp[3] = {-p[0],-p[1],-p[2]};

    //Note, FT conventions are  e^{-ipx} for forwards FFT into p-space with momentum +p
    //Hence x and x' have momentum +p and likewise y and y' have momentum -p

    a2a_prop.set_v_momentum(p);
    a2a_prop.set_wdag_momentum(p);
    a2a_prop.fft_vw();

    MesonField2 mf_left(a2a_prop,a2a_prop,mf_mat_left,src);
  
    a2a_prop.set_v_momentum(mp);
    a2a_prop.set_wdag_momentum(mp);
    a2a_prop.fft_vw();
 
    MesonField2 mf_right(a2a_prop,a2a_prop,mf_mat_right,src);
  
    CorrelationFunction result("result",1);
  
    MesonField2::contract(mf_left,mf_right,result);

    // if(!UniqueID()){
    //   printf("Result:\n");
    //   for(int t=0;t<GJP.TnodeSites()*GJP.Tnodes();++t){
    // 	printf("%d %.12le %.12le\n",t,result(0,t).real(),result(0,t).imag());
    //   }
    // }

    //Also do the same but with explicit source timeslice
    if(!UniqueID()) printf("Result with source and sink timeslice explicit:\n");

    for(int x4=0;x4<GJP.TnodeSites()*GJP.Tnodes();x4++){
      CorrelationFunction result("result",1);
      MesonField2::contract_specify_tsrc(mf_left,mf_right,0,x4,result);
    
      if(!UniqueID()){
	for(int y4=0;y4<GJP.TnodeSites()*GJP.Tnodes();++y4){
	  printf("%d %d %.12le %.12le\n",x4,y4,result(0,y4).real(),result(0,y4).imag());
	}
      }

    }
  }
  int t_size = GJP.Tnodes()*GJP.TnodeSites();
  std::complex<double> twopoint_results[t_size][t_size];
  if(0){
    //Do contraction using Daiqian's translationally covariant field source. Here we repeat the above but using the inbuilt code that makes the FFT propagators translationially covariant
    //\sum_{x,x',y,y'} e^{-ipx}e^{ipy}e^{-ipx'}e^{ipy'} 0.5 * 4 * tr{ [ w^dag(x') g5 \sigma_3 v(x) ] [ w^dag(y) g5 \sigma_3 v(y') ] } 
    //Note the extra factor of 4 which arises because previously we measured with   \sigma_3 (1 \pm \sigma_2) in each meson field  
    //and now we measure with 1/2(1 \mp \sigma_2)\sigma_3 1/2(1 \pm \sigma_2) = \sigma_3 1/2(1 \pm \sigma_2) in each meson field

    SpinColorFlavorMatrix s3g5(gamma5,sigma3);
    SpinColorFlavorMatrix lmat = s3g5*0.5;  //*4;
    SpinColorFlavorMatrix rmat = s3g5;

    MFqdpMatrix mf_mat_left(MFstructure::W, MFstructure::V, true, false, lmat);
    MFqdpMatrix mf_mat_right(MFstructure::W, MFstructure::V, true, false, rmat);
  
    MFBasicSource src(MFBasicSource::BoxSource, (double)(GJP.XnodeSites()*GJP.Xnodes()) );
  
    Float p[3] = {0.0,0.0,0.0};
    for(int i=0;i<3;++i) 
      if(GJP.Bc(i) == BND_CND_GPARITY)
	p[i] = 0.5; //units of pi/L
  
    const Float mp[3] = {-p[0],-p[1],-p[2]};

    //Note, FT conventions are  e^{-ipx} for forwards FFT into p-space with momentum +p
    //Hence x and x' have momentum +p and likewise y and y' have momentum -p

#ifdef USE_TRANSCONV_FIELDS
    a2a_prop.gparity_make_fields_transconv(true);
#endif

#ifdef IMPOSE_MOMENTUM
    a2a_prop.set_v_momentum(mp); //x momentum
    a2a_prop.set_wdag_momentum(mp);
#endif

    a2a_prop.fft_vw();

#ifdef USE_TRANSCONV_FIELDS
    a2a_prop.gparity_make_fields_transconv(false);
#endif

    MesonField2 mf_left(a2a_prop,a2a_prop,mf_mat_left,src);

#ifdef USE_TRANSCONV_FIELDS
    a2a_prop.gparity_make_fields_transconv(true);
#endif

#ifdef IMPOSE_MOMENTUM
    a2a_prop.set_v_momentum(p);
    a2a_prop.set_wdag_momentum(p);
#endif

    a2a_prop.fft_vw();

#ifdef USE_TRANSCONV_FIELDS
    a2a_prop.gparity_make_fields_transconv(false);
#endif

    MesonField2 mf_right(a2a_prop,a2a_prop,mf_mat_right,src);
  
    CorrelationFunction result("result",1);
  
    MesonField2::contract(mf_left,mf_right,result);

    if(!UniqueID()) printf("Result with inbuilt transconv code and source and sink timeslice explicit:\n");

    for(int x4=0;x4<GJP.TnodeSites()*GJP.Tnodes();x4++){
      CorrelationFunction result("result",1);
      MesonField2::contract_specify_tsrc(mf_left,mf_right,0,x4,result);
    
      for(int y4=0;y4<GJP.TnodeSites()*GJP.Tnodes();++y4)
	twopoint_results[x4][y4] = result(0,y4);

      if(!UniqueID()){
	for(int y4=0;y4<GJP.TnodeSites()*GJP.Tnodes();++y4){
	  printf("%d %d %.12le %.12le\n",x4,y4,result(0,y4).real(),result(0,y4).imag());
	}
      }

    }
  }


  if(0){
    //Test MesonField2::combine_mf_wv_wv:
    //\sum_{ t2 in 'range' }   [\sum_{\vec x} w_i^dag(t1) v_j(t1; t2)] [\sum_\vec y} w_j^dag(t2) v_k(t2; t3) ]
    //Where the stuff in square brackets are meson fields
    
    //Test it by using the function and taking the trace, which should give us the same as above:
    //\sum_{x,x',y,y'} e^{-ipx}e^{ipy}e^{-ipx'}e^{ipy'} 0.5 * 4 * tr{ [ w^dag(x') g5 \sigma_3 v(x) ] [ w^dag(y) g5 \sigma_3 v(y') ] } 

    //if we take t2=y_4  and  t3=t1=x_4

    //Conclusion: Appears to work although numerical errors are suprisingly large. Maybe this should be expected?

    SpinColorFlavorMatrix s3g5(gamma5,sigma3);
    SpinColorFlavorMatrix lmat = s3g5*0.5*4;
    SpinColorFlavorMatrix rmat = s3g5;

    MFqdpMatrix mf_mat_left(MFstructure::W, MFstructure::V, true, false, lmat);
    MFqdpMatrix mf_mat_right(MFstructure::W, MFstructure::V, true, false, rmat);
  
    MFBasicSource src(MFBasicSource::BoxSource, (double)(GJP.XnodeSites()*GJP.Xnodes()) );
  
    Float p[3] = {0.0,0.0,0.0};
    for(int i=0;i<3;++i) 
      if(GJP.Bc(i) == BND_CND_GPARITY)
	p[i] = 0.5; //units of pi/L
  
    const Float mp[3] = {-p[0],-p[1],-p[2]};

    //Note, FT conventions are  e^{-ipx} for forwards FFT into p-space with momentum +p
    //Hence x and x' have momentum +p and likewise y and y' have momentum -p

#ifdef USE_TRANSCONV_FIELDS
    a2a_prop.gparity_make_fields_transconv(true);
#endif

#ifdef IMPOSE_MOMENTUM
    a2a_prop.set_v_momentum(p);
    a2a_prop.set_wdag_momentum(p);
#endif 

    a2a_prop.fft_vw();
#ifdef USE_TRANSCONV_FIELDS
    a2a_prop.gparity_make_fields_transconv(false);
#endif

    MesonField2 mf_left(a2a_prop,a2a_prop,mf_mat_left,src);

#ifdef USE_TRANSCONV_FIELDS
    a2a_prop.gparity_make_fields_transconv(true);
#endif

#ifdef IMPOSE_MOMENTUM
    a2a_prop.set_v_momentum(mp);
    a2a_prop.set_wdag_momentum(mp);
#endif

    a2a_prop.fft_vw();
#ifdef USE_TRANSCONV_FIELDS
    a2a_prop.gparity_make_fields_transconv(false);
#endif
    MesonField2 mf_right(a2a_prop,a2a_prop,mf_mat_right,src);

    int t_size = GJP.Tnodes()*GJP.TnodeSites();

    std::complex<double> results[t_size][t_size];
    if(!UniqueID()) printf("\n");
    const int NA=-1;

    for(int y_4 = 0; y_4 < t_size; y_4++){
      RangeSpecificT range(y_4);
      MesonField2 prod;
      MesonField2::combine_mf_wv_wv(prod, mf_left, mf_right, range);
      printf("MesonField2::combine_mf_wv_wv y_4=%d combined meson fields of size  %d x %d  and  %d x %d into one of size %d x %d\n",y_4,mf_left.rows(),mf_left.cols(),mf_right.rows(),mf_right.cols(),prod.rows(),prod.cols());

      for(int x_4=0;x_4<t_size;x_4++){
	//Take trace
	std::complex<double> &result = results[x_4][y_4];
	result.real() = 0.0;
	result.imag() = 0.0;
	for(int J=0;J<prod.get_size(MesonField2::Right);J++){
	  int i,t;
	  prod.inv_idx(J,i,t,MesonField2::Right);
	  if(t == x_4 || t == NA)
	    result += *prod.mf_val(i,J,x_4);
	}
      }
    }

    for(int x_4=0;x_4<t_size;x_4++){
      for(int y_4 = 0; y_4 < t_size; y_4++){
      	if(!UniqueID()){
	  printf("%d %d %.12le %.12le, ratio wrt expectance %.12le %.12le\n",x_4,y_4,results[x_4][y_4].real(),results[x_4][y_4].imag(),  results[x_4][y_4].real()/twopoint_results[x_4][y_4].real(),  results[x_4][y_4].imag()/twopoint_results[x_4][y_4].imag() );
	}
      }
    }

  }


  if(0){
    //Test MesonField2::contract_vleft_wright by using it to form the 2 point contraction
    
    //Contract the meson field with the specified vector component of an A2APropbfm at a given site over the mode indices to the left and right, leaving a SpinColorFlavorMatrix
    //In terms of mode indices that include the time-slice dilution in the mode index, this is
    //  result_{a,b} = \sum_{I,J} v_{aI}(x) M_{IJ}(t) w^dag_{bJ}(z)
    //where M is this MesonField, assumed to be of the form w^dag \otimes v
    //In terms of explicit source time-slice indices for v (to right of semicolon) [implicit spin-color sum in square brackets]
    //  result_{a,b} = \sum_{i,j} v_{ai}(\vec x,x_4 ; t) [\sum_{\vec y} w_i^dag(\vec y , t) v_j(\vec y,t ; z_4)] w^dag_{bj}(\vec z,z_4)
    //Note, x,z are local coordinates,  t is a *global* coordinate

    //Want: \sum_{x,x',y,y'} e^{-ipx}e^{ipy}e^{-ipx'}e^{ipy'} 0.5 * 4 * tr{ [ w^dag(x') g5 \sigma_3 v(x) ] [ w^dag(y) g5 \sigma_3 v(y') ] } 

    SpinColorFlavorMatrix s3g5(gamma5,sigma3);
    SpinColorFlavorMatrix rmat = s3g5;

    SpinColorFlavorMatrix s1g5(gamma5,sigma1);
    static std::complex<Float> _i(0.0,1.0);

#ifdef USE_TRANSCONV_FIELDS
    SpinColorFlavorMatrix lmat = s3g5 - s1g5*_i;
#else
    SpinColorFlavorMatrix lmat = s3g5*2.0;
#endif

    MFqdpMatrix mf_mat_right(MFstructure::W, MFstructure::V, true, false, rmat);
  
    MFBasicSource src(MFBasicSource::BoxSource, (double)(GJP.XnodeSites()*GJP.Xnodes()) );
  
    Float p[3] = {0.0,0.0,0.0};
    for(int i=0;i<3;++i) 
      if(GJP.Bc(i) == BND_CND_GPARITY)
	p[i] = 0.5; //units of pi/L
  
    const Float mp[3] = {-p[0],-p[1],-p[2]};
#ifdef USE_TRANSCONV_FIELDS
    a2a_prop.gparity_make_fields_transconv(true);
#endif

#ifdef IMPOSE_MOMENTUM
    a2a_prop.set_v_momentum(mp);
    a2a_prop.set_wdag_momentum(mp);
#endif
    a2a_prop.fft_vw();

#ifdef USE_TRANSCONV_FIELDS
    a2a_prop.gparity_make_fields_transconv(false);
#endif

    MesonField2 mf_right(a2a_prop,a2a_prop,mf_mat_right,src);

    //Loop over source timeslice and coords x and x' on that timeslice
    int loc3vol = GJP.VolNodeSites()/GJP.TnodeSites();
    static const Float _mpi = -3.141592654;

    if(!UniqueID()) printf("Doing 2-pt function for testing MesonField2::contract_vleft_wright\n");

    //    for(int x_4=0;x_4<t_size;x_4++){
    for(int x_4=0;x_4<1;x_4++){
      CorrelationFunction result("result",1,CorrelationFunction::THREADED); //result as a function of y_4
      
      //Only the nodes with x_4 on-node contribute
      int x_4_loc = x_4 - GJP.TnodeCoor()*GJP.TnodeSites();
      if( x_4_loc >= 0 && x_4_loc < GJP.TnodeSites() ){
	
	for(int y_4=0;y_4<t_size;y_4++){
#pragma omp parallel for
	  for(int x3d_loc=0 ; x3d_loc<loc3vol ; x3d_loc++){
	    for(int xpr3d_loc=0 ; xpr3d_loc<loc3vol ; xpr3d_loc++){
	      int x4d_loc = x3d_loc + loc3vol*x_4_loc;
	      int xpr4d_loc = xpr3d_loc + loc3vol*x_4_loc;

	      SpinColorFlavorMatrix scf;
	      MesonField2::contract_vleft_wright(scf,a2a_prop,x4d_loc,a2a_prop,xpr4d_loc, mf_right, y_4);	      

#ifdef DO_GAUGE_FIX
	      //Need gauge fixing matrices too
	      const cps::Matrix* gfmatx_f[2] = { lattice->FixGaugeMatrix(x4d_loc,0), lattice->FixGaugeMatrix(x4d_loc,1) };
	      scf(0,0).LeftTimesEqual(*gfmatx_f[0]);
	      scf(0,1).LeftTimesEqual(*gfmatx_f[0]);
	      scf(1,0).LeftTimesEqual(*gfmatx_f[1]);
	      scf(1,1).LeftTimesEqual(*gfmatx_f[1]);
#endif

	      //Multiply by spin-flavor structure
	      scf.LeftTimesEqual(lmat);

#ifdef DO_GAUGE_FIX
	      //For wdag we need to left-multiply by V^dag where V is the gauge fixing matrix
	      cps::Matrix gfmatxpr_f[2];
	      gfmatxpr_f[0].Dagger( *(lattice->FixGaugeMatrix(xpr4d_loc,0)) );
	      gfmatxpr_f[1].Dagger( *(lattice->FixGaugeMatrix(xpr4d_loc,1)) );

	      scf(0,0).LeftTimesEqual(gfmatxpr_f[0]);
	      scf(0,1).LeftTimesEqual(gfmatxpr_f[0]);
	      scf(1,0).LeftTimesEqual(gfmatxpr_f[1]);
	      scf(1,1).LeftTimesEqual(gfmatxpr_f[1]);
#endif

	      //Now the phases
	      cps::Complex phases[2] = { cps::Complex(1.0,0.0), cps::Complex(1.0,0.0) };

#ifdef IMPOSE_MOMENTUM
	      int pos_loc[2] = { x3d_loc, xpr3d_loc };
	      for(int pp=0;pp<2;pp++){
		Float pdotx = 0.0;
		int rem(pos_loc[pp]);
		for(int d=0;d<3;d++){ 
		  int szd = GJP.NodeSites(d);
		  int pos_d_loc = rem % szd; rem /= szd;
		  int pos_d = pos_d_loc + GJP.NodeCoor(d)*GJP.NodeSites(d);
		  pdotx += p[d] * pos_d * _mpi / szd;
		}
		phases[pp] = cps::Complex( cos(pdotx), sin(pdotx) );
	      }
#endif

	      //Take the trace
	      result(omp_get_thread_num(),0,y_4) += 0.5*2* scf.Trace() * phases[0]*phases[1];
	    }
	  }
	}
      }
      cps::sync(); //some nodes will be waiting
      result.sumLattice();

      if(!UniqueID()){
	for(int y4=0;y4<t_size;++y4){
	  printf("%d %d %.12le %.12le\n",x_4,y4,result(0,y4).real(),result(0,y4).imag());
	}
      }
      
    }//end of source time loop

  }

  if(0){
    //TESTING
    ContractionTypeKtoPiPi ktopipi_args;
    setup_ktopipi_args(ktopipi_args,"a2aprop","a2aprop_s",'+','-');

    Gparity_KtoPiPi gpcon;
    gpcon.setup(ktopipi_args,*lattice);

    Testing::test_prod1(gpcon,ktopipi_args,*lattice);
    Testing::test_conLLLH(gpcon,ktopipi_args,*lattice);
  }
  if(0){
    //TESTING 2
    ContractionTypeKtoPiPi ktopipi_args;
    setup_ktopipi_args(ktopipi_args,"a2aprop","a2aprop_s",'+','-');

    Gparity_KtoPiPi gpcon;
    gpcon.setup(ktopipi_args,*lattice);
    Testing::type1_propHg5conj_test(gpcon,ktopipi_args,*lattice);
  }
  if(0){
    //TESTING 3
    ContractionTypeKtoPiPi ktopipi_args;
    setup_ktopipi_args(ktopipi_args,"a2aprop","a2aprop_s",'+','-');

    Gparity_KtoPiPi gpcon;
    gpcon.setup(ktopipi_args,*lattice);
    Testing::type4_test(gpcon,ktopipi_args,*lattice);
  }

  if(0){
    //Do the K->(pipi)_I=0 contractions
    //setup_ktopipi_args(ContractionTypeKtoPiPi &into, const char *prop_L, const char* prop_H){
    ContractionTypeKtoPiPi ktopipi_args;
    setup_ktopipi_args(ktopipi_args,"a2aprop","a2aprop_s",'+','-');

    Gparity_KtoPiPi gpcon;
    gpcon.setup(ktopipi_args,*lattice);
  
    int tsize = GJP.Tnodes()*GJP.TnodeSites();

    //Also write out in Daiqian's format
    FILE *fp = fopen("type1.ck","w");

    //Daiqian's output shows different results for each kaon timeslice
    for(int tk = 0; tk < tsize; ++tk){
      std::vector<int> tkv(1,tk);
      std::vector<CorrelationFunction> result;

      if(!UniqueID()) printf("Calculating type1 diagrams with tk=%d\n",tk);
      //gpcon.type1(2,1,result,&tkv);
      gpcon.type1_propHg5conj(2,1,result,&tkv);

      for(int c=0;c<result.size();c++){
	if(!UniqueID()) printf("Contraction %s\n",result[c].getLabel().c_str());
	for(int top=0;top<tsize;++top){
	  if(!UniqueID()) printf("%d %d %.12le %.12le\n",tk,top,result[c](0,top) );
	}
      }

      //Daiqian prints out
      //<t_K> <(t_op-t_K)%T> <C1(0,V,0,A)> <C1(0,A,0,V)> <C1(0,V,1,A)> <C1(0,A,1,V)> <C1(1,V,0,A)> <C1(1,A,0,V)> <C1(1,V,1,A)> <C1(1,A,1,V)> <C2.............
      //where each C contains both the real and imaginary parts
      //And prints 0 when  t_op - t_K outside the range  0 < t_op - t_K < delta    where delta is the pi-K separation (2 here)
      //const static std::string GammaNames[4] = { "M_{0,V}","M_{0,A}", "M_{1,V}", "M_{1,A}" };
      const int _0V(0);
      const int _0A(1);
      const int _1V(2);
      const int _1A(3);
      const int delta(2);
      const int Lt( GJP.Tnodes()*GJP.TnodeSites() );

      int ncon =  result.size()/4/4 ;

      //We want the second time index to be tt= (top-tk+Lt) % Lt in order
      std::vector< std::vector< std::vector<cps::Complex> > > ordered(tsize);  //[tt][c][gammacomb]
      for(int top=0;top<tsize;++top){
	int tt = (top-tk+Lt )%Lt;
	ordered[tt].resize(ncon);

	printf("tk=%d top=%d -> tt=%d\n",tk,top,tt);

	for(int c=0;c<ncon;c++){
	  std::vector<cps::Complex> &v = ordered[tt][c];
	  v.resize( 8,cps::Complex(0.0,0.0) );

	  //if( (top-tk+Lt)%Lt > 0 && (top-tk+Lt)%Lt < delta ){
	    v[0] = result[gpcon.result_map(c,_0V,_0A)](0,top);
	    v[1] = result[gpcon.result_map(c,_0A,_0V)](0,top);

	    v[2] = result[gpcon.result_map(c,_0V,_1A)](0,top);
	    v[3] = result[gpcon.result_map(c,_0A,_1V)](0,top);

	    v[4] = result[gpcon.result_map(c,_1V,_0A)](0,top);
	    v[5] = result[gpcon.result_map(c,_1A,_0V)](0,top);

	    v[6] = result[gpcon.result_map(c,_1V,_1A)](0,top);
	    v[7] = result[gpcon.result_map(c,_1A,_1V)](0,top);
	    //}
	  
	}
      }

      for(int tt=0;tt<tsize;++tt){
	fprintf(fp,"%d %d",tk,tt);
	for(int c=0;c<ncon;c++){
	  std::vector<cps::Complex> &v = ordered[tt][c];
	  for(int vv=0;vv<8;vv++) fprintf(fp," %.6le %.6le", std::real(v[vv]), std::imag(v[vv]) );	  
	}
	fprintf(fp,"\n");
      }


    }
    fclose(fp);
  }

  //Original test with Daiqian
  if(0){
    //For the other types we avoid the ambiguity in momentum directions / timeslices for the two pions by summing over all combinations:
    //pi1(t, +p), pi2(t+delta, -p)  +   pi1(t+delta, +p), pi2(t, -p)  +   pi1(t, -p), pi2(t+delta, +p)  +   pi1(t+delta, -p), pi2(t,+p) 
    ContractionTypeKtoPiPi ktopipi_args_pm,  ktopipi_args_mp;
    setup_ktopipi_args(ktopipi_args_pm,"a2aprop","a2aprop_s",'+','-');
    setup_ktopipi_args(ktopipi_args_mp,"a2aprop","a2aprop_s",'-','+');

    Gparity_KtoPiPi gpcon_pm;
    gpcon_pm.setup(ktopipi_args_pm,*lattice);
  
    Gparity_KtoPiPi gpcon_mp;
    gpcon_mp.setup(ktopipi_args_mp,*lattice);

    int tsize = GJP.Tnodes()*GJP.TnodeSites();

    int n_types = 6;
    bool do_type[] = { false, false, false, false, false, false };
    std::string names[] = { "type1.ck", "type2.ck", "type3.ck", "type4.ck", "psvertex_type3.ck", "psvertex_type4.ck" };
    do_type[4] = true; //pseudoscalar vertex of type4 form

    //Also write out in Daiqian's format
    FILE *fp[n_types]; 
    for(int i=0;i<n_types;i++){
      if(do_type[i]){
	fp[i] = fopen(names[i].c_str(),"w");
      }else fp[i] = NULL;
    }

    //Daiqian's output shows different results for each kaon timeslice
    for(int tk = 0; tk < tsize; ++tk){
      std::vector<int> tkv(1,tk);
      std::vector<CorrelationFunction> result_pm, result_mp;

      for(int type=0; type < n_types; ++type){
	if(!do_type[type]) continue;

	if(type == 0){
	  gpcon_pm.type1(2,1,result_pm,&tkv);
	  gpcon_mp.type1(2,1,result_mp,&tkv);
	}else if(type == 1){
	  gpcon_pm.type2(2,1,result_pm,&tkv);
	  gpcon_mp.type2(2,1,result_mp,&tkv);
	}else if(type == 2){
	  gpcon_pm.type3(2,1,result_pm,&tkv);
	  gpcon_mp.type3(2,1,result_mp,&tkv);
	}else if(type == 3){
	  gpcon_pm.type4(2,1,result_pm,&tkv);
	  gpcon_mp.type4(2,1,result_mp,&tkv);
	}else if(type == 4){
	  gpcon_pm.psvertex_type3(2,1,result_pm,&tkv);
	  gpcon_mp.psvertex_type3(2,1,result_mp,&tkv);
	}else if(type == 5){
	  gpcon_pm.psvertex_type4(2,1,result_pm,&tkv);
	  gpcon_mp.psvertex_type4(2,1,result_mp,&tkv);
	}else{
	  printf("Invalid contraction type\n");
	  exit(-1);
	}

	//Sum pm and mp
	std::vector<CorrelationFunction> sum(result_pm.size());
	for(int ii=0;ii<sum.size();ii++){
	  sum[ii].setNcontractions(1);
	  for(int t=0;t<tsize;t++)
	    sum[ii](0,t) = result_pm[ii](0,t) + result_mp[ii](0,t);
	}

	//Daiqian prints out
	//<t_K> <(t_op-t_K)%T> <C1(0,V,0,A)> <C1(0,A,0,V)> <C1(0,V,1,A)> <C1(0,A,1,V)> <C1(1,V,0,A)> <C1(1,A,0,V)> <C1(1,V,1,A)> <C1(1,A,1,V)> <C2.............
	//where each C contains both the real and imaginary parts
	//And prints 0 when  t_op - t_K outside the range  0 < t_op - t_K < delta    where delta is the pi-K separation (2 here)
	//const static std::string GammaNames[4] = { "M_{0,V}","M_{0,A}", "M_{1,V}", "M_{1,A}" };
	std::vector< std::vector< std::vector<cps::Complex> > > ordered;
	daiqian_order_result(ordered,sum,type+1,tk);

	for(int tt=0;tt<tsize;++tt){
	  fprintf(fp[type],"%d %d",tk,tt);
	  for(int c=0;c<ncon;c++){
	    std::vector<cps::Complex> &v = ordered[tt][c];
	    for(int vv=0;vv<v.size();vv++) fprintf(fp[type]," %.16le %.16le", std::real(v[vv]), std::imag(v[vv]) );	  
	  }
	  fprintf(fp[type],"\n");
	}
	
	
      }


    }
    for(int f=0;f<n_types;f++)
      if(fp[f]!=NULL) fclose(fp[f]);

  }

  if(1){
    //08/12/14  Test against Daiqian's results where each pion comprises two momentum combinations
    //Have to hand-specify the momenta
    int ngp = 0;
    for(int i=0;i<3;i++) if(GJP.Bc(i) == BND_CND_GPARITY) ++ngp;

    //We want to consider pions moving along the G-parity axis and perpendicular to it
    std::vector< std::pair< std::vector<int>, std::vector<int> > > mom_comb_diag;
    std::vector< std::pair< std::vector<int>, std::vector<int> > > mom_comb_perp; 
		
    if(ngp==3){
      //Two combinations for pion moving in +++ direction
      mom_comb_diag.resize(2);
      mom_comb_diag[0].first = three_vec(1,1,1);
      mom_comb_diag[0].second = three_vec(1,1,1);

      mom_comb_diag[1].first = three_vec(-1,-1,-1);
      mom_comb_diag[1].second = three_vec(3,3,3);
      
      //Two combinations for pion moving in -++
      mom_comb_perp.resize(2);
      mom_comb_perp[0].first = three_vec(1,1,1);
      mom_comb_perp[0].second = three_vec(-3,1,1);

      mom_comb_perp[1].first = three_vec(-1,-1,-1);
      mom_comb_perp[1].second = three_vec(-1,3,3);
    }else{
      printf("Momentum combinations for ngp = %d not specified\n",ngp); exit(-1);
    }
    
    int n_types = 7;
    bool do_type[] = {false,     false, false, false, false, false, false };  //first is dummy
    std::string names[] = {"",      "type1", "type2", "type3", "type4", "psvertex_type3", "psvertex_type4" };
    int ncon[] = {-1,   6,6,10,10,1,1};
    do_type[1] = true;

    //Because the problem is symmetric under rotations around the diagonal axis, we need only consider the diagonal axis and one off-diagonal
    //First consider   pi1(+++) pi2(---). There are 4 combinations (2x2), which we sum. Then do the parity-flipped version, then similar for the perpendicular momenta
    std::vector< std::pair< std::vector<int>, std::vector<int> > >* dir_mcomb[2] = { &mom_comb_diag, &mom_comb_perp };
    const char* pi1_mom_names[] = { "111", "_1_1_1", "_111", "1_1_1" }; //Daiqian's notation

    for(int dir = 0; dir < 2; dir++){ //dir=0 -> diag,  dir=1 -> perp
      const std::vector< std::pair< std::vector<int>, std::vector<int> > > &mom_comb = *dir_mcomb[dir];      
      for(int parity = 0; parity < 2; parity++){

	const char* pi1_mom_name = pi1_mom_names[parity + 2*dir]; //for filenames

	for(int type = 1; type < n_types; type++){
	    if(!do_type[type]) continue;

	    std::ostringstream fn; fn << "type" << type << "_deltat_2_sep_1_mom" << pi1_mom_name << ".ck";
	    FILE *fp = fopen(fn.str().c_str(),"w");

	    for(int tk = 0; tk < tsize; ++tk){

	      //Sum over the 4 combinations of quark momenta
	      std::vector<CorrelationFunction> sum(4*4*ncon[type]);   //(4*4 operator pairs) * (number of contractions of type)
	      for(int ii=0;ii<sum.size();ii++) sum[ii].setNcontractions(1); //zeroes

	      for(int pcomb1=0;pcomb1<2;pcomb1++){
		for(int pcomb2=0;pcomb2<2;pcomb2++){
		  int idx = pcomb2 + 2*pcomb[1];

		  const std::pair< std::vector<int>, std::vector<int> > & p_pi1 = mom_comb_diag[pcomb1];
		  const std::pair< std::vector<int>, std::vector<int> > & minus_p_pi2 = mom_comb_diag[pcomb2]; //is in +++ dir, we need --- so apply minus below

		  std::vector<int> p_pi1_q1 =  ( parity==0 ? p_pi1.first : minus_vec( p_pi1.first ) );
		  std::vector<int> p_pi1_q2 =  ( parity==0 ? p_pi1.second : minus_vec( p_pi1.second ) );
		  std::vector<int> p_pi2_q1 =  ( parity==0 ? minus_vec( p_pi2.first ) : p_pi1.first );
		  std::vector<int> p_pi2_q2 =  ( parity==0 ? minus_vec( p_pi2.second ) : p_pi1.second );

		  Gparity_KtoPiPi gpcon;
		  ContractionTypeKtoPiPi gpcon_args;

		  setup_ktopipi_multimom_args(gpcon_args,"a2aprop","a2aprop_s", 
					      p_pi1_q1, p_pi1_q1,
					      p_pi2_q1, p_pi2_q2);
		
		  gpcon.setup(gpcon_args,*lattice);
		  
		  std::vector<CorrelationFunction> tmp;
		  compute_type(tmp,type,tk,2,1,gpcon);
		
		  for(int ii=0;ii<sum.size();ii++)
		    for(int t=0;t<GJP.Tnodes()*GJP.TnodeSites();t++)
		      sum[ii](0,t) += tmp[ii](0,t);
	    
		}
	      }

	      //Order as in Daiqian's files
	      std::vector< std::vector< std::vector<cps::Complex> > > ordered;//[tt][c][gammacomb]
	      daiqian_order_result(ordered,sum,type,tk);

	      //Write this tk to file
	      for(int tt=0;tt<tsize;++tt){
		fprintf(fp,"%d %d",tk,tt);
		for(int c=0;c<ncon;c++){
		  std::vector<cps::Complex> &v = ordered[tt][c];
		  for(int vv=0;vv<v.size();vv++) fprintf(fp[type]," %.16le %.16le", std::real(v[vv]), std::imag(v[vv]) );	  
		}
		fprintf(fp[type],"\n");
	      }

	    }//end of tk loop
	    
	}//end of type
      }//parity loop
    }//dir loop


  }//end of test













  }

  if(0){
    //pipi scattering
    ContractionTypeKtoPiPi ktopipi_args;
    setup_ktopipi_args(ktopipi_args,"a2aprop","a2aprop_s",'+','-');

    Gparity_KtoPiPi gpcon;
    gpcon.setup(ktopipi_args,*lattice);
  
    int tsize = GJP.Tnodes()*GJP.TnodeSites();

    //Also write out in Daiqian's format
    std::string filenames[4] = {"pipi_C_+.ck", "pipi_D_+.ck", "pipi_R_+.ck", "pipi_V_+.ck"}; //+ indicates that p_pi1_src = p_pi1_snk

    FILE *fp[4];
    for(int i=0;i<4;i++)
      fp[i] = fopen(filenames[i].c_str(),"w");
    
    std::vector< std::vector< std::vector< std::complex<double> > > > reorg(t_size);//T * T * 4

    for(int t=0;t<t_size;t++){
      reorg[t].resize(t_size);
      for(int u=0;u<t_size;u++)
	reorg[t][u].resize(4);
    }

    //My output is x4, (y4-x4+T)%T   where x4 is the pi1 source timeslice and y4 the pi1 sink timeslice
    //For pipi sep 1 (as here),  r4=(x4+1)%T and s4=(y4+1)%T are the source and sink timeslices of pi2
    //Daiqian writes out r4 (y4-r4+T)%T so I reorganise my output to match
    //matches traj_0_FigureR_sep1_mom_1_1_1_mom111

    for(int x4=0;x4<tsize;x4++){
      std::vector<int> x4_rng(1,x4);
      std::vector<CorrelationFunction> result;
      gpcon.pipi(1,'+',result,&x4_rng);

      for(int tsep = 0; tsep<tsize; tsep++){
	int r4=(x4+1)%t_size;
	int y4=(tsep+x4)%t_size;
	int s4=(y4+1)%t_size;

	int dq_sep = (y4-r4+t_size)%t_size;

	for(int c=0;c<4;c++) reorg[x4][tsep][c] = result[c](0,tsep);
	//for(int c=0;c<4;c++) reorg[r4][dq_sep][c] = result[c](0,tsep);
      }
    }

    for(int tpi1_src=0;tpi1_src<tsize;tpi1_src++){
      //Output   <r4> <(y4-r4)%T> <real> <imag>
      for(int c=0;c<4;c++){
	for(int tsep = 0; tsep<tsize; tsep++){	
	  fprintf(fp[c],"%d %d %.16le %.16le\n",tpi1_src,tsep, reorg[tpi1_src][tsep][c].real(), reorg[tpi1_src][tsep][c].imag());
	}
      }
    }
    
    for(int i=0;i<4;i++)
      fclose(fp[i]);

  }


  if(0){
    //pipi scattering
    //In this version   pi1(-p)pi2(+p) -> pi1(+p)pi2(-p)

    ContractionTypeKtoPiPi ktopipi_args;
    setup_ktopipi_args(ktopipi_args,"a2aprop","a2aprop_s",'-','+');

    Gparity_KtoPiPi gpcon;
    gpcon.setup(ktopipi_args,*lattice);
  
    int tsize = GJP.Tnodes()*GJP.TnodeSites();

    //Also write out in Daiqian's format
    std::string filenames[4] = {"pipi_C_-+.ck", "pipi_D_-+.ck", "pipi_R_-+.ck", "pipi_V_-+.ck"}; 

    FILE *fp[4];
    for(int i=0;i<4;i++)
      fp[i] = fopen(filenames[i].c_str(),"w");
    
    std::vector< std::vector< std::vector< std::complex<double> > > > reorg(t_size);//T * T * 4

    for(int t=0;t<t_size;t++){
      reorg[t].resize(t_size);
      for(int u=0;u<t_size;u++)
	reorg[t][u].resize(4);
    }

    for(int x4=0;x4<tsize;x4++){
      std::vector<int> x4_rng(1,x4);
      std::vector<CorrelationFunction> result;
      gpcon.pipi(1,'-',result,&x4_rng); //pi1 at sink has opposite momentum to that at source

      for(int tsep = 0; tsep<tsize; tsep++){
	// int r4=(x4+1)%t_size;
	// int y4=(tsep+x4)%t_size;
	// int s4=(y4+1)%t_size;

	// int dq_sep = (y4-r4+t_size)%t_size;

	for(int c=0;c<4;c++) reorg[x4][tsep][c] = result[c](0,tsep);
	//for(int c=0;c<4;c++) reorg[r4][dq_sep][c] = result[c](0,tsep);
      }
    }

    for(int tpi1_src=0;tpi1_src<tsize;tpi1_src++){
      //Output   <r4> <(y4-r4)%T> <real> <imag>
      for(int c=0;c<4;c++){
	for(int tsep = 0; tsep<tsize; tsep++){	
	  fprintf(fp[c],"%d %d %.16le %.16le\n",tpi1_src,tsep, reorg[tpi1_src][tsep][c].real(), reorg[tpi1_src][tsep][c].imag());
	}
      }
    }
    
    for(int i=0;i<4;i++)
      fclose(fp[i]);

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
