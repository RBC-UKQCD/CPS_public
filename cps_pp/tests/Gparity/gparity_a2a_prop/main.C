//Test G-parity modified version of Daiqian's a2a code. This presupposes that the Lanczos works, so test this first (gparity_lanczos)
//NOTE: You will need to link against libfftw3 and libfftw3_threads

//Enable flavor dilution with -dilute_flavor command line option


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

#include<omp.h>

#ifdef HAVE_BFM
#include <chroma.h>
#endif

using namespace std;
USING_NAMESPACE_CPS

void setup_double_latt(Lattice &double_latt, cps::Matrix* orig_gfield, bool gparity_X, bool gparity_Y){
  //orig latt ( U_0 U_1 ) ( U_2 U_3 ) ( U_4 U_5 ) ( U_6 U_7 )
  //double tatt ( U_0 U_1 U_2 U_3 ) ( U_4 U_5 U_6 U_7 ) ( U_0* U_1* U_2* U_3* ) ( U_4* U_5* U_6* U_7* )

  cps::Matrix *dbl_gfield = double_latt.GaugeField();

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
void setup_double_matrixfield(cps::Matrix* double_mat, cps::Matrix* orig_mat, int nmat_per_site, bool gparity_X, bool gparity_Y){
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
  
void GaugeTransformU(cps::Matrix *gtrans, Lattice &lat);

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

void lanczos_arg(LancArg &into, const bool &precon){
  into.mass = 0.01;
  into.stop_rsd = 1e-06;
  into.qr_rsd = 1e-14; ///convergence of intermediate QR solves, defaults to 1e-14
  into.EigenOper = DDAGD;
  into.precon = precon; //also try this with true
  into.N_get = 10;///Want K converged vectors
  into.N_use = 11;///Dimension M Krylov space
  into.N_true_get = 10;//Actually number of eigen vectors you will get
  into.ch_ord = 3;///Order of Chebyshev polynomial
  into.ch_alpha = 5.5;///Spectral radius
  into.ch_beta = 0.5;///Spectral offset (ie. find eigenvalues of magnitude less than this)
  into.ch_sh = false;///Shifting or not
  into.ch_mu = 0;///Shift the peak
  into.lock = false;///Use locking transofrmation or not
  into.maxits =10000;///maxiterations
  into.fname = "Lanczos";
}

void a2a_arg(A2AArg &into, const int &flavor_dilution = 0, const RandomType &rand_type = UONE, const int &nl= 8){
  into.nl = nl;
  into.nhits = 1;
  into.rand_type = rand_type;
  into.src_width = 1;
  into.dilute_flavor = flavor_dilution;
  into.do_gauge_fix = true;
}

void create_eig(GwilsonFdwf* lattice, Lanczos_5d<double>* &eig, const bool &precon){
  //Run in 2f environment. Version for 1f is directly converted from this
  bfm_evo<double>* dwf = new bfm_evo<double>();
  bfmarg dwfa;
  setup_bfmargs(dwfa);

  dwf->init(dwfa);
  
  lattice->BondCond(); //Don't forget to apply the boundary conditions!
  Float* gauge = (Float*) lattice->GaugeField();
  dwf->cps_importGauge(gauge); 

  //Setup and run the lanczos algorithm
  LancArg lanc_arg; lanczos_arg(lanc_arg,precon);
  eig = new Lanczos_5d<double>(*dwf,lanc_arg);
  eig->Run();
  lattice->BondCond(); //Don't forget to un-apply the boundary conditions!
}
  
//Note, only the eigenvectors and eigenvalues of the Lanczos_5d instance are used, hence for a cleaner test I generate these on the 2f lattice
//and simply convert them to the 1f lattice
//The extraction of the eigenvectors to CPS format (needed for conversion) should be done BEFORE moving to the 1f environment
void eig_convert_cps(Lanczos_5d<double> &eig, Float** &eigv_2f_cps, const bool &precon){
  //EXECUTED IN 2F ENVIRONMENT
  long f_size = (long)2*24 * GJP.VolNodeSites() * GJP.SnodeSites();
  multi1d<bfm_fermion> &eigenvecs = eig.bq;
  
  eigv_2f_cps = (Float**)pmalloc(eigenvecs.size() * sizeof(Float*) );
  for(int i=0;i<eigenvecs.size();i++){
    eigv_2f_cps[i] = (Float *)pmalloc(sizeof(Float) * f_size); 
    if(precon){ eigenvecs[i][0] = eig.dop.allocCompactFermion();  eig.dop.set_zero(eigenvecs[i][0]); } //for precon ferms this is not allocated    
    eig.dop.cps_impexFermion(eigv_2f_cps[i],eigenvecs[i],0);
  }
}

void eig_convert_2f_1f(Lanczos_5d<double> &eig, Float** eigv_2f_cps, const std::vector<double> &eval_2f, const bool &precon, const bool &gparity_X, const bool &gparity_Y){
  //THIS SHOULD BE EXECUTED WITHIN THE 1F ENVIRONMENT
  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();
  multi1d<bfm_fermion> &eigenvecs = eig.bq;
  //temp storage
  Float *v_1f_cps = (Float *)pmalloc(sizeof(Float) * f_size); 

  for(int i=0;i<eigenvecs.size();i++){
    setup_double_5d_vector((Vector*)v_1f_cps,(Vector*)eigv_2f_cps[i], gparity_X,gparity_Y);
    if(precon){ eigenvecs[i][0] = eig.dop.allocCompactFermion(); } //for precon ferms this is not allocated 
    eig.dop.cps_impexFermion(v_1f_cps,eigenvecs[i],1); //import
    if(precon){ eig.dop.freeFermion(eigenvecs[i][0]); } //delete the even checkerboard that we allocated above
    pfree(eigv_2f_cps[i]);
  }
  eig.evals = eval_2f;

  pfree(v_1f_cps);
  pfree(eigv_2f_cps);
}

//Create eig but do not run, instead copy over evals and converted evecs from 2f
void create_eig_1f(GwilsonFdwf* lattice, Lanczos_5d<double>* &eig_1f, const Lanczos_5d<double> &eig_2f, Float** eigv_2f_cps, const bool &precon, const bool &gparity_X, const bool &gparity_Y){
  bfm_evo<double>* dwf = new bfm_evo<double>();
  bfmarg dwfa;
  setup_bfmargs(dwfa);

  dwf->init(dwfa);
  
  lattice->BondCond(); //Don't forget to apply the boundary conditions!
  Float* gauge = (Float*) lattice->GaugeField();
  dwf->cps_importGauge(gauge); 
  
  LancArg lanc_arg; lanczos_arg(lanc_arg,precon);
  eig_1f = new Lanczos_5d<double>(*dwf,lanc_arg);  
  eig_1f->do_init(); //normally called in Run
  eig_convert_2f_1f(*eig_1f, eigv_2f_cps, eig_2f.evals, precon, gparity_X, gparity_Y);
  lattice->BondCond(); //Don't forget to un-apply the boundary conditions!
}


CPS_START_NAMESPACE

class A2APropbfmTesting{
public:
  //Lanczos should be precomputed. prop pointer will be assigned to an a2a propagator
  static void a2a_prop_gen(A2APropbfm* &prop, Lattice* lattice, Lanczos_5d<double> &eig, const int &dilute_flavor = 0, const RandomType &rand_type = UONE, const int &nl= 8){
    A2AArg arg;
    a2a_arg(arg,dilute_flavor,rand_type,nl);
  
    bfm_evo<double> &dwf = eig.dop;

    CommonArg carg;
    prop = new A2APropbfm(*lattice,arg,carg,&eig);

    prop->allocate_vw();
    QDPIO::cout << "Computing A2a low modes component\n";
    prop->compute_vw_low(dwf);
    QDPIO::cout << "Computing A2a high modes component\n";
    prop->compute_vw_high(dwf);
    QDPIO::cout << "Completed A2A prop\n";
  }

  static void convert_2f_A2A_prop_to_1f(A2APropbfm &prop_2f, bool gparity_X, bool gparity_Y){
    //Called in 1f context!

    int four_vol = GJP.VolNodeSites(); //size includes factor of 2 for second G-parity flavour as we are in the 1f context
    //Vector **v; // v  size   [nvec][24*four_vol *2(gp)]
    for(int vec = 0; vec < prop_2f.get_nvec(); ++vec){
      int vsz = 24*four_vol*sizeof(Float);
      Float* new_v = (Float*)pmalloc(vsz);
      SingleToDoubleField doubler(gparity_X, gparity_Y, 24, (Float*)prop_2f.get_v(vec), new_v);
      doubler.Run();
      memcpy( (void*)prop_2f.get_v(vec), (void*)new_v, vsz);
      pfree(new_v);
    }
  
    //Vector **wl; // low modes w [A2AArg.nl][24*four_vol *2(gp)]

    for(int l = 0; l < prop_2f.get_nl(); ++l){
      int wlsz = 24*four_vol*sizeof(Float);
      Float* new_wl = (Float*)pmalloc(wlsz);
      SingleToDoubleField doubler(gparity_X, gparity_Y, 24, (Float*)prop_2f.get_wl(l), new_wl);
      doubler.Run();
      memcpy( (void*)prop_2f.get_wl(l), (void*)new_wl, wlsz);
      pfree(new_wl);
    }

    //Vector *wh; // high modes w, compressed  [2*four_vol*nh_site*2]   where nh_site = A2AArg.nhits

    //G-parity wh mapping  (re/im=[0,1], site=[0..four_vol-1],flav=[0,1],hit=[0..nh_site]) ->    re/im + 2*( site + flav*four_vol + hit*2*four_vol )
    //Standard or G-parity and flavor dilution (re/im=[0,1], site=[0..four_vol-1],hit=[0..nh_site]) ->    re/im + 2*( site + hit*four_vol )
    {
      if(!prop_2f.a2a.dilute_flavor){
	int hit_stride = 2*four_vol; //complex number * sites   (in double lattice context so four_vol is actually twice the original four_vol, and covers both flavours)
	int wh_bytes = hit_stride*prop_2f.get_nhits()*sizeof(Float);
	Float* new_wh = (Float*)pmalloc(wh_bytes);
	
	for(int hit = 0; hit < prop_2f.get_nhits(); hit++){
	  SingleToDoubleField doubler(gparity_X, gparity_Y, 2, (Float*)prop_2f.get_wh() + hit*hit_stride, new_wh + hit*hit_stride);
	  doubler.Run();
	}
	memcpy( (void*)prop_2f.get_wh(), (void*)new_wh, wh_bytes);
	pfree(new_wh);
      }else{
	//wh is the same for both flavours and is not stored in a double field. To move to a double lattice version it is simplest just
	//to copy wh onto a double-size vector and use the usual doubling routines

	int hit_stride_1f = 2*four_vol;
	int hit_stride_2f = four_vol; //complex number * number of sites for 1 flavour (half number of sites of double lattice)

	int wh_bytes = hit_stride_1f*prop_2f.get_nhits()*sizeof(Float);
	Float* new_wh = (Float*)pmalloc(wh_bytes);	
	Float* double_wh = (Float*)pmalloc(wh_bytes);

	for(int hit = 0; hit < prop_2f.get_nhits(); hit++){
	  Float* from = (Float*)prop_2f.get_wh() + hit*hit_stride_2f;
	  memcpy( (void*)double_wh, (void*)from, hit_stride_2f*sizeof(Float));
	  memcpy( (void*)(double_wh+hit_stride_2f), (void*)from, hit_stride_2f*sizeof(Float));

	  SingleToDoubleField doubler(gparity_X, gparity_Y, 2, double_wh, new_wh + hit*hit_stride_1f);
	  doubler.Run();
	}
	//the memory location for wh in the 2f version is a factor of 2 smaller, hence we must deallocate and reallocate
	sfree(prop_2f.wh);
	prop_2f.wh = (Vector*)new_wh;

	//memcpy( (void*)prop_2f.get_wh(), (void*)new_wh, wh_bytes);
	//pfree(new_wh);
	pfree(double_wh);
      }
    }

    //FFTW fields

    if(prop_2f.get_v_fftw(0)!=NULL){

      //Assumes just 1 set of v and w (i.e. no strange quark)
      //v_fftw   size a2a_prop->get_nvec()
      for(int i=0;i< prop_2f.get_nvec(); i++){
	//On each index is a CPS fermion
	int vsz = 24*four_vol*sizeof(Float);
	Float* new_v = (Float*)pmalloc(vsz);
	SingleToDoubleField doubler(gparity_X, gparity_Y, 24, (Float*)prop_2f.get_v_fftw(i), new_v);
	doubler.Run();
	memcpy( (void*)prop_2f.get_v_fftw(i), (void*)new_v, vsz);
	pfree(new_v);
      }
      //wl_fftw[0] (light part)  size a2a_prop->get_nl()
      for(int i=0;i< prop_2f.get_nl(); i++){
	//On each index is a CPS fermion
	int vsz = 24*four_vol*sizeof(Float);
	Float* new_v = (Float*)pmalloc(vsz);
	SingleToDoubleField doubler(gparity_X, gparity_Y, 24, (Float*)prop_2f.get_wl_fftw(i), new_v);
	doubler.Run();
	memcpy( (void*)prop_2f.get_wl_fftw(i), (void*)new_v, vsz);
	pfree(new_v);
      }
      //wh_fftw[0] is 12*a2a_prop->get_nhits() stacked fermion vectors (factor of 2 for flavor dilution)
      for(int i=0;i< 12*prop_2f.get_nhits() * (prop_2f.a2a.dilute_flavor ? 2 : 1) ; i++){
	int off = i*24*four_vol; //(float units)
	Float* fld = (Float*)prop_2f.get_wh_fftw() + off;

	int vsz = 24*four_vol*sizeof(Float);
	Float* new_v = (Float*)pmalloc(vsz);
	SingleToDoubleField doubler(gparity_X, gparity_Y, 24, fld, new_v);
	doubler.Run();
	memcpy( fld, (void*)new_v, vsz);
	pfree(new_v);
      }     

    }
  }


  //Compare 1f 2f A2A propagator components (v and w). 2f propagator should have been converted to 1f layout using the above
  static void compare_1f_2f_A2A_prop(A2APropbfm &prop_2f, A2APropbfm &prop_1f){
    //CALL IN 1F CONTEXT
    int four_vol = GJP.VolNodeSites(); 
    printf("Comparing 1f and 2f A2A propagators\n");
    //Compare wl
    for(int l = 0; l < prop_2f.get_nl(); ++l){
      printf("\n\n1f 2f A2A prop starting wl comparison, mode %d\n",l);
      int wlsz = 24*four_vol;
      Float* wl_2f = (Float*)prop_2f.get_wl(l);
      Float* wl_1f = (Float*)prop_1f.get_wl(l);

      bool fail = false;
      for(int i=0;i<wlsz;++i)
	if( fabs(wl_2f[i]-wl_1f[i]) > 1e-12 ){
	  printf("1f 2f A2A prop wl comparison fail mode %d, offset %d.  1f: %.13e     2f: %.13e\n",l,i,wl_1f[i],wl_2f[i]);
	  fail = true;
	}
      if(fail){
	printf("1f 2f A2A prop wl comparison failed on mode %d\n",l); exit(-1);
      }
    }
    //Compare wh
    for(int hit = 0; hit < prop_2f.get_nhits(); hit++){
      printf("\n\n1f 2f A2A prop starting wh comparison, hit %d\n",hit);
      int whsz = 2*four_vol;

      int hit_off = whsz*hit;

      Float* wh_2f = (Float*)prop_2f.get_wh() + hit_off;
      Float* wh_1f = (Float*)prop_1f.get_wh() + hit_off;

      bool fail = false;
      for(int i=0;i<whsz;++i){
	int rem = i;
	int reim = rem % 2; rem/=2;
	int pos[4];
	for(int xx=0;xx<4;xx++){ pos[xx] = rem % GJP.NodeSites(xx); rem /= GJP.NodeSites(xx); }
	
	if( fabs(wh_2f[i]-wh_1f[i]) > 1e-12 ){
	  printf("1f 2f A2A prop wh comparison fail hit %d, offset %d (reim %d, pos %d %d %d %d).  1f: %.13e     2f: %.13e\n",hit,i,reim,pos[0],pos[1],pos[2],pos[3],wh_1f[i],wh_2f[i]);
	  fail = true;
	}else printf("1f 2f A2A prop wh comparison pass hit %d, offset %d (reim %d, pos %d %d %d %d).  1f: %.13e     2f: %.13e\n",hit,i,reim,pos[0],pos[1],pos[2],pos[3],wh_1f[i],wh_2f[i]);
      }
      if(fail){
	printf("1f 2f A2A prop wh comparison failed on hit %d\n",hit); exit(-1);
      }
    }  
    //Compare v
    for(int vec = 0; vec < prop_2f.get_nvec(); ++vec){
      printf("\n\n1f 2f A2A prop starting v comparison, vec %d\n",vec);
      int vsz = 24*four_vol;  
      Float* v_2f = (Float*)prop_2f.get_v(vec);
      Float* v_1f = (Float*)prop_1f.get_v(vec);

      bool fail = false;
      for(int i=0;i<vsz;++i)
	if( fabs(v_2f[i]-v_1f[i]) > 1e-12 ){
	  printf("1f 2f A2A prop v comparison fail vec %d, offset %d.  1f: %.13e     2f: %.13e\n",vec,i,v_1f[i],v_2f[i]);
	  fail = true;
	}
      if(fail){
	printf("1f 2f A2A prop v comparison failed on vec %d\n",vec); exit(-1);
      }
    }
    printf("Successfully compared 1f and 2f A2A propagators\n");
  }

};

CPS_END_NAMESPACE

void setup_gfix_args(FixGaugeArg &r){
  r.fix_gauge_kind = FIX_GAUGE_COULOMB_T;
  r.hyperplane_start = 0;
  r.hyperplane_step = 1;
  r.hyperplane_num = GJP.TnodeSites();
  r.stop_cond = 1e-06;
  r.max_iter_num = 6000;
}

CPS_START_NAMESPACE

class MesonFieldTesting{
public:
  static void convert_mesonfield_2f_1f(MesonField &mf, const bool &gparity_X, const bool &gparity_Y){
    //EXECUTE IN 1F ENVIRONMENT

    int four_vol = GJP.VolNodeSites();
    //Assumes just 1 set of v and w (i.e. no strange quark)
    //v_fftw   size a2a_prop->get_nvec()
    for(int i=0;i< mf.a2a_prop->get_nvec(); i++){
      //On each index is a CPS fermion
      int vsz = 24*four_vol*sizeof(Float);
      Float* new_v = (Float*)pmalloc(vsz);
      SingleToDoubleField doubler(gparity_X, gparity_Y, 24, (Float*)mf.v_fftw[i], new_v);
      doubler.Run();
      memcpy( (void*)mf.v_fftw[i], (void*)new_v, vsz);
      pfree(new_v);
    }
    //wl_fftw[0] (light part)  size a2a_prop->get_nl()
    for(int i=0;i< mf.a2a_prop->get_nl(); i++){
      //On each index is a CPS fermion
      int vsz = 24*four_vol*sizeof(Float);
      Float* new_v = (Float*)pmalloc(vsz);
      SingleToDoubleField doubler(gparity_X, gparity_Y, 24, (Float*)mf.wl_fftw[0][i], new_v);
      doubler.Run();
      memcpy( (void*)mf.wl_fftw[0][i], (void*)new_v, vsz);
      pfree(new_v);
    }
    //wh_fftw[0] is 12*a2a_prop->get_nhits() stacked fermion vectors
    for(int i=0;i< 12*mf.a2a_prop->get_nhits(); i++){
      int off = i*24*four_vol; //(float units)
      Float* fld = (Float*)mf.wh_fftw[0] + off;

      int vsz = 24*four_vol*sizeof(Float);
      Float* new_v = (Float*)pmalloc(vsz);
      SingleToDoubleField doubler(gparity_X, gparity_Y, 24, fld, new_v);
      doubler.Run();
      memcpy( fld, (void*)new_v, vsz);
      pfree(new_v);
    }       
    //Note 1f and 2f mesonfields themselves have the same layout
  }
  static void compare_fftw_vecs(MesonField &_1f, MesonField &_2f){
    printf("Comparing meson fields\n");
    int four_vol = GJP.VolNodeSites();

    for(int i=0;i< _1f.a2a_prop->get_nl(); i++){
      printf("Comparing wl_fftw for vec idx %d\n",i);
      bool fail(false);
      for(int j=0;j<24*four_vol;j++){
	Float v1f = ((Float*)_1f.wl_fftw[0][i])[j];
	Float v2f = ((Float*)_2f.wl_fftw[0][i])[j];
	if( fabs( v1f  - v2f ) > 1e-12 ){ printf("Failed on wl_fftw mode idx %d at ferm offset %d: %.14e  %.14e\n",i,j,v1f,v2f); fail = true; }
      }
      if(fail){ printf("Comparison of wl_fftw for vec idx %d failed\n",i); exit(-1); }
    }
    //wh_fftw[0] is 12*a2a_prop->get_nhits() stacked fermion vectors
    for(int i=0;i< 12*_1f.a2a_prop->get_nhits(); i++){
      int off = i*24*four_vol; //(float units)
      Float* fl1f = (Float*)_1f.wh_fftw[0] + off;
      Float* fl2f = (Float*)_2f.wh_fftw[0] + off;

      printf("Comparing wh_fftw for hit/spin-color idx %d\n",i);
      bool fail(false);
      for(int j=0;j<24*four_vol;j++){
	Float v1f = fl1f[j];
	Float v2f = fl2f[j];
	if( fabs( v1f  - v2f ) > 1e-12 ){ printf("Failed on wh_fftw hit/spin-color idx %d at ferm offset %d: %.14e  %.14e\n",i,j,v1f,v2f); fail = true; }
      }
      if(fail){ printf("Comparison of wh_fftw for hit/spin-color idx %d failed\n",i); exit(-1); }
    }
    for(int i=0;i< _1f.a2a_prop->get_nvec(); i++){
      printf("Comparing v_fftw for vec idx %d\n",i);
      bool fail(false);
      for(int j=0;j<24*four_vol;j++){
	Float v1f = ((Float*)_1f.v_fftw[i])[j];
	Float v2f = ((Float*)_2f.v_fftw[i])[j];
	if( fabs( v1f  - v2f ) > 1e-12 ){ printf("Failed on v_fftw vec idx %d at ferm offset %d: %.14e  %.14e\n",i,j,v1f,v2f); fail = true; }
      }
      if(fail){ printf("Comparison of v_fftw for vec idx %d failed\n",i); exit(-1); }
    }


    printf("Passed comparison of meson fields\n");
  }
  //Ensure that the old MesonField FFT results are the same as those obtained from new A2APropBfm
  static void compare_fftw_fields(MesonField &mf, A2APropbfm &a2a){
    printf("Comparing FFTW fields between MesonField and A2APropbfm\n");
    int four_vol = GJP.VolNodeSites();
    int gpfac = (GJP.Gparity()?2:1);

    for(int i=0;i< a2a.get_nl(); i++){
      printf("Comparing wl_fftw for vec idx %d\n",i);
      bool fail(false);
      for(int j=0;j< gpfac*24*four_vol;j++){
	Float v1f = ((Float*)mf.wl_fftw[0][i])[j];
	Float v2f = ((Float*)a2a.get_wl_fftw(i))[j];
	if( fabs( v1f  - v2f ) > 1e-12 ){ printf("Failed on wl_fftw mode idx %d at ferm offset %d: %.14e  %.14e\n",i,j,v1f,v2f); fail = true; }
      }
      if(fail){ printf("Comparison of wl_fftw for vec idx %d failed\n",i); exit(-1); }
    }
    for(int i=0;i< 12* a2a.get_nhits(); i++){
      int off = i*gpfac*24*four_vol; //(float units)
      Float* fl1f = (Float*)mf.wh_fftw[0] + off;
      Float* fl2f = (Float*)a2a.get_wh_fftw() + off;

      printf("Comparing wh_fftw for hit/spin-color idx %d\n",i);
      bool fail(false);
      for(int j=0;j<gpfac*24*four_vol;j++){
	Float v1f = fl1f[j];
	Float v2f = fl2f[j];
	if( fabs( v1f  - v2f ) > 1e-12 ){ printf("Failed on wh_fftw hit/spin-color idx %d at ferm offset %d: %.14e  %.14e\n",i,j,v1f,v2f); fail = true; }
      }
      if(fail){ printf("Comparison of wh_fftw for hit/spin-color idx %d failed\n",i); exit(-1); }
    }
    for(int i=0;i< a2a.get_nvec(); i++){
      printf("Comparing v_fftw for vec idx %d\n",i);
      bool fail(false);
      for(int j=0;j<gpfac*24*four_vol;j++){
	Float v1f = ((Float*)mf.v_fftw[i])[j];
	Float v2f = ((Float*)a2a.get_v_fftw(i))[j];
	if( fabs( v1f  - v2f ) > 1e-12 ){ printf("Failed on v_fftw vec idx %d at ferm offset %d: %.14e  %.14e\n",i,j,v1f,v2f); fail = true; }
      }
      if(fail){ printf("Comparison of v_fftw for vec idx %d failed\n",i); exit(-1); }
    }
    printf("Passed comparison of FFTW fields between MesonField and A2APropbfm\n");
  }
  
  //Compare the 1f and 2f mesonfield object MesonField.mf. Has the same memory layout for the 2 implementations so no conversion is necessary
  static void compare_mf_ll(MesonField &mf_1f, MesonField &mf_2f){
    int t_size = GJP.TnodeSites()*GJP.Tnodes();
    int size = mf_1f.nvec[0]*(mf_1f.nl[0]+12*mf_1f.nhits[0])*t_size*2;
    bool fail(false);
    for(int i=0;i<size;i++){
      int rem = i;
      int reim = rem % 2; rem/=2;
      int glb_t = rem % t_size; rem/=t_size;
      int mat_j = rem % mf_1f.nvec[0]; rem/=mf_1f.nvec[0];
      int mat_i = rem;

      Float* _1f = (Float*)mf_1f.mf + i;
      Float* _2f = (Float*)mf_2f.mf + i;
      if( fabs(*_1f - *_2f) > 1e-12 ){ if(!UniqueID()) printf("Failed on mf_ll offset %d: %.14e   %.14e\n",i,*_1f,*_2f); fail = true; }      
    }
    if(fail){ printf("Failed mf_ll test\n"); exit(-1); }
    else printf("Passed mf_ll test\n");
  }

  static void compare_source_MesonField2_2f(MesonField &mf, MFsource &src){
    bool fail(false);
    for(int i=0;i<2*GJP.VolNodeSites()/GJP.TnodeSites();i++){
      Float* _orig = (Float*)mf.src + i;
      Float* _new = (Float*)src.src + i;
      if( fabs(*_orig - *_new) > 1e-12 ){ if(!UniqueID()) printf("Failed on source comparison offset %d: %.14e   %.14e\n",i,*_orig,*_new); fail = true; }
      else printf("Passed on source comparison offset %d: %.14e   %.14e\n",i,*_orig,*_new);
    }
    if(fail){ printf("Failed source comparison test\n"); exit(-1); }
    else printf("Passed source comparison test\n");
  }

  static void compare_mf_ll_MesonField2_2f(MesonField &mf_2f, MesonField2 &mf2_2f){
    int t_size = GJP.TnodeSites()*GJP.Tnodes();
    int size = mf_2f.nvec[0]*(mf_2f.nl[0]+12*mf_2f.nhits[0])*t_size*2;
    bool fail(false);
    for(int i=0;i<size;i++){
      int rem = i;
      int reim = rem % 2; rem/=2;
      int glb_t = rem % t_size; rem/=t_size;
      int mat_i = rem % mf_2f.nvec[0]; rem/=mf_2f.nvec[0];
      int mat_j = rem;

      Float* _orig = (Float*)mf_2f.mf + i;
      Float* _new = (Float*)mf2_2f.mf + i;
      if( fabs(*_orig - *_new) > 1e-12 ){ if(!UniqueID()) printf("Failed on mf_ll_MesonField2_2f offset %d (i %d j %d t %d reim %d): %.14e   %.14e\n",i,mat_i,mat_j,glb_t,reim,*_orig,*_new); fail = true; }      
      else if(!UniqueID()) printf("Passed on mf_ll_MesonField2_2f offset %d (i %d j %d t %d reim %d): %.14e   %.14e\n",i,mat_i,mat_j,glb_t,reim,*_orig,*_new);
    }
    if(fail){ printf("Failed mf_ll_MesonField2_2f test\n"); exit(-1); }
    else printf("Passed mf_ll_MesonField2_2f test\n");
  }
  static void compare_mf_sl_MesonField2_2f(MesonField &mf_2f, MesonField2 &mf2_2f){
    int t_size = GJP.TnodeSites()*GJP.Tnodes();
    int size = mf_2f.nvec[0]*(mf_2f.nl[1]+12*mf_2f.nhits[1])*t_size*2;
    bool fail(false);
    for(int i=0;i<size;i++){
      int rem = i;
      int reim = rem % 2; rem/=2;
      int glb_t = rem % t_size; rem/=t_size;
      int mat_i = rem % mf_2f.nvec[0]; rem/=mf_2f.nvec[0];
      int mat_j = rem;

      Float* _orig = (Float*)mf_2f.mf_sl + i;
      Float* _new = (Float*)mf2_2f.mf + i;
      if( fabs(*_orig - *_new) > 1e-12 ){ if(!UniqueID()) printf("Failed on mf_sl_MesonField2_2f offset %d (i %d j %d t %d reim %d): %.14e   %.14e\n",i,mat_i,mat_j,glb_t,reim,*_orig,*_new); fail = true; }      
      else if(!UniqueID()) printf("Passed on mf_sl_MesonField2_2f offset %d (i %d j %d t %d reim %d): %.14e   %.14e\n",i,mat_i,mat_j,glb_t,reim,*_orig,*_new);
    }
    if(fail){ printf("Failed mf_sl_MesonField2_2f test\n"); exit(-1); }
    else printf("Passed mf_sl_MesonField2_2f test\n");
  }
  static void compare_mf_ls_MesonField2_2f(MesonField &mf_2f, MesonField2 &mf2_2f){
    int t_size = GJP.TnodeSites()*GJP.Tnodes();
    int size = mf_2f.nvec[1]*(mf_2f.nl[0]+12*mf_2f.nhits[0])*t_size*2;
    bool fail(false);
    for(int i=0;i<size;i++){
      int rem = i;
      int reim = rem % 2; rem/=2;
      int glb_t = rem % t_size; rem/=t_size;
      int mat_i = rem % mf_2f.nvec[1]; rem/=mf_2f.nvec[1];
      int mat_j = rem;

      Float* _orig = (Float*)mf_2f.mf_ls + i;
      Float* _new = (Float*)mf2_2f.mf + i;
      if( fabs(*_orig - *_new) > 1e-12 ){ if(!UniqueID()) printf("Failed on mf_ls_MesonField2_2f offset %d (i %d j %d t %d reim %d): %.14e   %.14e\n",i,mat_i,mat_j,glb_t,reim,*_orig,*_new); fail = true; }      
      else if(!UniqueID()) printf("Passed on mf_ls_MesonField2_2f offset %d (i %d j %d t %d reim %d): %.14e   %.14e\n",i,mat_i,mat_j,glb_t,reim,*_orig,*_new);
    }
    if(fail){ printf("Failed mf_ls_MesonField2_2f test\n"); exit(-1); }
    else printf("Passed mf_ls_MesonField2_2f test\n");
  }


  static void compare_mf_ll_MesonField2_2f_1f(MesonField2 &mf2_2f, MesonField2 &mf2_1f){
    int t_size = GJP.TnodeSites()*GJP.Tnodes();
    int size = mf2_2f.nvec[0]*(mf2_2f.nl[0]+mf2_2f.dilute_size*mf2_2f.nhits[0])*t_size*2;
    bool fail(false);
    if(mf2_2f.dilute_flavor) printf("Flavour dilution\n");
    else printf("No flavour dilution\n");

    for(int i=0;i<size;i++){
      int rem = i;
      int reim = rem % 2; rem/=2;
      int glb_t = rem % t_size; rem/=t_size;
      int mat_i = rem % mf2_2f.nvec[0]; rem/=mf2_2f.nvec[0];
      int mat_j = rem;

      Float* _orig = (Float*)mf2_2f.mf + i;
      Float* _new = (Float*)mf2_1f.mf + i;
      if( fabs(*_orig - *_new) > 1e-12 ){ if(!UniqueID()) printf("Failed on mf_ll_MesonField2_2f_1f offset %d (i %d j %d t %d reim %d): %.14e   %.14e\n",i,mat_i,mat_j,glb_t,reim,*_orig,*_new); fail = true; }      
      else if(!UniqueID()) printf("Passed on mf_ll_MesonField2_2f_1f offset %d (i %d j %d t %d reim %d): %.14e   %.14e\n",i,mat_i,mat_j,glb_t,reim,*_orig,*_new);
    }
    if(fail){ printf("Failed mf_ll_MesonField2_2f_1f test\n"); exit(-1); }
    else printf("Passed mf_ll_MesonField2_2f_1f test\n");
  }

  static void compare_kaon_corr_MesonField2_2f(std::complex<double> *kaoncorr_orig, CorrelationFunction &kaoncorr_mf2){
    //std::complex<double> *kaoncorr_orig = (std::complex<double> *)pmalloc(sizeof(std::complex<double>)*GJP.Tnodes()*GJP.TnodeSites());
    bool fail(false);

    const char* thr = kaoncorr_mf2.threadType() == CorrelationFunction::THREADED ? "threaded" : "unthreaded";

    for(int t=0;t<GJP.Tnodes()*GJP.TnodeSites();++t){
      double *_orig = (double*) &kaoncorr_orig[t];
      double *_new = (double*) &kaoncorr_mf2(0,t);
      for(int i=0;i<2;i++){	
	if( fabs(_orig[i] - _new[i]) > 1e-12 ){ if(!UniqueID()) printf("Failed on compare_kaon_corr_MesonField2_2f t = %d i = %d: %.14e   %.14e\n",t,i,_orig[i],_new[i]); fail = true; }      
      }
    }
    if(fail){ printf("Failed compare_kaon_corr_MesonField2_2f (%s) test\n",thr); exit(-1); }
    else printf("Passed compare_kaon_corr_MesonField2_2f (%s) test\n",thr);
  }


#define SETUP_ARRAY(OBJ,ARRAYNAME,TYPE,SIZE)	\
  OBJ . ARRAYNAME . ARRAYNAME##_len = SIZE; \
  OBJ . ARRAYNAME . ARRAYNAME##_val = new TYPE [SIZE]

#define ELEM(OBJ,ARRAYNAME,IDX) OBJ . ARRAYNAME . ARRAYNAME##_val[IDX]

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

  static void wallsource_amplitude_nogparity(Lattice &lattice, Lanczos_5d<double> &eig){
    //Use a MesonField2 object with flavor dilution and NORAND high modes for this comparison
    //Source type must be a box source filling the entire lattice 3-volume
    //Compare to traditional wall source contraction

    // With zero low modes,    w_a^i(x) = \delta_a^{i}   \forall x   where i runs from 0 to 11

    // \sum_i w_a^i(x) [w_b^i(y)]*  =  \sum_{ij} w_a^i(x) \delta_{ij} [w_b^j(y)]*
    //                              =  \delta_a^{i}  \delta_{ij}  \delta_b^{j}
    //                              =  \delta_{ab}

    // for all x,y.

    // v_a^i(x) = \sum_y D^{-1}_{ab}(x,y) w_b^i(y) = \sum_y D^{-1}_{ai}(x,y)

    // which is just a wall source propagator.

    // \sum_i v_a^i(x) [w_b^i(y)]* = \sum_y D^{-1}_{ai}(x,y)\delta_b^{i} = \sum_y D^{-1}_{ab}(x,y)

    // is also a wall source propagator.

    // Using   \sum_i v_a^i(x) [w_b^i(y)]* = \sum_y D^{-1}_{ai}(x,y)\delta_b^{i} = \sum_y D^{-1}_{ab}(x,y)   we calculate a generic meson correlation function:

    // \sum_{x,x'}\sum_{y,y'} f(|x-x'|) D^{-1}_{ab}(x,y) M_{bc} f(|y-y'|) D^{-1}_{cd}(y',x') M_{da}
    // =
    // \sum_{x,x'}\sum_{y,y'} f(|x-x'|) [\sum_i v_a^i(x) [w_b^i(y)]* ]  M_{bc} f(|y-y'|)  [\sum_j v_c^j(y') [w_d^j(x')]* ] M_{da}
    // =
    // \sum_{ij} \sum_{x,x'}\sum_{y,y'} [ f(|x-x'|) [w_d^j(x')]* M_{da} v_a^i(x)  ]    [ f(|y-y'|) [w_b^i(y)]* M_{bc} v_c^j(y')  ]

    // then use

    // \sum_{x,x'} [ f(|x-x'|) [w_d^j(x')]* M_{da} v_a^i(x)  ]
    // =
    // \sum_{x,x',z} [ f(|z|)\delta(z-x'+x) [w_d^j(x')]* M_{da} v_a^i(x)  ]
    // =

    // \sum_p [ \sum_z f(|z|) e^{ipz} ] [\sum_x' e^{ipx'} w_d^j(x')]* M_{da}  [\sum_x v_a^i(x) e^{ipx} ]

    //Gauge fix lattice if not done already
    CommonArg c_arg;
    FixGaugeArg gfix_arg;
    setup_gfix_args(gfix_arg);
    AlgFixGauge fix_gauge(lattice,&c_arg,&gfix_arg);
    bool free_gfix(false);
    if(lattice.FixGaugeKind()==FIX_GAUGE_NONE){
      fix_gauge.run();
      free_gfix = true; //will be freed at end of method
    }
   
    bool test_with_gauge_fixing = true;

    int neig =0;

    //Generate the a2a prop
    A2APropbfm* prop_f;
    A2APropbfmTesting::a2a_prop_gen(prop_f, &lattice, eig, 1,NORAND,neig);
    if(!test_with_gauge_fixing) prop_f->do_gauge_fix(false);

    printf("Wall source prop with %d eigenvectors\n",prop_f->get_args().nl);

    {
      //Check wh is just a wall source on each time-slice
      //wh_offset = wh_id * nodesites + (x + Lx*(y+Ly*(z+Lz*t)))    (wh_id is a hit index)
      Float* wh = (Float*) prop_f->get_wh();
      bool fail(false);
      for(int i=0;i<GJP.VolNodeSites();i++){
	if( wh[2*i]!=1.0 || wh[2*i+1]!=0.0 ){ printf("Wh check fail site %d: %.14e %.14e\n",i,wh[2*i], wh[2*i+1]);  fail=true; }
      }
      if(fail){ printf("Wh check failed\n"); exit(-1); }
      else printf("Wh check passed\n"); 
    }


    //Calculate mesonfield with wall source, i.e. box source of size L^3
    MFqdpMatrix structure_wdagv(MFstructure::W, MFstructure::V, true, false,15); //gamma^5 spin
    MFBasicSource source(MFBasicSource::BoxSource,GJP.Xnodes()*GJP.XnodeSites());
    
    {
      //Check source structure is a delta function p=0 (Fourier transform of a uniform wall source)
      bool fail(false);
      for(int i=0;i<GJP.Xnodes()*GJP.XnodeSites()*GJP.Ynodes()*GJP.YnodeSites()*GJP.Znodes()*GJP.ZnodeSites();i++){
	std::complex<Float> val = source(i);	
	if(i==0 && (val.real()!=1.0 || val.imag()!=0.0)){ printf("Source failed %d : %.14e %.14e expected 1,0\n",i,val.real(),val.imag()); fail=true; }
	else if(i!=0 && (val.real()!=0.0 || val.imag()!=0.0)){ printf("Source failed %d : %.14e %.14e expected 0,0\n",i,val.real(),val.imag()); fail=true; }
      }
      if(fail){ printf("Source check failed\n"); exit(-1); }
      else printf("Source check passed\n");
    }

    MesonField2 mf(*prop_f,*prop_f, structure_wdagv, source);

    //Calculate correlation function
    CorrelationFunction corr_a2a("",1, CorrelationFunction::THREADED);

    int tsrc_con = 0;
    MesonField2::contract_specify_tsrc(mf,mf,0, tsrc_con, corr_a2a);
    corr_a2a.sumThreads();

    //Generate wall source propagators (yeah I know it should be a momentum source, but for this comparison with a single translationally non-invariant lattice it should not matter)
    PropManager::clear();

    int Lt = GJP.Tnodes()*GJP.TnodeSites();

    JobPropagatorArgs prop_args;
    SETUP_ARRAY(prop_args,props,PropagatorArg,Lt);

    std::string prop_names_srct[Lt];

    for(int t=0;t<Lt;t++){
      PropagatorArg &parg = prop_args.props.props_val[t];
    
      std::ostringstream os; os << "prop_"<<t;
      prop_names_srct[t] = os.str();

      parg.generics.tag = const_cast<char*>(prop_names_srct[t].c_str());
      parg.generics.type = QPROPW_TYPE;
      parg.generics.mass = 0.01; //needs to be the same as the mass in the a2a_arg function!
      parg.generics.bc[0] = GJP.Xbc();
      parg.generics.bc[1] = GJP.Ybc();
      parg.generics.bc[2] = GJP.Zbc();
      parg.generics.bc[3] = GJP.Tbc();

      //int sz = test_with_gauge_fixing ? 3 : 2;
      int sz = 2;

      SETUP_ARRAY(parg,attributes,AttributeContainer,sz);
      
      ELEM(parg,attributes,0).type = WALL_SOURCE_ATTR;
      WallSourceAttrArg &srcarg = ELEM(parg,attributes,0).AttributeContainer_u.wall_source_attr;
      srcarg.t = t;

      ELEM(parg,attributes,1).type = CG_ATTR;
      CGAttrArg &cgattr = ELEM(parg,attributes,1).AttributeContainer_u.cg_attr;
      cgattr.max_num_iter = 6000;
      cgattr.stop_rsd = 1e-10;
      cgattr.true_rsd = 1e-10;
      
      //In A2A the v vector is D^{-1} w  where w does not contain the gauge fixing matrices. Only after v is calculated are the gauge fixing matrices applied to w (after it is spin-color diluted)
      //and the intrinsic delta function in space ensures the meson field w*v is gauge invariant
      
      // if(test_with_gauge_fixing){
      // 	ELEM(parg,attributes,2).type = GAUGE_FIX_ATTR;
      // 	GaugeFixAttrArg &gfattr = ELEM(parg,attributes,2).AttributeContainer_u.gauge_fix_attr;
      // 	gfattr.gauge_fix_src = 1;
      // 	gfattr.gauge_fix_snk = 0;
      // }
    }

    PropManager::setup(prop_args);   
    PropManager::calcProps(lattice);



    for(int tsrc=0;tsrc<Lt;tsrc++){
      //Check propagators match. We have set up v and w such that v are just wall source propagators
      //There should be  Lt * 12 vectors. We want to pick out the one originating from timeslice tsrc
      //Mapping is
      //i-a2a.nl = (sc_id + 3*4/src_width * t_id +  3*4/src_width*Lt * wh_id)
      //where a2a.nl=0 here,  src_width=1  and nhits=1
      QPropWcontainer &pcon = PropManager::getProp(prop_names_srct[tsrc].c_str()).convert<QPropWcontainer>();

      bool fail = false;
      for(int x=0;x<GJP.VolNodeSites();++x){
	WilsonMatrix & qpw_site = pcon.getProp(lattice).SiteMatrix(x);

	for(int sc_id=0;sc_id<12;sc_id++){
	  int scol=sc_id/3;
	  int ccol=sc_id%3;

	  int i = sc_id + 12*tsrc;
	  Float* v_i = (Float*)prop_f->get_v(i) + 24*x;
	
	  //Here sc_id is the column spin-color index
	  for(int rowidx=0;rowidx<12;rowidx++){
	    int srow = rowidx/3;
	    int crow = rowidx %3;
	    Float* qpw_cmp = (Float*)& qpw_site(srow,crow,scol,ccol);
	    for(int reim=0;reim<2;reim++){
	      Float qpw_val = qpw_cmp[reim];
	      Float v_val = v_i[reim + 2*rowidx];
	    
	      //NOTE THERE IS A NORMALIZATION DIFFERENCE! QPW = (5-M5) A2A for DWF
	      qpw_val/= 3.2;

	      if( !ratio_diff_match(qpw_val,v_val,1e-3) ){
		printf("Wall source tsrc %d v and qpropw test fail x=%d sc_col=%d sc_row=%d reim=%d: %.14e %.14e, ratio %.14e\n",tsrc,x,sc_id,rowidx,reim,qpw_val,v_val,qpw_val/v_val);
		fail=true;
	      }
	    }
	  }
	}
      }
      if(fail){
	printf("Wall source tsrc %d v and qpropw test failed\n",tsrc); exit(-1);
      }else printf("Wall source tsrc %d v and qpropw test passed\n",tsrc);    
    }

    std::vector<Float> mom(3,0.0);
    FourierProp<WilsonMatrix> ftprop_calc;
    ftprop_calc.add_momentum(mom);
    if(!test_with_gauge_fixing) ftprop_calc.gaugeFixSink(false);

    std::vector<WilsonMatrix> ftprop_srct[Lt];

    for(int t=0;t<Lt;t++)
      ftprop_srct[t] = ftprop_calc.getFTProp(lattice,mom,prop_names_srct[t].c_str());

    if(!prop_f->fft_vw_computed()){ printf("FFTW vectors have not been computed, why?\n"); exit(-1); }
  
    {
      //Independently calculate zero momentum FTprop and confirm FourierProp result
      WilsonMatrix p0propft[Lt][Lt]; //src, snk
      for(int tsrc = 0; tsrc <Lt; tsrc++) for(int tsnk=0;tsnk<Lt;tsnk++) p0propft[tsrc][tsnk] = 0.0;

      for(int tsrc = 0; tsrc <Lt; tsrc++){
	for(int x=0;x<GJP.VolNodeSites();x++){
	  int tsnk = x/(GJP.VolNodeSites()/GJP.TnodeSites()) + GJP.TnodeSites()*GJP.TnodeCoor();
	  
	  WilsonMatrix site = PropManager::getProp(prop_names_srct[tsrc].c_str()).convert<QPropWcontainer>().getProp(lattice).SiteMatrix(x);
	  if(test_with_gauge_fixing) site.LeftTimesEqual(*lattice.FixGaugeMatrix(x));

	  p0propft[tsrc][tsnk] += site;
	}
      }
      for(int tsrc = 0; tsrc <Lt; tsrc++) 
	for(int tsnk=0;tsnk<Lt;tsnk++){
	  Float* w = (Float*)p0propft[tsrc][tsnk].ptr(); //returns Rcomplex*
	  static const int size = 2*12*12;
	  slice_sum(w, size, 99);
	}

      for(int tsrc=0;tsrc<Lt;tsrc++){
	bool fail = false;

	for(int tsnk=0;tsnk<Lt;tsnk++){
	  for(int sc_id=0;sc_id<12;sc_id++){
	    int scol=sc_id/3;
	    int ccol=sc_id%3;

	    //Here sc_id is the column spin-color index
	    for(int rowidx=0;rowidx<12;rowidx++){
	      int srow = rowidx/3;
	      int crow = rowidx %3;

	      Float* manualgf_ft = (Float*)pmalloc( 2 *sizeof(Float) );
	      for(int p=0;p<2;++p) manualgf_ft[p] = 0.0;

	      for(int x=0;x<GJP.VolNodeSites();x++){
		WilsonMatrix site = PropManager::getProp(prop_names_srct[tsrc].c_str()).convert<QPropWcontainer>().getProp(lattice).SiteMatrix(x);
		int t_glb = x/(GJP.VolNodeSites()/GJP.TnodeSites()) + GJP.TnodeCoor()*GJP.TnodeSites();
		if(t_glb!=tsnk) continue;

		if(test_with_gauge_fixing) site.LeftTimesEqual(*lattice.FixGaugeMatrix(x));
		
		manualgf_ft[0] += site(srow,crow,scol,ccol).real();
		manualgf_ft[1] += site(srow,crow,scol,ccol).imag();
	      }
	      slice_sum(manualgf_ft, 2*Lt, 99);

	      Float* qpw_cmp = (Float*)& ftprop_srct[tsrc][tsnk](srow,crow,scol,ccol);
	      Float* p0manual_cmp = (Float*)& p0propft[tsrc][tsnk](srow,crow,scol,ccol);

	      for(int reim=0;reim<2;reim++){
		Float qpw_val = qpw_cmp[reim];
		Float v_val = p0manual_cmp[reim];
		Float manual_val = manualgf_ft[reim];

		if( fabs(qpw_val-v_val) > 1e-12 ){
		  printf("FFT propw and manually computed p0 tsrc %d test fail tsnk=%d sc_col=%d sc_row=%d reim=%d: qpw %.14e manual %.14e, ratio %.14e\n",tsrc,tsnk,sc_id,rowidx,reim,qpw_val,v_val,qpw_val/v_val);
		  fail=true;
		}else if( fabs(qpw_val-manual_val)>1e-12 ){
		  printf("FFT propw and manually computed p0 tsrc %d component sum test fail tsnk=%d sc_col=%d sc_row=%d reim=%d: qpw %.14e manual %.14e, ratio %.14e\n",tsrc,tsnk,sc_id,rowidx,reim,qpw_val,manual_val,qpw_val/manual_val);
		  fail=true;
		}
	      }
	    }
	  }
	}
      
	if(fail){
	  printf("FFT propw and manually computed p0 tsrc %d test failed\n",tsrc); exit(-1);
	}else printf("FFT propw and manually computed p0 tsrc %d test passed\n",tsrc);   
      }
    }
  
    if(test_with_gauge_fixing){
      //Test A2A prop gauge fixing
      for(int tsrc=0;tsrc<Lt;tsrc++){
	for(int src_sc = 0; src_sc < 12 ; src_sc ++){
	  Float* a2agf = (Float*)pmalloc( 24*GJP.VolNodeSites() *sizeof(Float) );
	  Float* manualgf = (Float*)pmalloc( 24*GJP.VolNodeSites() *sizeof(Float) );
	  Float* ungf = (Float*)pmalloc( 24*GJP.VolNodeSites() *sizeof(Float) );
	  Float* ungf_ft = (Float*)pmalloc( 24*Lt *sizeof(Float) );
	  Float* manualgf_ft = (Float*)pmalloc( 24*Lt *sizeof(Float) );
	  for(int p=0;p<24*Lt;++p){
	    ungf_ft[p] = 0.0;
	    manualgf_ft[p] = 0.0;
	  }
	  
	  for(int x=0;x<GJP.VolNodeSites();x++){
	    WilsonMatrix site = PropManager::getProp(prop_names_srct[tsrc].c_str()).convert<QPropWcontainer>().getProp(lattice).SiteMatrix(x);

	    int t_glb = x/(GJP.VolNodeSites()/GJP.TnodeSites()) + GJP.TnodeCoor()*GJP.TnodeSites();

	    for(int snk_sc=0;snk_sc<12;snk_sc++){
	      int off = 24*x + 2*snk_sc;
	      *((Rcomplex*)(ungf+off)) = site(snk_sc/3,snk_sc%3, src_sc/3, src_sc%3);

	      ungf_ft[24*t_glb + 2*snk_sc] +=  ungf[off];
	      ungf_ft[24*t_glb + 2*snk_sc+1] +=  ungf[off+1];
	    }
	    int off = 24*x;
	    memcpy( (void*)(manualgf+off), (void*)(ungf+off), 24*sizeof(Float));
	    memcpy( (void*)(a2agf+off), (void*)(ungf+off), 24*sizeof(Float));
	    
	    const Matrix &V = *lattice.FixGaugeMatrix(x);
	    for(int s=0;s<4;s++){
	      Vector tmp = *((Vector*)(manualgf+off) + s);
	      ((Vector*)(manualgf+off) + s)->DotXEqual(V,tmp);
	    }

	    //IS THIS SOMEHOW DIFFERING FROM THE WILSON MATRIX MULTIPLY????? Yes it does appear so.... why??? OK. so it was a bug in WilsonMatrix that has been fixed in the latest CPS but clearly not in this version!
	    bool fail=false;
	    WilsonMatrix site_copy(site);

	    site_copy.LeftTimesEqual(V); //this one looks right! Check the other one
	    
	    Float * dotxeq_check_out = (Float*)pmalloc( 24 * sizeof(Float) );
	    Float * dotxeq_check_in = (Float*)pmalloc( 24 * sizeof(Float) );
	    memcpy( (void*)dotxeq_check_in, (void*)(ungf+off), 24*sizeof(Float));
	    for(int ii=0;ii<24;ii++) dotxeq_check_out[ii] = 0.0;

	    for(int s_i=0;s_i<4;s_i++){
	      Rcomplex* out_si = (Rcomplex*)(dotxeq_check_out + 6*s_i);
	      Rcomplex* in_si = (Rcomplex*)(dotxeq_check_in + 6*s_i);
	      for(int ii=0;ii<3;ii++){
		for(int jj=0;jj<3;jj++){
		  out_si[ii] += V(ii, jj) * in_si[jj];
		}
	      }
	    }
	    
	    for(int s=0;s<4;s++){
	      Vector *mgf = (Vector*)(manualgf+off) + s;
	      Vector *dotxcheck_s = (Vector*)dotxeq_check_out + s;
	      Vector wmm;
	      Float* wmmf = (Float*)&wmm;
	      for(int c=0;c<3;c++){
		wmmf[2*c] = site_copy(s,c,src_sc/3,src_sc%3).real();
		wmmf[2*c+1] = site_copy(s,c,src_sc/3,src_sc%3).imag();
	      }
	      for(int ss=0;ss<6;ss++){
		Float mf_v = *( (Float*)mgf + ss );
		Float wmm_v = *( wmmf + ss );
		Float dotxcheck_v = *( (Float*)dotxcheck_s + ss );

		if( fabs(mf_v-wmm_v) > 1e-12 ){ printf("GF matrix multiply check failed tsrc=%d src_sc=%d x=%d sink spin=%d col_reim_idx=%d: WMM %.14e  VM %.14e, ratio %.14e, dotXcheck %.14e\n",tsrc,src_sc,x,s,ss,wmm_v,mf_v,wmm_v/mf_v,dotxcheck_v); fail=true; }
	      }
	    }
	    if(fail){ printf("GF matrix multiply check failed\n"); exit(-1); }

	    for(int snk_sc=0;snk_sc<12;snk_sc++){
	      for(int reim=0;reim<2;reim++){
		int off2 = reim + 2*snk_sc + 24*x;
		manualgf_ft[24*t_glb + 2*snk_sc + reim] += manualgf[off2];
		
		printf("x=%d tsrc=%d src_sc=%d snk_sc=%d reim=%d gauge fixed vs not gauge-fixed: %.14e  %.14e\n",x,tsrc,src_sc,snk_sc,reim,ungf[off2],manualgf[off2]);
	      }
	    }
	  }
	  Float* tmp = (Float*)pmalloc( 24*GJP.VolNodeSites() *sizeof(Float) );
	  memcpy((void*)tmp,(void*)a2agf, 24*GJP.VolNodeSites() *sizeof(Float));

	  prop_f->gf_vec(  (Vector*)a2agf,  (Vector*)tmp );

	  bool fail=false;
	  for(int x=0;x<24*GJP.VolNodeSites();x++){
	    if(fabs(a2agf[x]-manualgf[x])>1e-12){ printf("A2A manual prop gfix test fail i=%d: got %.14e expect %.14e\n",x,a2agf[x],manualgf[x]); fail=true; }
	  }
	  if(fail){ printf("A2A manual prop gfix test fail\n"); exit(-1); }
	  else printf("A2A manual prop gfix test pass\n");

	  //Try using A2APropbfm FFT code to get p=0 FTprop and compare to FTprop calculated above
	  fftw_init_threads();
	  fftw_plan_with_nthreads(bfmarg::threads);
	  const int fft_dim[3] = { GJP.ZnodeSites() * GJP.Znodes(),
				   GJP.YnodeSites() * GJP.Ynodes(),
				   GJP.XnodeSites() * GJP.Xnodes()};
	  const int t_size = GJP.TnodeSites()*GJP.Tnodes();
	  const int size_4d = GJP.VolNodeSites();
	  const int size_3d_glb = fft_dim[0]*fft_dim[1]*fft_dim[2];
	  const int sc_size = 12;

	  int fftw_alloc_sz = size_3d_glb * t_size * sc_size;
	  fftw_complex *fft_mem = fftw_alloc_complex(fftw_alloc_sz);

	  prop_f->fft_vector( (Vector*)ungf, (Vector*) ungf, fft_dim, fft_mem);

	  prop_f->do_gauge_fix(false);
	  prop_f->fft_vector( (Vector*)manualgf, (Vector*) manualgf, fft_dim, fft_mem);
	  prop_f->do_gauge_fix(true);

	  fftw_free(fft_mem);
	  fftw_cleanup();
	  fftw_cleanup_threads();

	  if(GJP.XnodeCoor()==0 && GJP.YnodeCoor()==0 && GJP.ZnodeCoor()==0){ 	    //we want p=0
	    bool fail = false;
	    for(int tsnk=0; tsnk<Lt;tsnk++){
	      if(tsnk < GJP.TnodeCoor()*GJP.TnodeSites() || tsnk >= (GJP.TnodeCoor()+1)*GJP.TnodeSites() ) continue;
	      int tsnk_loc = tsnk - GJP.TnodeCoor()*GJP.TnodeSites();
	    
	      for(int snk_sc=0;snk_sc<12;snk_sc++){
		Float *a2a_val = ungf + 2*snk_sc + 24*GJP.VolNodeSites()/GJP.TnodeSites()*tsnk_loc;
		Float *qpw_val = (Float*)&ftprop_srct[tsrc][tsnk](snk_sc/3,snk_sc%3, src_sc/3, src_sc%3);
		Float *pregf_val = manualgf + 2*snk_sc + 24*GJP.VolNodeSites()/GJP.TnodeSites()*tsnk_loc;
		Float* mangf_val = manualgf_ft + 2*snk_sc + 24*tsnk;

		for(int reim=0;reim<2;reim++){
		  if( fabs(a2a_val[reim]-qpw_val[reim]) > 1e-07){ printf("FFTprop p=0 ftprop comp test fail tsrc=%d tsnk=%d sc_src=%d sc_snk=%d reim=%d: expect %.14e got %.14e   ratio %.14e  (FFT pre-gf version %.14e  manual FT gf version %.14e)\n",tsrc,tsnk,src_sc,snk_sc,reim,qpw_val[reim], a2a_val[reim], qpw_val[reim]/a2a_val[reim], pregf_val[reim], mangf_val[reim] ); fail=true; }

		}
		
	      }
	    }
	    if(fail){ printf("FFTprop p=0 ftprop comp test fail\n"); exit(-1); }
	    else printf("FFTprop p=0 ftprop comp test passed\n");
	  }

	}
      }
    }

    //Check FT props at p=0 agree
    for(int tsrc=0;tsrc<Lt;tsrc++){
      bool fail = false;

      for(int tsnk=0;tsnk<Lt;tsnk++){
	if( tsnk < GJP.TnodeCoor()*GJP.TnodeSites() || tsnk >= (GJP.TnodeCoor()+1)*GJP.TnodeSites() ) continue; //fftw vecs are redistributed locally after fft

	for(int sc_id=0;sc_id<12;sc_id++){
	  int scol=sc_id/3;
	  int ccol=sc_id%3;

	  int i = sc_id + 12*tsrc;
	  Float* v_i = (Float*)prop_f->get_v_fftw(i) + 24*GJP.VolNodeSites()/GJP.TnodeSites()*tsnk ; //p=0 at offset 0 spatial coord and t temporal
	
	  //Here sc_id is the column spin-color index
	  for(int rowidx=0;rowidx<12;rowidx++){
	    int srow = rowidx/3;
	    int crow = rowidx %3;
	    Float* qpw_cmp = (Float*)& ftprop_srct[tsrc][tsnk](srow,crow,scol,ccol);
	    for(int reim=0;reim<2;reim++){
	      Float qpw_val = qpw_cmp[reim];
	      Float v_val = v_i[reim + 2*rowidx];
	    
	      //NOTE THERE IS A NORMALIZATION DIFFERENCE! QPW = (5-M5) A2A for DWF
	      qpw_val/= 3.2;

	      if( !ratio_diff_match(qpw_val,v_val,1e-3) && fabs(qpw_val-v_val) > 1e-7 ){ //latter because we don't expect elements to be accurate to much more than 1e-07
		printf("FFTW v and qpropw tsrc %d test fail tsnk=%d sc_col=%d sc_row=%d reim=%d: qpw %.14e a2a %.14e, ratio %.14e\n",tsrc,tsnk,sc_id,rowidx,reim,qpw_val,v_val,qpw_val/v_val);
		fail=true;
	      }
	    }
	  }
	}
      }
      
      if(fail){
	printf("FFTW v and qpropw tsrc %d test failed\n",tsrc); exit(-1);
      }else printf("FFTW v and qpropw tsrc %d test passed\n",tsrc);   
    }

    int Vglob = 1;
    for(int i=0;i<3;i++) Vglob*= GJP.Nodes(i)*GJP.NodeSites(i);

    WilsonMatrix FTw[Lt];
    for(int t=0;t<Lt;t++) FTw[t] = 0.0;
    WilsonMatrix _one(0.0);
    for(int i=0;i<12;i++) _one(i/3,i%3, i/3, i%3) = std::complex<Float>(1.0,0.0);
    
    //If gauge fixing we need the FT of w (zero mom) for comparison
    for(int t=0;t<GJP.TnodeSites();t++){
      int t_glb = GJP.TnodeCoor()*GJP.TnodeSites() + t;
      for(int x3d=0;x3d<GJP.VolNodeSites()/GJP.TnodeSites();x3d++){
	int x4d = x3d + GJP.VolNodeSites()/GJP.TnodeSites()*t;
	const Matrix &V =  *lattice.FixGaugeMatrix(x4d);
	WilsonMatrix tmp(_one);
	tmp.LeftTimesEqual(V);

	FTw[t_glb] += tmp;
      }
    }
    slice_sum( (Float*)&FTw, 12*12*2*Lt, 99);
    {
      //Check FFT w field is  (no gauge fixing)  w_a^i = \sum_x e^{ipx} \delta_a^i = Vol3d \delta(p) \delta_a^i
      //                                               = \sum_x e^{ipx} V_ab(x,t) \delta_b^i 

      bool fail=false;
      for(int i=0;i<12;i++){
	Float* wh_fftw = (Float*)prop_f->get_wh_fftw() + 24*GJP.VolNodeSites()*i;
	for(int t=0;t<GJP.TnodeSites();t++){
	  for(int p=0;p<GJP.VolNodeSites()/GJP.TnodeSites();p++){
	    int t_glb = GJP.TnodeCoor()*GJP.TnodeSites() + t;
	    
	    int p_glb[3]; int rem=p;
	    for(int i=0;i<3;i++){ 
	      p_glb[i] = rem % GJP.NodeSites(i) + GJP.NodeCoor(i)*GJP.NodeSites(i); rem/=GJP.NodeSites(i); 
	    }

	    for(int sc=0;sc<12;sc++){
	      for(int reim=0;reim<2;reim++){
		Float* val = wh_fftw + reim + 2*sc + 24*( p+ GJP.VolNodeSites()/GJP.TnodeSites()*t );

		if(!test_with_gauge_fixing){
		  if( p_glb[0]==0 && p_glb[1]==0 && p_glb[2]==0 && sc ==i && reim==0 && *val != Vglob ){ printf("w fail at p=0 on t=%d, i=%d sc=%d: expect 1.0 got %.14e\n",t_glb,i,sc,*val); fail=true; }
		  else if( !(p_glb[0]==0 && p_glb[1]==0 && p_glb[2]==0 && sc ==i && reim==0) && *val != 0.0 ){ printf("w fail at p=0 on t=%d, i=%d sc=%d: expect 0.0 got %.14e\n",t_glb,i,sc,*val); fail=true; }
		}else{
		  //Just do p=0 part
		  if( p_glb[0]==0 && p_glb[1]==0 && p_glb[2]==0){
		    Float expect = ((Float*)&(FTw[t_glb](sc/3,sc%3, i/3,i%3)))[reim];
		    
		    if( fabs(*val - expect)>1e-7 ){ printf("w fail at p=0 on t=%d, i=%d sc=%d: expect %.14e got %.14e\n",t_glb,i,sc,expect,*val); fail=true; }
		  }
		}

	      }
	    }

	  }
	}
      }
      if(fail){ printf("w fail\n"); exit(-1); }
      else printf("w pass\n");
    }

      //Meson field is (no gauge fixing)
      //M^{ij}(t;t') = \sum_{p} \delta(p) [\sum_x e^{-ipx}[w_a^i(x,t)]* ] \gamma^5_{ab} [\sum_y e^{ipy} v_b^j(y,t; t')]   (wall source is zero momentum delta function - tested above)
      //Here i and j both run from 0..11 and the diluted source time-slice is t'
      //          = [\sum_x \delta_a^i] \gamma^5_{ab}  \sum_y v_b^j(y,t; t')
      //          = [\sum_x \delta_a^i] \gamma^5_{ab}  \sum_{y,y'}  D^{-1}_{bc}(y,t; y',t') w_c^i(y',t')
      //          = [\sum_x \delta_a^i] \gamma^5_{ab}  \sum_{y,y'}  D^{-1}_{bc}(y,t; y',t') \delta_c^j
      //          = V \gamma^5_{ib}  \sum_{y,y'}  D^{-1}_{bj}(y,t; y',t')
      //where V is the *spatial* three-volume 



    {
      //Check   [w_a^i]* \gamma^5_ab v_b^j contraction working
      WilsonMatrix g5(0.0);
      for(int s=0;s<4;s++) for(int c=0;c<4;c++) g5(s,c,s,c) = 1.0;
      g5.gr(-5);
      
      bool fail = false;
      int size_4d = GJP.VolNodeSites();
      int size_3d = size_4d/GJP.TnodeSites();
      
      for(int i = 0; i < 12; i++){ // nvec = nl + nhits * Lt * sc_size / width   for v generated independently for each hit, source timeslice and spin-color index
	for(int j = 0; j < 12*Lt; j++){
	  //We can also build the zero momentum part of the contraction from the FTprops about and check that too
	  int t_j = j/12; //v source time
	  int sc_j = j%12;

	  complex<double> mf_ij_repro[Lt]; //indexed by sink time

	  for(int tsnk = 0; tsnk<Lt; tsnk++){
	    WilsonMatrix tmp = ftprop_srct[t_j][tsnk];
	    tmp.gl(-5);
	    if(!test_with_gauge_fixing){
	      tmp *= Float(Vglob)/3.2; //normalise
	    }else{
	      WilsonMatrix wstar = FTw[tsnk];
	      wstar.hconj();

	      tmp.LeftTimesEqual(wstar);
	      tmp *= 1.0/3.2;
	    }

	    mf_ij_repro[tsnk] = tmp(i/3, i%3, sc_j/3, sc_j%3);
	  }
	  for(int x = 0; x < size_4d; x++){	
	    int x_3d = x % size_3d;
	    int t = x / size_3d;
	    
	    int t_glb = t + GJP.TnodeCoor()*GJP.TnodeSites();

	    int rem = x_3d;
	    int pos_glb[3];
	    for(int ii=0;ii<3;ii++){
	      pos_glb[ii] = rem % GJP.NodeSites(ii) + GJP.NodeCoor(ii)*GJP.Nodes(ii);
	      rem/=GJP.NodeSites(ii);
	    }

	    complex<double> *left_vec = prop_f->get_w_fftw(i,x,0);
	    complex<double> *right_vec = prop_f->get_v_fftw(j,x,0);
	    
	    complex<double> con = structure_wdagv.contract_internal_indices(left_vec,right_vec);

	    complex<double> expect(0.0,0.0);
	    for(int a=0;a<12;a++) for(int b=0;b<12;b++) expect +=  conj(left_vec[a]) * g5(a/3,a%3, b/3,b%3) * right_vec[b];

	    for(int reim=0;reim<2;reim++){
	      Float ex = ((Float*)&expect)[reim];
	      Float c = ((Float*)&con)[reim];

	      if( fabs(ex-c) > 1e-12 ){
		printf("Contraction check fail i=%d j=%d x=%d t=%d reim %d: expect %.14e go %.14e  ratio %.14e\n",i,j,x_3d,t,reim,((Float*)&expect)[reim], ((Float*)&con)[reim], ((Float*)&expect)[reim]/((Float*)&con)[reim] );
		fail = true;
	      }
	      if(pos_glb[0]==0 && pos_glb[1]==0 && pos_glb[2] == 0){
		//FOR MULTI-NODE NEED TO SUM OVER NODES FOR THIS COMPARISON
		Float r = ((Float*)&mf_ij_repro[t_glb])[reim];
		//Compare to FTprop result
		if( !ratio_diff_match(r,c,1e-4) && fabs(r-c) > 1e-07 ){
		  printf("Zero mom prop mf repro test fail i=%d j=%d t=%d reim %d: expect %.14e go %.14e  ratio %.14e\n",i,j,t_glb,reim,((Float*)&mf_ij_repro[t_glb])[reim], ((Float*)&con)[reim], ((Float*)&mf_ij_repro[t_glb])[reim]/((Float*)&con)[reim] );
		  fail = true;
		}
	      }

	    }
	  }
	}
      }
      if(fail){ printf("Contraction check fail\n"); exit(-1); }
      else printf("Contraction check pass\n");
	
    }


    {
      bool fail = false;
      for(int tsnk=0;tsnk<Lt;tsnk++){
	for(int tsrc=0;tsrc<Lt;tsrc++){
	  WilsonMatrix expect = ftprop_srct[tsrc][tsnk];
	  expect.gl(-5);

	  if(!test_with_gauge_fixing){
	     expect *= Float(Vglob)/3.2; //normalise
	  }else{
	    WilsonMatrix wstar = FTw[tsnk];
	    wstar.hconj();
	    
	    expect.LeftTimesEqual(wstar);
	    expect *= 1.0/3.2;
	  }


	  for(int i=0;i<12;i++){
	    for(int j=0;j<12;j++){
	      const static int NA = -1;
	      Float* mf_val = (Float*)&mf(i,j,tsnk,NA,tsrc);
	      Float* e_val = (Float*)& expect(i/3,i%3, j/3,j%3);
	      //Check idx function
	      int lidx_e = i;
	      int ridx_e = j+12*tsrc;
	      int lidx_mf = mf.idx(i,NA,MesonField2::Left);
	      int ridx_mf = mf.idx(j,tsrc,MesonField2::Right);
	      if( lidx_e != lidx_mf || ridx_e != ridx_mf ){ 
		printf("Idx check fail: for W expect %d got %d,  for V expect %d got %d\n",lidx_e,lidx_mf,ridx_e,ridx_mf); fail = true;
	      }

	      //Check operator() is pulling out the right term
	      Float* chk = (Float*)mf.mf + 2*tsnk + 2*GJP.Tnodes()*GJP.TnodeSites()*( (j+12*tsrc) + 12*Lt*i );
	      if( chk!=mf_val ){
		printf("Pointer check fail: expect %p, got %p. Contents %.14e and %.14e compare to expected val from FTprop %.14e\n",chk,mf_val,*chk,*mf_val,*e_val);
		fail = true;
	      }
	      
	      for(int reim=0;reim<2;reim++){
		if( !ratio_diff_match(mf_val[reim],e_val[reim],1e-4) && fabs(mf_val[reim]-e_val[reim])>5e-7 ){ printf("Mesonfield test failed tsrc=%d tsnk=%d i=%d j=%d reim=%d: got %.14e expect %.14e, ratio %.14e\n",tsrc,tsnk,i,j,reim,mf_val[reim],e_val[reim],e_val[reim]/mf_val[reim] ); fail=true; }
	      }	    
	    }
	  }
	}
      }
      if(fail){ printf("Mesonfield test failed\n"); exit(-1); }
      else printf("Mesonfield test passed\n");
    }


    //Wall source wall sink
    CorrelationFunction corr_qpw("",1,CorrelationFunction::UNTHREADED);
    for(int t=0;t<Lt;t++){
      WilsonMatrix g1 = ftprop_srct[tsrc_con][t];
      g1.gr(-5);
      g1.gl(-5);

      WilsonMatrix g2 = ftprop_srct[t][tsrc_con];

      if(!test_with_gauge_fixing){
	g1 *= Float(Vglob)/3.2; //normalise
	g2 *= Float(Vglob)/3.2; //normalise
      }else{
	WilsonMatrix wstar = FTw[t];
	wstar.hconj();
	
	g1.LeftTimesEqual(wstar);
	g1 *= 1.0/3.2;
	
	wstar = FTw[tsrc_con];
	wstar.hconj();
	
	g2.LeftTimesEqual(wstar);
	g2 *= 1.0/3.2;
      }


      int t_dis = (t-tsrc_con+Lt)% Lt;

      //Each FTprop needs to be normalised to Float(Vglob)/3.2

      corr_qpw(0,t_dis) += Trace(g1,g2);
    }
    
    //Compare 
    bool fail=false;
    for(int t=0;t<GJP.Tnodes()*GJP.TnodeSites();++t){
      if( fabs( corr_a2a(0,t).real() - corr_qpw(0,t).real() ) > 1e-7 ){
	if(!UniqueID()) printf("Wall source comparison test real part fail t=%d:  %.14e %.14e, ratio qpw/a2a %.14e\n",t,corr_a2a(0,t).real(),corr_qpw(0,t).real(),corr_qpw(0,t).real()/corr_a2a(0,t).real());
	fail=true;
      }else if(!UniqueID()) printf("Wall source comparison test real part pass t=%d:  %.14e %.14e\n",t,corr_a2a(0,t).real(),corr_qpw(0,t).real());
      
      if( fabs( corr_a2a(0,t).imag() - corr_qpw(0,t).imag() ) > 1e-7 ){
	if(!UniqueID()) printf("Wall source comparison test imag part fail t=%d:  %.14e %.14e, ratio qpw/a2a %.14e\n",t,corr_a2a(0,t).imag(),corr_qpw(0,t).imag(),corr_qpw(0,t).imag()/corr_a2a(0,t).imag());
	fail=true;
      }else if(!UniqueID()) printf("Wall source comparison imag real part pass t=%d:  %.14e %.14e\n",t,corr_a2a(0,t).imag(),corr_qpw(0,t).imag());
    }
    if(fail){
      if(!UniqueID()) printf("Wall source comparison test failed, exiting\n");
      exit(-1);
    }else if(!UniqueID()) printf("Wall source comparison test passed\n");

    if(free_gfix) fix_gauge.free();
  }




  static void wallsource_amplitude_2f(Lattice &lattice, Lanczos_5d<double> &eig){
    //Use a MesonField2 object with flavor dilution and NORAND high modes for this comparison
    //Source type must be a box source filling the entire lattice 3-volume
    //Compare to traditional wall source contraction

    //Gauge fix lattice if not done already
    CommonArg c_arg;
    FixGaugeArg gfix_arg;
    setup_gfix_args(gfix_arg);
    AlgFixGauge fix_gauge(lattice,&c_arg,&gfix_arg);
    bool free_gfix(false);
    if(lattice.FixGaugeKind()==FIX_GAUGE_NONE){
      fix_gauge.run();
      free_gfix = true; //will be freed at end of method
    }
   
    //Generate the a2a prop
    int neig = 0;

    A2APropbfm* prop_f;
    A2APropbfmTesting::a2a_prop_gen(prop_f, &lattice, eig, 1,NORAND,neig);

    printf("Wall source prop with %d eigenvectors\n",prop_f->get_args().nl);

    //Calculate mesonfield with wall source, i.e. box source of size L^3
    MFqdpMatrix structure_wdagv(MFstructure::W, MFstructure::V, true, false,15,sigma1); //gamma^5 spin, sigma1 flavour
    MFBasicSource source(MFBasicSource::BoxSource,GJP.Xnodes()*GJP.XnodeSites());
    
    MesonField2 mf(*prop_f,*prop_f, structure_wdagv, source);

    //Calculate correlation function
    CorrelationFunction corr_a2a("",1, CorrelationFunction::THREADED);
    CorrelationFunction corr_a2a_tsepsum("",1, CorrelationFunction::THREADED);

    int tsrc=0;

    MesonField2::contract_specify_tsrc(mf,mf,0, tsrc, corr_a2a);
    MesonField2::contract(mf,mf,0, corr_a2a_tsepsum);

    //Generate a wall source propagator (yeah I know it should be a momentum source, but for this comparison with a single translationally non-invariant lattice it should not matter)
    PropManager::clear();

    int Lt = GJP.Tnodes()*GJP.TnodeSites();

    JobPropagatorArgs prop_args;
    SETUP_ARRAY(prop_args,props,PropagatorArg,Lt);
    prop_args.lanczos.lanczos_len = 0;
    std::string prop_names_srct[Lt];

    for(int t=0;t<Lt;t++){
      PropagatorArg &parg = prop_args.props.props_val[t];
    
      std::ostringstream os; os << "prop_"<<t;
      prop_names_srct[t] = os.str();

      parg.generics.tag = const_cast<char*>(prop_names_srct[t].c_str());
      parg.generics.type = QPROPW_TYPE;
      parg.generics.mass = 0.01; //needs to be the same as the mass in the a2a_arg function!
      parg.generics.bc[0] = GJP.Xbc();
      parg.generics.bc[1] = GJP.Ybc();
      parg.generics.bc[2] = GJP.Zbc();
      parg.generics.bc[3] = GJP.Tbc();

      SETUP_ARRAY(parg,attributes,AttributeContainer,3);
    
      ELEM(parg,attributes,0).type = WALL_SOURCE_ATTR;
      WallSourceAttrArg &srcarg = ELEM(parg,attributes,0).AttributeContainer_u.wall_source_attr;
      srcarg.t = t;

      ELEM(parg,attributes,1).type = GPARITY_FLAVOR_ATTR;
      GparityFlavorAttrArg &gparg = ELEM(parg,attributes,1).AttributeContainer_u.gparity_flavor_attr;
      gparg.flavor = 0; //other flavor will be automatically generated using the flavour relation

      ELEM(parg,attributes,2).type = CG_ATTR;
      CGAttrArg &cgattr = ELEM(parg,attributes,2).AttributeContainer_u.cg_attr;
      cgattr.max_num_iter = 5000;
      cgattr.stop_rsd = 1e-08;
      cgattr.true_rsd = 1e-08;
    }

    PropManager::setup(prop_args);
    PropManager::calcProps(lattice);

    for(int tsrc=0;tsrc<Lt;tsrc++){
      //Check propagators match. We have set up v and w such that v are just wall source propagators
      //There should be  Lt * 24 vectors. We want to pick out the one originating from timeslice tsrc
      //Mapping is
      //i-a2a.nl = (sc_id + 12*flav_id + 24/src_width * t_id +  24/src_width*Lt * wh_id)
      //where a2a.nl=0 here,  src_width=1  and nhits=1
      QPropWcontainer &pcon = PropManager::getProp(prop_names_srct[tsrc].c_str()).convert<QPropWcontainer>();

      bool fail = false;
      for(int x=0;x<GJP.VolNodeSites();++x){
	SpinColorFlavorMatrix site_matrix(pcon,lattice,x);

	for(int scf_id=0;scf_id<24;scf_id++){
	  int sc_id = scf_id % 12,  flav_id = scf_id/12;

	  int scol=sc_id/3;
	  int ccol=sc_id%3;

	  int i = sc_id + 12*flav_id + 24*tsrc;
	  Float* v_i = (Float*)prop_f->get_v(i) + 24*x;
	
	  //Here scf_id is the column spin-color index
	  for(int rowidx=0;rowidx<24;rowidx++){
	    int f_row = rowidx/12, sc_row = rowidx%12;

	    int srow = sc_row/3;
	    int crow = sc_row %3;
	    Float* qpw_cmp = (Float*)& site_matrix(srow,crow,f_row, scol,ccol,flav_id);
	    for(int reim=0;reim<2;reim++){
	      Float qpw_val = qpw_cmp[reim];
	      Float v_val = v_i[reim + 2*sc_row + f_row*24*GJP.VolNodeSites() ];
	    
	      //NOTE THERE IS A NORMALIZATION DIFFERENCE! QPW = (5-M5) A2A for DWF
	      qpw_val/= 3.2;

	      if( !ratio_diff_match(qpw_val,v_val,1e-3) ){
		printf("Wall source tsrc %d v and qpropw test fail x=%d sc_col=%d f_col=%d sc_row=%d f_row=%d reim=%d: %.14e %.14e, ratio %.14e\n",tsrc,x,sc_id,flav_id,sc_row,f_row,reim,qpw_val,v_val,qpw_val/v_val);
		fail=true;
	      }
	    }
	  }
	}
      }
      if(fail){
	printf("Wall source tsrc %d v and qpropw test failed\n",tsrc); exit(-1);
      }else printf("Wall source tsrc %d v and qpropw test passed\n",tsrc);    
    }



    std::vector<Float> mom(3,0.0);
    FourierProp<SpinColorFlavorMatrix> ftprop_calc;
    ftprop_calc.add_momentum(mom);

    std::vector<SpinColorFlavorMatrix> ftprop_srct[Lt];

    for(int t=0;t<Lt;t++)
      ftprop_srct[t] = ftprop_calc.getFTProp(lattice,mom,prop_names_srct[t].c_str());

    SpinColorFlavorMatrix FTw[Lt];
    for(int t=0;t<Lt;t++) FTw[t] = 0.0;

    WilsonMatrix _onew(0.0);
    for(int i=0;i<12;i++) _onew(i/3,i%3, i/3, i%3) = std::complex<Float>(1.0,0.0);

    SpinColorFlavorMatrix _one(0.0);
    for(int f=0;f<2;f++) _one(f,f) = _onew;

    //If gauge fixing we need the FT of w (zero mom) for comparison
    for(int t=0;t<GJP.TnodeSites();t++){
      int t_glb = GJP.TnodeCoor()*GJP.TnodeSites() + t;
      for(int x3d=0;x3d<GJP.VolNodeSites()/GJP.TnodeSites();x3d++){
	int x4d = x3d + GJP.VolNodeSites()/GJP.TnodeSites()*t;
	const Matrix &V_f0 =  *lattice.FixGaugeMatrix(x4d,0);
	const Matrix &V_f1 =  *lattice.FixGaugeMatrix(x4d,1);

	SpinColorFlavorMatrix tmp(_one);
	for(int j=0;j<2;j++){
	  tmp(0,j).LeftTimesEqual(V_f0);
	  tmp(1,j).LeftTimesEqual(V_f1);
	}
	FTw[t_glb] += tmp;
      }
    }
    for(int t=0;t<Lt;t++)
      for(int f1=0;f1<2;f1++)
	for(int f2=0;f2<2;f2++)
	  slice_sum( (Float*)& FTw[t](f1,f2), 12*12*2, 99);

    {
      //Check FFT w field is  w_a^i = \sum_x e^{ipx} V_ab(x,t) \delta_b^i 

      bool fail=false;
      for(int i=0;i<24;i++){
	int i_f = i/12, i_sc=i%12;

	Float* wh_fftw = (Float*)prop_f->get_wh_fftw() + 24*2*GJP.VolNodeSites()*( 12*i_f + i_sc );
	for(int t=0;t<GJP.TnodeSites();t++){
	  for(int p=0;p<GJP.VolNodeSites()/GJP.TnodeSites();p++){
	    int t_glb = GJP.TnodeCoor()*GJP.TnodeSites() + t;
	    
	    int p_glb[3]; int rem=p;
	    for(int i=0;i<3;i++){ 
	      p_glb[i] = rem % GJP.NodeSites(i) + GJP.NodeCoor(i)*GJP.NodeSites(i); rem/=GJP.NodeSites(i); 
	    }

	    for(int scf=0;scf<24;scf++){
	      int sc = scf%12, f= scf/12;

	      for(int reim=0;reim<2;reim++){
		Float* val = wh_fftw + reim + 2*sc + 24*( p+ GJP.VolNodeSites()/GJP.TnodeSites()*t + f*GJP.VolNodeSites() );

		//Just do p=0 part
		if( p_glb[0]==0 && p_glb[1]==0 && p_glb[2]==0){
		  Float expect = ((Float*)&(FTw[t_glb](f,i_f)(sc/3,sc%3, i_sc/3,i_sc%3)))[reim];
		    
		  if( fabs(*val - expect)>1e-7 ){ printf("w fail at p=0 on t=%d, i=%d sc=%d f=%d: expect %.14e got %.14e, ratio %.14e\n",t_glb,i,sc,f,expect,*val, expect/(*val)); fail=true; }
		}
	      }
	    }

	  }
	}
      }
      if(fail){ printf("w fail\n"); exit(-1); }
      else printf("w pass\n");
    }
    CorrelationFunction corr_qpw("",1,CorrelationFunction::UNTHREADED);

    for(int t=0;t<Lt;t++){
      SpinColorFlavorMatrix g1 = ftprop_srct[tsrc][t];
      g1.gl(-5).pl(sigma1);

      SpinColorFlavorMatrix g2 = ftprop_srct[t][tsrc];
      g2.gl(-5).pl(sigma1);

      SpinColorFlavorMatrix wstar = FTw[t];
      wstar.hconj();
	
      g1.LeftTimesEqual(wstar);
      g1 *= 1.0/3.2;
	
      wstar = FTw[tsrc];
      wstar.hconj();
	
      g2.LeftTimesEqual(wstar);
      g2 *= 1.0/3.2;

      int t_dis = (t-tsrc+Lt)% Lt;

      //Each FTprop needs to be normalised to Float(Vglob)/3.2

      corr_qpw(0,t_dis) += Trace(g1,g2);
    }

    //Compare 
    bool fail(false);
    for(int t=0;t<GJP.Tnodes()*GJP.TnodeSites();++t){
      if( !ratio_diff_match(corr_a2a(0,t).real(),corr_qpw(0,t).real(),1e-3) && fabs( corr_a2a(0,t).real() - corr_qpw(0,t).real() ) > 1e-7 ){
	if(!UniqueID()) printf("Wall source comparison test real part fail t=%d:  %.14e %.14e, ratio qpw/a2a %.14e\n",t,corr_a2a(0,t).real(),corr_qpw(0,t).real(), corr_qpw(0,t).real()/corr_a2a(0,t).real() );
	fail=true;
      }else if(!UniqueID()) printf("Wall source comparison test real part pass t=%d:  %.14e %.14e\n",t,corr_a2a(0,t).real(),corr_qpw(0,t).real());
      
      if( !ratio_diff_match(corr_a2a(0,t).imag(),corr_qpw(0,t).imag(),1e-3) && fabs( corr_a2a(0,t).imag() - corr_qpw(0,t).imag() ) > 1e-7 ){
	if(!UniqueID()) printf("Wall source comparison test imag part fail t=%d:  %.14e %.14e, ratio qpw/a2a %.14e\n",t,corr_a2a(0,t).imag(),corr_qpw(0,t).imag(), corr_qpw(0,t).imag()/corr_a2a(0,t).imag());
	fail=true;
      }else if(!UniqueID()) printf("Wall source comparison imag real part pass t=%d:  %.14e %.14e\n",t,corr_a2a(0,t).imag(),corr_qpw(0,t).imag());
    }
    if(fail){
      if(!UniqueID()) printf("Wall source comparison test failed, exiting\n");
      exit(-1);
    }else if(!UniqueID()) printf("Wall source comparison test passed, exiting\n");


    //Test using alg_gparitycontract methods
    PropManager::clear();

    //Put the existing Lanczos object in PropManager for tag retrieval
    {
      LanczosContainerArg larg;
      larg.tag = "lanczos";
      larg.lanc_arg.fname = "";

      PropManager::addLanczos(larg).set_lanczos(&eig);
    }


    {
      JobPropagatorArgs prop_args_a2a;
      SETUP_ARRAY(prop_args_a2a,props,PropagatorArg,1);
      
      PropagatorArg &parg = prop_args_a2a.props.props_val[0];
    
      parg.generics.tag = "a2a_prop";
      parg.generics.type = A2A_PROP_TYPE;
      parg.generics.mass = 0.01; //needs to be the same as the mass in the a2a_arg function!
      parg.generics.bc[0] = GJP.Xbc();
      parg.generics.bc[1] = GJP.Ybc();
      parg.generics.bc[2] = GJP.Zbc();
      parg.generics.bc[3] = GJP.Tbc();

      SETUP_ARRAY(parg,attributes,AttributeContainer,1);
    
      ELEM(parg,attributes,0).type = A2A_ATTR;
      A2AAttrArg &a2a_arg = ELEM(parg,attributes,0).AttributeContainer_u.a2a_attr;

      a2a_arg  = prop_f->get_args();
      a2a_arg.lanczos_tag = "lanczos";

      PropManager::setup(prop_args_a2a);   
      PropManager::calcProps(lattice);
    }

    {
      ContractionTypeA2ABilinear con_args;
      con_args.prop_src_snk = "a2a_prop";
      con_args.prop_snk_src = "a2a_prop";
      con_args.source_smearing.type = BOX_3D_SMEARING;
      con_args.source_smearing.A2ASmearing_u.box_3d_smearing.side_length = GJP.Xnodes()*GJP.XnodeSites();
      con_args.sink_smearing.type = BOX_3D_SMEARING;
      con_args.sink_smearing.A2ASmearing_u.box_3d_smearing.side_length = GJP.Xnodes()*GJP.XnodeSites();
      SETUP_ARRAY(con_args,source_spin_matrix,MatIdxAndCoeff,1);
      SETUP_ARRAY(con_args,sink_spin_matrix,MatIdxAndCoeff,1);
      SETUP_ARRAY(con_args,source_flavor_matrix,MatIdxAndCoeff,1);
      SETUP_ARRAY(con_args,sink_flavor_matrix,MatIdxAndCoeff,1);

      con_args.source_spin_matrix.source_spin_matrix_val[0].idx = 15;
      con_args.source_spin_matrix.source_spin_matrix_val[0].coeff = 1.0;

      con_args.sink_spin_matrix.sink_spin_matrix_val[0].idx = 15;
      con_args.sink_spin_matrix.sink_spin_matrix_val[0].coeff = 1.0;

      con_args.source_flavor_matrix.source_flavor_matrix_val[0].idx = 1;
      con_args.source_flavor_matrix.source_flavor_matrix_val[0].coeff = 1.0;

      con_args.sink_flavor_matrix.sink_flavor_matrix_val[0].idx = 1;
      con_args.sink_flavor_matrix.sink_flavor_matrix_val[0].coeff = 1.0;

      con_args.file = "";
    
      CommonArg c_arg;
      AlgGparityContract con(lattice,c_arg);
      CorrelationFunction gpcon_result("",1, CorrelationFunction::THREADED);
      con.contract_a2a_bilinear(con_args, gpcon_result);

      for(int t=0;t<GJP.Tnodes()*GJP.TnodeSites();++t){
	if( fabs( corr_a2a_tsepsum(0,t).real() - gpcon_result(0,t).real() ) > 1e-12 ){
	  if(!UniqueID()) printf("AlgGparityContract comparison test real part fail t=%d:  %.14e %.14e, ratio gpcon/a2a %.14e\n",t,corr_a2a_tsepsum(0,t).real(),gpcon_result(0,t).real(), gpcon_result(0,t).real()/corr_a2a_tsepsum(0,t).real() );
	  fail=true;
	}else if(!UniqueID()) printf("AlgGparityContract comparison test real part pass t=%d:  %.14e %.14e\n",t,corr_a2a_tsepsum(0,t).real(),gpcon_result(0,t).real());
      
	if( fabs( corr_a2a_tsepsum(0,t).imag() - gpcon_result(0,t).imag() ) > 1e-12 ){
	  if(!UniqueID()) printf("AlgGparityContract comparison test imag part fail t=%d:  %.14e %.14e, ratio gpcon/a2a %.14e\n",t,corr_a2a_tsepsum(0,t).imag(),gpcon_result(0,t).imag(), gpcon_result(0,t).imag()/corr_a2a_tsepsum(0,t).imag());
	  fail=true;
	}else if(!UniqueID()) printf("AlgGparityContract comparison imag real part pass t=%d:  %.14e %.14e\n",t,corr_a2a_tsepsum(0,t).imag(),gpcon_result(0,t).imag());
      }
      if(fail){	if(!UniqueID()) printf("AlgGparityContract comparison test failed, exiting\n"); exit(-1); }
      else if(!UniqueID()) printf("AlgGparityContract comparison test passed, exiting\n");    
    }


    PropManager::getLanczos("lanczos").set_lanczos(NULL); //prevent future PropManager::clear() operations from destroying the Lanczos object

    if(free_gfix) fix_gauge.free();
  }
  

};
CPS_END_NAMESPACE

//(complex<double> *)mf + (j*nvec[0]+i)*t_size*n_flav2 + t_size*(2*g + f) +glb_t;
  // int mf_fac = GJP.Gparity() || GJP.Gparity1fX() ? 4 : 1; //flavour matrix indices also

  // mf = (Vector *)smalloc(cname, fname, "mf", sizeof(Float)*nvec[0]*(nl[0]+sc_size*nhits[0])*t_size*2 * mf_fac);

using namespace Chroma;
using namespace cps;

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

int cout_time(char *);
void ReadGaugeField(const MeasArg &meas_arg);
void bfm_init(bfm_evo<double> &dwf,double mq);

int main (int argc,char **argv )
{
  Start(&argc, &argv);

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
    printf("Doing G-parity test in X direction\n");
  }else if(arg0==1){
    printf("Doing G-parity test in X and Y directions\n");
    gparity_X = true;
    gparity_Y = true;
  }else{
    printf("Doing standard lattice test\n");
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
  bool verbose(false);
  bool skip_gparity_inversion(false);
  bool unit_gauge(false);

  int size[] = {2,2,2,2,2};

  int dilute_flavor = 0;

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
    }else if( strncmp(cmd,"-save_lrg",15) == 0){
      if(i==argc-1){ printf("-save_lrg requires an argument\n"); exit(-1); }
      save_lrg=true;
      save_lrg_file = argv[i+1];
      i+=2;  
    }else if( strncmp(cmd,"-verbose",15) == 0){
      verbose=true;
      i++;
    }else if( strncmp(cmd,"-dilute_flavor",15) == 0){
      printf("Running with flavor dilution\n");
      dilute_flavor = 1;
      i++;
    }else if( strncmp(cmd,"-skip_gparity_inversion",30) == 0){
      skip_gparity_inversion=true;
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

  if(gparity_X) do_arg.x_bc = BND_CND_GPARITY;
  if(gparity_Y) do_arg.y_bc = BND_CND_GPARITY;

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

  cps_qdp_init(&argc,&argv);
  
#if TARGET == BGQ
  omp_set_num_threads(64);
#else
  omp_set_num_threads(1);
#endif

  bool precon = true;
  //Backup LRG
  LatRanGen LRGbak = LRG;

  //Generate eigenvectors  
  Lanczos_5d<double>* eig_2f;
  create_eig(lattice,eig_2f,precon);

  //Different tests if no G-parity active
  if(!gparity_X && !gparity_Y){
    MesonFieldTesting::wallsource_amplitude_nogparity(*lattice, *eig_2f);
#ifdef HAVE_BFM
    Chroma::finalize();
#endif
    
    if(UniqueID()==0) printf("Main job complete\n"); 
    return 0;
  }

  //Generate A2A prop in 2f environment
  A2APropbfm* prop_2f;
  A2APropbfmTesting::a2a_prop_gen(prop_2f, lattice, *eig_2f, dilute_flavor);

  //Pull the eigenvectors out of eig into CPS canonical ordered fields for later 2f->1f conversion
  Float** eigv_2f_cps;
  eig_convert_cps(*eig_2f, eigv_2f_cps,precon);

  //Generate gauge fixing matrices
  CommonArg c_arg;
  FixGaugeArg gfix_arg;
  setup_gfix_args(gfix_arg);
  AlgFixGauge fix_gauge_2f(*lattice,&c_arg,&gfix_arg);
  
  //Generate MesonField object
  MesonField mesonfield_2f(*lattice, prop_2f, &fix_gauge_2f, &c_arg);
  mesonfield_2f.allocate_vw_fftw();
  
  //Do the FFTW on the A2APropbfm also and test this new code
  fix_gauge_2f.run();
  prop_2f->allocate_vw_fftw();
  prop_2f->fft_vw();
  fix_gauge_2f.free(); //free gauge fixing matrices

  mesonfield_2f.prepare_vw(); //re-performs gauge fixing
  if(!dilute_flavor) MesonFieldTesting::compare_fftw_fields(mesonfield_2f,*prop_2f); //I did not modify Daiqian's original code for flavor dilution

  //Calculate mesonfield with exponential source with radius 2
  int source_type = 2;   //exp 1 box 2
  double radius = 2;

  MFBasicSource::SourceType source_type2 = MFBasicSource::BoxSource; //MFBasicSource::ExponentialSource;
  if(!dilute_flavor) mesonfield_2f.cal_mf_ll(radius,source_type);

  //Try to duplicate mf_ll using MesonField2
  MFqdpMatrix structure_2f(MFstructure::W, MFstructure::V, true, false,15,sigma0); //gamma^5 spin, unit mat flavour
  MFBasicSource source_2f(source_type2,radius);
  
  if(!dilute_flavor) MesonFieldTesting::compare_source_MesonField2_2f(mesonfield_2f,source_2f);  
  MesonField2 mf2_2f(*prop_2f,*prop_2f, structure_2f, source_2f);
  if(!dilute_flavor) MesonFieldTesting::compare_mf_ll_MesonField2_2f(mesonfield_2f, mf2_2f);

  //Calculate a 'kaon' correlation function in both Daiqian's code and mine (both modified for G-parity)
  if(!dilute_flavor){
    AlgFixGauge fix_gauge_2f_ls(*lattice,&c_arg,&gfix_arg);
    MesonField mesonfield_2f_ls(*lattice, prop_2f, prop_2f, &fix_gauge_2f_ls, &c_arg); //pretend the propagator is also a strange quark for testing :)
    mesonfield_2f_ls.allocate_vw_fftw();
    mesonfield_2f_ls.prepare_vw();

    mesonfield_2f_ls.cal_mf_sl(radius,source_type);
    MesonFieldTesting::compare_mf_sl_MesonField2_2f(mesonfield_2f_ls,mf2_2f);

    mesonfield_2f_ls.cal_mf_ls(radius,source_type);
    MesonFieldTesting::compare_mf_ls_MesonField2_2f(mesonfield_2f_ls,mf2_2f);

    std::complex<double> *kaoncorr_orig = (std::complex<double> *)pmalloc(sizeof(std::complex<double>)*GJP.Tnodes()*GJP.TnodeSites());
    mesonfield_2f_ls.run_kaon(kaoncorr_orig);
    
    CorrelationFunction kaoncorr_mf2("",1);
    MesonField2::contract(mf2_2f, mf2_2f, kaoncorr_mf2);
    
    MesonFieldTesting::compare_kaon_corr_MesonField2_2f(kaoncorr_orig, kaoncorr_mf2);

    //Check threaded version
    CorrelationFunction kaoncorr_mf2_thr("",1, CorrelationFunction::THREADED);
    MesonField2::contract(mf2_2f, mf2_2f, kaoncorr_mf2_thr);
    kaoncorr_mf2_thr.sumThreads();

    MesonFieldTesting::compare_kaon_corr_MesonField2_2f(kaoncorr_orig, kaoncorr_mf2_thr);

  }
  
  //Test flavor dilution and everything else by generating a2a props with non-random high modes (1.0 on all sites) and with a wall source
  //and compare to traditional wall source
  if(dilute_flavor) MesonFieldTesting::wallsource_amplitude_2f(*lattice, *eig_2f);

  //Restore LRG backup to reset RNG for 1f section
  LRG = LRGbak;

  //Move to 1f environment
  if(UniqueID()==0) printf("Starting double lattice section\n");
  
  int array_size = 2*lattice->GsiteSize() * GJP.VolNodeSites() * sizeof(Float);
  cps::Matrix *orig_lattice = (cps::Matrix *) pmalloc(array_size);
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
    QDP::Layout::setLattSize(nrow);
    QDP::Layout::create();
  }
#endif
  lattice = new GwilsonFdwf;
  setup_double_latt(*lattice,orig_lattice,gparity_X,gparity_Y);
  setup_double_rng(gparity_X,gparity_Y);
   
  //Convert eigenvectors from 2f to 1f to use here (assumes 2f and 1f lanczos agree - this just ensures they agree perfectly, such that any numerical disparities are localised to the A2A code)
  Lanczos_5d<double>* eig_1f;
  create_eig_1f(lattice, eig_1f, *eig_2f, eigv_2f_cps, precon, gparity_X, gparity_Y);

  //Create 1f A2A prop
  A2APropbfm* prop_1f;
  A2APropbfmTesting::a2a_prop_gen(prop_1f, lattice, *eig_1f, dilute_flavor);

  //In-place convert 2f A2A prop to 1f format
  A2APropbfmTesting::convert_2f_A2A_prop_to_1f(*prop_2f, gparity_X, gparity_Y);

  //Compare props
  A2APropbfmTesting::compare_1f_2f_A2A_prop(*prop_2f, *prop_1f);

  //Fix the gauge
  AlgFixGauge fix_gauge_1f(*lattice,&c_arg,&gfix_arg);
  
  if(!dilute_flavor){
    //Generate MesonField object
    MesonField mesonfield_1f(*lattice, prop_1f, &fix_gauge_1f, &c_arg);
    mesonfield_1f.allocate_vw_fftw();
    mesonfield_1f.prepare_vw();
    
    //Calculate mesonfield with exponential source with radius 2
    mesonfield_1f.cal_mf_ll(radius,source_type);
    
    MesonFieldTesting::convert_mesonfield_2f_1f(mesonfield_2f,gparity_X,gparity_Y);
    
    //Compare meson fields
    MesonFieldTesting::compare_fftw_vecs(mesonfield_1f,mesonfield_2f);
    MesonFieldTesting::compare_mf_ll(mesonfield_1f, mesonfield_2f);
  }
  //Test MesonField2 1f code against its 2f version (its 2f version has already been compared to the original MesonField code above)
  fix_gauge_1f.run();
  prop_1f->allocate_vw_fftw();
  prop_1f->fft_vw();
  fix_gauge_1f.free();

  MFqdpMatrix structure_1f(MFstructure::W, MFstructure::V, true, false,15,sigma0); //gamma^5 spin, unit mat flavour
  MFBasicSource source_1f(source_type2,radius);
  
  MesonField2 mf2_1f(*prop_1f,*prop_1f, structure_1f, source_1f);
  MesonFieldTesting::compare_mf_ll_MesonField2_2f_1f(mf2_2f,mf2_1f);



#ifdef HAVE_BFM
  Chroma::finalize();
#endif

  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  
  return 0;
}
