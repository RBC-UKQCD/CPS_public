//Periodically calculate true residual of single precision multi-mass solve using same parameters as the 32^3 gparity run, only using the 16^3 mobius+DSDR configs on a 128 node machine
//Local volume will not be identical because 32^3 job has 8^4 local volume whereas 16^3 has 8x4x4x8
//However we might be able to determine the optimum single precision residual

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

#include <util/lat_cont.h>

#include <chroma.h>
#include <omp.h>
#include <pthread.h>
#include <sys/time.h>
#include <sstream>

#ifdef HAVE_BFM
#include <chroma.h>
#endif



//-------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

const char *cname = "";

DoArg do_arg;
RationalQuotientRemezArg rq_arg;

#define decode_vml(arg_name) \
  if(!UniqueID()) printf("Decoding %s.vml\n",#arg_name); fflush(stdout); \
  do{									\
        if ( ! arg_name.Decode(#arg_name".vml", #arg_name) )            \
            ERR.General(cname, fname, "Bad " #arg_name ".vml.\n");      \
    } while(0)

void decode_vml_all(void)
{
    char *fname = "decode_vml_all()";

    //decode_vml(b_arg);
    decode_vml(do_arg);
    decode_vml(rq_arg);
}

void setup_bfmargs(bfmarg &dwfa, const Float &residual, const int &nthreads){
  if(!UniqueID()) printf("Setting up bfmargs\n");


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
    if(!UniqueID()) printf("G-parity directions: ");
    for(int d=0;d<3;d++)
      if(GJP.Bc(d) == BND_CND_GPARITY){ dwfa.gparity_dir[d] = 1; if(!UniqueID()) printf("%d ",d); }
      else dwfa.gparity_dir[d] = 0;
    for(int d=0;d<4;d++){
      dwfa.nodes[d] = procs[d];
      dwfa.ncoor[d] = ncoor[d];
    }
    if(!UniqueID()) printf("\n");
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
      if(!UniqueID()) printf("Non-local comms in direction %d\n",mu);
    } else { 
      dwfa.local_comm[mu] = 1;
      if(!UniqueID()) printf("Local comms in direction %d\n",mu);
    }
  }

  dwfa.precon_5d = 0; //mobius uses 4d preconditioning
  dwfa.mobius_scale = 32.0/12.0;

  dwfa.Ls   = GJP.SnodeSites();
  dwfa.solver = HmCayleyTanh;
  dwfa.M5   = toDouble(GJP.DwfHeight());
  dwfa.mass = toDouble(0.0001);
  dwfa.Csw  = 0.0;
  dwfa.max_iter = 100000;
  dwfa.residual = residual;
  if(!UniqueID()) printf("Finished setting up bfmargs\n");
}



static int multi_shift_trueresid(Fermion_t psi[], 
				 Fermion_t src,
				 double    mass[],
				 double    alpha[],
				 int       nshift,
				 double mresidual[],
				 int single,
				 bfm_evo<float> &bfm_f,
				 bfm_evo<double> &bfm_d,
				 int report_freq = 100)
{
  //NOTE: Assumes bfm_d comms are active

  int me = bfm_d.thread_barrier();
 
  // Per shift fields
  Fermion_t ps [nshift]; // search directions 
  double    bs [nshift];
  double    rsq[nshift];
  double    z[nshift][2];
  int       converged[nshift];

  const int       primary =0;

  //Primary shift fields CG iteration
  double a,b,c,d;
  double cp,bp; //prev
  Fermion_t p = bfm_f.threadedAllocFermion(mem_slow);
  Fermion_t r = bfm_f.threadedAllocFermion(mem_slow);

  mixed_cg::switch_comm(bfm_f, bfm_d);

  //Single precision copy of src
  Fermion_t src_f = bfm_f.threadedAllocFermion(mem_slow);
  mixed_cg::threaded_convFermion(src_f, src, bfm_f, bfm_d);

  // Matrix mult fields in single and double prec
  Fermion_t tmp = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t tmp2 = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t mp  = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t mmp = bfm_f.threadedAllocFermion(mem_fast); 

  Fermion_t tmp_d = bfm_d.threadedAllocFermion(mem_slow); 
  Fermion_t mp_d  = bfm_d.threadedAllocFermion(mem_slow); 
  Fermion_t mmp_d = bfm_d.threadedAllocFermion(mem_slow); 
  Fermion_t r_d = bfm_d.threadedAllocFermion(mem_slow); 
  Fermion_t p_d = bfm_d.threadedAllocFermion(mem_slow); 

  // Check lightest mass
  for(int s=0;s<nshift;s++){
    if ( mass[s] < mass[primary] ) {
      printf("First shift not lightest - oops\n");
      exit(-1);
    }
  }

  Fermion_t psi_f [nshift];

  for(int i=0;i<nshift;i++){
    ps[i] = bfm_f.threadedAllocFermion(mem_slow);
    converged[i]=0;

    //init single prec solution
    psi_f[i] = bfm_f.threadedAllocFermion(mem_slow);
    mixed_cg::threaded_convFermion(psi_f[i], psi[i], bfm_f, bfm_d);
  }

  // Wire guess to zero
  // Residuals "r" are src
  // First search direction "p" is also src
  cp = bfm_f.norm(src_f);
  for(int s=0;s<nshift;s++){
    rsq[s] = cp * mresidual[s] * mresidual[s];
    bfm_f.copy(ps[s],src_f);
  }
  // r and p for primary
  bfm_f.copy(r,src_f);
  bfm_f.copy(p,src_f);
  

  d= bfm_f.Mprec(p,mp,tmp,DaggerNo,1); //mp = Mpc p,  what is the 'norm' of?? I think its |Mpc p|^2
  bfm_f.Mprec(mp,mmp,tmp,DaggerYes); //mmp = Mpc^dag mp = Mpc^dag Mpc p
  bfm_f.axpy(mmp,p,mmp,mass[0]); //mmp = p*mass[0]+mmp
  
  //Overall  mmp = (Mpc^dag Mpc + mass[0])*p
  double rn = bfm_f.norm(p);

  d += rn*mass[0];
  b = -cp /d;
  if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec_multi: b = -cp/d = -%le/%le = %le\n",cp,d,b);

  // Set up the various shift variables
  int       iz=0;

  z[0][1-iz] = 1.0;
  z[0][iz]   = 1.0;
  bs[0]      = b;
  for(int s=1;s<nshift;s++){
    z[s][1-iz] = 1.0;
    z[s][iz]   = 1.0/( 1.0 - b*(mass[s]-mass[0]));
    bs[s]      = b*z[s][iz]; // Sign relative to Mike - FIXME
  }

  // r += b[0] A.p[0]
  // c= norm(r)
  c= bfm_f.axpy_norm(r,mmp,r,b);

  if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec_multi: k=0 residual %le \n",c);

  for(int s=0;s<nshift;s++) {
    bfm_f.axpby(psi_f[s],src_f,src_f,0.,-bs[s]*alpha[s]);
  }
  

  // Iteration loop
  for (int k=1;k<=bfm_f.max_iter;k++){
    a = c /cp;

    // if((k%100)==0){
    //   mixed_cg::threaded_convFermion(r_d,r, bfm_d, bfm_f);
    //   mixed_cg::threaded_convFermion(p_d,p, bfm_d, bfm_f);
    // }

    bfm_f.axpy(p,p,r,a);

    for(int s=0;s<nshift;s++){
      if ( ! converged[s] ) {
        if (s==0){
	  //DEBUG
	  double prenorm;
	  static const int printfreq = 1;
	  
	  if(k%printfreq==0) prenorm = bfm_f.norm(ps[0]);
	  //DEBUG

          bfm_f.axpy(ps[s],ps[s],r,a);

	  //DEBUG
	  if(k%printfreq==0){
	    double modr = bfm_f.norm(r);
	    double modps = bfm_f.norm(ps[0]);
	    if ( bfm_d.isBoss() && (!me) ) printf("shift[0] : Search direction with orig norm2 %.12e multiplied by %.12e and shifted in direction of r vector with norm2 %.12e to give resulting norm2 %.12e\n",prenorm,a,modr,modps);
	  }
	  //DEBUG

        } else{
	  double as =a *z[s][iz]*bs[s] /(z[s][1-iz]*b);
	  bfm_f.axpby(ps[s],r,ps[s],z[s][iz],as);
        }
      }
    }

    cp=c;
    
    d= bfm_f.Mprec(p,mp,tmp,DaggerNo,1); 
    bfm_f.Mprec(mp,mmp,tmp,DaggerYes);
    bfm_f.axpy(mmp,p,mmp,mass[0]);
    double rn = bfm_f.norm(p);

    d += rn*mass[0];
    bp=b;
    b=-cp/d;
    
    c= bfm_f.axpy_norm(r,mmp,r,b);

    // if((k%100)==0){
    //   //Repeat above steps in double precision and see if there is a difference
    //   mixed_cg::switch_comm(bfm_d, bfm_f);
    //   bfm_d.axpy(p_d,p_d,r_d,a);
      
    //   double d_d = bfm_d.Mprec(p_d,mp_d,tmp_d,DaggerNo,1); 
    //   bfm_d.Mprec(mp_d,mmp_d,tmp_d,DaggerYes);
    //   bfm_d.axpy(mmp_d,p_d,mmp_d,mass[0]);
    //   double rn_d = bfm_d.norm(p_d);

    //   d_d += rn_d*mass[0];
    //   double b_d=-cp/d_d;
    
    //   double c_d= bfm_d.axpy_norm(r_d,mmp_d,r_d,b_d);

    //   if ( bfm_d.isBoss() && (!me) )printf("single prec c=%g,  double prec c=%g\n",c,c_d);

    //   mixed_cg::switch_comm(bfm_f, bfm_d);
    // }


    // Toggle the recurrence history
    bs[0] = b;
    iz = 1-iz;

    for(int s=1;s<nshift;s++){
      if(!converged[s]){
	double z0 = z[s][1-iz];
	double z1 = z[s][iz];
	z[s][iz] = z0*z1*bp
	  / (b*a*(z1-z0) + z1*bp*(1- (mass[s]-mass[0])*b)); 
	bs[s] = b*z[s][iz]/z0; // NB sign  rel to Mike
      }
    }
    
    for(int s=0;s<nshift;s++){
      int ss = s;
      if(!converged[s])
	bfm_f.axpy(psi_f[ss],ps[s],psi_f[ss],-bs[s]*alpha[s]);
    }

    // Convergence checks
    int all_converged = 1;
    if ( ((k%100)==0) && bfm_f.isBoss() && (!me) ) printf("bfmbase::CGNE_prec_multi: k=%d c=%g, shift in current dir for lightest pole %.12e\n",k,c,-bs[0]*alpha[0]);
    for(int s=0;s<nshift;s++){
      if (!converged[s]){
	double css  = c * z[s][iz]* z[s][iz];	
	if(css<rsq[s]) converged[s]=1;
	else all_converged=0;
	if(bfm_f.isBoss() && (!me) && converged[s]) printf("bfmbase::CGNE_prec_multi: Shift %d converged on iter %d: test cur %g, targ %g   [Stated true resid %g].\n",s,k,css,rsq[s],(css/rsq[s])*mresidual[s]);
	else if ( ((k%100)==0) && bfm_f.isBoss() && (!me) ) printf("bfmbase::CGNE_prec_multi: Shift %d convergence test cur %g, targ %g   [Stated true resid %g].\n",s,css,rsq[s],sqrt(css/rsq[s])*mresidual[s]);	
      }
    }

    if ( all_converged ){
      if ( bfm_f.isBoss() && (!me) )printf("bfmbase::CGNE_prec_multi: k=%d All shifts have converged\n",k);
      if ( bfm_f.isBoss() && (!me) )printf("bfmbase::CGNE_prec_multi: k=%d Checking solutions\n",k);
      // Check answers 

      for(int s=0; s < nshift; s++) {
	bfm_f.Mprec(psi_f[s],mp,tmp,DaggerNo);
	bfm_f.Mprec(mp,mmp,tmp,DaggerYes);
	bfm_f.axpy(tmp,psi_f[s],mmp,mass[s]);
	bfm_f.axpy(tmp2,tmp,src_f,-1);
	double rn = bfm_f.norm(tmp2);
	double cn = bfm_f.norm(src_f);
	if ( bfm_f.isBoss() && !me ) {
	  printf("single prec final: shift[%d] true residual %le \n",s,sqrt(rn/cn));
	}
      } 

      mixed_cg::switch_comm(bfm_d, bfm_f);

      for(int s=0; s < nshift; s++) {
	//Convert solution to double precision
	mixed_cg::threaded_convFermion(psi[s], psi_f[s], bfm_d, bfm_f);
	
	bfm_d.Mprec(psi[s],mp_d,tmp_d,DaggerNo);
	bfm_d.Mprec(mp_d,mmp_d,tmp_d,DaggerYes);
        bfm_d.axpy(tmp_d,psi[s],mmp_d,mass[s]);
	bfm_d.axpy(r_d,tmp_d,src,-1);
	double rn = bfm_d.norm(r_d);
	double cn = bfm_d.norm(src);
	if ( bfm_d.isBoss() && !me ) {
	  printf("double prec final: shift[%d] true residual %.12le \n",s,sqrt(rn/cn));
	}
      }

      if ( single ) {
	for(int s=1; s < nshift; s++) { 
	  bfm_d.axpy(psi[0],psi[s],psi[0],1.0);
	}      
      }

      bfm_f.threadedFreeFermion(src_f);
      bfm_f.threadedFreeFermion(tmp);
      bfm_f.threadedFreeFermion(tmp2);
      bfm_f.threadedFreeFermion(p);
      bfm_f.threadedFreeFermion(mp);
      bfm_f.threadedFreeFermion(mmp);
      bfm_f.threadedFreeFermion(r);
      for(int i=0;i<nshift;i++){
	bfm_f.threadedFreeFermion(ps[i]);
	bfm_f.threadedFreeFermion(psi_f[i]);
      }  

      bfm_d.threadedFreeFermion(tmp_d);
      bfm_d.threadedFreeFermion(mp_d);
      bfm_d.threadedFreeFermion(mmp_d);
      bfm_d.threadedFreeFermion(r_d);
      bfm_d.threadedFreeFermion(p_d);

      return k;

    }else if(k % report_freq == 0){
      if ( bfm_f.isBoss() && !me ) printf("Report: iter %d\n",k);

      int nreport = nshift; //1;//nshift;

      for(int s=0; s < nreport; s++) {
	double css  = c * z[s][iz]* z[s][iz];
	if ( bfm_f.isBoss() && !me ) printf("running: shift[%d] true residual %.12e\n",s,sqrt(css/rsq[s])*mresidual[s]);
      }

      //Periodically report double precision true residual as well as single prec residual
      //Single
      for(int s=0; s < nreport; s++) {
	bfm_f.Mprec(psi_f[s],mp,tmp,DaggerNo);
	bfm_f.Mprec(mp,mmp,tmp,DaggerYes);
	bfm_f.axpy(tmp,psi_f[s],mmp,mass[s]);
	bfm_f.axpy(tmp2,tmp,src_f,-1);
	double rn = bfm_f.norm(tmp2);
	double cn = bfm_f.norm(src_f);
	if ( bfm_f.isBoss() && !me ) {
	  printf("single prec: shift[%d] true residual %.12le [converged = %d]\n",s,sqrt(rn/cn),converged[s]);
	}
      }      

      //Double
      mixed_cg::switch_comm(bfm_d, bfm_f);
      for(int s=0; s < nreport; s++) {
	//Convert solution to double precision
	mixed_cg::threaded_convFermion(psi[s], psi_f[s], bfm_d, bfm_f);
	
	bfm_d.Mprec(psi[s],mp_d,tmp_d,DaggerNo);
	bfm_d.Mprec(mp_d,mmp_d,tmp_d,DaggerYes);
        bfm_d.axpy(tmp_d,psi[s],mmp_d,mass[s]);
	bfm_d.axpy(r_d,tmp_d,src,-1);
	double rn = bfm_d.norm(r_d);
	double cn = bfm_d.norm(src);
	if ( bfm_d.isBoss() && !me ) {
	  printf("double prec: shift[%d] true residual %.12le [converged = %d]\n",s,sqrt(rn/cn),converged[s]);
	}
      }
      mixed_cg::switch_comm(bfm_f, bfm_d);
    }
  }

  if ( bfm_d.isBoss() && !me ) printf("bfmbase::CGNE_prec_multi: CG not converged \n");

  bfm_f.threadedFreeFermion(src_f);
  bfm_f.threadedFreeFermion(tmp);
  bfm_f.threadedFreeFermion(tmp2);
  bfm_f.threadedFreeFermion(p);
  bfm_f.threadedFreeFermion(mp);
  bfm_f.threadedFreeFermion(mmp);
  bfm_f.threadedFreeFermion(r);
  for(int i=0;i<nshift;i++){
    bfm_f.threadedFreeFermion(ps[i]);
    bfm_f.threadedFreeFermion(psi_f[i]);
  }  

  bfm_d.threadedFreeFermion(tmp_d);
  bfm_d.threadedFreeFermion(mp_d);
  bfm_d.threadedFreeFermion(mmp_d);
  bfm_d.threadedFreeFermion(r_d);
  bfm_d.threadedFreeFermion(p_d);

  return -1;
}

void setup(int *argc, char ***argv)
{
    const char *fname = "setup()";

    Start(argc, argv);

    if(*argc < 2) {
        ERR.General(cname, fname, "Must provide VML directory.\n");
    }

    if(chdir( (*argv)[1]) != 0) {
      ERR.General(cname, fname, "Changing directory to %s failed.\n", (*argv)[1]);
    }

    decode_vml_all();
 
    VRB.Result(cname, fname, "Read VML files successfully.\n");

#if TARGET == BGQ
  LRG.setSerial();
#endif

    GJP.Initialize(do_arg);
    LRG.Initialize();

    VRB.Result(cname, fname, "VRB.Level(%d)\n", do_arg.verbose_level);
    VRB.Level(do_arg.verbose_level);

    Chroma::initialize(argc, argv);
    multi1d<int> nrow(Nd);

    for(int i = 0; i< Nd; ++i)
        nrow[i] = GJP.Sites(i);

    Layout::setLattSize(nrow);
    Layout::create();
}




static int multi_shift_dp_spmprec (Fermion_t psi[], 
				   Fermion_t src,
				   double    mass[],
				   double    alpha[],
				   int       nshift,
				   double mresidual[],
				   int single,
				   bfm_evo<float> &bfm_f,
				   bfm_evo<double> &bfm_d)
{
  //Assumes bfm_d comms active!
  int me = bfm_d.thread_barrier();
   
  // Per shift fields
  Fermion_t ps [nshift]; // search directions 
  double    bs [nshift];
  double    rsq[nshift];
  double    z[nshift][2];
  int       converged[nshift];

  const int       primary =0;

  //Primary shift fields CG iteration
  double a,b,c,d;
  double cp,bp; //prev
  Fermion_t p ;
  Fermion_t r ;
  Fermion_t p_f ;

  // Matrix mult fields
  Fermion_t tmp = bfm_d.threadedAllocFermion(mem_fast); 
  Fermion_t mp  = bfm_d.threadedAllocFermion(mem_fast); 
  Fermion_t mmp = bfm_d.threadedAllocFermion(mem_fast); 

  Fermion_t tmp_f = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t mp_f  = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t mmp_f = bfm_f.threadedAllocFermion(mem_fast); 

  // Check lightest mass
  for(int s=0;s<nshift;s++){
    if ( mass[s] < mass[primary] ) {
      printf("First shift not lightest - oops\n");
      exit(-1);
    }
  }

  //Allocate per shift fields
  r = bfm_d.threadedAllocFermion(mem_slow);
  p = bfm_d.threadedAllocFermion(mem_slow);
  p_f = bfm_f.threadedAllocFermion(mem_slow);

  for(int i=0;i<nshift;i++){
    ps[i] = bfm_d.threadedAllocFermion(mem_slow);
    converged[i]=0;
  }

  // Wire guess to zero
  // Residuals "r" are src
  // First search direction "p" is also src
  cp = bfm_d.norm(src);
  for(int s=0;s<nshift;s++){
    rsq[s] = cp * mresidual[s] * mresidual[s];
    if ( bfm_d.isBoss() && (!me) && bfm_d.verbose ) printf("bfmbase::CGNE_prec_multi: shift %d, shift amount %le, target resid %le, src norm %le, targ true resid %le\n",s,mass[s],rsq[s],cp,mresidual[s]);
    bfm_d.copy(ps[s],src);
  }
  // r and p for primary
  bfm_d.copy(r,src);
  bfm_d.copy(p,src);
  
  mixed_cg::switch_comm(bfm_f, bfm_d);
  mixed_cg::threaded_convFermion(p_f, p, bfm_f, bfm_d);
  
  d= bfm_f.Mprec(p_f,mp_f,tmp_f,DaggerNo,1); 
  bfm_f.Mprec(mp_f,mmp_f,tmp_f,DaggerYes); 
  bfm_f.axpy(mmp_f,p_f,mmp_f,mass[0]); 

  mixed_cg::threaded_convFermion(mmp, mmp_f, bfm_d, bfm_f);
  mixed_cg::switch_comm(bfm_d, bfm_f);

  double rn = bfm_d.norm(p);
  d += rn*mass[0];


  b = -cp /d;
  if ( bfm_d.isBoss() && !me ) printf("bfmbase::CGNE_prec_multi: b = -cp/d = -%le/%le = %le\n",cp,d,b);

  // Set up the various shift variables
  int       iz=0;

  z[0][1-iz] = 1.0;
  z[0][iz]   = 1.0;
  bs[0]      = b;
  for(int s=1;s<nshift;s++){
    z[s][1-iz] = 1.0;
    z[s][iz]   = 1.0/( 1.0 - b*(mass[s]-mass[0]));
    bs[s]      = b*z[s][iz]; // Sign relative to Mike - FIXME
  }

  // r += b[0] A.p[0]
  // c= norm(r)
  c=bfm_d.axpy_norm(r,mmp,r,b);

  if ( bfm_d.isBoss() && !me ) printf("bfmbase::CGNE_prec_multi: k=0 residual %le \n",c);

  // Linear algebra overhead can be minimised if summing results
  // psi-= b[0]p[0] = -b[0]src

  for(int s=0;s<nshift;s++) {
    bfm_d.axpby(psi[s],src,src,0.,-bs[s]*alpha[s]);
  }
  

  // Iteration loop
  for (int k=1;k<= bfm_d.max_iter;k++){
    //if ( isBoss() && !me ) printf("bfmbase::CGNE_prec_multi: k=%d residual %le \n",k,c);

    a = c /cp;
    bfm_d.axpy(p,p,r,a);

    // Note to self - direction ps is iterated seperately
    // for each shift. Does not appear to have any scope
    // for avoiding linear algebra in "single" case.
    // 
    // However SAME r is used. Could load "r" and update
    // ALL ps[s]. Bandwidth goes to that of "copy".
    // Could compare AXPY and COPY performance to estimate gain.
    for(int s=0;s<nshift;s++){
      if ( ! converged[s] ) {
        if (s==0){
          bfm_d.axpy(ps[s],ps[s],r,a);
        } else{
	  double as =a *z[s][iz]*bs[s] /(z[s][1-iz]*b);
	  bfm_d.axpby(ps[s],r,ps[s],z[s][iz],as);
        }
      }
    }

    cp=c;

    mixed_cg::switch_comm(bfm_f, bfm_d);
    mixed_cg::threaded_convFermion(p_f, p, bfm_f, bfm_d);

    d= bfm_f.Mprec(p_f,mp_f,tmp_f,DaggerNo,1); 
    bfm_f.Mprec(mp_f,mmp_f,tmp_f,DaggerYes);
    bfm_f.axpy(mmp_f,p_f,mmp_f,mass[0]);

    mixed_cg::threaded_convFermion(mmp, mmp_f, bfm_d, bfm_f);
    mixed_cg::switch_comm(bfm_d, bfm_f);

    double rn = bfm_d.norm(p);
    d += rn*mass[0];

    bp=b;
    b=-cp/d;
    
    c=bfm_d.axpy_norm(r,mmp,r,b);

    // Toggle the recurrence history
    bs[0] = b;
    iz = 1-iz;

    for(int s=1;s<nshift;s++){
      if(!converged[s]){
	double z0 = z[s][1-iz];
	double z1 = z[s][iz];
	z[s][iz] = z0*z1*bp
	  / (b*a*(z1-z0) + z1*bp*(1- (mass[s]-mass[0])*b)); 
	bs[s] = b*z[s][iz]/z0; // NB sign  rel to Mike
      }
    }
    
    // Fixme - Mike has a variable "at" that is missing here.
    // Understand this. NB sign of bs rel to Mike
    for(int s=0;s<nshift;s++){
      int ss = s;
      // Scope for optimisation here in case of "single".
      // Could load psi[0] and pull all ps[s] in.
      //      if ( single ) ss=primary;
      if(!converged[s])
	bfm_d.axpy(psi[ss],ps[s],psi[ss],-bs[s]*alpha[s]);
    }

    // Convergence checks
    int all_converged = 1;
    if ( ((k%100)==0) && bfm_d.isBoss() && (!me) ) printf("bfmbase::CGNE_prec_multi: k=%d c=%g\n",k,c);
    for(int s=0;s<nshift;s++){

      if (!converged[s]){

	double css  = c * z[s][iz]* z[s][iz];
	
	if ( ((k%100)==0) && bfm_d.isBoss() && (!me) ) printf("bfmbase::CGNE_prec_multi: Shift %d convergence test cur %g, targ %g   [True resid %g]\n",s,css,rsq[s],sqrt(css/rsq[s])*mresidual[s]);

	if(css<rsq[s]){
	  converged[s]=1;
	  //if ( isBoss() && (!me) ) printf("bfmbase::CGNE_prec_multi: Shift %d convergence test cur %g, targ %g\n",s,css,rsq[s]);
	  if ( bfm_d.isBoss() && (!me) ) printf("bfmbase::CGNE_prec_multi: k=%d Shift %d has converged\n",k,s);
	} else {
	  all_converged=0;
	}

      }
    }

    if ( all_converged ){
      if ( bfm_d.isBoss() && (!me) )printf("bfmbase::CGNE_prec_multi: k=%d All shifts have converged\n",k);
      if ( bfm_d.isBoss() && (!me) )printf("bfmbase::CGNE_prec_multi: k=%d Checking solutions\n",k);
      // Check answers 
      for(int s=0; s < nshift; s++) { 
	bfm_d.Mprec(psi[s],mp,tmp,DaggerNo);
	bfm_d.Mprec(mp,mmp,tmp,DaggerYes);
        bfm_d.axpy(tmp,psi[s],mmp,mass[s]);
	bfm_d.axpy(r,tmp,src,-1);
	double rn = bfm_d.norm(r);
	double cn = bfm_d.norm(src);
	if ( bfm_d.isBoss() && !me ) {
	  printf("bfmbase::CGNE_prec_multi: shift[%d] true residual %le \n",s,sqrt(rn/cn));
	}
      }

      if ( single ) {
	for(int s=1; s < nshift; s++) { 
	  bfm_d.axpy(psi[0],psi[s],psi[0],1.0);
	}      
      }

      bfm_d.threadedFreeFermion(tmp);
      bfm_d.threadedFreeFermion(p);
      bfm_d.threadedFreeFermion(mp);
      bfm_d.threadedFreeFermion(mmp);
      bfm_d.threadedFreeFermion(r);
      for(int i=0;i<nshift;i++){
	bfm_d.threadedFreeFermion(ps[i]);
      }  
      bfm_f.threadedFreeFermion(tmp_f);
      bfm_f.threadedFreeFermion(p_f);
      bfm_f.threadedFreeFermion(mp_f);
      bfm_f.threadedFreeFermion(mmp_f);
      return k;
    }
  }

  if ( bfm_d.isBoss() && !me ) printf("bfmbase::CGNE_prec_multi: CG not converged \n");

  bfm_d.threadedFreeFermion(tmp);
  bfm_d.threadedFreeFermion(p);
  bfm_d.threadedFreeFermion(mp);
  bfm_d.threadedFreeFermion(mmp);
  bfm_d.threadedFreeFermion(r);
  for(int i=0;i<nshift;i++){
    bfm_d.threadedFreeFermion(ps[i]);
  }  
  bfm_f.threadedFreeFermion(tmp_f);
  bfm_f.threadedFreeFermion(p_f);
  bfm_f.threadedFreeFermion(mp_f);
  bfm_f.threadedFreeFermion(mmp_f);

  return -1;
}





static int multi_shift_sp_dp_switch (Fermion_t psi[], 
				     Fermion_t src,
				     double    mass[],
				     double    alpha[],
				     int       nshift,
				     double mresidual[],
				     int single,
				     bfm_evo<float> &bfm_f,
				     bfm_evo<double> &bfm_d,
				     int switch_iter,
				     int report_freq = -1)
{
  //Assumes bfm_d comms active!
  int me = bfm_d.thread_barrier();
   
  // Per shift fields
  Fermion_t ps [nshift], ps_f [nshift]; // search directions 
  double    bs [nshift];
  double    rsq[nshift];
  double    z[nshift][2];
  int       converged[nshift];

  const int       primary =0;

  //Primary shift fields CG iteration
  double a,b,c,d;
  double cp,bp; //prev
  Fermion_t p, p_f ;
  Fermion_t r, r_f ;

  // Matrix mult fields
  Fermion_t tmp = bfm_d.threadedAllocFermion(mem_fast); 
  Fermion_t mp  = bfm_d.threadedAllocFermion(mem_fast); 
  Fermion_t mmp = bfm_d.threadedAllocFermion(mem_fast); 

  Fermion_t tmp_f = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t mp_f  = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t mmp_f = bfm_f.threadedAllocFermion(mem_fast); 

  // Check lightest mass
  for(int s=0;s<nshift;s++){
    if ( mass[s] < mass[primary] ) {
      printf("First shift not lightest - oops\n");
      exit(-1);
    }
  }

  //Allocate per shift fields
  r = bfm_d.threadedAllocFermion(mem_slow);
  p = bfm_d.threadedAllocFermion(mem_slow);

  r_f = bfm_f.threadedAllocFermion(mem_slow);
  p_f = bfm_f.threadedAllocFermion(mem_slow);

  Fermion_t psi_f [nshift]; //single prec solutions

  for(int i=0;i<nshift;i++){
    ps[i] = bfm_d.threadedAllocFermion(mem_slow);
    ps_f[i] = bfm_f.threadedAllocFermion(mem_slow);
    psi_f[i] = bfm_f.threadedAllocFermion(mem_slow);
    converged[i]=0;
  }

  Fermion_t src_f = bfm_f.threadedAllocFermion(mem_slow); 
  mixed_cg::threaded_convFermion(src_f, src, bfm_f, bfm_d);
  
  bool dp = (switch_iter == 0 ? true : false); //switch: double prec (true) single (false)

  //Initially have single prec active
  if(!dp) mixed_cg::switch_comm(bfm_f, bfm_d);

  cp = (dp ? bfm_d.norm(src) : bfm_f.norm(src_f));
  for(int s=0;s<nshift;s++){
    rsq[s] = cp * mresidual[s] * mresidual[s];
    if ( bfm_f.isBoss() && (!me) && bfm_d.verbose ) printf("bfmbase::CGNE_prec_multi: shift %d, shift amount %le, target resid %le, src norm %le, targ true resid %le\n",s,mass[s],rsq[s],cp,mresidual[s]);
    dp ? bfm_d.copy(ps[s],src) : bfm_f.copy(ps_f[s],src_f);
  }
  // r and p for primary
  double rn;

  if(dp){
    bfm_d.copy(r,src);
    bfm_d.copy(p,src);
    
    d= bfm_d.Mprec(p,mp,tmp,DaggerNo,1); 
    bfm_d.Mprec(mp,mmp,tmp,DaggerYes); 
    bfm_d.axpy(mmp,p,mmp,mass[0]);

    rn = bfm_d.norm(p); 
  }else{
    bfm_f.copy(r_f,src_f);
    bfm_f.copy(p_f,src_f);
    
    d= bfm_f.Mprec(p_f,mp_f,tmp_f,DaggerNo,1); 
    bfm_f.Mprec(mp_f,mmp_f,tmp_f,DaggerYes); 
    bfm_f.axpy(mmp_f,p_f,mmp_f,mass[0]); 

    rn = bfm_f.norm(p_f);
  }

  d += rn*mass[0];

  b = -cp /d;
  if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec_multi: b = -cp/d = -%le/%le = %le\n",cp,d,b);

  // Set up the various shift variables
  int       iz=0;

  z[0][1-iz] = 1.0;
  z[0][iz]   = 1.0;
  bs[0]      = b;
  for(int s=1;s<nshift;s++){
    z[s][1-iz] = 1.0;
    z[s][iz]   = 1.0/( 1.0 - b*(mass[s]-mass[0]));
    bs[s]      = b*z[s][iz]; // Sign relative to Mike - FIXME
  }

  // r += b[0] A.p[0]
  c= (dp ? bfm_d.axpy_norm(r,mmp,r,b) : bfm_f.axpy_norm(r_f,mmp_f,r_f,b));

  if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec_multi: k=0 residual %le \n",c);

  for(int s=0;s<nshift;s++) dp ? bfm_d.axpby(psi[s],src,src,0.,-bs[s]*alpha[s]) : bfm_f.axpby(psi_f[s],src_f,src_f,0.,-bs[s]*alpha[s]);
  
  // Iteration loop
  for (int k=1;k<= bfm_d.max_iter;k++){
    if( k == switch_iter ){
      if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec_multi: Switching to double precision on iteration %d\n",k);
      mixed_cg::threaded_convFermion(r, r_f, bfm_d, bfm_f);
      mixed_cg::threaded_convFermion(p, p_f, bfm_d, bfm_f);
      for(int s=0;s<nshift;s++){
	mixed_cg::threaded_convFermion(ps[s], ps_f[s], bfm_d, bfm_f);
      	mixed_cg::threaded_convFermion(psi[s], psi_f[s], bfm_d, bfm_f);
      }
      mixed_cg::switch_comm(bfm_d, bfm_f);
      dp = true;
    }

    a = c /cp;
    dp ? bfm_d.axpy(p,p,r,a) : bfm_f.axpy(p_f,p_f,r_f,a);

    for(int s=0;s<nshift;s++){
      if ( ! converged[s] ) {
        if (s==0){
          dp ? bfm_d.axpy(ps[s],ps[s],r,a) : bfm_f.axpy(ps_f[s],ps_f[s],r_f,a);
        } else{
	  double as =a *z[s][iz]*bs[s] /(z[s][1-iz]*b);
	  dp ? bfm_d.axpby(ps[s],r,ps[s],z[s][iz],as) : bfm_f.axpby(ps_f[s],r_f,ps_f[s],z[s][iz],as);
        }
      }
    }

    cp=c;

    double rn;
    
    if(dp){
      d= bfm_d.Mprec(p,mp,tmp,DaggerNo,1); 
      bfm_d.Mprec(mp,mmp,tmp,DaggerYes);
      bfm_d.axpy(mmp,p,mmp,mass[0]);
      rn = bfm_d.norm(p);
    }else{
      d= bfm_f.Mprec(p_f,mp_f,tmp_f,DaggerNo,1); 
      bfm_f.Mprec(mp_f,mmp_f,tmp_f,DaggerYes);
      bfm_f.axpy(mmp_f,p_f,mmp_f,mass[0]);
      rn = bfm_f.norm(p_f);
    }

    d += rn*mass[0];

    bp=b;
    b=-cp/d;
    
    c = (dp ? bfm_d.axpy_norm(r,mmp,r,b) : bfm_f.axpy_norm(r_f,mmp_f,r_f,b));

    // Toggle the recurrence history
    bs[0] = b;
    iz = 1-iz;

    for(int s=1;s<nshift;s++){
      if(!converged[s]){
	double z0 = z[s][1-iz];
	double z1 = z[s][iz];
	z[s][iz] = z0*z1*bp
	  / (b*a*(z1-z0) + z1*bp*(1- (mass[s]-mass[0])*b)); 
	bs[s] = b*z[s][iz]/z0; // NB sign  rel to Mike
      }
    }
    
    // Fixme - Mike has a variable "at" that is missing here.
    // Understand this. NB sign of bs rel to Mike
    for(int s=0;s<nshift;s++){
      int ss = s;
      if(!converged[s])
	dp ? bfm_d.axpy(psi[ss],ps[s],psi[ss],-bs[s]*alpha[s]) : bfm_f.axpy(psi_f[ss],ps_f[s],psi_f[ss],-bs[s]*alpha[s]);
    }

    // Convergence checks
    int all_converged = 1;
    if ( ((k%100)==0) && bfm_d.isBoss() && (!me) ) printf("bfmbase::CGNE_prec_multi: k=%d c=%g\n",k,c);
    for(int s=0;s<nshift;s++){

      if (!converged[s]){

	double css  = c * z[s][iz]* z[s][iz];
	
	if ( ((k%100)==0) && bfm_d.isBoss() && (!me) ) printf("bfmbase::CGNE_prec_multi: Shift %d convergence test cur %g, targ %g   [True resid %g]\n",s,css,rsq[s],sqrt(css/rsq[s])*mresidual[s]);

	if(css<rsq[s]){
	  converged[s]=1;
	  if ( bfm_d.isBoss() && (!me) ) printf("bfmbase::CGNE_prec_multi: k=%d Shift %d has converged\n",k,s);
	} else {
	  all_converged=0;
	}
      }
    }

    if ( all_converged ){
      if ( bfm_d.isBoss() && (!me) )printf("bfmbase::CGNE_prec_multi: k=%d All shifts have converged\n",k);
      if ( bfm_d.isBoss() && (!me) )printf("bfmbase::CGNE_prec_multi: k=%d Checking solutions\n",k);
      
      if(!dp){ //Didn't reach the switchover point before converging
	for(int s=0;s<nshift;s++) mixed_cg::threaded_convFermion(psi[s], psi_f[s], bfm_d, bfm_f);
	mixed_cg::switch_comm(bfm_d, bfm_f);
      }

      // Check answers 
      for(int s=0; s < nshift; s++) { 
	bfm_d.Mprec(psi[s],mp,tmp,DaggerNo);
	bfm_d.Mprec(mp,mmp,tmp,DaggerYes);
        bfm_d.axpy(tmp,psi[s],mmp,mass[s]);
	bfm_d.axpy(r,tmp,src,-1);
	double rn = bfm_d.norm(r);
	double cn = bfm_d.norm(src);
	if ( bfm_d.isBoss() && !me ) {
	  printf("bfmbase::CGNE_prec_multi: shift[%d] true residual %le \n",s,sqrt(rn/cn));
	}
      }

      if ( single ) {
	for(int s=1; s < nshift; s++) { 
	  bfm_d.axpy(psi[0],psi[s],psi[0],1.0);
	}      
      }

      bfm_d.threadedFreeFermion(tmp);
      bfm_d.threadedFreeFermion(p);
      bfm_d.threadedFreeFermion(mp);
      bfm_d.threadedFreeFermion(mmp);
      bfm_d.threadedFreeFermion(r);
      for(int i=0;i<nshift;i++){
	bfm_d.threadedFreeFermion(ps[i]);
	bfm_f.threadedFreeFermion(ps_f[i]);
	bfm_f.threadedFreeFermion(psi_f[i]);
      }  
      bfm_f.threadedFreeFermion(tmp_f);
      bfm_f.threadedFreeFermion(p_f);
      bfm_f.threadedFreeFermion(mp_f);
      bfm_f.threadedFreeFermion(mmp_f);
      return k;

    }else if(report_freq != -1 && k % report_freq == 0){
      if ( bfm_f.isBoss() && !me ) printf("Report: iter %d\n",k);

      int nreport = nshift; //1;//nshift;

      if(!dp){
	mixed_cg::switch_comm(bfm_d, bfm_f);
	for(int s=0; s < nreport; s++) mixed_cg::threaded_convFermion(psi[s], psi_f[s], bfm_d, bfm_f);
      }
      for(int s=0; s < nreport; s++) {
	bfm_d.Mprec(psi[s],mp,tmp,DaggerNo);
	bfm_d.Mprec(mp,mmp,tmp,DaggerYes);
        bfm_d.axpy(tmp,psi[s],mmp,mass[s]);
	bfm_d.axpy(mp,tmp,src,-1);
	double rn = bfm_d.norm(mp);
	double cn = bfm_d.norm(src);
	double css  = c * z[s][iz]* z[s][iz];
	if ( bfm_d.isBoss() && !me ) {
	  printf("Report: shift[%d] true residual %.12le, running resid %.12le\n",s,sqrt(rn/cn),sqrt(css/rsq[s])*mresidual[s]);
	}
      }
      if(!dp) mixed_cg::switch_comm(bfm_f, bfm_d);
    }
  }

  if ( bfm_d.isBoss() && !me ) printf("bfmbase::CGNE_prec_multi: CG not converged \n");

  bfm_d.threadedFreeFermion(tmp);
  bfm_d.threadedFreeFermion(p);
  bfm_d.threadedFreeFermion(mp);
  bfm_d.threadedFreeFermion(mmp);
  bfm_d.threadedFreeFermion(r);
  for(int i=0;i<nshift;i++){
    bfm_d.threadedFreeFermion(ps[i]);
    bfm_f.threadedFreeFermion(ps_f[i]);
    bfm_f.threadedFreeFermion(psi_f[i]);
  }  
  bfm_f.threadedFreeFermion(tmp_f);
  bfm_f.threadedFreeFermion(p_f);
  bfm_f.threadedFreeFermion(r_f);
  bfm_f.threadedFreeFermion(mp_f);
  bfm_f.threadedFreeFermion(mmp_f);

  return -1;
}


static int multi_shift_sp_reliable_update_dp_vects(Fermion_t psi[], 
						   Fermion_t src,
						   double    mass[],
						   double    alpha[],
						   int       nshift,
						   double mresidual[],
						   int single,
						   bfm_evo<float> &bfm_f,
						   bfm_evo<double> &bfm_d,
						   int update_freq = 100,
						   int report_freq = -1)
{
  //NOTE: Assumes bfm_d comms are active
  //update_freq is the frequency at which the reliable update step is performed
  //report_freq prints the double precision true residual when k % report_freq = 0. Use -1 to disable

  int me = bfm_d.thread_barrier();
 
  double    bs [nshift];
  double    rsq[nshift];
  double    z[nshift][2];
  int       converged[nshift];

  const int       primary =0;

  //Primary shift fields CG iteration
  double a,b,c,d;
  double cp,bp; //prev

  //Single precision fields
  Fermion_t r = bfm_f.threadedAllocFermion(mem_slow); //residual vector, single precision  
  Fermion_t tmp = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t p = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t mp  = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t mmp = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t src_f = bfm_f.threadedAllocFermion(mem_slow);
  mixed_cg::threaded_convFermion(src_f, src, bfm_f, bfm_d);

  //Double precision fields
  Fermion_t p_d = bfm_d.threadedAllocFermion(mem_fast); //search direction, double precision
  Fermion_t tmp_d = bfm_d.threadedAllocFermion(mem_fast); 
  Fermion_t mp_d  = bfm_d.threadedAllocFermion(mem_fast); 
  Fermion_t mmp_d = bfm_d.threadedAllocFermion(mem_fast); 
  Fermion_t ps_d [nshift]; // search directions (double precision)

  for(int i=0;i<nshift;i++){
    ps_d[i] = bfm_d.threadedAllocFermion(mem_slow);
    converged[i]=0;
  }

#define DEALLOCATE_ALL \
  bfm_f.threadedFreeFermion(r); \
  bfm_f.threadedFreeFermion(tmp); \
  bfm_f.threadedFreeFermion(p); \
  bfm_f.threadedFreeFermion(mp); \
  bfm_f.threadedFreeFermion(mmp); \
  bfm_f.threadedFreeFermion(src_f); \
  bfm_d.threadedFreeFermion(p_d);   \
  bfm_d.threadedFreeFermion(tmp_d); \
  bfm_d.threadedFreeFermion(mp_d); \
  bfm_d.threadedFreeFermion(mmp_d); \
  for(int s=0;s<nshift;s++) bfm_d.threadedFreeFermion(ps_d[s])

  // Check lightest mass
  for(int s=0;s<nshift;s++){
    if ( mass[s] < mass[primary] ) {
      printf("First shift not lightest - oops\n");
      exit(-1);
    }
  }

  cp = bfm_d.norm(src);
  for(int s=0;s<nshift;s++){
    rsq[s] = cp * mresidual[s] * mresidual[s];
    bfm_d.copy(ps_d[s],src);
  }
  // r and p for primary
  bfm_f.copy(r,src_f); //residual vector in single prec
  bfm_d.copy(p_d,src);

  double rn = cp; //norm of src = p_d
  
  mixed_cg::switch_comm(bfm_f, bfm_d);
  mixed_cg::threaded_convFermion(p, p_d, bfm_f, bfm_d);
  d= bfm_f.Mprec(p,mp,tmp,DaggerNo,1); //mp = Mpc p,  what is the 'norm' of?? I think its |Mpc p|^2
  bfm_f.Mprec(mp,mmp,tmp,DaggerYes); //mmp = Mpc^dag mp = Mpc^dag Mpc p
  bfm_f.axpy(mmp,p,mmp,mass[0]); //mmp = p*mass[0]+mmp
  
  d += rn*mass[0];
  b = -cp /d;
  if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec_multi: b = -cp/d = -%le/%le = %le\n",cp,d,b);

  // Set up the various shift variables
  int       iz=0;

  z[0][1-iz] = 1.0;
  z[0][iz]   = 1.0;
  bs[0]      = b;
  for(int s=1;s<nshift;s++){
    z[s][1-iz] = 1.0;
    z[s][iz]   = 1.0/( 1.0 - b*(mass[s]-mass[0]));
    bs[s]      = b*z[s][iz]; // Sign relative to Mike - FIXME
  }

  c= bfm_f.axpy_norm(r,mmp,r,b);
  if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec_multi: k=0 residual %le \n",c);

  for(int s=0;s<nshift;s++) {
    bfm_d.axpby(psi[s],src,src,0.,-bs[s]*alpha[s]); //initialize double prec solutions
  }
  
  // Iteration loop
  for (int k=1;k<=bfm_f.max_iter;k++){
    a = c /cp;

    mixed_cg::threaded_convFermion(tmp_d, r, bfm_d, bfm_f); //store double prec copy of r in tmp_d
    bfm_d.axpy(p_d,p_d,tmp_d,a);

    for(int s=0;s<nshift;s++){
      if ( ! converged[s] ) {
        if (s==0){
          bfm_d.axpy(ps_d[s],ps_d[s],tmp_d,a);
        } else{
	  double as =a *z[s][iz]*bs[s] /(z[s][1-iz]*b);
	  bfm_d.axpby(ps_d[s],tmp_d,ps_d[s],z[s][iz],as);
        }
      }
    }

    cp=c;
    
    mixed_cg::threaded_convFermion(p, p_d, bfm_f, bfm_d);
    d= bfm_f.Mprec(p,mp,tmp,DaggerNo,1); 
    bfm_f.Mprec(mp,mmp,tmp,DaggerYes);
    bfm_f.axpy(mmp,p,mmp,mass[0]);
    double rn = bfm_f.norm(p);

    d += rn*mass[0];
    bp=b;
    b=-cp/d;
    
    // Toggle the recurrence history
    bs[0] = b;
    iz = 1-iz;

    for(int s=1;s<nshift;s++){
      if(!converged[s]){
	double z0 = z[s][1-iz];
	double z1 = z[s][iz];
	z[s][iz] = z0*z1*bp
	  / (b*a*(z1-z0) + z1*bp*(1- (mass[s]-mass[0])*b)); 
	bs[s] = b*z[s][iz]/z0; // NB sign  rel to Mike
      }
    }
    
    for(int s=0;s<nshift;s++){
      int ss = s;
      if(!converged[s])
    	bfm_d.axpy(psi[ss],ps_d[s],psi[ss],-bs[s]*alpha[s]);
    }

    //Reliable update
    if(k % update_freq == 0){
      double c_sp = bfm_f.axpy_norm(r,mmp,r,b);

      //Replace r with true residual
      mixed_cg::switch_comm(bfm_d, bfm_f);
      bfm_d.Mprec(psi[0],mp_d,tmp_d,0,1);
      bfm_d.Mprec(mp_d,mmp_d,tmp_d,1);
      bfm_d.axpy(mmp_d,psi[0],mmp_d,mass[0]);

      c = bfm_d.axpy_norm(tmp_d,mmp_d,src,-1.0);

      if( bfm_d.isBoss() && !me) printf("bfmbase::CGNE_prec_multi: reliable update iter %d, replaced |r|^2 = %.12le with |r|^2 = %.12le\n",k,c_sp,c);

      mixed_cg::threaded_convFermion(r, tmp_d, bfm_f, bfm_d);
      mixed_cg::switch_comm(bfm_f, bfm_d);
    }else{
      c= bfm_f.axpy_norm(r,mmp,r,b);
    }

    // Convergence checks
    int all_converged = 1;
    if ( ((k%100)==0) && bfm_f.isBoss() && (!me) ) printf("bfmbase::CGNE_prec_multi: k=%d c=%g, shift in current dir for lightest pole %.12e\n",k,c,-bs[0]*alpha[0]);
    for(int s=0;s<nshift;s++){
      if (!converged[s]){
	double css  = c * z[s][iz]* z[s][iz];	
	if(css<rsq[s]) converged[s]=1;
	else all_converged=0;
	if(bfm_f.isBoss() && (!me) && converged[s]) printf("bfmbase::CGNE_prec_multi: Shift %d converged on iter %d: test cur %g, targ %g   [Stated true resid %g].\n",s,k,css,rsq[s],(css/rsq[s])*mresidual[s]);
	else if ( ((k%100)==0) && bfm_f.isBoss() && (!me) ) printf("bfmbase::CGNE_prec_multi: Shift %d convergence test cur %g, targ %g   [Stated true resid %g].\n",s,css,rsq[s],sqrt(css/rsq[s])*mresidual[s]);	
      }
    }

    if ( all_converged ){
      if ( bfm_f.isBoss() && (!me) )printf("bfmbase::CGNE_prec_multi: k=%d All shifts have converged\n",k);
      if ( bfm_f.isBoss() && (!me) )printf("bfmbase::CGNE_prec_multi: k=%d Checking solutions\n",k);

      // Check answers
      mixed_cg::switch_comm(bfm_d, bfm_f);

      for(int s=0; s < nshift; s++) {
	//Convert solution to double precision
	bfm_d.Mprec(psi[s],mp_d,tmp_d,DaggerNo);
	bfm_d.Mprec(mp_d,mmp_d,tmp_d,DaggerYes);
        bfm_d.axpy(tmp_d,psi[s],mmp_d,mass[s]);
	bfm_d.axpy(mp_d,tmp_d,src,-1);
	double rn = bfm_d.norm(mp_d);
	double cn = bfm_d.norm(src);
	if ( bfm_d.isBoss() && !me ) {
	  printf("double prec final: shift[%d] true residual %.12le \n",s,sqrt(rn/cn));
	}
      }

      if ( single ) {
	for(int s=1; s < nshift; s++) { 
	  bfm_d.axpy(psi[0],psi[s],psi[0],1.0);
	}      
      }

      DEALLOCATE_ALL;

      return k;

    }else if(report_freq != -1 && k % report_freq == 0){
      mixed_cg::switch_comm(bfm_d, bfm_f);
      for(int s=0; s < nshift; s++) {
	double css  = c * z[s][iz]* z[s][iz];
	bfm_d.Mprec(psi[s],mp_d,tmp_d,DaggerNo);
	bfm_d.Mprec(mp_d,mmp_d,tmp_d,DaggerYes);
        bfm_d.axpy(tmp_d,psi[s],mmp_d,mass[s]);
	bfm_d.axpy(mp_d,tmp_d,src,-1);
	double rn = bfm_d.norm(mp_d);
	double cn = bfm_d.norm(src);
	if ( bfm_d.isBoss() && !me ) {
	  printf("iter %d, double prec: shift[%d] true residual %.12le, running true residual %.12le [converged = %d]\n",k,s,sqrt(rn/cn),sqrt(css/rsq[s])*mresidual[s],converged[s]);
	}
      }
      mixed_cg::switch_comm(bfm_f, bfm_d);
    }
  }

  if ( bfm_d.isBoss() && !me ) printf("bfmbase::CGNE_prec_multi: CG not converged \n");

  DEALLOCATE_ALL;

  return -1;
}
#undef DEALLOCATE_ALL



static int CGNE_singleprec_trueresid(Fermion_t psi_d, Fermion_t src_d,
			       bfm_evo<float> &bfm_f,
			       bfm_evo<double> &bfm_d,
			       int report_freq = -1)
{

  double f;
  double cp,c,a,d,b;
  int me = bfm_f.thread_barrier();

  //Assumes bfm_d comms active
  mixed_cg::switch_comm(bfm_f, bfm_d);

  Fermion_t src = bfm_f.threadedAllocFermion(mem_slow); 
  mixed_cg::threaded_convFermion(src, src_d, bfm_f, bfm_d);

  Fermion_t psi = bfm_f.threadedAllocFermion(mem_slow); 
  mixed_cg::threaded_convFermion(psi, psi_d, bfm_f, bfm_d);

  Fermion_t p   = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t tmp = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t mp  = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t mmp = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t r   = bfm_f.threadedAllocFermion(mem_fast); 

  Fermion_t tmp_d = bfm_d.threadedAllocFermion(mem_fast); 
  Fermion_t mp_d  = bfm_d.threadedAllocFermion(mem_fast); 
  Fermion_t mmp_d = bfm_d.threadedAllocFermion(mem_fast); 

#define DEALLOCATE_ALL \
      bfm_f.threadedFreeFermion(tmp); \
      bfm_f.threadedFreeFermion(p); \
      bfm_f.threadedFreeFermion(mp); \
      bfm_f.threadedFreeFermion(mmp); \
      bfm_f.threadedFreeFermion(r); \
      bfm_f.threadedFreeFermion(psi); \
      bfm_f.threadedFreeFermion(src); \
      bfm_d.threadedFreeFermion(tmp_d); \
      bfm_d.threadedFreeFermion(mp_d); \
      bfm_d.threadedFreeFermion(mmp_d)

  //Initial residual computation & set up
  double guess = bfm_f.norm(psi);
  d= bfm_f.Mprec(psi,mp,tmp,DaggerNo);
  b= bfm_f.Mprec(mp,mmp,tmp,DaggerYes);

  bfm_f.axpy (r, mmp, src,-1.0);
  bfm_f.axpy (p, mmp, src,-1.0);

  a = bfm_f.norm(p);
  cp= bfm_f.norm(r);

  Float ssq =  bfm_f.norm(src);
  if ( bfm_f.isBoss() && !me ) {
    printf("bfmbase::CGNE_prec gues %le \n",guess);
    printf("bfmbase::CGNE_prec src  %le \n",ssq);
    printf("bfmbase::CGNE_prec  mp  %le \n",d);
    printf("bfmbase::CGNE_prec  mmp %le \n",b);
    printf("bfmbase::CGNE_prec   r  %le \n",cp);
    printf("bfmbase::CGNE_prec   p  %le \n",a);
  }
  Float rsq =  bfm_f.residual* bfm_f.residual*ssq;

  //Check if guess is really REALLY good :)
  if ( cp <= rsq ) {
    if ( bfm_f.isBoss() && !me ) {
      printf("bfmbase::CGNE_prec k=0 converged - suspiciously nice guess %le %le\n",cp,rsq);
    }
    mixed_cg::threaded_convFermion(psi_d, psi, bfm_d, bfm_f);
    mixed_cg::switch_comm(bfm_d, bfm_f);
    DEALLOCATE_ALL;
    
    return 0;
  }
  if ( bfm_f.isBoss() && !me ) 
    printf("bfmbase::CGNE_prec k=0 residual %le rsq %le\n",cp,rsq);

  struct timeval start,stop;
  gettimeofday(&start,NULL);

  for (int k=1;k<=bfm_f.max_iter;k++){
    c=cp; 
    d = bfm_f.Mprec(p,mp,tmp,0,1);
    a = c/d;

    double qq= bfm_f.Mprec(mp,mmp,tmp,1); 
    double b_pred = a*(a*qq-d)/c;

    cp = bfm_f.axpy_norm(r,mmp,r,-a);
    b = cp/c;

    bfm_f.axpy(psi,p,psi,a);
    bfm_f.axpy(p,p,r,b);
    
    if ( k%100 == 0 ){
      if ( bfm_f.isBoss() && !me ) {
	printf("bfmbase::CGNE_prec: k=%d r^2=%le %le\n",k,cp,sqrt(cp/ssq));
      }
    }

    // Stopping condition
    if ( cp <= rsq ) { 

      gettimeofday(&stop,NULL);
      struct timeval diff;
      timersub(&stop,&start,&diff);

      if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec converged in %d iterations\n",k);
      if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec converged in %d.%6.6d s\n",diff.tv_sec,diff.tv_usec);

      //Calculate double precision true resid
      mixed_cg::threaded_convFermion(psi_d, psi, bfm_d, bfm_f);
      mixed_cg::switch_comm(bfm_d, bfm_f);
      
      bfm_d.Mprec(psi_d,mp_d,tmp_d,0);
      bfm_d.Mprec(mp_d,mmp_d,tmp_d,1); 
      bfm_d.axpy(tmp_d,src_d,mmp_d,-1.0);
      double true_residual = sqrt(bfm_d.norm(tmp_d)/bfm_d.norm(src_d));
      if ( bfm_d.isBoss() && !me ) 
	printf("bfmbase::CGNE_prec: true residual is %le \n",true_residual);

      DEALLOCATE_ALL;

      return k;
    }else if(report_freq != -1 && k % report_freq == 0){
      mixed_cg::threaded_convFermion(psi_d, psi, bfm_d, bfm_f);
      mixed_cg::switch_comm(bfm_d, bfm_f);

      bfm_d.Mprec(psi_d,mp_d,tmp_d,0);
      bfm_d.Mprec(mp_d,mmp_d,tmp_d,1); 
      bfm_d.axpy(tmp_d,src_d,mmp_d,-1.0);
      double true_residual = sqrt(bfm_d.norm(tmp_d)/bfm_d.norm(src_d));
      if ( bfm_d.isBoss() && !me ) 
	printf("bfmbase::CGNE_prec: iter %d, double prec true residual is %.12le, running true residual is %.12le\n",k,true_residual,sqrt(cp/rsq)*bfm_f.residual);

      mixed_cg::switch_comm(bfm_f, bfm_d);
    }

  }
  if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec: CG not converged \n");
  DEALLOCATE_ALL;
  mixed_cg::switch_comm(bfm_d, bfm_f);

  return -1;
}
#undef DEALLOCATE_ALL


static int CGNE_singleprec_trueresid_reliable_update(Fermion_t psi_d, Fermion_t src_d,
						     bfm_evo<float> &bfm_f,
						     bfm_evo<double> &bfm_d,
						     int update_freq = 10,
						     int report_freq = -1)
{

  double f;
  double cp,c,a,d,b;
  int me = bfm_f.thread_barrier();

  //Assumes bfm_d comms active
  mixed_cg::switch_comm(bfm_f, bfm_d);

  Fermion_t src = bfm_f.threadedAllocFermion(mem_slow); 
  mixed_cg::threaded_convFermion(src, src_d, bfm_f, bfm_d);

  Fermion_t psi = bfm_f.threadedAllocFermion(mem_slow); 
  mixed_cg::threaded_convFermion(psi, psi_d, bfm_f, bfm_d);

  Fermion_t p   = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t tmp = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t mp  = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t mmp = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t r   = bfm_f.threadedAllocFermion(mem_fast); 

  Fermion_t tmp_d = bfm_d.threadedAllocFermion(mem_fast); 
  Fermion_t mp_d  = bfm_d.threadedAllocFermion(mem_fast); 
  Fermion_t mmp_d = bfm_d.threadedAllocFermion(mem_fast); 

#define DEALLOCATE_ALL \
      bfm_f.threadedFreeFermion(tmp); \
      bfm_f.threadedFreeFermion(p); \
      bfm_f.threadedFreeFermion(mp); \
      bfm_f.threadedFreeFermion(mmp); \
      bfm_f.threadedFreeFermion(r); \
      bfm_f.threadedFreeFermion(psi); \
      bfm_f.threadedFreeFermion(src); \
      bfm_d.threadedFreeFermion(tmp_d); \
      bfm_d.threadedFreeFermion(mp_d); \
      bfm_d.threadedFreeFermion(mmp_d)

  //Initial residual computation & set up
  double guess = bfm_f.norm(psi);

  d= bfm_f.Mprec(psi,mp,tmp,DaggerNo);
  b= bfm_f.Mprec(mp,mmp,tmp,DaggerYes);
  
  bfm_f.axpy (r, mmp, src,-1.0);
  bfm_f.axpy (p, mmp, src,-1.0);
  
  a = bfm_f.norm(p);
  cp= bfm_f.norm(r);

  Float ssq =  bfm_f.norm(src);
  if ( bfm_f.isBoss() && !me ) {
    printf("bfmbase::CGNE_prec gues %le \n",guess);
    printf("bfmbase::CGNE_prec src  %le \n",ssq);
    printf("bfmbase::CGNE_prec  mp  %le \n",d);
    printf("bfmbase::CGNE_prec  mmp %le \n",b);
    printf("bfmbase::CGNE_prec   r  %le \n",cp);
    printf("bfmbase::CGNE_prec   p  %le \n",a);
  }
  Float rsq =  bfm_f.residual* bfm_f.residual*ssq;

  //Check if guess is really REALLY good :)
  if ( cp <= rsq ) {
    if ( bfm_f.isBoss() && !me ) {
      printf("bfmbase::CGNE_prec k=0 converged - suspiciously nice guess %le %le\n",cp,rsq);
    }
    mixed_cg::threaded_convFermion(psi_d, psi, bfm_d, bfm_f);
    mixed_cg::switch_comm(bfm_d, bfm_f);
    DEALLOCATE_ALL;
    
    return 0;
  }
  if ( bfm_f.isBoss() && !me ) 
    printf("bfmbase::CGNE_prec k=0 residual %le rsq %le\n",cp,rsq);

  struct timeval start,stop;
  gettimeofday(&start,NULL);

  for (int k=1;k<=bfm_f.max_iter;k++){
    c=cp; 

    // if(k % update_freq == 0){
    //   //Replace r with true residual
    //   mixed_cg::switch_comm(bfm_d, bfm_f);
    //   mixed_cg::threaded_convFermion(psi_d, psi, bfm_d, bfm_f);
    //   bfm_d.Mprec(psi_d,mp_d,tmp_d,0,1);
    //   bfm_d.Mprec(mp_d,mmp_d,tmp_d,1);
    //   double c_old = c;
    //   c = bfm_d.axpy_norm(tmp_d,src_d,mmp_d,-1.0);
    //   printf("bfmbase::CGNE_prec: reliable update iter %d, replaced |r|^2 = %.12le with |r|^2 = %.12le\n",k,c_old,c);

    //   mixed_cg::threaded_convFermion(r, tmp_d, bfm_f, bfm_d);
    //   mixed_cg::switch_comm(bfm_f, bfm_d);
    // } 

    d = bfm_f.Mprec(p,mp,tmp,0,1);
    a = c/d;
    
    bfm_f.axpy(psi,p,psi,a);

    double qq= bfm_f.Mprec(mp,mmp,tmp,1); 
    double b_pred = a*(a*qq-d)/c;
    
    if(k % update_freq == 0){
      double cp_sp = bfm_f.axpy_norm(r,mmp,r,-a);

      //Replace r with true residual
      mixed_cg::switch_comm(bfm_d, bfm_f);
      mixed_cg::threaded_convFermion(psi_d, psi, bfm_d, bfm_f);
      bfm_d.Mprec(psi_d,mp_d,tmp_d,0,1);
      bfm_d.Mprec(mp_d,mmp_d,tmp_d,1);
      cp = bfm_d.axpy_norm(tmp_d,mmp_d,src_d,-1.0);

      if( bfm_d.isBoss() && !me) printf("bfmbase::CGNE_prec: reliable update iter %d, replaced |r|^2 = %.12le with |r|^2 = %.12le\n",k,cp_sp,cp);

      mixed_cg::threaded_convFermion(r, tmp_d, bfm_f, bfm_d);

      //mixed_cg::threaded_convFermion(p, tmp_d, bfm_f, bfm_d); //make next search direction equal to residual
      mixed_cg::switch_comm(bfm_f, bfm_d);

    }else{
      cp = bfm_f.axpy_norm(r,mmp,r,-a);
      //b = cp/c;    
      //bfm_f.axpy(p,p,r,b);
    }
    b = cp/c;    
    bfm_f.axpy(p,p,r,b);

    if ( k%100 == 0 ){
      if ( bfm_f.isBoss() && !me ) {
	printf("bfmbase::CGNE_prec: k=%d r^2=%le %le\n",k,cp,sqrt(cp/ssq));
      }
    }

    // Stopping condition
    if ( cp <= rsq ) { 

      gettimeofday(&stop,NULL);
      struct timeval diff;
      timersub(&stop,&start,&diff);

      if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec converged in %d iterations\n",k);
      if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec converged in %d.%6.6d s\n",diff.tv_sec,diff.tv_usec);

      //Calculate double precision true resid
      mixed_cg::threaded_convFermion(psi_d, psi, bfm_d, bfm_f);
      mixed_cg::switch_comm(bfm_d, bfm_f);
      
      bfm_d.Mprec(psi_d,mp_d,tmp_d,0);
      bfm_d.Mprec(mp_d,mmp_d,tmp_d,1); 
      bfm_d.axpy(tmp_d,src_d,mmp_d,-1.0);
      double true_residual = sqrt(bfm_d.norm(tmp_d)/bfm_d.norm(src_d));
      if ( bfm_d.isBoss() && !me ) 
	printf("bfmbase::CGNE_prec: true residual is %le \n",true_residual);

      DEALLOCATE_ALL;

      return k;
    }else if(report_freq != -1 && k % report_freq == 0){
      mixed_cg::threaded_convFermion(psi_d, psi, bfm_d, bfm_f);
      mixed_cg::switch_comm(bfm_d, bfm_f);

      bfm_d.Mprec(psi_d,mp_d,tmp_d,0);
      bfm_d.Mprec(mp_d,mmp_d,tmp_d,1); 
      bfm_d.axpy(tmp_d,src_d,mmp_d,-1.0);
      double true_residual = sqrt(bfm_d.norm(tmp_d)/bfm_d.norm(src_d));
      if ( bfm_d.isBoss() && !me ) 
	printf("bfmbase::CGNE_prec: iter %d, double prec true residual is %.12le, running true residual is %.12le  (|r|^2=%.12le)\n",k,true_residual,sqrt(cp/rsq)*bfm_f.residual,cp);

      mixed_cg::switch_comm(bfm_f, bfm_d);
    }

  }
  if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec: CG not converged \n");
  DEALLOCATE_ALL;
  mixed_cg::switch_comm(bfm_d, bfm_f);

  return -1;
}
#undef DEALLOCATE_ALL



static int CGNE_singleprec_trueresid_reliable_update_dp_vects(Fermion_t psi_d, Fermion_t src_d,
							      bfm_evo<float> &bfm_f,
							      bfm_evo<double> &bfm_d,
							      int update_freq = 10,
							      int report_freq = -1)
{
  //Maintain seach vector in double precision

  double f;
  double cp,c,a,d,b;
  int me = bfm_d.thread_barrier();

  //Assumes bfm_d comms active

  Fermion_t src = bfm_f.threadedAllocFermion(mem_slow); 
  mixed_cg::threaded_convFermion(src, src_d, bfm_f, bfm_d);

  Fermion_t psi = bfm_f.threadedAllocFermion(mem_slow); 
  mixed_cg::threaded_convFermion(psi, psi_d, bfm_f, bfm_d);

  Fermion_t p   = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t tmp = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t mp  = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t mmp = bfm_f.threadedAllocFermion(mem_fast); 
  Fermion_t r   = bfm_f.threadedAllocFermion(mem_fast); 

  Fermion_t p_d = bfm_d.threadedAllocFermion(mem_fast); 
  Fermion_t tmp_d = bfm_d.threadedAllocFermion(mem_fast); 
  Fermion_t mp_d  = bfm_d.threadedAllocFermion(mem_fast); 
  Fermion_t mmp_d = bfm_d.threadedAllocFermion(mem_fast); 

#define DEALLOCATE_ALL \
  bfm_f.threadedFreeFermion(p);     \
      bfm_f.threadedFreeFermion(tmp); \
      bfm_f.threadedFreeFermion(mp); \
      bfm_f.threadedFreeFermion(mmp); \
      bfm_f.threadedFreeFermion(r); \
      bfm_f.threadedFreeFermion(psi); \
      bfm_f.threadedFreeFermion(src); \
      bfm_d.threadedFreeFermion(p_d); \
      bfm_d.threadedFreeFermion(tmp_d); \
      bfm_d.threadedFreeFermion(mp_d); \
      bfm_d.threadedFreeFermion(mmp_d)

  //Initial residual computation & set up
  mixed_cg::switch_comm(bfm_f, bfm_d);
  double guess = bfm_f.norm(psi);

  d= bfm_f.Mprec(psi,mp,tmp,DaggerNo);
  b= bfm_f.Mprec(mp,mmp,tmp,DaggerYes);
  
  cp = bfm_f.axpy_norm (r, mmp, src,-1.0);

  mixed_cg::threaded_convFermion(tmp_d, mmp, bfm_d, bfm_f);
  mixed_cg::switch_comm(bfm_d, bfm_f);
  a = bfm_d.axpy_norm (p_d, tmp_d, src_d,-1.0);
  mixed_cg::switch_comm(bfm_f, bfm_d);

  Float ssq =  bfm_f.norm(src);
  if ( bfm_f.isBoss() && !me ) {
    printf("bfmbase::CGNE_prec gues %le \n",guess);
    printf("bfmbase::CGNE_prec src  %le \n",ssq);
    printf("bfmbase::CGNE_prec  mp  %le \n",d);
    printf("bfmbase::CGNE_prec  mmp %le \n",b);
    printf("bfmbase::CGNE_prec   r  %le \n",cp);
    printf("bfmbase::CGNE_prec   p  %le \n",a);
  }
  Float rsq =  bfm_f.residual* bfm_f.residual*ssq;

  //Check if guess is really REALLY good :)
  if ( cp <= rsq ) {
    if ( bfm_f.isBoss() && !me ) {
      printf("bfmbase::CGNE_prec k=0 converged - suspiciously nice guess %le %le\n",cp,rsq);
    }
    mixed_cg::switch_comm(bfm_d, bfm_f);
    DEALLOCATE_ALL;
    
    return 0;
  }
  if ( bfm_f.isBoss() && !me ) 
    printf("bfmbase::CGNE_prec k=0 residual %le rsq %le\n",cp,rsq);

  struct timeval start,stop;
  gettimeofday(&start,NULL);

  for (int k=1;k<=bfm_f.max_iter;k++){
    c=cp; 

    mixed_cg::threaded_convFermion(p, p_d, bfm_f, bfm_d);
    d = bfm_f.Mprec(p,mp,tmp,0,1);
    double qq= bfm_f.Mprec(mp,mmp,tmp,1); 
    a = c/d;
    double b_pred = a*(a*qq-d)/c;

    bfm_d.axpy(psi_d,p_d,psi_d,a);
    
    if(k % update_freq == 0){
      double cp_sp = bfm_f.axpy_norm(r,mmp,r,-a);

      //Replace r with true residual
      mixed_cg::switch_comm(bfm_d, bfm_f);
      bfm_d.Mprec(psi_d,mp_d,tmp_d,0,1);
      bfm_d.Mprec(mp_d,mmp_d,tmp_d,1);
      cp = bfm_d.axpy_norm(tmp_d,mmp_d,src_d,-1.0);

      if( bfm_d.isBoss() && !me) printf("bfmbase::CGNE_prec: reliable update iter %d, replaced |r|^2 = %.12le with |r|^2 = %.12le\n",k,cp_sp,cp);

      mixed_cg::threaded_convFermion(r, tmp_d, bfm_f, bfm_d);
      mixed_cg::switch_comm(bfm_f, bfm_d);
    }else{
      cp = bfm_f.axpy_norm(r,mmp,r,-a);
      mixed_cg::threaded_convFermion(tmp_d, r, bfm_d, bfm_f);
    }

    //Update search direction
    b = cp/c;
    bfm_d.axpy(p_d,p_d,tmp_d,b);

    if ( k%100 == 0 ){
      if ( bfm_f.isBoss() && !me ) {
	printf("bfmbase::CGNE_prec: k=%d r^2=%le %le\n",k,cp,sqrt(cp/ssq));
      }
    }

    // Stopping condition
    if ( cp <= rsq ) { 

      gettimeofday(&stop,NULL);
      struct timeval diff;
      timersub(&stop,&start,&diff);

      if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec converged in %d iterations\n",k);
      if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec converged in %d.%6.6d s\n",diff.tv_sec,diff.tv_usec);

      //Calculate double precision true resid
      mixed_cg::switch_comm(bfm_d, bfm_f);
      
      bfm_d.Mprec(psi_d,mp_d,tmp_d,0);
      bfm_d.Mprec(mp_d,mmp_d,tmp_d,1); 
      bfm_d.axpy(tmp_d,src_d,mmp_d,-1.0);
      double true_residual = sqrt(bfm_d.norm(tmp_d)/bfm_d.norm(src_d));
      if ( bfm_d.isBoss() && !me ) 
	printf("bfmbase::CGNE_prec: true residual is %le \n",true_residual);

      DEALLOCATE_ALL;

      return k;
    }else if(report_freq != -1 && k % report_freq == 0){
      mixed_cg::switch_comm(bfm_d, bfm_f);

      bfm_d.Mprec(psi_d,mp_d,tmp_d,0);
      bfm_d.Mprec(mp_d,mmp_d,tmp_d,1); 
      bfm_d.axpy(tmp_d,src_d,mmp_d,-1.0);
      double true_residual = sqrt(bfm_d.norm(tmp_d)/bfm_d.norm(src_d));
      if ( bfm_d.isBoss() && !me ) 
	printf("bfmbase::CGNE_prec: iter %d, double prec true residual is %.12le, running true residual is %.12le  (|r|^2=%.12le)\n",k,true_residual,sqrt(cp/rsq)*bfm_f.residual,cp);

      mixed_cg::switch_comm(bfm_f, bfm_d);
    }

  }
  if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec: CG not converged \n");
  DEALLOCATE_ALL;
  mixed_cg::switch_comm(bfm_d, bfm_f);

  return -1;
}
#undef DEALLOCATE_ALL





Float* rand_5d_canonical_fermion(Lattice &lat){
  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();
  if(GJP.Gparity()) f_size*=2;
  Float *v1 = (Float *)pmalloc(sizeof(Float) * f_size);
  if(!UniqueID()) printf("Making random gaussian 5d vector\n");
  lat.RandGaussVector((Vector*)v1, 0.5, 2, CANONICAL, FIVE_D);
  if(!UniqueID()) printf("Finished making random gaussian vector\n");
  return v1;
}

int main(int argc, char *argv[])
{
    const char *fname = "main()";

    if(argc<3){
      if(!UniqueID()) printf("Usage: <executable> <vml dir> <test number> <options>\n");
      exit(0);
    }
    int test_no;
    { std::stringstream ss; ss << argv[2]; ss >> test_no; }

    if(test_no == 9) do_arg.start_conf_kind = START_CONF_ORD;

    setup(&argc, &argv);

    Float residual = 1e-10;
    int report_freq = 100;
    bool skip_test1 = false;
    int switch_iter = 90;
    int rel_up_freq = 10;

    bool resid_from_file = false;
    char* resid_file;

   int nthreads = 1; 
#if TARGET == BGQ
   nthreads = 64;
#endif

    int i=3;
    while(i<argc){
      char* cmd = argv[i];
      if( strncmp(cmd,"-residual",15) == 0){
	std::stringstream ss; ss << argv[i+1]; ss >> residual;
	i+=2;
      }else if( strncmp(cmd,"-resid_file",15) == 0){
	resid_from_file = true; resid_file = argv[i+1];
	i+=2;
      }else if( strncmp(cmd,"-report_freq",20) == 0){
	std::stringstream ss; ss << argv[i+1]; ss >> report_freq;
	i+=2;
      }else if( strncmp(cmd,"-switch_iter",20) == 0){
	std::stringstream ss; ss << argv[i+1]; ss >> switch_iter;
	i+=2;
      }else if( strncmp(cmd,"-rel_up_freq",20) == 0){
	std::stringstream ss; ss << argv[i+1]; ss >> rel_up_freq;
	i+=2;
      }else if( strncmp(cmd,"-threads",20) == 0){
	std::stringstream ss; ss << argv[i+1]; ss >> nthreads;
	i+=2;
      }else{
	ERR.General("","main","Unknown argument %s\n",cmd);
      }
    }

    // VRB.Result( "CPS", "main", "omp_get_num_threads[1] -> %d", omp_get_num_threads() );
    
    // #pragma omp parallel
    // {
    //   if ( UniqueID() == 0 && omp_get_thread_num() == 0 ) {
    //     VRB.Result( "CPS", "main", "omp_get_num_threads[2] -> %d", omp_get_num_threads() );
    //   }
    // }

    // VRB.Result( "CPS", "main", "omp_get_num_threads[3] -> %d", omp_get_num_threads() );
    
    GnoneFdwf lattice;

    //Get a RemezArg from RationalQuotientRemezArg. We only use the first fermion MD here
    RemezArg &remez_arg = rq_arg.frm_md.frm_md_val[0];

    Float* src_cps = rand_5d_canonical_fermion(lattice);

    lattice.BondCond();
    Float* gauge = (Float*) lattice.GaugeField();
    
    bfmarg dwfa;
    setup_bfmargs(dwfa,residual,nthreads);
    
    if(!UniqueID()){ printf("sizeof float %d double %d\n",sizeof(float),sizeof(double)); fflush(stdout); }

    bfm_evo<double> bfm_d;
    bfm_d.init(dwfa);
    // sync();
    if(!UniqueID()){ printf("Importing gauge\n"); fflush(stdout); }
    bfm_d.cps_importGauge(gauge);

    //Init the single prec instance
    bfm_d.comm_end();

    bfm_evo<float> bfm_f;
    bfm_f.init(dwfa);
    bfm_f.cps_importGauge(gauge);
    bfm_f.comm_end();

    bfm_d.comm_init();

    // sync();
    if(!UniqueID()){ printf("Allocating fermions\n"); fflush(stdout); }
    Fermion_t src_bfm[2] = {bfm_d.allocFermion(), bfm_d.allocFermion()}; //odd/even
    
    // sync();
    if(!UniqueID()){ printf("Impexing random vector to bfm src vectors\n"); fflush(stdout); }
    bfm_d.cps_impexFermion(src_cps,src_bfm,1);

    //Allocate npoles Fermion vectors in double precision

    // sync();
    if(!UniqueID()){ printf("Allocating solution vectors\n"); fflush(stdout); }
    Fermion_t sol[remez_arg.degree];
    for(int i=0;i<remez_arg.degree;++i) sol[i] = bfm_d.allocFermion();

    Float alpha[remez_arg.degree];
    Float resid[remez_arg.degree];
    for(int i=0;i<remez_arg.degree;++i){
      alpha[i] = 1.0;
      resid[i] = residual;
    }
    if(resid_from_file){
      if(!UniqueID()){ printf("Loading residuals from ActionRationalQuotientArg %s\n",resid_file); fflush(stdout); }
      ActionRationalQuotientArg rq_resid_file;
      if (!rq_resid_file.Decode(resid_file, "rq_resid_file") ) ERR.General(cname, fname, "Bad file\n");      

      Float* resids_in = rq_resid_file.fermions.fermions_val[0].md_approx.stop_rsd.stop_rsd_val;

      for(int i=0;i<remez_arg.degree;++i) resid[i] = resids_in[i];
    }

    // if(!UniqueID()){ printf("Before initial CPS sync\n"); fflush(stdout); }
    // sync();
    // if(!UniqueID()){ printf("After initial CPS sync\n"); fflush(stdout); }   
    

//     if(test_no == 0){
// #pragma omp parallel
//       {
// 	for(int i=0;i<remez_arg.degree;++i) bfm_d.set_zero(sol[i]);
// 	multi_shift_trueresid(sol,src_bfm[0], remez_arg.pole_inv, alpha, remez_arg.degree, resid, 0, bfm_f, bfm_d, report_freq);
//       }
//     }    

//     if(test_no == 1){
//       if(!UniqueID()) printf("Running multi-shift double-prec, single-prec Mprec\n");
    
// #pragma omp parallel
//       {
// 	for(int i=0;i<remez_arg.degree;++i)
// 	  bfm_d.set_zero(sol[i]);

// 	struct timeval start, stop, diff_dpsp, diff_switch, diff_std;

// 	if(!UniqueID() && !omp_get_thread_num()) printf("Starting mixed multi-shift\n");

// 	gettimeofday(&start,NULL); 
// 	//multi_shift_dp_spmprec(sol,src_bfm[0], remez_arg.pole_inv, alpha, remez_arg.degree, resid, 0, bfm_f, bfm_d);
// 	gettimeofday(&stop,NULL); 
// 	timersub(&stop,&start,&diff_dpsp);

// 	for(int i=0;i<remez_arg.degree;++i)
// 	  bfm_d.set_zero(sol[i]);

// 	if(!UniqueID() && !omp_get_thread_num()) printf("Starting switched prec multi-shift\n");
      
// 	gettimeofday(&start,NULL); 
// 	multi_shift_sp_dp_switch(sol,src_bfm[0], remez_arg.pole_inv, alpha, remez_arg.degree, resid, 0, bfm_f, bfm_d, switch_iter,report_freq);
// 	gettimeofday(&stop,NULL); 
// 	timersub(&stop,&start,&diff_switch);

// 	for(int i=0;i<remez_arg.degree;++i)
// 	  bfm_d.set_zero(sol[i]);

// 	if(!UniqueID() && !omp_get_thread_num()) printf("Starting double prec multi-shift\n");

// 	gettimeofday(&start,NULL);
// 	bfm_d.CGNE_prec_MdagM_multi_shift(sol,src_bfm[0], remez_arg.pole_inv, alpha, remez_arg.degree, resid, 0);
// 	gettimeofday(&stop,NULL); 
// 	timersub(&stop,&start,&diff_std);
      
// 	if(!UniqueID() && !omp_get_thread_num()) printf("Mixed time %d.%6.6d s, switched time %d.%6.6d s, double prec time %d.%6.6d s\n",
// 							diff_dpsp.tv_sec,diff_dpsp.tv_usec,
// 							diff_switch.tv_sec,diff_switch.tv_usec,
// 							diff_std.tv_sec,diff_std.tv_usec); 
//       }
//     }

//     if(test_no == 2){
//       if(!UniqueID()) printf("Running single precision CG test\n");
// #pragma omp parallel
//       {
// 	bfm_d.set_zero(sol[0]);
// 	CGNE_singleprec_trueresid(sol[0],src_bfm[0],bfm_f,bfm_d,report_freq);
//       }
//     }	

//     if(test_no == 3){
//       if(!UniqueID()) printf("Running reliable-update single precision CG test\n");
// #pragma omp parallel
//       {
// 	CGNE_singleprec_trueresid_reliable_update(sol[0],src_bfm[0],bfm_f,bfm_d,rel_up_freq,report_freq);
//       }
//     }

//     if(test_no == 4){
//       if(!UniqueID()) printf("Running reliable-update single precision with double vects CG test\n");
// #pragma omp parallel
//       {
// 	struct timeval start, stop, diff;
// 	gettimeofday(&start,NULL);       
// 	CGNE_singleprec_trueresid_reliable_update_dp_vects(sol[0],src_bfm[0],bfm_f,bfm_d,rel_up_freq,report_freq);
// 	gettimeofday(&stop,NULL); 
// 	timersub(&stop,&start,&diff);
// 	if(!UniqueID() && !omp_get_thread_num()) printf("Reliable-update single precision with double vects time %d.%6.6d s\n", diff.tv_sec,diff.tv_usec);
//       }
//     }

//     if(test_no == 5){
//       if(!UniqueID()) printf("Running double precision CG test\n");
// #pragma omp parallel
//       {
// 	struct timeval start, stop, diff;
// 	gettimeofday(&start,NULL);
// 	bfm_d.CGNE_prec(sol[0],src_bfm[0]);
// 	gettimeofday(&stop,NULL); 
// 	timersub(&stop,&start,&diff);
// 	if(!UniqueID() && !omp_get_thread_num()) printf("Double precision CG time %d.%6.6d s\n", diff.tv_sec,diff.tv_usec);
//       }
//     }
    
//     if(test_no == 6){
//       if(!UniqueID()) printf("Running reliable-update single precision multi-shift with double vects\n");
// #pragma omp parallel
//       {
// 	struct timeval start, stop, diff;
// 	gettimeofday(&start,NULL); 
// 	multi_shift_sp_reliable_update_dp_vects(sol,src_bfm[0], remez_arg.pole_inv, alpha, remez_arg.degree, resid, 0, bfm_f, bfm_d, rel_up_freq, report_freq);
// 	gettimeofday(&stop,NULL); 
// 	timersub(&stop,&start,&diff);
// 	if(!UniqueID() && !omp_get_thread_num()) printf("Reliable-update single precision multi-mass with double vects time %d.%6.6d s\n", diff.tv_sec,diff.tv_usec);
//       }
//     }
    
    if(test_no == 7){
      //Test the defect corrected version of the above
      if(!UniqueID()){ printf("Testing reliable-update single precision multi-mass with double vects and defect correction\n"); fflush(stdout); }
#pragma omp parallel
      {
	struct timeval start, stop, diff;
	gettimeofday(&start,NULL); 
	if(!UniqueID() && !omp_get_thread_num()){ printf("Prior to first thread barrier\n"); fflush(stdout); }
	int me = bfm_d.thread_barrier();
	if(!UniqueID() && !omp_get_thread_num()){ printf("After first thread barrier\n"); fflush(stdout); }

	// if ( bfm_d.isBoss() && !me ){ printf("MAIN: 1\n"); fflush(stdout); }
	// //TESTING
	// bfm_d.thread_barrier();
	// if(!me) bfm_d.comm_spi_barrier();
	// bfm_d.thread_barrier();
	// //TESTING
	// if ( bfm_d.isBoss() && !me ){ printf("MAIN: 1.5\n"); fflush(stdout); }

	mixed_cg::threaded_cg_mixed_multi_shift_MdagM_sp_relup_dp_defect_correction(sol,src_bfm[0], remez_arg.pole_inv, alpha, remez_arg.degree, resid, 0, bfm_f, bfm_d, rel_up_freq, report_freq);
	
	gettimeofday(&stop,NULL); 
	timersub(&stop,&start,&diff);
	if(!UniqueID() && !omp_get_thread_num()) printf("Reliable-update single precision multi-mass with double vects time %d.%6.6d s\n", diff.tv_sec,diff.tv_usec);
      }
    }


    if(test_no == 8){
      if(!UniqueID()){ printf("Running double precision multi-shift test\n"); fflush(stdout); }
#pragma omp parallel
      {
	if(!UniqueID() && omp_get_thread_num()==0){ printf("Started DOUBLE_PREC multi-shift with %d threads\n",omp_get_num_threads()); fflush(stdout); }
	struct timeval start, stop, diff;
	gettimeofday(&start,NULL);
	// if(!UniqueID() && omp_get_thread_num()==0){ printf("Prior to initial thread barrier\n"); fflush(stdout); }
	// bfm_d.thread_barrier();
	// if(!UniqueID() && omp_get_thread_num()==0){ printf("After initial thread barrier\n"); fflush(stdout); }

	// if(!UniqueID() && omp_get_thread_num()==0){ printf("Prior to initial node sync\n"); fflush(stdout); }
	// double pooh(1.0);
	// bfm_d.comm_gsum(pooh);
	// if(!UniqueID() && omp_get_thread_num()==0){ printf("After initial node sync\n"); fflush(stdout); }

	bfm_d.CGNE_prec_MdagM_multi_shift(sol,src_bfm[0], remez_arg.pole_inv, alpha, remez_arg.degree, resid, 0);

	gettimeofday(&stop,NULL); 
	timersub(&stop,&start,&diff);
	if(!UniqueID() && !omp_get_thread_num()){ printf("Double precision multi-shift time %d.%6.6d s\n", diff.tv_sec,diff.tv_usec); fflush(stdout); }
      }
    }
    
//     if(test_no == 9){
//       if(!UniqueID()) printf("Running convfermion test\n");
      
//       Fermion_t r1 = bfm_f.allocFermion();
//       Fermion_t r2 = bfm_f.allocFermion();

//       Fermion_t diff = bfm_f.allocFermion();      
// #pragma omp parallel
//       {
// 	struct timeval start, stop, diff1, diff2;
// 	gettimeofday(&start,NULL);

// 	mixed_cg::threaded_convFermion(r1,src_bfm[0],bfm_f,bfm_d);

// 	gettimeofday(&stop,NULL); 
// 	timersub(&stop,&start,&diff1);	
// 	gettimeofday(&start,NULL);

// 	mixed_cg::threaded_convFermion_fast(r2,src_bfm[0],bfm_f,bfm_d);

// 	gettimeofday(&stop,NULL); 
// 	timersub(&stop,&start,&diff2);	

// 	double dnorm = bfm_f.axpy_norm(diff, r1, r2, -1.0);
// 	double rnorm = bfm_f.norm(r1);
// 	if(!UniqueID() && !omp_get_thread_num()) printf("Norm of difference %.12le, norm of vector %.12le, orig time %d.%6.6d s, fast time %d.%6.6d s\n",dnorm,rnorm,diff1.tv_sec,diff1.tv_usec,diff2.tv_sec,diff2.tv_usec);
//       }
//     }


    VRB.Result(cname, fname, "Program ended normally.\n");
    
    return 0;
}
