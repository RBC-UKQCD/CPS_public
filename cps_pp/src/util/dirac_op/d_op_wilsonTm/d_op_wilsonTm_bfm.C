#include <config.h>
#include <stdio.h>
#include <stdarg.h>
#ifdef USE_OMP
#include <omp.h>
#endif

#ifdef USE_BFM_TM
#ifdef USE_CHROMA
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#endif
#include <bfm.h>
#include <bfm_qdp.h>
//#include <qdp.h>
#include <util/gjp.h>
#include <comms/sysfunc_cps.h>
#include <util/error.h>
#include <util/verbose.h>
#include <util/smalloc.h>
#include <util/time_cps.h>
#include <util/dirac_op/d_op_dwf.h>
#include <util/wilson.h>
#include <util/lattice/bfm_evo.h>

//static int  Printf(char *format,...){}
#define Printf if ( !UniqueID() ) printf
//#define Printf printf


USING_NAMESPACE_CPS

void importGauge(
CPS_NAMESPACE::Lattice &Lat,
multi1d<LatticeColorMatrix> &U,
CPS_NAMESPACE::Float *gauge,
int dag);

void impexFermion(
int if_export,
Lattice &Lat,
multi1d<LatticeFermion> const &qdp,
CPS_NAMESPACE::Float *cps_p,
int even, int odd,int Ls=0, double fac_p=1.);


CPS_START_NAMESPACE
static void print_vec(void *p, char *name){
 Float *tmp_p = (Float *)p;
 if (!UniqueID()) printf("%s(%p) = %g\n",name, p, *tmp_p);
}
int cps_qdp_finalize();
int cps_qdp_init(int *argc, char ***argv);

int DiracOpWilsonTm::InvCg(Vector *out,
                   Vector *in,
		   Float src_norm_sq,
                   Float *true_res){
#if 0
   return DiracOp::InvCg(out,in,src_norm_sq,true_res);
#else

  const char *fname="InvCg(V*,V*,F,F*) [bfm version]";
  double cg_time = -dclock();
  size_t f_size_cb =  GJP.VolNodeSites() * lat.FsiteSize() / (lat.FchkbEvl()+1);
  if(GJP.Gparity()) f_size_cb*=2;
  
  src_norm_sq = in->NormSqNode(f_size_cb);
  DiracOpGlbSum(&src_norm_sq);

  VRB.Flow(cname,fname,"src_norm_sq=%e\n",src_norm_sq);
  
// Calculate stopping condition
//------------------------------------------------------------------
  Float stp_cnd = src_norm_sq * dirac_arg->stop_rsd * dirac_arg->stop_rsd;
  VRB.Result(cname,fname,
           "stp_cnd =%e\n", IFloat(stp_cnd));

  int iter=0;   

  {
     Float mass = dirac_arg->mass;
     Float epsilon = dirac_arg->epsilon;
     Float *sol_cps = (Float*) out;
     Float *src_cps = (Float*) in;
     Float residual = dirac_arg->stop_rsd;
     int max_iter = dirac_arg->max_num_iter;
  
    if (GJP.Snodes() != 1)
      ERR.General("",fname,"Snodes()(%d)!=1",GJP.Snodes());
  
    cps_qdp_init(GJP.argc_p(), GJP.argv_p());
    double M5 = GJP.DwfHeight();
  
    Printf("src[0]=%g\n",*src_cps);
    int lx = QDP::Layout::subgridLattSize()[0];
    int ly = QDP::Layout::subgridLattSize()[1];
    int lz = QDP::Layout::subgridLattSize()[2];
    int lt = QDP::Layout::subgridLattSize()[3];

    int nferm = 1;
    if(GJP.Gparity()) nferm = 2; //two flavours separately
    multi1d<LatticeFermion> src_qdp(nferm);
  
    /********************************************************
     * Setup Dirac operator
     ********************************************************
     */
    bfmarg wilsa;
    bfm_qdp<double>  wils;
  
    wilsa.node_latt[0]  = lx;
    wilsa.node_latt[1]  = ly;
    wilsa.node_latt[2]  = lz;
    wilsa.node_latt[3]  = lt;
    wilsa.verbose=0;
    wilsa.reproduce=0;

    bfmarg::Threads(GJP.Nthreads());
    bfmarg::Reproduce(0);
    bfmarg::ReproduceChecksum(0);
    bfmarg::ReproduceMasterCheck(0);
    bfmarg::Verbose(1);

    static int CGcount = 0;
    CGcount++;
    int test_freq = GJP.CGreprodFreq();
    if (test_freq && (CGcount % test_freq == 0)) {
      wilsa.reproduce=1;
      bfmarg::Reproduce(1);
    }
  
    multi1d<int> procs = QDP::Layout::logicalSize();
    Printf("%d dim machine\n\t", procs.size());
    for(int mu=0;mu<4;mu++){
      Printf("%d ", procs[mu]);
      if ( procs[mu]>1 ) {
        wilsa.local_comm[mu] = 0;
      } else { 
        wilsa.local_comm[mu] = 1;
      }
    }
    Printf("\nLocal comm = ");
    for(int mu=0;mu<4;mu++){
      Printf("%d ", wilsa.local_comm[mu]);
    }
    Printf("\n");
    
    multi1d<int> ncoor = QDP::Layout::nodeCoord();
  
#ifdef BFM_GPARITY
    if(GJP.Gparity()){
      wilsa.gparity = 1;
      Printf("G-parity directions: ");
      for(int d=0;d<3;d++)
	if(GJP.Bc(d) == BND_CND_GPARITY){ wilsa.gparity_dir[d] = 1; Printf("%d ",d); }
	else wilsa.gparity_dir[d] = 0;
      for(int d=0;d<4;d++){
	wilsa.nodes[d] = procs[d];
	wilsa.ncoor[d] = ncoor[d];
      }
      Printf("\n");
    }
#endif
  
    GJP.SetNthreads();
    
    Float *gauge = (Float*) lat.GaugeField();
    wilsa.precon_5d = 0;
    wilsa.Ls   = 1;
    wilsa.solver = WilsonTM;
    wilsa.M5   = 0.0;
    wilsa.mass = toDouble(mass);
    wilsa.twistedmass = toDouble(epsilon);
    wilsa.Csw  = 0.0;
#if 0
    wilsa.list_engine=0;
    wilsa.list_length=0;
#endif
    wilsa.max_iter = max_iter;
    wilsa.residual = toDouble(residual);

    for(int i = 0;i<1;i++){
      Printf("Initialising bfm operator\n");
      wils.init(wilsa);
    
      {
	int nlatt = Nd;
	if(GJP.Gparity()) nlatt*=2;

        multi1d<LatticeColorMatrix> U(nlatt);
        importGauge(lat,U,gauge,1);
        wils.importGauge(U);
      }
    
      Printf("Setup gauge field\n");
      /********************************************************
       * Gaussian source and result vectors
       ********************************************************
       */
    
      Fermion_t src_bfm = wils.allocFermion();
      Printf("src_bfm=%p \n",src_bfm);
      Fermion_t psi[1];
      psi[0] = wils.allocFermion();
    
      Printf("psi[0]=%p \n",psi[0]);
      double fac= 0.25/(kappa*kappa);

      for(int i = 0;i<1;i++){      
        /********************************************************
         * Import gauge field to BAGEL
         ********************************************************
         */
	if(!GJP.Gparity()){
	  impexFermion(0,lat,src_qdp,sol_cps,0,1,1,1./fac);
	  wils.importFermion(src_qdp[0],psi[0],1);
	  impexFermion(0,lat,src_qdp,src_cps,0,1,1);
	  wils.importFermion(src_qdp[0],src_bfm,1);
	}else{
	  //CK: src_qdp contains more than 2 entries corresponding to the 2 flavour fields
	  impexFermion(0,lat,src_qdp,sol_cps,0,1,1,1./fac);
	  //Printf("QDP Norms psi %g, %g\n", toDouble(norm2(src_qdp[0])), toDouble(norm2(src_qdp[1])) );
	  wils.importFermion(src_qdp,psi[0],1);

	  impexFermion(0,lat,src_qdp,src_cps,0,1,1);
	  //Printf("QDP Norms src %g, %g\n", toDouble(norm2(src_qdp[0])), toDouble(norm2(src_qdp[1])) );
	  wils.importFermion(src_qdp,src_bfm,1);
	}

        Float *tmp_p = (Float *)src_bfm;        
	Printf("src_bfm[0]=%g, norm %g\n",*tmp_p, wils.norm(src_bfm));

        Printf("Calling half cb inverter\n"); fflush(stdout);
        wils.inv_type=CG_PREC_MDAGM;
        wils.qdp_chi_h[0]=psi[0];
        wils.qdp_chi_h[1]=psi[0];
        wils.qdp_psi_h[0]=src_bfm;
        wils.qdp_psi_h[1]=src_bfm;
        bfm_spawn_cg(wils);
        tmp_p = (Float *)psi[0];
        Printf("psi[0]=%g\n",*(tmp_p));
	if(!GJP.Gparity()){
	  wils.exportFermion(src_qdp[0],psi[0],1);
	}else{
	  wils.exportFermion(src_qdp,psi[0],1);
	}
      }
      
      impexFermion(1,lat,src_qdp,sol_cps,0,1,1,fac);
      Printf("sol[0]=%g\n",*sol_cps);
      
      Printf("src_bfm=%p freed\n",src_bfm);
      wils.freeFermion(src_bfm);
      Printf("psi[0]=%p freed\n",psi[0]);
      wils.freeFermion(psi[0]);
      Printf("Done\n"); 
      wils.end();
    }
    cps_qdp_finalize();
  }

  //------------------------------------------------------------------
  // Done. Finish up and return
  //------------------------------------------------------------------
    // Calculate and set true residual: 
    // true_res = |src - MatPcDagMatPc * sol| / |src|
  {
    Vector *sol = (Vector *) out;
    Vector *src = (Vector *) in;
    Vector *mmp = (Vector *) smalloc(cname,fname,"mmp",f_size_cb * sizeof(Float));
    Vector *res = (Vector *) smalloc(cname,fname,"res",f_size_cb * sizeof(Float));
    MatPcDagMatPc(mmp, sol);
    res->CopyVec(src, f_size_cb);
    res->VecMinusEquVec(mmp, f_size_cb);
    Float res_norm_sq_cur = res->NormSqNode(f_size_cb);
    DiracOpGlbSum(&res_norm_sq_cur);
    Float tmp = res_norm_sq_cur / src_norm_sq;
    tmp = sqrt(tmp);
    if(true_res != 0){
      *true_res = tmp;
    }
    VRB.Result(cname,fname,
  	     "True |res| / |src| = %e, iter = %d\n", IFloat(tmp), iter);
    sfree(cname,fname,"mmp",mmp);
    sfree(cname,fname,"res",res);
  }
  cg_time +=dclock();
  print_flops(cname,fname,0,cg_time);
  return iter;
#endif
}




//Added by CK
int DiracOpWilsonTm::MInvCG(Vector **out, Vector *in, Float in_norm, Float *mass, 
		    int Nmass, int isz, Float *RsdCG,
		    MultiShiftSolveType type, Float *alpha)
{
  char *fname="MInvCG(V**,V*,F,F*,i,i,F*,t,F*)[bfm WilsonTm version]";

  VRB.Func(cname,fname);
  size_t f_size_cb =  GJP.VolNodeSites() * lat.FsiteSize() / 2;
  if(GJP.Gparity()) f_size_cb*=2;

  VRB.Result(cname,fname,"cps_qdp_init(GJP.argc_p(), GJP.argv_p())");
  cps_qdp_init(GJP.argc_p(), GJP.argv_p());
  int iter=0;

  double minv_time = -dclock();
  Float src_norm_sq = in->NormSqNode(f_size_cb);
  DiracOpGlbSum(&src_norm_sq);

  {
    if (( type != SINGLE ) && ( type != MULTI )){
      ERR.General(cname,fname,"( type != SINGLE ) && ( type != MULTI )\n");
    }

    /********************************************************
     * Setup QDP
     ********************************************************
     */

    int lx = QDP::Layout::subgridLattSize()[0];
    int ly = QDP::Layout::subgridLattSize()[1];
    int lz = QDP::Layout::subgridLattSize()[2];
    int lt = QDP::Layout::subgridLattSize()[3];

    /********************************************************
     * Setup WILS operator
     ********************************************************
     */
    bfmarg wilsa;
    GJP.SetNthreads();

    wilsa.solver = WilsonTM;
    wilsa.node_latt[0]  = lx;
    wilsa.node_latt[1]  = ly;
    wilsa.node_latt[2]  = lz;
    wilsa.node_latt[3]  = lt;
    wilsa.verbose=1;
    wilsa.reproduce=0;

    bfmarg::Threads(GJP.Nthreads());
    bfmarg::Reproduce(0);
    bfmarg::ReproduceChecksum(0);
    bfmarg::ReproduceMasterCheck(0);
    bfmarg::Verbose(1);
    bfmarg::onepluskappanorm = 0;

    multi1d<int> ncoor = QDP::Layout::nodeCoord();
    multi1d<int> procs = QDP::Layout::logicalSize();

#ifdef BFM_GPARITY
    if(GJP.Gparity()){
      wilsa.gparity = 1;
      Printf("G-parity directions: ");
      for(int d=0;d<3;d++)
	if(GJP.Bc(d) == BND_CND_GPARITY){ wilsa.gparity_dir[d] = 1; Printf("%d ",d); }
	else wilsa.gparity_dir[d] = 0;
      for(int d=0;d<4;d++){
	wilsa.nodes[d] = procs[d];
	wilsa.ncoor[d] = ncoor[d];
      }
      Printf("\n");
    }
#endif


    Printf("%d dim machine\n\t", procs.size());
    for(int mu=0;mu<4;mu++){
      Printf("%d ", procs[mu]);
      if ( procs[mu]>1 ) wilsa.local_comm[mu] = 0;
      else wilsa.local_comm[mu] = 1;
    }
    Printf("\nLocal comm = ");
    for(int mu=0;mu<4;mu++){
      Printf("%d ", wilsa.local_comm[mu]);
    }
    Printf("\n");

    double mq= dirac_arg->mass;
    double epsilon = dirac_arg->epsilon;
    Printf("mq=%g epsilon=%g\n",mq,epsilon);

    wilsa.precon_5d = 0;
    wilsa.Ls   = 1;
    wilsa.M5   = 0.0;
    wilsa.mass = toDouble(mq);
    wilsa.twistedmass = toDouble(epsilon);
    wilsa.Csw  = 0.0;
    wilsa.max_iter = dirac_arg->max_num_iter;
    wilsa.residual = RsdCG[0];
    Printf("Initialising bfm operator\n");

    bfm_qdp<double>  wils;
    wils.init(wilsa);

#define NMULTI (3)
    const int MaxShift = 20;
    if(Nmass>MaxShift)
      ERR.General(cname,fname,"Nmass(%d)>MaxShift(%d)\n",Nmass,MaxShift);

    Fermion_t sol_bfm[MaxShift];
    double masses[MaxShift];
    double alpha_bfm [MaxShift];
    double residuals[MaxShift];
   
    //In InvCg BFM solution vector has to be multiplied by 0.25/kappa^2 when converting to cps
    //   sol = (MdagM)^-1 src
    //Thus (MdagM)^-1 BFM is 4*kappa^2 larger than (MdagM)^-1 CPS
    //Thus MdagM_BFM = 0.25/kappa^2 * MdagM CPS
    //Here we invert  MdagM + shift
    //We therefore need to convert the shift to  shift*0.25/kappa^2  such that the normalization matches that of the MdagM operator

    //Residuals
    // true_res = |src - (MatPcDagMatPc + shift) * sol| / |src|
    //BFM calculation;  // true_res = |src - (MatPcDagMatPc + shift)_BFM * sol_BFM| / |src| = |src - (MatPcDagMatPc + shift)_CPS*0.25/kappa^2 * sol_BFM| / |src| = |src - (MatPcDagMatPc + shift)_CPS * sol_CPS| / |src|
    //Hence residuals should agree between BFM and CPS
    double fac= 0.25/(kappa*kappa);

    //For negative mass, CPS sets the shift of the inverted matrix to zero. Without this fix the BFM code has convergence issues
    //Fix this here
    #define CK_NEGATIVE_SHIFT_FIX 0
    
    int nshift_bfm = Nmass;
    if(mass[0] >= 0.0 || (mass[0]<0.0 && !CK_NEGATIVE_SHIFT_FIX) ){
      for(int i =0; i<nshift_bfm ;i++){
	sol_bfm[i] = wils.allocFermion();
	wils.master_fill(sol_bfm[i],0.0);
	masses[i] = fac*mass[i];
	Printf("sol_bfm[%d]=%p\n",i,sol_bfm[i]);
	alpha_bfm[i] = 1.;
	residuals[i]=RsdCG[i];
	
	Printf("%d: shift=%g alpha_bfm=%g RsdCG=%g\n",i,mass[i],alpha_bfm[i],RsdCG[i]);
      }
    }else{
      nshift_bfm = Nmass+1;

      for(int i =0; i<nshift_bfm;i++){
	sol_bfm[i] = wils.allocFermion();
	wils.master_fill(sol_bfm[i],0.0);

	if(i==0) masses[i] = 0.0;
	else masses[i] = fac*mass[i-1];

	Printf("sol_bfm[%d]=%p\n",i,sol_bfm[i]);
	alpha_bfm[i] = 1.;

	if(i==0) residuals[i] = RsdCG[i];
	else residuals[i]=RsdCG[i-1];
	
	if(i==0) Printf("%d: shift=%g alpha_bfm=%g RsdCG=%g\n",i,0.0,alpha_bfm[i],RsdCG[0]);
	else Printf("%d: shift=%g alpha_bfm=%g RsdCG=%g\n",i,mass[i-1],alpha_bfm[i],RsdCG[i-1]);
      }
    }



    /********************************************************
     * Import gauge field to BAGEL
     ********************************************************
     */

    Float *gauge = (Float*) lat.GaugeField();
    {
      int nlatt = Nd;
      if(GJP.Gparity()) nlatt*=2;
      
      multi1d<LatticeColorMatrix> U(nlatt);
      importGauge(lat,U,gauge,1);
      wils.importGauge(U);
    }

    Float *src_cps = (Float*) in;
    Fermion_t src_bfm = wils.allocFermion();
    wils.master_fill(src_bfm,0.0);

    int nferm = 1;
    if(GJP.Gparity()) nferm=2;
      
    {
      multi1d<LatticeFermion> src_qdp(nferm);
      impexFermion(0,lat,src_qdp,src_cps,0,1,1);
      wils.importFermion(src_qdp,src_bfm,1);
    }

    Printf("Calling multi-shift inverter\n");

    wils.inv_type=CG_PREC_MDAGM_MULTI;
    wils.qdp_psi_h[0]=src_bfm;
    wils.qdp_psi_h[1]=src_bfm;
    wils.qdp_chi_multi_h=sol_bfm;
    wils.qdp_chi_h[0]=sol_bfm[0];
    wils.qdp_chi_h[1]=sol_bfm[0];
    wils.shifts=masses;
    wils.alpha =alpha_bfm;
    wils.nshift=nshift_bfm;
    wils.mresidual=residuals;
    wils.single=0;

    double bfm_time = -dclock();
    bfm_spawn_cg(wils);
    bfm_time  += dclock();
    print_flops(fname,"bfm",0,bfm_time);

    print_vec(sol_bfm[0],"sol_bfm[0]");
    bfm_time = -dclock();
    {
      multi1d<LatticeFermion> sol_qdp(nferm);
      int n_sol=Nmass;

      VRB.Result(cname,fname,"type=%d\n",type);

      Float *sol_cps = (Float *)fmalloc(cname,fname,"sol_cps",f_size_cb*2*sizeof(Float));
      memset(sol_cps,0,f_size_cb*2*sizeof(Float));
      Printf("sol_cps=%p\n",sol_cps);
      Vector *sol = (Vector *) sol_cps;
      Vector *src = (Vector *) in;
      Vector *mmp = (Vector *) smalloc(cname,fname,"mmp",f_size_cb * sizeof(Float));
      Vector *res = (Vector *) smalloc(cname,fname,"res",f_size_cb * sizeof(Float));

      double fac= 0.25/(kappa*kappa);
      
      int solvec_offset = 0;
      if(mass[0]<0.0 && CK_NEGATIVE_SHIFT_FIX){
	//Discard the first solution
	Printf("sol_bfm[%d]=%p freed\n",0,sol_bfm[0]);
	wils.freeFermion(sol_bfm[0]);
	solvec_offset = 1;
      }

      for(int i=0;i<n_sol;i++){
	int bfm_idx = i+solvec_offset;

	wils.exportFermion(sol_qdp,sol_bfm[bfm_idx],1);
	print_vec(sol_bfm[bfm_idx],"sol_bfm");
	impexFermion(1,lat,sol_qdp,sol_cps,0,1,1,fac);
	print_vec(sol_cps,"sol_cps");
	Float coef = 1.;
	// calculating true residual. Agrees with BFM, so turned off.
	if (0) { //CK: enabled for test
	  MatPcDagMatPc(mmp, sol);
	  mmp -> FTimesV1PlusV2(mass[i],sol,mmp, f_size_cb);
	  res->CopyVec(src, f_size_cb);
	  res->VecMinusEquVec(mmp, f_size_cb);
	  Float res_norm_sq_cur = res->NormSqNode(f_size_cb);
	  DiracOpGlbSum(&res_norm_sq_cur);
	  Float tmp = res_norm_sq_cur / src_norm_sq;
	  tmp = sqrt(tmp);
	  VRB.Result(cname,fname,
  		     "True %d: |res| / |src| = %e, iter = %d\n", i,IFloat(tmp), iter);
	}

	if (alpha) coef=alpha[i];
	if ( type == SINGLE ){
	  Float *tmp_p = (Float *)out[0];
#pragma omp parallel for default(shared)
	  for(long j=0;j<f_size_cb;j++){
	    tmp_p[j] +=coef*sol_cps[j];
	  }
	} else {
	  Float *tmp_p = (Float *)out[i];
#pragma omp parallel for default(shared)
	  for(long j=0;j<f_size_cb;j++)
	    tmp_p[j] =coef*sol_cps[j];
			
	}
	Printf("sol_bfm[%d]=%p freed\n",bfm_idx,sol_bfm[bfm_idx]);
	wils.freeFermion(sol_bfm[bfm_idx]);
      }
      ffree(cname,fname,"sol_cps",sol_cps);
      ffree(cname,fname,"mmp",mmp);
      ffree(cname,fname,"res",res);
    }
    bfm_time  += dclock();
    print_flops(fname,"convert",0,bfm_time);
  
    Printf("src_bfm=%p freed\n",src_bfm);
    wils.freeFermion(src_bfm);
    Printf("Done\n"); 
    wils.end();
    iter= wils.iter;
  }
  minv_time += dclock();
  print_flops(cname,fname,0,minv_time);

  cps_qdp_finalize();
  return iter;

}


#if 0 //Cant seem to get working
void DiracOpWilsonTm::CalcHmdForceVecs(Vector *chi)
{
  const char *fname = "CalcHmdForceVecs(V*) [bfm version]" ;
  VRB.Func(cname,fname) ;

  if (f_out == 0)
    ERR.Pointer(cname, fname, "f_out") ;

  if (f_in == 0)
    ERR.Pointer(cname, fname, "f_in") ;

//------------------------------------------------------------------
// f_out stores (chi,rho) (the v of eqn B12): 
// rho = gamma_5(-theta) Dslash chi
// f_in stores (psi,sigma) (the w of eqn B11):
// psi = gamma_5(-theta) MatPc chi
// sigma = gamma_5(-theta) Dslash psi
//------------------------------------------------------------------

  Vector *chi_new, *rho, *psi, *sigma ;

  int vol =  ((Wilson *)wilson_lib_arg)->vol[0];
  size_t f_size_cb = 12 * GJP.VolNodeSites() ;
  if(GJP.Gparity()){ 
    vol*=2;
    f_size_cb *= 2; //Layout is   |   odd   |   even  |
                    //            | f0 | f1 | f0 | f1 |
                    //where for each checkerboard, each flavour field occupies one half-volume
  }

  chi_new = f_out ;
  rho = (Vector *)((Float *)f_out + f_size_cb) ;
  psi = f_in ;
  sigma = (Vector *)((Float *)f_in + f_size_cb) ;

  //Start BFM engines!
  bfmarg wilsa;
  GJP.SetNthreads();

  wilsa.solver = WilsonTM;
  for(int i=0;i<4;i++) wilsa.node_latt[i] = QDP::Layout::subgridLattSize()[i];
  wilsa.verbose=1;
  wilsa.reproduce=0;

  bfmarg::Threads(GJP.Nthreads());
  bfmarg::Reproduce(0);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(1);
  bfmarg::onepluskappanorm = 0;

  multi1d<int> ncoor = QDP::Layout::nodeCoord();
  multi1d<int> procs = QDP::Layout::logicalSize();

  if(GJP.Gparity()){
    wilsa.gparity = 1;
    for(int d=0;d<3;d++)
      if(GJP.Bc(d) == BND_CND_GPARITY){ wilsa.gparity_dir[d] = 1; Printf("%d ",d); }
      else wilsa.gparity_dir[d] = 0;
    for(int d=0;d<4;d++){
      wilsa.nodes[d] = procs[d];
      wilsa.ncoor[d] = ncoor[d];
    }
  }
  for(int mu=0;mu<4;mu++){
    if ( procs[mu]>1 ) wilsa.local_comm[mu] = 0;
    else wilsa.local_comm[mu] = 1;
  }

  double mq= dirac_arg->mass;
  double epsilon = dirac_arg->epsilon;

  wilsa.precon_5d = 0;
  wilsa.Ls   = 1;
  wilsa.M5   = 0.0;
  wilsa.mass = toDouble(mq);
  wilsa.twistedmass = toDouble(epsilon);
  wilsa.Csw  = 0.0;
  wilsa.max_iter = dirac_arg->max_num_iter;
  wilsa.residual = 1e-08;

  bfm_evo<double>  wils;
  wils.init(wilsa);

  //In InvCg BFM solution vector has to be multiplied by 0.25/kappa^2 when converting to cps
  //   sol = (MdagM)^-1 src
  //Thus (MdagM)^-1 BFM is 4*kappa^2 larger than (MdagM)^-1 CPS
  //Thus MdagM_BFM = 0.25/kappa^2 * MdagM CPS
  //Here we invert  MdagM + shift
  //We therefore need to convert the shift to  shift*0.25/kappa^2  such that the normalization matches that of the MdagM operator
    
  wils.cps_importGauge((Float*)lat.GaugeField());
  
  Fermion_t chi_bfm = wils.allocFermion();
  Fermion_t rho_bfm = wils.allocFermion();
  Fermion_t psi_bfm = wils.allocFermion();
  Fermion_t sigma_bfm = wils.allocFermion();
  Fermion_t tmp = wils.allocFermion();
  Fermion_t tmp2 = wils.allocFermion();

  //import fermions
  wils.cps_impexcbFermion((Float*)chi, chi_bfm, 1, 1); //odd
  wils.cps_impexcbFermion((Float*)rho, rho_bfm, 1, 0); //even
  wils.cps_impexcbFermion((Float*)psi, psi_bfm, 1, 1); //odd
  wils.cps_impexcbFermion((Float*)sigma, sigma_bfm, 1, 0); //even
  
  chi_new->CopyVec(chi, f_size_cb) ;
#pragma omp parallel    
  {
    /*BFM:  Mprec = Doo = Moo-MoeMee^{-1}Meo
     * WilsonTM : Mee = Moo = 4+m + i tm g5
     *            Meeinv    = (4+m - i tm g5)/ ( (4+m)^2 + tm^2 ) 
     *                      =    (4+m)
     *                        -----------------   .  ( 1 -i tm/(4+m) g5 )
     *                        ( (4+m)^2 + tm^2 )
     *
     * Meo = -1/2 Deo
     */
    /*CPS: 1_OO - kappa^2 * g5theta(ctheta,-stheta) * Dslash_0E * g5theta(ctheta,-stheta) * Dslash_E0*/
    //g5theta(ctheta,stheta) = (ctheta + i stheta gamma_5)M

    //Orig code:
    //MatPc(psi,chi) ;
    //g5theta(psi, vol, ctheta, stheta);
    //psi->VecTimesEquFloat(-kappa*kappa,f_size_cb) ;
    //Is
    //=  -kappa^2{ (ctheta + i stheta g5)(1 - kappa^2 * g5theta(ctheta,-stheta) * Dslash_0E * g5theta(ctheta,-stheta) * Dslash_E0) psi }
    //=  -kappa^2{ (ctheta + i stheta g5) - kappa^2 (ctheta + i stheta g5)(ctheta - i stheta g5) Dslash_0E * g5theta(ctheta,-stheta) * Dslash_E0) psi }
    //=  -kappa^2{ (ctheta + i stheta g5) - kappa^2 (ctheta^2 + stheta^2) Dslash_0E * g5theta(ctheta,-stheta) * Dslash_E0) psi  } 
    //=  -kappa^2{ (ctheta + i stheta g5) - kappa^2 ( (m+4)^2 + epsilon^2 )kappa^2 Dslash_0E * g5theta(ctheta,-stheta) * Dslash_E0) psi  }
    //=  -kappa^2{ (ctheta + i stheta g5) - 1/4 kappa^2 Dslash_0E * g5theta(ctheta,-stheta) * Dslash_E0) psi }

    wils.Mprec(chi_bfm, psi_bfm, tmp, 1, 0); //ODD_OUT DAG_NO

    //CPS kappa = 1/[2*sqrt( (4+m)^2 + epsilon^2 )]
    //    ctheta = (m+4)*kappa
    //    stheta = epsilon * kappa

    //(BFM  Mprec) = (4+m + i tm g5) - 1/4 D_oe [(4+m - i tm g5)/ ( (4+m)^2 + tm^2 )] D_eo
    //             = (4+m + i tm g5) - D_oe kappa^2 ( (4+m) -i tm g5 ) D_eo
    //kappa^2 (BFM  Mprec) = (ctheta + i stheta g5) - kappa^2 D_oe ( ctheta -i stheta g5 ) D_eo
    wils.scale(psi_bfm,-kappa*kappa*kappa*kappa);

    

    wils.dslash(chi_bfm,rho_bfm,0,0); //EVEN_OUT DAG_NO
    wils.dslash(psi_bfm,sigma_bfm,0,1); //EVEN_OUT DAG_YES
  }
  
  //export results
  wils.cps_impexcbFermion((Float*)rho, rho_bfm, 0, 0); //even
  wils.cps_impexcbFermion((Float*)psi, psi_bfm, 0, 1); //odd
  wils.cps_impexcbFermion((Float*)sigma, sigma_bfm, 0, 0); //even

  return ;
}
#endif


CPS_END_NAMESPACE
#endif
