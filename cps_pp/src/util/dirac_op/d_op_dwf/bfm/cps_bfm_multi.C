#include <config.h>
#include <stdio.h>
#include <stdarg.h>
#include <omp.h>
#ifdef USE_BFM

#include <chroma.h>
#include <bfm.h>
#include <bfm_qdp.h>
#ifdef UNIFORM_SEED_TESTING
#include "majority_vote.h"
#endif
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
//#include <qdp.h>
#include <util/gjp.h>
#include <comms/sysfunc_cps.h>
#include <util/error.h>
#include <util/verbose.h>
#include <util/smalloc.h>
#include <util/time_cps.h>
#include <util/dirac_op/d_op_dwf.h>
#ifdef USE_BFM_MINV

//static int  Printf(char *format,...){}
#define Printf if ( QMP_is_primary_node() ) printf
//#define Printf printf


using namespace Chroma;

//typedef LatticeFermion T;
//typedef multi1d<LatticeFermion> T5;
//typedef multi1d<LatticeColorMatrix> U;
//
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
int even, int odd, 
int Ls=0, double fac_t=1.);


CPS_START_NAMESPACE
static void print_vec(void *p, char *name){
 Float *tmp_p = (Float *)p;
// if (!UniqueID()) printf("%s(%p) = %g\n",name, p, *tmp_p);
}
int cps_qdp_init(int *argc, char ***argv);

int cps_qdp_finalize();

int DiracOpDwf::MInvCG(Vector **out, Vector *in, Float in_norm, Float *shift, 
             int Nshift, int isz, Float *RsdCG,
             MultiShiftSolveType type, Float *alpha)
{
        char *fname="MInvCG(V**,V*,F,F*,i,i,F*,t,F*)";
        VRB.Func(cname,fname);
  unsigned int f_size_cb =  GJP.VolNodeSites() * lat.FsiteSize() / 2;

	VRB.Result(cname,fname,"cps_qdp_init(GJP.argc_p(), GJP.argv_p())");
	cps_qdp_init(GJP.argc_p(), GJP.argv_p());
	int iter=0;

	double minv_time = -dclock();
	Float src_norm_sq = in->NormSqNode(f_size_cb);
	DiracOpGlbSum(&src_norm_sq);



{
//  Chroma::initialize(&argc,&argv);

	if (( type != SINGLE ) && ( type != MULTI )){
    ERR.General(cname,fname,"( type != SINGLE ) && ( type != MULTI )\n");
	}


  /********************************************************
   * Command line parsing
   ********************************************************
   */
#if 0
#undef COMMANDLINE
#ifdef COMMANDLINE
  if ( argc != 6 ) { 
   Printf("Usage: %s lx ly lz lt Ls\n All must be even\n",argv[0]);
   Printf("argc is %d\n",argc);
   for ( int i=0;i<argc;i++){
      Printf("%d %s\n",i,argv[i]);
   }
      exit
(-1);

  }
#endif
  /********************************************************
   * Setup QDP
   ********************************************************
   */
  multi1d<int> nrow(Nd);
#ifdef COMMANDLINE
  nrow[0] = atoi(argv[1]);
  nrow[1] = atoi(argv[2]);
  nrow[2] = atoi(argv[3]);
  nrow[3] = atoi(argv[4]);
  int Ls = atoi(argv[5]);
#else
  nrow[0] = 4;
  nrow[1] = 4;
  nrow[2] = 4;
  nrow[3] = 4;
  int Ls  = 4;
#endif
#endif

//  Layout::setLattSize(nrow);
//  Layout::create();

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  /********************************************************
   * Setup DWF operator
   ********************************************************
   */
  bfmarg dwfa;
  omp_set_num_threads(64);
#ifdef UNIFORM_SEED_TESTING
  majorityVote  dwf;
#else
  bfm_qdp<double>  dwf;
#endif
  dwfa.solver = DWF;
  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;
  dwfa.verbose=1;
  dwfa.reproduce=0;
  bfmarg::Threads(64);
  bfmarg::Reproduce(0);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(1);
  bfmarg::onepluskappanorm = 0;

  multi1d<int> procs = QDP::Layout::logicalSize();
  Printf("%d dim machine\n\t", procs.size());
  for(int mu=0;mu<4;mu++){
    Printf("%d ", procs[mu]);
    if ( procs[mu]>1 ) {
      dwfa.local_comm[mu] = 0;
    } else { 
      dwfa.local_comm[mu] = 1;
    }
  }
  Printf("\nLocal comm = ");
  for(int mu=0;mu<4;mu++){
    Printf("%d ", dwfa.local_comm[mu]);
  }
  Printf("\n");
  
  multi1d<int> ncoor = QDP::Layout::nodeCoord();

  int Ls = GJP.NodeSites(4);
  double M5 = GJP.DwfHeight();
  double mq= dirac_arg->mass;
  Printf("Ls=%d M5=%g mq=%g\n",Ls,M5,mq);

  dwfa.precon_5d = 1;
  dwfa.Ls   = Ls;
  dwfa.M5   = toDouble(M5);
  dwfa.mass = toDouble(mq);
  dwfa.Csw  = 0.0;
  dwfa.max_iter = dirac_arg->max_num_iter;
  dwfa.residual = 1.e-8;
  Printf("Initialising bfm operator\n");
  dwf.init(dwfa);

  /********************************************************
   * Gaussian gauge field
   ********************************************************
   */
#if 0

  multi1d<LatticeColorMatrix>  u(Nd);
  HotSt(u);
  multi1d<LatticeColorMatrixF> uf(Nd);
  for(int m=0; m < Nd; ++m){
    uf[m] = u [m];
    u [m] = uf[m];
  }
  Printf("Setup gauge field\n");
#endif
  /********************************************************
   * Gaussian source and result vectors
   ********************************************************
   */
  Fermion_t test= dwf.allocFermion();
  Printf("Filling with zeroes field\n");
//  fflush(stdout);
  dwf.master_fill(test,0.0);
  print_vec(test,"test");

#define NMULTI (3)
   const int MaxShift = 20;
  if(Nshift>MaxShift)
    ERR.General(cname,fname,"Nshift(%d)>MaxShift(%d)\n",Nshift,MaxShift);
   Fermion_t sol_bfm[MaxShift];
   double masses[MaxShift];
   double alpha_bfm [MaxShift];
   double residuals[MaxShift];
   
  double fac = (5-M5)*(5-M5);
   Printf("type=%d shift=%p alpha=%p RsdCG=%p\n",type,shift,alpha,RsdCG);
   for(int i =0; i<Nshift;i++){
        sol_bfm[i] = dwf.allocFermion();
        dwf.master_fill(sol_bfm[i],0.0);
        masses[i] = shift[i]*fac;
        Printf("sol_bfm[%d]=%p\n",i,sol_bfm[i]);
        alpha_bfm[i] = 1.;
        residuals[i]=RsdCG[i];
#if 0
	if (i==(Nshift-1)) residuals[i]*=0.1;
	if ((type==MULTI) || (alpha==NULL)){
        alpha_bfm[i] = 1.;
        residuals[i]=RsdCG[i];
	} else {
        alpha_bfm[i] = alpha[i];
        residuals[i]=RsdCG[i]*alpha[i];
	}
#endif

   Printf("%d: shift=%g alpha_bfm=%g RsdCG=%g\n",i,shift[i],alpha_bfm[i],RsdCG[i]);
   }


  /********************************************************
   * Import gauge field to BAGEL
   ********************************************************
   */

  Float *gauge = (Float*) lat.GaugeField();
{
  multi1d<LatticeColorMatrix> U(Nd);
  importGauge(lat,U,gauge,1);
  dwf.importGauge(U);
}

  Float *src_cps = (Float*) in;
  Fermion_t src_bfm = dwf.allocFermion();
  dwf.master_fill(src_bfm,0.0);
{
#if 0
  multi1d<LatticeFermion> source(Ls);
  for(int s=0;s<Ls;s++) gaussian(source[s]);
#endif
  multi1d<LatticeFermion> src_qdp(Ls);
  impexFermion(0,lat,src_qdp,src_cps,0,1);
  dwf.importFermion(src_qdp,src_bfm,1);
}

#if 0
//  Printf("Calling half cb inverter\n"); fflush(stdout);
  dwf.inv_type=CG_PREC_MDAGM;
#endif

  Printf("Calling multi-shift inverter\n");
//fflush(stdout);
  dwf.inv_type=CG_PREC_MDAGM_MULTI;
  dwf.qdp_psi_h[0]=src_bfm;
  dwf.qdp_psi_h[1]=src_bfm;
  dwf.qdp_chi_multi_h=sol_bfm;
  dwf.qdp_chi_h[0]=sol_bfm[0];
  dwf.qdp_chi_h[1]=sol_bfm[0];
  dwf.shifts=masses;
  dwf.alpha =alpha_bfm;
  dwf.nshift=Nshift;
  dwf.mresidual=residuals;
  dwf.single=0;
double bfm_time = -dclock();
  bfm_spawn_cg(dwf);
  bfm_time  += dclock();
  print_flops(fname,"bfm",0,bfm_time);

  print_vec(sol_bfm[0],"sol_bfm[0]");
  bfm_time = -dclock();
{
	multi1d<LatticeFermion> sol_qdp(Ls);
     int n_sol=Nshift;
//    if (dwf.single) n_sol=1;
	VRB.Result(cname,fname,"type=%d\n",type);

	Float *sol_cps =NULL;
	sol_cps = (Float *)fmalloc(cname,fname,"sol_cps",f_size_cb*2*sizeof(Float));
	memset(sol_cps,0,f_size_cb*2*sizeof(Float));
	Printf("sol_cps=%p\n",sol_cps);
	Vector *sol = (Vector *) sol_cps;
	Vector *src = (Vector *) in;
	Vector *mmp = (Vector *) smalloc(cname,fname,"mmp",f_size_cb * sizeof(Float));
	Vector *res = (Vector *) smalloc(cname,fname,"res",f_size_cb * sizeof(Float));

	double M5 = GJP.DwfHeight();
	double fac = (5-M5)*(5-M5);
	for(int i=0;i<n_sol;i++){
		dwf.exportFermion(sol_qdp,sol_bfm[i],1);
  		print_vec(sol_bfm[i],"sol_bfm");
		impexFermion(1,lat,sol_qdp,sol_cps,0,1,0,fac);
  		print_vec(sol_cps,"sol_cps");
		Float coef = 1.;
// calculating true residual. Agrees with BFM, so turned off.
		if (0) {
			Float shift_local=shift[i];
//			InvCgShift(sol,src,0.,NULL,&shift_local);
			MatPcDagMatPc(mmp, sol);
			mmp -> FTimesV1PlusV2(shift[i],sol,mmp, f_size_cb);
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
//			VRB.Result(cname,fname,"out[0]=%g\n",*tmp_p);
		} else {
			Float *tmp_p = (Float *)out[i];
#pragma omp parallel for default(shared)
			for(long j=0;j<f_size_cb;j++)
			tmp_p[j] =coef*sol_cps[j];
			
		}
		Printf("sol_bfm[%d]=%p freed\n",i,sol_bfm[i]);
		dwf.freeFermion(sol_bfm[i]);
	}
	ffree(cname,fname,"sol_cps",sol_cps);
	ffree(cname,fname,"mmp",mmp);
	ffree(cname,fname,"res",res);
}
  bfm_time  += dclock();
  print_flops(fname,"convert",0,bfm_time);
  
  //  dwf.CGNE_prec_MdagM_multi_shift(psi,
  //				  src,
  //				  masses,
  //				  alpha,
  //				  NMULTI,
  //				  residuals,
  //				  0);

  Printf("src_bfm=%p freed\n",src_bfm);
  dwf.freeFermion(src_bfm);
  Printf("test=%p freed\n",test);
  dwf.freeFermion(test);
  Printf("Done\n"); 
  dwf.end();
  iter= dwf.iter;

}
	minv_time += dclock();
	print_flops(cname,fname,0,minv_time);

	cps_qdp_finalize();
	return iter;

}
CPS_END_NAMESPACE
#endif
#endif
