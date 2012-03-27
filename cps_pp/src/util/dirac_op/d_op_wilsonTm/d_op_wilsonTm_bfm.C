#include <config.h>
#include <stdio.h>
#include <stdarg.h>
#ifdef USE_OMP
#include <omp.h>
#endif

#ifdef USE_BFM_TM
#include <chroma.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
//#include <qdp.h>
#include <util/gjp.h>
#include <comms/sysfunc_cps.h>
#include <util/error.h>
#include <util/verbose.h>
#include <util/smalloc.h>
#include <util/time_cps.h>
#include <util/dirac_op/d_op_dwf.h>

static int  Printf(char *format,...){}
//#define Printf if ( QMP_is_primary_node() ) printf
//#define Printf printf


using namespace Chroma;

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

   char *fname="InvCg(V*,V*,F,F*)";
   double cg_time = -dclock();
  unsigned int f_size_cb =  GJP.VolNodeSites() * lat.FsiteSize() / (lat.FchkbEvl()+1);
//  if(src_norm_sq == 0){
  src_norm_sq = in->NormSqNode(f_size_cb);
  DiracOpGlbSum(&src_norm_sq);
//  }
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
  //  if (!qdp_initted)
  //    ERR.General("",fname,"QDP not initialized!(%d)",0);
  
  
  //  int Ls = GJP.SnodeSites();
    double M5 = GJP.DwfHeight();
  //  int Nd = 4;
  
    Printf("src[0]=%g\n",*src_cps);
    int lx = QDP::Layout::subgridLattSize()[0];
    int ly = QDP::Layout::subgridLattSize()[1];
    int lz = QDP::Layout::subgridLattSize()[2];
    int lt = QDP::Layout::subgridLattSize()[3];
    multi1d<LatticeFermion> src_qdp(1);
//    multi1d<LatticeFermion> sol_qdp(1);
  
    /********************************************************
     * Setup DWF operator
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
    bfmarg::Threads(64);
    bfmarg::Reproduce(0);
    bfmarg::ReproduceChecksum(0);
    bfmarg::ReproduceMasterCheck(0);
    bfmarg::Verbose(0);

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
  
  //  Real M5(1.8);
  //  Real mq(0.1);
  
    omp_set_num_threads(64);
    Float *gauge = (Float*) lat.GaugeField();
    wilsa.precon_5d = 0;
    wilsa.Ls   = 1;
    wilsa.solver = WilsonTM;
    wilsa.M5   = 0.0;
    wilsa.mass = toDouble(mass);
    wilsa.twistedmass = toDouble(epsilon);
    wilsa.Csw  = 0.0;
    wilsa.list_engine=0;
    wilsa.list_length=0;
    wilsa.max_iter = max_iter;
    wilsa.residual = toDouble(residual);
  //OK
    for(int i = 0;i<1;i++){
      Printf("Initialising bfm operator\n");
      wils.init(wilsa);
    
  //OK
      {
        multi1d<LatticeColorMatrix> U(Nd);
        importGauge(lat,U,gauge,1);
        wils.importGauge(U);
      }
    
      Printf("Setup gauge field\n");
      /********************************************************
       * Gaussian source and result vectors
       ********************************************************
       */
  //  multi1d<LatticeFermion> source(Ls);
  //  for(int s=0;s<Ls;s++) gaussian(source[s]);
    
      Fermion_t src_bfm = wils.allocFermion();
  //  wils.master_fill(src_bfm,0.0);
      Printf("src_bfm=%p \n",src_bfm);
      Fermion_t psi[1];
      psi[0] = wils.allocFermion();
  //  wils.master_fill(psi[0],0.0);
    
      Printf("psi[0]=%p \n",psi[0]);
      double fac= 0.25/(kappa*kappa);
      for(int i = 0;i<1;i++){
      
      
        /********************************************************
         * Import gauge field to BAGEL
         ********************************************************
         */
      
        impexFermion(0,lat,src_qdp,sol_cps,0,1,1,1./fac);
        wils.importFermion(src_qdp[0],psi[0],1);
        impexFermion(0,lat,src_qdp,src_cps,0,1,1);
        wils.importFermion(src_qdp[0],src_bfm,1);
        Float *tmp_p = (Float *)src_bfm;
        Printf("src_bfm[0]=%g\n",*tmp_p);
      
        Printf("Calling half cb inverter\n"); fflush(stdout);
        wils.inv_type=CG_PREC_MDAGM;
        wils.qdp_chi_h[0]=psi[0];
        wils.qdp_chi_h[1]=psi[0];
        wils.qdp_psi_h[0]=src_bfm;
        wils.qdp_psi_h[1]=src_bfm;
        bfm_spawn_cg(wils);
        tmp_p = (Float *)psi[0];
        Printf("psi[0]=%g\n",*(tmp_p));
        wils.exportFermion(src_qdp[0],psi[0],1);
      }
      
      impexFermion(1,lat,src_qdp,sol_cps,0,1,1,fac);
      Printf("sol[0]=%g\n",*sol_cps);
      
      Printf("src_bfm=%p freed\n",src_bfm);
      wils.freeFermion(src_bfm);
      Printf("psi[0]=%p freed\n",psi[0]);
      wils.freeFermion(psi[0]);
      Printf("Done\n"); 
//OK
      wils.end();
    }
//NOT OK
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
}

CPS_END_NAMESPACE
#endif
