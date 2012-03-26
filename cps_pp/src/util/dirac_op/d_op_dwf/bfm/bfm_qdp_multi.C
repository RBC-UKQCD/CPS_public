#if 0
//#ifdef USE_BFM
#include <qdp.h>
#include <bfm.h>
#include <bfm_qdp.h>

template<class Float>
void * bfm_thread_cg_routine_CJ(void * arg)
{
  bfm_qdp<Float> *dwf = (bfm_qdp<Float> *)arg;

  int me = dwf->thread_barrier();

  uint64_t cfg_pf_usr = L1P_CFG_PF_USR_dfetch_depth(2)
    | L1P_CFG_PF_USR_dfetch_max_footprint(3)   
    | L1P_CFG_PF_USR_ifetch_depth(0)       
    | L1P_CFG_PF_USR_ifetch_max_footprint(1)   
    | L1P_CFG_PF_USR_pf_stream_est_on_dcbt 
    | L1P_CFG_PF_USR_pf_stream_establish_enable
    | L1P_CFG_PF_USR_pf_adaptive_throttle(0xF) ;

  out64_sync((uint64_t *)L1P_CFG_PF_USR,cfg_pf_usr);
  out64_sync((uint64_t *)L1P_CFG_CLK_GATE,0x7F);

  int CGNE_MdagM(Fermion_t sol_guess[2], Fermion_t source[2]);//Not implemented
  int CGNE(Fermion_t sol_guess[2], Fermion_t source[2]);// deprecated
  int CGNE_prec(Fermion_t sol_guess, Fermion_t source); // Internal

  if ( dwf->inv_type == CG_PREC_MDAGM ) {
    dwf->iter = dwf->CGNE_prec_MdagM(dwf->qdp_chi_h[1],dwf->qdp_psi_h[1]);
  } else if ( dwf->inv_type == CG_PREC_M ) {
    dwf->iter = dwf->CGNE_M(dwf->qdp_chi_h,dwf->qdp_psi_h);
  } else if ( dwf->inv_type == CG_PREC_MDAG ) {
    dwf->iter = dwf->CGNE_Mdag(dwf->qdp_chi_h,dwf->qdp_psi_h);
  } else if ( dwf->inv_type == CG_PREC_MDAGM_MULTI ) {
    if (!me) printf("dwf->iter = dwf->CGNE_prec_MdagM_multi_shift_CJ\n");
    if (!me) printf("dwf->verbose=%d me=%d\n",dwf->verbose,me);
    dwf->iter = dwf->CGNE_prec_MdagM_multi_shift(dwf->qdp_chi_multi_h,
						 dwf->qdp_psi_h[1],// odd parity
						 dwf->shifts,
						 dwf->alpha,
						 dwf->nshift,
						 dwf->mresidual,
						 dwf->single); // sum single result
    if (!me) printf("dwf->iter = dwf->CGNE_prec_MdagM_multi_shift_CJ done \n");
  } else if ( dwf->inv_type == M_UNPREC ) {
    dwf->Munprec(dwf->qdp_psi_h,dwf->qdp_chi_h,dwf->mtmp,DaggerNo);
  } else if ( dwf->inv_type == MDAG_UNPREC ) {
    dwf->Munprec(dwf->qdp_psi_h,dwf->qdp_chi_h,dwf->mtmp,DaggerYes);
  } else if ( (dwf->inv_type == CG_UNPREC_M) 
	   || (dwf->inv_type == CG_UNPREC_MDAG)
	   || (dwf->inv_type == CG_UNPREC_MDAGM) ){
    QDP_error_exit("Not implmented UNPREC cg as threaded yet");
  }

  fflush(stdout);
  return (void *)dwf;
}

template<class Float> void   bfm_spawn_cg_CJ (bfm_qdp<Float> &dwf);

template<class Float>
void bfm_spawn_cg_CJ (bfm_qdp<Float> &dwf)
{
    int THREADS= dwf.threads;
     bfmarg::Threads(16);
     bfmarg::Reproduce(0);
    bfmarg::ReproduceChecksum(0);
     bfmarg::ReproduceMasterCheck(0);
    bfmarg::Verbose(1);
    printf("THREADS=%d\n",THREADS);
    pthread_t handles[THREADS];
    for(int t =0; t<THREADS; t++ ) {
      pthread_create(&handles[t],(pthread_attr_t *)NULL,bfm_thread_cg_routine_CJ<Float>,(void *)&dwf);
    }
    for(int t =0; t<THREADS; t++ ) {
      pthread_join(handles[t],NULL);
    }
}


template <class Float>
int dwf_CG_precMdagM_oo_multi(multi2d<LatticeFermion> &sol_qdp,
			      multi1d<LatticeFermion> &src_qdp,
			      multi1d<LatticeColorMatrix> &U,
                              double    shifts[],
                              double    alpha[],
                              int       nshift,
                              double mresidual_p[],
                              int single,
                              int Ls,
                              Real mass,
                              Real M5,
                              Real residual, int max_iter)
{
  // Set up BAGEL object
  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  multi1d<int> procs = QDP::Layout::logicalSize();

  bfmarg dwfa;

  //Physics parameters
  dwfa.Ls           = Ls;
  dwfa.M5           = toDouble(M5);
  dwfa.mass         = toDouble(mass);
  dwfa.precon_5d    = 1;
  dwfa.list_engine  = 0;
  dwfa.list_length  = 8192;
  dwfa.max_iter     = max_iter;
  dwfa.residual     = toDouble(residual);

  //Geometry
  dwfa.node_latt[0] = lx;
  dwfa.node_latt[1] = ly;
  dwfa.node_latt[2] = lz;
  dwfa.node_latt[3] = lt;
  for(int mu=0;mu<4;mu++){
    if (procs[mu]>1) dwfa.local_comm[mu] = 0;
    else             dwfa.local_comm[mu] = 1;
  }

//        bfmarg::Threads(16);
 //       bfmarg::Reproduce(0);
  //      bfmarg::ReproduceChecksum(0);
   //     bfmarg::ReproduceMasterCheck(0);
    //    bfmarg::Verbose(1);
  dwfa.verbose=1;
  dwfa.threads=64;


  QDPIO::cout << "bfm_qdp:: Initialising BAGEL-2 solver "<<endl;

  bfm_qdp<Float> dwf; dwf.init(dwfa);

  dwf.importGauge(U);

  Fermion_t src;
  double *mresidual = new double[nshift];
  for(int i =0;i<nshift;i++){
    mresidual[i] = mresidual_p[i];
    printf("mresidual[%d]=%e\n",i,mresidual[i]);
  }
  src = dwf.allocFermion();
  dwf.importFermion(src_qdp,src,1);
  dwf.qdp_psi_h[0] = src;
  dwf.qdp_psi_h[1] = src;
  dwf.shifts=shifts;
  dwf.alpha =alpha;
  dwf.nshift=nshift;
  dwf.mresidual=mresidual;
  dwf.single=0;

  Fermion_t single_sols[nshift];
  dwf.qdp_chi_multi_h=single_sols;
//  if ( single ) {
    for(int s=0;s<nshift;s++) {
      single_sols[s]=dwf.allocFermion();
      if (single_sols[s]==NULL){
        printf("bfm_qdp: bad allocate\n");
        exit(-1);
      }
    }
//  }


  int iter=0;
  dwf.inv_type = CG_PREC_MDAGM_MULTI;
//  if ( dwfa.threads == 1 ) {
//    bfm_thread_cg_routine<Float>((void *)&dwf);
//  } else { 
//    bfm_spawn_cg_CJ(dwf);
    bfm_spawn_cg(dwf);
//  }
    printf("bfm_spawn_cg(dwf) done\n");

  if ( single ) {
    Fermion_t sol_guess;
    sol_guess = dwf.allocFermion();
    dwf.master_fill(sol_guess,0);
    for(int s=0;s<nshift;s++){
      dwf.axpy(sol_guess,single_sols[s],sol_guess,alpha[s]);
      dwf.freeFermion(single_sols[s]);
    }
    multi1d<LatticeFermion> sol_qdp_tmp = sol_qdp[0];
    dwf.exportFermion(sol_qdp_tmp,sol_guess,1);
    dwf.freeFermion(sol_guess);
  } else {
    for(int s=0;s<nshift;s++){
    multi1d<LatticeFermion> sol_qdp_tmp = sol_qdp[s];
    dwf.exportFermion(sol_qdp_tmp,single_sols[s],1);
    dwf.freeFermion(single_sols[s]);
    }
  }

  delete[] mresidual;

 

#if 0
  /*******************************************************/
  /* verify the solution                                 */
  /*******************************************************/
  QDPIO::cout << "bfm_qdp:: Verifying solution" <<endl;

  /*
   * Compute the residual according to CHROMA
   */
   multi1d<LatticeFermion> regress(Ls);
   multi1d<int> bcs(Nd);
   bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;

   Handle< FermBC<T,G,G> > fbc(new SimpleFermBC< T, G, G >(bcs));
   Handle<CreateFermState<T,G,G> > cfs( new CreateSimpleFermState<T,G,G>(fbc));

   UnprecDWFermActArray  S_f(cfs, M5, mass, Ls);
   Handle< FermState<T,G,G> > fs( S_f.createState(U) );

   Handle< UnprecLinearOperatorArray<T,G,G> > M(S_f.unprecLinOp(fs,mass));
  
   // Check the result
   PlusMinus pm = PLUS;

   (*M)(regress,sol,pm);

   regress = regress - src;
   QDPIO::cout << "bfm_qdp:: QDP regression check :  |M sol - src| = " << norm2(regress) << endl;
   if ( toDouble( norm2(regress) / norm2(src) ) > 1.0e-5 ) { 
     QDPIO::cout << "bfm_qdp:: QDP regression check : This is worryingly large - PANIC"<< endl;
     exit(-1);
   }
#endif
    printf("dwf_CG_precMdagM_oo_multi done\n");
   return iter;

}

template int dwf_CG_precMdagM_oo_multi<float>( multi2d<LatticeFermion> &sol,
			      multi1d<LatticeFermion> &src,
			      multi1d<LatticeColorMatrix> &U,
                              double    shifts[],
                              double    alpha[],
                              int       nshift,
                              double mresidual[],
                              int single,
                              int Ls,
                              Real mass,
                              Real M5,
Real residual, int max_iter);
 
template int dwf_CG_precMdagM_oo_multi<double>( multi2d<LatticeFermion> &sol,
			      multi1d<LatticeFermion> &src,
			      multi1d<LatticeColorMatrix> &U,
                              double    shifts[],
                              double    alpha[],
                              int       nshift,
                              double mresidual[],
                              int single,
                              int Ls,
                              Real mass,
                              Real M5,
                              Real residual, int max_iter);
#endif
