#ifndef _MAIN_CK_H
#define _MAIN_CK_H

#include <util/time_cps.h>
#include <alg/a2a/grid_lanczos.h>
#include <alg/ktopipi_jobparams.h>

//Useful functions for main programs
CPS_START_NAMESPACE

void ReadGaugeField(const MeasArg &meas_arg, bool double_latt = false){
  double time = -dclock();
  const char *cname = "main";
  const char *fname = "ReadGaugeField";

  GwilsonFdwf lat;
  std::ostringstream os;
  os << meas_arg.GaugeStem << '.' << meas_arg.TrajCur;
  std::string lat_file = os.str();

  ReadLatticeParallel rl;
  if(double_latt) rl.disableGparityReconstructUstarField();

  rl.read(lat,lat_file.c_str());
  if(!rl.good())ERR.General(cname,fname,"Failed read lattice %s",lat_file.c_str());

  time += dclock();
  print_time(cname,fname,time);
}

void ReadRngFile(const MeasArg &meas_arg, bool double_latt = false){
  double time = -dclock();
  const char *cname = "main";
  const char *fname = "ReadRngFile";

  std::ostringstream os;
  os << meas_arg.RNGStem << '.' << meas_arg.TrajCur;
  std::string rng_file = os.str();

  if(!LRG.Read(rng_file.c_str())) ERR.General(cname,fname,"Failed read rng file %s",rng_file.c_str());
  time += dclock();
  print_time(cname,fname,time);
}

#ifdef USE_BFM
template<typename mf_Float>
void setMass(bfm_evo<mf_Float> &dwf, const double &mass){
  dwf.mass = mass;
  dwf.GeneralisedFiveDimEnd(); // reinitialising since using a new mass
  dwf.GeneralisedFiveDimInit();
}

void setup_bfmargs(bfmarg &dwfa, int nthread, const BfmSolver &solver = HmCayleyTanh, const double mobius_scale = 1.){
  if(!UniqueID()) printf("Setting up bfmargs\n");

  omp_set_num_threads(nthread);

  dwfa.node_latt[0] = GJP.XnodeSites();
  dwfa.node_latt[1] = GJP.YnodeSites();
  dwfa.node_latt[2] = GJP.ZnodeSites();
  dwfa.node_latt[3] = GJP.TnodeSites();
  multi1d<int> ncoor(4);
  multi1d<int> procs(4);
  for(int i=0;i<4;i++){ ncoor[i] = GJP.NodeCoor(i); procs[i] = GJP.Nodes(i); }

  if(GJP.Gparity()){
    dwfa.gparity = 1;
    if(!UniqueID()) printf("G-parity directions: ");
    for(int d=0;d<3;d++)
      if(GJP.Bc(d) == BND_CND_GPARITY){ dwfa.gparity_dir[d] = 1; printf("%d ",d); }
      else dwfa.gparity_dir[d] = 0;
    for(int d=0;d<4;d++){
      dwfa.nodes[d] = procs[d];
      dwfa.ncoor[d] = ncoor[d];
    }
    if(!UniqueID()) printf("\n");
  }

  dwfa.verbose=1;
  dwfa.reproduce=0;
  bfmarg::Threads(nthread);
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

  dwfa.precon_5d = 1;
  if(solver == HmCayleyTanh){
    dwfa.precon_5d = 0; //mobius uses 4d preconditioning
    dwfa.mobius_scale = mobius_scale;
  }
  dwfa.Ls = GJP.SnodeSites();
  dwfa.solver = solver;
  dwfa.M5 = toDouble(GJP.DwfHeight());
  dwfa.mass = toDouble(0.01);
  dwfa.Csw = 0.0;
  dwfa.max_iter = 20000;
  dwfa.residual = 1e-08;
  if(!UniqueID()) printf("Finished setting up bfmargs\n");
}

 
void test_eigenvectors(BFM_Krylov::Lanczos_5d<double> &eig, bfm_evo<double> & dwf, bool singleprec_evecs){
  const int len = 24 * dwf.node_cbvol * (1 + dwf.gparity) * dwf.cbLs;
  omp_set_num_threads(bfmarg::threads);	

  Fermion_t bq_tmp = singleprec_evecs ? dwf.allocCompactFermion() : dwf.allocFermion(); 
  Fermion_t tmp1 = dwf.allocFermion();
  Fermion_t tmp2 = dwf.allocFermion();
  Fermion_t tmp3 = dwf.allocFermion();

  if(!UniqueID()) printf("Computing eigenvector residuals\n");

  for(int i=0;i<eig.get;i++){
    if(singleprec_evecs){ // eig->bq is in single precision
#pragma omp parallel for  //Bet I could reduce the threading overheads by parallelizing this entire method
      for(int j = 0; j < len; j++) {
	((double*)bq_tmp)[j] = ((float*)(eig.bq[i][1]))[j];
      }
    }else{
#pragma omp parallel
      {
	dwf.axpy(bq_tmp, eig.bq[i][1], eig.bq[i][1], 0.);
      }
    }
    
    double nrm_boss;
#pragma omp parallel
    {
      dwf.Mprec(bq_tmp,tmp1,tmp3, 0);
      dwf.Mprec(tmp1, tmp2, tmp3, 1); //tmp2 = M M^dag v
      
      //M M^dag v = lambda v
      dwf.set_zero(tmp1);	
      dwf.axpy(tmp3, bq_tmp, tmp1, eig.evals[i]); //tmp3 = lambda v
      
      double nrm = dwf.axpy_norm(tmp1, tmp2, tmp3, -1.); //tmp1 = tmp3 - tmp2
      if(dwf.isBoss()) nrm_boss = sqrt(nrm); //includes global sum
    }
    if(!UniqueID()) printf("%d %g\n",i,nrm_boss);
  }
  
  dwf.freeFermion(bq_tmp);
  dwf.freeFermion(tmp1);
  dwf.freeFermion(tmp2);
  dwf.freeFermion(tmp3);
}
#endif


#if defined(USE_GRID_LANCZOS)
template<typename GridPolicies>
void test_eigenvectors(const std::vector<typename GridPolicies::GridFermionField> &evec, const std::vector<Grid::RealD> &eval, const double mass, typename GridPolicies::FgridGFclass &lattice){
  typedef typename GridPolicies::GridFermionField GridFermionField;
  typedef typename GridPolicies::FgridFclass FgridFclass;
  typedef typename GridPolicies::GridDirac GridDirac;
  
  Grid::GridCartesian *UGrid = lattice.getUGrid();
  Grid::GridRedBlackCartesian *UrbGrid = lattice.getUrbGrid();
  Grid::GridCartesian *FGrid = lattice.getFGrid();
  Grid::GridRedBlackCartesian *FrbGrid = lattice.getFrbGrid();
  Grid::QCD::LatticeGaugeFieldD *Umu = lattice.getUmu();
  double mob_b = lattice.get_mob_b();
  double mob_c = mob_b - 1.;   //b-c = 1
  double M5 = GJP.DwfHeight();

  typename GridDirac::ImplParams params;
  lattice.SetParams(params);

  GridDirac Ddwf(*Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,mob_b,mob_c, params);
  Grid::SchurDiagMooeeOperator<GridDirac, GridFermionField> HermOp(Ddwf);
  
  GridFermionField tmp1(FrbGrid);
  GridFermionField tmp2(FrbGrid);
  GridFermionField tmp3(FrbGrid);
  
  for(int i=0;i<evec.size();i++){
    HermOp.Mpc(evec[i], tmp1);
    HermOp.MpcDag(tmp1, tmp2); //tmp2 = M^dag M v

    tmp3 = eval[i] * evec[i]; //tmp3 = lambda v

    double nrm = sqrt(axpy_norm(tmp1, -1., tmp2, tmp3)); //tmp1 = tmp3 - tmp2
    
    if(!UniqueID()) printf("%d %g\n",i,nrm);
  }

}
#endif


//Keep code clean by wrapping BFM or Grid underlay
struct LatticeSolvers{
#if !defined(USE_BFM_LANCZOS) && !defined(USE_BFM_A2A)

  LatticeSolvers(const JobParams &jp, const int nthreads){
    omp_set_num_threads(nthreads);
  }

#else
  bfm_evo<double> dwf_d;
  bfm_evo<float> dwf_f;
  bfmarg dwfa;

  LatticeSolvers(const JobParams &jp, const int nthreads){
    //Initialize both a double and single precision instance of BFM
    BfmSolver solver;
    switch(jp.solver){
    case BFM_DWF:
      solver = DWF; break;
    case BFM_HmCayleyTanh:
      solver = HmCayleyTanh; break;
    default:
      ERR.General("LatticeSolvers","constructor","Unknown solver\n");
    }
    setup_bfmargs(dwfa,nthreads,solver,jp.mobius_scale);
    dwf_d.init(dwfa);
    dwf_d.comm_end(); dwf_f.init(dwfa); dwf_f.comm_end(); dwf_d.comm_init();
  }

  ~LatticeSolvers(){
    dwf_d.end();
    dwf_f.end();
  }

#endif

};

template<typename LattType>
struct LatticeSetup{
# if defined(USE_BFM_LANCZOS) || defined(USE_BFM_A2A)
  static void importBFMlattice(Lattice *lat, LatticeSolvers &solvers){
    lat->BondCond(); //Apply the boundary conditions!
    Float* gauge = (Float*) lat->GaugeField();
    solvers.dwf_d.cps_importGauge(gauge);
    solvers.dwf_d.comm_end(); 
    solvers.dwf_f.comm_init(); solvers.dwf_f.cps_importGauge(gauge); solvers.dwf_f.comm_end(); 
    solvers.dwf_d.comm_init();
    lat->BondCond(); //Un-apply the boundary conditions! 
  }
#endif

  typedef LattType LatticeType;
  LatticeType *lat;
  
  //Grid or Grid/BFM mixed
#if defined(USE_GRID_LANCZOS) || defined(USE_GRID_A2A)

  LatticeSetup(const JobParams &jp, LatticeSolvers &solvers){
    assert(jp.solver == BFM_HmCayleyTanh);
    FgridParams grid_params; 
    grid_params.mobius_scale = jp.mobius_scale;
    lat = new LatticeType(grid_params);
    //lat->ImportGauge(); //lattice -> Grid  (applied APRD - signs internally then reverses)

    NullObject null_obj;
    lat->BondCond();
    CPSfield<cps::ComplexD,4*9,FourDpolicy,OneFlavorPolicy> cps_gauge((cps::ComplexD*)lat->GaugeField(),null_obj);
    cps_gauge.exportGridField(*lat->getUmu());
    lat->BondCond();
    
# if defined(USE_BFM_LANCZOS) || defined(USE_BFM_A2A)
    importBFMlattice(lat,solvers);
# endif
  }

#else
  //BFM only


  LatticeSetup(const JobParams &jp, LatticeSolvers &solvers){
    lat = new LatticeType; //doesn't actually matter
    importBFMlattice(lat,solvers);
  }

#endif

  LatticeType & getLattice(){ return *lat; }

  ~LatticeSetup(){
    delete lat;
  }

};

//Generates and stores evecs and evals
template<typename GridPolicies = void>
struct Lanczos{
#if defined(USE_GRID_LANCZOS)
  std::vector<Grid::RealD> eval; 
  std::vector<typename GridPolicies::GridFermionField> evec;
  std::vector<typename GridPolicies::GridFermionFieldF> evec_f;
  double mass;
  double resid;

  //For precision change
  Grid::GridCartesian *UGrid_f;
  Grid::GridRedBlackCartesian *UrbGrid_f;
  Grid::GridCartesian *FGrid_f;
  Grid::GridRedBlackCartesian *FrbGrid_f;
  
  Lanczos(): UGrid_f(NULL), UrbGrid_f(NULL), FGrid_f(NULL), FrbGrid_f(NULL){}
  
  void compute(const LancArg &lanc_arg, LatticeSolvers &solvers, typename GridPolicies::FgridGFclass &lat){
    mass = lanc_arg.mass;
    resid = lanc_arg.stop_rsd;
    
#ifdef A2A_LANCZOS_SINGLE
    //Make single precision Grids
    int Ls = GJP.Snodes()*GJP.SnodeSites();
    std::vector<int> nodes(4);
    std::vector<int> vol(4);
    for(int i=0;i<4;i++){
      vol[i]= GJP.NodeSites(i)*GJP.Nodes(i);;
      nodes[i]= GJP.Nodes(i);
    }
    std::vector<int> simd_layout = Grid::GridDefaultSimd(Grid::QCD::Nd,Grid::vComplexF::Nsimd());
    if(!UniqueID()) printf("Created single-prec Grids: nodes (%d,%d,%d,%d) vol (%d,%d,%d,%d) and SIMD layout (%d,%d,%d,%d)\n",nodes[0],nodes[1],nodes[2],nodes[3],vol[0],vol[1],vol[2],vol[3],simd_layout[0],simd_layout[1],simd_layout[2],simd_layout[3]);
    
    UGrid_f = Grid::QCD::SpaceTimeGrid::makeFourDimGrid(vol,simd_layout,nodes);
    UrbGrid_f = Grid::QCD::SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_f);
    FGrid_f = Grid::QCD::SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_f);
    FrbGrid_f = Grid::QCD::SpaceTimeGrid::makeFiveDimRedBlackGrid(GJP.SnodeSites()*GJP.Snodes(),UGrid_f); 

    gridSinglePrecLanczos<GridPolicies>(eval,evec_f,lanc_arg,lat,UGrid_f,UrbGrid_f,FGrid_f,FrbGrid_f);
#else    
    gridLanczos<GridPolicies>(eval,evec,lanc_arg,lat);
#  ifndef MEMTEST_MODE
    test_eigenvectors<GridPolicies>(evec,eval,lanc_arg.mass,lat);
#  endif

#endif
  }
  void toSingle(){
    typedef typename GridPolicies::GridFermionField GridFermionField;
    typedef typename GridPolicies::GridFermionFieldF GridFermionFieldF;
    
    //Make a single precision 5D checkerboarded Grid
    std::vector<int> nodes(4);
    std::vector<int> vol(4);
    for(int i=0;i<4;i++){
      vol[i]= GJP.NodeSites(i)*GJP.Nodes(i);;
      nodes[i]= GJP.Nodes(i);
    }
    UGrid_f = Grid::QCD::SpaceTimeGrid::makeFourDimGrid(vol,Grid::GridDefaultSimd(Grid::QCD::Nd,Grid::vComplexF::Nsimd()),nodes);
    FrbGrid_f = Grid::QCD::SpaceTimeGrid::makeFiveDimRedBlackGrid(GJP.SnodeSites()*GJP.Snodes(),UGrid_f);

    int nev = evec.size();
    for(int i=0;i<nev;i++){      
      GridFermionFieldF tmp_f(FrbGrid_f);
#ifndef MEMTEST_MODE
      precisionChange(tmp_f, evec.back());
#endif
      evec.pop_back();
      evec_f.push_back(std::move(tmp_f));
    }
    //These are in reverse order!
    std::reverse(evec_f.begin(), evec_f.end());
  }

  void freeEvecs(){
    std::vector<typename GridPolicies::GridFermionField>().swap(evec); //evec.clear();
    std::vector<typename GridPolicies::GridFermionFieldF>().swap(evec_f);
    if(UGrid_f != NULL) delete UGrid_f;
    if(UrbGrid_f != NULL) delete UrbGrid_f;
    if(FGrid_f != NULL) delete FGrid_f;
    if(FrbGrid_f != NULL) delete FrbGrid_f;
  }

#else
  
  BFM_Krylov::Lanczos_5d<double> *eig;

  Lanczos(): eig(NULL){}
  
  void compute(const LancArg &lanc_arg, LatticeSolvers &solvers, Lattice &lat){
    eig = new BFM_Krylov::Lanczos_5d<double>(solvers.dwf_d,const_cast<LancArg&>(lanc_arg)); //sets up the mass of dwf_d correctly
    eig->Run();

    solvers.dwf_f.mass = solvers.dwf_d.mass; //keep the single-prec solver in sync
    solvers.dwf_f.GeneralisedFiveDimEnd(); // reinitialising since using a new mass
    solvers.dwf_f.GeneralisedFiveDimInit();

    test_eigenvectors(*eig,solvers.dwf_d,false);
  }

  void toSingle(){
    eig->toSingle(); 
    //Test the single-prec converted eigenvectors to make sure we haven't dropped too much precision
    test_eigenvectors(*eig,eig->dop,true);
  }

  void freeEvecs(){
    eig->free_bq();
  }

  ~Lanczos(){
    if(eig != NULL)
      delete eig;
  }
#endif

};


template<typename mf_Policies, typename LanczosPolicies>
struct computeA2Avectors{
  static void compute(A2AvectorV<mf_Policies> &V, A2AvectorW<mf_Policies> &W, bool mixed_solve, bool evecs_single_prec, Lattice &lat, Lanczos<LanczosPolicies> &eig, LatticeSolvers &solvers){
#ifdef USE_BFM_LANCZOS
    W.computeVW(V, lat, *eig.eig, evecs_single_prec, solvers.dwf_d, mixed_solve ? & solvers.dwf_f : NULL);
#else
    if(evecs_single_prec)
      W.computeVW(V, lat, eig.evec_f, eig.eval, eig.mass, eig.resid, 10000);
    else
      W.computeVW(V, lat, eig.evec, eig.eval, eig.mass, eig.resid, 10000);
#endif
  }
};

template<typename ComplexType>
void setupFieldParams(cps::NullObject &n){}

#ifdef USE_GRID
template<typename ComplexType>
void setupFieldParams(typename FourDSIMDPolicy::ParamType &p){
  int nsimd = ComplexType::Nsimd();
  FourDSIMDPolicy::SIMDdefaultLayout(p,nsimd,2); //only divide over spatial directions
  
  printf("4D field params: Nsimd = %d, SIMD dimensions:\n", nsimd);
  for(int i=0;i<4;i++)
    printf("%d ", p[i]);
  printf("\n");
}
template<typename ComplexType>
void setupFieldParams(typename ThreeDSIMDPolicy::ParamType &p){
  int nsimd = ComplexType::Nsimd();
  ThreeDSIMDPolicy::SIMDdefaultLayout(p,nsimd);
  
  printf("3D field params: Nsimd = %d, SIMD dimensions:\n", nsimd);
  for(int i=0;i<3;i++)
    printf("%d ", p[i]);
  printf("\n");
}
#endif


CPS_END_NAMESPACE


#endif
