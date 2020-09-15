#ifndef _A2A_GRID_LANCZOS_H
#define _A2A_GRID_LANCZOS_H

#ifdef USE_GRID
#include<util/lattice/fgrid.h>

CPS_START_NAMESPACE

template<typename GridFermionField, typename GridGaugeField, typename GridDirac>
void gridLanczos(std::vector<Grid::RealD> &eval, std::vector<GridFermionField> &evec, const LancArg &lanc_arg,		 
		 GridDirac &Ddwf, const GridGaugeField& Umu,
		 Grid::GridCartesian *UGrid, Grid::GridRedBlackCartesian *UrbGrid,
		 Grid::GridCartesian *FGrid, Grid::GridRedBlackCartesian *FrbGrid){
  
  if(lanc_arg.N_true_get == 0){
    std::vector<Grid::RealD>().swap(eval); 	std::vector<GridFermionField>().swap(evec);      
    //eval.clear(); evec.clear();
    if(!UniqueID()) printf("gridLanczos skipping because N_true_get = 0\n");
    return;
  }

  assert(lanc_arg.precon);
  Grid::SchurDiagMooeeOperator<GridDirac, GridFermionField> HermOp(Ddwf);

    // int Nstop;   // Number of evecs checked for convergence
    // int Nk;      // Number of converged sought
    // int Np;      // Np -- Number of spare vecs in kryloc space
    // int Nm;      // Nm -- total number of vectors

  const int Nstop = lanc_arg.N_true_get;
  const int Nk = lanc_arg.N_get;
  const int Np = lanc_arg.N_use - lanc_arg.N_get;
  const int Nm = lanc_arg.N_use;
  const int MaxIt= lanc_arg.maxits;
  Grid::RealD resid = lanc_arg.stop_rsd;

  double lo = lanc_arg.ch_beta * lanc_arg.ch_beta;
  double hi = lanc_arg.ch_alpha * lanc_arg.ch_alpha;
  int ord = lanc_arg.ch_ord + 1; //different conventions

  if(!UniqueID()) printf("Chebyshev lo=%g hi=%g ord=%d\n",lo,hi,ord);
  
  Grid::Chebyshev<GridFermionField> Cheb(lo,hi,ord);
  Grid::ImplicitlyRestartedLanczos<GridFermionField> IRL(HermOp,Cheb,Nstop,Nk,Nm,resid,MaxIt);

  if(lanc_arg.lock) IRL.lock = 1;

  eval.resize(Nm);
  evec.reserve(Nm);
  evec.resize(Nm, FrbGrid);
 
  for(int i=0;i<Nm;i++){
    evec[i].checkerboard = Grid::Odd;
  }
  GridFermionField src(FrbGrid);
  
#ifndef MEMTEST_MODE
  if(!UniqueID()) printf("Starting Grid RNG seeding for Lanczos\n");
  double time = -dclock();
  
# if 0
  std::vector<int> seeds5({5,6,7,8});
  Grid::GridParallelRNG RNG5rb(FrbGrid);  RNG5rb.SeedFixedIntegers(seeds5);

  print_time("gridLanczos","RNG seeding",time+dclock());
  time = -dclock();

  if(!UniqueID()) printf("Initializing Gaussian src\n");
  gaussian(RNG5rb,src);
  src.checkerboard = Grid::Odd;
# else
  {
    CPSfermion5D<cps::ComplexD> tmp;
    tmp.setGaussianRandom();
    
    GridFermionField src_all(FGrid);
    tmp.exportGridField(src_all);
    pickCheckerboard(Grid::Odd,src,src_all);
  }
# endif

  print_time("gridLanczos","Gaussian src",time+dclock());
  time = -dclock();

  
  if(!UniqueID()) printf("Starting Lanczos algorithm\n");
  int Nconv;
  IRL.calc(eval,evec,
	   src,
	   Nconv);

  print_time("gridLanczos","Algorithm",time+dclock());
#endif
}

  

template<typename GridPolicies>
void gridLanczos(std::vector<Grid::RealD> &eval, std::vector<typename GridPolicies::GridFermionField> &evec, const LancArg &lanc_arg, typename GridPolicies::FgridGFclass &lattice){
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
  if(!UniqueID()) printf("Grid b=%g c=%g b+c=%g\n",mob_b,mob_c,mob_b+mob_c);

  typename GridDirac::ImplParams params;
  lattice.SetParams(params);

  GridDirac Ddwf(*Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,lanc_arg.mass,M5,mob_b,mob_c, params);

  gridLanczos(eval, evec, lanc_arg, Ddwf, *Umu, UGrid, UrbGrid, FGrid, FrbGrid);
}

template<typename GridPolicies>
void gridSinglePrecLanczos(std::vector<Grid::RealD> &eval, std::vector<typename GridPolicies::GridFermionFieldF> &evec, const LancArg &lanc_arg, typename GridPolicies::FgridGFclass &lattice,
			   Grid::GridCartesian *UGrid_f, Grid::GridRedBlackCartesian *UrbGrid_f,
			   Grid::GridCartesian *FGrid_f, Grid::GridRedBlackCartesian *FrbGrid_f
			   ){
  typedef typename GridPolicies::GridFermionFieldF GridFermionFieldF;
  typedef typename GridPolicies::FgridFclass FgridFclass;
  typedef typename GridPolicies::GridDiracF GridDiracF;
  
  Grid::QCD::LatticeGaugeFieldD *Umu = lattice.getUmu();

  Grid::QCD::LatticeGaugeFieldF Umu_f(UGrid_f);
  //Grid::precisionChange(Umu_f,*Umu);
  
  NullObject null_obj;
  lattice.BondCond();
  CPSfield<cps::ComplexD,4*9,FourDpolicy,OneFlavorPolicy> cps_gauge((cps::ComplexD*)lattice.GaugeField(),null_obj);
  cps_gauge.exportGridField(Umu_f);
  lattice.BondCond();

  double mob_b = lattice.get_mob_b();
  double mob_c = mob_b - 1.;   //b-c = 1
  double M5 = GJP.DwfHeight();
  if(!UniqueID()) printf("Grid b=%g c=%g b+c=%g\n",mob_b,mob_c,mob_b+mob_c);

  typename GridDiracF::ImplParams params;
  lattice.SetParams(params);

  GridDiracF Ddwf(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,lanc_arg.mass,M5,mob_b,mob_c, params);

  gridLanczos(eval, evec, lanc_arg, Ddwf, Umu_f, UGrid_f, UrbGrid_f, FGrid_f, FrbGrid_f);
}



CPS_END_NAMESPACE

#endif
#endif
