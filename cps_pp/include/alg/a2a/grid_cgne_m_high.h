#ifndef _GRID_CGNE_M_HIGH_H
#define _GRID_CGNE_M_HIGH_H

#ifdef USE_GRID
#include<util/lattice/fgrid.h>
CPS_START_NAMESPACE

//nLowMode is the number of modes we actually use to deflate. This must be <= evals.size(). The full set of computed eigenvectors is used to improve the guess.
template<typename GridPolicies>
inline void Grid_CGNE_M_high(typename GridPolicies::GridFermionField &solution, const typename GridPolicies::GridFermionField &source, double resid, int max_iters,
			     EvecInterface<GridPolicies> &evecs, int nLowMode, 
			     typename GridPolicies::FgridFclass &latg, typename GridPolicies::GridDirac &Ddwf, Grid::GridCartesian *FGrid, Grid::GridRedBlackCartesian *FrbGrid){
  typedef typename GridPolicies::GridFermionField GridFermionField;
  typedef typename GridPolicies::FgridFclass FgridFclass;
  typedef typename GridPolicies::GridDirac GridDirac;
  
  double f = norm2(source);
  if (!UniqueID()) printf("Grid_CGNE_M_high: Source norm is %le\n",f);
  f = norm2(solution);
  if (!UniqueID()) printf("Grid_CGNE_M_high: Guess norm is %le\n",f);

  Grid::SchurDiagMooeeOperator<GridDirac, GridFermionField> linop(Ddwf);

  GridFermionField tmp_cb1(FrbGrid);
  GridFermionField tmp_cb2(FrbGrid);
  GridFermionField tmp_cb3(FrbGrid);
  GridFermionField tmp_cb4(FrbGrid);

  GridFermionField tmp_full(FGrid);

  // src_o = Mprecdag * (source_o - Moe MeeInv source_e)  , cf Daiqian's thesis page 60
  GridFermionField src_o(FrbGrid);

  pickCheckerboard(Grid::Even,tmp_cb1,source);  //tmp_cb1 = source_e
  pickCheckerboard(Grid::Odd,tmp_cb2,source);   //tmp_cb2 = source_o

  Ddwf.MooeeInv(tmp_cb1,tmp_cb3);
  Ddwf.Meooe     (tmp_cb3,tmp_cb4); //tmp_cb4 = Moe MeeInv source_e       (tmp_cb3 free)
  axpy    (tmp_cb3,-1.0,tmp_cb4, tmp_cb2); //tmp_cb3 = (source_o - Moe MeeInv source_e)    (tmp_cb4 free)
  linop.MpcDag(tmp_cb3, src_o); //src_o = Mprecdag * (source_o - Moe MeeInv source_e)    (tmp_cb3, tmp_cb4 free)

  //Compute low-mode projection and CG guess
  int Nev = evecs.nEvecs();

  GridFermionField lsol_full(FrbGrid); //full low-mode part (all evecs)
  lsol_full = Grid::zero;

  GridFermionField lsol_defl(FrbGrid); //low-mode part for subset of evecs with index < nLowMode
  lsol_defl = Grid::zero;
  lsol_defl.checkerboard = Grid::Odd;
  
  GridFermionField sol_o(FrbGrid); //CG solution
  sol_o = Grid::zero;

  if(Nev < nLowMode)
    ERR.General("","Grid_CGNE_M_High","Number of low eigen modes to do deflation is smaller than number of low modes to be substracted!\n");

  if(Nev > 0){
    if (!UniqueID()) printf("Grid_CGNE_M_High: deflating with %d evecs\n",Nev);

    for(int n = 0; n < Nev; n++){
      double eval = evecs.getEvec(tmp_cb1,n);
      Grid::ComplexD cn = innerProduct(tmp_cb1, src_o);	
      axpy(lsol_full, cn / eval, tmp_cb1, lsol_full);

      if(n == nLowMode - 1) lsol_defl = lsol_full;
    }
    sol_o = lsol_full; //sol_o = lsol   Set guess equal to low mode projection 
  }

  f = norm2(src_o);
  if (!UniqueID()) printf("Grid_CGNE_M_high: CGNE_prec_MdagM src norm %le\n",f);
  f = norm2(sol_o);
  if (!UniqueID()) printf("Grid_CGNE_M_high: CGNE_prec_MdagM guess norm %le\n",f);

  //MdagM inverse controlled by evec interface
#ifndef MEMTEST_MODE
  evecs.CGNE_MdagM(linop, sol_o, src_o, resid, max_iters);
#endif
  
  f = norm2(sol_o);
  if (!UniqueID()) printf("Grid_CGNE_M_high: CGNE_prec_MdagM sol norm %le\n",f);


  //Pull low-mode part out of solution
  axpy(sol_o, -1.0, lsol_defl, sol_o);

  f = norm2(sol_o);
  if (!UniqueID()) printf("Grid_CGNE_M_high: sol norm after subtracting low-mode part %le\n",f);

  assert(sol_o.checkerboard == Grid::Odd);
  setCheckerboard(solution, sol_o);
  
  // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
  pickCheckerboard(Grid::Even,tmp_cb1,source);  //tmp_cb1 = src_e
  
  Ddwf.Meooe(sol_o,tmp_cb2); //tmp_cb2 = Meo sol_o
  assert(tmp_cb2.checkerboard == Grid::Even);

  axpy(tmp_cb1, -1.0, tmp_cb2, tmp_cb1); //tmp_cb1 = (-Meo sol_o + src_e)   (tmp_cb2 free)
  
  Ddwf.MooeeInv(tmp_cb1,tmp_cb2);  //tmp_cb2 = Mee^-1(-Meo sol_o + src_e)   (tmp_cb1 free)

  f = norm2(tmp_cb2);
  if (!UniqueID()) printf("Grid_CGNE_M_high: even checkerboard of sol %le\n",f);

  assert(tmp_cb2.checkerboard == Grid::Even);
  setCheckerboard(solution, tmp_cb2);

  f = norm2(solution);
  if (!UniqueID()) printf("Grid_CGNE_M_high: unprec sol norm is %le\n",f);
}

CPS_END_NAMESPACE
#endif
#endif
