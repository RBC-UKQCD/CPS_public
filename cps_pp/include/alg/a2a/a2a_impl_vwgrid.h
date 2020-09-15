//Apply 1/2(1+-g5) to field. In Grid conventions this just zeroes the lower/upper spin components
template<typename FermionField>
void chiralProject(FermionField &out, const FermionField &in, const char sgn){
  int base; //where to start zeroing
  switch(sgn){
  case '+':
    base = 2;
    break;
  case '-':
    base = 0;
    break;
  default:
    assert(0);
  }
  
  out.checkerboard = in.checkerboard;
  conformable(in,out);

  const int Ns = 4;
  Grid::GridBase *grid=in._grid;

  //decltype(Grid::QCD::peekSpin(static_cast<const Grid::Lattice<typename FermionField::vector_object>&>(in),0)) zero_spn(in._grid);
  decltype(Grid::PeekIndex<SpinIndex>(in,0)) zero_spn(in._grid);
  Grid::zeroit(zero_spn);

  out = in;
  Grid::PokeIndex<SpinIndex>(out, zero_spn, base);
  Grid::PokeIndex<SpinIndex>(out, zero_spn, base+1);
}

//Convert a 5D field to a 4D field, with the upper 2 spin components taken from s-slice 's_u' and the lower 2 from 's_l'
template<typename FermionField>
void DomainWallFiveToFour(FermionField &out, const FermionField &in, int s_u, int s_l){
  assert(out._grid->Nd() == 4 && in._grid->Nd() == 5);

  FermionField tmp1_4d(out._grid);
  FermionField tmp2_4d(out._grid);
  FermionField tmp3_4d(out._grid);
  ExtractSlice(tmp1_4d,const_cast<FermionField&>(in),s_u, 0); //Note Grid conventions, s-dimension is index 0!
  chiralProject(tmp2_4d, tmp1_4d, '+'); // 1/2(1+g5)  zeroes lower spin components
  
  ExtractSlice(tmp1_4d,const_cast<FermionField&>(in),s_l, 0); 
  chiralProject(tmp3_4d, tmp1_4d, '-'); // 1/2(1-g5)  zeroes upper spin components

  out = tmp2_4d + tmp3_4d;
}
template<typename FermionField>
void DomainWallFourToFive(FermionField &out, const FermionField &in, int s_u, int s_l){
  assert(out._grid->Nd() == 5 && in._grid->Nd() == 4);

  zeroit(out);
  FermionField tmp1_4d(in._grid);
  chiralProject(tmp1_4d, in, '+'); // 1/2(1+g5)  zeroes lower spin components
  InsertSlice(tmp1_4d, out,s_u, 0);

  chiralProject(tmp1_4d, in, '-'); // 1/2(1-g5)  zeroes upper spin components
  InsertSlice(tmp1_4d, out,s_l, 0);
}

//Randomization of wh fields must have handled with care to ensure order preservation
template<typename ComplexFieldType>
struct _set_wh_random_impl<ComplexFieldType, grid_vector_complex_mark>{
  static void doit(std::vector<PtrWrapper<ComplexFieldType> > &wh, const RandomType &type, const int nhits){
    typedef typename Grid::GridTypeMapper<typename ComplexFieldType::FieldSiteType>::scalar_type ScalarComplexType;
    
    
    typedef CPSfield<ScalarComplexType, ComplexFieldType::FieldSiteSize,
		     typename ComplexFieldType::FieldDimensionPolicy::EquivalentScalarPolicy, typename ComplexFieldType::FieldFlavorPolicy, typename ComplexFieldType::FieldAllocPolicy>
      ScalarComplexFieldType;

    NullObject null_obj;
    
    //Use scalar generation code and import
    std::vector<PtrWrapper<ScalarComplexFieldType> > wh_scalar(nhits); for(int i=0;i<nhits;i++) wh_scalar[i].set(new ScalarComplexFieldType(null_obj));
    _set_wh_random_impl<ScalarComplexFieldType, complex_double_or_float_mark>::doit(wh_scalar,type,nhits);
    for(int i=0;i<nhits;i++) wh[i]->importField(*wh_scalar[i]);
  }
};

//Main implementations with generic interface
template< typename mf_Policies>
void A2AvectorW<mf_Policies>::computeVWlow(A2AvectorV<mf_Policies> &V, Lattice &lat, EvecInterface<mf_Policies> &evecs, const Float mass){
  if(!UniqueID()) printf("Computing VWlow using Grid\n");
  typedef typename mf_Policies::GridFermionField GridFermionField;
  typedef typename mf_Policies::FgridFclass FgridFclass;
  typedef typename mf_Policies::GridDirac GridDirac;
  
  const char *fname = "computeVQlow(....)";

  int ngp = 0;
  for(int i=0;i<3;i++) if(GJP.Bc(i) == BND_CND_GPARITY) ++ngp;

#ifdef USE_GRID_GPARITY
  if(ngp == 0) ERR.General("A2AvectorW","computeVWlow","Fgrid is currently compiled for G-parity\n");
#else
  if(ngp != 0) ERR.General("A2AvectorW","computeVWlow","Fgrid is not currently compiled for G-parity\n");
#endif

  assert(lat.Fclass() == mf_Policies::FGRID_CLASS_NAME);
  FgridFclass &latg = dynamic_cast<FgridFclass&>(lat);

  //Grids and gauge field
  Grid::GridCartesian *UGrid = latg.getUGrid();
  Grid::GridRedBlackCartesian *UrbGrid = latg.getUrbGrid();
  Grid::GridCartesian *FGrid = latg.getFGrid();
  Grid::GridRedBlackCartesian *FrbGrid = latg.getFrbGrid();
  Grid::QCD::LatticeGaugeFieldD *Umu = latg.getUmu();
  
  //Mobius parameters
  const double mob_b = latg.get_mob_b();
  const double mob_c = mob_b - 1.;   //b-c = 1
  const double M5 = GJP.DwfHeight();
  printf("Grid b=%g c=%g b+c=%g\n",mob_b,mob_c,mob_b+mob_c);

  const int gparity = GJP.Gparity();

  //Double precision temp fields
  CPSfermion4D<ComplexD> afield;
  CPSfermion5D<ComplexD> bfield;

  const int glb_ls = GJP.SnodeSites() * GJP.Snodes();

  //Setup Grid Dirac operator
  typename GridDirac::ImplParams params;
  latg.SetParams(params);

  GridDirac Ddwf(*Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,mob_b,mob_c, params);
  Grid::SchurDiagMooeeOperator<GridDirac,GridFermionField> linop(Ddwf);

  //Eigenvectors exist on odd checkerboard
  GridFermionField bq_tmp(FrbGrid);
  GridFermionField tmp(FrbGrid);
  GridFermionField tmp2(FrbGrid);
  GridFermionField tmp3(FrbGrid);

  GridFermionField tmp_full(FGrid);
  GridFermionField tmp_full2(FGrid);

  GridFermionField tmp_full_4d(UGrid);
  
  //The general method is described by page 60 of Daiqian's thesis
#ifndef MEMTEST_MODE
  for(int i = 0; i < nl; i++) {
    //Step 1) Compute V
    Float eval = evecs.getEvec(bq_tmp,i);
    assert(bq_tmp.checkerboard == Grid::Odd);

    //Compute  [ -(Mee)^-1 Meo bq_tmp, bg_tmp ]
    Ddwf.Meooe(bq_tmp,tmp2);	//tmp2 = Meo bq_tmp 
    Ddwf.MooeeInv(tmp2,tmp);   //tmp = (Mee)^-1 Meo bq_tmp
    tmp = -tmp; //even checkerboard
    
    assert(tmp.checkerboard == Grid::Even);
    
    setCheckerboard(tmp_full, tmp); //even checkerboard
    setCheckerboard(tmp_full, bq_tmp); //odd checkerboard

    //Get 4D part and poke into a
    //Recall that D^{-1} = <w^\dagger v> = <q \bar q>.  v therefore transforms like a conjugate spinor. For conjugate spinors \bar\psi(x) = P_R \bar\psi(x,Ls-1) + P_L \bar\psi(x,0),  i.e. s_u=Ls-1 and s_l=0 for CPS gamma5

    DomainWallFiveToFour(tmp_full_4d, tmp_full, glb_ls-1,0);
    tmp_full_4d = Grid::RealD(1./eval) * tmp_full_4d;
    V.getVl(i).importGridField(tmp_full_4d); //Multiply by 1/lambda[i] and copy into v (with precision change if necessary)
    
    
    //Step 2) Compute Wl

    //Do tmp = [ -[Mee^-1]^dag [Meo]^dag Doo bq_tmp,  Doo bq_tmp ]    (Note that for the Moe^dag in Daiqian's thesis, the dagger also implies a transpose of the spatial indices, hence the Meo^dag in the code)
    linop.Mpc(bq_tmp,tmp2);  //tmp2 = Doo bq_tmp
    
    Ddwf.MeooeDag(tmp2,tmp3); //tmp3 = Meo^dag Doo bq_tmp
    Ddwf.MooeeInvDag(tmp3,tmp); //tmp = [Mee^-1]^dag Meo^dag Doo bq_tmp
    tmp = -tmp;
    
    assert(tmp.checkerboard == Grid::Even);
    assert(tmp2.checkerboard == Grid::Odd);

    setCheckerboard(tmp_full, tmp);
    setCheckerboard(tmp_full, tmp2);

    //Left-multiply by D-^dag.  D- = (1-c*DW)
    Ddwf.DW(tmp_full, tmp_full2, 1);
    axpy(tmp_full, -mob_c, tmp_full2, tmp_full); 

    //Get 4D part, poke onto a then copy into wl
    //Recall that D^{-1} = <w^\dagger v> = <q \bar q>.  w (and w^\dagger) therefore transforms like a spinor. For spinors \psi(x) = P_R \bar\psi(x,0) + P_L \bar\psi(x,Ls-1),  i.e. s_u=0 and s_l=Ls-1 for CPS gamma5
    DomainWallFiveToFour(tmp_full_4d, tmp_full, 0, glb_ls-1);
    wl[i]->importGridField(tmp_full_4d);
  }
#endif
}



//Compute the high mode parts of V and W. 
//singleprec_evecs specifies whether the input eigenvectors are stored in single preciison
//You can optionally pass a single precision bfm instance, which if given will cause the underlying CG to be performed in mixed precision.
//WARNING: if using the mixed precision solve, the eigenvectors *MUST* be in single precision (there is a runtime check)
template< typename mf_Policies>
void A2AvectorW<mf_Policies>::computeVWhigh(A2AvectorV<mf_Policies> &V, Lattice &lat, EvecInterface<mf_Policies> &evecs, const Float mass, const Float residual, const int max_iter){
  typedef typename mf_Policies::GridFermionField GridFermionField;
  typedef typename mf_Policies::FgridFclass FgridFclass;
  typedef typename mf_Policies::GridDirac GridDirac;
  
  const char *fname = "computeVWhigh(....)";

  int ngp = 0;
  for(int i=0;i<3;i++) if(GJP.Bc(i) == BND_CND_GPARITY) ++ngp;

#ifdef USE_GRID_GPARITY
  if(ngp == 0) ERR.General("A2AvectorW","computeVWlow","Fgrid is currently compiled for G-parity\n");
#else
  if(ngp != 0) ERR.General("A2AvectorW","computeVWlow","Fgrid is not currently compiled for G-parity\n");
#endif

  assert(lat.Fclass() == mf_Policies::FGRID_CLASS_NAME);
  FgridFclass &latg = dynamic_cast<FgridFclass&>(lat);

  //Grids and gauge field
  Grid::GridCartesian *UGrid = latg.getUGrid();
  Grid::GridRedBlackCartesian *UrbGrid = latg.getUrbGrid();
  Grid::GridCartesian *FGrid = latg.getFGrid();
  Grid::GridRedBlackCartesian *FrbGrid = latg.getFrbGrid();
  Grid::QCD::LatticeGaugeFieldD *Umu = latg.getUmu();
  
  //Mobius parameters
  const double mob_b = latg.get_mob_b();
  const double mob_c = mob_b - 1.;   //b-c = 1
  const double M5 = GJP.DwfHeight();
  printf("Grid b=%g c=%g b+c=%g\n",mob_b,mob_c,mob_b+mob_c);

  const int gparity = GJP.Gparity();

  //Setup Grid Dirac operator
  typename GridDirac::ImplParams params;
  latg.SetParams(params);

  GridDirac Ddwf(*Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,mob_b,mob_c, params);
  Grid::SchurDiagMooeeOperator<GridDirac, GridFermionField> linop(Ddwf);

  VRB.Result("A2AvectorW", fname, "Start computing high modes using Grid.\n");
    
  //Generate the compact random sources for the high modes
#ifndef MEMTEST_MODE
  setWhRandom(args.rand_type);
#endif
  
  //Allocate temp *double precision* storage for fermions
  CPSfermion4D<typename mf_Policies::ComplexTypeD,typename mf_Policies::FermionFieldType::FieldDimensionPolicy, typename mf_Policies::FermionFieldType::FieldFlavorPolicy, typename mf_Policies::FermionFieldType::FieldAllocPolicy> v4dfield(wh[0]->getDimPolParams());
  
  const int glb_ls = GJP.SnodeSites() * GJP.Snodes();

  GridFermionField gtmp(FrbGrid);
  GridFermionField gtmp2(FrbGrid);
  GridFermionField gtmp3(FrbGrid);

  GridFermionField gsrc(FGrid);
  GridFermionField gtmp_full(FGrid);
  GridFermionField gtmp_full2(FGrid);

  GridFermionField tmp_full_4d(UGrid);

  //Details of this process can be found in Daiqian's thesis, page 60
#ifndef MEMTEST_MODE
  for(int i=0; i<nh; i++){
    //Step 1) Get the diluted W vector to invert upon
    getDilutedSource(v4dfield, i);

    //Step 2) Solve V
    v4dfield.exportGridField(tmp_full_4d);
    DomainWallFourToFive(gsrc, tmp_full_4d, 0, glb_ls-1);

    //Left-multiply by D-.  D- = (1-c*DW)
    Ddwf.DW(gsrc, gtmp_full, Grid::QCD::DaggerNo);
    axpy(gsrc, -mob_c, gtmp_full, gsrc); 

    //We can re-use previously computed solutions to speed up the calculation if rerunning for a second mass by using them as a guess
    //If no previously computed solutions this wastes a few flops, but not enough to care about
    //V vectors default to zero, so this is a zero guess if not reusing existing solutions
    V.getVh(i).exportGridField(tmp_full_4d);
    DomainWallFourToFive(gtmp_full, tmp_full_4d, 0, glb_ls-1);

    Ddwf.DW(gtmp_full, gtmp_full2, Grid::QCD::DaggerNo);
    axpy(gtmp_full, -mob_c, gtmp_full2, gtmp_full); 

    //Do the CG
    Grid_CGNE_M_high<mf_Policies>(gtmp_full, gsrc, residual, max_iter, evecs, nl, latg, Ddwf, FGrid, FrbGrid);
    
    //CPSify the solution, including 1/nhit for the hit average
    DomainWallFiveToFour(tmp_full_4d, gtmp_full, glb_ls-1,0);
    tmp_full_4d = Grid::RealD(1. / nhits) * tmp_full_4d;
    V.getVh(i).importGridField(tmp_full_4d);
  }
#endif
}


//Wrappers for generic interface

//BFM evecs
#ifdef USE_BFM_LANCZOS

//Compute the low mode part of the W and V vectors. In the Lanczos class you can choose to store the vectors in single precision (despite the overall precision, which is fixed to double here)
//Set 'singleprec_evecs' if this has been done
template< typename mf_Policies>
void A2AvectorW<mf_Policies>::computeVWlow(A2AvectorV<mf_Policies> &V, Lattice &lat, BFM_Krylov::Lanczos_5d<double> &eig, bfm_evo<double> &dwf, bool singleprec_evecs){
  EvecInterfaceBFM<mf_Policies> ev(eig,dwf,lat,singleprec_evecs);
  return computeVWlow(V,lat,ev,dwf.mass);
}

template< typename mf_Policies>
void A2AvectorW<mf_Policies>::computeVWhigh(A2AvectorV<mf_Policies> &V, BFM_Krylov::Lanczos_5d<double> &eig, bool singleprec_evecs, Lattice &lat, bfm_evo<double> &dwf_d, bfm_evo<float> *dwf_fp){
  bool mixed_prec_cg = dwf_fp != NULL; 
  if(mixed_prec_cg){
    //NOT IMPLEMENTED YET
    ERR.General("A2AvectorW","computeVWhigh","No grid implementation of mixed precision CG with BFM evecs\n");
  }

  if(mixed_prec_cg && !singleprec_evecs){ ERR.General("A2AvectorW","computeVWhigh","If using mixed precision CG, input eigenvectors must be stored in single precision"); }

  EvecInterfaceBFM<mf_Policies> ev(eig,dwf_d,lat,singleprec_evecs);
  return computeVWhigh(V,lat,ev,dwf_d.mass,dwf_d.residual,dwf_d.max_iter);
}

#endif



//Grid evecs
#ifdef USE_GRID_LANCZOS

template< typename mf_Policies>
void A2AvectorW<mf_Policies>::computeVWlow(A2AvectorV<mf_Policies> &V, Lattice &lat, const std::vector<typename mf_Policies::GridFermionField> &evec, const std::vector<Grid::RealD> &eval, const double mass){
  EvecInterfaceGrid<mf_Policies> ev(evec,eval);
  return computeVWlow(V,lat,ev,mass);
}

template< typename mf_Policies>
void A2AvectorW<mf_Policies>::computeVWhigh(A2AvectorV<mf_Policies> &V, Lattice &lat, const std::vector<typename mf_Policies::GridFermionField> &evec, const std::vector<Grid::RealD> &eval, const double mass, const Float residual, const int max_iter){
  EvecInterfaceGrid<mf_Policies> ev(evec,eval);
  return computeVWhigh(V,lat,ev,mass,residual,max_iter);
}

template< typename mf_Policies>
void A2AvectorW<mf_Policies>::computeVWlow(A2AvectorV<mf_Policies> &V, Lattice &lat, const std::vector<typename mf_Policies::GridFermionFieldF> &evec, const std::vector<Grid::RealD> &eval, const double mass){
  EvecInterfaceGridSinglePrec<mf_Policies> ev(evec,eval,lat,mass);
  return computeVWlow(V,lat,ev,mass);
}

template< typename mf_Policies>
void A2AvectorW<mf_Policies>::computeVWhigh(A2AvectorV<mf_Policies> &V, Lattice &lat, const std::vector<typename mf_Policies::GridFermionFieldF> &evec, const std::vector<Grid::RealD> &eval, const double mass, const Float residual, const int max_iter){
  EvecInterfaceGridSinglePrec<mf_Policies> ev(evec,eval,lat,mass);
  return computeVWhigh(V,lat,ev,mass,residual,max_iter);
}


#endif
