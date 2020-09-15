int IRL (Vector ** f_out, Vector * f_in, CgArg * cg_arg,
	 int Nk, int Np, int Nstop, Float low, Float high, int order)
{
  const char *fname ("IRL()");
  SetMass (cg_arg->mass);
  SetEpsilon (cg_arg->epsilon);
  VRB.Result (cname, fname, "mass=%0.14g epsilon=%g \n", mass, eps);
//      VRB.Result(cname,fname,"max_iter=%d\n",cg_arg->max_num_iter);
  RealD M5 = GJP.DwfHeight ();

  ImportGauge ();
  std::vector < int >twists = SetTwist ();

  DIRAC::ImplParams params;
  SetParams (params);
  DIRAC DdwfD (*Umu, FIVE_GRID * UGridD, *UrbGridD, mass MOB PARAMS);
  Grid::QCD::LatticeGaugeFieldF UmuF (UGridF);
  precisionChange (UmuF, *Umu);
  DIRAC_F DdwfF (UmuF, FIVE_GRID_F * UGridF, *UrbGridF, mass MOB PARAMS);
#ifndef TWOKAPPA
  Grid::SchurDiagMooeeOperator < DIRAC, LATTICE_FERMION > HermOp (DdwfD);
  Grid::SchurDiagMooeeOperator < DIRAC_F, LATTICE_FERMION_F > HermOpF (DdwfF);
#else
  Grid::SchurDiagTwoOperator < DIRAC, LATTICE_FERMION > HermOp (DdwfD);
  Grid::SchurDiagTwoOperator < DIRAC_F, LATTICE_FERMION_F > HermOpF (DdwfF);
#endif

  LATTICE_FERMION grid_in (FERM_GRID), psi (FERM_GRID);
//      std::vector<LATTICE_FERMION> grid_rb_out(Nshift,F_RB_GRID);
  LATTICE_FERMION grid_rb_in (F_RB_GRID), grid_out (FERM_GRID);
  LATTICE_FERMION_F grid_rb_in_f (F_RB_GRID_F), grid_out_f (FERM_GRID_F);
  LatVector cps_temp (FsiteSize () / 6 * n_gp, GJP.VolNodeSites () / 2);
  VRB.Debug (cname, fname, "cps_temp.Size()=%d\n", cps_temp.Size ());
//      VRB.Debug(cname,fname,"Nshift=%d alpha=%p\n",Nshift,alpha);
  Vector *f_temp = cps_temp.Vec ();

  double fac = 1.;
  ImportFermion (grid_in, f_in, Odd);

  pickCheckerboard (Odd, grid_rb_in, grid_in);
  const int Nm = Nk + Np;
  const int MaxIt = 10000;
  RealD resid = 1.0e-8;
  RealD resid_f = 1.0e-6;

  std::vector < double >Coeffs { 0., -1.};
  Grid::Polynomial < LATTICE_FERMION > PolyX (Coeffs);
  Grid::Chebyshev < LATTICE_FERMION_F > Cheb (low, high, order + 1);
//  ChebyshevLanczos<LatticeFermion> Cheb(9.,1.,0.,20);
//  Cheb.csv(std::cout);
//  exit(-24);
//  ImplicitlyRestartedLanczos<FermionField> IRL(HermOp,Cheb,Nstop,Nk,Nm,resid,MaxIt);
  Grid::FunctionHermOp< LATTICE_FERMION_F > OpCheby(Cheb,HermOpF);
  Grid::PlainHermOp< LATTICE_FERMION_F > Op(HermOpF);
  Grid::ImplicitlyRestartedLanczos<LATTICE_FERMION_F> IRL(OpCheby,Op,Nstop,Nk,Nm,resid,MaxIt);
//  Grid::ImplicitlyRestartedLanczosCJ < LATTICE_FERMION_F > IRL (HermOpF, Cheb, Nstop, Nk, Nm, resid_f, MaxIt);
  precisionChange (grid_rb_in_f, grid_rb_in);

  std::vector < Grid::RealD > eval (Nm);
  std::vector < LATTICE_FERMION > evec (Nm, F_RB_GRID);
  std::vector < LATTICE_FERMION_F > evec_f (Nm, F_RB_GRID_F);

  int Nconv;
//  IRL.calc(eval,evec,grid_rb_in,Nconv);
  IRL.calc (eval, evec_f, grid_rb_in_f, Nconv);

  return Nconv;

}

#define SINGLE

int FeigSolv (Vector ** f_eigenv, Float * lambda,
	      LanczosArg * eig_arg, CnvFrmType cnv_frm)
{
  const char *fname ("FeigSolv()");
  if (cnv_frm == CNV_FRM_YES)
    ERR.General (cname, fname, " only has support for checkerboarded\n");
  int Ncb = 1; //hard coded for CNV_FRM_NO, Ncb=1
  SetMass (eig_arg->mass);
  SetEpsilon (0.);
  VRB.Result (cname, fname, "mass=%0.14g epsilon=%g \n", mass, eps);
//      VRB.Result(cname,fname,"max_iter=%d\n",cg_arg->max_num_iter);
  RealD M5 = GJP.DwfHeight ();

  ImportGauge ();
  std::vector < int >twists = SetTwist ();
  int mem_save = eig_arg->mem_save;
//  if (mem_save) ERR.General (cname, fname, " mem_save does not work yet!\n");
  int precision = eig_arg->precision;
  DIRAC::ImplParams params;
  SetParams (params);

//      Staying with double precision, at least until gaining more confidence in single precision
  DIRAC DdwfD (*Umu, FIVE_GRID * UGridD, *UrbGridD, mass MOB PARAMS);
#ifndef TWOKAPPA
  Grid::SchurDiagMooeeOperator < DIRAC, LATTICE_FERMION > HermOp (DdwfD);
#else
  Grid::SchurDiagTwoOperator < DIRAC, LATTICE_FERMION > HermOp (DdwfD);
//  Grid::SchurDiagTwoKappaOperator < DIRAC, LATTICE_FERMION > HermOp (DdwfD);
#endif

#ifdef SINGLE
  Grid::QCD::LatticeGaugeFieldF UmuF (UGridF);
  precisionChange (UmuF, *Umu);
  DIRAC_F DdwfF (UmuF, FIVE_GRID_F * UGridF, *UrbGridF, mass MOB PARAMS);
#ifndef TWOKAPPA
  VRB.Result (cname, fname, "Using Asym preconditioner\n");
  Grid::SchurDiagMooeeOperator < DIRAC_F, LATTICE_FERMION_F > HermOpF (DdwfF);
#else
  VRB.Result (cname, fname, "Using SYM2 preconditioner\n");
  Grid::SchurDiagTwoOperator < DIRAC_F, LATTICE_FERMION_F > HermOpF (DdwfF);
#endif
#endif

  LATTICE_FERMION grid_in (FERM_GRID), psi (FERM_GRID);
  LATTICE_FERMION grid_rb_in (F_RB_GRID), grid_out (FERM_GRID);
#ifdef SINGLE
  LATTICE_FERMION_F grid_rb_in_f (F_RB_GRID_F), grid_out_f (FERM_GRID_F);
#endif
  LatVector cps_temp (FsiteSize () / 6 * n_gp, GJP.VolNodeSites () / 2);
  VRB.Debug (cname, fname, "cps_temp.Size()=%d\n", cps_temp.Size ());
//      VRB.Debug(cname,fname,"Nshift=%d alpha=%p\n",Nshift,alpha);
  Vector *f_temp = cps_temp.Vec ();
  RandGaussVector (f_temp, 0.5, 1);

//      double fac = 1.;


  ImportFermion (grid_in, f_temp, Odd);
  pickCheckerboard (Odd, grid_rb_in, grid_in);

  const int Np = eig_arg->np_lanczos_vectors;
  const int Nk = eig_arg->nk_lanczos_vectors;
  const int Nm = Np + Nk;
  const int Nstop = eig_arg->nt_lanczos_vectors;
  const int MaxIt = eig_arg->maxiters;
  RealD resid = eig_arg->stop_residual;

  Float low = eig_arg->matpoly_arg.params.params_val[1];
  low = low * low;
  Float high = eig_arg->matpoly_arg.params.params_val[0];
  high = high * high;
  int order = eig_arg->matpoly_arg.Npol + 1;	//convention difference
  VRB.Result(cname,fname,"low=%g high=%g order=%d\n",low,high,order);

  if(mem_save==1)
    for (int i = Nstop-1; i >= 0  ; i--) {
    if (f_eigenv[i] != NULL) sfree(f_eigenv[i]);
	f_eigenv[i]=NULL;
  }
#ifndef SINGLE
  Grid::Chebyshev < LATTICE_FERMION > Cheb (low, high, order);
  Grid::ImplicitlyRestartedLanczosCJ < LATTICE_FERMION > IRL (HermOp, Cheb, Nstop,
							    Nk, Nm, resid,
							    MaxIt);
#else
  Grid::Chebyshev < LATTICE_FERMION_F > Cheb (low, high, order);
//  Grid::ImplicitlyRestartedLanczosCJ <LATTICE_FERMION_F> IRL(HermOpF,Cheb,Nstop,Nk,Nm,resid,MaxIt);
  Grid::FunctionHermOp< LATTICE_FERMION_F > OpCheby(Cheb,HermOpF);
  Grid::PlainHermOp< LATTICE_FERMION_F > Op(HermOpF);
  Grid::ImplicitlyRestartedLanczos<LATTICE_FERMION_F> IRL(OpCheby,Op,Nstop,Nk,Nm,resid,MaxIt);
  precisionChange(grid_rb_in_f,grid_rb_in);
#endif

  std::vector < Grid::RealD > eval (Nm);

  int Nconv;
#ifndef SINGLE
  std::vector < LATTICE_FERMION > evec (Nm, F_RB_GRID);
  for (int i = 0; i < Nm; i++)
    evec[i].checkerboard = Grid::Odd;
  IRL.calc (eval, evec, grid_rb_in, Nconv);
#else
  std::vector < LATTICE_FERMION_F > evec_f (Nm, F_RB_GRID_F);
  for (int i = 0; i < Nm; i++) evec_f[i].checkerboard = Grid::Odd;
  IRL.calc(eval,evec_f,grid_rb_in_f,Nconv);
#endif

	uint64_t f_size = (FsiteSize () * n_gp) *  (GJP.VolNodeSites () / 2) * Ncb;
	if(precision=PREC_SINGLE) f_size *=sizeof(float);
	else				 f_size *=sizeof(double);

  for (int i = Nstop-1; i >= 0  ; i--) {
    if( (mem_save==0) && (!f_eigenv[i] ) )
      ERR.General (cname, fname, "f_eigenv[%d] not allocated\n", i);
    lambda[i] = eval[i];

    if(mem_save==1){
	f_eigenv[i] = (Vector *) smalloc(cname,fname,"f_eigenv[i]",f_size);
    } 
#ifndef SINGLE
    setCheckerboard (grid_out, evec[i]);
    ImportFermion (f_eigenv[i], grid_out, Odd);
	evec.resize(i,F_RB_GRID);
#else
    setCheckerboard (grid_out_f, evec_f[i]);
//    ImportFermion (f_eigenv[i], grid_out_f, Odd);
// forcing it to be float, although it is ostensively double
    ImpexFermion < float,LATTICE_FERMION_F, SITE_FERMION_F > (f_eigenv[i], grid_out_f, 0, Odd, NULL);
    VRB.Result(cname,fname,"f_eigenv[%d]=%p\n",i,f_eigenv[i]);
	evec_f.resize(i,F_RB_GRID_F);
#endif
    
  }

  return Nconv;

}

#if 0 // SimpleLanczos only present in Lanczos branch
#ifndef IF_TM
#undef SINGLE
int SimpleLanczos (std::vector<Float>  &lambda, LanczosArg * eig_arg)
{
  const char *fname ("FeigSolv()");
  int Ncb = 1; //hard coded for CNV_FRM_NO, Ncb=1
  SetMass (eig_arg->mass);
  SetEpsilon (0.);
  VRB.Result (cname, fname, "mass=%0.14g epsilon=%g \n", mass, eps);
//      VRB.Result(cname,fname,"max_iter=%d\n",cg_arg->max_num_iter);
  RealD M5 = GJP.DwfHeight ();

  ImportGauge ();
  std::vector < int >twists = SetTwist ();
  int mem_save = eig_arg->mem_save;
  int precision = eig_arg->precision;
  DIRAC::ImplParams params;
  SetParams (params);

//      Staying with double precision, at least until gaining more confidence in single precision
  DIRAC DdwfD (*Umu, FIVE_GRID * UGridD, *UrbGridD, mass MOB PARAMS);
#ifndef TWOKAPPA
  Grid::SchurDiagMooeeOperator < DIRAC, LATTICE_FERMION > HermOp (DdwfD);
#else
  Grid::SchurDiagTwoOperator < DIRAC, LATTICE_FERMION > HermOp (DdwfD);
#endif

#ifdef SINGLE
  Grid::QCD::LatticeGaugeFieldF UmuF (UGridF);
  precisionChange (UmuF, *Umu);
  DIRAC_F DdwfF (UmuF, FIVE_GRID_F * UGridF, *UrbGridF, mass MOB PARAMS);
#ifndef TWOKAPPA
  VRB.Result (cname, fname, "Using Asym preconditioner\n");
  Grid::SchurDiagMooeeOperator < DIRAC_F, LATTICE_FERMION_F > HermOpF (DdwfF);
#else
  VRB.Result (cname, fname, "Using SYM2 preconditioner\n");
  Grid::SchurDiagTwoOperator < DIRAC_F, LATTICE_FERMION_F > HermOpF (DdwfF);
#endif
#endif

  LATTICE_FERMION grid_in (FERM_GRID), psi (FERM_GRID);
  LATTICE_FERMION grid_rb_in (F_RB_GRID), grid_out (FERM_GRID);
#ifdef SINGLE
  LATTICE_FERMION_F grid_rb_in_f (F_RB_GRID_F), grid_out_f (FERM_GRID_F);
#endif
  LatVector cps_temp (FsiteSize () / 6 * n_gp, GJP.VolNodeSites () / 2);
  VRB.Debug (cname, fname, "cps_temp.Size()=%d\n", cps_temp.Size ());
  Vector *f_temp = cps_temp.Vec ();
  RandGaussVector (f_temp, 0.5, 1);

  ImportFermion (grid_in, f_temp, Odd);
  pickCheckerboard (Odd, grid_rb_in, grid_in);
  std:: cout<<"norm2(grid_rb_n)= " << norm2(grid_rb_in)<<std::endl;

  const int Np = eig_arg->np_lanczos_vectors;
  const int Nk = eig_arg->nk_lanczos_vectors;
  const int Nm = Np + Nk;
  const int Nstop = eig_arg->nt_lanczos_vectors;
  const int MaxIt = eig_arg->maxiters;
  RealD resid = eig_arg->stop_residual;

  Float low = eig_arg->matpoly_arg.params.params_val[1];
  low = low * low;
  Float high = eig_arg->matpoly_arg.params.params_val[0];
  high = high * high;
  int order = eig_arg->matpoly_arg.Npol + 1;	//convention difference
  std::vector < Grid::RealD > eval (Nm);

  int Nconv;

  std::vector < double >Coeffs { 0., -1.};

#ifndef SINGLE
  Grid::Polynomial < LATTICE_FERMION > PolyX (Coeffs);
  Grid::Chebyshev < LATTICE_FERMION > Cheb (low, high, order);
  if(order >2){
  Grid::SimpleLanczos<LATTICE_FERMION> IRL(HermOp,Cheb,Nstop,Nk,Nm,resid,MaxIt);
  IRL.calc(eval,grid_rb_in,Nconv);
  } else {
  Grid::SimpleLanczos<LATTICE_FERMION> IRL(HermOp,PolyX,Nstop,Nk,Nm,resid,MaxIt);
  IRL.calc(eval,grid_rb_in,Nconv);
  }
#else
  Grid::Polynomial < LATTICE_FERMION_F > PolyX (Coeffs);
  Grid::Chebyshev < LATTICE_FERMION_F > Cheb (low, high, order);
  if(order >2){
  Grid::SimpleLanczos<LATTICE_FERMION_F> IRL(HermOpF,Cheb,Nstop,Nk,Nm,resid,MaxIt);
  precisionChange(grid_rb_in_f,grid_rb_in);
  IRL.calc(eval,grid_rb_in_f,Nconv);
  } else {
  Grid::SimpleLanczos<LATTICE_FERMION_F> IRL(HermOpF,PolyX,Nstop,Nk,Nm,resid,MaxIt);
  precisionChange(grid_rb_in_f,grid_rb_in);
  IRL.calc(eval,grid_rb_in_f,Nconv);
  }
#endif

  lambda.resize(Nconv);
  for(int i=0;i<Nconv;i++) eval[i] = lambda[i];


  return Nconv;

}
#endif //#ifndef IF_TM
#endif //#if 0
#undef SINGLE
