void setZmobius (std::vector < std::complex < double > >&omega, cps::Float * bs,
		 int ls)
{
  omega.clear ();
  for (int i = 0; i < ls; i++) {
    std::complex < double >temp2 (bs[2 * i], bs[2 * i + 1]);
    std::complex < double >temp = 1. / (2. * temp2 - 1.);
    VRB.Result ("FgridParams", "setZmobius", "bs[%d]=%g %g, omega=%g %g\n",
		i, bs[2 * i], bs[2 * i + 1], i, temp.real (), temp.imag ());
    omega.push_back (temp);
  }
}

void Schur (LATTICE_FERMION_F & grid_out_l, const LATTICE_FERMION_F & grid_in_l,
	    DIRAC_F & DdwfF, Float rsd, int maxinner)
{
  Grid::ConjugateGradient < LATTICE_FERMION_F > CG (rsd, maxinner, false);
  Grid::SchurRedBlackDiagTwoSolve < LATTICE_FERMION_F > SchurSolver (CG);
  SchurSolver (DdwfF, grid_in_l, grid_out_l);
}

void MixedSchur (LATTICE_FERMION & grid_out, LATTICE_FERMION & grid_in,
		 DIRAC & Ddwf, DIRAC_F & DdwfF, Float rsd, int maxinner,
		 int maxouter)
{
  Grid::QCD::SchurDifferentiableDiagTwo < IMPL > MdagMD (Ddwf);
  Grid::QCD::SchurDifferentiableDiagTwo < IMPL_F > MdagMF (DdwfF);
  Grid::MixedPrecisionConjugateGradient < LATTICE_FERMION, LATTICE_FERMION_F >
    CG_m (rsd, maxinner, maxouter, DdwfF.FermionRedBlackGrid (), MdagMF,
	  MdagMD);
  Grid::SchurRedBlackDiagTwoMixed < LATTICE_FERMION > SchurSolver (CG_m);
  SchurSolver (Ddwf, grid_in, grid_out);
}

void CheckFermion (Vector * cps_f, LATTICE_FERMION & grid_f)
{
  ImportFermion (cps_f, grid_f);
  Float *ftmp = (Float *) cps_f;
  printf ("CheckFermion: %g %g\n", ftmp[0], ftmp[1]);
}


//Copied from Christoph's GPT, except src and dst reversed
template < typename VV >
  void SpinProject (VV & _dst, const VV & _src, int inverse)
     const
     {
       _dst = _src;

       int Ls = _dst._grid->_ldimensions[0];
       int offset = inverse ? (Ls - 1) : (1);

       for (int i = 0; i < Ls; i++)
       {
	 Grid::QCD::axpby_ssp_pminus (_dst, 0., _src, 1., _src, i, i);
	 Grid::QCD::axpby_ssp_pplus (_dst, 1., _dst, 1., _src, i /* out */ ,
				     (i + offset) % Ls /*in */ );
       }
     }



int FmatInv (Vector * f_out,
	     Vector * f_in,
	     MobiusArg * mob_l,
	     MobiusArg * mob_s,
	     Float * true_res, CnvFrmType cnv_frm, PreserveType prs_f_in)
{
  const char *fname ("FmatInv()");
  RealD M5 = GJP.DwfHeight ();
  ImportGauge ();
  std::vector < int >twists = SetTwist ();

  DIRAC::ImplParams params;
  SetParams (params);

  // n_restart and stop_rsd is taken from mob_l
  int n_restart = mob_l->rsd_vec.rsd_vec_len;
  Float *rsd_vec = mob_l->rsd_vec.rsd_vec_val;

  CgArg _cg_arg_l = mob_l->cg;
  CgArg _cg_arg_s = mob_s->cg;

  CgArg *cg_arg_l = &_cg_arg_l;
  CgArg *cg_arg_s = &_cg_arg_s;


  std::vector < std::complex < double >>omegas_l;
  std::vector < std::complex < double >>omegas_s;

  setZmobius (omegas_l, mob_l->zmobius_b_coeff.zmobius_b_coeff_val,
	      mob_l->zmobius_b_coeff.zmobius_b_coeff_len);
  setZmobius (omegas_s, mob_s->zmobius_b_coeff.zmobius_b_coeff_val,
	      mob_s->zmobius_b_coeff.zmobius_b_coeff_len);

  int mob_l_ls = mob_l->ls;
  Grid::GridCartesian * FGridD_l;
  Grid::GridCartesian * FGridF_l;
  Grid::GridRedBlackCartesian * FrbGridF_l;
  Grid::GridRedBlackCartesian * FrbGridD_l;
  FGridD_l = Grid::QCD::SpaceTimeGrid::makeFiveDimGrid (mob_l_ls, UGridD);
  FGridF_l = Grid::QCD::SpaceTimeGrid::makeFiveDimGrid (mob_l_ls, UGridF);
  FrbGridD_l =
    Grid::QCD::SpaceTimeGrid::makeFiveDimRedBlackGrid (mob_l_ls, UGridD);
  FrbGridF_l =
    Grid::QCD::SpaceTimeGrid::makeFiveDimRedBlackGrid (mob_l_ls, UGridF);

  LATTICE_FERMION grid_in_l (FGridD_l), grid_out_l (FGridD_l);
  LATTICE_FERMION grid_src_l (FGridD_l);
  LATTICE_FERMION grid_sol_l (FGridD_l);

  int mob_s_ls = mob_s->ls;
  Grid::GridCartesian * FGridD_s;
  Grid::GridCartesian * FGridF_s;
  Grid::GridRedBlackCartesian * FrbGridF_s;
  Grid::GridRedBlackCartesian * FrbGridD_s;
  FGridD_s = Grid::QCD::SpaceTimeGrid::makeFiveDimGrid (mob_s_ls, UGridD);
  FGridF_s = Grid::QCD::SpaceTimeGrid::makeFiveDimGrid (mob_s_ls, UGridF);
  FrbGridD_s =
    Grid::QCD::SpaceTimeGrid::makeFiveDimRedBlackGrid (mob_s_ls, UGridD);
  FrbGridF_s =
    Grid::QCD::SpaceTimeGrid::makeFiveDimRedBlackGrid (mob_s_ls, UGridF);

  char fname_eig_root_bc[1024];
  snprintf (fname_eig_root_bc, 1024, "%s", cg_arg_s->fname_eigen);
  VRB.Result (cname, fname, "fname=%s\n", fname_eig_root_bc);
  int n_fields = GJP.SnodeSites ();
  const size_t f_size_per_site = (size_t) FsiteSize () / GJP.SnodeSites () / 2;
  int neig = 0;
  if_defl = false;
  if ((cg_arg_s->Inverter == CG_LOWMODE_DEFL)
      || (cg_arg_s->Inverter == LOWMODEAPPROX)) {
    if_defl = true;
    neig = cg_arg_s->neig;
  }

  std::cout << "neig= " << neig << std::endl;
  std::vector < Float > eval_s (neig);
  std::vector < LATTICE_FERMION_F > evec_f_s (neig, FrbGridF_s);
  Grid::Guesser < Float, LATTICE_FERMION_F > guesser (neig, eval_s, evec_f_s);
  if (if_defl) {
    LATTICE_FERMION_F grid_f (FGridF_s), grid_f_rb (FrbGridF_s);
    // search for eigen cache
    EigenCache *ecache;
    if ((ecache = EigenCacheListSearch ((char *) fname_eig_root_bc, neig)) == 0) {
      ERR.General (cname, fname,
		   "Eigenvector cache does not exist: neig %d name %s \n",
		   neig, fname_eig_root_bc);
    }
    Float mass = cg_arg_s->mass;
    EigenContainer eigcon (*this, (char *) fname_eig_root_bc, neig,
			   f_size_per_site / 2, n_fields, ecache);
    Float *eval_cps = eigcon.load_eval ();
    for (int i = 0; i < neig; i++) {
      eval_s[i] = eval_cps[i];
      VRB.Result (cname, fname, "eval[%d]=%0.14e\n", i, eval_s[i]);
      Vector *evec = eigcon.nev_load (i);
      ImpexFermion < float, LATTICE_FERMION_F, SITE_FERMION_F > (evec,
								 grid_f,
								 1, Odd, NULL);
      pickCheckerboard (Odd, evec_f_s[i], grid_f);
    }
  }


  LATTICE_FERMION grid_tmp2_l (FGridD_l);
  LATTICE_FERMION_F grid_fin_s (FGridF_s), grid_fout_s (FGridF_s);

  mass = cg_arg_l->mass;
  VRB.Result (cname, fname, "mass=%0.14g \n", mass);

  Grid::QCD::LatticeGaugeFieldF UmuF (UGridF);
  precisionChange (UmuF, *Umu);
  DIRAC Ddwf_l (*Umu, *FGridD_l, *FrbGridD_l, *UGridD, *UrbGridD, mass, M5,
		omegas_l, 1., 0. PARAMS);
  DIRAC_F DdwfF_l (UmuF, *FGridF_l, *FrbGridF_l, *UGridF, *UrbGridF, mass, M5,
		   omegas_l, 1., 0. PARAMS);
  DIRAC Ddwf_l_PV (*Umu, *FGridD_l, *FrbGridD_l, *UGridD, *UrbGridD, 1., M5,
		   omegas_l, 1., 0. PARAMS);
  DIRAC_F DdwfF_l_PV (UmuF, *FGridF_l, *FrbGridF_l, *UGridF, *UrbGridF, 1., M5,
		      omegas_l, 1., 0. PARAMS);

  DIRAC Ddwf_s (*Umu, *FGridD_s, *FrbGridD_s, *UGridD, *UrbGridD, mass, M5,
		omegas_s, 1., 0. PARAMS);
  DIRAC_F DdwfF_s (UmuF, *FGridF_s, *FrbGridF_s, *UGridF, *UrbGridF, mass, M5,
		   omegas_s, 1., 0. PARAMS);
  DIRAC Ddwf_s_PV (*Umu, *FGridD_s, *FrbGridD_s, *UGridD, *UrbGridD, 1., M5,
		   omegas_s, 1., 0. PARAMS);
  DIRAC_F DdwfF_s_PV (UmuF, *FGridF_s, *FrbGridF_s, *UGridF, *UrbGridF, 1., M5,
		      omegas_s, 1., 0. PARAMS);



  ImportFermion (grid_out_l, f_in);
  Ddwf_l.Dminus (grid_out_l, grid_in_l);	//left to right 
  grid_src_l = grid_in_l;	//dminus_in
  MixedSchur (grid_out_l, grid_in_l, Ddwf_l, DdwfF_l, rsd_vec[0],
	      cg_arg_l->max_num_iter, 1);
  Float stp_cnd = cg_arg_l->stop_rsd * cg_arg_l->stop_rsd * norm2 (grid_src_l);
  VRB.Result (cname, fname, "stp_cnd = %e\n", stp_cnd);
//  CheckFermion (f_out, grid_out_l); exit(-1);
  for (int i = 0; i < n_restart; i++) {
    if (0) {
      LATTICE_FERMION grid_tmp_l (FGridD_l);
      Ddwf_l.M (grid_out_l, grid_tmp_l);	//left to right
      grid_tmp_l = grid_in_l - grid_tmp_l;
      VRB.Result (cname, fname, "inner iter %d: residual %e\n",
		  i, norm2 (grid_tmp_l));
    }
    if (i == 0)
      grid_sol_l = grid_out_l;
    else
      grid_sol_l += grid_out_l;
    Ddwf_l.M (grid_sol_l, grid_out_l);	//left to right
    grid_in_l = grid_src_l - grid_out_l;	//deficit
    VRB.Result (cname, fname, "inner iter %d: residual %e\n",
		i, norm2 (grid_in_l));
    if (norm2 (grid_in_l) < stp_cnd)
      break;
    grid_out_l = 0.;
// 1st PV
    MixedSchur (grid_out_l, grid_in_l, Ddwf_l_PV, DdwfF_l_PV,
		1e-5, cg_arg_l->max_num_iter, 5);
    {
      LATTICE_FERMION grid_tmp_s (FGridD_s);
      LATTICE_FERMION grid_tmp2_s (FGridD_s);
      LATTICE_FERMION _c0 (UGridD);
      SpinProject (grid_tmp2_l, grid_out_l, 1);
      ExtractSlice (_c0, grid_tmp2_l, 0, 0);	// grid_tmp2_l should be kept
      grid_tmp_s = 0.;
      InsertSlice (_c0, grid_tmp_s, 0, 0);	//left to right
      SpinProject (grid_tmp2_s, grid_tmp_s, 0);
      Ddwf_s_PV.M (grid_tmp2_s, grid_tmp_s);	// left to right
      precisionChange (grid_fin_s, grid_tmp_s);
    }

//    Schur (grid_fout_s, grid_fin_s, DdwfF_s, rsd_vec[i], cg_arg_s->max_num_iter);
//void Schur (LATTICE_FERMION_F & grid_out_l, const LATTICE_FERMION_F & grid_in_l,
//          DIRAC_F & DdwfF, Float rsd, int maxinner)
    {
      Grid::ConjugateGradient < LATTICE_FERMION_F >
	CG (rsd_vec[i], cg_arg_s->max_num_iter, false);
      if (if_defl) guesser(grid_fin_s,grid_fout_s);
      Grid::SchurRedBlackDiagTwoSolve < LATTICE_FERMION_F > SchurSolver (CG);
      SchurSolver (DdwfF_s, grid_fin_s, grid_fout_s);
    }
    {
      LATTICE_FERMION grid_tmp_s (FGridD_s);
      LATTICE_FERMION grid_tmp_l (FGridD_l);
      LATTICE_FERMION _yprime (FGridD_s);
      LATTICE_FERMION _y0prime (UGridD);
      LATTICE_FERMION _ci (UGridD);
      precisionChange (grid_tmp_s, grid_fout_s);
      SpinProject (_yprime, grid_tmp_s, 1);
      ExtractSlice (_y0prime, _yprime, 0, 0);
      grid_tmp_l = 0.;
      InsertSlice (_y0prime, grid_tmp_l, 0, 0);
      grid_tmp_l = (-1.) * grid_tmp_l;
      for (int i = 1; i < mob_l_ls; i++) {
	ExtractSlice (_ci, grid_tmp2_l, i, 0);
	InsertSlice (_ci, grid_tmp_l, i, 0);
      }
      SpinProject (grid_tmp2_l, grid_tmp_l, 0);
      Ddwf_l.M (grid_tmp2_l, grid_tmp_l);	// left to right
      MixedSchur (grid_tmp2_l, grid_tmp_l, Ddwf_l_PV, DdwfF_l_PV, 1e-5,
		  cg_arg_l->max_num_iter, 5);
      SpinProject (grid_tmp_l, grid_tmp2_l, 1);
      InsertSlice (_y0prime, grid_tmp_l, 0, 0);	// left to right
      SpinProject (grid_out_l, grid_tmp_l, 0);
    }

  }

  {
    LATTICE_FERMION temp (FGridD_l);
    Ddwf_l.M (grid_sol_l, temp);
    temp = temp - grid_src_l;
    Float temp_res = std::sqrt (norm2 (temp) / norm2 (grid_src_l));
    VRB.Result (cname, fname, "true_res before final solve =%g\n", temp_res);
    if (true_res)
      *true_res = temp_res;
  }

#if 1
  MixedSchur (grid_sol_l, grid_src_l, Ddwf_l, DdwfF_l, cg_arg_l->stop_rsd,
	      cg_arg_l->max_num_iter, 5);

#else
//  ImportFermion (grid_out_l, f_out);
  {
    int maxinner = cg_arg_l->max_num_iter;
    int maxouter = 10;
    Grid::QCD::SchurDifferentiableDiagTwo < IMPL > MdagMD (Ddwf_l);
    Grid::QCD::SchurDifferentiableDiagTwo < IMPL_F > MdagMF (DdwfF_l);
    Grid::MixedPrecisionConjugateGradient < LATTICE_FERMION,
      LATTICE_FERMION_F > CG_m (stp_cnd, maxinner, maxouter,
				F_RB_GRID_F, MdagMF, MdagMD);
    Grid::SchurRedBlackDiagTwoMixed < LATTICE_FERMION > SchurSolver (CG_m);
    SchurSolver (Ddwf_l, grid_in_l, grid_out_l);
  }
#endif
  {
    LATTICE_FERMION temp (FGridD_l);
    Ddwf_l.M (grid_sol_l, temp);
    temp = temp - grid_src_l;
    Float temp_res = std::sqrt (norm2 (temp) / norm2 (grid_src_l));
    VRB.Result (cname, fname, "true_res=%g\n", temp_res);
    if (true_res)
      *true_res = temp_res;
  }
  ImportFermion (f_out, grid_sol_l);
// probably doesn't need to create and destroy *_l ones, or even both. 
  delete (FGridD_l);
  delete (FGridF_l);
  delete (FrbGridD_l);
  delete (FrbGridF_l);
  delete (FGridD_s);
  delete (FGridF_s);
  delete (FrbGridD_s);
  delete (FrbGridF_s);

  return 0;
}
