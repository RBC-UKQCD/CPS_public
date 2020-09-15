
#ifdef IF_FIVE_D
#define FIVE_GRID *FGridD,*FrbGridD,
#define FIVE_GRID_F *FGridF,*FrbGridF,
#define FERM_GRID FGridD
#define FERM_GRID_F FGridF
#define F_RB_GRID FrbGridD
#define F_RB_GRID_F FrbGridF
#else
#define FIVE_GRID
#define FIVE_GRID_F
#define FERM_GRID UGridD
#define FERM_GRID_F UGridF
#define F_RB_GRID UrbGridD
#define F_RB_GRID_F UrbGridF
#endif

#define LATTICE_FERMION DIRAC ::FermionField
#define LATTICE_FERMION_F DIRAC_F ::FermionField

#ifdef TWOKAPPA
#include<Grid/qcd/action/fermion/SchurDiagTwoKappa.h>
#endif


CPS_START_NAMESPACE class FGRID:public virtual Lattice, public FgridParams,
  public virtual FgridBase
{

  using RealD = Grid::RealD;
  using RealF = Grid::RealF;
public:
  const char *cname;

private:
    template < typename CPSfloat, typename vobj, typename sobj >
    void ImpexFermion (Vector * cps_f, vobj & grid_f, int cps2grid, EvenOdd eo,
		       Float * fac)
  {
    const char *fname = "ImpexFermion()";
      VRB.Debug (cname, fname, "%p %p %d %d F5D()=%d\n", cps_f, &grid_f,
		 cps2grid, eo, this->F5D ());
      fflush (stdout);
      fflush (stderr);
    if (fac)
        VRB.Result (cname, fname, "fac=%g\n", *fac);
    if (!cps_f)
        ERR.General (cname, fname, "CPS field not allocated\n");
    unsigned long vol;
      Grid::GridBase * grid = grid_f._grid;
    unsigned long fourvol = GJP.VolNodeSites ();
    int ncb = 2;
    if (eo != All)
    {
      fourvol /= 2;
      ncb = 1;
    }
    vol = fourvol * (2 / ncb);
    if (F5D ())
      vol *= GJP.SnodeSites ();
    if (((grid_f._grid)->lSites ()) != vol)
      ERR.General (cname, fname,
		   "numbers of grid(%d) and GJP(%d) does not match\n",
		   grid->lSites (), vol);
#pragma omp parallel for
    for (int site = 0; site < vol; site++) {
//        VRB.Result(cname,fname,"thread=%d of %d\n",omp_get_thread_num(),omp_get_num_threads());
      std::vector < int >grid_coor;
      sobj siteGrid;
      Grid::Lexicographic::CoorFromIndex (grid_coor, site, grid->_ldimensions);
      int pos[4], offset = 0;
      if (F5D ())
	offset = 1;
      for (int i = 0; i < 4; i++)
	pos[i] = grid_coor[i + offset];
      if ((ncb > 1) || (pos[0] + pos[1] + pos[2] + pos[3]) % 2 == eo)
	if (cps2grid) {
	  //                    printf("%d: %s:%s: site=%d grid_coor=%d %d %d %d %d 4dindex = %d \n", UniqueID(),cname,fname, site,grid_coor[0],grid_coor[1],grid_coor[2],grid_coor[3],grid_coor[4], (FsiteOffset(pos)/(2/ncb)) );fflush(stdout);
	  for (int gp = 0; gp < n_gp; gp++)
	    for (int s = 0; s < Ns; s++)
	      for (int i = 0; i < Nc; i++) {
		int i_gp = 0;
		if (F5D ())
		  i_gp = grid_coor[0];	// s coor when 5D
		int index =
		  (FsiteOffset (pos) / (2 / ncb)) + fourvol * (gp +
							       n_gp * i_gp);
//              Float *cps = (Float *) cps_f;
		CPSfloat *cps = (CPSfloat *) cps_f;
		cps += 2 * (i + Nc * (s + Ns * (index)));
		std::complex < double >elem (*cps, *(cps + 1));
//                              if (norm(elem)>0.01) printf("%d %d %d %d: cps[%d][%d][%d][%d] = %g %g\n", GJP.NodeCoor(0),GJP.NodeCoor(1),GJP.NodeCoor(2),GJP.NodeCoor(3), site,gp,s,i,elem.real(),elem.imag());
		if (fac)
		  siteGrid (GP) (s) (i) = elem * (*fac);
		else
		  siteGrid (GP) (s) (i) = elem;
	      }
	  pokeLocalSite (siteGrid, grid_f, grid_coor);
	} else {
	  Grid::Lexicographic::CoorFromIndex (grid_coor, site,
					      grid->_ldimensions);
	  peekLocalSite (siteGrid, grid_f, grid_coor);
	  for (int gp = 0; gp < n_gp; gp++)
	    for (int s = 0; s < Ns; s++)
	      for (int i = 0; i < Nc; i++) {
		int i_gp = 0;
		if (F5D ())
		  i_gp = grid_coor[0];	// s coor when 5D
		std::complex < double >elem;
		elem = siteGrid (GP) (s) (i);
//                              if (norm(elem)>0.01) printf("%d %d %d %d %d: grid[%d][%d][%d][%d] = %g %g\n",
//grid_coor[0], grid_coor[1], grid_coor[2], grid_coor[3], grid_coor[4],
//site,gp,s,i,elem.real(),elem.imag());
		int index =
		  (FsiteOffset (pos) / (2 / ncb)) + fourvol * (gp +
							       n_gp * i_gp);
//              Float *cps = (Float *) cps_f;
		CPSfloat *cps = (CPSfloat *) cps_f;
		cps += 2 * (i + Nc * (s + Ns * (index)));
		if (fac) {
		  *cps = elem.real () * (*fac);
		  *(cps + 1) = elem.imag () * (*fac);
		} else {
		  *cps = elem.real ();
		  *(cps + 1) = elem.imag ();
		}
	      }
	}
    }
  }

public:
  bool if_defl;
  void PrintMass (const char *fname)
  {
#ifdef IF_TM
    VRB.Result (cname, fname, "mass=%0.14g epsilon=%g \n", mass, this->eps);
#else
    VRB.Result (cname, fname, "mass=%0.14g\n", mass);
#endif
  }


FGRID (FgridParams & params):cname (XSTR (FGRID)), FgridBase (params),
    if_defl (false) {
#ifdef GRID_GPARITY
    if (!GJP.Gparity ())
      ERR.General (cname, XSTR (FGRID),
		   "Trying to instantiate Fgrid class with Gparity on non-Gparity lattice\n");
#endif
  }
  void ImportFermion (Vector * cps_f, LATTICE_FERMION & grid_f, EvenOdd eo =
		      All, Float * fac = NULL) {
    ImpexFermion < Float, LATTICE_FERMION, SITE_FERMION > (cps_f, grid_f, 0, eo,
							   fac);
  }

  void ImportFermion (LATTICE_FERMION & grid_f, Vector * cps_f, EvenOdd eo =
		      All, Float * fac = NULL) {
    ImpexFermion < Float, LATTICE_FERMION, SITE_FERMION > (cps_f, grid_f, 1, eo,
							   fac);
  }

  void ImportFermion (Vector * cps_f, LATTICE_FERMION_F & grid_f, EvenOdd eo =
		      All, Float * fac = NULL) {
    ImpexFermion < Float, LATTICE_FERMION_F, SITE_FERMION_F > (cps_f, grid_f, 0,
							       eo, fac);
  }


  void ImportFermion (LATTICE_FERMION_F & grid_f, Vector * cps_f, EvenOdd eo =
		      All, Float * fac = NULL) {
    ImpexFermion < Float, LATTICE_FERMION_F, SITE_FERMION_F > (cps_f, grid_f, 1,
							       eo, fac);
  }

  void SetParams (DIRAC::ImplParams & params)
  {
#ifdef GRID_GPARITY
    std::vector < int >twists = SetTwist ();
    params.twists = twists;
#endif
//#ifdef GRID_ZMOB
#if 0
//      this->setZmobius(GJP.zmobius_b);
// assumes b=1 c=0
    this->omegas.clear ();
    for (int i = 0; i < GJP.zmobius_b.size (); i++) {
      std::complex < double >temp = 1. / (2. * GJP.zmobius_b[i] - 1.);
      VRB.Result ("FgridParams", "setZmobius", "bs[%d]=%g %g, omega=%g %g\n",
		  i, GJP.zmobius_b[i].real (), GJP.zmobius_b[i].imag (), i,
		  temp.real (), temp.imag ());
      this->omegas.push_back (temp);
    }

#endif
  }



  void Fdslash (Vector * f_out, Vector * f_in, CgArg * cg_arg,
		CnvFrmType cnv_frm, int dir_flag)
  {
    const char *fname ("Fdslash()");
#if 1
//      mass = cg_arg->mass;
    SetMass (cg_arg->mass);
    SetEpsilon (cg_arg->epsilon);
    PrintMass (fname);
    RealD M5 = GJP.DwfHeight ();
    ImportGauge ();
    std::vector < int >twists = SetTwist ();


    LATTICE_FERMION grid_in (FERM_GRID), grid_out (FERM_GRID);
    ImportFermion (grid_in, f_in);
    DIRAC::ImplParams params;
    SetParams (params);
    DIRAC Ddwf (*Umu, FIVE_GRID * UGridD, *UrbGridD, mass MOB PARAMS);
    Ddwf.M (grid_in, grid_out);
#ifdef IF_TM
    LATTICE_FERMION grid_tmp (FERM_GRID);
    grid_tmp = Grid::QCD::Gamma (Grid::QCD::Gamma::Algebra::Gamma5) * grid_in;
    grid_out += std::complex < double >(0, 1.) * grid_tmp;
#endif
    ImportFermion (f_out, grid_out);
#endif
  }

  // It returns the type of fermion class
  FclassType Fclass () const
  {
    return CLASS_NAME;
  }

  int FsiteSize () const
  {
#ifdef IF_FIVE_D
    int size = 24 * Ls;
#else
    int size = 24;
#endif
      return size;
  }
  // Returns the number of fermion field 
  // components (including real/imaginary) on a
  // site of the 4-D lattice.

  int FchkbEvl () const
  {
    return 1;
  }
  // Returns 0 => If no checkerboard is used for the evolution
  //      or the CG that inverts the evolution matrix.

  // It calculates f_out where A * f_out = f_in and
  // A is the preconditioned fermion matrix that appears
  // in the HMC evolution (even/odd preconditioning 
  // of [Dirac^dag Dirac]). The inversion is done
  // with the conjugate gradient. cg_arg is the structure
  // that contains all the control parameters, f_in is the
  // fermion field source vector, f_out should be set to be
  // the initial guess and on return is the solution.
  // f_in and f_out are defined on a checkerboard.
  // If true_res !=0 the value of the true residual is returned
  // in true_res.
  // *true_res = |src - MatPcDagMatPc * sol| / |src|
  // The function returns the total number of CG iterations.
  int FmatEvlInv (Vector * f_out, Vector * f_in,
		  CgArg * cg_arg,
		  Float * true_res, CnvFrmType cnv_frm = CNV_FRM_YES) {
    const char *fname ("FmatEvlInv()");

    SetMass (cg_arg->mass);
    SetEpsilon (cg_arg->epsilon);
    PrintMass (fname);
    VRB.Result (cname, fname, "max_iter=%d\n", cg_arg->max_num_iter);
    RealD M5 = GJP.DwfHeight ();

    Float time_total, time1, time2;
    time_total = time1 = dclock ();
    int total; //total iteration count

    ImportGauge ();
    std::vector < int >twists = SetTwist ();
    int threads = Grid::GridThread::GetThreads ();
    std::
      cout << Grid::GridLogMessage << "Grid is setup to use " << threads <<
      " threads" << std::endl;
    if (0) {
      Float max_dev = 1;
      Float max_diff = 1;
      CheckUnitarity (max_dev, max_diff);
      VRB.Result (cname, fname, "CheckUnitarity():dev=%g max_diff=%g\n",
		  max_dev, max_diff);
    }


    DIRAC::ImplParams params;
    SetParams (params);
    DIRAC DdwfD (*Umu, FIVE_GRID * UGridD, *UrbGridD, mass MOB PARAMS);
    Grid::QCD::LatticeGaugeFieldF UmuF (UGridF);
    precisionChange (UmuF, *Umu);
    DIRAC_F DdwfF (UmuF, FIVE_GRID_F * UGridF, *UrbGridF, mass MOB PARAMS);

    LATTICE_FERMION grid_in (FERM_GRID), psi (FERM_GRID);
    LATTICE_FERMION grid_rb_out (F_RB_GRID);
    LATTICE_FERMION grid_rb_in (F_RB_GRID), grid_out (FERM_GRID);

    LatVector cps_temp (FsiteSize () / 6 * n_gp, GJP.VolNodeSites () / 2);
    VRB.Debug (cname, fname, "cps_temp.Size()=%d\n", cps_temp.Size ());
//    VRB.Debug (cname, fname, "Nshift=%d alpha=%p\n", Nshift, alpha);
    Vector *f_temp = cps_temp.Vec ();

#if (defined IF_TM ) && ( defined USE_F_CLASS_WILSON_TM)
    Float kappa = cg_arg[0]->mass + 4.;
    kappa = kappa * kappa + cg_arg[0]->epsilon * cg_arg[0]->epsilon;
    kappa = 0.5 / sqrt (kappa);
    double fac = 0.25 / (kappa * kappa);
    VRB.Debug (cname, fname, "kappa=%0.14g fac=%g \n", kappa, fac);
#else
    double fac = 1.;
#endif

    double inner_fac = 1.;

    ImportFermion (grid_in, f_in, Odd);

    pickCheckerboard (Odd, grid_rb_in, grid_in);
#ifndef TWOKAPPA
    Grid::QCD::SchurDifferentiableOperator < IMPL > MdagMD (DdwfD);
    Grid::QCD::SchurDifferentiableOperator < IMPL_F > MdagMF (DdwfF);
#else
    Grid::QCD::SchurDifferentiableDiagTwo < IMPL > MdagMD (DdwfD);
    Grid::QCD::SchurDifferentiableDiagTwo < IMPL_F > MdagMF (DdwfF);
#endif
    VRB.Result (cname, fname, "MultiShiftController.getMode()=%d\n",
		MultiShiftController.getMode ());
    if (MultiShiftController.getMode () == MultiShiftCGcontroller::DOUBLE_PREC) {
      VRB.Result (cname, fname, "running full double\n");
      time1 = dclock ();
      Float stp_cnd = cg_arg->stop_rsd;
      Grid::ConjugateGradient < LATTICE_FERMION > CG_d (stp_cnd,
							cg_arg->max_num_iter);
      CG_d (MdagMD, grid_rb_in, grid_rb_out);
      time2 = dclock ();

#if 0
      Float total_flops =
	(Float) GJP.VolNodeSites () * n_gp * (2 * 1440) * CG_d.iter;
      if (F5D ())
	total_flops *= Ls;
      print_flops (fname, "CG(D)", total_flops, time2 - time1);
#endif
      time1 = time2;
    } else {
      VRB.Result (cname, fname,
		  "running single precision Multimass followed by shifted mixed precision CG\n");
      int maxouter = 10;
      int maxinner = cg_arg->max_num_iter;
      int maxmulti = cg_arg->max_num_iter;
      LATTICE_FERMION_F grid_in_f (FERM_GRID_F), psi_f (FERM_GRID_F);
      LATTICE_FERMION_F grid_rb_out_f (F_RB_GRID_F);
      LATTICE_FERMION_F grid_rb_in_f (F_RB_GRID_F);
      Grid::MixedPrecisionConjugateGradient < LATTICE_FERMION,
	LATTICE_FERMION_F > CG_m (cg_arg->stop_rsd, maxinner, maxouter,
				  F_RB_GRID_F, MdagMF, MdagMD);
      CG_m (grid_rb_in, grid_rb_out);
      total = CG_m.TotalInnerIterations;
      time2 = dclock ();
      print_time (fname, "MixedCG", time2 - time1);
      time1 = time2;
    }
    setCheckerboard (grid_out, grid_rb_out);
    ImportFermion (f_out, grid_out, Odd);


    time2 = dclock ();
    print_time (fname, "Acc/Export", time2 - time1);
    time1 = time2;
    // print_flops(cname,fname,time2,-time_total,
   return total;
  }

  int FmatEvlMInv (Vector ** f_out, Vector * f_in, Float * shift,
		   int Nshift, int isz, CgArg ** cg_arg,
		   CnvFrmType cnv_frm, MultiShiftSolveType type, Float * alpha,
		   Vector ** f_out_d)
  {
    const char *fname ("FmatEvlMInv()");
//      mass = cg_arg[0]->mass;
    SetMass (cg_arg[0]->mass);
    SetEpsilon (cg_arg[0]->epsilon);
    PrintMass (fname);
    VRB.Result (cname, fname, "max_iter=%d\n", cg_arg[0]->max_num_iter);
    RealD M5 = GJP.DwfHeight ();

    Float time_total, time1, time2;
    time_total = time1 = dclock ();

    ImportGauge ();
    std::vector < int >twists = SetTwist ();
#if 1
    int threads = Grid::GridThread::GetThreads ();
    std::cout << Grid::GridLogMessage << "Grid is setup to use " << threads <<
      " threads" << std::endl;
#endif


    DIRAC::ImplParams params;
    SetParams (params);
    DIRAC DdwfD (*Umu, FIVE_GRID * UGridD, *UrbGridD, mass MOB PARAMS);
    Grid::QCD::LatticeGaugeFieldF UmuF (UGridF);
    precisionChange (UmuF, *Umu);
    DIRAC_F DdwfF (UmuF, FIVE_GRID_F * UGridF, *UrbGridF, mass MOB PARAMS);

    LATTICE_FERMION grid_in (FERM_GRID), psi (FERM_GRID);
    std::vector < LATTICE_FERMION > grid_rb_out (Nshift, F_RB_GRID);
    LATTICE_FERMION grid_rb_in (F_RB_GRID), grid_out (FERM_GRID);

    LatVector cps_temp (FsiteSize () / 6 * n_gp, GJP.VolNodeSites () / 2);
    VRB.Debug (cname, fname, "cps_temp.Size()=%d\n", cps_temp.Size ());
    VRB.Debug (cname, fname, "Nshift=%d alpha=%p\n", Nshift, alpha);
    Vector *f_temp = cps_temp.Vec ();

    Grid::MultiShiftFunction shifts;
    Grid::MultiShiftFunction shiftsF;
    shifts.order = Nshift;
    shifts.poles.resize (Nshift);
    shifts.tolerances.resize (Nshift);
    shifts.residues.resize (Nshift);
    shiftsF.order = Nshift;
    shiftsF.poles.resize (Nshift);
    shiftsF.tolerances.resize (Nshift);
    shiftsF.residues.resize (Nshift);
//#ifdef IF_TM This is needed if F_CLASS_WILSON_TM is used
#if (defined IF_TM ) && ( defined USE_F_CLASS_WILSON_TM)
    Float kappa = cg_arg[0]->mass + 4.;
    kappa = kappa * kappa + cg_arg[0]->epsilon * cg_arg[0]->epsilon;
    kappa = 0.5 / sqrt (kappa);
    double fac = 0.25 / (kappa * kappa);
    VRB.Debug (cname, fname, "kappa=%0.14g fac=%g \n", kappa, fac);
#else
    double fac = 1.;
#endif
    double inner_fac = 1.;
    for (int i = 0; i < Nshift; i++) {
      shiftsF.poles[i] = shifts.poles[i] = shift[i] * fac;
      if (type == SINGLE)
	shifts.residues[i] = alpha[i];
      if (type == SINGLE)
	shiftsF.residues[i] = alpha[i];
      shiftsF.tolerances[i] = shifts.tolerances[i] = cg_arg[i]->stop_rsd;
      if (fabs (shiftsF.tolerances[i]) < 1e-8)
	shiftsF.tolerances[i] = 1e-8;
    }

    ImportFermion (grid_in, f_in, Odd);

    pickCheckerboard (Odd, grid_rb_in, grid_in);
#ifndef TWOKAPPA
    Grid::QCD::SchurDifferentiableOperator < IMPL > MdagMD (DdwfD);
    Grid::QCD::SchurDifferentiableOperator < IMPL_F > MdagMF (DdwfF);
#else
    Grid::QCD::SchurDifferentiableDiagTwo < IMPL > MdagMD (DdwfD);
    Grid::QCD::SchurDifferentiableDiagTwo < IMPL_F > MdagMF (DdwfF);
#endif
    VRB.Result (cname, fname, "MultiShiftController.getMode()=%d\n",
		MultiShiftController.getMode ());
    if (MultiShiftController.getMode () == MultiShiftCGcontroller::DOUBLE_PREC)
//if(0)
    {
      VRB.Result (cname, fname, "running full double\n");
      time1 = dclock ();
      Grid::ConjugateGradientMultiShift < LATTICE_FERMION >
	MSCG (cg_arg[0]->max_num_iter, shifts);
      MSCG (MdagMD, grid_rb_in, grid_rb_out);
      time2 = dclock ();
#if 0
      Float total_flops =
	(Float) GJP.VolNodeSites () * n_gp * (2 * 1440) * MSCG.iter;
      if (F5D ())
	total_flops *= Ls;
      print_flops (fname, "MultiCG(D)", total_flops, time2 - time1);
#endif
      time1 = time2;
    } else {
      VRB.Result (cname, fname,
		  "running single precision Multimass followed by shifted mixed precision CG\n");
      int maxouter = 10;
      int maxinner = cg_arg[0]->max_num_iter;
      int maxmulti = cg_arg[0]->max_num_iter;
      LATTICE_FERMION_F grid_in_f (FERM_GRID_F), psi_f (FERM_GRID_F);
      std::vector < LATTICE_FERMION_F > grid_rb_out_f (Nshift, F_RB_GRID_F);
      LATTICE_FERMION_F grid_rb_in_f (F_RB_GRID_F);
//#ifdef IF_TM
// reliable update not working! reverting back to single precision for now
#if 1
      maxouter = cg_arg[0]->max_num_iter;
      maxinner = cg_arg[0]->max_num_iter;
      maxmulti = cg_arg[0]->max_num_iter;

      precisionChange (grid_rb_in_f, grid_rb_in);

      time2 = dclock ();
      print_time (fname, "setup()", time2 - time1);
      time1 = time2;
      Grid::ConjugateGradientMultiShift < LATTICE_FERMION_F > MSCG (maxmulti,
								    shiftsF);
      MSCG (MdagMF, grid_rb_in_f, grid_rb_out_f);
      time2 = dclock ();

#if 0
      Float total_flops =
	(Float) GJP.VolNodeSites () * n_gp * (2 * 1440) * MSCG.iter;
      if (F5D ())
	total_flops *= Ls;
      print_flops (fname, "MultiCG(F)", total_flops, time2 - time1);
#endif
      time1 = time2;
      for (int i = 0; i < Nshift; i++)
	precisionChange (grid_rb_out[i], grid_rb_out_f[i]);
#else
// Doesn't work for WilsonTM. Why??
      time2 = dclock ();
      print_time (fname, "setup()", time2 - time1);
      time1 = time2;
      Grid::MixedPrecisionConjugateGradientMultiShift < LATTICE_FERMION,
	LATTICE_FERMION_F > MSCG (F_RB_GRID_F, MdagMF, MdagMD, 20, shifts);
      MSCG.MaxOuterIterations = cg_arg[0]->max_num_iter;
      MSCG (grid_rb_in, grid_rb_out);
      time2 = dclock ();
      Float total_flops =
	(Float) GJP.VolNodeSites () * n_gp * (2 * 1440) * MSCG.iter;
      if (F5D ())
	total_flops *= Ls;
      print_flops (fname, "MultiCG(Mixed)", total_flops, time2 - time1);
      time1 = time2;
#endif



      for (int i = 0; i < Nshift; i++) {
//      precisionChange(grid_rb_out[i],grid_rb_out_f[i]);
#if 0
	Grid::MixedPrecisionConjugateGradient < LATTICE_FERMION,
	  LATTICE_FERMION_F > CG_m (shifts.tolerances[i], maxinner, maxouter,
				    F_RB_GRID_F, MdagMF, MdagMD);
	RealD shift = shifts.poles[i];
	CG_m (grid_rb_in, grid_rb_out[i], &shift);
#else
#endif
      }
      time2 = dclock ();
      print_time (fname, "ShiftedCG", time2 - time1);
      time1 = time2;
    }

    zeroit (grid_rb_in);	//using it as the accumulated solution for SINGLE
    for (int i = 0; i < Nshift; i++) {
      if (type == MULTI) {
	setCheckerboard (grid_out, grid_rb_out[i]);
#if (defined IF_TM ) && ( defined USE_F_CLASS_WILSON_TM)
	ImportFermion (f_temp, grid_out, Odd);
	f_out[i]->FTimesV1PlusV2 (fac - 1., f_temp, f_temp, cps_temp.Size ());
#else
	ImportFermion (f_out[i], grid_out, Odd);
#endif
      } else {			// type == SINGLE
#define CPS_VAXPY
#ifdef CPS_VAXPY
	setCheckerboard (grid_out, grid_rb_out[i]);
	ImportFermion (f_temp, grid_out, Odd);
	f_out[0]->FTimesV1PlusV2 ((alpha[i] * fac), f_temp, f_out[0],
				  cps_temp.Size ());
#else
	Grid::axpy (grid_rb_in, (alpha[i] * fac), grid_rb_out[i], grid_rb_in);
#endif
      }
    }
#ifndef CPS_VAXPY
    if (type == SINGLE) {
      setCheckerboard (grid_out, grid_rb_in);
      ImportFermion (f_out[0], grid_out, Odd);
    }
#endif
    time2 = dclock ();
    print_time (fname, "Acc/Export", time2 - time1);
    time1 = time2;
    // print_flops(cname,fname,time2,-time_total,

#if 0
    if (1) {
      LATTICE_FERMION temp (FERM_GRID);
      Ddwf.M (grid_out, temp);
      temp = temp - grid_in;
      true_res = std::sqrt (norm2 (temp) / norm2 (grid_in));
      VRB.Debug (cname, fname, "true_res=%g\n", true_res);
    }
    ImportFermion (f_out, grid_out, Odd);
#endif
   return 0;
  }

  void FminResExt (Vector * sol, Vector * source, Vector ** sol_old,
		   Vector ** vm, int degree, CgArg * cg_arg, CnvFrmType cnv_frm)
  {
    size_t f_size = (size_t) GJP.VolNodeSites () * FsiteSize () / 2;
    if (GJP.Gparity ())
      f_size *= 2;
    sol->VecZero (f_size);

//    ERR.NotImplemented (cname, "FminResExt");
  }



  int FmatInv (Vector * f_out, Vector * f_in,
	       CgArg * cg_arg,
	       Float * true_res,
	       CnvFrmType cnv_frm, PreserveType prs_f_in, int if_dminus)
  {
    const char *fname ("FmatInv()");

    static int verbose = 1;		// turning off debug output after first time
    mass = cg_arg->mass;
    VRB.Result (cname, fname, "mass=%0.14g dminus=%d\n", mass, if_dminus);
    RealD M5 = GJP.DwfHeight ();
    ImportGauge ();
    std::vector < int >twists = SetTwist ();

    LATTICE_FERMION grid_in (FERM_GRID), grid_out (FERM_GRID);
    DIRAC::ImplParams params;
    SetParams (params);
    VRB.Debug (cname, fname,
	       "DIRAC Ddwf (*Umu(%p), FIVE_GRID (%p %p) * UGridD(%p), *UrbGridD(%p), mass(%g), MOB PARAMS)\n",
	       Umu, FGridD, FrbGridD, UGridD, UrbGridD, mass);

    DIRAC Ddwf (*Umu, FIVE_GRID * UGridD, *UrbGridD, mass MOB PARAMS);
    Grid::QCD::LatticeGaugeFieldF UmuF (UGridF);
    precisionChange (UmuF, *Umu);
    DIRAC_F DdwfF (UmuF, FIVE_GRID_F * UGridF, *UrbGridF, mass MOB PARAMS);

#ifndef IF_TM
    if (if_dminus == 1) {
#if 1
//was broken in Grid.
      ImportFermion (grid_out, f_in);
      Ddwf.Dminus (grid_out, grid_in);
#else
// Obviously only for mobius, needs more work for zMobius
      ImportFermion (grid_in, f_in);
      Ddwf.DW (grid_in, grid_out, Grid::QCD::DaggerNo);
      double mob_c = mob_b - 1.;
      axpy (grid_in, -mob_c, grid_out, grid_in);
#endif
    } else if (if_dminus == -1) {
      ImportFermion (grid_out, f_in);
      Ddwf.Dminus (grid_out, grid_in);
      ImportFermion (f_out, grid_in);
      return 0;
    } else
#endif
      ImportFermion (grid_in, f_in);

    char fname_eig_root_bc[1024];
    snprintf (fname_eig_root_bc, 1024, "%s", cg_arg->fname_eigen);
    VRB.Result (cname, fname, "fname=%s\n", fname_eig_root_bc);
    int n_fields = GJP.SnodeSites ();
    const size_t f_size_per_site = (size_t) FsiteSize () / GJP.SnodeSites () / 2;

    int neig = 1;
    if_defl = false;
    if ((cg_arg->Inverter == CG_LOWMODE_DEFL) ||
	(cg_arg->Inverter == LOWMODEAPPROX)) {
      if_defl = true;
      neig = cg_arg->neig;
    }

    std::cout << "neig= " << neig << std::endl;
    std::vector < Float > eval (neig);
    std::vector < LATTICE_FERMION_F > evec_f (neig, F_RB_GRID_F);
    Grid::Guesser < Float, LATTICE_FERMION_F > guesser (neig, eval, evec_f);

    if (if_defl) {
      LATTICE_FERMION_F grid_f (FERM_GRID_F), grid_f_rb (F_RB_GRID_F);
      // search for eigen cache
      EigenCache *ecache;
      if ((ecache =
	   EigenCacheListSearch ((char *) fname_eig_root_bc, neig)) == 0) {
	ERR.General (cname, fname,
		     "Eigenvector cache does not exist: neig %d name %s \n",
		     neig, fname_eig_root_bc);
      }
      Float mass = cg_arg->mass;
      EigenContainer eigcon (*this, (char *) fname_eig_root_bc, neig,
			     f_size_per_site / 2, n_fields, ecache);

      Float *eval_cps = eigcon.load_eval ();
     VRB.Result (cname, fname, "eval_cps=%p\n", eval_cps);
      for (int i = 0; i < neig; i++) {
	eval[i] = eval_cps[i];
	Vector *evec = eigcon.nev_load (i);
	VRB.Debug (cname, fname, "eval[%d]=%0.14e evec[%d]=%p\n", i, eval[i],i,evec);
	ImpexFermion < float, LATTICE_FERMION_F, SITE_FERMION_F > (evec, grid_f,
								   1, Odd,
								   NULL);
	pickCheckerboard (Odd, evec_f[i], grid_f);
      }


    }

    ImportFermion (grid_out, f_out);

    Float stp_cnd = cg_arg->stop_rsd;
//      stp_cnd = stp_cnd*stp_cnd*norm2(grid_in);

// should be fixed to use evecs gnenerated from Lanczos. Making it compile for now.
#ifndef TWOKAPPA
    {
      Grid::ConjugateGradient < LATTICE_FERMION > CG (stp_cnd,
						      cg_arg->max_num_iter);
      Grid::SchurRedBlackDiagMooeeSolve < LATTICE_FERMION > SchurSolver (CG);
      SchurSolver (Ddwf, grid_in, grid_out);
    }
#else
    {
      int maxinner = cg_arg->max_num_iter;
      int maxouter = 10;
      Grid::QCD::SchurDifferentiableDiagTwo < IMPL > MdagMD (Ddwf);
      Grid::QCD::SchurDifferentiableDiagTwo < IMPL_F > MdagMF (DdwfF);
      Grid::MixedPrecisionConjugateGradient < LATTICE_FERMION,
	LATTICE_FERMION_F > CG_m (stp_cnd, maxinner, maxouter, F_RB_GRID_F,
				  MdagMF, MdagMD);
	VRB.Result (cname, fname, "if_defl=%d verbose=%d\n",if_defl,verbose);
      if (if_defl) {
	if (verbose) {
	  for (int i = 0; i < neig; i++) {
	    LATTICE_FERMION_F tmp (F_RB_GRID_F);
	    MdagMF.HermOp (evec_f[i], tmp);
	    Grid::ComplexD zalph = innerProduct (evec_f[i], tmp);
//              std::cout <<" <lambda|A|lambda>= " <<zalph<<std::endl;
	    Grid::RealD alph = real (zalph);
	    tmp = tmp - eval[i] * evec_f[i];
	    VRB.Result (cname, fname, "evec %d: norm alpha eval residual= %e %e %e %e\n",
			i,norm2 (evec_f[i]), alph, eval[i],norm2 (tmp));
	  }
	  verbose=0;
        }
	CG_m.useGuesser (guesser);
      }
      Grid::SchurRedBlackDiagTwoMixed < LATTICE_FERMION > SchurSolver (CG_m);
      SchurSolver (Ddwf, grid_in, grid_out);
    }
#endif
    if (true_res) {
      LATTICE_FERMION temp (FERM_GRID);
      Ddwf.M (grid_out, temp);
      temp = temp - grid_in;
      *true_res = std::sqrt (norm2 (temp) / norm2 (grid_in));
      VRB.Result (cname, fname, "true_res=%g\n", *true_res);
    }
    ImportFermion (f_out, grid_out);
    return 0;
  }
  int FmatInv (Vector * f_out, Vector * f_in,
	       CgArg * cg_arg,
	       Float * true_res,
	       CnvFrmType cnv_frm = CNV_FRM_YES,
	       PreserveType prs_f_in = PRESERVE_YES) {
    return FmatInv (f_out, f_in, cg_arg, true_res, cnv_frm, prs_f_in, 1);
  }

  int FmatInvTest (Vector * f_out, Vector * f_in,
		   CgArg * cg_arg,
		   Float * true_res,
		   CnvFrmType cnv_frm = CNV_FRM_YES,
		   PreserveType prs_f_in = PRESERVE_YES) {
    return FmatInv (f_out, f_in, cg_arg, true_res, cnv_frm, prs_f_in, 0);
  }


  int FmatInv (Vector * f_out, Vector * f_in,
	       CgArg * cg_arg,
	       CnvFrmType cnv_frm = CNV_FRM_YES,
	       PreserveType prs_f_in = PRESERVE_YES) {
    return FmatInv (f_out, f_in, cg_arg, NULL, cnv_frm, prs_f_in, 1);
  }

  int FeigSolv (Vector ** f_eigenv, Float * lambda,
		Float * chirality, int *valid_eig,
		Float ** hsum,
		EigArg * eig_arg, CnvFrmType cnv_frm = CNV_FRM_YES) {
    ERR.NotImplemented (cname, "FeigSolv");
    return 0;
  }

  // It sets the pseudofermion field phi from frm1, frm2.
  Float SetPhi (Vector * phi, Vector * frm1, Vector * frm2,
		Float mass, DagType dag)
  {
    return SetPhi (phi, frm1, frm2, mass, -12345, dag);
  }

  Float SetPhi (Vector * phi, Vector * frm1, Vector * frm2,
		Float mass, Float epsilon, DagType dag)
  {
    const char *fname = "SetPhi(V*,V*,V*,F)";

    if (phi == 0)
      ERR.Pointer (cname, fname, "phi");

    if (frm1 == 0)
      ERR.Pointer (cname, fname, "frm1");

    MatPc (phi, frm1, mass, epsilon, dag);
    return FhamiltonNode (frm1, frm1);
  }


//Seem to have a minus sign relative to Fbfm. Putting in for now.
  void MatPc (Vector * f_out, Vector * f_in, Float mass, Float epsilon,
	      DagType dag)
  {
    const char *fname ("MatPc()");
#if 1
//      mass = cg_arg->mass;
    SetMass (mass);
    SetEpsilon (epsilon);
    PrintMass (fname);
    RealD M5 = GJP.DwfHeight ();
#ifdef IF_TM
    if (epsilon < -10000)
      ERR.General (cname, fname,
		   "Non-Gparity version called with Gparity lattice class\n");
#endif
    ImportGauge ();
    std::vector < int >twists = SetTwist ();


    LATTICE_FERMION grid_in (FERM_GRID), grid_out (FERM_GRID);
    LATTICE_FERMION grid_rb_in (F_RB_GRID), grid_rb_out (F_RB_GRID);

//      ImportFermion(f_out, grid_out);

    ImportFermion (grid_in, f_in, Odd);
    pickCheckerboard (Odd, grid_rb_in, grid_in);
    DIRAC::ImplParams params;
    SetParams (params);
    DIRAC Ddwf (*Umu, FIVE_GRID * UGridD, *UrbGridD, mass MOB PARAMS);
#ifndef TWOKAPPA
    Grid::QCD::SchurDifferentiableOperator < IMPL > Mpc (Ddwf);
#else
    Grid::QCD::SchurDifferentiableDiagTwo < IMPL > Mpc (Ddwf);
#endif
    if (dag == DAG_NO)
      Mpc.Mpc (grid_rb_in, grid_rb_out);
    else
      Mpc.MpcDag (grid_rb_in, grid_rb_out);
    setCheckerboard (grid_out, grid_rb_out);
    ImportFermion (f_out, grid_out, Odd);
#endif
  }

  void MatPc (Vector * out, Vector * in, Float mass, DagType dag)
  {
    MatPc (out, in, mass, -12345, dag);
  }

  ForceArg EvolveMomFforce (Matrix * mom, Vector * frm,
			    Float mass, Float epsilon, Float step_size)
  {
    LatVector cps_temp (FsiteSize () / 6 * n_gp, GJP.VolNodeSites () / 2);
    MatPc (cps_temp.Vec (), frm, mass, epsilon, DAG_NO);
    return EvolveMomFforce (mom, cps_temp.Vec (), frm, mass, epsilon,
			    -step_size);
  }

  ForceArg EvolveMomFforce (Matrix * mom, Vector * frm,
			    Float mass, Float step_size)
  {
    return EvolveMomFforce (mom, frm, mass, -12345, step_size);
  }
  // It evolves the canonical momentum mom by step_size
  // using the fermion force.

  ForceArg EvolveMomFforce (Matrix * mom, Vector * phi, Vector * eta,
			    Float mass, Float epsilon, Float step_size)
  {
    return EvolveMomFforceBase (mom, phi, eta, mass, epsilon, -step_size);
  }
  ForceArg EvolveMomFforce (Matrix * mom, Vector * phi, Vector * eta,
			    Float mass, Float step_size)
  {
    return EvolveMomFforceBase (mom, phi, eta, mass, -12345, -step_size);
  }

  // It evolves the canonical Momemtum mom:
  // mom += coef * (phi1^\dag e_i(M) \phi2 + \phi2^\dag e_i(M^\dag) \phi1)
  // note: this function does not exist in the base Lattice class.
  //CK: This function is not used, so I have not modified it for WilsonTM

  ForceArg EvolveMomFforceBase (Matrix * mom,
				Vector * phi1,
				Vector * phi2,
				Float mass_, Float epsilon, Float coef)
  {
    const char *fname ("EvolveMomFforceBase()");
    mass = mass_;
    VRB.Result (cname, fname, "mass=%0.14g\n", mass);
    RealD M5 = GJP.DwfHeight ();
    ForceArg Fdt;



    ImportGauge ();
    SetMass (mass);
    SetEpsilon (epsilon);
    PrintMass (fname);

    Grid::QCD::LatticeGaugeFieldD grid_mom (UGridD), dSdU (UGridD),
      force (UGridD);
    ImportGauge (&grid_mom, mom);

//      std::vector<int> twists = SetTwist();
    DIRAC::ImplParams params;
    SetParams (params);

    LATTICE_FERMION grid_phi (FERM_GRID), grid_Y (FERM_GRID);
    LATTICE_FERMION X (F_RB_GRID), Y (F_RB_GRID);
    ImportFermion (grid_phi, phi1, Odd);
    ImportFermion (grid_Y, phi2, Odd);
    pickCheckerboard (Odd, X, grid_phi);
    pickCheckerboard (Odd, Y, grid_Y);

    DIRAC DenOp (*Umu, FIVE_GRID * UGridD, *UrbGridD, mass MOB PARAMS);
    DIRAC NumOp (*Umu, FIVE_GRID * UGridD, *UrbGridD, 1. MOB PARAMS);

// not really used. Just to fill invoke TwoFlavourEvenOddRatioPseudoFermionAction
    Grid::ConjugateGradient < LATTICE_FERMION > CG (1e-8, 10000);

    Grid::QCD::TwoFlavourEvenOddRatioPseudoFermionAction < IMPL > quo (NumOp,
								       DenOp,
								       CG, CG);

    {
      DenOp.ImportGauge (*Umu);

      Grid::QCD::SchurDifferentiableOperator < IMPL > Mpc (DenOp);
#ifdef TWOKAPPA
      ERR.General (cname, fname,
		   "Not implemnted for DiagTwo Preconditioning yet\n");
//      MpDeriv is written yet.
//      Grid::QCD::SchurDifferentiableDiagTwo < IMPL > Mpc (DenOp);
#endif

      force = 0.;
      Mpc.MpcDeriv (force, X, Y);
      dSdU = force;
      force = 0.;
      Mpc.MpcDagDeriv (force, Y, X);
      dSdU = dSdU + force;


    }
    grid_mom += -coef * 2. * Ta (dSdU);

    ImportGauge (mom, &grid_mom);

//#ifdef IF_TM
#if 0
    {
      int pos[4], mu;
      for (pos[0] = 0; pos[0] < GJP.XnodeSites (); pos[0]++)
	for (pos[1] = 0; pos[1] < GJP.YnodeSites (); pos[1]++)
	  for (pos[2] = 0; pos[2] < GJP.ZnodeSites (); pos[2]++)
	    for (pos[3] = 0; pos[3] < GJP.TnodeSites (); pos[3]++)
	      for (mu = 0; mu < 4; mu++)
		if (GJP.NodeBc (mu) == BND_CND_GPARITY)
		  if (pos[mu] == GJP.NodeSites (mu) - 1)
		    if ((pos[0] + pos[1] + pos[2] + pos[3]) % 2 == 0) {
		      Matrix *tmp_p = GaugeField () + GsiteOffset (pos) + mu;
		      Matrix mtmp = *tmp_p;
		      mtmp *= -1.;
		      *tmp_p = mtmp;
		    }
    }
#endif

    return Fdt;
  }

  // It evolves the canonical Momemtum mom:
  // mom += coef * (phi1^\dag e_i(M) \phi2 + \phi2^\dag e_i(M^\dag) \phi1)
  // note: this function does not exist in the base Lattice class.
  ForceArg EvolveMomFforceBase (Matrix * mom,
				Vector * phi1,
				Vector * phi2, Float mass, Float coef)
  {
    return EvolveMomFforceBase (mom, phi1, phi2, mass, -12345, coef);
  }

  // It evolve the canonical momentum mom  by step_size
  // using the bosonic quotient force.

  ForceArg RHMC_EvolveMomFforce (Matrix * mom, Vector ** sol, int degree,
				 int isz, Float * alpha, Float mass,
				 Float epsilon, Float dt, Vector ** sol_d,
				 ForceMeasure force_measure)
  {
    const char *fname ("RHMC_EvolveMomFforce()");
    char *force_label = NULL;

    int g_size = GJP.VolNodeSites () * GsiteSize ();
    if (GJP.Gparity ())
      g_size *= 2;

    Matrix *mom_tmp;

    if (force_measure == FORCE_MEASURE_YES) {
      mom_tmp =
	(Matrix *) smalloc (g_size * sizeof (Float), cname, fname, "mom_tmp");
      ((Vector *) mom_tmp)->VecZero (g_size);
      force_label = new char[100];
    } else {
      mom_tmp = mom;
    }

    if (!UniqueID ()) {
      Float pvals[4];
      for (int ii = 0; ii < 4; ii++) {
	int off = 18 * ii + 2;
	pvals[ii] = ((Float *) mom)[off];
      }
      VRB.Debug (cname, fname,
		 "initial mom Px(0) = %e, Py(0) = %e, Pz(0) = %e, Pt(0) = %e\n",
		 pvals[0], pvals[1], pvals[2], pvals[3]);
    }

    for (int i = 0; i < degree; i++) {
      Float *sol_f = (Float *) sol[i];
      VRB.Debug (cname, fname,
		 "sol[%d] = %g mass=%g epsilon=%g alpha[%d]=%g dt=%g\n", i,
		 *sol_f, mass, epsilon, i, alpha[i], dt);
      ForceArg Fdt =
	EvolveMomFforce (mom_tmp, sol[i], mass, epsilon, alpha[i] * dt);

      if (force_measure == FORCE_MEASURE_YES) {
	sprintf (force_label, "Rational, mass = %e, pole = %d:", mass, i + isz);
	Fdt.print (dt, force_label);
      }

      if (!UniqueID ()) {
	Float pvals[4];
	for (int ii = 0; ii < 4; ii++) {
	  int off = 18 * ii + 2;
	  pvals[ii] = ((Float *) mom_tmp)[off];
	}
	VRB.Debug (cname, fname,
		   "mom_tmp after pole %d:  Px(0) = %e, Py(0) = %e, Pz(0) = %e, Pt(0) = %e\n",
		   i, pvals[0], pvals[1], pvals[2], pvals[3]);
      }
    }

    ForceArg ret;

    // If measuring the force, need to measure and then sum to mom
    if (force_measure == FORCE_MEASURE_YES) {
      ret.measure (mom_tmp);
      ret.glb_reduce ();

      fTimesV1PlusV2 ((IFloat *) mom, 1.0, (IFloat *) mom_tmp, (IFloat *) mom,
		      g_size);

      delete[]force_label;
      sfree (mom_tmp, cname, fname, "mom_tmp");
    }

    return ret;
  }

  ForceArg RHMC_EvolveMomFforce (Matrix * mom, Vector ** sol, int degree,
				 int isz, Float * alpha, Float mass, Float dt,
				 Vector ** sol_d, ForceMeasure measure)
  {
    return RHMC_EvolveMomFforce (mom, sol, degree, isz, alpha, mass, -12345, dt,
				 sol_d, measure);
  }

  // The fermion Hamiltonian of the node sublattice.
  // chi must be the solution of Cg with source phi.           
  Float FhamiltonNode (Vector * phi, Vector * chi)
  {
    const char *fname = "FhamiltonNode(V*, V*)";

    if (phi == 0)
      ERR.Pointer (cname, fname, "phi");
    if (chi == 0)
      ERR.Pointer (cname, fname, "chi");

    size_t f_size = (size_t) GJP.VolNodeSites () * FsiteSize () / 2;
    if (GJP.Gparity ())
      f_size *= 2;

    // Sum accross s nodes is not necessary for MDWF since the library
    // does not allow lattice splitting in s direction.
    return phi->ReDotProductNode (chi, f_size);
  }

#if 0
  // doing nothing for now
  // Convert fermion field f_field from -> to
  void Fconvert (Vector * f_field, StrOrdType to, StrOrdType from)
  {
  }
#endif

  Float BhamiltonNode (Vector * boson, Float mass, Float epsilon)
  {
    ERR.NotImplemented (cname, "BhamiltonNode");
    return 0.;
  }
  Float BhamiltonNode (Vector * boson, Float mass)
  {
    ERR.NotImplemented (cname, "BhamiltonNode");
    return 0.;
  }
  // The boson Hamiltonian of the node sublattice
  int SpinComponents () const
  {
    return 4;
  }

  int ExactFlavors () const
  {
    return 2;
  }

  //!< Method to ensure bosonic force works (does nothing for Wilson
  //!< theories.
  void BforceVector (Vector * in, CgArg * cg_arg)
  {
    ERR.NotImplemented (cname, "BforceVec");
  }

#if 0
  int F5D () const
  {
#ifdef IF_FIVE_D
    return 1;
#else
    return 0;
#endif
  }
#endif

#ifndef IF_TM
#include "fgrid_lanczos.h.inc"
#endif

#ifdef GRID_MADWF
#include "fgrid_madwf.h.inc"
#endif


};
#if 1


//------------------------------------------------------------------
#define GFCLASS_NAME XSTR(GFCLASS(Gnone))
class GFCLASS (Gnone)
:public virtual Lattice,
  public virtual Gnone, public virtual FGRID, public virtual FgridBase
{
private:
  const char *cname;		// Class name.

public:
GFCLASS (Gnone) (FgridParams & params):cname (GFCLASS_NAME), FGRID (params),
    FgridBase (params) {
    const char *fname = GFCLASS_NAME "()";
    VRB.Func (cname, fname);
  }

  ~GFCLASS (Gnone) () {
    const char *fname = "~" GFCLASS_NAME "()";
    VRB.Func (cname, fname);
  }

};

#undef GFCLASS_NAME
#endif
CPS_END_NAMESPACE
#undef LATTICE_FERMION
#undef LATTICE_FERMION_F
#undef FIVE_GRID
#undef FIVE_GRID_F
#undef FERM_GRID
#undef FERM_GRID_F
#undef F_RB_GRID
#undef F_RB_GRID_F
#undef GRID_GPARITY
#undef IF_FIVE_D
#undef IF_TM
#undef FGRID
#undef CLASS_NAME
#undef DIRAC
#undef DIRAC_F
#undef MOB
#undef SITE_FERMION
#undef SITE_FERMION_F
#undef IMPL
#undef IMPL_F
#undef PARAMS
#undef GP
#undef GRID_GPARITY
