#ifndef INCLUDED_FGRID_H
#define INCLUDED_FGRID_H
#include<stdlib.h>
#include<config.h>
#include<assert.h>
#ifdef USE_GRID
#include<util/lattice.h>
#include<util/time_cps.h>
#ifdef USE_BFM
#include<util/lattice/bfm_mixed_solver.h>
#endif
#include<util/multi_cg_controller.h>
#include<util/eigen_container.h>
#include<Grid/Grid.h>
//using namespace Grid;
//using namespace Grid::QCD;
#undef HAVE_HANDOPT

namespace Grid
{

  template < typename Float, class Field >
    class Guesser:public LinearFunction < Field >
  {
  public:
    int neig;
      std::vector < Float > &eval;
      std::vector < Field > &evec;
      Guesser (int n, std::vector < Float > &_eval,
	       std::vector < Field > &_evec)
    : neig (n), eval (_eval), evec (_evec)
    {
      assert (eval.size () >= neig);
      assert (evec.size () >= neig);
    }
    void operator  () (const Field & in, Field & out)
    {
      std::cout << GridLogMessage << "guesser() called " << std::endl;
      out = 0.;
      for (int i = 0; i < neig; i++) {
	Grid::ComplexD coef = innerProduct (evec[i], in);
	coef = coef / eval[i];
//	if (cps::VRB.IsActivated(cps::VERBOSE_DEBUG_LEVEL))
          std::cout<<GridLogMessage <<"eval coef norm(evec) "<<i<<" : "<<eval[i]<<" "<<coef<<" "<<norm2(evec[i])<< std::endl;
	out += coef * evec[i];
      }
      std::
	cout << GridLogMessage << "norm(out)  : " << norm2 (out) << std::endl;
    }
  };

}

CPS_START_NAMESPACE class FgridParams
{
public:
  Float mobius_scale;
  Float mobius_bmc;
    std::vector < std::complex < double >>omega;	//ZMobius
  Float epsilon;		//WilsonTM
    FgridParams ():mobius_scale (1.), mobius_bmc (1.)
  {
  }
   ~FgridParams ()
  {
  }
  void setZmobius (cps::Complex * bs, int ls)
  {
    omega.clear ();
    for (int i = 0; i < ls; i++) {
      std::complex < double >temp = 1. / (2. * bs[i] - 1.);
      VRB.Result ("FgridParams", "setZmobius", "bs[%d]=%g %g, omega=%g %g\n",
		  i, bs[i].real (), bs[i].imag (), i, temp.real (),
		  temp.imag ());
      omega.push_back (temp);
    }
  }
};

class FgridBase:public virtual Lattice, public virtual FgridParams,
  public virtual FwilsonTypes
{

  using RealD = Grid::RealD;
  using RealF = Grid::RealF;
public:
//  typedef enum EvenOdd
//  { Even, Odd, All } EvenOdd;
  const char *cname;
  static bool grid_initted;

protected:
  const int Nc = Grid::QCD::Nc;
  const int Nd = Grid::QCD::Nd;
  const int Ns = Grid::QCD::Ns;
  int n_gp;
    Grid::GridCartesian * UGridD;
    Grid::GridCartesian * UGridF;
    Grid::GridRedBlackCartesian * UrbGridD;
    Grid::GridRedBlackCartesian * UrbGridF;
    Grid::GridCartesian * FGridD;
    Grid::GridCartesian * FGridF;
    Grid::GridRedBlackCartesian * FrbGridF;
    Grid::GridRedBlackCartesian * FrbGridD;
    Grid::QCD::LatticeGaugeFieldD * Umu;
//      Grid::QCD::LatticeGaugeFieldF *Umu_f;
  int threads;
    std::vector < int >vol;	// global volume
    std::vector < int >nodes;
  RealD mass;
  RealD mob_b;			//Mobius
  RealD mob_c;			//Mobius
  RealD eps;			// WilsonTM 
    std::vector < std::complex < double >>omegas;	//ZMobius
  int Ls;

public:
  double get_mob_b ()
  {
    return mob_b;
  };

  Grid::GridCartesian * getFGrid () {
    return FGridD;
  }
  Grid::GridRedBlackCartesian * getFrbGrid () {
    return FrbGridD;
  }
  Grid::GridCartesian * getUGrid () {
    return UGridD;
  }
  Grid::GridCartesian * getUGridF () {
    return UGridF;
  }
  Grid::GridRedBlackCartesian * getUrbGrid () {
    return UrbGridD;
  }
  Grid::QCD::LatticeGaugeFieldD * getUmu () {
    return Umu;
  }
//      Grid::QCD::LatticeGaugeFieldF *getUmu_f(){return Umu_f;}
FgridBase (FgridParams & params):cname ("FgridBase"), vol (4, 1), nodes (4, 1), mass (1.),
    Ls (1) {
//,epsilon(0.),
    const char *fname ("FgridBase()");
    if (!grid_initted)
      Grid::Grid_init (GJP.argc_p (), GJP.argv_p ());
    grid_initted = true;
    *((FgridParams *) this) = params;
    eps = params.epsilon;
    omegas = params.omega;

//              VRB.Debug(cname,fname,"mobius_scale=%g\n",mobius_scale);
    mob_b = 0.5 * (mobius_scale + mobius_bmc);
    mob_c = mob_b - mobius_bmc;
    VRB.Func (cname, fname);
    if (!GJP.Gparity ()) {
//              ERR.General(cname,fname,"Only implemented for Grid with Gparity at the moment\n");
      n_gp = 1;
    } else
      n_gp = 2;
    VRB.Debug (cname, fname, "Grid initted\n");
    threads = Grid::GridThread::GetThreads ();
    for (int i = 0; i < 4; i++)
      vol[i] = GJP.NodeSites (i) * GJP.Nodes (i);;
    for (int i = 0; i < 4; i++)
      nodes[i] = GJP.Nodes (i);
    VRB.Result (cname, fname,
		"vol nodes Nd=%d Grid::vComplexD::Nsimd()=%d threads=%d omp_get_max_threads()=%d\n",
		Nd, Grid::vComplexD::Nsimd (), threads, omp_get_max_threads ());
    for (int i = 0; i < 4; i++)
      VRB.Debug (cname, fname, "%d %d \n", vol[i], nodes[i]);
    UGridD =
      Grid::QCD::SpaceTimeGrid::makeFourDimGrid (vol,
						 Grid::GridDefaultSimd (Nd, Grid:: vComplexD:: Nsimd ()),
						 nodes);
    UGridF =
      Grid::QCD::SpaceTimeGrid::makeFourDimGrid (vol,
						 Grid::GridDefaultSimd (Nd, Grid:: vComplexF:: Nsimd ()),
						 nodes);
    VRB.Debug (cname, fname, "UGridD=%p UGridF=%p\n", UGridD, UGridF);
    bool fail = false;
    for (int i = 0; i < 4; i++)
      if (GJP.NodeCoor (i) != UGridD->_processor_coor[i])
	fail = true;
    if (fail)
      for (int i = 0; i < 4; i++) {
	printf ("CPS: %d  pos[%d]=%d Grid: %d pos[%d]=%d\n", UniqueID (), i,
		GJP.NodeCoor (i), UGridD->_processor, i,
		UGridD->_processor_coor[i]);
      }
#ifdef HAVE_HANDOPT
    if (GJP.Gparity ())
      Grid::QCD::WilsonKernelsStatic::HandOpt = 0;	//Doesn't seem to be working with Gparity
    else
      Grid::QCD::WilsonKernelsStatic::HandOpt = 1;
#endif
    VRB.Debug (cname, fname, "UGrid.lSites()=%d\n", UGridD->lSites ());
    SetLs (GJP.SnodeSites ());
    UrbGridD = Grid::QCD::SpaceTimeGrid::makeFourDimRedBlackGrid (UGridD);
    UrbGridF = Grid::QCD::SpaceTimeGrid::makeFourDimRedBlackGrid (UGridF);
    VRB.Debug (cname, fname, "UrbGridD=%p UrbGridF=%p\n", UrbGridD, UrbGridF);
    FGridD = Grid::QCD::SpaceTimeGrid::makeFiveDimGrid (Ls, UGridD);
    FGridF = Grid::QCD::SpaceTimeGrid::makeFiveDimGrid (Ls, UGridF);
    VRB.Debug (cname, fname, "FGridD=%p FGridF=%p\n", FGridD, FGridF);
    VRB.Debug (cname, fname, "FGridD.lSites()=%d\n", FGridD->lSites ());
    FrbGridD = Grid::QCD::SpaceTimeGrid::makeFiveDimRedBlackGrid (Ls, UGridD);
    FrbGridF = Grid::QCD::SpaceTimeGrid::makeFiveDimRedBlackGrid (Ls, UGridF);
    VRB.Debug (cname, fname, "FrbGridD=%p FrbGridF=%p\n", FrbGridD, FrbGridF);
    Umu = new Grid::QCD::LatticeGaugeFieldD (UGridD);
//  Umu_f = new Grid::QCD::LatticeGaugeFieldF(UGrid_f);
    grid_initted = true;
    VRB.FuncEnd (cname, fname);
//#ifdef USE_QMP
#if 0
    fail = false;
    for (int t = 0; t < GJP.Tnodes (); t++)
      for (int z = 0; z < GJP.Znodes (); z++)
	for (int y = 0; y < GJP.Ynodes (); y++)
	  for (int x = 0; x < GJP.Xnodes (); x++) {
	    std::vector < int >node
	    {
	    x, y, z, t};
	    int cps_rank = QMP_get_node_number_from (&node[0]);	//is a MPI_COMM_WORLD rank
	    int grid_rank = UGridD->RankFromProcessorCoor (node);	//is an MPI_Cart rank. However this MPI_Cart is drawn from MPI_COMM_WORLD and so the rank mapping to physical processors should be the same. However check below
	    int fail = 0;
	    if (UGridD->_processor == grid_rank) {
	      int world_rank;
	      MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);
	      if (world_rank != UGridD->_processor)
		fail = 1;
	    }
	    QMP_status_t ierr = QMP_sum_int (&fail);
	    if (ierr != QMP_SUCCESS)
	      ERR.General ("FgridBase", "FgridBase", "Rank check sum failed\n");
	    if (fail != 0)
	      ERR.General ("FgridBase", "FgridBase",
			   "Grid MPI_Cart rank does not align with MPI_COMM_WORLD rank\n");

	    if (cps_rank != grid_rank) {
	      if (!UniqueID ()) {
		std::cout << "Error in FgridBase constructor: node (" << node[0]
		  << "," << node[1] << "," << node[2] << "," << node[3] <<
		  ") maps to different MPI ranks for Grid " << grid_rank <<
		  " and CPS " << cps_rank << std::endl;
		std::cout.flush ();
	      }
	      fail = true;
	    }
	  }
    if (fail)
      exit (0);
#endif

  }
  void ResetParams (FgridParams & params)
  {
    this->mobius_scale = params.mobius_scale;
    this->mobius_bmc = params.mobius_bmc;
    this->eps = params.epsilon;
    this->omegas = params.omega;
  }
  virtual ~ FgridBase (void)
  {
    if (Umu)
      delete Umu;
//              if(Umu_f) delete Umu_f;
    delete UGridD;
    delete UGridF;
    delete UrbGridD;
    delete UrbGridF;
    delete FGridD;
    delete FGridF;
    delete FrbGridD;
    delete FrbGridF;
//              BondCond();
//              Grid_finalize();
  }
  int SetLs (int _Ls)
  {
    Ls = _Ls;
    return Ls;
  }
  Float SetMass (Float _mass)
  {
    mass = _mass;
    return mass;
  }
  Float SetEpsilon (Float _epsilon)
  {
    eps = _epsilon;
    return eps;
  }


  void ImportGauge (Grid::QCD::LatticeGaugeFieldD * grid_lat, Matrix * mom)
  {
    ImpexGauge (grid_lat, NULL, mom, 1);
  }

  void ImportGauge (Matrix * mom, Grid::QCD::LatticeGaugeFieldD * grid_lat)
  {
    ImpexGauge (grid_lat, NULL, mom, 0);
  }

  void ImportGauge ()
  {
    ImpexGauge (Umu, NULL, NULL, 1);
  }

  void ImpexGauge (Grid::QCD::LatticeGaugeFieldD * grid_lat,
		   Grid::QCD::LatticeGaugeFieldF * grid_lat_f, Matrix * mom,
		   int cps2grid)
  {

    BondCond ();
    Float *gauge = (Float *) mom;
    if (!mom)
      gauge = (Float *) GaugeField ();
//              if (!grid_lat)  grid_lat = Umu;
//              if (!grid_lat_f && cps2grid )  grid_lat_f = Umu_f;
    unsigned long vol;
    const char *fname = "ImpexGauge()";
    Grid::GridBase * grid = grid_lat->_grid;
    if (grid_lat->_grid->lSites () != (vol = GJP.VolNodeSites ()))
      ERR.General (cname, fname,
		   "numbers of grid(%d) and GJP(%d) does not match\n",
		   grid_lat->_grid->lSites (), vol);
    std::vector < int >grid_coor;
    Grid::QCD::LorentzColourMatrixD siteGrid;
    Grid::QCD::LorentzColourMatrixF siteGrid_f;
    for (int site = 0; site < vol; site++)
      for (int mu = 0; mu < 4; mu++) {
	if (cps2grid) {
	  for (int i = 0; i < Nc; i++)
	    for (int j = 0; j < Nc; j++) {
	      Float *cps = gauge + 18 * (site * 4 + mu) + 6 * j + 2 * i;
	      std::complex < double >elem (*cps, *(cps + 1));
	      siteGrid (mu) ()(j, i) = elem;
	      siteGrid_f (mu) ()(j, i) = elem;
//                              if (norm(elem)>0.01) printf("gauge[%d][%d][%d][%d] = %g %g\n",site,mu,i,j,elem.real(),elem.imag());
	    }
	  Grid::Lexicographic::CoorFromIndex (grid_coor, site,
					      grid->_ldimensions);
	  pokeLocalSite (siteGrid, *grid_lat, grid_coor);
	  if (grid_lat_f)
	    pokeLocalSite (siteGrid_f, *grid_lat_f, grid_coor);
	} else {
	  Grid::Lexicographic::CoorFromIndex (grid_coor, site,
					      grid->_ldimensions);
	  peekLocalSite (siteGrid, *grid_lat, grid_coor);
	  for (int i = 0; i < Nc; i++)
	    for (int j = 0; j < Nc; j++) {
	      std::complex < double >elem;
	      elem = siteGrid (mu) ()(j, i);
	      Float *cps = gauge + 18 * (site * 4 + mu) + 6 * j + 2 * i;
	      *cps = elem.real ();
	      *(cps + 1) = elem.imag ();
//i	if (norm(elem)>0.01) printf("gauge[%d][%d][%d][%d] = %g %g\n",site,mu,i,j,elem.real(),elem.imag());
	    }
	}
      }
    Float *f_tmp = (Float *) gauge;
    VRB.Debug (cname, fname, "mom=(%g %g)(%g %g)(%g %g)\n",
	       *f_tmp, *(f_tmp + 1), *(f_tmp + 2),
	       *(f_tmp + 3), *(f_tmp + 4), *(f_tmp + 5));
    BondCond ();
  }
  std::vector < int >SetTwist ()
  {
    std::vector < int >twists (Nd, 0);
    for (int i = 0; i < 3; i++) {
      twists[i] = (GJP.Bc (i) == BND_CND_GPARITY) ? 1 : 0;
      if (twists[i])
	VRB.Debug (cname, "SetTwist()", "gparity[%d]=1\n", i);
    }
    return twists;
  }



  FclassType Fclass () const
  {
    return F_CLASS_GRID;
  }
  // It returns the type of fermion class

  //! Multiplication of a lattice spin-colour vector by gamma_5.
//  void Gamma5(Vector *v_out, Vector *v_in, int num_sites);

#if 0
  int FsiteOffsetChkb (const int *x) const
  {
    ERR.NotImplemented (cname, "FsiteOffsetChkb");
  }
  // Sets the offsets for the fermion fields on a 
  // checkerboard. The fermion field storage order
  // is not the canonical one but it is particular
  // to the Dwf fermion type. x[i] is the 
  // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
//  int FsiteOffset(const int *x) const;
#endif


#if 1
  // doing nothing for now
  // Convert fermion field f_field from -> to
  void Fconvert (Vector * f_field, StrOrdType to, StrOrdType from)
  {
  }
#endif

  // The boson Hamiltonian of the node sublattice
  int SpinComponents () const
  {
    return 4;
  }

  int ExactFlavors () const
  {
    return 2;
  }
  virtual int FsiteSize () const = 0;

#if 1
  void Ffour2five (Vector * five, Vector * four, int s_u, int s_l, int Ncb = 2);
  //!< Transforms a 4-dimensional fermion field into a 5-dimensional field.
  /* The 5d field is zero */
  // The 5d field is zero
  // except for the upper two components (right chirality)
  // at s = s_u which are equal to the ones of the 4d field
  // and the lower two components (left chirality) 
  // at s_l, which are equal to the ones of the 4d field
  // For spread-out DWF s_u, s_l refer to the global
  // s coordinate i.e. their range is from 
  // 0 to [GJP.Snodes() * GJP.SnodeSites() - 1]

  void Ffive2four (Vector * four, Vector * five, int s_u, int s_l, int Ncb = 2);
  //!< Transforms a 5-dimensional fermion field into a 4-dimensional field.
  //The 4d field has
  // the upper two components (right chirality) equal to the
  // ones of the 5d field at s = s_u and the lower two 
  // components (left chirality) equal to the
  // ones of the 5d field at s = s_l, where s is the 
  // coordinate in the 5th direction.
  // For spread-out DWF s_u, s_l refer to the global
  // s coordinate i.e. their range is from 
  // 0 to [GJP.Snodes() * GJP.SnodeSites() - 1]
  // The same 4D field is generarted in all s node slices.
#endif

};

CPS_END_NAMESPACE
#define PASTER(x,y) x ## y
#define EVALUATOR(x,y) PASTER(x,y)
#define GFCLASS(class) EVALUATOR(class, FGRID )
#define XSTR(s) STR(s)
#define STR(s) #s
// Match to Fbfm twisted wilson instead of FwilsonTM
#undef USE_F_CLASS_WILSON_TM
#undef TWOKAPPA
#define GRID_GPARITY
#define IF_FIVE_D
#undef IF_TM
#define FGRID FgridGparityMobius
#define CLASS_NAME F_CLASS_GRID_GPARITY_MOBIUS
#define DIRAC Grid::QCD::GparityMobiusFermionD
#define DIRAC_F Grid::QCD::GparityMobiusFermionF
#define MOB	,M5,mob_b,mob_b-1.
#define IMPL Grid::QCD::GparityWilsonImplD
#define IMPL_F Grid::QCD::GparityWilsonImplF
#define SITE_FERMION Grid::QCD::iGparitySpinColourVector<Grid::ComplexD>
#define SITE_FERMION_F Grid::QCD::iGparitySpinColourVector<Grid::ComplexF>
#define PARAMS	,params
#define GP gp
#include "fgrid.h.inc"
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
#undef PARAMS
#undef GP
#define GRID_GPARITY
#undef IF_FIVE_D
#define IF_TM
#define FGRID FgridGparityWilsonTM
#define CLASS_NAME F_CLASS_GRID_GPARITY_WILSON_TM
#define DIRAC Grid::QCD::GparityWilsonTMFermionD
#define DIRAC_F Grid::QCD::GparityWilsonTMFermionF
#define MOB  ,eps
#define IMPL Grid::QCD::GparityWilsonImplD
#define IMPL_F Grid::QCD::GparityWilsonImplF
#define SITE_FERMION Grid::QCD::iGparitySpinColourVector<Grid::ComplexD>
#define SITE_FERMION_F Grid::QCD::iGparitySpinColourVector<Grid::ComplexF>
#define PARAMS	,params
#define GP gp
#include "fgrid.h.inc"
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
#undef PARAMS
#undef GP
#undef GRID_GPARITY
#define IF_FIVE_D
#undef IF_TM
#define FGRID FgridMobius
#define CLASS_NAME F_CLASS_GRID_MOBIUS
#define DIRAC Grid::QCD::MobiusFermionD
#define DIRAC_F Grid::QCD::MobiusFermionF
#define MOB	,M5,mob_b,mob_b-1.
#define SITE_FERMION Grid::QCD::iSpinColourVector<Grid::ComplexD>
#define SITE_FERMION_F Grid::QCD::iSpinColourVector<Grid::ComplexF>
#define IMPL Grid::QCD::WilsonImplD
#define IMPL_F Grid::QCD::WilsonImplF
#define PARAMS
#define GP
#include "fgrid.h.inc"
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
#undef PARAMS
#undef GP
#undef GRID_GPARITY
#define IF_FIVE_D
#define GRID_ZMOB
#define FGRID FgridZmobius
#define CLASS_NAME F_CLASS_GRID_ZMOBIUS
#define DIRAC Grid::QCD::ZMobiusFermionD
#define DIRAC_F Grid::QCD::ZMobiusFermionF
#define MOB	,M5,omegas,1.,0.
#define SITE_FERMION Grid::QCD::iSpinColourVector<Grid::ComplexD>
#define SITE_FERMION_F Grid::QCD::iSpinColourVector<Grid::ComplexF>
#define IMPL Grid::QCD::ZWilsonImplD
#define IMPL_F Grid::QCD::ZWilsonImplF
#define PARAMS
#define GP
// Using TwoKappa only for zMobius for now, really SYM2
//#define TWOKAPPA
//#define GRID_MADWF
#include "fgrid.h.inc"
#undef GRID_GPARITY
#undef IF_FIVE_D
#undef GRID_ZMOB
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
#undef TWOKAPPA
#undef GRID_MADWF
#undef IF_FIVE_D
#define IF_TM
#define FGRID FgridWilsonTM
#define CLASS_NAME F_CLASS_GRID_WILSON_TM
#define DIRAC Grid::QCD::WilsonTMFermionD
#define DIRAC_F Grid::QCD::WilsonTMFermionF
#define MOB  ,eps
#define IMPL Grid::QCD::WilsonImplD
#define IMPL_F Grid::QCD::WilsonImplF
#define SITE_FERMION Grid::QCD::iSpinColourVector<Grid::ComplexD>
#define SITE_FERMION_F Grid::QCD::iSpinColourVector<Grid::ComplexF>
#define PARAMS
#define GP
#include "fgrid.h.inc"
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
#undef PARAMS
#undef GP
#endif //#ifdef USE_GRID
#endif
