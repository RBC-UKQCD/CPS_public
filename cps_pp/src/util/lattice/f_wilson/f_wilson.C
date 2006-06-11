#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of Fwilson class.

  $Id: f_wilson.C,v 1.22 2006-06-11 05:35:06 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_wilson/f_wilson.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_wilson.C
//
// Fwilson is derived from FwilsonTypes and is relevant to
// wilson fermions
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/wilson.h>
#include <util/verbose.h>
#include <util/gjp.h>
#include <util/error.h>
#include <comms/scu.h>
#include <comms/glb.h>

#define BENCHMARK
#ifdef BENCHMARK
#include <util/qcdio.h>
#include <sys/time.h>
unsigned long WfmFlops;
#ifndef timersub
#define timersub(a, b, result)                                                \
  do {                                                                        \
    (result)->tv_sec = (a)->tv_sec - (b)->tv_sec;                             \
    (result)->tv_usec = (a)->tv_usec - (b)->tv_usec;                          \
    if ((result)->tv_usec < 0) {                                              \
      --(result)->tv_sec;                                                     \
      (result)->tv_usec += 1000000;                                           \
    }                                                                         \
  } while (0)
#endif
#endif

CPS_START_NAMESPACE

//------------------------------------------------------------------
// Initialize static variables.
//------------------------------------------------------------------

//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
Fwilson::Fwilson()
: FwilsonTypes()
{
  cname = "Fwilson";
  char *fname = "Fwilson()";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Check if anisotropy is present and exit since Fwilson has
  // not been tested for anisotropic lattices.
  //----------------------------------------------------------------
  if(GJP.XiBare() != 1 ||
     GJP.XiV()    != 1 ||
     GJP.XiVXi()  != 1   ){
    ERR.General(cname,fname,
    "XiBare=%g, XiV=%g, XiVXi=%g : Fwilson has not been tested with anisotropy\n",
                GJP.XiBare(), GJP.XiV(), GJP.XiVXi());
  }

  //----------------------------------------------------------------
  // Do initializations before the wilson library can be used
  // Initialization involve memory allocation.
  //----------------------------------------------------------------

  //  static Wilson wilson_struct;
  //  f_dirac_op_init_ptr = &wilson_struct;
  //  wilson_init((Wilson *) f_dirac_op_init_ptr);

  static Wilson wilson_struct;
  f_dirac_op_init_ptr = &wilson_struct;
  wilson_init((Wilson*) f_dirac_op_init_ptr);
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
Fwilson::~Fwilson()
{
  char *fname = "~Fwilson()";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Un-initialize the wilson library. Memory is set free here.
  //----------------------------------------------------------------
  wilson_end((Wilson *) f_dirac_op_init_ptr);
}


//------------------------------------------------------------------
// FclassType Fclass():
// It returns the type of fermion class.
//------------------------------------------------------------------
FclassType Fwilson::Fclass() const{
  return F_CLASS_WILSON;
}

int Fwilson::FsiteSize() const
{
  return 2 * Colors() * SpinComponents();  
  // re/im * colors * spin_components
}



//------------------------------------------------------------------
// int FchkbEvl() :
// returns 1 => The fermion fields in the evolution
//      or the CG that inverts the evolution matrix
//      are defined on a single checkerboard (half the 
//      lattice).
//------------------------------------------------------------------
int Fwilson::FchkbEvl() const
{
  return 1;
}


//------------------------------------------------------------------
// int FmatEvlInv(Vector *f_out, Vector *f_in, 
//                CgArg *cg_arg, 
//                Float *true_res,
//		  CnvFrmType cnv_frm = CNV_FRM_YES):
// It calculates f_out where A * f_out = f_in and
// A is the preconditioned fermion matrix that appears
// in the HMC evolution (even/odd  preconditioning 
// of [Dirac^dag Dirac]. The inversion is done
// with the conjugate gradient. cg_arg is the structure
// that contains all the control parameters, f_in is the
// fermion field source vector, f_out should be set to be
// the initial guess and on return is the solution.
// f_in and f_out are defined on a checkerboard.
// If true_res !=0 the value of the true residual is returned
// in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int Fwilson::FmatEvlInv(Vector *f_out, Vector *f_in, 
			CgArg *cg_arg, 
			Float *true_res,
			CnvFrmType cnv_frm)
{
  int iter;
  char *fname = "FmatEvlInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

  DiracOpWilson wilson(*this, f_out, f_in, cg_arg, cnv_frm);

  WfmFlops = 0;
  struct timeval t_start, t_stop;
  gettimeofday(&t_start,NULL);
  
  iter = wilson.InvCg(&(cg_arg->true_rsd));
  if (true_res) *true_res = cg_arg ->true_rsd;

  gettimeofday(&t_stop,NULL);
  timersub(&t_stop,&t_start,&t_start);
  double flops= (double)WfmFlops;
  double secs = t_start.tv_sec + 1.E-6 *t_start.tv_usec;
//  printf("Wilson solve: %d iteratations %d flops %f Mflops per node\n",
//	 iter,WfmFlops,flops/(secs*1000000) );

  
  // Return the number of iterations
  return iter;
}

//------------------------------------------------------------------
// int FmatEvlMInv(Vector *f_out, Vector *f_in, 
//                Float shift[], int Nshift, 
//                CgArg **cg_arg, Float *true_res,
//		  CnvFrmType cnv_frm = CNV_FRM_YES):
// It calculates f_out where (A + shift)* f_out = f_in and
// A is the fermion matrix that appears in the HMC 
// evolution ([Dirac^dag Dirac]) and shift is a real shift of the 
// fermion matrix, with Nshift such shifts. The inversion is done 
// with the multishift conjugate gradient. cg_arg is the structure
// that contains all the control parameters, f_in is the
// fermion field source vector, f_out is the array of solution 
// vectors, f_in and f_out are defined on a checkerboard.
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int Fwilson::FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
			 int Nshift, int isz, CgArg **cg_arg,
			 CnvFrmType cnv_frm, MultiShiftSolveType type, 
			 Float *alpha, Vector **f_out_d)
{
  char *fname = "FmatMInv(V*, V*, .....)";
  VRB.Func(cname,fname);

  int f_size = GJP.VolNodeSites() * FsiteSize() / (FchkbEvl()+1);
  Float dot = f_in -> NormSqGlbSum4D(f_size);

  Float *RsdCG = new Float[Nshift];
  for (int s=0; s<Nshift; s++) RsdCG[s] = cg_arg[s]->stop_rsd;

  //Fake the constructor
  DiracOpWilson wilson(*this, f_out[0], f_in, cg_arg[0], cnv_frm);

  int return_value = wilson.MInvCG(f_out,f_in,dot,shift,Nshift,isz,RsdCG,type,alpha);  

  for (int s=0; s<Nshift; s++) cg_arg[s]->true_rsd = RsdCG[s];
  delete[] RsdCG;
  return return_value;
}

//------------------------------------------------------------------
// Lattice class api to the chronological inverter
//------------------------------------------------------------------
void Fwilson::FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
			 Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm)
{

  char *fname = "FminResExt(V*, V*, V**, V**, int, CgArg *, CnvFrmType)";
  VRB.Func(cname,fname);
  
  DiracOpWilson wilson(*this, sol, source, cg_arg, cnv_frm);
  wilson.MinResExt(sol,source,sol_old,vm,degree);
  
}

//------------------------------------------------------------------
// int FmatInv(Vector *f_out, Vector *f_in, 
//             CgArg *cg_arg, 
//             Float *true_res,
//             CnvFrmType cnv_frm = CNV_FRM_YES,
//             PreserveType prs_f_in = PRESERVE_YES):
// It calculates f_out where A * f_out = f_in and
// A is the fermion matrix (Dirac operator). The inversion
// is done with the conjugate gradient. cg_arg is the 
// structure that contains all the control parameters, f_in 
// is the fermion field source vector, f_out should be set 
// to be the initial guess and on return is the solution.
// f_in and f_out are defined on the whole lattice.
// If true_res !=0 the value of the true residual is returned
// in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// cnv_frm is used to specify if f_in should be converted 
// from canonical to fermion order and f_out from fermion 
// to canonical. 
// prs_f_in is used to specify if the source
// f_in should be preserved or not. If not the memory usage
// is less by half the size of a fermion vector.
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int Fwilson::FmatInv(Vector *f_out, Vector *f_in, 
		     CgArg *cg_arg, 
		     Float *true_res,
		     CnvFrmType cnv_frm,
		     PreserveType prs_f_in)
{
  int iter;
  char *fname = "FmatInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

  DiracOpWilson wilson(*this, f_out, f_in, cg_arg, cnv_frm);
  
  iter = wilson.MatInv(true_res, prs_f_in);
  
  // Return the number of iterations
  return iter;
}



//------------------------------------------------------------------
// int FeigSolv(Vector **f_eigenv, Float *lambda, int *valid_eig,
//              EigArg *eig_arg, 
//              CnvFrmType cnv_frm = CNV_FRM_YES):
// It solve  A * f_eigenv = lambda * f_eigenv where
// A is the fermion matrix (Dirac operator). The solution
// is done with the Ritz algorithm. eig_arg is the 
// structure that contains all the control parameters, f_eigenv
// is the fermion field eigenvectors, lambda are the
// returned eigenvalues.
// f_eigenv is defined on the whole lattice.
// The function returns the total number of Ritz iterations.
//------------------------------------------------------------------
int Fwilson::FeigSolv(Vector **f_eigenv, Float *lambda,
		      Float *chirality, int *valid_eig,
		      Float **hsum,
		      EigArg *eig_arg, 
		      CnvFrmType cnv_frm)
{
  int iter;
  char *fname = "FeigSolv(EigArg*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);
  CgArg cg_arg;
  cg_arg.mass = eig_arg->mass;
  cg_arg.RitzMatOper = eig_arg->RitzMatOper;
  int N_eig = eig_arg->N_eig;
  int i;

  //=========================
  // convert fermion field
  //=========================


  for(i=0; i < N_eig; ++i)
    Fconvert(f_eigenv[i], WILSON, CANONICAL);



  //------------------------------------------------------------------
  //  we want both the eigenvalues of D_{hermitian} and
  //  D^{+}D.  To not change the arguments passed to RitzEig,
  //  we pass a float pointer which points to 2 * N_eig values
  //  and return both lambda and lambda^2 from RitzEig
  //------------------------------------------------------------------

  Float * lambda2 = (Float * ) smalloc (N_eig*2*sizeof(Float));
  if ( lambda2 == 0 ) ERR.Pointer(cname,fname, "lambda2");
  
  {
    DiracOpWilson wilson(*this, (Vector*) 0 , (Vector*) 0, &cg_arg, CNV_FRM_NO);
    iter = wilson.RitzEig(f_eigenv, lambda2, valid_eig, eig_arg);
  }


  for(i=0; i < N_eig; ++i)
    {
      Fconvert(f_eigenv[i], CANONICAL, WILSON);
    }

  /*
    the call to RitzEig returns a negative number if either the KS or CG maxes
    out, we wish to cope with this in alg_eig, so "pass it up". Clean up the
    storage order first in case we still want to use the eigenvectors as a
    guess.
  */
  if ( iter < 0 ) { return iter ; }


  // Compute chirality
  int f_size = (GJP.VolNodeSites() * FsiteSize());

  Vector* v1 = (Vector *)smalloc(f_size*sizeof(Float));
  if (v1 == 0)
    ERR.Pointer(cname, fname, "v1");
  VRB.Smalloc(cname, fname, "v1", v1, f_size*sizeof(Float));

  for(i=0; i < N_eig; ++i)
  {
    Gamma5(v1, f_eigenv[i], GJP.VolNodeSites());
    chirality[i] = f_eigenv[i]->ReDotProductGlbSum4D(v1, f_size);
  }

  VRB.Sfree(cname, fname, "v1", v1);
  sfree(v1);


  // rescale wilson eigenvalues to the convention  m + Dslash(U)
  Float factor = 4.0 + eig_arg->mass;
    
  
  FILE* fp=Fopen(eig_arg->fname,"a");
  for(i=0; i<N_eig; ++i)
    {
      lambda2[i] *= factor;	 		 //rescale eigenvalue
      lambda2[N_eig + i] *= ( factor * factor ); //rescale squared evalue
      lambda[i]=lambda2[i];                      //copy back
      
      //print out eigenvalue, eigenvalue^2, chirality 
      Fprintf(fp,"%d %g %g %g %d\n",i,
              (float)lambda2[i],
              (float)lambda2[N_eig + i],
	      (float)chirality[i],valid_eig[i]);
    }
  Fclose(fp);
  sfree(lambda2); 


  // Slice-sum the eigenvector density to make a 1D vector
  if (eig_arg->print_hsum)
    for(i=0; i < N_eig; ++i)
      f_eigenv[i]->NormSqArraySliceSum(hsum[i], FsiteSize(), eig_arg->hsum_dir);


  // The remaining part in QCDSP version are all about "downloading
  // eigenvectors", supposedly not applicable here.

  // Return the number of iterations
  return iter;
}


//------------------------------------------------------------------
// SetPhi(Vector *phi, Vector *frm1, Vector *frm2, Float mass,
//        DagType dag):
// It sets the pseudofermion field phi from frm1, frm2.
// Note that frm2 is not used.
// Modified - now returns the (trivial) value of the action
//------------------------------------------------------------------
Float Fwilson::SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
		      Float mass, DagType dag){
  char *fname = "SetPhi(V*,V*,V*,F)";
  VRB.Func(cname,fname);
  CgArg cg_arg;
  cg_arg.mass = mass;

  if (phi == 0)
    ERR.Pointer(cname,fname,"phi") ;

  if (frm1 == 0)
    ERR.Pointer(cname,fname,"frm1") ;

  DiracOpWilson wilson(*this, frm1, frm2, &cg_arg, CNV_FRM_NO) ;
  
  if (dag == DAG_YES) wilson.MatPcDag(phi, frm1) ;
  else wilson.MatPc(phi, frm1) ;

  return FhamiltonNode(frm1, frm1);
}


//------------------------------------------------------------------
// EvolveMomFforce(Matrix *mom, Vector *frm, Float mass, 
//                 Float dt):
// It evolves the canonical momentum mom by dt
// using the fermion force.
//------------------------------------------------------------------
ForceArg Fwilson::EvolveMomFforce(Matrix *mom, Vector *chi, 
			      Float mass, Float dt)
{
  char *fname = "EvolveMomFforce(M*,V*,F,F,F)";
  VRB.Func(cname,fname);

  Matrix *gauge = GaugeField() ;

  if (Colors() != 3)
    ERR.General(cname,fname,"Wrong nbr of colors.") ;

  if (SpinComponents() != 4)
    ERR.General(cname,fname,"Wrong nbr of spin comp.") ;

  if (mom == 0)
    ERR.Pointer(cname,fname,"mom") ;

  if (chi == 0)
    ERR.Pointer(cname,fname,"chi") ;

//------------------------------------------------------------------
// allocate space for two CANONICAL fermion fields.
//------------------------------------------------------------------
  int f_size = FsiteSize() * GJP.VolNodeSites() ;

  char *str_v1 = "v1" ;
  Vector *v1 = (Vector *)smalloc(f_size*sizeof(Float)) ;
  if (v1 == 0)
    ERR.Pointer(cname, fname, str_v1) ;
  VRB.Smalloc(cname, fname, str_v1, v1, f_size*sizeof(Float)) ;

  char *str_v2 = "v2" ;
  Vector *v2 = (Vector *)smalloc(f_size*sizeof(Float)) ;
  if (v2 == 0)
    ERR.Pointer(cname, fname, str_v2) ;
  VRB.Smalloc(cname, fname, str_v2, v2, f_size*sizeof(Float)) ;

//------------------------------------------------------------------
// allocate space for two CANONICAL fermion field on a site.
//------------------------------------------------------------------

  char *str_site_v1 = "site_v1";
  Float *site_v1 = (Float *)smalloc(FsiteSize()*sizeof(Float));
  if (site_v1 == 0)
    ERR.Pointer(cname, fname, str_site_v1) ;
  VRB.Smalloc(cname, fname, str_site_v1, site_v1,
    FsiteSize()*sizeof(Float)) ;

  char *str_site_v2 = "site_v2";
  Float *site_v2 = (Float *)smalloc(FsiteSize()*sizeof(Float));
  if (site_v2 == 0)
    ERR.Pointer(cname, fname, str_site_v2) ;
  VRB.Smalloc(cname, fname, str_site_v2, site_v2,
    FsiteSize()*sizeof(Float)) ;

  {
    CgArg cg_arg ;
    cg_arg.mass = mass ;

    DiracOpWilson wilson(*this, v1, v2, &cg_arg, CNV_FRM_YES) ;
    wilson.CalcHmdForceVecs(chi) ;
  }

  int x, y, z, t, lx, ly, lz, lt ;

  lx = GJP.XnodeSites() ;
  ly = GJP.YnodeSites() ;
  lz = GJP.ZnodeSites() ;
  lt = GJP.TnodeSites() ;

//------------------------------------------------------------------
// start by summing first over direction (mu) and then over site
// to allow SCU transfers to happen face-by-face in the outermost
// loop.
//------------------------------------------------------------------

  int mu ;

  Matrix tmp, f ;

  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;

  for (mu=0; mu<4; mu++) {
    for (t=0; t<lt; t++)
    for (z=0; z<lz; z++)
    for (y=0; y<ly; y++)
    for (x=0; x<lx; x++) {
      int gauge_offset = x+lx*(y+ly*(z+lz*t)) ;
      int vec_offset = FsiteSize()*gauge_offset ;
      gauge_offset = mu+4*gauge_offset ;

      Float *v1_plus_mu ;
      Float *v2_plus_mu ;
      int vec_plus_mu_offset = FsiteSize() ;

      Float coeff = -2.0 * dt ;

      switch (mu) {
        case 0 :
          vec_plus_mu_offset *= (x+1)%lx+lx*(y+ly*(z+lz*t)) ;
          if ((x+1) == lx) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;                        
            v2_plus_mu = site_v2 ;                        
            if (GJP.XnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
          break ;
        case 1 :
          vec_plus_mu_offset *= x+lx*((y+1)%ly+ly*(z+lz*t)) ;
          if ((y+1) == ly) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;                        
            v2_plus_mu = site_v2 ;                        
            if (GJP.YnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
          break ;
        case 2 :
          vec_plus_mu_offset *= x+lx*(y+ly*((z+1)%lz+lz*t)) ;
          if ((z+1) == lz) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;
            v2_plus_mu = site_v2 ;
            if (GJP.ZnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
          break ;
        case 3 :
          vec_plus_mu_offset *= x+lx*(y+ly*(z+lz*((t+1)%lt))) ;
          if ((t+1) == lt) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;
            v2_plus_mu = site_v2 ;
            if (GJP.TnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
      } // end switch mu

      sproj_tr[mu](   (IFloat *)&tmp,
                      (IFloat *)v1_plus_mu,
                      (IFloat *)v2+vec_offset, 1, 0, 0);

      sproj_tr[mu+4]( (IFloat *)&f,
                      (IFloat *)v2_plus_mu,
                      (IFloat *)v1+vec_offset, 1, 0, 0);

      tmp += f ;

      f.DotMEqual(*(gauge+gauge_offset), tmp) ;

      tmp.Dagger(f) ;

      f.TrLessAntiHermMatrix(tmp) ;

      f *= coeff ;

      *(mom+gauge_offset) += f ;
      Float norm = f.norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);
    }
  }

//------------------------------------------------------------------
// deallocate space for two CANONICAL fermion fields on a site.
//------------------------------------------------------------------
  VRB.Sfree(cname, fname, str_site_v2, site_v2) ;
  sfree(site_v2) ;

  VRB.Sfree(cname, fname, str_site_v1, site_v1) ;
  sfree(site_v1) ;

//------------------------------------------------------------------
// deallocate space for two CANONICAL fermion fields.
//------------------------------------------------------------------
  VRB.Sfree(cname, fname, str_v2, v2) ;
  sfree(v2) ;

  VRB.Sfree(cname, fname, str_v1, v1) ;
  sfree(v1) ;

  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();

  return ForceArg(L1, sqrt(L2), Linf);
}

ForceArg Fwilson::RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
				    int isz, Float *alpha, Float mass, 
				    Float dt, Vector **sol_d, 
				    ForceMeasure force_measure) {
  char *fname = "RHMC_EvolveMomFforce";
  char *force_label;

  ForceArg Fdt;
  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;

  int g_size = GJP.VolNodeSites() * GsiteSize();

  Matrix *mom_tmp;

  if (force_measure == FORCE_MEASURE_YES) {
    mom_tmp = (Matrix*)smalloc(g_size*sizeof(Float),cname, fname, "mom_tmp");
    ((Vector*)mom_tmp) -> VecZero(g_size);
    force_label = new char[100];
  } else {
    mom_tmp = mom;
  }

  for (int i=0; i<degree; i++) {
    ForceArg Fdt = EvolveMomFforce(mom_tmp,sol[i],mass,dt*alpha[i]);
    if (force_measure == FORCE_MEASURE_YES) {
      sprintf(force_label, "Rational, mass = %e, pole = %d:", mass, i+isz);
      Fdt.print(dt, force_label);
    }
  }

  // If measuring the force, need to measure and then sum to mom
  if (force_measure == FORCE_MEASURE_YES) {
    for (int i=0; i<g_size/18; i++) {
      Float norm = (mom_tmp+i)->norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);
    }
    glb_sum(&L1);
    glb_sum(&L2);
    glb_max(&Linf);

    L1 /= 4.0*GJP.VolSites();
    L2 /= 4.0*GJP.VolSites();

    fTimesV1PlusV2((IFloat*)mom, 1.0, (IFloat*)mom_tmp, (IFloat*)mom, g_size);

    delete[] force_label;
    sfree(mom_tmp, cname, fname, "mom_tmp");
  }

  return ForceArg(L1, sqrt(L2), Linf);
}

//------------------------------------------------------------------
// Float BhamiltonNode(Vector *boson, Float mass):
// The boson Hamiltonian of the node sublattice.
//------------------------------------------------------------------
Float Fwilson::BhamiltonNode(Vector *boson, Float mass){
  char *fname = "BhamiltonNode(V*,F)";
  VRB.Func(cname,fname);
  CgArg cg_arg;
  cg_arg.mass = mass;

  if (boson == 0)
    ERR.Pointer(cname,fname,"boson");

  int f_size = (GJP.VolNodeSites() * FsiteSize()) >> 1 ;

  Vector *bsn_tmp = (Vector *)
    smalloc(f_size*sizeof(Float));

  char *str_tmp = "bsn_tmp" ;

  if (bsn_tmp == 0)
    ERR.Pointer(cname,fname,str_tmp) ;

  VRB.Smalloc(cname,fname,str_tmp,bsn_tmp,f_size*sizeof(Float));

  DiracOpWilson wilson(*this, boson, bsn_tmp, &cg_arg, CNV_FRM_NO) ;

  wilson.MatPc(bsn_tmp,boson);

  Float ret_val = bsn_tmp->NormSqNode(f_size) ;

  VRB.Sfree(cname,fname,str_tmp,bsn_tmp);

  sfree(bsn_tmp) ;

  return ret_val;
}

ForceArg Fwilson::EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
		      Float mass, Float dt) {
  char *fname = "EvolveMomFforce(M*,V*,V*,F,F)";
  ERR.General(cname,fname,"Not Implemented\n");
  return ForceArg(0.0,0.0,0.0);
}

CPS_END_NAMESPACE
