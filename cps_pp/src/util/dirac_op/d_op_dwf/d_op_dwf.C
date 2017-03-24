#include <config.h>
#include <stdio.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOpDwf class methods.

*/
//------------------------------------------------------------------
//
// d_op_dwf.C
//
// DiracOpDwf is derived from the DiracOp base class. 
// DiracOpDwf is the front end for a library that contains
// all Dirac operators associated with Dwf fermions.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/wilson.h>
#include <util/time_cps.h>
#include <util/dwf.h>
//#include <mem/p2v.h>
#include <comms/glb.h>
#ifdef USE_BLAS
#include <util/qblas_extend.h>
#endif

#ifdef USE_CG_DWF_WRAPPER
#include "cps_cg_dwf.h"
int CgDwfWrapper::lat_allocated=0;
const char *CgDwfWrapper::cname = "CPSCgDwf";
LatMatrix *CgDwfWrapper::Plus[4];
LatMatrix *CgDwfWrapper::Minus[4];
#endif

CPS_START_NAMESPACE



//------------------------------------------------------------------
/*!
  Only one instance of this class is allowed to be in existence at
  any time.
  \param latt The lattice on which this Dirac operator is defined
  \param f_field_out A (pointer to) a spin-colour field (optionally). 
  \param f_field_in A (pointer to) a spin-colour field (optionally).
  \param arg Parameters for the solver.
  \param cnv_frm_flag Whether the lattice fields should be converted to
  to a new storage order appropriate for the type of fermion action.
  If this is ::CNV_FRM_NO, then just the gauge field is converted.
  If this is ::CNV_FRM_YES, then the fields \a f_field_out and \a f_field_in
  are also converted: This assumes they are initially in the same order as
  the gauge field.
 */
//------------------------------------------------------------------
static  Matrix *new_gauge_field;
static  Matrix *old_gauge_field;
DiracOpDwf::DiracOpDwf(Lattice & latt,
			     Vector *f_field_out,
			     Vector *f_field_in,
			     CgArg *arg,
			     CnvFrmType cnv_frm_flg) :
			     DiracOpWilsonTypes(latt, 
						f_field_out,
						f_field_in, 
						arg,
						cnv_frm_flg)
{
  cname = "DiracOpDwf";
  char *fname = "DiracOpDwf(L&,V*,V*,CgArg*,CnvFrmType)";
  VRB.Func(cname,fname);


  //----------------------------------------------------------------
  // Do the necessary conversions
  //----------------------------------------------------------------
#define PROFILE
#ifdef PROFILE
  Float time = -dclock();
#endif
  if(cnv_frm == CNV_FRM_YES)
    lat.Convert(WILSON, f_out, f_in);
  else if(cnv_frm == CNV_FRM_NO)
    lat.Convert(WILSON);
#ifdef PROFILE
  time += dclock();
  print_flops("lattice","Convert()",0,time);
#endif

  //----------------------------------------------------------------
  // Initialize parameters
  //----------------------------------------------------------------
  DiracArg(dirac_arg);

  //----------------------------------------------------------------
  // Initialize the pointer to the initialized Dwf structure
  // (the structure has been initialized by the Lattice::Fdwf
  // constructor.
  //----------------------------------------------------------------
  dwf_lib_arg = lat.FdiracOpInitPtr();

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // call to dwf_init but is done here again in case the user
  // has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the Fdwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );
//  printf("new_gauge_field=%p size = %x \n",new_gauge_field,sizeof(Matrix)*GJP.VolNodeSites()*4);

}


//------------------------------------------------------------------
/*!
  If the storage order of any fields was changed by the constructor
  then they are changed to the canonical order by the destructor.
*/
//------------------------------------------------------------------
DiracOpDwf::~DiracOpDwf() {
  char *fname = "~DiracOpDwf()";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Do the necessary conversions
  //----------------------------------------------------------------
#define PROFILE
#ifdef PROFILE
  Float time = -dclock();
#endif
  if(cnv_frm == CNV_FRM_YES)
    lat.Convert(CANONICAL, f_out, f_in);
  else if(cnv_frm == CNV_FRM_NO)
    lat.Convert(CANONICAL);
#ifdef PROFILE
  time += dclock();
  print_flops("lattice","Convert()",0,time);
#endif

}


//------------------------------------------------------------------
// DiracArg(CgArg *arg):
// It sets the dirac_arg pointer to arg and initializes
// mass.
//------------------------------------------------------------------
void DiracOpDwf::DiracArg(CgArg *arg){
  dirac_arg = arg;
  mass = dirac_arg->mass;
}


//------------------------------------------------------------------
// MatPcDagMatPc :
// MatPcDagMatPc is the fermion matrix that appears in the HMC 
// evolution. It is a Hermitian matrix.
// The in, out fields are defined on the checkerboard lattice.
// If dot_prd is not 0 then the dot product (on node)
// <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.
//------------------------------------------------------------------
void DiracOpDwf::MatPcDagMatPc(Vector *out, 
			       Vector *in, 
			       Float *dot_prd){

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

//  printf("kappa:%e \n",dwf_arg->dwf_kappa);
  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  dwf_mdagm(out, 
	    gauge_field, 
	    in, 
	    dot_prd,
	    mass,
	    (Dwf *) dwf_lib_arg);

}


//------------------------------------------------------------------
// Dslash(Vector *out, Vector *in, ChkbType cb, DagType dag) :
// Dslash is the derivative part of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
// cb is not used.
//------------------------------------------------------------------
void DiracOpDwf::Dslash(Vector *out, 
			Vector *in, 
			ChkbType cb, 
			DagType dag) {

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  dwf_dslash(out, 
	     gauge_field, 
	     in, 
	     mass,
	     cb,
	     dag,
	     (Dwf *) dwf_lib_arg);
}

//------------------------------------------------------------------
/*!
  The preconditioned matrix connects sites of odd parity.
  \param out The resulting vector.
  \param in The vector to be multiplied.
*/
//------------------------------------------------------------------
void DiracOpDwf::MatPc(Vector *out, Vector *in) {  

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  dwf_m(out, 
	gauge_field, 
	in, 
	mass,
	(Dwf *) dwf_lib_arg);
}

//------------------------------------------------------------------
/*!
  The preconditioned matrix connects sites of odd parity.
  \param out The resulting vector.
  \param in The vector to be multiplied.
*/
//------------------------------------------------------------------
void DiracOpDwf::MatPcDag(Vector *out, Vector *in) {

  char *fname = "MatPcDag(*V,*V)";
  VRB.Func(fname,cname);
  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  IFloat *tmp = (IFloat *)in;
  dwf_mdag(out, 
	   gauge_field, 
	   in, 
	   mass,
	   (Dwf *) dwf_lib_arg);
  tmp = (IFloat *)out;
}

#if defined (USE_CG_DWF_WRAPPER)|| !defined (USE_CG_DWF)  || (TARGET != NOARCH)
//------------------------------------------------------------------
// int MatInv(Vector *out, Vector *in, 
//            Float *true_res, PreserveType prs_in);
// The inverse of the unconditioned Dirac Operator 
// using Conjugate gradient.
// If true_res !=0 the value of the true residual is returned
// in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// prs_in is used to specify if the source
// in should be preserved or not. If not the memory usage
// is less by half the size of a fermion vector.
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int DiracOpDwf::MatInv(Vector *out, 
		       Vector *in, 
		       Float *true_res,
		       PreserveType prs_in) {
  char *fname = "MatInv(V*,V*,F*)";
  VRB.Func(cname,fname);
//  VRB.Result(cname,fname,"Not using cg-dwf");

#define PROFILE
#ifdef PROFILE
  Float time = -dclock();
#endif
  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  Vector *temp2;
  unsigned long long temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;


  Vector *temp = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp == 0) ERR.Pointer(cname, fname, "temp");
  VRB.Smalloc(cname,fname, "temp", temp, temp_size * sizeof(Float));

  if(prs_in == PRESERVE_YES){
    temp2 = (Vector *) smalloc(temp_size * sizeof(Float));
    if (temp2 == 0) ERR.Pointer(cname, fname, "temp2");
    VRB.Smalloc(cname,fname, "temp2", temp2, temp_size * sizeof(Float));
  }

  // points to the even part of fermion source 
  Vector *even_in = (Vector *) ( (IFloat *) in + temp_size );

  // points to the even part of fermion solution
  Vector *even_out = (Vector *) ( (IFloat *) out + temp_size );

//  printf("MatInv : even : %e %e\n",even_in->NormSqNode(temp_size),even_out->NormSqNode(temp_size));

	
  Dslash(temp, even_in, CHKB_EVEN, DAG_NO);
//  printf("MatInv : even : Dslash : temp:%e even:%e\n",temp->NormSqNode(temp_size),even_in->NormSqNode(temp_size));

  fTimesV1PlusV2((IFloat *)temp, (IFloat) dwf_arg->dwf_kappa, (IFloat *)temp,
    (IFloat *)in, temp_size);

//  printf("MatInv : even : Dslash : temp:%e \n",temp->NormSqNode(temp_size));

  int iter;


  // save source
  if(prs_in == PRESERVE_YES){
    moveMem((IFloat *)temp2, (IFloat *)in, 
		temp_size * sizeof(IFloat) / sizeof(char));
  }

  switch (dirac_arg->Inverter) {
  case CG:
#ifdef USE_CG_DWF_WRAPPER
	{
  		CgDwfWrapper cg_dwf;
		cg_dwf.Init(lat);
		printf("lat.FsiteSize()=%d\n",lat.FsiteSize());
		cg_dwf.Inv(&lat,out,in,dirac_arg,temp_size);
	}
#endif
    MatPcDag(in, temp);
#ifdef PROFILE
  time += dclock();
  print_flops(fname,"Before InvCg()",0,time);
  time = -dclock();
#endif
#ifdef USE_QUDA
    iter = QudaInvert(out, in, true_res, 1);
#else
    iter = InvCg(out,in,true_res);
#endif
#ifdef PROFILE
  time += dclock();
  print_flops(fname,"MatPcDat+InvCg()",0,time);
  time = -dclock();
#endif
    break;
  case BICGSTAB:
#ifdef USE_QUDA
    iter = QudaInvert(out, temp, true_res, 0);
#else
    iter = BiCGstab(out,temp,0.0,dirac_arg->bicgstab_n,true_res);
#endif
    break;
  case LOWMODEAPPROX :
    MatPcDag(in, temp);
    iter = InvLowModeApprox(out,in, dirac_arg->fname_eigen, dirac_arg->neig, true_res );
    break;
  case CG_LOWMODE_DEFL :
    MatPcDag(in, temp);
    InvLowModeApprox(out,in, dirac_arg->fname_eigen, dirac_arg->neig, true_res );
    iter = InvCg(out,in,true_res);
    break;
  default:
    ERR.General(cname,fname,"InverterType %d not implemented\n",
		dirac_arg->Inverter);
  }

  // restore source
  if(prs_in == PRESERVE_YES){
    moveMem((IFloat *)in, (IFloat *)temp2, 
		temp_size * sizeof(IFloat) / sizeof(char));
  }

  Dslash(temp, out, CHKB_ODD, DAG_NO);

  fTimesV1PlusV2((IFloat *)even_out, (IFloat) dwf_arg->dwf_kappa, (IFloat *)temp,
    (IFloat *) even_in, temp_size);

  VRB.Sfree(cname, fname, "temp", temp);
  sfree(temp);

  if(prs_in == PRESERVE_YES){
    VRB.Sfree(cname, fname, "temp2", temp2);
    sfree(temp2);
  }
#ifdef PROFILE
  time += dclock();
  print_flops(fname,"After InvCg()",0,time);
#endif

  return iter;
}
#endif


//------------------------------------------------------------------
// Overloaded function is same as original 
// but true_res=0.
//------------------------------------------------------------------
int DiracOpDwf::MatInv(Vector *out, Vector *in, PreserveType prs_in)
{ return MatInv(out, in, 0, prs_in); }


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in and out = f_out.
//------------------------------------------------------------------
int DiracOpDwf::MatInv(Float *true_res, PreserveType prs_in)
{ return MatInv(f_out, f_in, true_res, prs_in); }


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in, out = f_out, true_res=0.
//------------------------------------------------------------------
int DiracOpDwf::MatInv(PreserveType prs_in)
{ return MatInv(f_out, f_in, 0, prs_in); }


//------------------------------------------------------------------
// Mat(Vector *out, Vector *in) :
// Mat is the unpreconditioned fermion matrix.  
// Mat works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpDwf::Mat(Vector *out, Vector *in) {  
  char *fname = "Mat(V*,V*)";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  // points to the even part of fermion source 
  Vector *even_in = (Vector *) ( (IFloat *) in + temp_size );
  // points to the even part of fermion solution
  Vector *even_out = (Vector *) ( (IFloat *) out + temp_size );

  Dslash(out, even_in, CHKB_EVEN, DAG_NO);

  fTimesV1PlusV2((IFloat *)out, -((Dwf *)dwf_lib_arg)->dwf_kappa, 
		 (IFloat *)out, (IFloat *)in, temp_size);

  Dslash(even_out, in, CHKB_ODD, DAG_NO);
  
  fTimesV1PlusV2((IFloat *)even_out, -((Dwf *)dwf_lib_arg)->dwf_kappa, 
		 (IFloat *)even_out, (IFloat *)even_in, temp_size);

}


//------------------------------------------------------------------
// MatDag(Vector *out, Vector *in) :
// MatDag is the unpreconditioned fermion matrix.  
// MatDag works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpDwf::MatDag(Vector *out, Vector *in) {
  char *fname = "MatDag(V*,V*)";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  // points to the even part of fermion source 
  Vector *even_in = (Vector *) ( (IFloat *) in + temp_size );
  // points to the even part of fermion solution
  Vector *even_out = (Vector *) ( (IFloat *) out + temp_size );

  Dslash(out, even_in, CHKB_EVEN, DAG_YES);

  fTimesV1PlusV2((IFloat *)out, -((Dwf *) dwf_lib_arg)->dwf_kappa, 
		 (IFloat *)out, (IFloat *)in, temp_size);

  Dslash(even_out, in, CHKB_ODD, DAG_YES);
  
  fTimesV1PlusV2((IFloat *)even_out, -((Dwf *) dwf_lib_arg)->dwf_kappa,
		 (IFloat *)even_out, (IFloat *)even_in, temp_size);

}


//------------------------------------------------------------------
// MatHerm(Vector *out, Vector *in) :
// MatHerm is gamma5*R*Mat.
// MatHerm works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpDwf::MatHerm(Vector *out, Vector *in) {
  char *fname = "MatHerm(V*,V*)";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize();
  Vector *temp = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp == 0) 
    ERR.Pointer(cname, fname, "temp");
  VRB.Smalloc(cname,fname, "temp", temp, temp_size);

  Mat(out, in);
  lat.Freflex(temp, out);
  MultGamma(out, temp, 15, GJP.VolNodeSites()*GJP.SnodeSites());
  
  VRB.Sfree(cname, fname, "temp", temp);
  sfree(temp);

}


//------------------------------------------------------------------
/*!
  \pre This method is to be used when the instance of this object has been
  created with \a f_field_out and \a f_field_in pointers to spin-colour
  vectors defined over the whole lattice. 

  \param chi A spin-colour vector defined on odd parity lattice sites.

  \post The vector \a f_field_out is \f$ (1+D)\chi \f$

  and the vector \a f_field_in is \f$ (D^\dagger-\kappa^2 M)\chi \f$

  where \e M is the odd-even preconditioned fermion matrix connecting odd to
  odd parity sites and \e D is the hopping term connecting odd to
  even parity sites. Recall that \a chi is defined on odd sites only.
  The new vectors are in odd-even order.
*/
//------------------------------------------------------------------

void DiracOpDwf::CalcHmdForceVecs(Vector *chi)
{
  char *fname = "CalcHmdForceVecs(V*)" ;
  VRB.Func(cname,fname) ;

  if (f_out == 0)
    ERR.Pointer(cname, fname, "f_out") ;

  if (f_in == 0)
    ERR.Pointer(cname, fname, "f_in") ;

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------

//------------------------------------------------------------------
// f_out stores (chi,rho), f_in stores (psi,sigma)
//------------------------------------------------------------------

  Vector *chi_new, *rho, *psi, *sigma ;

  int f_size_cb = 12 * GJP.VolNodeSites() * GJP.SnodeSites() ;

  chi_new = f_out ;

  chi_new->CopyVec(chi, f_size_cb) ;

  psi = f_in ;

//  fprintf(stderr,"psi=%p chi=%p rho=%p sigma=%p\n",psi,chi,rho,sigma);
  MatPc(psi,chi) ;
//  fprintf(stderr,"MatPc\n");

  {
    Float kappa = ((Dwf *)dwf_lib_arg)->dwf_kappa ;
    psi->VecTimesEquFloat(-kappa*kappa,f_size_cb) ;
  }

  rho = (Vector *)((Float *)f_out + f_size_cb) ;

  Dslash(rho, chi, CHKB_ODD, DAG_NO) ;
//  fprintf(stderr,"Dslash\n");

  sigma = (Vector *)((Float *)f_in + f_size_cb) ;

  Dslash(sigma, psi, CHKB_ODD, DAG_YES) ;
//  fprintf(stderr,"Dslash\n");

  return ;
}

//------------------------------------------------------------------
// DiracOpGlbSum(Float *): 
// The global sum used by InvCg. If s_nodes = 1
// it is the usual global sum. If s_nodes > 1 it
// is the 5-dimensional globals sum glb_sum_five.
//------------------------------------------------------------------
void DiracOpDwf::DiracOpGlbSum(Float *float_p) {
//  if(GJP.Snodes() == 1) {
//    glb_sum(float_p);
//  }
//  else {
    glb_sum_five(float_p);
//  }
}


//#undef DEBUG_DWF_DSLASH
//#define DEBUG_DWF_DSLASH
#ifdef  DEBUG_DWF_DSLASH
#undef DEBUG_DWF_DSLASH
#define DEBUG_DWF_DSLASH(msg,a ...) do		\
    if( UniqueID()%100==0 )			\
      printf("[%05d] " msg, UniqueID() ,##a);	\
  while(0);

#else
#define time_elapse() 0
#define DEBUG_DWF_DSLASH(msg,a ...) {}
#endif


#ifndef USE_BLAS
#define MOVE_FLOAT( pa, pb, n )  moveFloat(pa, pb, n) 
#define VEC_TIMESEQU_FLOAT(py, fact, n ) vecTimesEquFloat( py, fact, n)
#define AXPY(n, fact, px, py)  fTimesV1PlusV2(py, fact, px, py, n)
#else
#define MOVE_FLOAT( pa, pb, n )  cblas_dcopy(n, pb, 1, pa, 1)
#define VEC_TIMESEQU_FLOAT(py, fact, n ) cblas_dscal( n,  fact, py,1 )
#define AXPY(n, fact, px, py)  cblas_daxpy(n, fact, px,1,py,1)
#endif


//------------------------------------------------------------------
// RitzMat(Vector *out, Vector *in) :
// RitzMat is the base operator used in in Ritz.
// RitzMat works on the full or half lattice.
// The in, out fields are defined on the full or half lattice.
//------------------------------------------------------------------
void DiracOpDwf::RitzMat(Vector *out, Vector *in) {
  char *fname = "RitzMat(V*,V*)";
  VRB.Func(cname,fname);

  //printf("single ritzmat %d\n",dirac_arg->RitzMatOper);
  switch(dirac_arg->RitzMatOper)
    {
      // Now always call MatPcDagMatPcShift even for MATPCDAG_MATPC case
      //
      //case MATPCDAG_MATPC:
    case MATPCDAG_MATPC:
      MatPcDagMatPc(out, in);
      break;
    //case MATPCDAG_MATPC_SHIFT:
      //MatPcDagMatPcShift(out, in);
      //break;
    case MATPC_HERM :
      MatPcHerm(out, in);
      break;

    case MAT_HERM:
    case MATDAG_MAT:
      MatDagMat(out, in);
      break;
    case NEG_MATPCDAG_MATPC:
      MatPcDagMatPc(out, in);
      out -> VecNegative(out, RitzLatSize());
      break;    
      
    case NEG_MATDAG_MAT:
      MatDagMat(out, in);
      out->VecNegative(out, RitzLatSize());
      break;
      
    default:
      ERR.General(cname,fname,"RitzMatOper %d not implemented\n",
		  dirac_arg->RitzMatOper);
    }

#if 0
  //debug
  const int size = RitzLatSize(); // this is number of Float, f_size .
  Complex deb = out->CompDotProductGlbSum(in, size);
  Float d_deb = in->NormSqGlbSum(size);
  printf("single Ritz %e %e %e\n", deb.real()/d_deb, deb.imag()/d_deb,d_deb);
  //debug
#endif
}



// PolynomialAccerelation
//
//  Q = [ -2 RitzMat + (alpha + beta) ] / [ alpha - beta ]
//
//  Output:  out =  T_n(Q) in
//
//  T_0 = 1,    T_1 = Q
//   T_{n+1}(Q) =  2 Q T_n(Q)  - T_{n-1}(Q)
//  
// Calling virtual RitzMat(V*,V*)
//
//   alpha = param[0]^2
//   beta = ( param[1] + fabs(eigend_shift) )^2
//
void DiracOpDwf::RitzMat(Vector *out, Vector *in,
			 MatrixPolynomialArg* cheby_arg) {
  char *fname = "RitzMat(V*,V*,I,F,F,F,V*,V*)";
  VRB.Func(cname,fname);

  //MatPcDagMatPcShift(out, in);
  //return;
  
  //DEBUG_DWF_DSLASH( "dummy %e\n", time_elapse() );
  double time_start=dclock();
    
  //printf("cheby ritzmat %d\n", dirac_arg->RitzMatOper);


  //debug const
  const Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  const Float shift = dirac_arg->eigen_shift;
  //printf("shift =%e\n",shift);exit(1);


  const int Npol = cheby_arg-> Npol;
  const int size = RitzLatSize(); // this is number of Float, f_size .
  //
  // Q = 2 / (alpha-beta)  RitzMat  -   (alpha+beta)/(alpha-beta)
  // 2 Q =   c1  (  c0 Ddag D  -   1 )
  //  c1 = 2 (alpha+beta)/(alpha-beta),   c0 =  2 / (alpha+beta)
  //
  const Float alpha = pow( cheby_arg-> params.params_val[0], 2);
  const Float beta  = pow( cheby_arg-> params.params_val[1] + fabs(shift), 2);
  //printf("alpha=%e beta=%e\n", alpha,beta);
  
  const Float c1 =   2.0*(alpha+beta)/(alpha-beta);
  const Float c0 =   2.0/(alpha+beta);

  Vector *tmp  = (Vector*)cheby_arg->tmp1;
  Vector *tmp2 = (Vector*)cheby_arg->tmp2;
  

  //  tmp2 =  T_0 v = v = in
  //tmp2 -> CopyVec(in, size);
  MOVE_FLOAT( (Float*)tmp2, (Float*)in, size );
  //  tmp =  T_1 v = Q v = Q in
  //  QV = 0.5* (2Q)V = 0.5 c1 ( c0 Ddag D - 1)

  RitzMat(tmp, in);
  
#if 0
  tmp->VecTimesEquFloat(c0, size);
  tmp->VecMinusEquVec(in,size);
  tmp->VecTimesEquFloat(0.5*c1, size);
  //  tmp =  0.5 c1 ( c0 DdagD in - in )
#else
  VEC_TIMESEQU_FLOAT((Float*)tmp,0.5*c1*c0, size);
  AXPY( size, -0.5*c1, (Float*)in, (Float*)tmp );
#endif
  
  // debug
  //out->CopyVec(tmp,size);
  //printf("cheby %f %f\n", alpha,beta);
  
  // loop over
  for(int i=2; i<=Npol; ++i){
#if 0
    // out = 2 Q tmp
    RitzMat(out, tmp);

    out->VecTimesEquFloat(c0, size);
    out->VecMinusEquVec(tmp,size);
    out->VecTimesEquFloat(c1, size);

    // out = out - tmp2
    out->VecMinusEquVec(tmp2, size);
#else
    // out = c1 (  c0 DagD tmp - tmp) -tmp2

    RitzMat(out, tmp);

    VEC_TIMESEQU_FLOAT((Float*)out, c1*c0, size);
    AXPY( size , - c1, (Float*)tmp, (Float*)out);
    AXPY( size,  -1.0, (Float*)tmp2, (Float*)out);
#endif

    
    if( i!=Npol) {
#if 0
      // tmp2 = tmp
      tmp2->CopyVec(tmp, size);
      // tmp = out
      tmp->CopyVec(out, size);
#else
      // tmp2 = tmp
      Vector* swap_tmp2 = tmp2;
      tmp2 = tmp;
      tmp = swap_tmp2;
      // tmp = out
      tmp->CopyVec(out, size);
#endif
    }
  }
 
  
  if(!UniqueID()) printf("mpoly tot = %e\n", dclock() - time_start);
  
#if 0
  //debug
  Complex deb = out->CompDotProductGlbSum(in, size);
  printf("debug %e %e %e\n", deb.real(), deb.imag(),in->NormSqGlbSum(size));
  //debug
#endif
  
}


//!! N.B. This overwrites contents of  dwf_arg->frm_tmp2
void DiracOpDwf::MatPcHerm(Vector *out, Vector *in) {
  char *fname = "MatPcHerm(V*,V*)";
  VRB.Func(cname,fname);
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  Vector* vtmp = (Vector*)(dwf_arg->frm_tmp1);
  
  MatPc(vtmp,in);
  ReflectAndMultGamma5( out, vtmp,  dwf_arg->vol_4d/2, dwf_arg->ls);
  
}

CPS_END_NAMESPACE
