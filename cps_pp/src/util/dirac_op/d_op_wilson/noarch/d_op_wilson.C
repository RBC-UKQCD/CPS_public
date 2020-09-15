#include <config.h>
#if 0
#include "../sse/d_op_wilson.C"
#else
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOpWilson class methods.

*/
//------------------------------------------------------------------
//
// d_op_wilson.C
//
// DiracOpWilson is derived from the DiracOp base class. 
// DiracOpWilson is the front end for a library that contains
// all Dirac operators associated with Wilson fermions.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/wilson.h>
#include <comms/glb.h>
CPS_START_NAMESPACE

bool DiracOpWilson::use_bfm = false;
#ifdef USE_BFM
bfmarg DiracOpWilson::bfm_arg;
#endif

//------------------------------------------------------------------
/*!
  Only one instance of this class is allowed to be in existence at
  any time.
  \param latt The lattice on which this Dirac operator is defined
  \param f_field_out A (pointer to) a spin-colour field (optionally). 
  \param f_field_in A (pointer to) a spin-colour field (optionally).
  \param arg Parameters for the solver.
  \param convert Whether the lattice fields should be converted to
  to a new storage order appropriate for the type of fermion action.
  If this is ::CNV_FRM_NO, then just the gauge field is converted.
  If this is ::CNV_FRM_YES, then the fields \a f_field_out and \a f_field_in
  are also converted: This assumes they are initially in the same order as
  the gauge field.
 */
//------------------------------------------------------------------
DiracOpWilson::DiracOpWilson(Lattice & latt,
			     Vector *f_field_out,
			     Vector *f_field_in,
			     CgArg *arg,
			     CnvFrmType convert) :
			     DiracOpWilsonTypes(latt, 
						f_field_out,
						f_field_in, 
						arg,
						convert)
{
  cname = "DiracOpWilson";
  char *fname = "DiracOpWilson(L&,V*,V*,CgArg*,CnvFrmType)";
  VRB.Func(cname,fname);

#ifdef USE_BFM
  if(use_bfm) {
    assert(bfm_arg.solver == WilsonFermion || bfm_arg.solver == WilsonTM);
    bevo.init(bfm_arg);

    assert(lat.StrOrd() == CANONICAL); //required by BondCond()
    Float *gauge = (Float *)(lat.GaugeField());
    latt.BondCond();
    bevo.cps_importGauge(gauge);
    latt.BondCond();
  }
#endif

  //----------------------------------------------------------------
  // Do the necessary conversions
  //----------------------------------------------------------------
  if(cnv_frm == CNV_FRM_YES)
    lat.Convert(WILSON, f_out, f_in);
  else if(cnv_frm == CNV_FRM_NO)
    lat.Convert(WILSON);

  //----------------------------------------------------------------
  // Initialize parameters
  //----------------------------------------------------------------
  DiracArg(dirac_arg);

  //----------------------------------------------------------------
  // Initialize the pointer to the initialized Wilson structure
  // (the structure has been initialized by the Lattice::Fwilson
  // constructor.
  //----------------------------------------------------------------
  wilson_lib_arg = lat.FdiracOpInitPtr();

}


//------------------------------------------------------------------
/*!
  If the storage order of any fields was changed by the constructor
  then they are changed to the canonical order by the destructor.
*/
//------------------------------------------------------------------
DiracOpWilson::~DiracOpWilson() {
  char *fname = "~DiracOpWilson()";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Do the necessary conversions
  //----------------------------------------------------------------
  if(cnv_frm == CNV_FRM_YES)
    lat.Convert(CANONICAL, f_out, f_in);
  else if(cnv_frm == CNV_FRM_NO)
    lat.Convert(CANONICAL);

#ifdef USE_BFM
  if(use_bfm) bevo.end();
#endif
}


//------------------------------------------------------------------
// DiracArg(CgArg *arg):
// It sets the dirac_arg pointer to arg and initializes
// kappa.
//------------------------------------------------------------------
  void DiracOpWilson::DiracArg(CgArg *arg){
    dirac_arg = arg;
    kappa = 1.0 / (2.0 * (dirac_arg->mass + 4.0));
  }


//------------------------------------------------------------------
// MatPcDagMatPc :
// MatPcDagMatPc is the fermion matrix that appears in the HMC 
// evolution. It is a Hermitian matrix where M is
// the even/odd preconditioned Dirac Operator matrix.        
// MatPcDagMatPc connects only odd-->odd sites.
// The in, out fields are defined on the odd checkerboard.
// If dot_prd is not 0 then the dot product (on node)
// <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.
//------------------------------------------------------------------
void DiracOpWilson::MatPcDagMatPc(Vector *out, 
					 Vector *in, 
					 Float *dot_prd){
  const char* fname = "MatPcDagMatPc()";

#ifdef USE_BFM
  if(use_bfm) {
    MatPcDagMatPc_BFM(out, in, dot_prd);
    return; 
  }
#endif

  wilson_mdagm((IFloat *)out, 
	       (IFloat *)gauge_field, 
	       (IFloat *)in, 
	       (IFloat *)dot_prd,
	       IFloat(kappa),
	       (Wilson *)wilson_lib_arg);
}


//------------------------------------------------------------------
// Dslash(Vector *out, Vector *in, ChkbType cb, DagType dag) :
// Dslash is the derivative part of the fermion matrix. 
// Dslash conects only odd-->even or even-->odd sites.
// The in, out fields are defined on a checkerboard.
// cb refers to the checkerboard of the in field.
// cb = 0/1 <--> even/odd checkerboard of in field.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//------------------------------------------------------------------
void DiracOpWilson::Dslash(Vector *out, 
			   Vector *in, 
			   ChkbType cb, 
			   DagType dag) {
  const char* fname = "Dslash()";
#ifdef USE_BFM
  if(use_bfm) ERR.NotImplemented(cname, fname);
#endif

  wilson_dslash((IFloat *)out, 
		(IFloat *)gauge_field, 
		(IFloat *)in, 
		int(cb),
		int(dag),
		(Wilson *)wilson_lib_arg);
}

//------------------------------------------------------------------
/*!
  The vectors are defined on odd parity lattice sites.
  \param out The resulting vector.
  \param in The vector to be multiplied.
*/
//------------------------------------------------------------------
void DiracOpWilson::MatPc(Vector *out, Vector *in) {  

#ifdef USE_BFM
  if(use_bfm) {
    MatPc_BFM(out, in, DAG_NO);
    return;
  }
#endif

  wilson_m((IFloat *)out, 
	   (IFloat *)gauge_field, 
	   (IFloat *)in, 
	   IFloat(kappa),
	   (Wilson *)wilson_lib_arg);
}

//------------------------------------------------------------------
/*!
  Multiplication by the odd-even preconditioned fermion matrix.
  The vectors are defined on odd parity lattice sites.
  \param out The resulting vector.
  \param in The vector to be multiplied.
*/
//------------------------------------------------------------------
void DiracOpWilson::MatPcDag(Vector *out, Vector *in) {
  
#ifdef USE_BFM
  if(use_bfm) {
    MatPc_BFM(out, in, DAG_YES);
    return;
  }
#endif
  
  wilson_mdag((IFloat *)out, 
	      (IFloat *)gauge_field, 
	      (IFloat *)in, 
	      IFloat(kappa),
	      (Wilson *)wilson_lib_arg);
}

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
int DiracOpWilson::MatInv(Vector *out, 
			  Vector *in, 
			  Float *true_res,
			  PreserveType prs_in) {
  char *fname = "MatInv(V*,V*,F*)";
  VRB.Func(cname,fname);
  
#ifdef USE_BFM
  if(use_bfm) ERR.NotImplemented(cname, fname);
#endif
  
  Vector *temp2 = NULL;

  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  if(GJP.Gparity()) temp_size *= 2;

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

  Dslash(temp, even_in, CHKB_EVEN, DAG_NO);

  fTimesV1PlusV2((IFloat *)temp, (IFloat) kappa, (IFloat *)temp,
    (IFloat *)in, temp_size);

  // save source
  if(prs_in == PRESERVE_YES){
    moveMem((IFloat *)temp2, (IFloat *)in, 
	    temp_size * sizeof(IFloat) / sizeof(char));
  }

  int iter;
  switch (dirac_arg->Inverter) {
  case CG:
    MatPcDag(in, temp);
    iter = InvCg(out,in,true_res);
    break;
  case BICGSTAB:
    iter = BiCGstab(out,temp,0.0,dirac_arg->bicgstab_n,true_res);
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

  fTimesV1PlusV2((IFloat *)even_out, (IFloat) kappa, (IFloat *)temp,
    (IFloat *) even_in, temp_size);

  VRB.Sfree(cname, fname, "temp", temp);
  sfree(temp);

  if(prs_in == PRESERVE_YES){
    VRB.Sfree(cname, fname, "temp2", temp2);
    sfree(temp2);
  }

  return iter;
}

//------------------------------------------------------------------
// Overloaded function is same as original 
// but true_res=0.
//------------------------------------------------------------------
int DiracOpWilson::MatInv(Vector *out, Vector *in, PreserveType prs_in)
{ return MatInv(out, in, 0, prs_in); }


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in and out = f_out.
//------------------------------------------------------------------
int DiracOpWilson::MatInv(Float *true_res, PreserveType prs_in)
{ return MatInv(f_out, f_in, true_res, prs_in); }


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in, out = f_out, true_res=0.
//------------------------------------------------------------------
int DiracOpWilson::MatInv(PreserveType prs_in)
{ return MatInv(f_out, f_in, 0, prs_in); }


//------------------------------------------------------------------

void DiracOpWilson::Mat(Vector *out, Vector *in) {  
  char *fname = "Mat(V*,V*)";
  VRB.Func(cname,fname);

#ifdef USE_BFM
  if(use_bfm) {
    Mat_BFM(out, in, DAG_NO);
    return;
  }
#endif

  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  if(GJP.Gparity()) temp_size *= 2;

  // points to the even part of fermion source 
  Vector *even_in = (Vector *) ( (IFloat *) in + temp_size );
  // points to the even part of fermion solution
  Vector *even_out = (Vector *) ( (IFloat *) out + temp_size );

  Dslash(out, even_in, CHKB_EVEN, DAG_NO);

  fTimesV1PlusV2((IFloat *)out, -(IFloat)kappa, (IFloat *)out,
    (IFloat *)in, temp_size);

  Dslash(even_out, in, CHKB_ODD, DAG_NO);
  
  fTimesV1PlusV2((IFloat *)even_out, -(IFloat)kappa, (IFloat *)even_out,
    (IFloat *)even_in, temp_size);
}

//------------------------------------------------------------------

void DiracOpWilson::MatDag(Vector *out, Vector *in) {
  char *fname = "MatDag(V*,V*)";
  VRB.Func(cname,fname);
  
#ifdef USE_BFM
  if(use_bfm) {
    Mat_BFM(out, in, DAG_YES);
    return;
  }
#endif

  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  if(GJP.Gparity()) temp_size *= 2;

  // points to the even part of fermion source 
  Vector *even_in = (Vector *) ( (IFloat *) in + temp_size );
  // points to the even part of fermion solution
  Vector *even_out = (Vector *) ( (IFloat *) out + temp_size );

  Dslash(out, even_in, CHKB_EVEN, DAG_YES);

  fTimesV1PlusV2((IFloat *)out, -(IFloat) kappa, (IFloat *)out,
    (IFloat *)in, temp_size);

  Dslash(even_out, in, CHKB_ODD, DAG_YES);
  
  fTimesV1PlusV2((IFloat *)even_out, -(IFloat) kappa, (IFloat *)even_out,
    (IFloat *)even_in, temp_size);
}

//------------------------------------------------------------------

void DiracOpWilson::MatHerm(Vector *out, Vector *in) {
  char *fname = "MatHerm(V*,V*)";
  VRB.Func(cname,fname);

  int temp_size = GJP.VolNodeSites() * lat.FsiteSize();
  if(GJP.Gparity()) temp_size *= 2;

  Vector *temp = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp == 0) 
    ERR.Pointer(cname, fname, "temp");
  VRB.Smalloc(cname,fname, "temp", temp, temp_size);

  Mat(temp, in);
  MultGamma(out, temp, 15, GJP.VolNodeSites());
  
  VRB.Sfree(cname, fname, "temp", temp);
  sfree(temp);
}

//------------------------------------------------------------------
/*!
  \pre This method is to be used when the instance of this object has been
  created with \a f_field_out and \a f_field_in pointers to spin-colour
  vectors defined over the whole lattice. 

  \param chi A spin-colour vector defined on odd parity lattice sites.

  \post The vector \a f_field_out is \f$ \chi \f$ on odd sites and
  \f$ D\chi \f$ on even sites.

  The vector \a f_field_in is \f$ -\kappa^2 M\chi \f$ on odd sites and
  \f$ -\kappa^2 D M\chi \f$ on even sites

  where \e M is the odd-even preconditioned fermion matrix connecting odd to
  odd parity sites and \e D is the hopping term connecting odd to
  even parity sites. Recall that \a chi is defined on odd sites only.
  The new vectors are in odd-even order.
*/
//------------------------------------------------------------------

void DiracOpWilson::CalcHmdForceVecs(Vector *chi)
{
  const char *fname = "CalcHmdForceVecs(V*)" ;
  VRB.Func(cname,fname) ;
  
#ifdef USE_BFM
  if(use_bfm) ERR.NotImplemented(cname, fname);
#endif

  if (f_out == 0)
    ERR.Pointer(cname, fname, "f_out") ;

  if (f_in == 0)
    ERR.Pointer(cname, fname, "f_in") ;

//------------------------------------------------------------------
// f_out stores (chi,rho), f_in stores (psi,sigma)
//------------------------------------------------------------------

  Vector *chi_new, *rho, *psi, *sigma ;

  size_t f_size_cb = 12 * GJP.VolNodeSites() ;
  if(GJP.Gparity()) f_size_cb *= 2;

  chi_new = f_out ;

  chi_new->CopyVec(chi, f_size_cb) ;

  psi = f_in ;

  MatPc(psi,chi) ;

  psi->VecTimesEquFloat(-kappa*kappa,f_size_cb) ;

  rho = (Vector *)((Float *)f_out + f_size_cb) ;

  Dslash(rho, chi, CHKB_ODD, DAG_NO) ;

  sigma = (Vector *)((Float *)f_in + f_size_cb) ;

  Dslash(sigma, psi, CHKB_ODD, DAG_YES) ;

  return ;
}




#ifdef USE_BFM
enum {
  ImportToBfm = 1,
  ExportFromBfm = 0
};

void DiracOpWilson::Mat_BFM(Vector *out, Vector *in, DagType dag) {  
  const char* fname = "Mat_BFM(V*,V*,DagType)";

  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  // points to the even part of fermion source 
  Vector *even_in = (Vector *) ( (IFloat *) in + temp_size );
  // points to the even part of fermion solution
  Vector *even_out = (Vector *) ( (IFloat *) out + temp_size );
  
  Vector *odd_in = in;
  Vector *odd_out = out;
  
  Fermion_t bfm_in[2]; // = {even, odd}
  Fermion_t bfm_out[2]; // = {even, odd}
  bfm_in[0] = bevo.allocFermion();
  bfm_in[1] = bevo.allocFermion();
  bfm_out[0] = bevo.allocFermion();
  bfm_out[1] = bevo.allocFermion();
  Fermion_t bfm_tmp = bevo.allocFermion();

  bevo.mass = dirac_arg->mass;

  bevo.cps_impexcbFermion((Float *)even_in , bfm_in[0], ImportToBfm, Even);
  bevo.cps_impexcbFermion((Float *)odd_in , bfm_in[1], ImportToBfm, Odd);
#pragma omp parallel
  {
    bevo.Munprec(bfm_in, bfm_out, bfm_tmp, (dag == DAG_YES ? DaggerYes : DaggerNo));
  }
  bevo.cps_impexcbFermion((Float *)even_out , bfm_out[0], ExportFromBfm, Even);
  bevo.cps_impexcbFermion((Float *)odd_out , bfm_out[1], ExportFromBfm, Odd);

  bevo.freeFermion(bfm_in[0]);
  bevo.freeFermion(bfm_in[1]);
  bevo.freeFermion(bfm_out[0]);
  bevo.freeFermion(bfm_out[1]);
  bevo.freeFermion(bfm_tmp);
}

void DiracOpWilson::MatPc_BFM(Vector *out, Vector *in, DagType dag) {  
  const char* fname = "MatPc_BFM(V*,V*,DagType)";

  Fermion_t bfm_in = bevo.allocFermion();
  Fermion_t bfm_out = bevo.allocFermion();
  Fermion_t bfm_tmp = bevo.allocFermion();

  bevo.mass = dirac_arg->mass;

  bevo.cps_impexcbFermion((Float *)in , bfm_in, ImportToBfm, Odd);
#pragma omp parallel
  {
    bevo.Mprec(bfm_in, bfm_out, bfm_tmp, (dag == DAG_YES ? DaggerYes : DaggerNo));
  }
  bevo.cps_impexcbFermion((Float *)out, bfm_out, ExportFromBfm, Odd);

  bevo.freeFermion(bfm_in);
  bevo.freeFermion(bfm_out);
  bevo.freeFermion(bfm_tmp);
}

void DiracOpWilson::MatPcDagMatPc_BFM(Vector *out, Vector *in, Float *dot_prd) {  
  const char* fname = "MatPc_BFM(V*,V*,DagType)";

  Fermion_t bfm_in = bevo.allocFermion();
  Fermion_t bfm_out = bevo.allocFermion();
  Fermion_t bfm_tmp = bevo.allocFermion();
  Fermion_t bfm_intermediate = bevo.allocFermion();

  bevo.mass = dirac_arg->mass;
  int do_nrm = (dot_prd != NULL);

  bevo.cps_impexcbFermion((Float *)in , bfm_in, ImportToBfm, Odd);

#pragma omp parallel
  {
    Float norm = bevo.Mprec(bfm_in, bfm_intermediate, bfm_tmp, DaggerNo, do_nrm);
    if(do_nrm) *dot_prd = norm;
    bevo.Mprec(bfm_intermediate, bfm_out, bfm_tmp, DaggerYes);
  }
  bevo.cps_impexcbFermion((Float *)out, bfm_out, ExportFromBfm, Odd);

  bevo.freeFermion(bfm_in);
  bevo.freeFermion(bfm_out);
  bevo.freeFermion(bfm_tmp);
  bevo.freeFermion(bfm_intermediate);
}

int DiracOpWilson::RitzEig_BFM(Vector **eigenv, Float lambda[], int valid_eig[], EigArg *eig_arg)
{
  const char* fname = "RitzEig_BFM";
  if(eig_arg->N_eig != 1) ERR.NotImplemented(cname, fname);
  if(eig_arg->RitzMatOper != MATPCDAG_MATPC &&
     eig_arg->RitzMatOper != NEG_MATPCDAG_MATPC) ERR.NotImplemented(cname, fname);

  bool compute_min = (eig_arg->RitzMatOper == MATPCDAG_MATPC);
  bevo.residual = 1e-8;
  bevo.max_iter = 1000000;

  Fermion_t bfm_eigenv = bevo.allocFermion();
  Fermion_t bfm_temp1 = bevo.allocFermion();
  Fermion_t bfm_temp2 = bevo.allocFermion();
  Fermion_t bfm_throwaway = bevo.allocFermion();
  
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  int do_nrm;
    
  Float bfm_norm_sq_eigenv, bfm_vAv, bfm_temp2_dot_eigenv, bfm_normsq_temp2;

  bevo.mass = dirac_arg->mass;

  bevo.cps_impexcbFermion((Float *)eigenv[0], bfm_eigenv, ImportToBfm, Odd);
#pragma omp parallel
  {
    lambda[0] = bevo.ritz(bfm_eigenv, compute_min);

    bfm_norm_sq_eigenv = bevo.norm(bfm_eigenv);

    do_nrm = true;
    bfm_vAv = bevo.Mprec(bfm_eigenv, bfm_temp1, bfm_throwaway, DaggerNo, do_nrm);
    bfm_normsq_temp2 = bevo.Mprec(bfm_temp1, bfm_temp2, bfm_throwaway, DaggerYes, do_nrm); 

    bfm_temp2_dot_eigenv = bevo.inner_real(bfm_temp2, bfm_eigenv);
  }
  bevo.cps_impexcbFermion((Float *)eigenv[0], bfm_eigenv, ExportFromBfm, Odd);
  
  bevo.freeFermion(bfm_eigenv);
  bevo.freeFermion(bfm_temp1);
  bevo.freeFermion(bfm_temp2);
  bevo.freeFermion(bfm_throwaway);
  
  valid_eig[0] = 1;
  
   


  Vector *temp = (Vector *)smalloc(temp_size * sizeof(Float));
  Vector *temp2 = (Vector *)smalloc(temp_size * sizeof(Float));
  MatPc(temp, eigenv[0]);  
  Float vAv = temp->NormSqGlbSum(temp_size);
  Float norm_sq_eigenv = eigenv[0]->NormSqGlbSum(temp_size);
  
  MatPcDag(temp2, temp);
  Float temp2_dot_eigenv = temp2->ReDotProductGlbSum(eigenv[0], temp_size);

  Float normsq_temp2 = temp2->NormSqGlbSum(temp_size);
  
  Fermion_t bfm_import_test = bevo.allocFermion();

  bevo.cps_impexcbFermion((Float *)eigenv[0], bfm_import_test, ImportToBfm, Odd);
  eigenv[0]->VecZero(temp_size);
  bevo.cps_impexcbFermion((Float *)eigenv[0], bfm_import_test, ExportFromBfm, Odd);

  bevo.freeFermion(bfm_import_test);
}
#endif

int DiracOpWilson::RitzEig(Vector **eigenv, Float lambda[], int valid_eig[], EigArg *eig_arg)
{
#ifdef USE_BFM
  if(use_bfm) {
    return RitzEig_BFM(eigenv, lambda, valid_eig, eig_arg);
  } else 
#endif
  {
    return DiracOpWilsonTypes::RitzEig(eigenv, lambda, valid_eig, eig_arg);
  }
}

//------------------------------------------------------------------


CPS_END_NAMESPACE
#endif
