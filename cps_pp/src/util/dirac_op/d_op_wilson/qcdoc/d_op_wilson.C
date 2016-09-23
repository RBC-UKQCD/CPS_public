/*!\file
  Wilson Dirac operator code for QCDOC

  $Id: d_op_wilson.C,v 1.10 2006/02/21 21:14:09 chulwoo Exp $
*/
//------------------------------------------------------------------
//
// d_op_wilson.C
//
// DiracOpWilsonQCDOCa is derived from the DiracOp base class. 
// DiracOpWilson is the front end for a library that contains
// all Dirac operators associated with Wilson fermions.
//
// PAB. Call the BAGEL wfm_* routines for QCDOC. 
//
// The Wilson struct is different on qcdoc, and I hide it
// in a Wilson arg and wfm class. This is wrapped in "C" linkage wrappers
// with a static instance of wfm used. The scope lock in the CPS dirac
// operators leaves this static instance safe to use. Two static instances
// will be used for the DWF case.
//
// Since (annoyingly) the wilson_init(struct Wilson *) is exposed all the way up to the
// f_wilson class, I create a "C" linkage wilson_init(struct Wilson *) wrapper in here that calls my
// wfm_init(struct WilsonArg). May wish to think about what to do for DWF here.
//
//------------------------------------------------------------------

#include <math.h>
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/wilson.h>
#include <util/wfm.h>
#include <comms/glb.h>

CPS_START_NAMESPACE


//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
DiracOpWilson::DiracOpWilson(Lattice & latt,
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
  cname = "DiracOpWilson";
  char *fname = "DiracOpWilson(L&,V*,V*,CgArg*,CnvFrmType)";
  VRB.Func(cname,fname);


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
// Destructor
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
					 Float *dot_prd)
{
  wfm_mdagm((Float *)out,
	    (Float *)gauge_field,
	    (Float *)in,
	    (Float *)dot_prd,
	    (Float)kappa);
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
			   DagType dag) 
{

    wfm_dslash((Float *)out,
	       (Float *)gauge_field,
	       (Float *)in,
	       int(cb),
	       int(dag));

}

//------------------------------------------------------------------
// MatPc(Vector *out, Vector *in) :
// MatPc is the preconditioned fermion matrix.  
// MatPc connects only odd-->odd sites.
// The in, out fields are defined on the odd checkerboard.
//------------------------------------------------------------------
void DiracOpWilson::MatPc(Vector *out, Vector *in) 
{  
    wfm_m((Float *)out,
	  (Float *)gauge_field,
	  (Float *)in,
	  (Float)kappa);
}

//------------------------------------------------------------------
// MatPcDag(Vector *out, Vector *in) :
// MatPcDag is the dagger of the preconditioned fermion matrix. 
// MatPcDag connects only odd-->odd sites.
// The in, out fields are defined on the odd checkerboard.
//------------------------------------------------------------------
void DiracOpWilson::MatPcDag(Vector *out, Vector *in) 
{
//    printf("kappa=%e\n",kappa);
    wfm_mdag((Float *)out,
	     (Float *)gauge_field,
	     (Float *)in,
             (Float)kappa);
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
  Vector *temp2;

  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  // check out if converted
  //for (int ii = 0; ii < 2 * temp_size; ii++) {
  //  VRB.Result(cname, fname, "in[%d] = %e\n", ii, 
  //  *((Float *)in + ii));
  //  VRB.Result(cname, fname, "out[%d] = %e\n", ii, 
  //  *((Float *)out + ii));
  //}

  Vector *temp = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp == 0) ERR.Pointer(cname, fname, "temp");
  VRB.Smalloc(cname,fname, "temp", temp, temp_size * sizeof(Float));

  if(prs_in == PRESERVE_YES){
    temp2 = (Vector *) smalloc(2*temp_size * sizeof(Float));
    if (temp2 == 0) ERR.Pointer(cname, fname, "temp2");
    VRB.Smalloc(cname,fname, "temp2", temp2, temp_size * sizeof(Float));
  }
  // save source
  if(prs_in == PRESERVE_YES){
    moveMem((Float *)temp2, (Float *)in, 2*temp_size*sizeof(Float));
  }

#if 0
{
  printf("in(before)=\n");
  IFloat *temp_p = (IFloat *)in;
  for(int ii = 0; ii< GJP.VolNodeSites();ii++){
    for(int jj = 0; jj< lat.FsiteSize();jj++){
      if (fabs(*temp_p)>1e-7){
        printf("i=%d j=%d\n",ii,jj);
        printf("%e\n",*(temp_p));
      }
      temp_p++;
    }
  }
}
#endif

  // points to the even part of fermion source 
  Vector *even_in = (Vector *) ( (Float *) in + temp_size );

  // points to the even part of fermion solution
  Vector *even_out = (Vector *) ( (Float *) out + temp_size );

  Dslash(temp, even_in, CHKB_EVEN, DAG_NO);

  fTimesV1PlusV2((Float *)temp, (Float) kappa, (Float *)temp,
    (Float *)in, temp_size);

#if 0
{
  printf("temp(before)=\n");
  IFloat *temp_p = (IFloat *)temp;
  for(int ii = 0; ii< GJP.VolNodeSites();ii++){
    for(int jj = 0; jj< lat.FsiteSize();jj++){
      if (fabs(*temp_p)>1e-7){
        printf("i=%d j=%d\n",ii,jj);
        printf("%e\n",*(temp_p));
      }
      temp_p++;
    }
  }
}
#endif

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

  Dslash(temp, out, CHKB_ODD, DAG_NO);

  fTimesV1PlusV2((Float *)even_out, (Float) kappa, (Float *)temp,
    (Float *) even_in, temp_size);

  VRB.Sfree(cname, fname, "temp", temp);
  sfree(temp);


  // restore source
  if(prs_in == PRESERVE_YES){
    moveMem((Float *)in, (Float *)temp2, 2*temp_size*sizeof(Float));
  }

#if 0
{
  printf("in(after)=\n");
  IFloat *temp_p = (IFloat *)in;
  for(int ii = 0; ii< GJP.VolNodeSites();ii++){
    for(int jj = 0; jj< lat.FsiteSize();jj++){
      if (fabs(*temp_p)>1e-7){
        printf("i=%d j=%d\n",ii,jj);
        printf("%e\n",*(temp_p));
      }
      temp_p++;
    }
  }
}
#endif

#if 0
{
  printf("temp2(after)=\n");
  IFloat *temp_p = (IFloat *)temp2;
  for(int ii = 0; ii< GJP.VolNodeSites();ii++){
    for(int jj = 0; jj< lat.FsiteSize();jj++){
      if (fabs(*temp_p)>1e-7){
        printf("i=%d j=%d\n",ii,jj);
        printf("%e\n",*(temp_p));
      }
      temp_p++;
    }
  }
}
#endif

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
{
  return MatInv(out, in, 0, prs_in); 
}


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in and out = f_out.
//------------------------------------------------------------------
int DiracOpWilson::MatInv(Float *true_res, PreserveType prs_in)
{
  return MatInv(f_out, f_in, true_res, prs_in); 
}


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in, out = f_out, true_res=0.
//------------------------------------------------------------------
int DiracOpWilson::MatInv(PreserveType prs_in)
{
  return MatInv(f_out, f_in, 0, prs_in); 
}


//------------------------------------------------------------------
// Mat(Vector *out, Vector *in) :
// Mat is the unpreconditioned fermion matrix.  
// Mat works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpWilson::Mat(Vector *out, Vector *in) 
{  
  char *fname = "Mat(V*,V*)";
  VRB.Func(cname,fname);

  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  // points to the even part of fermion source 
  Vector *even_in = (Vector *) ( (Float *) in + temp_size );
  // points to the even part of fermion solution
  Vector *even_out = (Vector *) ( (Float *) out + temp_size );

  Dslash(out, even_in, CHKB_EVEN, DAG_NO);

  fTimesV1PlusV2((Float *)out, -(Float)kappa, (Float *)out,
    (Float *)in, temp_size);

  Dslash(even_out, in, CHKB_ODD, DAG_NO);
  
  fTimesV1PlusV2((Float *)even_out, -(Float)kappa, (Float *)even_out,
    (Float *)even_in, temp_size);
}

//------------------------------------------------------------------
// MatDag(Vector *out, Vector *in) :
// MatDag is the unpreconditioned fermion matrix.  
// MatDag works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpWilson::MatDag(Vector *out, Vector *in) 
{
  char *fname = "MatDag(V*,V*)";
  VRB.Func(cname,fname);

  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  // points to the even part of fermion source 
  Vector *even_in = (Vector *) ( (Float *) in + temp_size );
  // points to the even part of fermion solution
  Vector *even_out = (Vector *) ( (Float *) out + temp_size );

  Dslash(out, even_in, CHKB_EVEN, DAG_YES);

  fTimesV1PlusV2((Float *)out, -(Float) kappa, (Float *)out,
    (Float *)in, temp_size);

  Dslash(even_out, in, CHKB_ODD, DAG_YES);
  
  fTimesV1PlusV2((Float *)even_out, -(Float) kappa, (Float *)even_out,
    (Float *)even_in, temp_size);
}

//------------------------------------------------------------------
// MatHerm(Vector *out, Vector *in) :
// MatHerm is gamma5*Mat.
// MatHerm works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpWilson::MatHerm(Vector *out, Vector *in) 
{
  char *fname = "MatHerm(V*,V*)";
  VRB.Func(cname,fname);

  int temp_size = GJP.VolNodeSites() * lat.FsiteSize();
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
// GRF
// chi is the solution to MatPcInv.  The user passes two full size
// CANONICAL fermion vectors with conversion enabled to the
// constructor.  Using chi, the function fills these vectors;
// the result may be used to compute the HMD fermion force.
//------------------------------------------------------------------

void DiracOpWilson::CalcHmdForceVecs(Vector *chi)
{
  char *fname = "CalcHmdForceVecs(V*)" ;
  VRB.Func(cname,fname) ;

  if (f_out == 0)
    ERR.Pointer(cname, fname, "f_out") ;

  if (f_in == 0)
    ERR.Pointer(cname, fname, "f_in") ;

//------------------------------------------------------------------
// f_out stores (chi,rho), f_in stores (psi,sigma)
//------------------------------------------------------------------

  Vector *chi_new, *rho, *psi, *sigma ;

  int f_size_cb = 12 * GJP.VolNodeSites() ;

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

CPS_END_NAMESPACE
