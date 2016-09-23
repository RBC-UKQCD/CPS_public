#include <config.h>
#include <stdio.h>
#include <alg/no_arg.h>
#include <alg/common_arg.h>
#include <alg/alg_plaq.h>
#include <util/data_types.h>

CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOpMobius class methods.

  $Id: d_op_mobius.C,v 1.6 2013-06-07 19:26:34 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-06-07 19:26:34 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/d_op_mobius.C,v 1.6 2013-06-07 19:26:34 chulwoo Exp $
//  $Id: d_op_mobius.C,v 1.6 2013-06-07 19:26:34 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: d_op_mobius.C,v $
//  $Revision: 1.6 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/d_op_mobius.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// d_op_mobius.C
//
// DiracOpMobius is derived from the DiracOp base class. 
// DiracOpMobius is the front end for a library that contains
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
#include <util/zmobius.h>
//#include <mem/p2v.h>
#include <comms/glb.h>

//#define USE_BLAS
#ifdef USE_BLAS
#include "noarch/blas-subs.h"
#endif

CPS_START_NAMESPACE

void vecTimesEquComplex(Complex *a, Complex b, int len)
{
//#pragma omp parallel for
    for(int i = 0; i < len/2; ++i) {
    	a[i] *= b;
    }
}



void zTimesV1PlusV2(Complex *a, Complex b, const Complex *c,
		    const Complex *d, int len)
{
//#pragma omp parallel for
    for(int i = 0; i < len/2; ++i) {
    	a[i] = b * c[i] + d[i];
    }
}


void vecEqualsVecTimesEquComplex(Complex *a, Complex *b, Complex c, int len)
{
//  VRB.Result("", "vecEqualsVecTimesEquComplex()", "(%p %p %g %g %d)\n", a, b, c.real(),c.imag(),len);
  for (int i=0; i<len/2; i++) {
    *a++ = c * *b++;
  }

}

//  temp[s]  <-  kappa_b[s] *temp2[s]  +  in[s]  s=0..ls-1,
// ls_stride is the number of Floats stored in one s-slice
// s_node_coor is the coordinate of the current node in s-direction 
void zmobius_zvectTimesV1PlusV2 (Vector* temp, Complex* kappa_b,  Vector* temp2,
				 Vector* in, int local_ls, int ls_stride, int s_node_coor )
{
  for(int s=0; s<local_ls;++s){
    int glb_s = s + local_ls*s_node_coor;
    int idx = s*ls_stride/2;// "/2" is for complex
    zTimesV1PlusV2((Complex*) temp+idx, kappa_b[glb_s], (Complex*) temp2+idx,
		   (Complex*)in+idx, ls_stride);
  }
}

//  temp[s]  <-  kappa_b[s] *temp[s]  s=0..ls-1,
// ls_stride is the number of Floats stored in one s-slice
// s_node_coor is the coordinate of the current node in s-direction 
  void zmobius_zvectTimesEquComplex(Vector* temp, Complex* kappa_b, 
				  int local_ls, int ls_stride, int s_node_coor )
{
  for(int s=0; s<local_ls;++s){
    int glb_s = s + local_ls*s_node_coor;
    int idx = s*ls_stride/2;// "/2" is for complex
    vecTimesEquComplex((Complex*)temp+idx,  kappa_b[glb_s], ls_stride);
  }
}



//B matrix in Andrewes' note, MIT
//  ( b + C dslash5 )  or its dagger
// FIXME:  this is slow,  for experimental purpose only
void zmobius_B_MIT( Vector* out, Float mass,
		    int dag, Zmobus *mobius_lib_arg,
		    Complex* b, Complex*c )
{
  const char* cname="";
  const char* fname="zmobius_B_MIT(...)";
  
  const int ls = GJP.SnodeSites();
  const int vol_4d_cb = mobius_lib_arg->vol_4d / 2;
  const int f_size = 24 * mobius_lib_arg->ls * vol_4d_cb;
  const int ls_stride = 24 * vol_4d_cb;
  
  Vector *temp = (Vector *) smalloc( cname,fname, "temp", f_size * sizeof(Float));

  Complex* fact = new Complex [ls]; 

  for(int i=0;i<ls;++i)
    fact[i] = GJP.ZMobius_c()[i]/ GJP.ZMobius_b()[i];
  
  if(dag==0){
    moveFloat((IFloat*)temp, (IFloat*) out, f_size);
    zmobius_kappa_dslash_5_plus_cmplx(out,temp,
				      mass,dag, 
				      mobius_lib_arg,
				      fact);
    for(int i=0;i<ls;++i){
      Complex b = GJP.ZMobius_b()[i];  //b=1.0;
      int idx=i*ls_stride/2; // 2 is for complex
      vecEqualsVecTimesEquComplex((Complex*)out +idx, (Complex*)out+idx,
				  b, ls_stride);
    }
    
  } else { // dag==1

    moveFloat((IFloat*)temp, (IFloat*) out, f_size);
    for(int i=0;i<ls;++i){
      Complex b = GJP.ZMobius_b()[i]; //b=1.0;
      int idx=i*ls_stride/2; // 2 is for complex
      vecEqualsVecTimesEquComplex((Complex*)temp +idx, (Complex*)temp+idx,
				  conj(b), ls_stride);
    }
    moveFloat((IFloat*)out, (IFloat*) temp, f_size);
    zmobius_kappa_dslash_5_plus_cmplx(out,temp,
				      mass,dag, 
				      mobius_lib_arg,
				      fact);
  }


  

  sfree(temp);
  VRB.Sfree(cname, fname, "temp", temp);
  delete [] fact;
}
  
//B matrix in Andrewes' note, MIT, slow
//  ( b + C dslash5 )inverse  or its dagger
void zmobius_Binv_MIT( Vector* out, Float mass,
		    int dag, Zmobus *mobius_lib_arg,
		    Complex* b, Complex*c )
{
  const char* cname="";
  const char* fname="zmobius_Binv_MIT(...)";

  const int ls = GJP.SnodeSites();
  const int vol_4d_cb = mobius_lib_arg->vol_4d / 2;
  const int f_size = 24 * mobius_lib_arg->ls * vol_4d_cb;
  const int ls_stride = 24 * vol_4d_cb;
  


  Complex* fact = new Complex [ls]; 
  for(int i=0;i<ls;++i)
    fact[i] = GJP.ZMobius_c()[i]/ GJP.ZMobius_b()[i];

  if(dag==0) {
    for(int i=0;i<ls;++i){
      Complex inv_b = 1.0/ GJP.ZMobius_b()[i]; //inv_b=1.0;
      int idx=i*ls_stride/2; // 2 is for complex
      vecEqualsVecTimesEquComplex((Complex*)out +idx, (Complex*)out+idx,
				  inv_b, ls_stride);
    }
    zmobius_m5inv(out, mass, dag, mobius_lib_arg, fact);
    
  } else {//dag==1

    for(int i=0;i<ls;++i)
      fact[i] = GJP.ZMobius_c()[i]/ GJP.ZMobius_b()[i];
    zmobius_m5inv(out, mass, dag, mobius_lib_arg, fact);
    
    for(int i=0;i<ls;++i){
      Complex inv_b = 1.0/ conj(GJP.ZMobius_b()[i]);  //inv_b=1.0;
      int idx=i*ls_stride/2; // 2 is for complex
      vecEqualsVecTimesEquComplex((Complex*)out +idx, (Complex*)out+idx,
				  inv_b, ls_stride);
    }

  }


  delete [] fact;
}




//----------------------------------------------------------------
// Initialize kappa and ls. This has already been done by the Fmobius
// call to mobius_init but is done here again in case the user
// has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
// GJP.SnodeSites() while in the scope of the Fmobius object.
//----------------------------------------------------------------

void reset_mobius_arg(void* mobius_lib_arg){
  Zmobus *mobius_arg = (Zmobus *) mobius_lib_arg;
  mobius_arg->ls = GJP.SnodeSites();

  mobius_arg->pc_type = GJP.ZMobius_PC_Type();

  if(mobius_arg->zmobius_kappa_b){
    delete [] mobius_arg->zmobius_kappa_b;
  }
  mobius_arg->zmobius_kappa_b = (Complex*) new  Complex [mobius_arg->ls];
  if(mobius_arg->zmobius_kappa_c){
    delete [] mobius_arg->zmobius_kappa_c;
  }
  mobius_arg->zmobius_kappa_c = (Complex*) new  Complex [mobius_arg->ls];
  if(mobius_arg->zmobius_kappa_ratio){
    delete [] mobius_arg->zmobius_kappa_ratio;
  }
  mobius_arg->zmobius_kappa_ratio = (Complex*) new  Complex [mobius_arg->ls];


  //  printf("TIZB Height %e %e\n", GJP.DwfHeight(), GJP.DwfA5Inv());
  for(int i=0;i<mobius_arg->ls;++i){
    mobius_arg->zmobius_kappa_b[i] = 
      1.0 / ( 2 * (GJP.ZMobius_b()[i] *(4 - GJP.DwfHeight()) + GJP.DwfA5Inv()) );
    mobius_arg->zmobius_kappa_c[i] = 
      1.0 / ( 2 * (GJP.ZMobius_c()[i] *(4 - GJP.DwfHeight()) - GJP.DwfA5Inv()) );
    mobius_arg->zmobius_kappa_ratio[i] =
      mobius_arg->zmobius_kappa_b[i]  /      mobius_arg->zmobius_kappa_c[i];

#if 0
    printf("TIZB: kappas  %d %e %e %e %e  %e %e  %e %e  %e %e\n",
	   i,
	   GJP.ZMobius_b()[i].real(),	   GJP.ZMobius_b()[i].imag(),
	   GJP.ZMobius_c()[i].real(),	   GJP.ZMobius_c()[i].imag(),
	   mobius_arg->zmobius_kappa_b[i].real(),
	   mobius_arg->zmobius_kappa_b[i].imag(),
	   mobius_arg->zmobius_kappa_c[i].real(),
	   mobius_arg->zmobius_kappa_c[i].imag(),
	   mobius_arg->zmobius_kappa_ratio[i].real(),
	   mobius_arg->zmobius_kappa_ratio[i].imag());
   
#endif
  }
}

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
DiracOpZMobius::DiracOpZMobius(Lattice & latt,
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
  cname = "DiracOpZMobius";
  char *fname = "DiracOpZMobius(L&,V*,V*,CgArg*,CnvFrmType)";
  VRB.Func(cname,fname);


  //----------------------------------------------------------------
  // Do the necessary conversions
  //----------------------------------------------------------------
#undef PROFILE
#ifdef PROFILE
  Float time = -dclock();
#endif
  if(cnv_frm == CNV_FRM_YES)
    lat.Convert( DWF_4D_EOPREC_EE, f_out, f_in);
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
  // Initialize the pointer to the initialized Mobius structure
  // (the structure has been initialized by the Lattice::Fmobius
  // constructor.
  //----------------------------------------------------------------
  mobius_lib_arg = lat.FdiracOpInitPtr();
  
  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // call to mobius_init but is done here again in case the user
  // has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the Fmobius object.
  //----------------------------------------------------------------
  reset_mobius_arg( mobius_lib_arg);


}


//------------------------------------------------------------------
/*!
  If the storage order of any fields was changed by the constructor
  then they are changed to the canonical order by the destructor.
*/
//------------------------------------------------------------------
DiracOpZMobius::~DiracOpZMobius() {
  char *fname = "~DiracOpZMobius()";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Do the necessary conversions
  //----------------------------------------------------------------
#undef PROFILE
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
void DiracOpZMobius::DiracArg(CgArg *arg){
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
void DiracOpZMobius::MatPcDagMatPc(Vector *out, 
				  Vector *in, 
				  Float *dot_prd){

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpZMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpZMobius object.
  //----------------------------------------------------------------
  reset_mobius_arg( mobius_lib_arg);

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  zmobius_mdagm(out, 
	       gauge_field, 
	       in, 
	       dot_prd,
	       mass,
	       (Zmobus *) mobius_lib_arg);

}

//------------------------------------------------------------------
// MatPcDagMatPc :
// MatPcDagMatPc is the fermion matrix that appears in the HMC 
// evolution. It is a Hermitian matrix.
// The in, out fields are defined on the checkerboard lattice.
// If dot_prd is not 0 then the dot product (on node)
// <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.
//
// When dirac_arg->eigen_shift is non-zero, it shift the spectrum of matrix:
//    MatPcDagMatPc = H^2  ->  (H-shift)(H-shift)
// where H = Gamma_5 MatPc
//
// For other fermions, one could also implement similar shifts.
// For wilson, H = gamma_5 MatPC .
//------------------------------------------------------------------
void DiracOpZMobius::MatPcDagMatPcShift(Vector *out, 
				       Vector *in, 
				       Float *dot_prd){

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpZMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpZMobius object.
  //----------------------------------------------------------------

  const Float shift = dirac_arg->eigen_shift;
  reset_mobius_arg( mobius_lib_arg);

  

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  
  // we still check if shift is really needed
  if( shift == 0.0 )  {
    zmobius_mdagm(out, 
		 gauge_field, 
		 in, 
		 dot_prd,
		 mass,
		 (Zmobus *) mobius_lib_arg);
  } else {
    
    //mobius_arg->eigen_shift = dirac_arg->eigen_shift;
    
    zmobius_mdagm_shift(out, 
		       gauge_field, 
		       in, 
		       dot_prd,
		       mass,
		       (Zmobus *) mobius_lib_arg,
		       shift);
  }
  
}


//------------------------------------------------------------------
// Dslash(Vector *out, Vector *in, ChkbType cb, DagType dag) :
// Dslash is the derivative part of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
// cb is not used.
//------------------------------------------------------------------
void DiracOpZMobius::Dslash(Vector *out, 
			   Vector *in, 
			   ChkbType cb, 
			   DagType dag) {
  
  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpZMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpZMobius object.
  //----------------------------------------------------------------
  reset_mobius_arg( mobius_lib_arg);

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  zmobius_unprec(out, 
		gauge_field, 
		in, 
		mass,
//		cb,
		dag,
		(Zmobus *) mobius_lib_arg);
}

//------------------------------------------------------------------
/*!
  The preconditioned matrix connects sites of even parity.
  \param out The resulting vector.
  \param in The vector to be multiplied.
*/
//------------------------------------------------------------------
void DiracOpZMobius::MatPc(Vector *out, Vector *in) {  

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpZMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpZMobius object.
  //----------------------------------------------------------------
  reset_mobius_arg( mobius_lib_arg);

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  zmobius_m(out, 
	   gauge_field, 
	   in, 
	   mass,
	   (Zmobus *) mobius_lib_arg);
}

//------------------------------------------------------------------
/*!
  The preconditioned matrix connects sites of even parity.
  \param out The resulting vector.
  \param in The vector to be multiplied.
*/
//------------------------------------------------------------------
void DiracOpZMobius::MatPcDag(Vector *out, Vector *in) {

  char *fname = "MatPcDag(*V,*V)";
  VRB.Func(fname,cname);
  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpZMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpZMobius object.
  //----------------------------------------------------------------
  reset_mobius_arg( mobius_lib_arg);

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  IFloat *tmp = (IFloat *)in;
  zmobius_mdag(out, 
	      gauge_field, 
	      in, 
	      mass,
	      (Zmobus *) mobius_lib_arg);
  tmp = (IFloat *)out;
}



 //----------------------------------
void test_m5inv(void* mobius_lib_arg, Lattice &lat, Float mass, int dag )
{
  
  const char* cname=""; const char* fname="";
  Zmobus *mobius_arg = (Zmobus *) mobius_lib_arg;

  //Float minus_kappa_b = -mobius_arg->mobius_kappa_b;
  //Float kappa_b = - minus_kappa_b;
  Float norm;

  //printf("KAPPA_B %g\n",kappa_b); exit(0);

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  Vector *temp2;
  Vector *temp3;
  Vector *save_in;

  const int local_ls = GJP.SnodeSites();
  const int s_nodes = GJP.Snodes();
  const int global_ls = local_ls * s_nodes;
  const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = ((Zmobus*)mobius_lib_arg)->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;

  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  Vector *temp = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp == 0) ERR.Pointer(cname, fname, "temp");
  VRB.Smalloc(cname,fname, "temp", temp, temp_size * sizeof(Float));

  temp2 = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp2 == 0) ERR.Pointer(cname, fname, "temp2");
  VRB.Smalloc(cname,fname, "temp2", temp2, temp_size * sizeof(Float));


  
  if(!UniqueID())printf("TIZB mass=%e\n",mass);
  for(int TID=0; TID<=12; TID+=12){
    Complex mat_orig[100][100];
    Complex mat_inv[100][100];
    for(int iturn=0;iturn<2;++iturn){
      
      for (int s1=0;s1<local_ls;++s1){
	//set src
	for(int s=0;s<local_ls;++s){
	  for(int i=0;i< vol_4d_cb;++i){
	    for(int id=0;id<24;++id){
	      int idx=id+24*(i+vol_4d_cb*s);
	      Float *fp=(Float*)temp;
	      fp[idx]=0.0;
	      if(s==s1 && i==0 &&  id==TID)
		fp[idx]=1.0;
	    }}}
	
	if(iturn==0)
	  {
	    zmobius_m5inv(temp2, temp, mass, dag, mobius_arg,
			  mobius_arg->zmobius_kappa_ratio);  	  
	  }
	else { // this is the M5
	  moveFloat( (IFloat*)temp2, (IFloat*)temp, temp_size );
	  Complex* kappa_ratio= mobius_arg->zmobius_kappa_ratio;
	  zmobius_kappa_dslash_5_plus_cmplx(temp2,
					    temp,
					    mass,
					    dag,
					    mobius_arg,
					    kappa_ratio);
	}
	
	
        {
	  Float *fp=(Float*)temp2;
	  int id=TID;int i=0;
	  for(int s2=0;s2<local_ls;++s2){
	    int idx=id+24*(i+vol_4d_cb*s2);
	    if(iturn==0)
	      mat_inv[s2][s1]=Complex(fp[idx], fp[idx+1]);
	    else
	      mat_orig[s2][s1]=Complex(fp[idx], fp[idx+1]);
	  }
	}
	
      }
    }
    
    Complex mat_prod[100][100];
    for(int s1=0;s1<local_ls;++s1) for(int s2=0;s2<local_ls;++s2)
			      mat_prod[s1][s2]=0.0;
    
    for(int s1=0;s1<local_ls;++s1)
      for(int s2=0;s2<local_ls;++s2)
	for(int s=0;s<local_ls;++s)
	  mat_prod[s1][s2] += mat_orig[s1][s]*mat_inv[s][s2];
    
    printf("TID=\%d\n",TID);
#if 0
    printf("===  mat_orig ===\n");
    for(int s1=0;s1<local_ls;++s1){
      for(int s2=0;s2<local_ls;++s2) printf("%+.2e ", mat_orig[s1][s2].real());
      printf("\n");
    }
    printf("===  mat_inv ===\n");
    for(int s1=0;s1<local_ls;++s1){
      for(int s2=0;s2<local_ls;++s2) printf("%+.2e ", mat_inv[s1][s2].real());
      printf("\n");
    }
#endif
    printf("===  mat_prod ===\n");
    for(int s1=0;s1<local_ls;++s1){
      for(int s2=0;s2<local_ls;++s2) 
	printf("(%+.2e,%+.2e) ", mat_prod[s1][s2].real(),mat_prod[s1][s2].imag());
      printf("\n");
    }
  }
  
}

//----------------------------------
void test_m5inv_norm(void* mobius_lib_arg, Lattice &lat, Float mass, int dag )
{
  
  const char* cname=""; const char* fname="test_m5inv_norm";
  Zmobus *mobius_arg = (Zmobus *) mobius_lib_arg;

  //Float minus_kappa_b = -mobius_arg->mobius_kappa_b;
  //Float kappa_b = - minus_kappa_b;
  Float norm;

  //printf("KAPPA_B %g\n",kappa_b); exit(0);

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  Vector *temp2;
  Vector *temp3;
  Vector *save_in;

  const int local_ls = GJP.SnodeSites();
  const int s_nodes = GJP.Snodes();
  const int global_ls = local_ls * s_nodes;
  const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = ((Zmobus*)mobius_lib_arg)->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;

  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  Vector *temp = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp == 0) ERR.Pointer(cname, fname, "temp");
  VRB.Smalloc(cname,fname, "temp", temp, temp_size * sizeof(Float));

  temp2 = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp2 == 0) ERR.Pointer(cname, fname, "temp2");
  VRB.Smalloc(cname,fname, "temp2", temp2, temp_size * sizeof(Float));

  temp3 = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp2 == 0) ERR.Pointer(cname, fname, "temp3");
  VRB.Smalloc(cname,fname, "temp3", temp3, temp_size * sizeof(Float));

  srand48(19700215);
  
  for(int s=0;s<local_ls;++s){
    for(int i=0;i< vol_4d_cb;++i){
      for(int id=0;id<24;++id){
	int idx=id+24*(i+vol_4d_cb*s);
	Float *fp=(Float*)temp;
	fp[idx]=drand48();
      }
    }
  }
  moveFloat((IFloat*)temp3,(IFloat*)temp,temp_size);
  
  {
    Float norm = temp-> NormSqGlbSum(temp_size);
    if(!UniqueID()) printf("%s in Norm %.16e\n", fname, norm);
  }

  
  zmobius_m5inv(temp2, temp, mass, dag, mobius_arg,
		mobius_arg->zmobius_kappa_ratio);  	  
  {
    moveFloat( (IFloat*)temp, (IFloat*)temp2, temp_size );
    Complex* kappa_ratio= mobius_arg->zmobius_kappa_ratio;
    zmobius_kappa_dslash_5_plus_cmplx(temp,temp2,mass,dag,
				      mobius_arg, kappa_ratio);
  }

  {
    Float norm = temp-> NormSqGlbSum(temp_size);
    if(!UniqueID()) printf("%s check Norm %.16e\n", fname, norm);
  }

  IFloat *fp=(IFloat*)temp;
  IFloat *fp3=(IFloat*)temp3;
  Float max=-100;
  Float max_rel=-100;
  for(int s=0;s<local_ls;++s){
    for(int i=0;i< vol_4d_cb;++i){
      for(int id=0;id<24;++id){
	int idx=id+24*(i+vol_4d_cb*s);
	Float diff = fabs( fp[idx]-fp3[idx] );
	if(diff > max) max=diff;
	Float diff_rel = 2.0*fabs( fp[idx]-fp3[idx] )/fabs( fp[idx]+fp3[idx] );
	if(diff_rel > max_rel) max_rel=diff_rel;
      }
    }
  }
  glb_max(&max);
  glb_max(&max_rel);
  printf("max diff %e max rel diff %e\n", max, max_rel);
  
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
int DiracOpZMobius::MatInv(Vector *out, 
			  Vector *in, 
			  Float *true_res,
			  PreserveType prs_in) {
  char *fname = "MatInv(V*,V*,F*)";
  VRB.Func(cname,fname);
  
  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpZMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpZMobius object.
  //----------------------------------------------------------------
  //printf("KAPPA_B %g\n",((Zmobus*)mobius_lib_arg)->mobius_kappa_b); exit(0);

  reset_mobius_arg( mobius_lib_arg);

  Zmobus *mobius_arg = (Zmobus *) mobius_lib_arg;

  //Float minus_kappa_b = -mobius_arg->mobius_kappa_b;
  //Float kappa_b = - minus_kappa_b;
  Float norm;

  //printf("KAPPA_B %g\n",kappa_b); exit(0);

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  Vector *temp2;
  Vector *temp3;
  Vector *save_in;

  const int local_ls = GJP.SnodeSites();
  const int s_nodes = GJP.Snodes();
  const int global_ls = local_ls * s_nodes;
  const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = ((Zmobus*)mobius_lib_arg)->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;

  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  Vector *temp = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp == 0) ERR.Pointer(cname, fname, "temp");
  VRB.Smalloc(cname,fname, "temp", temp, temp_size * sizeof(Float));

  temp2 = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp2 == 0) ERR.Pointer(cname, fname, "temp2");
  VRB.Smalloc(cname,fname, "temp2", temp2, temp_size * sizeof(Float));

  // points to the even part of fermion source 
  Vector *odd_in = (Vector *) ( (IFloat *) in + temp_size );

  // points to the even part of fermion solution
  Vector *odd_out = (Vector *) ( (IFloat *) out + temp_size );

  // prepare source
  // mult by Dminus to compare with Hantao
  // do outside in f_mobius class instead
#if 0
  Dminus(temp3,in);
  //moveFloat((IFloat *)in, (IFloat *)temp3, 2*temp_size);
  //VRB.Sfree(cname, fname, "temp3", temp3);
  //sfree(temp3);
#endif
  //  DEBTIZB("insrc", (Vector*) in, 2*temp_size);
  
    
  // save source
  if(prs_in == PRESERVE_YES){
    temp3 = (Vector *) smalloc(2*temp_size * sizeof(Float));
    if (temp3 == 0) ERR.Pointer(cname, fname, "temp2");
    VRB.Smalloc(cname,fname, "temp3", temp3, 2*temp_size * sizeof(Float));
    moveMem((IFloat *)temp3, (IFloat *)in, 2*temp_size * sizeof(IFloat));
  }


#if 0
  test_m5inv(mobius_lib_arg, lat,  mass,0);
  test_m5inv(mobius_lib_arg, lat,  mass,1);
  test_m5inv_norm(mobius_lib_arg, lat,  mass,0);
  test_m5inv_norm(mobius_lib_arg, lat,  mass,1);
  exit(1);
#endif
  
#if 0
  //----------------------------------
  // check for m5inv
  zmobius_m5inv(temp, odd_in, mass, DAG_NO, mobius_arg,
		mobius_arg->zmobius_kappa_ratio);  

  norm = odd_in->NormSqGlbSum(temp_size);
  for(int i=0;i<24;++i) printf("%e ", *((IFloat*)odd_in+i)); printf("\n");
  if(!UniqueID()) printf("TIZB M5 Norm odd_in %.14e\n",norm);
  norm = temp->NormSqGlbSum(temp_size);
  if(!UniqueID()) printf("M5 Norm temp  %.14e\n",norm);
  for(int i=0;i<24;++i) printf("%e ", *((IFloat*)temp+i)); printf("\n");
  {
    moveFloat( (IFloat*)temp2, (IFloat*)temp, temp_size );
    Complex* kappa_ratio= mobius_arg->zmobius_kappa_ratio;
    zmobius_kappa_dslash_5_plus_cmplx(temp2,
				    temp,
				    mass,
				    DAG_NO,
				    mobius_arg,
				      kappa_ratio);
    norm = temp2->NormSqGlbSum(temp_size);
    if(!UniqueID()) printf("M5 Norm temp2  %.14e\n",norm);
  }
  exit(1);
  //check end
  //----------------------------------
#endif

  //-------------------
  // Preconditioning source and guess vector
  //-------------------

  // for source
  switch (mobius_arg->pc_type){
  case ZMOB_PC_ORIG:
  case ZMOB_PC_SYM2:
  case ZMOB_PC_SYM3:
    zmobius_m5inv(temp, odd_in, mass, DAG_NO, mobius_arg,
		  mobius_arg->zmobius_kappa_ratio);  

    zmobius_dslash_4(temp2, gauge_field, temp, CHKB_ODD, DAG_NO,
		     mobius_arg, mass);
    zmobius_zvectTimesV1PlusV2 (temp, mobius_arg->zmobius_kappa_b,  temp2, in,
				local_ls, ls_stride, s_node_coor );
    break;
  default:
    ERR.NotImplemented(cname,fname);
    break;
  }
  // for guess vector
  switch(mobius_arg->pc_type){
  case ZMOB_PC_ORIG:
    // need nothing
    break;
  case ZMOB_PC_SYM2:
  case ZMOB_PC_SYM3:
    // Apply M5 to out for sym2 preconditioning
    moveFloat((IFloat *)temp2, (IFloat *)out, temp_size );
    zmobius_kappa_dslash_5_plus_cmplx(out, temp2, mass, DAG_NO, mobius_arg,
				      mobius_arg->zmobius_kappa_ratio);
    break;
    default:
      ERR.NotImplemented(cname,fname);
      break;
  }
    
  int iter;
  switch (dirac_arg->Inverter) {
  case CG:
    MatPcDag(in, temp);
    iter = InvCg(out,in,true_res);
    break;
  case BICGSTAB:
    iter = BiCGstab(out,temp,0.0,dirac_arg->bicgstab_n,true_res);
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
    break;
  }


#if 1
  // check solution
  norm = out->NormSqGlbSum(temp_size);
  if(!UniqueID()) printf("Norm out %.14e\n",norm);
  norm = in->NormSqGlbSum(temp_size);
  if(!UniqueID()) printf("Norm in %.14e\n",norm);
  MatPcDagMatPc(temp,out);  
  norm = temp->NormSqGlbSum(temp_size);
  if(!UniqueID()) printf("Norm MatPcDagMatPc*out %.14e\n",norm);
  //exit(0);
#endif

  // restore source
  if(prs_in == PRESERVE_YES){
    moveMem((IFloat *)in, (IFloat *)temp3,
            2*temp_size * sizeof(IFloat) / sizeof(char));
  }

  
  //----------------------------------------------
  // post preconditioning aka  unpreconditioning 
  //----------------------------------------------
  //  odd_out = M5inv kappa_b M4 out + M5inv odd_in

  // For ZMOBIUS_PC_SYM2
  //  out_final  = M5inv out
  //  out_out  = M5inv kappa_b M4 M5inv out  + M5inv odd_in
  //           = M5inv kappa_b M4 out_final  + M5inv odd_in

  
  
  switch(mobius_arg->pc_type){
  case   ZMOB_PC_SYM3:
  case   ZMOB_PC_SYM2: {
    zmobius_m5inv(temp, out, mass, DAG_NO, mobius_arg,
		  mobius_arg->zmobius_kappa_ratio);
    moveFloat((IFloat *)out, (IFloat *)temp, temp_size );

  // Below is the same original postconditioning (dare to write again for clarity)
    zmobius_dslash_4(temp, gauge_field, out, CHKB_EVEN, DAG_NO, mobius_arg, mass);
    zmobius_zvectTimesEquComplex(temp, mobius_arg->zmobius_kappa_b,
				   local_ls, ls_stride, s_node_coor);
    zmobius_m5inv(odd_out, temp, mass, DAG_NO, mobius_arg,
		  mobius_arg->zmobius_kappa_ratio);
    zmobius_m5inv(temp, odd_in, mass, DAG_NO, mobius_arg,
		  mobius_arg->zmobius_kappa_ratio);
    odd_out->VecAddEquVec(temp, temp_size); 
    break;
  }
  case ZMOB_PC_ORIG: {
    zmobius_dslash_4(temp, gauge_field, out, CHKB_EVEN, DAG_NO, mobius_arg, mass);
    zmobius_zvectTimesEquComplex(temp, mobius_arg->zmobius_kappa_b,
				   local_ls, ls_stride, s_node_coor);
    zmobius_m5inv(odd_out, temp, mass, DAG_NO, mobius_arg,
		  mobius_arg->zmobius_kappa_ratio);
    zmobius_m5inv(temp, odd_in, mass, DAG_NO, mobius_arg,
		  mobius_arg->zmobius_kappa_ratio);
    odd_out->VecAddEquVec(temp, temp_size); 
    break;}
  default:
    ERR.NotImplemented(cname,fname);
    break;
  }

  
  VRB.Sfree(cname, fname, "temp2", temp2);
  sfree(temp2);
  VRB.Sfree(cname, fname, "temp", temp);
  sfree(temp);

 
  return iter;
}


//------------------------------------------------------------------
// Overloaded function is same as original 
// but true_res=0.
//------------------------------------------------------------------
int DiracOpZMobius::MatInv(Vector *out, Vector *in, PreserveType prs_in)
{ return MatInv(out, in, 0, prs_in); }


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in and out = f_out.
//------------------------------------------------------------------
int DiracOpZMobius::MatInv(Float *true_res, PreserveType prs_in)
{ return MatInv(f_out, f_in, true_res, prs_in); }


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in, out = f_out, true_res=0.
//------------------------------------------------------------------
int DiracOpZMobius::MatInv(PreserveType prs_in)
{ return MatInv(f_out, f_in, 0, prs_in); }


//------------------------------------------------------------------
// Mat(Vector *out, Vector *in) :
// Mat is the unpreconditioned fermion matrix.  
// Mat works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpZMobius::Mat(Vector *out, Vector *in) {  
  char *fname = "Mat(V*,V*)";
  VRB.Func(cname,fname);
  
  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpZMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpZMobius object.
  //----------------------------------------------------------------
  reset_mobius_arg( mobius_lib_arg);
  
  Zmobus *mobius_arg = (Zmobus *) mobius_lib_arg;
  
  //Float kappa = mobius_arg->mobius_kappa_b;
  //Float minus_kappa = -kappa;
  //Float kappa_ratio = mobius_arg->mobius_kappa_b/ mobius_arg->mobius_kappa_c;

  const int local_ls = GJP.SnodeSites();
  const int s_nodes = GJP.Snodes();
  const int global_ls = local_ls * s_nodes;
  const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = ((Zmobus*)mobius_lib_arg)->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  // points to the even part of fermion source 
  Vector *odd_in = (Vector *) ( (IFloat *) in + temp_size );
  // points to the even part of fermion solution
  Vector *odd_out = (Vector *) ( (IFloat *) out + temp_size );
  // temp
  Vector *frm_tmp2 = (Vector *) mobius_arg->frm_tmp2;

  //odd part
  //mobius_dslash_4(out, gauge_field, odd_in, CHKB_EVEN, DAG_NO, mobius_arg, mass);
  zmobius_dslash_4(out, gauge_field, odd_in, CHKB_ODD, DAG_NO, mobius_arg, mass);

#if 0
  out->VecTimesEquFloat(minus_kappa, temp_size); 
#else
  for(int s=0; s<local_ls;++s){
    const Complex* kappa_b= mobius_arg->zmobius_kappa_b;
    int glb_s = s + local_ls*s_node_coor;
    int idx = s*ls_stride/2;// "/2" is for complex
    vecTimesEquComplex((Complex*) out+idx, -kappa_b[glb_s], ls_stride);
  }

#endif
  
#if 0
  // intialize to zero since using the "plus-equal version"
  for(int i=0;i<temp_size;i++){
    *((IFloat*)frm_tmp2+i)=0.0;
  }
  zmobius_dslash_5_plus(frm_tmp2, in, mass, 0, mobius_arg);
  fTimesV1PlusV2((IFloat*)frm_tmp2, kappa_ratio, (IFloat*)frm_tmp2, (IFloat *)in, temp_size);
  out->VecAddEquVec(frm_tmp2, temp_size); 
#else
  {
    Complex* kappa_ratio= mobius_arg->zmobius_kappa_ratio;
    out->VecAddEquVec(in, temp_size); 
    zmobius_kappa_dslash_5_plus_cmplx(out, in, mass, 0, mobius_arg,
				      kappa_ratio);
  }
#endif
  
  //even part
  zmobius_dslash_4(odd_out, gauge_field, in, CHKB_EVEN, DAG_NO, mobius_arg, mass);

#if 0
  odd_out->VecTimesEquFloat(minus_kappa, temp_size); 
#else
  for(int s=0; s<local_ls;++s){
    const Complex* kappa_b= mobius_arg->zmobius_kappa_b;
    int glb_s = s + local_ls*s_node_coor;
    int idx = s*ls_stride/2;// "/2" is for complex
    vecTimesEquComplex((Complex*) odd_out+idx, -kappa_b[glb_s], ls_stride);
  }
#endif

#if 0
  // intialize to zero since using the "plus-equal version"
  for(int i=0;i<temp_size;i++){
    *((IFloat*)frm_tmp2+i)=0.0;
  }
  zmobius_dslash_5_plus(frm_tmp2, odd_in, mass, 0, mobius_arg);
  fTimesV1PlusV2((IFloat*)frm_tmp2, kappa_ratio, (IFloat*)frm_tmp2, (IFloat *)odd_in, temp_size);
  odd_out->VecAddEquVec(frm_tmp2, temp_size);
#else
  {
    Complex* kappa_ratio= mobius_arg->zmobius_kappa_ratio;
    odd_out->VecAddEquVec(odd_in, temp_size); 
    zmobius_kappa_dslash_5_plus_cmplx(odd_out, odd_in, mass, 0, mobius_arg,
				kappa_ratio);
  }

#endif
  
}


void DiracOpZMobius::Dminus(Vector *out, Vector *in) {  
  char *fname = "Dminus(V*,V*)";
  VRB.Func(cname,fname);

  Zmobus *mobius_arg = (Zmobus *) mobius_lib_arg;
  //Float kappa_c_inv_div2 = 0.5*( 2 * (GJP.Mobius_c() *(4 - GJP.DwfHeight()) - GJP.DwfA5Inv()) );
  const int local_ls = GJP.SnodeSites();
  const int s_nodes = GJP.Snodes();
  const int global_ls = local_ls * s_nodes;
  const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = ((Zmobus*)mobius_lib_arg)->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;

  
  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  // points to the odd part of fermion source 
  Vector *odd_in = (Vector *) ( (IFloat *) in + temp_size );
  // points to the odd part of fermion solution
  Vector *odd_out = (Vector *) ( (IFloat *) out + temp_size );

  zmobius_dminus(out, gauge_field, odd_in, CHKB_ODD, DAG_NO, mobius_arg);
  zmobius_dminus(odd_out, gauge_field, in, CHKB_EVEN, DAG_NO, mobius_arg);
  // out = -(c*D_W-1)*in (= 1 for DWF)

#if 0
  fTimesV1PlusV2((IFloat*)out, kappa_c_inv_div2, (IFloat*)in, (IFloat *)out, 2*temp_size);
#else
  for(int s=0; s<local_ls;++s){
    int glb_s = s + local_ls*s_node_coor;
    Complex kappa_c_inv_div2 =
      GJP.ZMobius_c()[glb_s] *(4 - GJP.DwfHeight()) - GJP.DwfA5Inv();
       
    int idx = s*ls_stride/2;// "/2" is for complex
    zTimesV1PlusV2((Complex*)out+idx, kappa_c_inv_div2, (Complex*)in+idx,
		   (Complex*)out+idx, ls_stride);
  }
  for(int s=0; s<local_ls;++s){
    int glb_s = s + local_ls*s_node_coor;
    Complex kappa_c_inv_div2 =
      GJP.ZMobius_c()[glb_s] *(4 - GJP.DwfHeight()) - GJP.DwfA5Inv();
       
    int idx = s*ls_stride/2;// "/2" is for complex
    zTimesV1PlusV2((Complex*)odd_out+idx, kappa_c_inv_div2, (Complex*)odd_in+idx,
		   (Complex*)odd_out+idx, ls_stride);
  }
#endif

  out->VecTimesEquFloat(-1.0, 2*temp_size); 

}


//------------------------------------------------------------------
// MatDag(Vector *out, Vector *in) :
// MatDag is the unpreconditioned fermion matrix.  
// MatDag works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpZMobius::MatDag(Vector *out, Vector *in) {
  char *fname = "MatDag(V*,V*)";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpZMobius constructors but is done again in case the
  // user has modified the GJP.MobiusA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpZMobius object.
  //----------------------------------------------------------------
  Zmobus *mobius_arg = (Zmobus *) mobius_lib_arg;

  reset_mobius_arg( mobius_lib_arg);

  //  Float kappa = mobius_arg->mobius_kappa_b;
  //Float minus_kappa = -kappa;
  //Float kappa_ratio = mobius_arg->mobius_kappa_b/mobius_arg->mobius_kappa_c;


  const int local_ls = GJP.SnodeSites();
  const int s_nodes = GJP.Snodes();
  const int global_ls = local_ls * s_nodes;
  const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = ((Zmobus*)mobius_lib_arg)->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  // points to the even part of fermion source 
  Vector *odd_in = (Vector *) ( (IFloat *) in + temp_size );
  // points to the even part of fermion solution
  Vector *odd_out = (Vector *) ( (IFloat *) out + temp_size );
  // temp
  Vector *frm_tmp2 = (Vector *) mobius_arg->frm_tmp2;

  //odd part

#if 0
  frm_tmp2->VecTimesEquFloat(kappa, temp_size); 
#else
  moveFloat((IFloat*)frm_tmp2,(IFloat*)odd_in, temp_size);
  for(int s=0; s<local_ls;++s){
    const Complex* kappa_b= mobius_arg->zmobius_kappa_b;
    int glb_s = s + local_ls*s_node_coor;
    int idx = s*ls_stride/2;// "/2" is for complex
    vecTimesEquComplex((Complex*) frm_tmp2+idx, -conj(kappa_b[glb_s]), ls_stride);
  }
#endif

  zmobius_dslash_4(out, gauge_field, frm_tmp2, CHKB_ODD, DAG_YES, mobius_arg, mass);

  
#if 0
  // intialize to zero since using the "plus-equal version"
  for(int i=0;i<temp_size;i++){
    *((IFloat*)frm_tmp2+i)=0.0;
  }
  zmobius_dslash_5_plus(frm_tmp2, in, mass, DAG_YES, mobius_arg);
  fTimesV1PlusV2((IFloat*)frm_tmp2, kappa_ratio, (IFloat*)frm_tmp2, (IFloat *)in, temp_size);
  out->VecAddEquVec(frm_tmp2, temp_size); 
#else
  {
    Complex* kappa_ratio= mobius_arg->zmobius_kappa_ratio;
    out->VecAddEquVec(in, temp_size); 
    zmobius_kappa_dslash_5_plus_cmplx(out, in, mass, DAG_YES, mobius_arg,
				kappa_ratio);
  }
#endif
  
  //even part
#if 0
  odd_out->VecTimesEquFloat(kappa, temp_size);
#else
  moveFloat((IFloat*)frm_tmp2,(IFloat*)in, temp_size);
  for(int s=0; s<local_ls;++s){
    const Complex* kappa_b= mobius_arg->zmobius_kappa_b;
    int glb_s = s + local_ls*s_node_coor;
    int idx = s*ls_stride/2;// "/2" is for complex
    vecTimesEquComplex((Complex*) frm_tmp2+idx, -conj(kappa_b[glb_s]), ls_stride);
  }
#endif

  
  //mobius_dslash_4(odd_out, gauge_field, in, CHKB_ODD, DAG_YES, mobius_arg, mass);
  zmobius_dslash_4(odd_out, gauge_field, frm_tmp2, CHKB_EVEN, DAG_YES, mobius_arg, mass);

#if 0
  // intialize to zero since using the "plus-equal version"
  for(int i=0;i<temp_size;i++){
    *((IFloat*)frm_tmp2+i)=0.0;
  }
  zmobius_dslash_5_plus(frm_tmp2, odd_in, mass, DAG_YES, mobius_arg);
  fTimesV1PlusV2((IFloat*)frm_tmp2, kappa_ratio, (IFloat*)frm_tmp2, (IFloat *)odd_in, temp_size);
  odd_out->VecAddEquVec(frm_tmp2, temp_size); 
#else
    {
    Complex* kappa_ratio= mobius_arg->zmobius_kappa_ratio;
    odd_out->VecAddEquVec(odd_in, temp_size); 
    zmobius_kappa_dslash_5_plus_cmplx(odd_out, odd_in, mass, DAG_YES, mobius_arg,
				kappa_ratio);
  }
#endif

}


//------------------------------------------------------------------
// MatHerm(Vector *out, Vector *in) :
// MatHerm is gamma5*R*Mat.
// MatHerm works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpZMobius::MatHerm(Vector *out, Vector *in) {
  char *fname = "MatHerm(V*,V*)";
  VRB.Func(cname,fname);

  ERR.NotImplemented(cname,fname);
  
  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpZMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpZMobius object.
  //----------------------------------------------------------------
  Zmobus *mobius_arg = (Zmobus *) mobius_lib_arg;
  reset_mobius_arg( mobius_lib_arg);  
  
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

void DiracOpZMobius::CalcHmdForceVecs(Vector *chi)
{
  char *fname = "CalcHmdForceVecs(V*)" ;
  ERR.NotImplemented(cname,fname); // TIZB, didn't check for now
#if 0
  VRB.Func(cname,fname) ;

  
  if (f_out == 0)
    ERR.Pointer(cname, fname, "f_out") ;

  if (f_in == 0)
    ERR.Pointer(cname, fname, "f_in") ;

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpZMobius constructors but is done again in case the
  // user has modified the GJP.MobiusA5Inv(), GJP.MobiusHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpZMobius object.
  //----------------------------------------------------------------
  Zmobus *mobius_arg = (Zmobus *) mobius_lib_arg;
  reset_mobius_arg( mobius_lib_arg);

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
    Float kappa = ((Zmobus *)mobius_lib_arg)->mobius_kappa_b ;
    psi->VecTimesEquFloat(-kappa*kappa,f_size_cb) ;
  }

  rho = (Vector *)((Float *)f_out + f_size_cb) ;

  Dslash(rho, chi, CHKB_ODD, DAG_NO) ;
//  fprintf(stderr,"Dslash\n");

  sigma = (Vector *)((Float *)f_in + f_size_cb) ;

  Dslash(sigma, psi, CHKB_ODD, DAG_YES) ;
//  fprintf(stderr,"Dslash\n");

  return ;
#endif
}

//------------------------------------------------------------------
// DiracOpGlbSum(Float *): 
// The global sum used by InvCg. If s_nodes = 1
// it is the usual global sum. If s_nodes > 1 it
// is the 5-dimensional globals sum glb_sum_five.
//------------------------------------------------------------------
void DiracOpZMobius::DiracOpGlbSum(Float *float_p) {
//  if(GJP.Snodes() == 1) {
//    glb_sum(float_p);
//  }
//  else {
    glb_sum_five(float_p);
//  }
}



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
void DiracOpZMobius::RitzMat(Vector *out, Vector *in) {
  char *fname = "RitzMat(V*,V*)";
  VRB.Func(cname,fname);

  //printf("single ritzmat %d\n",dirac_arg->RitzMatOper);
  switch(dirac_arg->RitzMatOper)
    {
    case MATPCDAG_MATPC:
      MatPcDagMatPc(out, in);
      break;
    case MATPCDAG_MATPC_SHIFT:
      MatPcDagMatPcShift(out, in);
      break;
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
void DiracOpZMobius::RitzMat(Vector *out, Vector *in,
			 MatrixPolynomialArg* cheby_arg) {
  char *fname = "RitzMat(V*,V*,MatrixPolyArg*)";
  VRB.Func(cname,fname);

  double time_start=dclock();
    
  //debug const
  const Zmobus *mobius_arg = (Zmobus *) mobius_lib_arg;
  const Float shift = dirac_arg->eigen_shift;

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


//!! N.B. This overwrites contents of  mobius_arg->frm_tmp2
void DiracOpZMobius::MatPcHerm(Vector *out, Vector *in) {
  char *fname = "MatPcHerm(V*,V*)";
  ERR.NotImplemented(cname,fname);
  
  VRB.Func(cname,fname);
  Zmobus *mobius_arg = (Zmobus *) mobius_lib_arg;
  Vector* vtmp = (Vector*)(mobius_arg->frm_tmp1);
  
  MatPc(vtmp,in);
  ReflectAndMultGamma5( out, vtmp,  mobius_arg->vol_4d/2, mobius_arg->ls);
  
}
#if 0
// specific to dwf 
void ReflectAndMultGamma5( Vector *out, const Vector *in,  int nodevol, int ls)
{
  char *fname = "MultGamma5(V*,V*,i)";
  ERR.NotImplemented("d_op_zmobius.C",fname);
  VRB.Func("",fname);
  for(int s=0; s< ls; ++s) { 
    IFloat *p = (IFloat *)out + 24*nodevol*s;
    IFloat *q = (IFloat *)in + 24*nodevol*(ls-1-s);
    for(int n = 0; n < nodevol; ++n)
      {
	int i;
	for(i = 0; i < 12; ++i)
	  *p++ = *q++;
	
	for(i = 0; i < 12; ++i)
	  *p++ = - *q++;
      }
  }

}

void HermicianDWF_ee( Vector* vtmp, Vector* evec, Float mass, Lattice* lattice, Vector* Apsi )
{
  ERR.NotImplemented("d_op_zmobius.C","HermicianDWF_ee");
  CgArg cg_arg;
  cg_arg.mass = mass;
  cg_arg.RitzMatOper = MATPC_HERM; // could be MATPCDAG_MATPC;
  DiracOpZMobius dop( *lattice, 0, 0, &cg_arg, CNV_FRM_NO );
  
  dop. MatPc(Apsi, evec);
  ReflectAndMultGamma5( vtmp, Apsi,  
			GJP.VolNodeSites()/2, GJP.SnodeSites() );
}
#endif

CPS_END_NAMESPACE
