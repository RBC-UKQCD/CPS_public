#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
//
// d_op_wilsonTm.C
//
// DiracOpWilsonTm is derived from the DiracOpWilson base class. 
// DiracOpWilson is  contains the modifications 
// associated with twisted-mass Wilson fermions.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/dirac_op.h>
#include <util/lattice.h>	// lattice includes util/enum.h; util/data_types.h;
							// vector.h
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/wilson.h>
#include <comms/glb.h>
CPS_START_NAMESPACE

//! Access to the elements of the \e SU(3) matrix
/*!
  Gets the element of the \e SU(3) matrix \e u with row \e row,
  column \e col and complex component \e d
*/
#define U(r,row,col,d,n,cb) *(u+(r+2*(row+3*(col+3*(d+4*(n+vol[0]*(cb)))))))
//! Access to the elements of a spinor vector.
/*!
  Gets the element of the spinor \e psi with spin \e s,
  colour \e c and complex component \e r
  PSIA is the spinor address; PSIV is the value at PSIA
*/
#define PSIA(f,r,c,s,n)     (f)+(r+2*(c+3*(s+4*(n))))
#define PSIV(f,r,c,s,n)     *((f)+(r+2*(c+3*(s+4*(n)))))

//
// "extern 'C'" is used here, as in many other places, to hide 
// within a C routine type conversion from Vector * to Ifloat *
// that is not permitted in C++
//
extern "C" {
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*!
  \param re The real part of the scalar complex factor
  \param im The imaginary part of the scalar complex factor 
  \param a an address of a (real, imaginary) pair
 */
inline void cTimesC(IFloat *a, IFloat re, IFloat im)
{
        // a points to real part
	IFloat t;               
	t = (*a);                         // save real part
	*a = re * (*a) - im * *(a+1);     // real part
	*(a+1)   = re * *(a+1) + im * t;  // imag part
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*! 
	 multiply input vector by gamma_5(theta)
    input vector starts at in and has vol # of color spinors 
    gamma_5(theta) = ctheta + i stheta gamma_5 
    gamma_5 = diag(1,1,-1,-1)						    
    Loop over sites							    
*/
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void g5theta(Vector *in, int vol, IFloat ctheta, IFloat stheta) {

  int n, c;
  for(n=0;n<vol;n++){
     for(c=0;c<3;c++){
		cTimesC(PSIA(((IFloat *)in),0,c,0,n), ctheta, stheta);
		cTimesC(PSIA(((IFloat *)in),0,c,1,n), ctheta, stheta);
		cTimesC(PSIA(((IFloat *)in),0,c,2,n), ctheta, -stheta);
		cTimesC(PSIA(((IFloat *)in),0,c,3,n), ctheta, -stheta);
		}		// for (c=0; . . . 		
	}			// for (vol=0; . . . 
};				// void g5theta . . . 

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
//
// DiracOpWilson constructor does most necessary initialization up
//
DiracOpWilsonTm::DiracOpWilsonTm(Lattice & latt,
			     Vector *f_field_out,
			     Vector *f_field_in,
			     CgArg *arg,
			     CnvFrmType convert) :
			     DiracOpWilson(latt, 
						f_field_out,
						f_field_in, 
						arg,
						convert)
{
  cname = "DiracOpWilsonTm";
  char *fname = "DiracOpWilsonTm(L&,V*,V*,CgArg*,CnvFrmType)";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Re-initializes parameters to include epsilon that 
  // were first initialized in DiracOpWilson
  //----------------------------------------------------------------
  DiracArg(dirac_arg);
  // ~~temporary message
  /* printf("===>  created - DiracOpWilsonTm \n"); */
}

//
// DiracOpWilson destructor does necessary cleaning up
//
DiracOpWilsonTm::~DiracOpWilsonTm() {
  char *fname = "~DiracOpWilsonTm()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// DiracArg(CgArg *arg): It sets the dirac_arg pointer to arg and 
// initializes kappa as the absolute value of the imaginary mass
// and the gamma_5(theta) parameters ctheta and sthera
// mass is represented in vml files as mass + \imath epsilon, and in 
// this function in polar form as \kappa (cos (theta) + \gamma_5 sin (theta))
//------------------------------------------------------------------
  void DiracOpWilsonTm::DiracArg(CgArg *arg){

    dirac_arg = arg;
	 epsilon = dirac_arg->epsilon;
	 kappa = (dirac_arg->mass + 4.0) * (dirac_arg->mass + 4.0) + epsilon * epsilon;
	 kappa = 1.0 / sqrt(kappa);
	 ctheta = (dirac_arg->mass + 4.0) * kappa;
	 stheta = epsilon * kappa;
    kappa = (1.0 / 2.0) * kappa;
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
//
// Calls up to standard Dslash from DiracOpWilson
//
void DiracOpWilsonTm::Dslash(Vector *out, 
			   Vector *in, 
			   ChkbType cb, 
			   DagType dag) {

/*	    wilson_dslash((IFloat *)out, 
			(IFloat *)gauge_field, 
			(IFloat *)in, 
			int(cb),
			int(dag),
			(Wilson *)wilson_lib_arg); */ 

	DiracOpWilson::Dslash(out, in, cb, dag); 

}

//
// Dslash_tm is Dslash pre or post multiplied by 
// gamma_5(theta) according to whether DAG is YES
// or NO, respectively
//
void DiracOpWilsonTm::Dslash_tm(Vector *out, 
			   Vector *in, 
			   ChkbType cb, 
			   DagType dag) {

/*--------------------------------------------------------------------------*/
/* Initializations                                                          */
/*--------------------------------------------------------------------------*/
 int vol =  ((Wilson *)wilson_lib_arg)->vol[0];

/*--------------------------------------------------------------------------*/
/* if DAG_YES:  *out <-- Dslash [gamma_5(theta) * [*in]]			    */
/*--------------------------------------------------------------------------*/
if (dag == DAG_YES) g5theta(in, vol, ctheta, stheta);

Dslash(out, in, cb, dag); 

/*--------------------------------------------------------------------------*/
/* the following redundant construction is to remind that -                 */
/* for DAG_NO:  *out <-- gamma_5(-theta) * [Dslash [*in]]                   */
/* but must restore *in calling parameter to its original values, so        */
/* for DAG_YES, *in <-- gamma_5(-theta) * gamma_5(theta) * [*in]            */
/*---------------------------------------------------------------------------*/
if (dag == DAG_NO) 
	g5theta(out, vol, ctheta, -stheta);
else
	g5theta(in, vol, ctheta, -stheta);

}
//------------------------------------------------------------------
/*!
  The vectors are defined on odd parity lattice sites.
  \param out The resulting vector.
  \param in The vector to be multiplied.
   Calls Dslash_tm not Dslash otherwise the same as 
   DiracOpWilson::MatPc
*/
//------------------------------------------------------------------
void DiracOpWilsonTm::MatPc(Vector *out, Vector *in) {  
  
  Vector *tmp1;
  int vol;
  int r, c, s, n;

/*--------------------------------------------------------------------------*/
/* Initializations                                                          */
/*--------------------------------------------------------------------------*/
  vol =  ((Wilson *)wilson_lib_arg)->vol[0];
  tmp1 = (Vector *)(((Wilson *)wilson_lib_arg)->af[0]);
//printf("Yep, been here, matpc\n");
//int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;
//tmp1 = (Vector *) smalloc(temp_size * sizeof(Float));

  
/*--------------------------------------------------------------------------*/
/* Dslash_E0                                                                */
/*--------------------------------------------------------------------------*/
  Dslash_tm(tmp1, in, CHKB_ODD, DAG_NO);

/*--------------------------------------------------------------------------*/
/* Dslash_0E                                                                */
/*--------------------------------------------------------------------------*/
  Dslash_tm(out, tmp1, CHKB_EVEN, DAG_NO);

/*--------------------------------------------------------------------------*/
/* 1_OO - kappa^2 * Dslash_0E * Dslash_E0                                   */
/*--------------------------------------------------------------------------*/
  for(n=0;n<vol;n++){
    for(s=0;s<4;s++){
      for(c=0;c<3;c++){
	for(r=0;r<2;r++){
		PSIV(((IFloat *)out),r,c,s,n) = 
		PSIV(((IFloat *)in),r,c,s,n) - kappa*kappa * PSIV(((IFloat *)out),r,c,s,n);
	}
      }
    }
  }
//sfree(tmp1);
}

//------------------------------------------------------------------
/*!
  Multiplication by the odd-even preconditioned fermion matrix.
  The vectors are defined on odd parity lattice sites.
  \param out The resulting vector.
  \param in The vector to be multiplied.
  Calls Dslash_tm not Dslash otherwise the same as 
   DiracOpWilson::MatPcDag
*/
//------------------------------------------------------------------
void DiracOpWilsonTm::MatPcDag(Vector *out, Vector *in) {

  Vector *tmp1;
  int vol;
  int r, c, s, n;

/*--------------------------------------------------------------------------*/
/* Initializations                                                          */
/*--------------------------------------------------------------------------*/
  vol =  ((Wilson *)wilson_lib_arg)->vol[0];
  tmp1 = (Vector *)(((Wilson *)wilson_lib_arg)->af[0]);
//printf("Yep, been here, matpcdag\n");
//int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;
//tmp1 = (Vector *) smalloc(temp_size * sizeof(Float));
  
/*--------------------------------------------------------------------------*/
/* Dslash_E0                                                                */
/*--------------------------------------------------------------------------*/
  Dslash_tm(tmp1, in, CHKB_ODD, DAG_YES);

/*--------------------------------------------------------------------------*/
/* Dslash_0E                                                                */
/*--------------------------------------------------------------------------*/
  Dslash_tm(out, tmp1, CHKB_EVEN, DAG_YES);

/*--------------------------------------------------------------------------*/
/* 1_OO - kappa^2 * Dslash_0E * Dslash_E0                                   */
/*--------------------------------------------------------------------------*/
  for(n=0;n<vol;n++){
    for(s=0;s<4;s++){
      for(c=0;c<3;c++){
	for(r=0;r<2;r++){
		PSIV(((IFloat *)out),r,c,s,n) = 
		PSIV(((IFloat *)in),r,c,s,n) - kappa*kappa * PSIV(((IFloat *)out),r,c,s,n);
	}
      }
    }
  }
//sfree(tmp1);
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
// Calls Dslash_tm not Dslash otherwise the same as 
// DiracOpWilson::MatPcDagMatPC
//------------------------------------------------------------------
void DiracOpWilsonTm::MatPcDagMatPc(Vector *out, 
					 Vector *in, 
					 Float *dot_prd){

  Vector *tmp1;
  Vector *tmp2;
  IFloat sum;
  int vol;
  int r, c, s, n;


/*--------------------------------------------------------------------------*/
/* Initializations                                                          */
/*--------------------------------------------------------------------------*/
//printf("Yep, been here\n");
  vol =  ((Wilson *)wilson_lib_arg)->vol[0];
  tmp1 = (Vector *)(((Wilson *)wilson_lib_arg)->af[0]);
  tmp2 = (Vector *)(((Wilson *)wilson_lib_arg)->af[1]);

//int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;
//tmp1 = (Vector *) smalloc(temp_size * sizeof(Float));
//tmp2 = (Vector *) smalloc(temp_size * sizeof(Float));

/*--------------------------------------------------------------------------*/
/* Dslash_E0                                                                */
/*--------------------------------------------------------------------------*/
  Dslash_tm(tmp1, in, CHKB_ODD, DAG_NO);

/*--------------------------------------------------------------------------*/
/* Dslash_0E                                                                */
/*--------------------------------------------------------------------------*/
  Dslash_tm(tmp2, tmp1, CHKB_EVEN, DAG_NO);

/*--------------------------------------------------------------------------*/
/* 1_OO - kappa^2 * Dslash_0E * Dslash_E0                                   */
/* and compute dot product                                                  */
/*--------------------------------------------------------------------------*/
  sum = 0.0;
  for(n=0;n<vol;n++){
    for(s=0;s<4;s++){
      for(c=0;c<3;c++){
	for(r=0;r<2;r++){
	   PSIV(((IFloat *)tmp2),r,c,s,n) = 
	   PSIV(((IFloat *)in),r,c,s,n) - kappa*kappa * PSIV(((IFloat *)tmp2),r,c,s,n);
	   if(dot_prd != 0)
		sum = sum + PSIV(((IFloat *)tmp2),r,c,s,n) * PSIV(((IFloat *)tmp2),r,c,s,n);
	}
      }
    }
  }
  if(dot_prd != 0) *dot_prd = sum;

/*--------------------------------------------------------------------------*/
/* DslashDag_E0                                                             */
/*--------------------------------------------------------------------------*/
  Dslash_tm(tmp1, tmp2, CHKB_ODD, DAG_YES);

/*--------------------------------------------------------------------------*/
/* DslashDag_0E                                                             */
/*--------------------------------------------------------------------------*/
  Dslash_tm(out, tmp1, CHKB_EVEN, DAG_YES);

/*--------------------------------------------------------------------------*/
/* [1_OO - kappa * DslashDag_0E * DslashDag_E0] *                           */
/*                                 [1_OO - kappa^2 * Dslash_0E * Dslash_E0] */
/*--------------------------------------------------------------------------*/
  for(n=0;n<vol;n++){
    for(s=0;s<4;s++){
      for(c=0;c<3;c++){
	for(r=0;r<2;r++){
		PSIV(((IFloat *)out),r,c,s,n) = 
		PSIV(((IFloat *)tmp2),r,c,s,n) - kappa*kappa * PSIV(((IFloat *)out),r,c,s,n);
	}
      }
    }
  }
//sfree(tmp2);
//sfree(tmp1);

}


//------------------------------------------------------------------
// int MatInv(Vector *out, Vector *in, 
//            Float *true_res, PreserveType prs_in);
// The inverse of the unconditioned Dirac Operator 
// using Conjugate gradient.
//------------------------------------------------------------------
//
// Function MatInv is not (yet) implemented for twisted mass.
// It calls Dslash itself but no calls have been found.  Correct twisted 
// mass modification depends on calling context.
// 
//
int DiracOpWilsonTm::MatInv(Vector *out, 
			  Vector *in, 
			  Float *true_res,
			  PreserveType prs_in) {
  char *fname = "MatInv(V*,V*,F*)";
  VRB.Func(cname,fname);

  ERR.General(cname,fname,"d_op_wilsonTm: MatInv not implemented\n");

  return 0;
}

//------------------------------------------------------------------
// Overloaded function is same as original 
// but true_res=0.
//------------------------------------------------------------------
int DiracOpWilsonTm::MatInv(Vector *out, Vector *in, PreserveType prs_in)
{ return MatInv(out, in, 0, prs_in); }

//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in and out = f_out.
//------------------------------------------------------------------
int DiracOpWilsonTm::MatInv(Float *true_res, PreserveType prs_in)
{ return MatInv(f_out, f_in, true_res, prs_in); }

//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in, out = f_out, true_res=0.
//------------------------------------------------------------------
int DiracOpWilsonTm::MatInv(PreserveType prs_in)
{ return MatInv(f_out, f_in, 0, prs_in); }

//------------------------------------------------------------------

//
// Function Mat is not (yet) implemented for twisted mass.
// It calls Dslash itself but no calls have been found.  Correct twisted 
// mass modification depends on calling context.
//
void DiracOpWilsonTm::Mat(Vector *out, Vector *in) {  
  char *fname = "Mat(V*,V*)";
  VRB.Func(cname,fname);

  ERR.General(cname,fname,"d_op_wilsonTm: Mat not implemented\n");
}

//
// Function MatDag is not (yet) implemented for twisted mass.
// It calls Dslash itself but no calls have been found.  Correct twisted 
// mass modification depends on calling context.
//

void DiracOpWilsonTm::MatDag(Vector *out, Vector *in) {
  char *fname = "MatDag(V*,V*)";
  VRB.Func(cname,fname);

  ERR.General(cname,fname,"d_op_wilsonTm: MatDag not implemented\n");
}

//------------------------------------------------------------------

//
// Function MatHerm is not (yet) implemented for twisted mass.
// It calls Dslash itself but no calls have been found.  Correct twisted 
// mass modification depends on calling context.
//
void DiracOpWilsonTm::MatHerm(Vector *out, Vector *in) {
  char *fname = "MatHerm(V*,V*)";
  VRB.Func(cname,fname);

  ERR.General(cname,fname,"d_op_wilsonTm: MatHerm not implemented\n");
}

//------------------------------------------------------------------
/*!
  this function calculates force vectors for fermion force.
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

void DiracOpWilsonTm::CalcHmdForceVecs(Vector *chi)
{
  char *fname = "CalcHmdForceVecs(V*)" ;
  VRB.Func(cname,fname) ;

  if (f_out == 0)
    ERR.Pointer(cname, fname, "f_out") ;

  if (f_in == 0)
    ERR.Pointer(cname, fname, "f_in") ;

//------------------------------------------------------------------
// f_out stores (chi,rho) (the v of eqn B12): 
// rho = gamma_5(-theta) Dslash chi
// f_in stores (psi,sigma) (the w of eqn B11):
// psi = gamma_5(-theta) MatPc chi
// sigma = gamma_5(-theta) Dslash psi
//------------------------------------------------------------------

  Vector *chi_new, *rho, *psi, *sigma ;

  int vol =  ((Wilson *)wilson_lib_arg)->vol[0];
  int f_size_cb = 12 * GJP.VolNodeSites() ;

  chi_new = f_out ;

  chi_new->CopyVec(chi, f_size_cb) ;

  psi = f_in ;

  MatPc(psi,chi) ;
  g5theta(psi, vol, ctheta, stheta);

  psi->VecTimesEquFloat(-kappa*kappa,f_size_cb) ;

  rho = (Vector *)((Float *)f_out + f_size_cb) ;

  Dslash(rho, chi, CHKB_ODD, DAG_NO) ;
  g5theta(rho, vol, ctheta, -stheta);

  sigma = (Vector *)((Float *)f_in + f_size_cb) ;

  Dslash(sigma, psi, CHKB_ODD, DAG_YES) ;
  g5theta(sigma, vol, ctheta, stheta);

  return ;
}

//
// alternative form with sensible variable names
/*
void DiracOpWilsonTm::CalcHmdForceVecs(Vector *chi)
{
  char *fname = "CalcHmdForceVecs(V*)" ;
  VRB.Func(cname,fname) ;

  if (f_out == 0) ERR.Pointer(cname, fname, "f_out") ;
  if (f_in == 0) ERR.Pointer(cname, fname, "f_in") ;

  // v1, v2 correspond to usage in f_wilsonTm:  v1 = f_out & v2 = f_in
  // even, odd checkerboard parts of v1, v2
  Vector *v1even, *v2even,  *v1odd, *v2odd;

  int vol =  ((Wilson *)wilson_lib_arg)->vol[0];
  // size in Floats
  int f_size_cb = 12 * GJP.VolNodeSites() ;

//---------------------------------------------------------
// define v1, v2 in terms of DiracOp( ..., f_out, f_in, ... )
  v1odd = f_out ;
  v2odd = f_in ;
  // cast v1 to Float; add f_size_cb in Floats; cast to Vector*
  v1even = (Vector *)((Float *)v1odd + f_size_cb) ;
  v2even = (Vector *)((Float *)v2odd + f_size_cb) ;

//---------------------------------------------------------
// find odd, even cb's of v1  
  v1odd->CopyVec(chi, f_size_cb) ; // ?? uses size in IFLOATS

  Dslash(v1even, v1odd, vol, CHKB_ODD, DAG_NO) ;
  g5theta(v1even, vol, ctheta, -stheta);

//---------------------------------------------------------
// find odd, even cb's of v2
  MatPc(v2odd, v1odd)
  g5theta(v2odd, vol, ctheta, stheta);

  v2odd->VecTimesEquFloat(-kappa*kappa,f_size_cb) ; // ?? uses size in IFLOATS

  Dslash(v2even, v2odd, CHKB_ODD, DAG_YES) ;
  g5theta(v2even, vol, ctheta, stheta);

  return ;
}
*/


//------------------------------------------------------------------
/*!
  this function calculates force vectors for boson force
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

void DiracOpWilsonTm::CalcBsnForceVecs(Vector *chi, Vector *phi)
{
  char *fname = "CalcBsnForceVecs(V*)" ;
  VRB.Func(cname,fname) ;

  if (f_out == 0) ERR.Pointer(cname, fname, "f_out") ;
  if (f_in == 0) ERR.Pointer(cname, fname, "f_in") ;

  /* printf("\n  ~~~~~~~~~~d_op_wilsonTm::CalcBsnForceVecs:  kappa = %e\n", kappa); */ 

  // v1, v2 correspond to usage in f_wilsonTm:  v1 = f_out & v2 = f_in
  // even, odd checkerboard parts of v1, v2
  Vector *v1even, *v2even,  *v1odd, *v2odd;

  int vol =  ((Wilson *)wilson_lib_arg)->vol[0];
  // size in Floats
  int f_size_cb = 12 * GJP.VolNodeSites() ;

//---------------------------------------------------------
// define v1, v2 in terms of DiracOp( ..., f_out, f_in, ... )

  v1odd = f_out ;
  v2odd = f_in ;
  
  // cast v1 to Float; add f_size_cb in Floats; cast to Vector*
  v1even = (Vector *)((Float *)v1odd + f_size_cb) ;
  v2even = (Vector *)((Float *)v2odd + f_size_cb) ;

//---------------------------------------------------------
// find odd, even cb's of v1
  v1odd->CopyVec(chi, f_size_cb) ; // ?? uses size in IFLOATS

  Dslash(v1even, v1odd, CHKB_ODD, DAG_NO) ;
  g5theta(v1even, vol, ctheta, -stheta);

//---------------------------------------------------------
// find odd, even cb's of v2
  v2odd->CopyVec(phi, f_size_cb) ; // ?? uses size in IFLOATS
  g5theta(v2odd, vol, ctheta, stheta);

  v2odd->VecTimesEquFloat(-kappa*kappa,f_size_cb) ; // ?? uses size in IFLOATS

  Dslash(v2even, v2odd, CHKB_ODD, DAG_YES) ;
  g5theta(v2even, vol, ctheta, stheta);

  return ;
}


CPS_END_NAMESPACE

