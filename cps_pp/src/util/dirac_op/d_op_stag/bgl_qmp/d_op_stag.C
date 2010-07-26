//------------------------------------------------------------------
//
// d_op_stag.C
//
// DiracOpStag is derived from the DiracOpStagTypes class.
// DiracOpStag is the front end for a library that contains
// all Dirac operators associated with Staggered fermions.
//
//------------------------------------------------------------------

#include <stdio.h>
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/vector.h>
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/time_cps.h>
#include <util/stag.h>
#include <comms/cbuf.h>
#include <comms/glb.h>
#include <comms/scu.h>
#if TARGET == QCDOC
#include <qalloc.h>
#include <ppc_lib.h>
#endif
CPS_START_NAMESPACE

extern "C"{
  void vaxmy(Float *scale,Vector *mult,Vector *sub,int ncvec);
  void vaxmy_vxdot(Float *scale, Vector *mult, Vector *sub, int ncvec, Float *norm);
}

//extern "C" void dirac_comm_assert(void);
//extern SCUDirArgIR *SCUarg;
//extern SCUDirArgIR *SCUarg_1;

//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
DiracOpStag::DiracOpStag(Lattice & latt,
			 Vector *f_field_out,
			 Vector *f_field_in,
			 CgArg *arg,
			 CnvFrmType cnv_frm_flg) :
			 DiracOpStagTypes(latt, 
					  f_field_out,
					  f_field_in, 
					  arg,
					  cnv_frm_flg)
{
  cname = "DiracOpStag";
  char *fname = "DiracOpStag(L&,V*,V*,CgArg*,CnvFrmType)";
  VRB.Func(cname,fname);


  //----------------------------------------------------------------
  // Do the necessary conversions
  //----------------------------------------------------------------
  if(cnv_frm == CNV_FRM_YES)
    lat.Convert(STAG, f_out, f_in);
  else if(cnv_frm == CNV_FRM_NO)
    lat.Convert(STAG);

  //----------------------------------------------------------------
  // Make a copy of gauge fields and rearrange them
  // so that it is suitable for Dirac operation. Here we
  // Assume that only one DiracStagOp instance exists at any time,
  // otherwise we need a static data member to check initializtions
  //----------------------------------------------------------------
#if 0
  stag_dirac_init(latt.GaugeField());
#endif
  stag_dirac_init_g();
  stag_dirac_init_with_lat (lat);

  //----------------------------------------------------------------
  // Set the node checkerboard size of the fermion field
  //----------------------------------------------------------------
  f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / 2;


  //----------------------------------------------------------------
  // Allocate memory for the temporary fermion vector frm_tmp.
  //----------------------------------------------------------------
#if 1
  frm_tmp = (Vector *) fmalloc(cname,fname,"frm_tmp",f_size_cb * sizeof(Float));
#else
  frm_tmp = (Vector *) qalloc(QCOMMS|QFAST,f_size_cb * sizeof(Float));
  if(frm_tmp == 0){
    frm_tmp = (Vector *) qalloc(QCOMMS,f_size_cb * sizeof(Float));
    printf("frm_tmp is allocated int DDR (%p)\n",frm_tmp);
  }
  if(frm_tmp == 0)
    ERR.Pointer(cname,fname, "frm_tmp");
  VRB.Smalloc(cname,fname, "frm_tmp", 
	      frm_tmp, f_size_cb * sizeof(Float));
#endif
  


  //----------------------------------------------------------------
  // Initialize parameters
  //----------------------------------------------------------------
  DiracArg(dirac_arg);

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
DiracOpStag::~DiracOpStag() {
  char *fname = "~DiracOpStag()";
  VRB.Func(cname,fname);

#if 0
  stag_destroy_dirac_buf();
#endif
  stag_destroy_dirac_buf_g();

  if(cnv_frm == CNV_FRM_YES)
    lat.Convert(CANONICAL, f_out, f_in);
  else if(cnv_frm == CNV_FRM_NO)
    lat.Convert(CANONICAL);

  //----------------------------------------------------------------
  // Free memory
  //----------------------------------------------------------------
  VRB.Sfree(cname,fname, "frm_tmp", frm_tmp);
#if TARGET == QCDOC
  qfree(frm_tmp);
#else
  sfree(frm_tmp);
#endif
}


//------------------------------------------------------------------
// DiracArg(CgArg *arg):
// It sets the dirac_arg pointer to arg and initializes
// mass_sq = 4 * mass^2.
//------------------------------------------------------------------
  void DiracOpStag::DiracArg(CgArg *arg){
    dirac_arg = arg;

    // Added for anisotropic lattices
    //------------------------------------------------------------------
    mass_rs = dirac_arg->mass * GJP.XiBare()/GJP.XiV();
    mass_sq = 4 * mass_rs * mass_rs;
    // End modification

  }


//------------------------------------------------------------------
// MatPcDagMatPc :
// MatPcDagMatPc is the fermion matrix that appears in the HMC 
// evolution. It is a Hermitian matrix where M is
// the even/odd preconditioned Dirac Operator matrix.        
// MatPcDagMatPc connects only even-->even sites.
// The in, out fields are defined on the even checkerboard.
// If dot_prd is not 0 then the dot product (on node)
// <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.
//------------------------------------------------------------------
void DiracOpStag::MatPcDagMatPc(Vector *out, 
			       Vector *in, 
			       Float *dot_prd){
 static long nflops = (4416)*GJP.VolNodeSites();

#undef PROFILE
#ifdef PROFILE
  struct timeval start,end;
  gettimeofday(&start,NULL);
  CGflops=0;
#endif
  stag_dirac(frm_tmp, in, 0, 0);
  stag_dirac(out, frm_tmp, 1, 0);

  if( dot_prd !=0 ){
	vaxmy_vxdot(&mass_sq,in,out,f_size_cb/6,dot_prd);
    CGflops +=f_size_cb*4;
    nflops +=f_size_cb*4;
  } else {
    CGflops +=f_size_cb*2;
    nflops +=f_size_cb*2;
    vaxmy(&mass_sq,in,out,f_size_cb/6);
  }

#ifdef PROFILE
  gettimeofday(&end,NULL);
//  printf("DiracOpStag::MatPcDagMatPc:: ");
  print_flops(cname,fname,CGflops,&start,&end);
#endif
}


//------------------------------------------------------------------
// Dslash(Vector *out, Vector *in, ChkbType cb, DagType dag) :
// Dslash is the derivative part of the fermion matrix. 
// Dslash conects only odd-->even or even-->odd sites.
// The in, out fields are defined on a checkerboard.
// cb refers to the checkerboard of the in field.
//------------------------------------------------------------------
void DiracOpStag::Dslash(Vector *out, 
				  Vector *in, 
				  ChkbType cb, 
				  DagType dag) {

  stag_dirac(out,in, int(cb),int(dag));
}

//------------------------------------------------------------------
// Dslash(Vector *out, Vector *in, ChkbType cb, DagType dag, int dir_flag) :
// Dslash is the derivative part of the fermion matrix. 
// Dslash conects only odd-->even or even-->odd sites.
// The in, out fields are defined on a checkerboard.
// cb refers to the checkerboard of the in field.
// dir_flag is flag which takes value 0 when all direction contribute to D,
// 1 - when only the special anisotropic direction contributes to D,
// 2 - when all  except the special anisotropic direction. 
//------------------------------------------------------------------
void DiracOpStag::Dslash(Vector *out, 
			 Vector *in, 
			 ChkbType cb, 
			 DagType dag,
			 int dir_flag) {
  // FIXME: Is it wise to simply ignore dir_flag?  -- Xiao-Yong Jin
  // if (dir_flag != 0)
  //   {
  //     ERR.NotImplemented (cname, "Dslash (.. dir_flag != 0)");
  //   }
  stag_dirac (out, in, (int) cb, (int) dag);
}

#if 0
//------------------------------------------------------------------
// Dslash(Vector *out, Vector *in, ChkbType cb, DagType dag, int dir_flag) :
// Dslash is the derivative part of the fermion matrix. 
// Dslash conects only odd-->even or even-->odd sites.
// The in, out fields are defined on a checkerboard.
// cb refers to the checkerboard of the in field.
// dir_flag is flag which takes value 0 when all direction contribute to D,
// 1 - when only the special anisotropic direction contributes to D,
// 2 - when all  except the special anisotropic direction. 
//------------------------------------------------------------------
void DiracOpStag::Dslash(Vector *out, 
			 Vector *in, 
			 ChkbType cb, 
			 DagType dag,
			 int dir_flag) {

  const int VECT_LEN=6;

  int nx[4];
  int nb[4];
  int xv[3];
  
  //-----------------------------------------------------------
  //  nx[4]
  //-----------------------------------------------------------
  nx[0] = GJP.XnodeSites();
  nx[1] = GJP.YnodeSites();
  nx[2] = GJP.ZnodeSites();
  nx[3] = GJP.TnodeSites();
  
  //-----------------------------------------------------------
  //  nb[4]
  //-----------------------------------------------------------
  nb[0] = 4;
  nb[1] = nb[0]*nx[0];
  nb[2] = nb[1]*nx[1];
  nb[3] = nb[2]*nx[2];
  
  //-----------------------------------------------------------
  //  xv[3]
  //-----------------------------------------------------------
  xv[0] = nx[3]/2;
  xv[1] = (nx[3]*nx[0])/2;
  xv[2] = (nx[3]*nx[0]*nx[1])/2;
  
  int x[4], offset, nu=GJP.XiDir();

  const Matrix *uoff;
  const Matrix *mp0;

  Vector vtmp0, vtmp1;
  Vector *outoff;  
  Vector *vp0;
  Vector *vp1;
  
  for(x[0] = 0; x[0] < nx[0]; ++x[0]) {
    for(x[1] = 0; x[1] < nx[1]; ++x[1]) {
      for(x[2] = 0; x[2] < nx[2]; ++x[2]) {
	for(x[3] = 0; x[3] < nx[3]; ++x[3]) {
	  
	  int parity = (x[0]+x[1]+x[2]+x[3])%2;
	  if((cb && parity==0) || (cb==0 && parity)){
	    
	    uoff = (Matrix *) gauge_field +
	      nb[0]*x[0]+nb[1]*x[1]+nb[2]*x[2]+nb[3]*x[3];
	    outoff = out+(x[3]>>1)+xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	    
	    int iter = 0;
	    for (int mu = 0; mu < 4; ++mu) {
	      if ( (dir_flag==0) || 
		   (dir_flag==1 && mu==nu) || 
		   (dir_flag==2 && mu!=nu))  {
		
		//-------------------------------------------
		//  calculate U^dag_u(x) in(x+u)
		//-------------------------------------------
		if(x[mu] == nx[mu]-1) { 	// x+mu off node
		  x[mu] = 0;
		  offset=(x[3]>>1)+xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
		  getPlusData((IFloat *)&vtmp0, (IFloat *) (in+offset),
			      VECT_LEN, mu);
		  x[mu] = nx[mu]-1;
		  vp0 = &vtmp0;
		  
		} else { 			// x+mu on node
		  x[mu]++;
		  offset=(x[3]>>1)+xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
		  vp0 = in+offset;
		  x[mu]--;
		}
		
		mp0 = uoff+mu;
		
		if(iter == 0 && dag == 0)
		  uDagDotXEqual((IFloat *)outoff, (const IFloat *)mp0,
				(const IFloat *)vp0);
		else
		  uDagDotXPlus((IFloat *)outoff,(const IFloat *)mp0,
			       (const IFloat *)vp0);
		
		
		
		//-------------------------------------------
		//  calculate U_u(x-u) in(x-u)
		//-------------------------------------------
		if(x[mu] == 0) { 		// x-mu off node
		  x[mu] = nx[mu]-1;
		  mp0 = uoff+x[mu]*nb[mu]+mu;
		  offset=(x[3]>>1)+xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
		  vp1 = in+offset;
		  x[mu] = 0;
		  
		  uDotXEqual((IFloat *)&vtmp0, (const IFloat *)mp0,
			     (const IFloat *)vp1);
		  
		  getMinusData((IFloat *)&vtmp1, (IFloat *)&vtmp0,
			       VECT_LEN, mu);
		  
		  *outoff -= vtmp1;
		  
		} else { 			// x-mu on node
		  
		  x[mu]--;
		  mp0 = uoff-nb[mu]+mu;
		  offset=(x[3]>>1)+xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
		  vp0 = in+offset;
		  x[mu]++;
		  
		  uDotXMinus((IFloat *)outoff,(const IFloat *)mp0,
			     (const IFloat *)vp0);
		}
		
		iter++;
		
	      }
	    }
	  }
	}
      }
    }
  }
}

#endif

//------------------------------------------------------------------
// dMdmu(Vector *out, Vector *in, ChkbType cb, DagType dag, int order) :
// dMdmu is the derivative of the fermion matrix with respect to the 
// chemical potential.
// dMdmu conects only odd-->even or even-->odd sites.
// The in, out fields are defined on a checkerboard.
// cb refers to the checkerboard of the in field.
// order refers to the order of the derivative.
//------------------------------------------------------------------
// void DiracOpStag::dMdmu(Vector *out, 
// 				  Vector *in, 
// 				  ChkbType cb, 
// 				  DagType dag,
// 		                  int order) {

//   stag_dMdmu(out,in, int(cb),int(dag),int(order));
// }

//------------------------------------------------------------------
// int MatInv(Vector *out, Vector *in, 
//            Float *true_res, PreserveType prs_in);
// The inverse of the Dirac Operator (D+m)
// using Conjugate gradient.
// Assume: the vector in contains both even and odd src.
//               even part is the 1st part.
// Return: the vector out contains both even and odd solutions.
//       the even solution is the 1st part.
// If true_res !=0 the value of the true residual is returned
// in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// prs_in is not used. The source in is always preserved.
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int DiracOpStag::MatInv(Vector *out, 
			Vector *in, 
			Float *true_res,
			PreserveType prs_in) {
  char *fname = "MatInv(V*,V*,F*)";
  VRB.Func(cname,fname);
  timeval start,end;
  int vol = GJP.VolNodeSites();

  Vector *k_e = in;
  Vector *k_o = k_e + vol/2;

  Vector *tmp = (Vector *) smalloc(f_size_cb * sizeof(Float));
  if(tmp == 0)
    ERR.Pointer(cname,fname, "tmp");
  VRB.Smalloc(cname,fname, "tmp", 
	      tmp, f_size_cb * sizeof(Float));

  // tmp = (2m - D)k

  stag_dirac(tmp, k_o, 1, 0);
  fTimesV1MinusV2((IFloat *)tmp, 2.*mass_rs, (IFloat *)k_e,
  	(IFloat *)tmp, f_size_cb);

#define PROFILE
#ifdef PROFILE
  gettimeofday(&start,NULL);
#endif
  int iter = InvCg(out, tmp, true_res);
#ifdef PROFILE
  gettimeofday(&end,NULL);
//  printf("DiracOpStag::InvCg:: ");
  print_flops(cname,fname,DiracOp::CGflops,&start,&end);
#endif

  // calculate odd solution
  Vector *x_e = out;
  Vector *x_o = x_e+vol/2;
  moveMem((IFloat *)x_o, (IFloat *)k_o, f_size_cb*sizeof(Float));
  stag_dirac(tmp, x_e, 0, 0);
  vecMinusEquVec((IFloat *)x_o, (IFloat *)tmp, f_size_cb);
  vecTimesEquFloat((IFloat *)x_o, 0.5/mass_rs, f_size_cb);

  sfree(tmp);

  return iter;
}


//------------------------------------------------------------------
// Overloaded function is same as original 
// but true_res=0.
//------------------------------------------------------------------
int DiracOpStag::MatInv(Vector *out, Vector *in, PreserveType prs_in)
{ return MatInv(out, in, 0, prs_in); }


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in and out = f_out.
//------------------------------------------------------------------
int DiracOpStag::MatInv(Float *true_res, PreserveType prs_in)
{ return MatInv(f_out, f_in, true_res, prs_in); }


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in, out = f_out, true_res=0.
//------------------------------------------------------------------
int DiracOpStag::MatInv(PreserveType prs_in)
{ return MatInv(f_out, f_in, 0, prs_in); }

//------------------------------------------------------------------
// RitzMat(Vector *out, Vector *in) :
// RitzMat is the base operator used in in Ritz.
// RitzMat works on the full or half lattice.
// The in, out fields are defined on the full or half lattice.
//------------------------------------------------------------------
void DiracOpStag::RitzMat(Vector *out, Vector *in) {
  char *fname = "RitzMat(V*,V*)";
  VRB.Func(cname,fname);
  Float *dot=0;

  Float mass = dirac_arg->mass;
  Float c = 1.0/(64.0 + 4.0*mass*mass);

  switch(dirac_arg->RitzMatOper)
    {
    case MATDAG_MAT:
      MatPcDagMatPc(out, in, dot);
      break;

    case MATPCDAG_MATPC:
      MatPcDagMatPc(out, in, dot);
      break;
      
    case NEG_MATPCDAG_MATPC:
      MatPcDagMatPc(out, in, dot);
      out->VecNegative(out, RitzLatSize());
      break;
      
    case NEG_MATDAG_MAT:
      MatPcDagMatPc(out, in, dot);
      out->VecNegative(out, RitzLatSize());
      break;
      
    case MATDAG_MAT_NORM:
      MatPcDagMatPc(out, in, dot);
      out->VecTimesEquFloat(c,RitzLatSize());
      break;

    case NEG_MATDAG_MAT_NORM:
      MatPcDagMatPc(out, in, dot);
      out->VecTimesEquFloat(-c,RitzLatSize());
      break;

    default:
      ERR.General(cname,fname,"RitzMatOper %d not implemented",
		  dirac_arg->RitzMatOper);
    }
  
}

//------------------------------------------------------------------
// RitzEigMat(Vector *out, Vector *in) :
// RitzEigMat is the base operator used in in RitzEig.
// RitzEigMat works on the full or half lattice.
// The in, out fields are defined on the full or half lattice.
//------------------------------------------------------------------
void DiracOpStag::RitzEigMat(Vector *out, Vector *in) {
  ERR.NotImplemented(cname,"RitzEigMat");
}
CPS_END_NAMESPACE
