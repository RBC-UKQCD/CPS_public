#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of GimprOLSym class methods.

  $Id: g_impr_OLSym.C,v 1.12 2013-04-05 17:51:14 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:51:14 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_impr_OLSym/g_impr_OLSym.C,v 1.12 2013-04-05 17:51:14 chulwoo Exp $
//  $Id: g_impr_OLSym.C,v 1.12 2013-04-05 17:51:14 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.12 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_impr_OLSym/g_impr_OLSym.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//-----------------------------------------------------------------------------
//
// g_impr_OLSym.C
//
// GimprOLSym is derived from Lattice. It implements the
// One Loop Symanzik improved gauge action.
//
//-----------------------------------------------------------------------------

CPS_END_NAMESPACE
#include <math.h>
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/gw_hb.h>
#include <util/error.h>
//#include <comms/nga_reg.h>
#include <comms/glb.h>
#include <comms/cbuf.h>
CPS_START_NAMESPACE

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
enum { MATRIX_SIZE = 18 };


//-----------------------------------------------------------------------------
// static variables used only inside this file
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//  CRAM temp buffer
//-----------------------------------------------------------------------------
#ifdef _TARTAN
static Matrix *mp0 = (Matrix *)CRAM_SCRATCH_ADDR;	// ihdot
static Matrix *mp1 = mp0 + 1;
static Matrix *mp2 = mp1 + 1;
//static Matrix *mp3 = mp2 + 1;
#else
static Matrix mt0;
static Matrix mt1;
static Matrix mt2;
//static Matrix mt3;
static Matrix *mp0 = &mt0;		// ihdot
static Matrix *mp1 = &mt1;
static Matrix *mp2 = &mt2;
//static Matrix *mp3 = &mt3;
#endif 


//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------
GimprOLSym::GimprOLSym()
{
  cname = "GimprOLSym";
  char *fname = "GimprOLSym()";
  VRB.Func(cname,fname);

  // This action cannot be used with anisotropic lattices
  if(GJP.XiBare()!=1)
      ERR.NotImplemented(cname,fname,"Anisotropic version not implemented\n") ;

  // See hep-lat/9507010 page 5 eq. 5 through 9
  Float u0(GJP.u0()) ;
  Float Alpha_s = -4.0*log(u0)/3.06839 ;
  // The plaquett comes with coefficient 1.0 
  rect_coeff = -( 1.0 + 0.4805*Alpha_s) /(20.0*u0*u0) ;
  cube_coeff = -(      0.03325*Alpha_s) /(   u0*u0) ;
  minus_beta_over_3 = -GJP.Beta()/3.0 ;
}


//-----------------------------------------------------------------------------
// Destructor
//-----------------------------------------------------------------------------
GimprOLSym::~GimprOLSym()
{
  char *fname = "~GimprOLSym()";
  VRB.Func(cname,fname);
}


//-----------------------------------------------------------------------------
// GclassType Gclass(void):
// It returns the type of gauge class
//-----------------------------------------------------------------------------
GclassType GimprOLSym::Gclass(void){

    return G_CLASS_IMPR_OLSYM; // WARNING! Define somewhere G_CLASS_IMPR_OLSYM

}

const unsigned CBUF_MODE4 = 0xcca52112;


//-----------------------------------------------------------------------------
// GforceSite(Matrix& force, int *x, int mu):
// It calculates the gauge force at site x and direction mu.
//-----------------------------------------------------------------------------
void GimprOLSym::GforceSite(Matrix& force, int *x, int mu)
{
  char *fname = "GforceSite(M&,i*,i)";
  VRB.Func(cname,fname);

  setCbufCntrlReg(4, CBUF_MODE4);

  Matrix *u_off = GaugeField()+GsiteOffset(x)+mu;

  //---------------------------------------------------------------------------
  //  get generalized staple
  //     mp1 = staple
  //---------------------------------------------------------------------------
  AllStaple(*mp1, x, mu);	

  //---------------------------------------------------------------------------
  // mp2 = U_mu(x)
  //---------------------------------------------------------------------------
  moveMem((IFloat *)mp2, (IFloat *)u_off,
          MATRIX_SIZE*sizeof(IFloat)) ;

  //---------------------------------------------------------------------------
  // force = -beta/3.0*U_mu(x)*gen_stap
  //---------------------------------------------------------------------------
  mDotMEqual((IFloat *)&force, (const IFloat *)mp2, (const IFloat *)mp1);

  vecTimesEquFloat((IFloat *)&force, minus_beta_over_3, MATRIX_SIZE);

  mp1->Dagger((IFloat *)&force);
  force.TrLessAntiHermMatrix(*mp1);
}

#if TARGET !=QCDOC
//-----------------------------------------------------------------------------
// EvolveMomGforce(Matrix *mom, Float step_size):
// It evolves the canonical momentum mom by step_size
// using the pure gauge force.
//-----------------------------------------------------------------------------
ForceArg GimprOLSym::EvolveMomGforce(Matrix *mom, Float dt){
  char *fname = "EvolveMomGforce(M*,F)";
  VRB.Func(cname,fname);

  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;
  
  setCbufCntrlReg(4, CBUF_MODE4);

  int x[4];
  
  for(x[0] = 0; x[0] < GJP.XnodeSites(); ++x[0])
  for(x[1] = 0; x[1] < GJP.YnodeSites(); ++x[1])
  for(x[2] = 0; x[2] < GJP.ZnodeSites(); ++x[2])
  for(x[3] = 0; x[3] < GJP.TnodeSites(); ++x[3]) {

    int uoff = GsiteOffset(x);

    for (int mu = 0; mu < 4; ++mu) {
      GforceSite(*mp0, x, mu);

      IFloat *ihp = (IFloat *)(mom+uoff+mu);
      IFloat *dotp = (IFloat *)mp0;
      fTimesV1PlusV2(ihp, dt, dotp, ihp, 18);
      Float norm = ((Matrix*)dotp)->norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);
    }
  }

  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();

  return ForceArg(dt*L1, dt*sqrt(L2), dt*Linf);

}
#endif


//-----------------------------------------------------------------------------
// Float GhamiltonNode(void):
// The pure gauge Hamiltonian of the node sublattice.
//-----------------------------------------------------------------------------
Float GimprOLSym::GhamiltonNode(void){
  char *fname = "GhamiltonNode()";
  VRB.Func(cname,fname);

  Float sum ;
  sum  =              SumReTrPlaqNode() ;
  sum += rect_coeff * SumReTrRectNode() ;
  sum += cube_coeff * SumReTrCubeNode() ;

  sum *= minus_beta_over_3 ;

  return sum ;
}


//-----------------------------------------------------------------------------
// void GactionGradient(Matrix &grad, int *x, int mu)
//   Calculates the partial derivative of the gauge action
//   w.r.t. the link U_mu(x).
//-----------------------------------------------------------------------------
void GimprOLSym::GactionGradient(Matrix &grad, int *x, int mu)
{
  char *fname = "GactionGradient(M&,i*,i)" ;
  VRB.Func(cname, fname) ;

  //---------------------------------------------------------------------------
  // get generalized staple
  //---------------------------------------------------------------------------
  AllStaple(grad, x, mu) ;

  //---------------------------------------------------------------------------
  // grad = -beta/3.0 * generalized_staple
  //---------------------------------------------------------------------------
  vecTimesEquFloat((IFloat *)&grad, minus_beta_over_3, MATRIX_SIZE) ;
}


//------------------------------------------------------------------------
// GimprOLSym::AllStaple() gets all the staples relavent to the heat bath
//    in this case it's 
//   plaq_staple+rect_coeff*rect_staple+cube_coeff*chair_staple
/*!
    The staple sum around the link \f$ U_\mu(x) \f$ is
  \f[
\sum_{\pm \nu, |\nu|\neq \mu} \left\{\right.
  U_\nu(x+\mu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]\f[
  -\frac{ 1 + 0.4805 \alpha_s}{20 u_0^2} \left[\right.
 U_\mu(x+\mu) U_\nu(x+2\mu) U^\dagger_\mu(x+\mu+\nu)
 U^\dagger_\mu(x+\nu)  U^\dagger_\nu(x) \f]\f[

 + U_\nu(x+\mu)    U^\dagger_\mu(x+\nu)  U^\dagger_\mu(x-\mu+\nu)
 U^\dagger_\nu(x-\mu) U_\mu(x-\mu)  \f]\f[
 
 + U_\nu(x+\mu)    U_\nu(x+\mu+\nu)   U^\dagger_\mu(x+2\nu)
 U^\dagger_\nu(x+\nu)  U^\dagger_\nu(x)
 \left.\right]\f]\f[

 -\frac{ 0.03325 \alpha_s}{u_0^2}
 \sum_{\pm \rho, |\rho|\neq \nu, |\rho|\neq \mu}

 U_\nu(x+\mu) U_\rho(x+\mu+\nu) U_{-\mu}(x+\mu+\nu+\rho) U_{-\nu}(x+\nu+\rho) U_{-\rho}(x+\rho)
 
 \left.\right\}
 \f]
 where \f$u_0\f$ is the tadpole coefficient
 and \f$\alpha_s = -4\log(u_0)/3.06839\f$.
 
  \param stap The computed staple sum.
  \param x The coordinates of the lattice site 
  \param mu The link direction 
  
 */
//------------------------------------------------------------------------
void GimprOLSym::AllStaple(Matrix &stap, const int *x, int mu){
  char * fname = "AllStaple()"; 
  VRB.Func(cname, fname);

  Matrix mat;

  // The 3-staple
  //Staple(mat, x, mu );  
  BufferedStaple(mat, x, mu );  

  // The flat 5-staple (flat rectangle)
  BufferedRectStaple(stap, x, mu); 
  vecTimesEquFloat((IFloat*)&stap, rect_coeff, MATRIX_SIZE);

  // add them and store the result in stap
  vecAddEquVec((IFloat*)&stap, (IFloat*)&mat, MATRIX_SIZE);

  // the 3D 5-staple (twisted chair)
  BufferedCubeStaple(mat, x, mu );  
  vecTimesEquFloat((IFloat*)&mat, cube_coeff, MATRIX_SIZE);

  // add the result to stap
  vecAddEquVec((IFloat*)&stap, (IFloat*)&mat, MATRIX_SIZE); 
}


CPS_END_NAMESPACE
