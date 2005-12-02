#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of GimprRect class.

  $Id: g_impr_rect.C,v 1.8 2005-12-02 16:27:40 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2005-12-02 16:27:40 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_impr_rect/g_impr_rect.C,v 1.8 2005-12-02 16:27:40 chulwoo Exp $
//  $Id: g_impr_rect.C,v 1.8 2005-12-02 16:27:40 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.8 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_impr_rect/g_impr_rect.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdlib.h>
#include <math.h>
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/gw_hb.h>
#include <util/time.h>
#include <comms/nga_reg.h>
#include <comms/glb.h>
#include <comms/cbuf.h>
CPS_START_NAMESPACE
#ifdef _TARTAN
CPS_END_NAMESPACE
#include <stdlib.h> // exit()
CPS_START_NAMESPACE
#endif


//------------------------------------------------------------------------------
enum { MATRIX_SIZE = 18 };


//------------------------------------------------------------------------------
// static variables used only inside this file
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//  CRAM temp buffer
//------------------------------------------------------------------------------
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


//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
GimprRect::GimprRect()
{
  cname = "GimprRect";
  char *fname = "GimprRect()";
  VRB.Func(cname,fname);

// This action cannot be used with anisotropic lattices
  if(GJP.XiBare()!=1)
    {
      VRB.Warn(cname,fname,"Anisotropic version not implemented\n") ;
      exit(-1) ;
    }

  plaq_coeff = - GJP.Beta() * ( 1.0 - 8.0 * GJP.C1() ) / 3.0 ;
  rect_coeff = - GJP.Beta() * (             GJP.C1() ) / 3.0 ;
}


//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
GimprRect::~GimprRect()
{
  char *fname = "~GimprRect()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------------------
// GclassType Gclass(void):
// It returns the type of gauge class
//------------------------------------------------------------------------------
GclassType GimprRect::Gclass(void){
  return G_CLASS_IMPR_RECT;
}

unsigned GimprRect::CBUF_MODE4 = 0xcca52112;


//------------------------------------------------------------------------------
// GforceSite(Matrix& force, int *x, int mu):
// It calculates the gauge force at site x and direction mu.
//------------------------------------------------------------------------------
void GimprRect::GforceSite(Matrix& force, int *x, int mu)
{
  char *fname = "GforceSite(M&,i*,i)";
  VRB.Func(cname,fname);

  Float tmp ;

  setCbufCntrlReg(4, CBUF_MODE4);

  Matrix *u_off = GaugeField()+GsiteOffset(x)+mu;

  //----------------------------------------------------------------------------
  //  get staple
  //     mp1 = staple
  //----------------------------------------------------------------------------
  Staple(*mp1, x, mu);	
  ForceFlops += 198*3*3+12+216*3;

  //----------------------------------------------------------------------------
  // mp2 = U_mu(x)
  //----------------------------------------------------------------------------
  moveMem((IFloat *)mp2, (IFloat *)u_off+BANK4_BASE+BANK_SIZE,
          MATRIX_SIZE*sizeof(IFloat)) ;

  //----------------------------------------------------------------------------
  // force = -(beta*(1-8*c_1)/3)*U_mu(x)*stap
  //----------------------------------------------------------------------------
  mDotMEqual((IFloat *)&force, (const IFloat *)mp2, (const IFloat *)mp1);

  tmp = plaq_coeff ;
  vecTimesEquFloat((IFloat *)&force, tmp, MATRIX_SIZE);

  //----------------------------------------------------------------------------
  //  get rectangle staple
  //     mp1 = rect_stap
  //----------------------------------------------------------------------------
  RectStaple(*mp1, x, mu) ;
  ForceFlops += 198*3*18+216*3*6;

  //----------------------------------------------------------------------------
  // mp2 = -(beta*c_1/3)*U_mu(x)
  //----------------------------------------------------------------------------
  moveMem((IFloat *)mp2, (IFloat *)u_off+BANK4_BASE+BANK_SIZE,
          MATRIX_SIZE*sizeof(IFloat));

  tmp = rect_coeff ;
  vecTimesEquFloat((IFloat *)mp2, tmp, MATRIX_SIZE) ;
  ForceFlops +=234;

  //----------------------------------------------------------------------------
  // force += -(beta*c_1/3)*U_mu(x)*rect_stap
  //----------------------------------------------------------------------------
  mDotMPlus((IFloat *)&force, (const IFloat *)mp2, (const IFloat *)mp1);

  mp1->Dagger((IFloat *)&force);
  force.TrLessAntiHermMatrix(*mp1);
  ForceFlops +=198+24;
}

#if TARGET != QCDOC
#define PROFILE
//------------------------------------------------------------------------------
// EvolveMomGforce(Matrix *mom, Float step_size):
// It evolves the canonical momentum mom by step_size
// using the pure gauge force.
//------------------------------------------------------------------------------
Float GimprRect::EvolveMomGforce(Matrix *mom, Float step_size){
  char *fname = "EvolveMomGforce(M*,F)";
  VRB.Func(cname,fname);

  Float Fdt = 0.0;

#ifdef PROFILE
  Float time = -dclock();
  ForceFlops = 0;
#endif
  
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
      fTimesV1PlusV2(ihp, step_size, dotp, ihp+BANK4_BASE, 18);
      Fdt += step_size*step_size*dotProduct(dotp, dotp, 18);
    }
  }
  ForceFlops +=GJP.VolNodeSites()*4*18*2;
#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,ForceFlops,time);
#endif

  glb_sum(&Fdt);

  return sqrt(Fdt);

}
#endif

//------------------------------------------------------------------------------
// Float GhamiltonNode(void):
// The pure gauge Hamiltonian of the node sublattice.
//------------------------------------------------------------------------------
Float GimprRect::GhamiltonNode(void){
  char *fname = "GhamiltonNode()";
  VRB.Func(cname,fname);

  Float sum ;
  sum  = plaq_coeff * SumReTrPlaqNode() ;
  sum += rect_coeff * SumReTrRectNode() ;
  return sum ;
}


//------------------------------------------------------------------------------
// void GactionGradient(Matrix &grad, int *x, int mu)
//   Calculates the partial derivative of the gauge action
//   w.r.t. the link U_mu(x).
//------------------------------------------------------------------------------
void GimprRect::GactionGradient(Matrix &grad, int *x, int mu)
{
  char *fname = "GactionGradient(M&,i*,i)" ;
  VRB.Func(cname, fname) ;

  //----------------------------------------------------------------------------
  // get plaq_staple
  //----------------------------------------------------------------------------
  Staple(*mp1, x, mu) ;

  //----------------------------------------------------------------------------
  // grad = (-beta*(1-8*c_1)/3) * plaq_staple
  //----------------------------------------------------------------------------
  vecTimesEquFloat((IFloat *)&grad, plaq_coeff, MATRIX_SIZE) ;

  //----------------------------------------------------------------------------
  // get rect_staple
  //----------------------------------------------------------------------------
  RectStaple(*mp1, x, mu) ;

  //----------------------------------------------------------------------------
  // mp1 = (-beta*c_1/3) * rect_staple
  //----------------------------------------------------------------------------
  vecTimesEquFloat((IFloat *)mp1, rect_coeff, MATRIX_SIZE) ;

  //----------------------------------------------------------------------------
  // grad = (-beta*(1-8*c_1)/3)*plaq_staple + (-beta*c_1/3)*rect_staple
  //----------------------------------------------------------------------------
  vecAddEquVec((IFloat *)&grad, (const IFloat *)mp1, MATRIX_SIZE) ;
}


//------------------------------------------------------------------------
// GimprRect::AllStaple() gets all the staples relavent to the heat bath
//    in this case it's (1-8*c_1)*plaq_staple + c_1*rect_staple 
/*!
  The sum of the staples around the link \f$ U_\mu(x) \f$ is:
\f[
(1-8c_1)\sum_{\nu \neq \mu}[
              U_\nu(x+\mu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
 +  U^\dagger_\nu(x+\mu-\nu) U^\dagger_\mu(x-\nu) U_\nu(x-\nu)] \f]\f[
 +  c_1\sum_{\nu \neq \mu}\left[\right.     
 U_\mu(x+\mu) U_\nu(x+2\mu) U^\dagger_\mu(x+\mu+\nu)
 U^\dagger_\mu(x+\nu)  U^\dagger_\nu(x) \f]\f[
 + U_\mu(x+\mu)    U^\dagger_\nu(x+2\mu-\nu) U^\dagger_\mu(x+\mu-\nu) 
 U^\dagger_\mu(x-\nu)  U_\nu(x-\nu) \f]\f[
 + U_\nu(x+\mu)    U^\dagger_\mu(x+\nu)  U^\dagger_\mu(x-\mu+\nu)
 U^\dagger_\nu(x-\mu) U_\mu(x-\mu)  \f]\f[
 + U^\dagger_\nu(x+\mu-\nu) U^\dagger_\mu(x-\nu)  U^\dagger_\mu(x-\mu-\nu)
 U_\nu(x-\mu-\nu) U_\mu(x-\mu) \f]\f[
 + U_\nu(x+\mu)    U_\nu(x+\mu+\nu)   U^\dagger_\mu(x+2\nu)
 U^\dagger_\nu(x+\nu)  U^\dagger_\nu(x) \f]\f[
 + U^\dagger_\nu(x+\mu-\nu) U^\dagger_\nu(x+\mu-2\nu) U^\dagger_\mu(x-2\nu)
 U_\nu(x-2\nu)  U_\nu(x-\nu)       
 \left.\right]
 \f]

  \param x The coordinates of the lattice site 
  \param mu The link direction 
  \param stap The computed staple sum.
 */
//------------------------------------------------------------------------
void GimprRect::
AllStaple(Matrix &stap, const int *x, int mu){
  char * fname = "AllStaple()"; 
  VRB.Func(cname, fname);

  Matrix mat;
  BufferedStaple(mat, x, mu );  
  vecTimesEquFloat((IFloat*)&mat, 1-8*GJP.C1(), MATRIX_SIZE);
  BufferedRectStaple(stap, x, mu); 
  vecTimesEquFloat((IFloat*)&stap, GJP.C1(), MATRIX_SIZE);
  vecAddEquVec((IFloat*)&stap, (IFloat*)&mat, MATRIX_SIZE); 
}


CPS_END_NAMESPACE
