#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of GimprOLSym class methods.

  $Id: g_impr_OLSym.C,v 1.2 2003-07-24 16:53:54 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:54 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_impr_OLSym/g_impr_OLSym.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Id: g_impr_OLSym.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.5  2001/08/16 10:50:36  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.3  2001/08/01 15:17:31  anj
//  Minor changes to ensure painless compilarion.Anj
//
//  Revision 1.2  2001/06/19 18:13:26  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:09  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: g_impr_OLSym.C,v $
//  $Revision: 1.2 $
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
#include <stdlib.h> // exit()
#include <math.h>
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/gw_hb.h>
#include <comms/nga_reg.h>
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
    {
      VRB.Warn(cname,fname,"Anisotropic version not implemented\n") ;
      exit(-1) ;
    }

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
  moveMem((IFloat *)mp2, (IFloat *)u_off+BANK4_BASE+BANK_SIZE,
          MATRIX_SIZE*sizeof(IFloat)) ;

  //---------------------------------------------------------------------------
  // force = -beta/3.0*U_mu(x)*gen_stap
  //---------------------------------------------------------------------------
  mDotMEqual((IFloat *)&force, (const IFloat *)mp2, (const IFloat *)mp1);

  vecTimesEquFloat((IFloat *)&force, minus_beta_over_3, MATRIX_SIZE);

  mp1->Dagger((IFloat *)&force);
  force.TrLessAntiHermMatrix(*mp1);
}


//-----------------------------------------------------------------------------
// EvolveMomGforce(Matrix *mom, Float step_size):
// It evolves the canonical momentum mom by step_size
// using the pure gauge force.
//-----------------------------------------------------------------------------
void GimprOLSym::EvolveMomGforce(Matrix *mom, Float step_size){
  char *fname = "EvolveMomGforce(M*,F)";
  VRB.Func(cname,fname);
  
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
    }
  }
}


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
    The staple sum around the link U_\mu(x) is
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
