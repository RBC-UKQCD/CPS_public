#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of GpowerRect class.

  $Id: g_power_rect.C,v 1.8 2013-04-05 17:51:14 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:51:14 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_power_rect/g_power_rect.C,v 1.8 2013-04-05 17:51:14 chulwoo Exp $
//  $Id: g_power_rect.C,v 1.8 2013-04-05 17:51:14 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.8 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_power_rect/g_power_rect.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------------------
//
// g_power_rect.C
//
// GpowerRect is derived from Lattice and is relevant to the
// action which contains the standard Wilson plaquette operator
// plus the second order rectangle operator plus a power
// plaquette term plus a power rectangle term.
//
// The full action is:
//
// (  Sum_p [ c_0*beta*{ -Tr[U_p]/3} + ( {1 - Tr[U_p]/3} / c_p )^k_p ]
//  + Sum_r [ c_1*beta*{ -Tr[U_r]/3} + ( {1 - Tr[U_r]/3} / c_r )^k_r ] )
// with c_p = GJP.PowerPlaqCutoff(), k_p = GJP.PowerPlaqExponent(),
//      c_r = GJP.PowerRectCutoff(), k_r = GJP.PowerRectExponent(),
// c_0 = 1 - 8 * c_1, c_1 = GJP.C1()
// This action supresses plaquettes with {1 - ReTr[U_p]/3} > c_p
// and rectangles with {1 - ReTr[U_r]/3} > c_r
// and threfore reduces lattice dislocations.
//
//------------------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/gw_hb.h>
//#include <comms/nga_reg.h>
#include <comms/glb.h>
#include <comms/scu.h>
#include <comms/cbuf.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
enum { MATRIX_SIZE = 18 };

//------------------------------------------------------------------
// DRAM temp buffer used for scu transfer
//------------------------------------------------------------------
static Matrix m_tmp1, m_tmp2;

//------------------------------------------------------------------------------
// static variables used only inside this file
//------------------------------------------------------------------------------
static IFloat invs3 = -1./3.;

//------------------------------------------------------------------------------
//  CRAM temp buffer
//------------------------------------------------------------------------------
#ifdef _TARTAN
static Matrix *mp0 = (Matrix *)CRAM_SCRATCH_ADDR;	// ihdot
static Matrix *mp1 = mp0 + 1;
static Matrix *mp2 = mp1 + 1;
static Matrix *mp3 = mp2 + 1;
static Matrix *mp4 = mp3 + 1;
#else
static Matrix mt0;
static Matrix mt1;
static Matrix mt2;
static Matrix mt3;
static Matrix mt4;
static Matrix *mp0 = &mt0;		// ihdot
static Matrix *mp1 = &mt1;
static Matrix *mp2 = &mt2;
static Matrix *mp3 = &mt3;
static Matrix *mp4 = &mt4;
#endif 


//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
GpowerRect::GpowerRect()
{
  cname = "GpowerRect";
  char *fname = "GpowerRect()";
  VRB.Func(cname,fname);

  plaq_coeff = - GJP.Beta() * ( 1.0 - 8.0 * GJP.C1() ) / 3.0 ;
  rect_coeff = - GJP.Beta() * (             GJP.C1() ) / 3.0 ;
}


//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
GpowerRect::~GpowerRect()
{
  char *fname = "~GpowerRect()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------------------
// GclassType Gclass(void):
// It returns the type of gauge class
//------------------------------------------------------------------------------
GclassType GpowerRect::Gclass(void){
  return G_CLASS_POWER_RECT;
}

const unsigned CBUF_MODE2 = 0xcca52112;
const unsigned CBUF_MODE4 = 0xcca52112;


//------------------------------------------------------------------------------
// GforceSite(Matrix& force, int *x, int mu):
// It calculates the gauge force at site x and direction mu.
//------------------------------------------------------------------------------
void GpowerRect::GforceSite(Matrix& force, int *x, int mu)
{
  char *fname = "GforceSite(M&,i*,i)";
  VRB.Func(cname,fname);

  Float tmp ;
  Float coeff = invs3;

  setCbufCntrlReg(4, CBUF_MODE4);

  Matrix *u_off = GaugeField()+GsiteOffset(x)+mu;

  //----------------------------------------------------------------------------
  //  get power staple
  //     mp1 = power staple
  //----------------------------------------------------------------------------
  PowerStaple(*mp1, x, mu);	

  //----------------------------------------------------------------------------
  // mp2 = U_mu(x)
  //----------------------------------------------------------------------------
  moveMem((IFloat *)mp2, (IFloat *)u_off,
          MATRIX_SIZE*sizeof(IFloat)) ;

  //----------------------------------------------------------------------------
  // force = -(1/3)*U_mu(x)*power_stap
  //----------------------------------------------------------------------------
  mDotMEqual((IFloat *)&force, (const IFloat *)mp2, (const IFloat *)mp1);

  tmp = coeff ;
  vecTimesEquFloat((IFloat *)&force, tmp, MATRIX_SIZE);

  //----------------------------------------------------------------------------
  //  get power rectangle staple
  //     mp1 = power_rect_stap
  //----------------------------------------------------------------------------
  PowerRectStaple(*mp1, x, mu) ;

  //----------------------------------------------------------------------------
  // mp2 = -(1/3)*U_mu(x)
  //----------------------------------------------------------------------------
  moveMem((IFloat *)mp2, (IFloat *)u_off,
          MATRIX_SIZE*sizeof(IFloat));

  tmp = coeff ;
  vecTimesEquFloat((IFloat *)mp2, tmp, MATRIX_SIZE) ;

  //----------------------------------------------------------------------------
  // force += -(1/3)*U_mu(x)*power_rect_stap
  //----------------------------------------------------------------------------
  mDotMPlus((IFloat *)&force, (const IFloat *)mp2, (const IFloat *)mp1);

  mp1->Dagger((IFloat *)&force);
  force.TrLessAntiHermMatrix(*mp1);
}


//------------------------------------------------------------------------------
// EvolveMomGforce(Matrix *mom, Float step_size):
// It evolves the canonical momentum mom by step_size
// using the pure gauge force.
//------------------------------------------------------------------------------
ForceArg GpowerRect::EvolveMomGforce(Matrix *mom, Float dt){
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


//------------------------------------------------------------------------------
// Float GhamiltonNode(void):
// The pure gauge Hamiltonian of the node sublattice.
//------------------------------------------------------------------------------
Float GpowerRect::GhamiltonNode(void){
  char *fname = "GhamiltonNode()";
  VRB.Func(cname,fname);

  Float sum ;
  sum  = plaq_coeff * SumReTrPlaqNode() ;
  sum += SumPowerPlaqNode();
  sum += rect_coeff * SumReTrRectNode() ;
  sum += SumPowerRectNode();
  return sum ;
}


//------------------------------------------------------------------------------
// void GactionGradient(Matrix &grad, int *x, int mu)
//   Calculates the partial derivative of the gauge action
//   w.r.t. the link U_mu(x).
//------------------------------------------------------------------------------
void GpowerRect::GactionGradient(Matrix &grad, int *x, int mu)
{
  char *fname = "GactionGradient(M&,i*,i)" ;
  VRB.Func(cname, fname) ;

  Float coeff = invs3;

  //----------------------------------------------------------------------------
  // get power_plaq_staple
  //----------------------------------------------------------------------------
  PowerStaple(*mp1, x, mu) ;

  //----------------------------------------------------------------------------
  // grad = (-1/3) * power_plaq_staple
  //----------------------------------------------------------------------------
  vecTimesEquFloat((IFloat *)&grad, coeff, MATRIX_SIZE) ;

  //----------------------------------------------------------------------------
  // get power_rect_staple
  //----------------------------------------------------------------------------
  PowerRectStaple(*mp1, x, mu) ;

  //----------------------------------------------------------------------------
  // mp1 = (-1/3) * power_rect_staple
  //----------------------------------------------------------------------------
  vecTimesEquFloat((IFloat *)mp1, coeff, MATRIX_SIZE) ;

  //----------------------------------------------------------------------------
  // grad = (-1/3)*[power_plaq_staple + power_rect_staple]
  //----------------------------------------------------------------------------
  vecAddEquVec((IFloat *)&grad, (const IFloat *)mp1, MATRIX_SIZE) ;
}


//------------------------------------------------------------------
/*!
  The power plaquette at site \a x in the \f$ \mu-\nu \f$ plane with
  \f$ \mu < \nu \f$ is
  \f[
  \{(1 - Re Tr[U_p(x,\mu,\nu)]/3) / c \}^k
  \f]
  where the \f$ U_p \f$ is the plaquette and \a c and \a k are the coefficients
  appearing in the action.

  \param x the coordinates of the lattice site at the start of the plaquette
  \param mu The first plaquette direction.
  \param nu The second plaquette direction; should be different from \a mu.
  \return  The real part of the trace of the power plaquette.
*/
//------------------------------------------------------------------
Float GpowerRect::PowerPlaq(int *x, int mu, int nu) const
{
  //  char *fname = "PowerPlaq()";
  //  VRB.Func(cname,fname);
  Float plaq = ReTrPlaq(x, mu, nu);
  Float tmp = (1.0 + invs3 * plaq) / GJP.PowerPlaqCutoff();
  Float p_plaq = 1.0;
  for(int i=0; i < GJP.PowerPlaqExponent(); i++){
    p_plaq *= tmp;
  }

  return p_plaq;
}


//------------------------------------------------------------------
/*!
  The power plaquette at site \a x in the \f$ \mu-\nu \f$ plane with
  \f$ \mu < \nu \f$ is
  \f[
  \{(1 - Re Tr[U_p(x,\mu,\nu)]/3) / c \}^k
  \f]
  where the \f$ U_p \f$ is the plaquette and \a c and \a k are the coefficients
  appearing in the action.

  \return  The real part of the trace of the power plaquette summed over all
  local lattice sites and all six planes.
*/
//------------------------------------------------------------------
Float GpowerRect::SumPowerPlaqNode(void) const
{
  char *fname = "SumPowerPlaqNode() const";
  VRB.Func(cname,fname);
  
  Float sum = 0.0;
  int x[4];
  
  for(x[0] = 0; x[0] < node_sites[0]; ++x[0]) {
    for(x[1] = 0; x[1] < node_sites[1]; ++x[1]) {
      for(x[2] = 0; x[2] < node_sites[2]; ++x[2]) {
	for(x[3] = 0; x[3] < node_sites[3]; ++x[3]) {
	  
	  for (int mu = 0; mu < 3; ++mu) {
	    for(int nu = mu+1; nu < 4; ++nu) {
	      sum += PowerPlaq(x,mu,nu);
	    }
	  }
	}
      }
    }
  }
  
  return sum;
}


//------------------------------------------------------------------
/*!
  The power plaquette at site \a x in the \f$ \mu-\nu \f$ plane with
  \f$ \mu < \nu \f$ is
  \f[
  \{(1 - Re Tr[U_p(x,\mu,\nu)]/3) / c \}^k
  \f]
  where the \f$ U_p \f$ is the plaquette and \a c and \a k are the coefficients
  appearing in the action.

  \return  The real part of the trace of the power plaquette summed over all
  lattice sites and all six planes.
*/
//------------------------------------------------------------------
Float GpowerRect::SumPowerPlaq(void) const
{
  char *fname = "SumPowerPlaq() const";
  VRB.Func(cname,fname);

  Float sum = SumPowerPlaqNode();
  glb_sum(&sum);
  return sum;
}


//------------------------------------------------------------------
/*!
  The staple sum around the link \f$ U_\mu(x) \f$ is
\f[
  \sum_{v \neq u} [
      ps(x,\mu,\nu)     U_\nu(x+\mu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)     
    + ps(x-\nu,\mu,\nu)   U^\dagger_\nu(x+\mu-\nu) U^\dagger_\mu(x-\nu) U_\nu(x-\nu) ]          
\f]
  where
\f[
 ps(x,\mu,\nu) = 
 \beta (1-8c_1) + k/c \{(1 - ReTr[U_p(x,\mu,\nu)]/3) / c \}^{k-1}
 \f]

 where the real parameters \f$\beta,\f$ \a c and \a k are those in the
 definition of the power plaquette action and \f$ U_p \f$ is the standard
 plaquette.

  \param x The coordinates of the lattice site 
  \param mu The link direction 
  \param pstap The computed staple sum.
*/
//------------------------------------------------------------------
void GpowerRect::PowerStaple(Matrix& pstap, int *x, int mu)
{
  char *fname = "PowerStaple(M&,i*,i)";
  VRB.Func(cname,fname);

  Float plaq;
  Float tmp;
  Float ps;
  Float inv_cutoff = 1.0 / GJP.PowerPlaqCutoff();
  int exponent = GJP.PowerPlaqExponent();
  Float coeff = Float(exponent) * inv_cutoff;
  Float beta_c0 = GJP.Beta() * ( 1.0 - 8.0 * GJP.C1() );

  // set cbuf
  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(4, CBUF_MODE4);

  const Matrix *p1;
  int offset_x = GsiteOffset(x);
  Matrix *g_offset = GaugeField()+offset_x;

  for(int nu = 0; nu < 4; ++nu) {
    if(nu != mu) {

      //----------------------------------------------------------
      // mp3 = U_u(x+v)~
      //----------------------------------------------------------
      p1 = GetLinkOld(g_offset, x, nu, mu);
      mp3->Dagger((IFloat *)p1);

      //----------------------------------------------------------
      // p1 = &U_v(x+u)
      //----------------------------------------------------------
      p1 = GetLinkOld(g_offset, x, mu, nu);

      //----------------------------------------------------------
      // mp2 = U_v(x+u) U_u(x+v)~
      //----------------------------------------------------------
      mDotMEqual((IFloat *)mp2, 
		 (const IFloat *)p1,
		 (const IFloat *)mp3);
      
      //----------------------------------------------------------
      //  mp3 = U_v(x)~
      //----------------------------------------------------------
      mp3->Dagger((IFloat *)(g_offset+nu));
      
      //----------------------------------------------------------
      // m_tmp1 = mp2*mp3 = m_tmp1*U_v(x+u)*U_u(x+v)~*U_v(x)~
      //----------------------------------------------------------
      mDotMEqual((IFloat *)&m_tmp1, 
		 (const IFloat *)mp2,
		 (const IFloat *)mp3);

      //----------------------------------------------------------
      // m_tmp2 = U_u(x)*m_tmp1 = U_u(x)*U_v(x+u)*U_u(x+v)~*U_v(x)~
      //----------------------------------------------------------
      mDotMEqual((IFloat *)&m_tmp2, 
		 (const IFloat *)(g_offset+mu),
		 (const IFloat *)&m_tmp1);

      //----------------------------------------------------------
      // Calculate power plaquette staple coefficient
      //----------------------------------------------------------
      plaq = m_tmp2.ReTr();
      tmp = (1.0 + invs3 * plaq) * inv_cutoff;
      ps = 1.0;
      for(int i=0; i < exponent-1; i++){
	ps *= tmp;
      }
      ps = beta_c0 + coeff * ps;

      //----------------------------------------------------------
      // m_tmp1 = m_tmp1*ps = [ U_v(x+u)*U_u(x+v)*~U_v(x)~ ]*ps
      //----------------------------------------------------------
      vecTimesEquFloat((IFloat *)&m_tmp1, ps, MATRIX_SIZE);

      //----------------------------------------------------------
      // Add to the staple
      //----------------------------------------------------------
      if( nu == 0  ||  (mu==0 && nu==1) )
	moveMem((IFloat *)&pstap, 
		(const IFloat *)&m_tmp1, 
		MATRIX_SIZE*sizeof(Float));
      else
	vecAddEquVec((IFloat *)&pstap, 
		     (const IFloat *)&m_tmp1,
		     MATRIX_SIZE);

      

      //----------------------------------------------------------
      //  calculate U_v(x+u-v)~*U_u(x-v)~*U_v(x-v)
      //----------------------------------------------------------
      int off_pv = (x[nu] == 0) ?
              (node_sites[nu]-1)*g_dir_offset[nu]
	    : -g_dir_offset[nu];

      Matrix *g_offpv = g_offset+off_pv;

      //----------------------------------------------------------
      // p1 = U_v(x+u-v)
      // mp3 = U_v(x+u-v)
      //----------------------------------------------------------
      p1 = GetLinkOld(g_offpv, x, mu, nu);
      moveMem((IFloat *)mp3, (IFloat *)p1, 
	      MATRIX_SIZE * sizeof(IFloat));

      //----------------------------------------------------------
      // mp2 = U_u(x-v)*U_v(x+u-v)
      //----------------------------------------------------------
      mDotMEqual((IFloat *)mp2, 
		 (const IFloat *)(g_offpv+mu),
		 (const IFloat *)mp3);
      
      //----------------------------------------------------------
      // mp3 = U_v(x+u-v)~*U_u(x-v)~ = mp2~
      //----------------------------------------------------------
      mp3->Dagger((IFloat *)mp2);


      //----------------------------------------------------------
      // mp2 = U_v(x-v)
      //----------------------------------------------------------
      moveMem((IFloat *)mp2, 
	      (const IFloat *)(g_offpv+nu),
	      MATRIX_SIZE * sizeof(IFloat));
      
      //----------------------------------------------------------
      // m_tmp1 = mp3*mp2 = U_v(x+u-v)~*U_u(x-v)~*U_v(x-v)
      //----------------------------------------------------------
      if(x[nu] == 0) {	// x-v off node
	mDotMEqual((IFloat *)&m_tmp2, 
		   (const IFloat *)mp3,
		   (const IFloat *)mp2);

	// m_tmp1 = U_v(x+u-v)*U_u(x-v)*U_v(x-v)^dag
	getMinusData((IFloat *)&m_tmp1, 
		     (IFloat *)&m_tmp2,
		     MATRIX_SIZE, nu);	
      } else {
	mDotMEqual((IFloat *)&m_tmp1, 
		   (const IFloat *)mp3,
		   (const IFloat *)mp2);
      }

      //----------------------------------------------------------
      // m_tmp2 = U_m(x)*m_tmp1 = U_v(x+u-v)~*U_u(x-v)~*U_v(x-v)*U_u(x)
      //----------------------------------------------------------
      mDotMEqual((IFloat *)&m_tmp2, 
		 (const IFloat *)(g_offset+mu),
		 (const IFloat *)&m_tmp1);

      //----------------------------------------------------------
      // Calculate power plaquette staple coefficient
      //----------------------------------------------------------
      plaq = m_tmp2.ReTr();
      tmp = (1.0 + invs3 * plaq) * inv_cutoff;
      ps = 1.0;
      for(int j=0; j < exponent-1; j++){
	ps *= tmp;
      }
      ps = beta_c0 + coeff * ps;

      //----------------------------------------------------------
      // m_tmp1 = m_tmp1*ps = [ U_v(x+u-v)~*U_u(x-v)~*U_v(x-v) ]*ps
      //----------------------------------------------------------
      vecTimesEquFloat((IFloat *)&m_tmp1, ps, MATRIX_SIZE);

      //----------------------------------------------------------
      // Add to the staple
      //----------------------------------------------------------
      vecAddEquVec((IFloat *)&pstap, 
		   (IFloat *)&m_tmp1,
		   MATRIX_SIZE);

    }
  }
}


//------------------------------------------------------------------
/*!
  The staple sum around the link \f$ U_\mu(x) \f$ is

  \f[
 \sum_{\nu \neq \mu} \left[\right.                                     
R(x,\mu,\nu)  U_\mu(x+\mu)  U_\nu(x+2\mu)  
 U^\dagger_\mu(x+\mu+\nu) U^\dagger_\mu(x+\nu)  U^\dagger_\nu(x)   \f]\f[
 + R(x-\nu,\mu,\nu) U_\mu(x+\mu)  U^\dagger_\nu(x+2\mu-\nu) 
 U^\dagger_\mu(x+\mu-\nu) U^\dagger_\mu(x-\nu)  U_\nu(x-\nu)  \f]\f[
 + R(x-\mu,\mu,\nu) U_\nu(x+\mu)  U^\dagger_\mu(x+\nu) U^\dagger_\mu(x-\mu+\nu)
 U^\dagger_\nu(x-\mu) U_\mu(x-\mu)      \f]\f[
 + R(x-\mu-\nu,\mu,\nu) U^\dagger_\nu(x+\mu-\nu) U^\dagger_\mu(x-\nu)  
 U^\dagger_\mu(x-\mu-\nu) U_\nu(x-\mu-\nu) U_\mu(x-\mu)  \f]\f[
 + R(x+\mu, \nu,\mu) U_\nu(x+\mu)  U_\nu(x+\mu+\nu)  U^\dagger_\mu(x+2\nu)
 U^\dagger_\nu(x+\nu)  U^\dagger_\nu(x) \f]\f[
 + R(x+\mu-2\nu, \nu,\mu) U^\dagger_\nu(x+\mu-\nu) U^\dagger_\nu(x+\mu-2\nu)
  U^\dagger_\mu(x-2\nu) U_\nu(x-2\nu)  U_\nu(x-\nu)
 \left.\right]
 \f]
 where
 \f[
 R(x,\mu, \nu) = c_1+k/c \{(1 - ReTr[U_r(x,\mu,\nu)]/3) / c \}^{k-1}
 \f]
 where \f$U_r\f$ is the 6-link rectangle loop defined in the action.
 
  \param x The coordinates of the lattice site 
  \param mu The link direction 
  \param rect The computed staple sum.
*/
//------------------------------------------------------------------
void GpowerRect::PowerRectStaple(Matrix& rect, int *x, int mu)
{
  char *fname = "PowerRectStaple(M&,i*,i)";
  VRB.Func(cname,fname);

  int link_site[4] ;

  // set CBUF
  setCbufCntrlReg(4, CBUF_MODE4) ;

  //----------------------------------------------------------------------------
  // do a dummy read from the DRAM image controlled by CBUF mode ctrl reg 0
  // to guarantee CBUF will start a new proces  on the next read from
  // a DRAM image controlled by CBUF mode ctrl reg 4.
  //----------------------------------------------------------------------------
#ifdef _TARTAN
  *((unsigned *)mp2) = *((unsigned *)0x2000) ;
#endif

  const Matrix *p1;
  int offset_x = GsiteOffset(x);
  Matrix *g_offset = GaugeField()+offset_x;

  rect.ZeroMatrix() ;

  int i;

  for (i=0; i<4; ++i) link_site[i] = x[i] ;

  Float re_tr_rect;
  Float ps;
  Float tmp;
  Float inv_cutoff = 1.0 / GJP.PowerRectCutoff();
  int exponent = GJP.PowerRectExponent();
  Float coeff = Float(exponent) * inv_cutoff;
  Float beta_c1 = GJP.Beta() * GJP.C1();

  for(int nu = 0; nu < 4; ++nu)
  if(nu != mu) {

    //----------------------------------------------------------
    // mp4 = U_v(x)~
    //----------------------------------------------------------
    mp4->Dagger((IFloat *)GetLink(link_site, nu)) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+v)~
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    mp3->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp2 = U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x+u+v)~
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu));

    //----------------------------------------------------------
    // mp3 = U_u(x+u+v)~ U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp4, (const IFloat *)mp2);

    //----------------------------------------------------------
    // p1 = &U_v(x+2u)
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    --(link_site[nu]) ;
    p1 = GetLink(link_site, nu) ;

    //----------------------------------------------------------
    // mp4 = U_v(x+2u) U_u(x+u+v)~ U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp4, (const IFloat *)p1,
               (const IFloat *)mp3);

    //----------------------------------------------------------
    // mp2 = U_u(x+u)
    //----------------------------------------------------------
    --(link_site[mu]) ;
    moveMem((IFloat *)mp2, (const IFloat *)GetLink(link_site, mu),
            MATRIX_SIZE * sizeof(IFloat)) ;

//PMV

    //----------------------------------------------------------
    // m_tmp1 = mp2*mp4 = 
    // U_u(x+u) U_v(x+2u) U_u(x+u+v)~ U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)&m_tmp1, 
	       (const IFloat *)mp2,
	       (const IFloat *)mp4);

    //----------------------------------------------------------
    // m_tmp2 = U_u(x)*m_tmp1 
    //----------------------------------------------------------
    mDotMEqual((IFloat *)&m_tmp2, 
	       (const IFloat *)(g_offset+mu),
	       (const IFloat *)&m_tmp1);

    //----------------------------------------------------------
    // Calculate power rectangle staple coefficient
    //----------------------------------------------------------
    re_tr_rect = m_tmp2.ReTr();
    tmp = (1.0 + invs3 * re_tr_rect) * inv_cutoff;
    ps = 1.0;
    for(i=0; i < exponent-1; i++){
      ps *= tmp;
    }
    ps = beta_c1 + coeff * ps;

    //----------------------------------------------------------
    // m_tmp1 = m_tmp1*ps
    //----------------------------------------------------------
    vecTimesEquFloat((IFloat *)&m_tmp1, ps, MATRIX_SIZE);

    //----------------------------------------------------------
    // Add to the staple
    //----------------------------------------------------------
    vecAddEquVec((IFloat *)&rect, 
		 (const IFloat *)&m_tmp1,
		 MATRIX_SIZE);

//PMV

    //----------------------------------------------------------
    // mp4 = U_v(x+2u-v)~
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    --(link_site[nu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, nu)) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+u) U_v(x+2u-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x+u-v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp2 = U_u(x+u) U_v(x+2u-v)~ U_u(x+u-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x-v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+u) U_v(x+2u-v)~ U_u(x+u-v)~ U_u(x-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp2 = U_v(x-v)
    //----------------------------------------------------------
    moveMem((IFloat *)mp2, (const IFloat *)GetLink(link_site, nu),
            MATRIX_SIZE * sizeof(IFloat)) ;

//PMV

    //----------------------------------------------------------
    // m_tmp1 = mp3*mp2 = 
    // U_u(x+u) U_v(x+2u-v)~ U_u(x+u-v)~ U_u(x-v)~ U_v(x-v)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)&m_tmp1, 
	       (const IFloat *)mp3,
	       (const IFloat *)mp2);

    //----------------------------------------------------------
    // m_tmp2 = U_u(x)*m_tmp1 
    //----------------------------------------------------------
    mDotMEqual((IFloat *)&m_tmp2, 
	       (const IFloat *)(g_offset+mu),
	       (const IFloat *)&m_tmp1);

    //----------------------------------------------------------
    // Calculate power rectangle staple coefficient
    //----------------------------------------------------------
    re_tr_rect = m_tmp2.ReTr();
    tmp = (1.0 + invs3 * re_tr_rect) * inv_cutoff;
    ps = 1.0;
    for(i=0; i < exponent-1; i++){
      ps *= tmp;
    }
    ps = beta_c1 + coeff * ps;

    //----------------------------------------------------------
    // m_tmp1 = m_tmp1*ps
    //----------------------------------------------------------
    vecTimesEquFloat((IFloat *)&m_tmp1, ps, MATRIX_SIZE);

    //----------------------------------------------------------
    // Add to the staple
    //----------------------------------------------------------
    vecAddEquVec((IFloat *)&rect, 
		 (const IFloat *)&m_tmp1,
		 MATRIX_SIZE);

//PMV

    //----------------------------------------------------------
    // p1 = &U_v(x-2v)
    //----------------------------------------------------------
    --(link_site[nu]) ;
    p1 = GetLink(link_site, nu) ;

    //----------------------------------------------------------
    // mp3 = U_v(x-2v) U_v(x-v)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)p1,
               (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // mp4 = U_u(x-2v)~
    //----------------------------------------------------------
    mp4->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp2 = U_u(x-2v)~ U_v(x-2v) U_v(x-v)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp4, (const IFloat *)mp3) ;

    //----------------------------------------------------------
    // mp4 = U_v(x+u-2v)~
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, nu)) ;

    //----------------------------------------------------------
    // mp3 = U_v(x+u-2v)~ U_u(x-2v)~ U_v(x-2v) U_v(x-v)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp4, (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // mp2 = U_v(x+u-v)~
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    mp2->Dagger((IFloat *)GetLink(link_site, nu)) ;

//PMV

    //----------------------------------------------------------
    // m_tmp1 = mp2*mp3 = 
    // U_v(x+u-v)~ U_v(x+u-2v)~ U_u(x-2v)~ U_v(x-2v) U_v(x-v)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)&m_tmp1, 
	       (const IFloat *)mp2,
	       (const IFloat *)mp3);

    //----------------------------------------------------------
    // m_tmp2 = U_u(x)*m_tmp1 
    //----------------------------------------------------------
    mDotMEqual((IFloat *)&m_tmp2, 
	       (const IFloat *)(g_offset+mu),
	       (const IFloat *)&m_tmp1);

    //----------------------------------------------------------
    // Calculate power rectangle staple coefficient
    //----------------------------------------------------------
    re_tr_rect = m_tmp2.ReTr();
    tmp = (1.0 + invs3 * re_tr_rect) * inv_cutoff;
    ps = 1.0;
    for(i=0; i < exponent-1; i++){
      ps *= tmp;
    }
    ps = beta_c1 + coeff * ps;

    //----------------------------------------------------------
    // m_tmp1 = m_tmp1*ps
    //----------------------------------------------------------
    vecTimesEquFloat((IFloat *)&m_tmp1, ps, MATRIX_SIZE);

    //----------------------------------------------------------
    // Add to the staple
    //----------------------------------------------------------
    vecAddEquVec((IFloat *)&rect, 
		 (const IFloat *)&m_tmp1,
		 MATRIX_SIZE);

//PMV

    //----------------------------------------------------------
    // mp4 = U_u(x-v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp3 =  U_v(x+u-v)~ U_u(x-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x-u-v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp2 = U_v(x+u-v)~ U_u(x-v)~ U_u(x-u-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp4) ;

//  GRF: this code may cause a hang
//
//  //----------------------------------------------------------
//  // p1 = &U_v(x-u-v)
//  //----------------------------------------------------------
//  p1 = GetLink(link_site, nu) ;
//
//  //----------------------------------------------------------
//  // mp3 = U_v(x+u-v)~ U_u(x-v)~ U_u(x-u-v)~ U_v(x-u-v)
//  //----------------------------------------------------------
//  mDotMEqual((IFloat *)mp3, (const IFloat *)mp2,
//             (const IFloat *)p1+BANK4_BASE+BANK_SIZE) ;
//
//  GRF: here is a work-around

    //----------------------------------------------------------
    // mp4 = U_v(x-u-v)
    //----------------------------------------------------------
    moveMem((IFloat *)mp4,
            (const IFloat *)GetLink(link_site, nu),
            MATRIX_SIZE*sizeof(IFloat)) ;

    //----------------------------------------------------------
    // mp3 = U_v(x+u-v)~ U_u(x-v)~ U_u(x-u-v)~ U_v(x-u-v)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

//  GRF: end work-around

    //----------------------------------------------------------
    // mp2 = U_u(x-u)
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    moveMem((IFloat *)mp2, (const IFloat *)GetLink(link_site, mu),
            MATRIX_SIZE * sizeof(IFloat)) ;

//PMV

    //----------------------------------------------------------
    // m_tmp1 = mp3*mp2 = 
    // U_v(x+u-v)~ U_u(x-v)~ U_u(x-u-v)~ U_v(x-u-v) U_u(x-u)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)&m_tmp1, 
	       (const IFloat *)mp3,
	       (const IFloat *)mp2);

    //----------------------------------------------------------
    // m_tmp2 = U_u(x)*m_tmp1 
    //----------------------------------------------------------
    mDotMEqual((IFloat *)&m_tmp2, 
	       (const IFloat *)(g_offset+mu),
	       (const IFloat *)&m_tmp1);

    //----------------------------------------------------------
    // Calculate power rectangle staple coefficient
    //----------------------------------------------------------
    re_tr_rect = m_tmp2.ReTr();
    tmp = (1.0 + invs3 * re_tr_rect) * inv_cutoff;
    ps = 1.0;
    for(i=0; i < exponent-1; i++){
      ps *= tmp;
    }
    ps = beta_c1 + coeff * ps;

    //----------------------------------------------------------
    // m_tmp1 = m_tmp1*ps
    //----------------------------------------------------------
    vecTimesEquFloat((IFloat *)&m_tmp1, ps, MATRIX_SIZE);

    //----------------------------------------------------------
    // Add to the staple
    //----------------------------------------------------------
    vecAddEquVec((IFloat *)&rect, 
		 (const IFloat *)&m_tmp1,
		 MATRIX_SIZE);

//PMV

    //----------------------------------------------------------
    // mp4 = U_v(x-u)~
    //----------------------------------------------------------
    mp4->Dagger((IFloat *)GetLink(link_site, nu)) ;

    //----------------------------------------------------------
    // mp3 = U_v(x-u)~ U_u(x-u)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp4, (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // mp4 = U_u(x-u+v)~
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp2 = U_u(x-u+v)~ U_v(x-u)~ U_u(x-u)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp4, (const IFloat *)mp3) ;

    //----------------------------------------------------------
    // mp4 = U_u(x+v)~
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+v)~ U_u(x-u+v)~ U_v(x-u)~ U_u(x-u)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp4, (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // mp2 = U_v(x+u)
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    --(link_site[nu]) ;
    moveMem((IFloat *)mp2, (const IFloat *)GetLink(link_site, nu),
            MATRIX_SIZE * sizeof(IFloat)) ;


//PMV

    //----------------------------------------------------------
    // m_tmp1 = mp2*mp3 = 
    // U_v(x+u) U_u(x+v)~ U_u(x-u+v)~ U_v(x-u)~ U_u(x-u)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)&m_tmp1, 
	       (const IFloat *)mp2,
	       (const IFloat *)mp3);

    //----------------------------------------------------------
    // m_tmp2 = U_u(x)*m_tmp1 
    //----------------------------------------------------------
    mDotMEqual((IFloat *)&m_tmp2, 
	       (const IFloat *)(g_offset+mu),
	       (const IFloat *)&m_tmp1);

    //----------------------------------------------------------
    // Calculate power rectangle staple coefficient
    //----------------------------------------------------------
    re_tr_rect = m_tmp2.ReTr();
    tmp = (1.0 + invs3 * re_tr_rect) * inv_cutoff;
    ps = 1.0;
    for(i=0; i < exponent-1; i++){
      ps *= tmp;
    }
    ps = beta_c1 + coeff * ps;

    //----------------------------------------------------------
    // m_tmp1 = m_tmp1*ps
    //----------------------------------------------------------
    vecTimesEquFloat((IFloat *)&m_tmp1, ps, MATRIX_SIZE);

    //----------------------------------------------------------
    // Add to the staple
    //----------------------------------------------------------
    vecAddEquVec((IFloat *)&rect, 
		 (const IFloat *)&m_tmp1,
		 MATRIX_SIZE);

//PMV


//  GRF: this code may cause a hang
//
//  //----------------------------------------------------------
//  // p1 = &U_v(x+u+v)
//  //----------------------------------------------------------
//  ++(link_site[nu]) ;
//  p1 = GetLink(link_site, nu) ;
//
//  //----------------------------------------------------------
//  // mp3 = U_v(x+u) U_v(x+u+v)
//  //----------------------------------------------------------
//  mDotMEqual((IFloat *)mp3, (const IFloat *)mp2,
//             (const IFloat *)p1+BANK4_BASE+BANK_SIZE) ;
//
//  GRF: here is a work-around

    //----------------------------------------------------------
    // mp4 = U_v(x+u+v)
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    moveMem((IFloat *)mp4,
            (const IFloat *)GetLink(link_site, nu),
            MATRIX_SIZE*sizeof(IFloat)) ;

    //----------------------------------------------------------
    // mp3 = U_v(x+u) U_v(x+u+v)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

//  GRF: end work-around

    //----------------------------------------------------------
    // mp4 = U_u(x+2v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    ++(link_site[nu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)) ;

    //----------------------------------------------------------
    // mp2 = U_v(x+u) U_v(x+u+v) U_u(x+2v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_v(x+v)~
    //----------------------------------------------------------
    --(link_site[nu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, nu)) ;

    //----------------------------------------------------------
    // mp3 = U_v(x+u) U_v(x+u+v) U_u(x+2v)~ U_v(x+v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp2 = U_v(x)~
    //----------------------------------------------------------
    --(link_site[nu]) ;
    mp2->Dagger((IFloat *)GetLink(link_site, nu)) ;


//PMV

    //----------------------------------------------------------
    // m_tmp1 = mp3*mp2 = 
    // U_v(x+u) U_v(x+u+v) U_u(x+2v)~ U_v(x+v)~  U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)&m_tmp1, 
	       (const IFloat *)mp3,
	       (const IFloat *)mp2);

    //----------------------------------------------------------
    // m_tmp2 = U_u(x)*m_tmp1 
    //----------------------------------------------------------
    mDotMEqual((IFloat *)&m_tmp2, 
	       (const IFloat *)(g_offset+mu),
	       (const IFloat *)&m_tmp1);

    //----------------------------------------------------------
    // Calculate power rectangle staple coefficient
    //----------------------------------------------------------
    re_tr_rect = m_tmp2.ReTr();
    tmp = (1.0 + invs3 * re_tr_rect) * inv_cutoff;
    ps = 1.0;
    for(i=0; i < exponent-1; i++){
      ps *= tmp;
    }
    ps = beta_c1 + coeff * ps;

    //----------------------------------------------------------
    // m_tmp1 = m_tmp1*ps
    //----------------------------------------------------------
    vecTimesEquFloat((IFloat *)&m_tmp1, ps, MATRIX_SIZE);

    //----------------------------------------------------------
    // Add to the staple
    //----------------------------------------------------------
    vecAddEquVec((IFloat *)&rect, 
		 (const IFloat *)&m_tmp1,
		 MATRIX_SIZE);

//PMV

    //----------------------------------------------------------
    // dummy read to switch CBUF banks for looping
    //----------------------------------------------------------
    *((IFloat *)mp4) = *((IFloat *)p1) ;
  }

//VRB.FuncEnd(cname, fname) ;
}


//------------------------------------------------------------------
/*!
  The power rectangle at site \a x in the \f$ \mu-\nu \f$ plane with
  \f$ \mu < \nu \f$ is
  \f[
  \{(1 - Re Tr[U_r(x,\mu,\nu)]/3) / c \}^k
  \f]
  where the \f$ U_r \f$ is the rectangle and \a c and \a k are the coefficients
  appearing in the action.

  \param x the coordinates of the lattice site at the start of the plaquette
  \param mu The first loop direction.
  \param nu The second loop direction; should be different from \a mu.
  \return  The real part of the trace of the power rectangle.
*/
//------------------------------------------------------------------
Float GpowerRect::PowerRect(int *x, int mu, int nu) const
{
//char *fname = "PowerRect(i*,i,i)";
//VRB.Func(cname,fname);


  Float rect = ReTrRect(x, mu, nu);
  Float tmp = (1.0 + invs3 * rect) / GJP.PowerRectCutoff();
  Float p_rect = 1.0;
  for(int i=0; i < GJP.PowerRectExponent(); i++){
    p_rect *= tmp;
  }

  return p_rect;
}


//------------------------------------------------------------------
/*!
  The power rectangle at site \a x in the \f$ \mu-\nu \f$ plane with
  \f$ \mu < \nu \f$ is
  \f[
  \{(1 - Re Tr[U_r(x,\mu,\nu)]/3) / c \}^k
  \f]
  where the \f$ U_r \f$ is the rectangle and \a c and \a k are the coefficients
  appearing in the action.

  \return  The real part of the trace of the power plaquette summed over all
  local lattice sites and all six planes.
*/
//------------------------------------------------------------------
Float GpowerRect::SumPowerRectNode(void) const
{
  char *fname = "SumPowerRectNode()";
  VRB.Func(cname,fname);

  Float sum = 0.0 ;
  int x[4] ;

  for(x[0] = 0; x[0] < node_sites[0]; ++x[0])
  for(x[1] = 0; x[1] < node_sites[1]; ++x[1])
  for(x[2] = 0; x[2] < node_sites[2]; ++x[2])
  for(x[3] = 0; x[3] < node_sites[3]; ++x[3])
  for(int mu = 0; mu < 4; ++mu)
  for(int nu = 0; nu < 4; ++nu) {
    if (mu != nu) {
      sum += PowerRect(x,mu,nu);
    }
  }
  return sum;
}

//------------------------------------------------------------------
/*!
  The power rectangle at site \a x in the \f$ \mu-\nu \f$ plane with
  \f$ \mu < \nu \f$ is
  \f[
  \{(1 - Re Tr[U_r(x,\mu,\nu)]/3) / c \}^k
  \f]
  where the \f$ U_r \f$ is the rectangle and \a c and \a k are the coefficients
  appearing in the action.

  \return  The real part of the trace of the power rectangle summed over all
  lattice sites and all six planes.
*/
//------------------------------------------------------------------
Float GpowerRect::SumPowerRect(void) const
{
  char *fname = "SumPowerRect()";
  VRB.Func(cname,fname);

  Float sum = SumPowerRectNode() ;
  glb_sum(&sum) ;
  return sum ;
}

void GpowerRect::AllStaple(Matrix & stap, const int *x, int mu)
{
  ERR.NotImplemented(cname, "AllStaple(M&,i*,i)") ;
}

CPS_END_NAMESPACE
