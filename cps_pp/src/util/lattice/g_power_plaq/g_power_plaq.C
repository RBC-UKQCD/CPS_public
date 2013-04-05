#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of GpowerPlaq class.

  $Id: g_power_plaq.C,v 1.8 2013-04-05 17:51:14 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:51:14 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_power_plaq/g_power_plaq.C,v 1.8 2013-04-05 17:51:14 chulwoo Exp $
//  $Id: g_power_plaq.C,v 1.8 2013-04-05 17:51:14 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.8 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_power_plaq/g_power_plaq.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// g_power_plaq.C
//
// GpowerPlaq is derived from Lattice and is relevant to the 
// power plaquette action. This action is the same as
// the standard Wilson action with the irrelevant power plaquette
// term added to it. The full action is:
//
// Sum_p [ beta * { -Tr[U_p]/3} + ( {1 - Tr[U_p]/3} / c )^k ]
//
// with c = GJP.PowerPlaqCutoff() and k = GJP.PowerPlaqExponent()
//
// This action supresses plaquettes with {1 - Tr[U_p]/3} > c 
// and threfore reduces lattice dislocations.
//
//------------------------------------------------------------------

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


//------------------------------------------------------------------
//! The number of floating point numbers in an SU(3) matrix.
enum { MATRIX_SIZE = 18 };
//------------------------------------------------------------------

//------------------------------------------------------------------
// DRAM temp buffer used for scu transfer
//------------------------------------------------------------------
static Matrix m_tmp1, m_tmp2;

//------------------------------------------------------------------
// static variables used only inside this file
//------------------------------------------------------------------
static IFloat invs3 = -1./3.;




//  CRAM temp buffer
#ifdef _TARTAN
static Matrix *mp0 = (Matrix *)CRAM_SCRATCH_ADDR;	// ihdot
static Matrix *mp1 = mp0 + 1;
static Matrix *mp2 = mp1 + 1;
static Matrix *mp3 = mp2 + 1;
#else
static Matrix mt0;
static Matrix mt1;
static Matrix mt2;
static Matrix mt3;
static Matrix *mp0 = &mt0;		// ihdot
static Matrix *mp1 = &mt1;
static Matrix *mp2 = &mt2;
static Matrix *mp3 = &mt3;
#endif 


//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
GpowerPlaq::GpowerPlaq()
{
  cname = "GpowerPlaq";
  char *fname = "GpowerPlaq()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
GpowerPlaq::~GpowerPlaq()
{
  char *fname = "~GpowerPlaq()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// GclassType Gclass(void):
// It returns the type of gauge class
//------------------------------------------------------------------
GclassType GpowerPlaq::Gclass(void){
  return G_CLASS_POWER_PLAQ;
}

const unsigned CBUF_MODE2 = 0xcca52112;
const unsigned CBUF_MODE4 = 0xcca52112;


//------------------------------------------------------------------
/*!
  \param force The computed force from the gauge action.
  \param x the lattice site coordinates.
  \param mu The direction mu.
  \todo Could this not be be a virtual Lattice method, if only so that I could
  avoid rewriting this?
*/
//------------------------------------------------------------------
void GpowerPlaq::GforceSite(Matrix& force, int *x, int mu)
{
  char *fname = "GforceSite(M&,i*,i)";
  VRB.Func(cname,fname);

  setCbufCntrlReg(4, CBUF_MODE4);

  Matrix *u_off = GaugeField()+GsiteOffset(x)+mu;


  //----------------------------------------
  //  get power staple
  //     mp1 = power staple
  //----------------------------------------
  PowerStaple(*mp1, x, mu);	
  

  //----------------------------------------
  // mp2 = U_mu(x)
  //----------------------------------------
  moveMem((IFloat *)mp2, (IFloat *)u_off,
  	MATRIX_SIZE * sizeof(IFloat));

  
  //----------------------------------------
  // force = -1/3*U_mu(x)*stap
  // "-" because no staggered phase
  //----------------------------------------
  mDotMEqual((IFloat *)&force, (const IFloat *)mp2, (const IFloat *)mp1);

  Float tmp = invs3;
  vecTimesEquFloat((IFloat *)&force, tmp, MATRIX_SIZE);

  
  // mp1 and mp2 free

  mp1->Dagger((IFloat *)&force);
  force.TrLessAntiHermMatrix(*mp1);
}


//------------------------------------------------------------------
// EvolveMomGforce(Matrix *mom, Float dt):
// It evolves the canonical momentum mom by dt
// using the pure gauge force.
//------------------------------------------------------------------
ForceArg GpowerPlaq::EvolveMomGforce(Matrix *mom, Float dt){
  char *fname = "EvolveMomGforce(M*,F)";
  VRB.Func(cname,fname);
  
  Float L1=0.0;
  Float L2=0.0;
  Float Linf=0.0;

  setCbufCntrlReg(4, CBUF_MODE4);

  int x[4];
  
  for(x[0] = 0; x[0] < GJP.XnodeSites(); ++x[0]) {
    for(x[1] = 0; x[1] < GJP.YnodeSites(); ++x[1]) {
      for(x[2] = 0; x[2] < GJP.ZnodeSites(); ++x[2]) {
	for(x[3] = 0; x[3] < GJP.TnodeSites(); ++x[3]) {
	  
	  int uoff = GsiteOffset(x);
	  
	  for (int mu = 0; mu < 4; ++mu) {
	    GforceSite(*mp0, x, mu);
	    
	    IFloat *ihp = (IFloat *)(mom+uoff+mu);
	    IFloat *dotp = (IFloat *)mp0;
	    fTimesV1PlusV2(ihp, dt, dotp, ihp, 
			   MATRIX_SIZE);
	    Float norm = ((Matrix*)dotp)->norm();
	    Float tmp = sqrt(norm);
	    L1 += tmp;
	    L2 += norm;
	    Linf = (tmp>Linf ? tmp : Linf);
	  }
	}
      }
    }
  }

  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();

  return ForceArg(dt*L1, dt*sqrt(L2), dt*Linf);

}

//------------------------------------------------------------------
// Float GhamiltonNode(void):
// The pure gauge Hamiltonian of the node sublattice.
//------------------------------------------------------------------
Float GpowerPlaq::GhamiltonNode(void){
  char *fname = "GhamiltonNode()";
  VRB.Func(cname,fname);

  Float tmp = GJP.Beta()*invs3;
  Float sum = SumReTrPlaqNode();
  sum *= tmp;
  sum += SumPowerPlaqNode();

  return sum;

}


//------------------------------------------------------------------
/*!
  The power plaquette is
  \f[
  ( (1 - 1/3 Re Tr [U_\mu(x) U_\nu(x+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)])/c )^k
\f]
  where the real parameters \a c and \a k are those in the definition of
  the power plaquette action.

  \param x the coordinates of the lattice site at the start of the plaquette
  \param mu The first plaquette direction.
  \param nu The second plaquette direction; should be different from \a mu.
  \return The computed power plaquette term.  
*/
//------------------------------------------------------------------
Float GpowerPlaq::PowerPlaq(int *x, int mu, int nu) const
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
  At each site \a x and different directions \a mu and \a nu,
  the power plaquette is

\f[
    ( (1 - 1/3 Re Tr [U_\mu(x) U_\nu(x+\nu)
                      U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)])/c )^k
\f]

  where the real parameters \a c and \a k are those in the definition of
  the power plaquette action. This computes the
  sum of the power plaquette over all local lattice sites \a x and all six
  \f$ \mu-\nu \f$ planes. 
 
  \return The locally summed power plaquette term.
*/
//------------------------------------------------------------------
Float GpowerPlaq::SumPowerPlaqNode(void) const
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
  At each site \a x and different directions \a mu and \a nu,
  the power plaquette is
  \f[
  \{ (1 - 1/3 Re Tr [U_\mu(x) U_\nu(x+\nu)
                    U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)])/c  \}^k
  \f]
  where the real parameters \a c and \a k are those in the definition of
  the power plaquette action. This computes the
  sum of the power plaquette over all lattice sites \a x and all six
  \f$ \mu-\nu \f$ planes. 
 
  \return The globally summed power plaquette term.
*/
//------------------------------------------------------------------
Float GpowerPlaq::SumPowerPlaq(void) const
{
  char *fname = "SumPowerPlaq() const";
  VRB.Func(cname,fname);

  Float sum = SumPowerPlaqNode();
  glb_sum(&sum);
  return sum;
}


//------------------------------------------------------------------
/*!
  The staple is:
\f[
  \sum_{\nu \neq \mu} \{                                                 
  ps(x,\mu,v)
    [ U_\nu(x+\mu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x) ]
  + ps(x-\nu,\mu,v)
    [ U^\dagger_\nu(x+\mu-\nu) U^\dagger_\mu(x-\nu) U_\nu(x-\nu) ] \}
\f]
  where
\f[
 ps(x,\mu,\nu) = 
 \beta + k/c ((1 - ReTr[U_p(x,\mu,\nu)]/3) / c )^{k-1}
 \f]

 where the real parameters \f$\beta\f$, \a c and \a k are those in the
 definition of the power plaquette action and \f$ U_p \f$ is the standard
 plaquette.
*/
//------------------------------------------------------------------
void GpowerPlaq::PowerStaple(Matrix& pstap, int *x, int mu)
{
  char *fname = "PowerStaple(M&,i*,i)";
  VRB.Func(cname,fname);

  Float plaq;
  Float tmp;
  Float ps;
  Float inv_cutoff = 1.0 / GJP.PowerPlaqCutoff();
  int exponent = GJP.PowerPlaqExponent();
  Float coeff = Float(exponent) * inv_cutoff;
  Float beta = GJP.Beta();

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
      ps = beta + coeff * ps;

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
      ps = beta + coeff * ps;

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


//------------------------------------------------------------------------------
// void GactionGradient(Matrix &grad, int *x, int mu)
//   Calculates the partial derivative of the gauge action
//   w.r.t. the link U_mu(x).
//------------------------------------------------------------------------------
void GpowerPlaq::GactionGradient(Matrix &grad, int *x, int mu)
{
  ERR.NotImplemented(cname, "GactionGradient(M&,i*,i)") ;
}

void GpowerPlaq::AllStaple(Matrix & stap, const int *x, int mu)
{
  ERR.NotImplemented(cname, "AllStaple(M&,i*,i)") ;
}

CPS_END_NAMESPACE
