#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_power_rect/g_power_rect.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: g_power_rect.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:37  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:28  anj
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
//  Revision 1.2  2001/05/25 06:16:10  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: g_power_rect.C,v $
//  $Revision: 1.1.1.1 $
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
#include<util/lattice.h>
#include<util/verbose.h>
#include<util/vector.h>
#include<util/gjp.h>
#include<util/gw_hb.h>
#include<comms/nga_reg.h>
#include<comms/glb.h>
#include<comms/scu.h>
#include<comms/cbuf.h>
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
  moveMem((IFloat *)mp2, (IFloat *)u_off+BANK4_BASE+BANK_SIZE,
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
  moveMem((IFloat *)mp2, (IFloat *)u_off+BANK4_BASE+BANK_SIZE,
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
void GpowerRect::EvolveMomGforce(Matrix *mom, Float step_size){
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
// Float PowerPlaq(int *x, int mu, int nu) const:
// It calculates the power plaquette 
// field at site x, mu, nu with mu < nu.
// The power plaquette field is:
//
// pp(x,u,v) = { (1 - ReTr[U_p(x,u,v)]/3) / c }^k
//
// with c = GJP.PowerPlaqCutoff() and
//      k = GJP.PowerPlaqExponent()
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
// Float SumPowerPlaqNode(void) const:
// It calculates the sum of the power plaquette  
// field at each site of the node sublattice.
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
// Float SumPowerPlaq(void) const:
// It calculates the sum of the power plaquette 
// field at each site of the whole lattice
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
// void PowerStaple(Matrix& pstap, int *x, int mu):
// It calculates the staple field at x, mu.
// The staple field is:
//
// V_u(x) = \sum_v(!=u) {
//      ps(x,u,v)   * [ U_v(x+u) U_u(x+v)~ U_v(x)~     ]
//    + ps(x-v,u,v) * [ U_v(x+u-v)~ U_u(x-v)~ U_v(x-v) ] }
//
// where
//
// ps(x,u,v) = 
//      beta*c_0 + {k/c} * { (1 - ReTr[U_p(x,u,v)]/3) / c }^(k-1)
//
// with c = GJP.PowerPlaqCutoff() and
//      k = GJP.PowerPlaqExponent()
//    c_0 = ( 1.0 - 8.0 * GJP.C1() )
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
      mp3->Dagger((IFloat *)p1+BANK4_BASE+BANK_SIZE);

      //----------------------------------------------------------
      // p1 = &U_v(x+u)
      //----------------------------------------------------------
      p1 = GetLinkOld(g_offset, x, mu, nu);

      //----------------------------------------------------------
      // mp2 = U_v(x+u) U_u(x+v)~
      //----------------------------------------------------------
      mDotMEqual((IFloat *)mp2, 
		 (const IFloat *)p1+BANK2_BASE,
		 (const IFloat *)mp3);
      
      //----------------------------------------------------------
      //  mp3 = U_v(x)~
      //----------------------------------------------------------
      mp3->Dagger((IFloat *)(g_offset+nu)+BANK4_BASE);
      
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
		 (const IFloat *)(g_offset+mu)+BANK2_BASE,
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
		(const IFloat *)&m_tmp1+BANK4_BASE, 
		MATRIX_SIZE*sizeof(Float));
      else
	vecAddEquVec((IFloat *)&pstap, 
		     (const IFloat *)&m_tmp1+BANK4_BASE,
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
      moveMem((IFloat *)mp3, (IFloat *)p1+BANK2_BASE, 
	      MATRIX_SIZE * sizeof(IFloat));

      //----------------------------------------------------------
      // mp2 = U_u(x-v)*U_v(x+u-v)
      //----------------------------------------------------------
      mDotMEqual((IFloat *)mp2, 
		 (const IFloat *)(g_offpv+mu)+BANK4_BASE,
		 (const IFloat *)mp3);
      
      //----------------------------------------------------------
      // mp3 = U_v(x+u-v)~*U_u(x-v)~ = mp2~
      //----------------------------------------------------------
      mp3->Dagger((IFloat *)mp2);


      //----------------------------------------------------------
      // mp2 = U_v(x-v)
      //----------------------------------------------------------
      moveMem((IFloat *)mp2, 
	      (const IFloat *)(g_offpv+nu)+BANK2_BASE,
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
		 (const IFloat *)(g_offset+mu)+BANK4_BASE,
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
		   (IFloat *)&m_tmp1+BANK2_BASE,
		   MATRIX_SIZE);

    }
  }
}


//------------------------------------------------------------------
// It calculates the rectangle staple field at x, mu.
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
    mp4->Dagger((IFloat *)GetLink(link_site, nu)+BANK4_BASE) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+v)~
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    mp3->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE+BANK_SIZE) ;

    //----------------------------------------------------------
    // mp2 = U_u(x+v)~ U_v(x)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x+u+v)~
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE);

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
    mDotMEqual((IFloat *)mp4, (const IFloat *)p1+BANK4_BASE+BANK_SIZE,
               (const IFloat *)mp3);

    //----------------------------------------------------------
    // mp2 = U_u(x+u)
    //----------------------------------------------------------
    --(link_site[mu]) ;
    moveMem((IFloat *)mp2, (const IFloat *)GetLink(link_site, mu)+BANK4_BASE,
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
	       (const IFloat *)(g_offset+mu)+BANK4_BASE+BANK_SIZE,
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
		 (const IFloat *)&m_tmp1+BANK4_BASE,
		 MATRIX_SIZE);

//PMV

    //----------------------------------------------------------
    // mp4 = U_v(x+2u-v)~
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    --(link_site[nu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, nu)+BANK4_BASE+BANK_SIZE) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+u) U_v(x+2u-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x+u-v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE) ;

    //----------------------------------------------------------
    // mp2 = U_u(x+u) U_v(x+2u-v)~ U_u(x+u-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x-v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE+BANK_SIZE) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+u) U_v(x+2u-v)~ U_u(x+u-v)~ U_u(x-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp2 = U_v(x-v)
    //----------------------------------------------------------
    moveMem((IFloat *)mp2, (const IFloat *)GetLink(link_site, nu)+BANK4_BASE,
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
	       (const IFloat *)(g_offset+mu)+BANK4_BASE+BANK_SIZE,
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
		 (const IFloat *)&m_tmp1+BANK4_BASE,
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
    mDotMEqual((IFloat *)mp3, (const IFloat *)p1+BANK4_BASE+BANK_SIZE,
               (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // mp4 = U_u(x-2v)~
    //----------------------------------------------------------
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE) ;

    //----------------------------------------------------------
    // mp2 = U_u(x-2v)~ U_v(x-2v) U_v(x-v)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp4, (const IFloat *)mp3) ;

    //----------------------------------------------------------
    // mp4 = U_v(x+u-2v)~
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, nu)+BANK4_BASE+BANK_SIZE) ;

    //----------------------------------------------------------
    // mp3 = U_v(x+u-2v)~ U_u(x-2v)~ U_v(x-2v) U_v(x-v)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp4, (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // mp2 = U_v(x+u-v)~
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    mp2->Dagger((IFloat *)GetLink(link_site, nu)+BANK4_BASE) ;

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
	       (const IFloat *)(g_offset+mu)+BANK4_BASE+BANK_SIZE,
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
		 (const IFloat *)&m_tmp1+BANK4_BASE,
		 MATRIX_SIZE);

//PMV

    //----------------------------------------------------------
    // mp4 = U_u(x-v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE+BANK_SIZE) ;

    //----------------------------------------------------------
    // mp3 =  U_v(x+u-v)~ U_u(x-v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_u(x-u-v)~
    //----------------------------------------------------------
    --(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE) ;

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
            (const IFloat *)GetLink(link_site, nu)+BANK4_BASE+BANK_SIZE,
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
    moveMem((IFloat *)mp2, (const IFloat *)GetLink(link_site, mu)+BANK4_BASE,
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
	       (const IFloat *)(g_offset+mu)+BANK4_BASE+BANK_SIZE,
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
		 (const IFloat *)&m_tmp1+BANK4_BASE,
		 MATRIX_SIZE);

//PMV

    //----------------------------------------------------------
    // mp4 = U_v(x-u)~
    //----------------------------------------------------------
    mp4->Dagger((IFloat *)GetLink(link_site, nu)+BANK4_BASE+BANK_SIZE) ;

    //----------------------------------------------------------
    // mp3 = U_v(x-u)~ U_u(x-u)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp4, (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // mp4 = U_u(x-u+v)~
    //----------------------------------------------------------
    ++(link_site[nu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE) ;

    //----------------------------------------------------------
    // mp2 = U_u(x-u+v)~ U_v(x-u)~ U_u(x-u)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp4, (const IFloat *)mp3) ;

    //----------------------------------------------------------
    // mp4 = U_u(x+v)~
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE+BANK_SIZE) ;

    //----------------------------------------------------------
    // mp3 = U_u(x+v)~ U_u(x-u+v)~ U_v(x-u)~ U_u(x-u)
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp4, (const IFloat *)mp2) ;

    //----------------------------------------------------------
    // mp2 = U_v(x+u)
    //----------------------------------------------------------
    ++(link_site[mu]) ;
    --(link_site[nu]) ;
    moveMem((IFloat *)mp2, (const IFloat *)GetLink(link_site, nu)+BANK4_BASE,
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
	       (const IFloat *)(g_offset+mu)+BANK4_BASE+BANK_SIZE,
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
		 (const IFloat *)&m_tmp1+BANK4_BASE,
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
            (const IFloat *)GetLink(link_site, nu)+BANK4_BASE+BANK_SIZE,
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
    mp4->Dagger((IFloat *)GetLink(link_site, mu)+BANK4_BASE) ;

    //----------------------------------------------------------
    // mp2 = U_v(x+u) U_v(x+u+v) U_u(x+2v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp2, (const IFloat *)mp3, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp4 = U_v(x+v)~
    //----------------------------------------------------------
    --(link_site[nu]) ;
    mp4->Dagger((IFloat *)GetLink(link_site, nu)+BANK4_BASE+BANK_SIZE) ;

    //----------------------------------------------------------
    // mp3 = U_v(x+u) U_v(x+u+v) U_u(x+2v)~ U_v(x+v)~
    //----------------------------------------------------------
    mDotMEqual((IFloat *)mp3, (const IFloat *)mp2, (const IFloat *)mp4) ;

    //----------------------------------------------------------
    // mp2 = U_v(x)~
    //----------------------------------------------------------
    --(link_site[nu]) ;
    mp2->Dagger((IFloat *)GetLink(link_site, nu)+BANK4_BASE) ;


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
	       (const IFloat *)(g_offset+mu)+BANK4_BASE+BANK_SIZE,
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
		 (const IFloat *)&m_tmp1+BANK4_BASE,
		 MATRIX_SIZE);

//PMV

    //----------------------------------------------------------
    // dummy read to switch CBUF banks for looping
    //----------------------------------------------------------
    *((IFloat *)mp4) = *((IFloat *)p1+BANK4_BASE+BANK_SIZE) ;
  }

//VRB.FuncEnd(cname, fname) ;
}


//------------------------------------------------------------------
// It calculates the power rectangle 
// field at site x, mu, nu with mu < nu.
// The power plaquette rectangle is:
//
// pp(x,u,v) = { (1 - ReTr[U_r(x,u,v)]/3) / c }^k
//
// with c = GJP.PowerRectCutoff() and
//      k = GJP.PowerRectExponent()
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
// It calculates the sum of the power rectangle  
// field at each site of the node sublattice.
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
// It calculates the sum of the power rectangle 
// field at each site of the whole lattice
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
