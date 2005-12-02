#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of Gwilson class.

  $Id: g_wilson.C,v 1.8 2005-12-02 16:32:51 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2005-12-02 16:32:51 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_wilson/g_wilson.C,v 1.8 2005-12-02 16:32:51 chulwoo Exp $
//  $Id: g_wilson.C,v 1.8 2005-12-02 16:32:51 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.8 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_wilson/g_wilson.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// g_wilson.C
//
// Gwilson is derived from Lattice and is relevant to the 
// standard Wilson single plaquette action.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
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


//------------------------------------------------------------------
enum { MATRIX_SIZE = 18 };
//------------------------------------------------------------------

//------------------------------------------------------------------
// static variables used only inside this file
//------------------------------------------------------------------
static IFloat invs3 = -1./3.;




//  CRAM temp buffer
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


//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
Gwilson::Gwilson()
{
  cname = "Gwilson";
  char *fname = "Gwilson()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
Gwilson::~Gwilson()
{
  char *fname = "~Gwilson()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// GclassType Gclass(void):
// It returns the type of gauge class
//------------------------------------------------------------------
GclassType Gwilson::Gclass(void){
  return G_CLASS_WILSON;
}

const unsigned CBUF_MODE4 = 0xcca52112;


//------------------------------------------------------------------
/*!
  \param force The computed force from the gauge action.
  \param x the lattice site coordinates.
  \param mu The direction mu.
  \todo Could this not be be a virtual Lattice method?
*/
//------------------------------------------------------------------
void Gwilson::GforceSite(Matrix& force, int *x, int mu)
{
  char *fname = "GforceSite(M&,i*,i)";
//  VRB.Func(cname,fname);

  setCbufCntrlReg(4, CBUF_MODE4);

  Matrix *u_off = GaugeField()+GsiteOffset(x)+mu;


  //----------------------------------------
  //  get staple
  //     mp1 = staple
  //----------------------------------------
  Staple(*mp1, x, mu);	
  ForceFlops += 198*3*3+12+216*3;
  

  //----------------------------------------
  // mp2 = U_mu(x)
  //----------------------------------------
  moveMem((IFloat *)mp2, (IFloat *)u_off+BANK4_BASE+BANK_SIZE,
  	MATRIX_SIZE * sizeof(IFloat));

  
  //----------------------------------------
  // force = -beta/3*U_mu(x)*stap
  // "-" because no staggered phase
  //----------------------------------------
  mDotMEqual((IFloat *)&force, (const IFloat *)mp2, (const IFloat *)mp1);

  Float tmp = GJP.Beta()*invs3;
  vecTimesEquFloat((IFloat *)&force, tmp, MATRIX_SIZE);

  
  // mp1 and mp2 free

  mp1->Dagger((IFloat *)&force);
  force.TrLessAntiHermMatrix(*mp1);
  ForceFlops += 198+18+24;
}

#if 0
#define PROFILE
//------------------------------------------------------------------
// EvolveMomGforce(Matrix *mom, Float step_size):
// It evolves the canonical momentum mom by step_size
// using the pure gauge force.
//------------------------------------------------------------------
void Gwilson::EvolveMomGforce(Matrix *mom, Float step_size){
  char *fname = "EvolveMomGforce(M*,F)";
  VRB.Func(cname,fname);
#ifdef PROFILE
  Float time = -dclock();
  ForceFlops=0;
#endif
  
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
	    fTimesV1PlusV2(ihp, step_size, dotp, ihp+BANK4_BASE, 18);
	  }
	}
      }
    }
  }
#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,ForceFlops,time);
#endif
}
#endif

//------------------------------------------------------------------
// Float GhamiltonNode(void):
// The pure gauge Hamiltonian of the node sublattice.
//------------------------------------------------------------------
Float Gwilson::GhamiltonNode(void){
  char *fname = "GhamiltonNode()";
  VRB.Func(cname,fname);

  Float tmp = GJP.Beta()*invs3;
  Float sum = SumReTrPlaqNode();
  sum *= tmp;
  return sum;

}

//------------------------------------------------------------------------------
// void GactionGradient(Matrix &grad, int *x, int mu)
//   Calculates the partial derivative of the gauge action
//   w.r.t. the link U_mu(x).
//------------------------------------------------------------------------------
void Gwilson::GactionGradient(Matrix &grad, int *x, int mu)
{
  char *fname = "GactionGradient(M&,I*,I)" ;
  VRB.Func(cname, fname) ;

  //----------------------------------------------------------------------------
  // get staple
  //----------------------------------------------------------------------------
  Staple(grad, x, mu) ;

  //----------------------------------------------------------------------------
  // grad = - (beta/3) * staple
  //   N.B. invs3 should be -1/3
  //----------------------------------------------------------------------------
  Float tmp = GJP.Beta() * invs3 ;
  vecTimesEquFloat((IFloat *)&grad, tmp, MATRIX_SIZE) ;
}

/*!
  The staple sum around the link \f$ U_\mu(x) \f$ is
\f[
  \sum_{\nu \neq \mu} [
           U_\nu(x+\mu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)                  
        +  U^\dagger_\nu(x+\mu-\nu) U^\dagger_\mu(x-\nu) U_\nu(x-\nu)  ]
\f]
  \param x The coordinates of the lattice site 
  \param mu The link direction
  \param stap The computed staple sum.
*/
void Gwilson::AllStaple(Matrix & stap, const int *x, int mu){
  char * fname = "AllStaple()";
  VRB.Func(cname, fname);
  BufferedStaple(stap, x, mu); 
}

CPS_END_NAMESPACE
