#include<config.h>
CPS_START_NAMESPACE
#ifndef INCLUDED_W_BARYON
#define INCLUDED_W_BARYON

CPS_END_NAMESPACE
#include <util/data_types.h>        // Float, Complex
#include <alg/w_spect_arg.h>        // SourceKind
#include <alg/w_ginfo.h>
#include <alg/w_momenta.h>
CPS_START_NAMESPACE

//---------------------------------------------------------------------------
// Forward Declarations  -- classes defined in other translation units
//---------------------------------------------------------------------------
class Lattice;                      // defined in util/include/lattice.h
class CgArg;                        // defined in  alg/include/cg_arg.h
class CommonArg;                    // defined in  alg/include/common_arg.h
class AlgWspect;                    // defined in  alg/include/alg_w_spect.h
class WspectGinfo;                      
class WspectHyperRectangle;                     
class WspectQuark;                     

//---------------------------------------------------------------------------
// class WspectBaryon 
//---------------------------------------------------------------------------
class WspectBaryon : public WspectGinfo {

public:

  //-------------------------------------------------------------------------
  // WspectBaryon::Constitutes
  //-------------------------------------------------------------------------
  enum Constitutes { ABC,       // three flavours [could be degenerate]
		     AAA,       // one   flavour
		     AAC,       // two   flavours [could be degenerate]
		     DELTA_CONSTI   = AAA,
		     NUCLEON_CONSTI = AAC
  };


  //-------------------------------------------------------------------------
  // WspectBaryon::TwoDiracMats
  //-------------------------------------------------------------------------
  enum TwoDiracMats { Gamma5Gamma5,   // B = q1 [q2 C Gamma5 q3]
		      GammaXGammaX,   // B = q1 [q2 C GammaX q3]
		      GammaYGammaY,   // B = q1 [q2 C GammaY q3]
		      GammaZGammaZ,   // B = q1 [q2 C GammaZ q3]
		      GammaTGammaT,   // B = q1 [q2 C GammaT q3]
                      UnitUnit,       // B = q1 [q2 C Unit   q3]
		      DELTAX_DIRAC  = GammaXGammaX,  
		      DELTAY_DIRAC  = GammaYGammaY,
		      DELTAZ_DIRAC  = GammaZGammaZ,
		      DELTAT_DIRAC  = GammaTGammaT,
		      NUCLEON_DIRAC = Gamma5Gamma5 
  };
  
  // CTOR
  //-------------------------------------------------------------------------
  WspectBaryon(const WspectQuark & Q1,        // Baryon = q1 [q2 GAMMA q3]
	       const WspectQuark & Q2,
	       const WspectQuark & Q3,
	       const WspectHyperRectangle & whr,
	       WspectBaryon::Constitutes con,
	       WspectBaryon::TwoDiracMats dmats);

  // CTOR
  //-------------------------------------------------------------------------
  WspectBaryon(const IFloat * Q1, 
	       const IFloat * Q2,
	       const IFloat * Q3,
	       const WspectHyperRectangle & whr,
	       WspectBaryon::Constitutes con,
	       WspectBaryon::TwoDiracMats dmats);


  // DTOR
  //-------------------------------------------------------------------------
  ~WspectBaryon();


  // ACCESSORS
  //-------------------------------------------------------------------------
  void print(char *filename, WbaryonFold = BARYON_FOLD) const;  


  // PRIVATEs
  //-------------------------------------------------------------------------
private:
  // static data member
  static char *                d_class_name;

  // the correlators after Dirac projection with Tr(Unit...) etc
  Float *                      d_unit_p;      // Tr (Unit ...)
  Float *                      d_gamma_p;     // Tr (Gamma ...)
  Float *                      d_up_p;        // Tr (1+Gamma...)
  Float *                      d_down_p;      // Tr (1-Gamma...)
  Float *                      d_buffer_p;    // holds dynamical memory   
  int                          d_buffer_size;  

  // data members to specify the properties of this baryon
  Constitutes                  d_flavors;     // the quark constitutes
  TwoDiracMats                 d_diracs;      // the Dirac matrices

  const IFloat *                d_quark1_p;    // only reference
  const IFloat *                d_quark2_p;    // only reference
  const IFloat *                d_quark3_p;    // only reference
  const WspectHyperRectangle & d_whr;         // only reference

  // The following data members are here for computation efficiency.
  int                          d_prop_dir;  
  int                          d_lclMin[LORENTZs];
  int                          d_lclMax[LORENTZs];


  // The following two member functions are supporting DiracAlgebra(..)

  // Note:
  //    This is the ONLY term if q1, q2 and q3 are distinguishable, 
  //      OR ONE of the terms if q1, q2 (and possible q3) are in-dist...
  // Explictly:
  //    Q1(Dx, Dy) Q2(alpha, beta) Q3(rho, phi) * Gamma(alpha,rho) 
  //                                            * Gamma(beta, phi)
  void DiracAlgebraABC(int Dx, int Dy,
		       int a_local_site_offset_for_quark_prop,
		       Complex &result) ;    // result += ...

  // Note:
  //   This is the ADDITIONAL term if q1, q2 (and perhaps q3) 
  //   are indistinguishable
  // Explictly:
  //    Q1(Dx, beta) Q2(alpha, Dy) Q3(rho, phi) * Gamma(alpha,rho) 
  //                                            * Gamma(beta, phi)
  void DiracAlgebraAAC(int Dx, int Dy,
		       int a_local_site_offset_for_quark_prop,
		       Complex &result) ;    // result += ...


  //-------------------------------------------------------------------------
  // CALCULATION PROCEDURE -- BARYONs
  //   step 1: color algebra 
  //   step 2: Dirac Algebra
  //   step 3: Momentum Projection
  //   step 4: Dirac Projection
  //   finish: Everything()
  //-------------------------------------------------------------------------

  // step 1
  void ColorAlgebra(int D1x, int D2x, int D3x, 
		    int D1y, int D2y, int D3y, 
		    int a_local_site_offset_for_quark_prop,
		    Complex &result);           // result += ...

  // step 2
  void DiracAlgebra(int Dx, int Dy,
		    int a_local_site_offset_for_quark_prop,
		    Complex &result);           // result += ...
  
  // step 3
  Complex ZeroMomProject(int Dx, int Dy,
			 int a_local_time_slice);

  // step 4
  Float DiracProjectUnit(int a_local_time_slice) ;  
  Float DiracProjectGamma(int a_local_time_slice) ;  

  // finish
  void  Everything();  
}; 

#endif // ! _INCLUDED_W_BARYON



CPS_END_NAMESPACE
