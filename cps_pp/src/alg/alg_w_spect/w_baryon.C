#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:40 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_w_spect/w_baryon.C,v 1.9 2004-08-18 11:57:40 zs Exp $
//  $Id: w_baryon.C,v 1.9 2004-08-18 11:57:40 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: w_baryon.C,v $
//  $Revision: 1.9 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_w_spect/w_baryon.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <alg/w_all.h>
#include <util/qcdio.h>
#include <util/error.h>
#include <util/verbose.h>
#include <util/lattice.h>
#include <util/vector.h>
#include <comms/glb.h>                 // glb_sum
#include <alg/alg_w_spect.h>         // AlgWspect::GetCounter()
CPS_START_NAMESPACE

//---------------------------------------------------------------------------
// For the purpose of debugging during code upgrade
//---------------------------------------------------------------------------
//#define DEBUG_W_BARYON

#ifdef  DEBUG_W_BARYON
//  #define DEBUG_W_BARYON_PROJECT
//  #define DEBUG_W_BARYON_PROJECT_IMAG
//  #define DEBUG_W_BARYON_DIRAC
//  #define DEBUG_W_BARYON_COLOR
CPS_END_NAMESPACE
#include <util/gjp.h>
CPS_START_NAMESPACE
#endif

//---------------------------------------------------------------------------
// static data members
//---------------------------------------------------------------------------
char * WspectBaryon::d_class_name = "WspectBaryon";


//---------------------------------------------------------------------------
// WspectBaryon::WspectBaryon()
//---------------------------------------------------------------------------
WspectBaryon::WspectBaryon(const WspectQuark &q1, 
			   const WspectQuark &q2, 
			   const WspectQuark &q3,
			   const WspectHyperRectangle & whr,
			   WspectBaryon::Constitutes con,
			   WspectBaryon::TwoDiracMats dmats)
  : d_quark1_p(q1.Data()),
    d_quark2_p(q2.Data()),
    d_quark3_p(q3.Data()),
    d_flavors(con),
    d_whr(whr),
    d_diracs(dmats)
{
  Everything();
}



//---------------------------------------------------------------------------
// WspectBaryon::WspectBaryon()
//---------------------------------------------------------------------------
WspectBaryon::WspectBaryon(const IFloat *q1, 
			   const IFloat *q2, 
			   const IFloat *q3,
			   const WspectHyperRectangle & whr,
			   WspectBaryon::Constitutes con,
			   WspectBaryon::TwoDiracMats dmats) 
  : d_quark1_p(q1),
    d_quark2_p(q2),
    d_flavors(con),
    d_quark3_p(q3),
    d_whr(whr),
    d_diracs(dmats)
{
  Everything();
}


//---------------------------------------------------------------------------
// void WspectBaryon::Everything()
//---------------------------------------------------------------------------
void
WspectBaryon::Everything() 
{
  VRB.Func(d_class_name, ctor_str);

  // consistency check -- both quark prop's are from the same source
  //-------------------------------------------------------------------------
  // Ping -- later on, we may want to check for the consistency 
  //         of WspectHyperRectangle of each quark propagator. ???
  //   ERR.General(d_class_name, ctor_str, "quarks w/ diff dir or src\n");



  // calculate d_lclMin and d_lclMax, whose [prop_dir]th elements are 
  // meant to be modified.
  //-------------------------------------------------------------------------
  {
    d_prop_dir = d_whr.dir();  
    const int *low  = d_whr.lclMin();
    const int *high = d_whr.lclMax();
    for (int i = 0; i < LORENTZs; ++i) {
      d_lclMin[i] = low[i];
      d_lclMax[i] = high[i];
    }
  }
  
  // number of sites along the propagation direction
  //-------------------------------------------------------------------------
  int glb_walls = glb_sites[d_prop_dir];

  // dynamically allocate space and zero them out.
  //-------------------------------------------------------------------------
  {
    d_buffer_size = glb_walls * 4;
    d_buffer_p = (Float *)smalloc(d_buffer_size*sizeof(Float));

    if (!d_buffer_p) 
      ERR.Pointer(d_class_name, ctor_str, empty_str);
    VRB.Smalloc(d_class_name,ctor_str, empty_str, 
		d_buffer_p, d_buffer_size*sizeof(IFloat));

    for (int i = 0; i < d_buffer_size; ++i) {
      *d_buffer_p++ = 0.0;      
    }
    d_buffer_p -= d_buffer_size;
    
    d_unit_p  = d_buffer_p;
    d_gamma_p = d_unit_p  + glb_walls;    
    d_up_p    = d_gamma_p + glb_walls;    
    d_down_p  = d_up_p    + glb_walls;    
  }
  

  // Take care of local parts of data
  //-------------------------------------------------------------------------
  // Note: each node takes care of itself simultaneously.
  // 
  {
    int lcl_walls      = lcl_sites[d_whr.dir()];
    Float *lcl_unit_p  = d_unit_p  + lcl2glb_offset[d_whr.dir()];
    Float *lcl_gamma_p = d_gamma_p + lcl2glb_offset[d_whr.dir()];
    for (int lclW = 0; lclW < lcl_walls; ++lclW) {
      *lcl_unit_p++  = DiracProjectUnit(lclW);
      *lcl_gamma_p++ = DiracProjectGamma(lclW);
    }
  }
  
  // Global sum over all walls
  //-------------------------------------------------------------------------
  {
    for (int glbW = 0; glbW < glb_walls; ++glbW) {
      glb_sum((Float *)(d_unit_p + glbW));	
      glb_sum((Float *)(d_gamma_p + glbW));
      d_up_p  [glbW] = d_unit_p[glbW] + d_gamma_p[glbW];
      d_down_p[glbW] = d_unit_p[glbW] - d_gamma_p[glbW];
    }
  }


  
#ifdef DEBUG_W_BARYON_DIRAC
  {
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
    FILE *fp = Fopen("baryon.dirac.dat", "a");
    Complex answer;
    int lcl[LORENTZs];  
    for (lcl[0] = 0; lcl[0] < GJP.XnodeSites(); lcl[0]++) {
      for (lcl[1] = 0; lcl[1] < GJP.YnodeSites(); lcl[1]++) {
	for (lcl[2] = 0; lcl[2] < GJP.YnodeSites(); lcl[2]++) {
	  for (lcl[3] = 0; lcl[3] < GJP.YnodeSites(); lcl[3]++) {
	    for (int Dx = 0; Dx < DIRACs; ++Dx) {
	      for (int Dy = 0; Dy < DIRACs; ++Dy) {
		answer = DiracAlgebra(Dx, Dy, lcl);
		Fprintf(fp, 
			"site[%d %d %d %d] spin[%d%d]: [%g %g]\n",
			lcl[0], lcl[1], lcl[2], lcl[3], Dx, Dy,
			answer.real(), answer.imag());
	      }
	    }
	  }
	}
      }
    }
    Fclose(fp);
  }
#endif  // #ifdef DEBUG_W_BARYON_DIRAC

#ifdef DEBUG_W_BARYON_COLOR
  {
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
    FILE *fp = Fopen("baryon.color.dat", "a");
    Complex answer;
    int lcl[LORENTZs];  
    for (lcl[0] = 0; lcl[0] < GJP.XnodeSites(); lcl[0]++) {
      for (lcl[1] = 0; lcl[1] < GJP.YnodeSites(); lcl[1]++) {
	for (lcl[2] = 0; lcl[2] < GJP.YnodeSites(); lcl[2]++) {
	  for (lcl[3] = 0; lcl[3] < GJP.YnodeSites(); lcl[3]++) {
	    for (int D1x = 0; D1x < DIRACs; ++D1x) {
	      for (int D1y = 0; D1y < DIRACs; ++D1y) {
		for (int D2x = 0; D2x < DIRACs; ++D2x) {
		  for (int D2y = 0; D2y < DIRACs; ++D2y) {
		    for (int D3x = 0; D3x < DIRACs; ++D3x) {
		      for (int D3y = 0; D3y < DIRACs; ++D3y) {
			answer = ColorAlgebra(D1x, D2x, D3x,
					      D1y, D2y, D3y, lcl);
			Fprintf(fp, 
				"site[%d %d %d %d] spin[%d%d][%d%d][%d%d]: [%g %g]\n",
				lcl[0], lcl[1], lcl[2], lcl[3], 
				D1x, D1y, D2x, D2y, D3x, D3y,
				answer.real(), answer.imag());
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    Fclose(fp);
  }

#endif   //  #ifdef DEBUG_W_BARYON_COLOR


}



//---------------------------------------------------------------------------
// WspectBaryon::~WspectBaryon()
//---------------------------------------------------------------------------
WspectBaryon::~WspectBaryon()
{
  VRB.Func(d_class_name, dtor_str);
  VRB.Sfree(d_class_name, dtor_str, empty_str, d_buffer_p);
  sfree(d_buffer_p);
}



//---------------------------------------------------------------------------
// void
// WspectBaryon::ColorAlgebra(int...)
//--------------------------------------------------------------------------- 
// Two color tensors involved: 3! x 3! = 36 terms
//--------------------------------------------------------------------------- 
void
WspectBaryon::ColorAlgebra(int D1x, int D2x, int D3x,
			   int D1y, int D2y, int D3y,
			   int local_site_offset, Complex &answer)
{
  int off_Cy = WspectQuark::weightSrcColor();
  int off_Dy = WspectQuark::weightSrcDirac();
  int off_Cy2 = off_Cy / COMPLEXs;  

  // pointer to quark prop data for the specified local site and dirac index 
  const Complex *q1 = (const Complex *)(d_quark1_p + off_Dy * D1y + 
					local_site_offset + 
					(COMPLEXs*COLORs) * D1x);
  const Complex *q2 = (const Complex *)(d_quark2_p + off_Dy * D2y + 
					local_site_offset + 
					(COMPLEXs*COLORs) * D2x); 
  const Complex *q3 = (const Complex *)(d_quark3_p + off_Dy * D3y + 
					local_site_offset + 
					(COMPLEXs*COLORs) * D3x);

  // 6 color indexes:    sink x (q1, q2, q3) ; source y (q1, q2, q3)
  answer += q1[0        ] * q2[1+off_Cy2] * q3[2+off_Cy ]; //012012
  answer += q1[  off_Cy ] * q2[1        ] * q3[2+off_Cy2]; //012201
  answer += q1[  off_Cy2] * q2[1+off_Cy ] * q3[2        ]; //012120

  answer -= q1[  off_Cy ] * q2[1+off_Cy2] * q3[2        ]; //012210
  answer -= q1[0        ] * q2[1+off_Cy ] * q3[2+off_Cy2]; //012021
  answer -= q1[  off_Cy2] * q2[1        ] * q3[2+off_Cy ]; //012102

  answer += q1[2        ] * q2[  off_Cy2] * q3[1+off_Cy ]; //201012
  answer += q1[2+off_Cy ] * q2[0        ] * q3[1+off_Cy2]; //201201
  answer += q1[2+off_Cy2] * q2[  off_Cy ] * q3[1        ]; //201120

  answer -= q1[2+off_Cy ] * q2[  off_Cy2] * q3[1        ]; //201210
  answer -= q1[2        ] * q2[  off_Cy ] * q3[1+off_Cy2]; //201021
  answer -= q1[2+off_Cy2] * q2[0        ] * q3[1+off_Cy ]; //201102

  answer += q1[1        ] * q2[2+off_Cy2] * q3[  off_Cy ]; //120012
  answer += q1[1+off_Cy ] * q2[2        ] * q3[  off_Cy2]; //120201
  answer += q1[1+off_Cy2] * q2[2+off_Cy ] * q3[0        ]; //120120

  answer -= q1[1+off_Cy ] * q2[2+off_Cy2] * q3[0        ]; //120210
  answer -= q1[1        ] * q2[2+off_Cy ] * q3[  off_Cy2]; //120021
  answer -= q1[1+off_Cy2] * q2[2        ] * q3[  off_Cy ]; //120102

  answer -= q1[2        ] * q2[1+off_Cy2] * q3[  off_Cy ]; //210012
  answer -= q1[2+off_Cy ] * q2[1        ] * q3[  off_Cy2]; //210201
  answer -= q1[2+off_Cy2] * q2[1+off_Cy ] * q3[0        ]; //210120

  answer += q1[2+off_Cy ] * q2[1+off_Cy2] * q3[0        ]; //210210
  answer += q1[2        ] * q2[1+off_Cy ] * q3[  off_Cy2]; //210021
  answer += q1[2+off_Cy2] * q2[1        ] * q3[  off_Cy ]; //210102

  answer -= q1[0        ] * q2[2+off_Cy2] * q3[1+off_Cy ]; //021012
  answer -= q1[  off_Cy ] * q2[2        ] * q3[1+off_Cy2]; //021201
  answer -= q1[  off_Cy2] * q2[2+off_Cy ] * q3[1        ]; //021120

  answer += q1[  off_Cy ] * q2[2+off_Cy2] * q3[1        ]; //021210
  answer += q1[0        ] * q2[2+off_Cy ] * q3[1+off_Cy2]; //021021
  answer += q1[  off_Cy2] * q2[2        ] * q3[1+off_Cy ]; //021102

  answer -= q1[1        ] * q2[  off_Cy2] * q3[2+off_Cy ]; //102012
  answer -= q1[1+off_Cy ] * q2[0        ] * q3[2+off_Cy2]; //102201
  answer -= q1[1+off_Cy2] * q2[  off_Cy ] * q3[2        ]; //102120

  answer += q1[1+off_Cy ] * q2[  off_Cy2] * q3[2        ]; //102210
  answer += q1[1        ] * q2[  off_Cy ] * q3[2+off_Cy2]; //102021
  answer += q1[1+off_Cy2] * q2[0        ] * q3[2+off_Cy ]; //102102

}



//---------------------------------------------------------------------------
// void
// WspectBaryon::DiracAlgebraABC(int...)
//--------------------------------------------------------------------------- 
// Note:
//    An additional minus should be here for the mu = prop_dir 
//    case, however, since we are allowing the prop_dir to be 
//    different than t, it is better for us to taken it into account
//    in the end.
//
// Conventions of Gamma Matrices:
// 
//      gamma(X)      gamma(Y)      gamma(Z)      gamma(T)      gamma(5)
//    0  0  0  i    0  0  0 -1    0  0  i  0    0  0  1  0    1  0  0  0 
//    0  0  i  0    0  0  1  0    0  0  0 -i    0  0  0  1    0  1  0  0
//    0 -i  0  0    0  1  0  0   -i  0  0  0    1  0  0  0    0  0 -1  0
//   -i  0  0  0   -1  0  0  0    0  i  0  0    0  1  0  0    0  0  0 -1
//
//    C = gamma(Y) * gamma(T) = 0  1  0  0
//                             -1  0  0  0
//                              0  0  0 -1
//                              0  0  1  0
//
// Form of Gamma Matrices:                   Calculation here:
// 
//    C GammaX =  0  0  i  0                 (C GammaX)...(C GammaX)^(c.c.)
//                0  0  0 -i
//                i  0  0  0
//                0 -i  0  0
//
//    C GammaY =  0  0  1  0                 (C GammaY)...(C GammaY)^(c.c.)
//                0  0  0  1
//                1  0  0  0
//                0  1  0  0
//
//    C GammaZ =  0  0  0 -i                 (C GammaZ)...(C GammaZ)^(c.c.)
//                0  0 -i  0
//                0 -i  0  0
//               -i  0  0  0
//   
//    C GammaT =  0  0  0  1                 (C GammaT)...(C GammaT)^(c.c.)
//                0  0 -1  0
//                0 -1  0  0
//                1  0  0  0
//
//    C Gamma5 =  0  1  0  0                 (C Gamma5)...(C Gamma5)
//               -1  0  0  0
//                0  0  0  1
//                0  0 -1  1
// 
//--------------------------------------------------------------------------- 
void
WspectBaryon::DiracAlgebraABC(int Dx, int Dy, int local_site_offset, 
			      Complex &answer) 
{
  Complex subtract(0.0, 0.0);
  
  switch (d_diracs) {
  case GammaXGammaX:       

    ColorAlgebra(Dx, 0, 2, Dy, 0, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 3, Dy, 0, 2, local_site_offset, subtract);
    ColorAlgebra(Dx, 2, 0, Dy, 0, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 1, Dy, 0, 2, local_site_offset, subtract);
  
    ColorAlgebra(Dx, 0, 2, Dy, 1, 3, local_site_offset, subtract);
    ColorAlgebra(Dx, 1, 3, Dy, 1, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 0, Dy, 1, 3, local_site_offset, subtract);
    ColorAlgebra(Dx, 3, 1, Dy, 1, 3, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 2, Dy, 2, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 3, Dy, 2, 0, local_site_offset, subtract);
    ColorAlgebra(Dx, 2, 0, Dy, 2, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 1, Dy, 2, 0, local_site_offset, subtract);
  
    ColorAlgebra(Dx, 0, 2, Dy, 3, 1, local_site_offset, subtract);
    ColorAlgebra(Dx, 1, 3, Dy, 3, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 0, Dy, 3, 1, local_site_offset, subtract);
    ColorAlgebra(Dx, 3, 1, Dy, 3, 1, local_site_offset, answer);

    break;

  case GammaYGammaY:

    ColorAlgebra(Dx, 0, 2, Dy, 0, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 3, Dy, 0, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 0, Dy, 0, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 1, Dy, 0, 2, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 2, Dy, 1, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 3, Dy, 1, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 0, Dy, 1, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 1, Dy, 1, 3, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 2, Dy, 2, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 3, Dy, 2, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 0, Dy, 2, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 1, Dy, 2, 0, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 2, Dy, 3, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 3, Dy, 3, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 0, Dy, 3, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 1, Dy, 3, 1, local_site_offset, answer);
    
    break;
    
  case GammaZGammaZ:  

    ColorAlgebra(Dx, 0, 3, Dy, 0, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 2, Dy, 0, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 1, Dy, 0, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 0, Dy, 0, 3, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 3, Dy, 1, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 2, Dy, 1, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 1, Dy, 1, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 0, Dy, 1, 2, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 3, Dy, 2, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 2, Dy, 2, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 1, Dy, 2, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 0, Dy, 2, 1, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 3, Dy, 3, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 2, Dy, 3, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 1, Dy, 3, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 0, Dy, 3, 0, local_site_offset, answer);

    break;

  case GammaTGammaT:

    ColorAlgebra(Dx, 0, 3, Dy, 0, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 2, Dy, 0, 3, local_site_offset, subtract);
    ColorAlgebra(Dx, 2, 1, Dy, 0, 3, local_site_offset, subtract);
    ColorAlgebra(Dx, 3, 0, Dy, 0, 3, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 3, Dy, 1, 2, local_site_offset, subtract);
    ColorAlgebra(Dx, 1, 2, Dy, 1, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 1, Dy, 1, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 0, Dy, 1, 2, local_site_offset, subtract);
  
    ColorAlgebra(Dx, 0, 3, Dy, 2, 1, local_site_offset, subtract);
    ColorAlgebra(Dx, 1, 2, Dy, 2, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 1, Dy, 2, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 0, Dy, 2, 1, local_site_offset, subtract);
  
    ColorAlgebra(Dx, 0, 3, Dy, 3, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 2, Dy, 3, 0, local_site_offset, subtract);
    ColorAlgebra(Dx, 2, 1, Dy, 3, 0, local_site_offset, subtract);
    ColorAlgebra(Dx, 3, 0, Dy, 3, 0, local_site_offset, answer);

    break;    

  case Gamma5Gamma5:

    ColorAlgebra(Dx, 0, 1, Dy, 0, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 0, Dy, 0, 1, local_site_offset, subtract);
    ColorAlgebra(Dx, 2, 3, Dy, 0, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 2, Dy, 0, 1, local_site_offset, subtract);
  
    ColorAlgebra(Dx, 0, 1, Dy, 1, 0, local_site_offset, subtract);
    ColorAlgebra(Dx, 1, 0, Dy, 1, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 3, Dy, 1, 0, local_site_offset, subtract);
    ColorAlgebra(Dx, 3, 2, Dy, 1, 0, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 1, Dy, 2, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 0, Dy, 2, 3, local_site_offset, subtract);
    ColorAlgebra(Dx, 2, 3, Dy, 2, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 2, Dy, 2, 3, local_site_offset, subtract);
  
    ColorAlgebra(Dx, 0, 1, Dy, 3, 2, local_site_offset, subtract);
    ColorAlgebra(Dx, 1, 0, Dy, 3, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 3, Dy, 3, 2, local_site_offset, subtract);
    ColorAlgebra(Dx, 3, 2, Dy, 3, 2, local_site_offset, answer);
    
    break;    

  case UnitUnit:

    ColorAlgebra(Dx, 0, 1, Dy, 0, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 0, Dy, 0, 1, local_site_offset, subtract);
    ColorAlgebra(Dx, 2, 3, Dy, 0, 1, local_site_offset, subtract);
    ColorAlgebra(Dx, 3, 2, Dy, 0, 1, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 1, Dy, 1, 0, local_site_offset, subtract);
    ColorAlgebra(Dx, 1, 0, Dy, 1, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 3, Dy, 1, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 2, Dy, 1, 0, local_site_offset, subtract);
  
    ColorAlgebra(Dx, 0, 1, Dy, 2, 3, local_site_offset, subtract);
    ColorAlgebra(Dx, 1, 0, Dy, 2, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 3, Dy, 2, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 2, Dy, 2, 3, local_site_offset, subtract);
  
    ColorAlgebra(Dx, 0, 1, Dy, 3, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 0, Dy, 3, 2, local_site_offset, subtract);
    ColorAlgebra(Dx, 2, 3, Dy, 3, 2, local_site_offset, subtract);
    ColorAlgebra(Dx, 3, 2, Dy, 3, 2, local_site_offset, answer);
    
    break;    
  }

  answer -= subtract;  
}



//---------------------------------------------------------------------------
// void
// WspectBaryonW::DiracAlgebraAAC(int...)
//--------------------------------------------------------------------------- 
void
WspectBaryon::DiracAlgebraAAC(int Dx, int Dy, int local_site_offset,
			      Complex &answer)
{

  Complex subtract(0.0, 0.0);
  
  switch (d_diracs) {

  case GammaXGammaX:

    ColorAlgebra(Dx, 0, 2, 0, Dy, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 3, 0, Dy, 2, local_site_offset, subtract);
    ColorAlgebra(Dx, 2, 0, 0, Dy, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 1, 0, Dy, 2, local_site_offset, subtract);
  
    ColorAlgebra(Dx, 0, 2, 1, Dy, 3, local_site_offset, subtract);
    ColorAlgebra(Dx, 1, 3, 1, Dy, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 0, 1, Dy, 3, local_site_offset, subtract);
    ColorAlgebra(Dx, 3, 1, 1, Dy, 3, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 2, 2, Dy, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 3, 2, Dy, 0, local_site_offset, subtract);
    ColorAlgebra(Dx, 2, 0, 2, Dy, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 1, 2, Dy, 0, local_site_offset, subtract);
  
    ColorAlgebra(Dx, 0, 2, 3, Dy, 1, local_site_offset, subtract);
    ColorAlgebra(Dx, 1, 3, 3, Dy, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 0, 3, Dy, 1, local_site_offset, subtract);
    ColorAlgebra(Dx, 3, 1, 3, Dy, 1, local_site_offset, answer);

    break;    

  case GammaYGammaY:

    ColorAlgebra(Dx, 0, 2, 0, Dy, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 3, 0, Dy, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 0, 0, Dy, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 1, 0, Dy, 2, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 2, 1, Dy, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 3, 1, Dy, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 0, 1, Dy, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 1, 1, Dy, 3, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 2, 2, Dy, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 3, 2, Dy, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 0, 2, Dy, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 1, 2, Dy, 0, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 2, 3, Dy, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 3, 3, Dy, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 0, 3, Dy, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 1, 3, Dy, 1, local_site_offset, answer);

    break;
    
  case GammaZGammaZ:

    ColorAlgebra(Dx, 0, 3, 0, Dy, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 2, 0, Dy, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 1, 0, Dy, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 0, 0, Dy, 3, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 3, 1, Dy, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 2, 1, Dy, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 1, 1, Dy, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 0, 1, Dy, 2, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 3, 2, Dy, 1, local_site_offset, answer);      
    ColorAlgebra(Dx, 1, 2, 2, Dy, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 1, 2, Dy, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 0, 2, Dy, 1, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 3, 3, Dy, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 2, 3, Dy, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 1, 3, Dy, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 0, 3, Dy, 0, local_site_offset, answer);

    break;
    
  case GammaTGammaT:

    ColorAlgebra(Dx, 0, 3, 0, Dy, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 2, 0, Dy, 3, local_site_offset, subtract);
    ColorAlgebra(Dx, 2, 1, 0, Dy, 3, local_site_offset, subtract);
    ColorAlgebra(Dx, 3, 0, 0, Dy, 3, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 3, 1, Dy, 2, local_site_offset, subtract);
    ColorAlgebra(Dx, 1, 2, 1, Dy, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 1, 1, Dy, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 0, 1, Dy, 2, local_site_offset, subtract);
  
    ColorAlgebra(Dx, 0, 3, 2, Dy, 1, local_site_offset, subtract);
    ColorAlgebra(Dx, 1, 2, 2, Dy, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 1, 2, Dy, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 0, 2, Dy, 1, local_site_offset, subtract);
  
    ColorAlgebra(Dx, 0, 3, 3, Dy, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 2, 3, Dy, 0, local_site_offset, subtract);
    ColorAlgebra(Dx, 2, 1, 3, Dy, 0, local_site_offset, subtract);
    ColorAlgebra(Dx, 3, 0, 3, Dy, 0, local_site_offset, answer);

    break;
    
  case Gamma5Gamma5:

    ColorAlgebra(Dx, 0, 1, 0, Dy, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 0, 0, Dy, 1, local_site_offset, subtract);
    ColorAlgebra(Dx, 2, 3, 0, Dy, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 2, 0, Dy, 1, local_site_offset, subtract);
  
    ColorAlgebra(Dx, 0, 1, 1, Dy, 0, local_site_offset, subtract);
    ColorAlgebra(Dx, 1, 0, 1, Dy, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 3, 1, Dy, 0, local_site_offset, subtract);
    ColorAlgebra(Dx, 3, 2, 1, Dy, 0, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 1, 2, Dy, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 0, 2, Dy, 3, local_site_offset, subtract);
    ColorAlgebra(Dx, 2, 3, 2, Dy, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 2, 2, Dy, 3, local_site_offset, subtract);
  
    ColorAlgebra(Dx, 0, 1, 3, Dy, 2, local_site_offset, subtract);
    ColorAlgebra(Dx, 1, 0, 3, Dy, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 3, 3, Dy, 2, local_site_offset, subtract);
    ColorAlgebra(Dx, 3, 2, 3, Dy, 2, local_site_offset, answer);

    break;

  case UnitUnit:

    ColorAlgebra(Dx, 0, 1, 0, Dy, 1, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 0, 0, Dy, 1, local_site_offset, subtract);
    ColorAlgebra(Dx, 2, 3, 0, Dy, 1, local_site_offset, subtract);
    ColorAlgebra(Dx, 3, 2, 0, Dy, 1, local_site_offset, answer);
  
    ColorAlgebra(Dx, 0, 1, 1, Dy, 0, local_site_offset, subtract);
    ColorAlgebra(Dx, 1, 0, 1, Dy, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 3, 1, Dy, 0, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 2, 1, Dy, 0, local_site_offset, subtract);
  
    ColorAlgebra(Dx, 0, 1, 2, Dy, 3, local_site_offset, subtract);
    ColorAlgebra(Dx, 1, 0, 2, Dy, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 2, 3, 2, Dy, 3, local_site_offset, answer);
    ColorAlgebra(Dx, 3, 2, 2, Dy, 3, local_site_offset, subtract);
  
    ColorAlgebra(Dx, 0, 1, 3, Dy, 2, local_site_offset, answer);
    ColorAlgebra(Dx, 1, 0, 3, Dy, 2, local_site_offset, subtract);
    ColorAlgebra(Dx, 2, 3, 3, Dy, 2, local_site_offset, subtract);
    ColorAlgebra(Dx, 3, 2, 3, Dy, 2, local_site_offset, answer);

    break;
  }

  answer -= subtract;  
}



//---------------------------------------------------------------------------
// void
// WspectBaryonW::DiracAlgebra(int...)
//--------------------------------------------------------------------------- 
void
WspectBaryon::DiracAlgebra(int Dx, int Dy, int local_site_offset,
			   Complex &answer) 
{
  char *fname = "DiracAlgebra(i,i,i,C&)";
  // answer += DiracAlgebraABC
  DiracAlgebraABC(Dx, Dy, local_site_offset, answer);    

  switch (d_flavors) {
  case AAC:                  // answer += DiracAlgebraAAC
    DiracAlgebraAAC(Dx, Dy, local_site_offset, answer);  
    break;
  case AAA:                  // answer += 2 DiracAlgebraAAC   
    Complex tmp(0, 0);
    DiracAlgebraAAC(Dx, Dy, local_site_offset, tmp);
    answer += 2 * tmp;  
    break;
  }
}




  

//---------------------------------------------------------------------------
// void
// WspectBaryon::ZeroMomProject(...)
//---------------------------------------------------------------------------
Complex
WspectBaryon::ZeroMomProject(int Dx, int Dy, int lcl_time_slice)
{
  d_lclMin[d_prop_dir] = d_lclMax[d_prop_dir] = lcl_time_slice;

  Complex answer(0.0, 0.0);

  int lcl[LORENTZs];  
  for (lcl[0] = d_lclMin[0]; lcl[0] <= d_lclMax[0]; lcl[0]++) {
    for (lcl[1] = d_lclMin[1]; lcl[1] <= d_lclMax[1]; lcl[1]++) {
      for (lcl[2] = d_lclMin[2]; lcl[2] <= d_lclMax[2]; lcl[2]++) {
	for (lcl[3] = d_lclMin[3]; lcl[3] <= d_lclMax[3]; lcl[3]++) {
	  DiracAlgebra(Dx, Dy, siteOffset(lcl)*SPINORs, answer);
#ifdef  DEBUG_W_BARYON_PROJECT_IMAG
	  printf("site[%d %d %d %d] spin[%d%d]: [%g %g]\n",
		 lcl[0], lcl[1], lcl[2], lcl[3], dirac, dirac,
		 answer.real(), answer.imag());
#endif
	}
      }
    }
  }

  return answer;  
}


//---------------------------------------------------------------------------
// Float
// WspectBaryon::DiracProject(...)
//---------------------------------------------------------------------------
Float
WspectBaryon::DiracProjectUnit(int a_local_time_slice) 
{
  Complex answer(0.0, 0.0);
  
  for (int dirac = 0; dirac < DIRACs; ++dirac) {
    answer += ZeroMomProject(dirac, dirac, a_local_time_slice);
  }
  
  // Note: a complier bug ? The following command screw up things.
  //    return answer.real();  
  // For example in one of our complex [3.27444 -0.0132897], 
  //    3.76367e-09 is passed back instead of 3.27444.
  // So one has to assign answer.real() to some other temporary 
  //    IFloat before return.
  //-------------------------------------------------------------------------
  Float ret = answer.real();  
  return ret;  
}
 

//---------------------------------------------------------------------------
// IFloat
// WspectBaryon::DiracProject(...)
//---------------------------------------------------------------------------
Float
WspectBaryon::DiracProjectGamma(int a_local_time_slice) 
{
  Float answer = 0.0;

  switch (d_whr.dir()) {
  case 0:
    answer -= ZeroMomProject(3, 0, a_local_time_slice).imag();
    answer -= ZeroMomProject(2, 1, a_local_time_slice).imag();
    answer += ZeroMomProject(1, 2, a_local_time_slice).imag();
    answer += ZeroMomProject(0, 3, a_local_time_slice).imag();
    break; 

  case 1:
    answer -= ZeroMomProject(3, 0, a_local_time_slice).real();
    answer += ZeroMomProject(2, 1, a_local_time_slice).real();
    answer += ZeroMomProject(1, 2, a_local_time_slice).real();
    answer -= ZeroMomProject(0, 3, a_local_time_slice).real();
    break; 

  case 2:
    answer -= ZeroMomProject(2, 0, a_local_time_slice).imag();
    answer += ZeroMomProject(3, 1, a_local_time_slice).imag();
    answer += ZeroMomProject(0, 2, a_local_time_slice).imag();
    answer -= ZeroMomProject(1, 3, a_local_time_slice).imag();
    break;

  case 3:
    answer += ZeroMomProject(2, 0, a_local_time_slice).real();
    answer += ZeroMomProject(3, 1, a_local_time_slice).real();
    answer += ZeroMomProject(0, 2, a_local_time_slice).real();
    answer += ZeroMomProject(1, 3, a_local_time_slice).real();
    break;
  } 

  return answer;  
}
 

//---------------------------------------------------------------------------
// void WspectBaryon::print(const CommonArg *common_arg) const
//---------------------------------------------------------------------------
void
WspectBaryon::print(char *filename, WbaryonFold fold) const
{

#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  char *fname = "print";
  VRB.Func(d_class_name, fname);

  // Open the output file
  //------------------------------------------------------------------
  FILE *fp=NULL;
  if (!filename || !(fp = Fopen(filename, "a"))) 
    ERR.FileA(d_class_name,fname, filename);


  // Right now d_unit_p  = Tr (Unit    Baryon Baryon^bar)
  //           d_gamma_p = Tr (Gamma   Baryon Baryon^bar)
  //           d_up_p    = Tr (1+Gamma Baryon Baryon^bar)
  //           d_down_p  = Tr (1-Gamma Baryon Baryon^bar)
  //-------------------------------------------------------------------------
  int glb_walls = glb_sites[d_whr.dir()];
  int src_wall  = d_whr.glbCoord();  
  int bound_cnd   = bnd_cnd[d_whr.dir()];
  Float f1 = 1.0;
  Float f2 = 1.0;


  // Fold the time slices
  //-------------------------------------------------------------------------
  int w;
  switch (fold) {
  case BARYON_PAST:
    Fprintf(fp, "%d 0 0 %e\n", AlgWspect::GetCounter(), 
	    IFloat(d_unit_p[src_wall]));
    for (w = 1; w < glb_walls/2; ++w){
      if (bound_cnd == -1) {
	f1 = (src_wall + w)/glb_walls             ? -1 : 1;
	f2 = (src_wall + glb_walls - w)/glb_walls ? -1 : 1;
      }
      Fprintf(fp, "%d %d 0 %e\n", AlgWspect::GetCounter(), w,
	      IFloat((f1 * d_unit_p[(src_wall + w)%glb_walls] +
		     f2 * bound_cnd * 
		     d_unit_p[(src_wall + glb_walls - w)%glb_walls])/2) );
    }
    Fprintf(fp,"%d %d 0 %e\n", AlgWspect::GetCounter(), w,
	    IFloat(d_unit_p[(src_wall+w)%glb_walls]));
    break;
    
  case BARYON_RAW:
    for (w = 0; w < glb_walls; ++w) {
      Fprintf(fp, "%d %d 0 %e %e\n", AlgWspect::GetCounter(), w,
      	      IFloat(d_up_p  [(src_wall + w)%glb_walls]),
	      IFloat(d_down_p[(src_wall + w)%glb_walls]));
    }
    break;
    
  case BARYON_FOLD:
  default:
    for (w = 0; w < glb_walls; ++w) {
      if (bound_cnd == -1) {
	f1 = (src_wall + w)/glb_walls             ? -1 : 1;
	f2 = (src_wall + glb_walls - w)/glb_walls ? -1 : 1;
      }
      Fprintf(fp, "%d %d 0 %e \n", AlgWspect::GetCounter(), w,
	      IFloat((f1 * d_up_p[(src_wall + w)%glb_walls] + 
		     f2 * bound_cnd *
		     d_down_p[(src_wall + glb_walls - w)%glb_walls])/2.0));
    }
    break;
  }
  
  Fclose(fp);
}




CPS_END_NAMESPACE
