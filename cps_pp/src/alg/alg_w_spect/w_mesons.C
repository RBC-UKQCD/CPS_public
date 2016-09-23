#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008/02/08 18:35:06 $
//  $Header: /space/cvs/cps/cps++/src/alg/alg_w_spect/w_mesons.C,v 1.12 2008/02/08 18:35:06 chulwoo Exp $
//  $Id: w_mesons.C,v 1.12 2008/02/08 18:35:06 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: w_mesons.C,v $
//  $Revision: 1.12 $
//  $Source: /space/cvs/cps/cps++/src/alg/alg_w_spect/w_mesons.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <alg/w_all.h>
#include <alg/w_gamma_mat.h>
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif

#include <util/error.h>                // ERR
#include <util/verbose.h>              // VRB
#include <util/vector.h>               // dotProduct
#include <util/qcdio.h> 
#include <comms/glb.h>                   // glb_sum
#include <alg/alg_w_spect.h>          // AlgWspect::GetCounter()
CPS_START_NAMESPACE

//---------------------------------------------------------------------------
// For the purpose of debugging or timing during code upgrade
//---------------------------------------------------------------------------
//#define DEBUG_W_MESON
//#define TIME_W_MESON

#ifdef  TIME_W_MESON
CPS_END_NAMESPACE
#include <time.h>                 // clock()
CPS_START_NAMESPACE
#endif


//---------------------------------------------------------------------------
// static data members
//---------------------------------------------------------------------------
char *      WspectMesons::d_class_name = "WspectMesons";


//---------------------------------------------------------------------------
// WspectMesons::WspectMesons(...)
//--------------------------------------------------------------------------- 
WspectMesons::WspectMesons(const IFloat *quark1, const IFloat *quark2, 
			   const WspectHyperRectangle &whr,
			   const WspectMomenta & mom)
  : d_quark1_p(quark1),
    d_quark2_p(quark2),
    d_whr(whr),
    d_mom(mom)
{
  Everything();
}




//---------------------------------------------------------------------------
// WspectMesons::WspectMesons(...)
//--------------------------------------------------------------------------- 
WspectMesons::WspectMesons(const WspectQuark &quark1, 
			   const WspectQuark &quark2,
			   const WspectHyperRectangle &whr,
			   const WspectMomenta & mom)
  : d_quark1_p(quark1.Data()),
    d_quark2_p(quark2.Data()),
    d_whr(whr),
    d_mom(mom)
{
  Everything();
}



//---------------------------------------------------------------------------
// WspectMesons::Everything()
//--------------------------------------------------------------------------- 
void
WspectMesons::Everything()
{
  VRB.Func(d_class_name, ctor_str);

#ifdef TIME_W_MESON
  int beforeClock = clock();
#endif

  //  Length in the propagation direction
  //-------------------------------------------------------------------------
  d_prop_dir  = d_whr.dir();  
  d_glb_walls = glb_sites[d_prop_dir];
  {
    const int *low  = d_whr.lclMin();
    const int *high = d_whr.lclMax();
    for (int i = 0; i < LORENTZs; ++i) {
      d_lclMin[i] = low[i];
      d_lclMax[i] = high[i];
    }
  }


  // d_num_mom
  //-------------------------------------------------------------------------
  d_num_mom = 1;
  if (d_mom)
    d_num_mom += d_mom.numNonZeroMom();


  // allocate space and clear them
  //-------------------------------------------------------------------------
  d_size   = d_num_mom * d_glb_walls * MESONs * COMPLEXs;
  {
    int fsize = d_size * sizeof(IFloat);
    d_data_p = (Complex *)smalloc(fsize);
    if (!d_data_p)
      ERR.Pointer(d_class_name, ctor_str, empty_str);
    VRB.Smalloc(d_class_name, ctor_str, empty_str, d_data_p, fsize);

    IFloat *flt_p = (IFloat *)d_data_p;
    for (int i = 0; i < d_size; ++i)
      *flt_p++ = 0.0;
  }

  // Set the on-node part of d_data_p 
  //-------------------------------------------------------------------------
  {
    int lcl_walls = lcl_sites[d_prop_dir];

    for (int lclW = 0; lclW < lcl_walls; ++lclW) {
      DiracAlgebra(lclW);
    }

  }


#ifdef DEBUG_W_MESON
  dumpData("before.gsum.meson.dat");
#endif
  

  // Global sum over all data
  //-------------------------------------------------------------------------
  {
    Float *Flt_p = (Float *)d_data_p;
    for (int i = 0; i < d_size; ++i)    glb_sum(Flt_p++);
  }


#ifdef DEBUG_W_MESON
  dumpData("after.gsum.meson.Mdat");  
#endif


#ifdef TIME_W_MESON
  int afterClock = clock();
  int clockTime = afterClock - beforeClock;  
  printf("Meson %d = [%d - %d]", clockTime, afterClock, beforeClock);  
#endif
  
}

//---------------------------------------------------------------------------
// WspectMesons::~WspectMesons()
//--------------------------------------------------------------------------- 
// Purpose:
//    Free all the memory on heap.
//--------------------------------------------------------------------------- 
WspectMesons::~WspectMesons()
{
  VRB.Func(d_class_name, dtor_str);
  VRB.Sfree(d_class_name, dtor_str, empty_str, d_data_p);
  sfree(d_data_p);
}


//---------------------------------------------------------------------------
// void WspectMesons::ColorAlgebra(...)
//---------------------------------------------------------------------------
// The quark propagator data is stored as
//        IFloat[Dy][Cy][T][Z][Y][X][Dx][Cx][2]
// if we only need the real part, we could use
//     dotProduct(q1_p, q2_p, COMPLEXs*COLORs);
// which is much more efficient (by a rough factor of 6).
//---------------------------------------------------------------------------
void
WspectMesons::ColorAlgebra(int D1x, int D2x, int D1y, int D2y,
			   const int lcl[LORENTZs], Complex &result) const
{
  int local_site_offset = siteOffset(lcl) * SPINORs;  

  const Complex *q1_p = (const Complex *)(d_quark1_p + 
					  WspectQuark::weightSrcDirac()*D1y + 
					  local_site_offset + 
					  (COMPLEXs*COLORs) * D1x);
  const Complex *q2_p = (const Complex *)(d_quark2_p + 
					  WspectQuark::weightSrcDirac()*D2y + 
					  local_site_offset + 
					  (COMPLEXs*COLORs) * D2x);

  int off_Cy_in_complex = WspectQuark::weightSrcColor() / COMPLEXs;

  for (int Cy = 0; Cy < COLORs; Cy++) {
    result += q1_p[0] * conj(q2_p[0]);
    result += q1_p[1] * conj(q2_p[1]);
    result += q1_p[2] * conj(q2_p[2]);
    q1_p += off_Cy_in_complex;    
    q2_p += off_Cy_in_complex;    
  }
}


//---------------------------------------------------------------------------
// void WspectMesons::MomProject(...)
//---------------------------------------------------------------------------
// res[D1x][D1y][D2x][D2y] 
//---------------------------------------------------------------------------
void
WspectMesons::MomProject(int D1x, int D2x, int D1y, int D2y, int lclW)
{
  VRB.Func(d_class_name, "MomProject");
  
  // set the coordinate limit in the prop_dir to avoid looping in that dir.
  //-------------------------------------------------------------------------
  d_lclMin[d_prop_dir] = d_lclMax[d_prop_dir] = lclW; 
  
  // loop over local sites [first for the efficiency in calculation of
  // local_site_offset and non-zero momentum projection] and Dirac indexes.
  //-------------------------------------------------------------------------
  int lcl[LORENTZs];                         
  for (lcl[0] = d_lclMin[0]; lcl[0] <= d_lclMax[0]; lcl[0]++) {
    for (lcl[1] = d_lclMin[1]; lcl[1] <= d_lclMax[1]; lcl[1]++) {	
      for (lcl[2] = d_lclMin[2]; lcl[2] <= d_lclMax[2]; lcl[2]++) {
	for (lcl[3] = d_lclMin[3]; lcl[3] <= d_lclMax[3]; lcl[3]++) {
	  
	  Complex colorSum(0.0, 0.0);		
	  ColorAlgebra(D1x, D2x, D1y, D2y, lcl, colorSum);

	  // zero momentum projection
	  d_mom_proj[0][D1x][D1y][D2x][D2y] += colorSum;  

	  // non-zero meomentum projection
	  if (d_mom) {	  
	    const Complex *mom = d_mom[lcl];
	    for (int iMom = 1; iMom < d_num_mom; ++iMom) {
	      d_mom_proj[iMom][D1x][D1y][D2x][D2y] += mom[iMom-1] * colorSum;  
	    }
	  }
	}       // for (lcl[3] = d_lclMin[3];
      }         // for (lcl[2] = d_lclMin[2]
    }           // for (lcl[1] = d_lclMin[1];
  }             // for (lcl[0] = d_lclMin[0];
}


//Used in DiracAlgebra

void WspectMesons::traceDirac(int sign, IFloat *gam1, IFloat *gam2,int mom,
				Complex *result_p){
  char *fname="traceDirac";
   //mesonId gives the offset into correlator result
  //coor_data_p[glbwall][
  //pass the offset from outside directly?
  //int corr_offset=mesonId+(lclw+lcl2glb_offset[d_prop_dir])*EXTMESONs;
 
  //check for errors
  if(sign!=1 && sign !=-1) {
    ERR.General(d_class_name, fname,"sign must be 1 or -1");
  }
  if(DIRACs!=4) ERR.General(d_class_name,fname,"DIRACs must be 4!");
 
  if(gam1==0 || gam2==0 || result_p==0) {
    ERR.General(d_class_name,fname,"null pointer passed");
  }

  for(int d1=0;d1<DIRACs;d1++){
    for(int d2=0;d2<DIRACs;d2++){
      for(int d3=0;d3<DIRACs;d3++){
        for(int d4=0;d4<DIRACs;d4++){

          Complex gam1_element(gam1[d1*4*2+d2*2],gam1[d1*4*2+d2*2+1]);
          Complex gam2_element(gam2[d3*4*2+d4*2],gam2[d3*4*2+d4*2+1]);
          //Complex tmp;
          //tmp=sign*gam1_element*d_zero_mom_proj[d2][d3][d1][d4]*gam2_element;
          //if(tmp.real()!=0 && tmp.imag()!=0) 
   	  //printf("Non-zero term %d,%d,%d,%d\n",d2,d3,d1,d4);

          if(sign==1){
            *(result_p) += 
 		gam1_element*gam2_element*d_mom_proj[mom][d2][d3][d1][d4];
          }else{
            *(result_p) -= 
		gam1_element*gam2_element*d_mom_proj[mom][d2][d3][d1][d4];
          }
        }//d1
      }//d2
    }//d3
  }//d4
 
  //done
}


//---------------------------------------------------------------------------
// void WspectMesons::DiracAlgebra(...)
//---------------------------------------------------------------------------
void
WspectMesons::DiracAlgebra(int lclW) 
{ 
  VRB.Func(d_class_name, "DiracAlgebra");

  // Loop over Dirac indexes:    Momentum Projection
  //
  // Ping -- we seem to be able to approve that why only the following
  //         combinations of spinor indexes will ever be used in 
  //         current calcuations. Can we prove it for non-degenerate
  //         quark masses??????
  //-----------------------------------------------------------------------
  {
    // clear space for result -- more effecient this way
    IFloat *p = (IFloat *)d_mom_proj;
    for (int  i = 0; i < sizeof(d_mom_proj)/sizeof(IFloat); ++i)
      *p++ = 0.0;

    // loop
    int D1x, D2x, D1y, D2y;  
    for (D1x = 0; D1x < DIRACs; D1x++) {    
      for (D1y = 0; D1y < DIRACs; D1y++) {      
	for (D2x = 0; D2x < DIRACs; D2x++) {
	  if (D1x==D1y) 
	    D2y=D2x;            // D1x==D1y  D2x==D2y
	  else if (D1x==D2x) 
	    D2y=D1y;            // D1x==D2x  D1y==D2y
	  else if (D1y==D2x) 
	    D2y=D1x;            // D1x==D2y  D1y==D2x
	  else {                    // D1x!=D2y!=D1y!=D2x
	    for (D2y = 0; D2y == D1x || D2y == D1y || D2y == D2x; ++D2y);
	  } 
	  MomProject(D1x, D2x, D1y, D2y, lclW);                   
	}
      }
    }
  }
  
 
#ifdef DEBUG_W_MESON
  {
    static int calls = 0; calls++;
    FILE *fp = Fopen("meson.dirac.dat", "a");    
    Fprintf(fp, "Call %d\n", calls);
    for (int i0 = 0; i0 < DIRACs; ++i0) {
      for (int i1 = 0; i1 < DIRACs; ++i1) {
	for (int i2 = 0; i2 < DIRACs; ++i2) {
	  for (int i3 = 0; i3 < DIRACs; ++i3) {
	    for (int iMom = 0; iMom < d_num_mom; ++iMom) {
	      const Complex & tmp = d_mom_proj[iMom][i0][i1][i2][i3];
	      Fprintf(fp, "mom %d %d %d %d %d [%g %g]\n", 
		      iMom, i0, i1, i2, i3, tmp.real(), tmp.imag());
	    }
	  }
	}
      }
    }
    Fclose(fp);
  }
#endif
 

  // point to the data for this global time slice
  //-----------------------------------------------------------------------
  Complex *mesonsAddr = d_data_p + (lclW + lcl2glb_offset[d_prop_dir]
				    ) * MESONs;

  // mesonsAddr increment for each momentum loop:
  //     each momentum:       d_glb_walls * MESONs
  //     within a iMom loop:  MESONs - 1   increments already
  //           subtract:     (d_glb_walls-1) * MESONs + 1;
  //-----------------------------------------------------------------------
  int numComplexToNextMom = (d_glb_walls-1) * MESONs + 1;

  // for each momentum and each meson
  //-----------------------------------------------------------------------
  for (int iMom = 0; iMom < d_num_mom; ++iMom) {    
    // Meson 0   GAMMA = [gamma_t^0][gamma_z^0][gy^0][gx^0]
    //traceDirac(1,(IFloat *)WG5, (IFloat *)WGamma5, iMom, mesonsAddr);

    *mesonsAddr    = d_mom_proj[iMom][0][0][0][0];
    *mesonsAddr   += d_mom_proj[iMom][0][1][0][1];
    *mesonsAddr   -= d_mom_proj[iMom][0][2][0][2];
    *mesonsAddr   -= d_mom_proj[iMom][0][3][0][3];
    *mesonsAddr   += d_mom_proj[iMom][1][0][1][0];
    *mesonsAddr   += d_mom_proj[iMom][1][1][1][1];
    *mesonsAddr   -= d_mom_proj[iMom][1][2][1][2];
    *mesonsAddr   -= d_mom_proj[iMom][1][3][1][3];
    *mesonsAddr   -= d_mom_proj[iMom][2][0][2][0];
    *mesonsAddr   -= d_mom_proj[iMom][2][1][2][1];
    *mesonsAddr   += d_mom_proj[iMom][2][2][2][2];
    *mesonsAddr   += d_mom_proj[iMom][2][3][2][3];
    *mesonsAddr   -= d_mom_proj[iMom][3][0][3][0];
    *mesonsAddr   -= d_mom_proj[iMom][3][1][3][1];
    *mesonsAddr   += d_mom_proj[iMom][3][2][3][2];
    *mesonsAddr++ += d_mom_proj[iMom][3][3][3][3]; 
    // Meson 1  GAMMA = [gt^0][gz^0][gy^0][gx^1]   Catalin's 8
    //mesonsAddr++;
    //traceDirac(-1,(IFloat *)WGamma1Gamma5, (IFloat *)WGamma1Gamma5, iMom, mesonsAddr);

    *mesonsAddr    = d_mom_proj[iMom][0][0][3][3];
    *mesonsAddr   += d_mom_proj[iMom][0][1][3][2];
    *mesonsAddr   += d_mom_proj[iMom][0][2][3][1];
    *mesonsAddr   += d_mom_proj[iMom][0][3][3][0];
    *mesonsAddr   += d_mom_proj[iMom][1][0][2][3];
    *mesonsAddr   += d_mom_proj[iMom][1][1][2][2];
    *mesonsAddr   += d_mom_proj[iMom][1][2][2][1];
    *mesonsAddr   += d_mom_proj[iMom][1][3][2][0];
    *mesonsAddr   += d_mom_proj[iMom][2][0][1][3];
    *mesonsAddr   += d_mom_proj[iMom][2][1][1][2];
    *mesonsAddr   += d_mom_proj[iMom][2][2][1][1];
    *mesonsAddr   += d_mom_proj[iMom][2][3][1][0];
    *mesonsAddr   += d_mom_proj[iMom][3][0][0][3];
    *mesonsAddr   += d_mom_proj[iMom][3][1][0][2];
    *mesonsAddr   += d_mom_proj[iMom][3][2][0][1];
    *mesonsAddr++ += d_mom_proj[iMom][3][3][0][0];     

    // Meson 2  GAMMA = [gt^0][gz^0][gy^1][gx^0]   Catalin's 1
    //mesonsAddr++;
    //traceDirac(-1,(IFloat *)WGamma2Gamma5, (IFloat *)WGamma2Gamma5, iMom, mesonsAddr);

    *mesonsAddr    = d_mom_proj[iMom][0][0][3][3];
    *mesonsAddr   -= d_mom_proj[iMom][0][1][3][2];
    *mesonsAddr   += d_mom_proj[iMom][0][2][3][1];
    *mesonsAddr   -= d_mom_proj[iMom][0][3][3][0];
    *mesonsAddr   -= d_mom_proj[iMom][1][0][2][3];
    *mesonsAddr   += d_mom_proj[iMom][1][1][2][2];
    *mesonsAddr   -= d_mom_proj[iMom][1][2][2][1];
    *mesonsAddr   += d_mom_proj[iMom][1][3][2][0];
    *mesonsAddr   += d_mom_proj[iMom][2][0][1][3];
    *mesonsAddr   -= d_mom_proj[iMom][2][1][1][2];
    *mesonsAddr   += d_mom_proj[iMom][2][2][1][1];
    *mesonsAddr   -= d_mom_proj[iMom][2][3][1][0];
    *mesonsAddr   -= d_mom_proj[iMom][3][0][0][3];
    *mesonsAddr   += d_mom_proj[iMom][3][1][0][2];
    *mesonsAddr   -= d_mom_proj[iMom][3][2][0][1];
    *mesonsAddr++ += d_mom_proj[iMom][3][3][0][0];     


    // Meson 3:   GAMMA =  [gt^0][gz^0][gy^1][gx^1]  Catalin's 9
    //mesonsAddr++;
    //traceDirac(-1,(IFloat *)WGamma3Gamma4, (IFloat *)WGamma3Gamma4, iMom,mesonsAddr);

    *mesonsAddr    = d_mom_proj[iMom][0][1][0][1];
    *mesonsAddr   -= d_mom_proj[iMom][0][0][0][0];
    *mesonsAddr   += d_mom_proj[iMom][0][2][0][2];
    *mesonsAddr   -= d_mom_proj[iMom][0][3][0][3];
    *mesonsAddr   += d_mom_proj[iMom][1][0][1][0];
    *mesonsAddr   -= d_mom_proj[iMom][1][1][1][1];
    *mesonsAddr   -= d_mom_proj[iMom][1][2][1][2];
    *mesonsAddr   += d_mom_proj[iMom][1][3][1][3];
    *mesonsAddr   += d_mom_proj[iMom][2][0][2][0];
    *mesonsAddr   -= d_mom_proj[iMom][2][1][2][1];
    *mesonsAddr   -= d_mom_proj[iMom][2][2][2][2];
    *mesonsAddr   += d_mom_proj[iMom][2][3][2][3];
    *mesonsAddr   -= d_mom_proj[iMom][3][0][3][0];
    *mesonsAddr   += d_mom_proj[iMom][3][1][3][1];
    *mesonsAddr   += d_mom_proj[iMom][3][2][3][2];
    *mesonsAddr++ -= d_mom_proj[iMom][3][3][3][3];     

    // Meson 4:     GAMMA = [gt^0][gz^1][gy^0][gx^0]   Catalin's 2
    //mesonsAddr++;
    //traceDirac(-1,(IFloat *)WGamma3Gamma5, (IFloat *)WGamma3Gamma5, iMom,mesonsAddr);

    *mesonsAddr    = d_mom_proj[iMom][0][0][2][2];
    *mesonsAddr   -= d_mom_proj[iMom][0][1][2][3];
    *mesonsAddr   += d_mom_proj[iMom][0][2][2][0];
    *mesonsAddr   -= d_mom_proj[iMom][0][3][2][1];
    *mesonsAddr   -= d_mom_proj[iMom][1][0][3][2];
    *mesonsAddr   += d_mom_proj[iMom][1][1][3][3];
    *mesonsAddr   -= d_mom_proj[iMom][1][2][3][0];
    *mesonsAddr   += d_mom_proj[iMom][1][3][3][1];
    *mesonsAddr   += d_mom_proj[iMom][2][0][0][2];
    *mesonsAddr   -= d_mom_proj[iMom][2][1][0][3];
    *mesonsAddr   += d_mom_proj[iMom][2][2][0][0];
    *mesonsAddr   -= d_mom_proj[iMom][2][3][0][1];
    *mesonsAddr   -= d_mom_proj[iMom][3][0][1][2];
    *mesonsAddr   += d_mom_proj[iMom][3][1][1][3];
    *mesonsAddr   -= d_mom_proj[iMom][3][2][1][0];
    *mesonsAddr++ += d_mom_proj[iMom][3][3][1][1]; 

    // Meson 5:   GAMMA = [gt^0][gz^1][gy^0][gx^1]     Catalin's 10
    //mesonsAddr++;
    //traceDirac(-1,(IFloat *)WGamma2Gamma4, (IFloat *)WGamma2Gamma4, iMom,mesonsAddr);

    *mesonsAddr    = d_mom_proj[iMom][0][1][1][0];
    *mesonsAddr   -= d_mom_proj[iMom][0][0][1][1];
    *mesonsAddr   += d_mom_proj[iMom][0][2][1][3];
    *mesonsAddr   -= d_mom_proj[iMom][0][3][1][2];
    *mesonsAddr   += d_mom_proj[iMom][1][0][0][1];
    *mesonsAddr   -= d_mom_proj[iMom][1][1][0][0];
    *mesonsAddr   -= d_mom_proj[iMom][1][2][0][3];
    *mesonsAddr   += d_mom_proj[iMom][1][3][0][2];
    *mesonsAddr   += d_mom_proj[iMom][2][0][3][1];
    *mesonsAddr   -= d_mom_proj[iMom][2][1][3][0];
    *mesonsAddr   -= d_mom_proj[iMom][2][2][3][3];
    *mesonsAddr   += d_mom_proj[iMom][2][3][3][2];
    *mesonsAddr   -= d_mom_proj[iMom][3][0][2][1];
    *mesonsAddr   += d_mom_proj[iMom][3][1][2][0];
    *mesonsAddr   += d_mom_proj[iMom][3][2][2][3];
    *mesonsAddr++ -= d_mom_proj[iMom][3][3][2][2]; 

    // Meson 6:     GAMMA = [gt^0][gz^0][gy^1][gx^1]  Catalin's 3
    //mesonsAddr++;
    //traceDirac(-1,(IFloat *)WGamma1Gamma4, (IFloat *)WGamma1Gamma4, iMom,mesonsAddr);

    *mesonsAddr    = d_mom_proj[iMom][0][2][1][3];
    *mesonsAddr   -= d_mom_proj[iMom][0][0][1][1];
    *mesonsAddr   -= d_mom_proj[iMom][0][1][1][0];
    *mesonsAddr   += d_mom_proj[iMom][0][3][1][2];
    *mesonsAddr   -= d_mom_proj[iMom][1][0][0][1];
    *mesonsAddr   -= d_mom_proj[iMom][1][1][0][0];
    *mesonsAddr   += d_mom_proj[iMom][1][2][0][3];
    *mesonsAddr   += d_mom_proj[iMom][1][3][0][2];
    *mesonsAddr   += d_mom_proj[iMom][2][0][3][1];
    *mesonsAddr   += d_mom_proj[iMom][2][1][3][0];
    *mesonsAddr   -= d_mom_proj[iMom][2][2][3][3];
    *mesonsAddr   -= d_mom_proj[iMom][2][3][3][2];
    *mesonsAddr   += d_mom_proj[iMom][3][0][2][1];
    *mesonsAddr   += d_mom_proj[iMom][3][1][2][0];
    *mesonsAddr   -= d_mom_proj[iMom][3][2][2][3];
    *mesonsAddr++ -= d_mom_proj[iMom][3][3][2][2];     

    // Meson 7:   GAMMA = [gt^0][gz^1][gy^1][gx^1]  Catalin's 11
    //mesonsAddr++;
    //traceDirac(-1,(IFloat *)WGamma4, (IFloat *)WGamma4, iMom, mesonsAddr);

    *mesonsAddr    = 0;    
    *mesonsAddr   -= d_mom_proj[iMom][0][0][2][2];
    *mesonsAddr   -= d_mom_proj[iMom][0][1][2][3];
    *mesonsAddr   -= d_mom_proj[iMom][0][2][2][0];
    *mesonsAddr   -= d_mom_proj[iMom][0][3][2][1];
    *mesonsAddr   -= d_mom_proj[iMom][1][0][3][2];
    *mesonsAddr   -= d_mom_proj[iMom][1][1][3][3];
    *mesonsAddr   -= d_mom_proj[iMom][1][2][3][0];
    *mesonsAddr   -= d_mom_proj[iMom][1][3][3][1];
    *mesonsAddr   -= d_mom_proj[iMom][2][0][0][2];
    *mesonsAddr   -= d_mom_proj[iMom][2][1][0][3];
    *mesonsAddr   -= d_mom_proj[iMom][2][2][0][0];
    *mesonsAddr   -= d_mom_proj[iMom][2][3][0][1];
    *mesonsAddr   -= d_mom_proj[iMom][3][0][1][2];
    *mesonsAddr   -= d_mom_proj[iMom][3][1][1][3];
    *mesonsAddr   -= d_mom_proj[iMom][3][2][1][0];
    *mesonsAddr++ -= d_mom_proj[iMom][3][3][1][1];     

    // Meson 8:     GAMMA = [gt^1][gz^0][gy^0][gx^0] Catalin's 4
    //mesonsAddr++;
    //traceDirac(1,(IFloat *)WGamma4Gamma5, (IFloat *)WGamma4Gamma5, iMom,mesonsAddr);

    *mesonsAddr    = d_mom_proj[iMom][0][0][2][2];
    *mesonsAddr   += d_mom_proj[iMom][0][1][2][3];
    *mesonsAddr   -= d_mom_proj[iMom][0][2][2][0];
    *mesonsAddr   -= d_mom_proj[iMom][0][3][2][1];
    *mesonsAddr   += d_mom_proj[iMom][1][0][3][2];
    *mesonsAddr   += d_mom_proj[iMom][1][1][3][3];
    *mesonsAddr   -= d_mom_proj[iMom][1][2][3][0];
    *mesonsAddr   -= d_mom_proj[iMom][1][3][3][1];
    *mesonsAddr   -= d_mom_proj[iMom][2][0][0][2];
    *mesonsAddr   -= d_mom_proj[iMom][2][1][0][3];
    *mesonsAddr   += d_mom_proj[iMom][2][2][0][0];
    *mesonsAddr   += d_mom_proj[iMom][2][3][0][1];
    *mesonsAddr   -= d_mom_proj[iMom][3][0][1][2];
    *mesonsAddr   -= d_mom_proj[iMom][3][1][1][3];
    *mesonsAddr   += d_mom_proj[iMom][3][2][1][0];
    *mesonsAddr++ += d_mom_proj[iMom][3][3][1][1]; 

    // Meson 9:   GAMMA = [gt^1][gz^0][gy^0][gx^1]   Catalin's 12
    //mesonsAddr++;
    //traceDirac(1,(IFloat *)WGamma1Gamma4Gamma5, (IFloat *)WGamma1Gamma4Gamma5,iMom, mesonsAddr);

    *mesonsAddr    = 0;    
    *mesonsAddr   -= d_mom_proj[iMom][0][0][1][1];
    *mesonsAddr   -= d_mom_proj[iMom][0][1][1][0];
    *mesonsAddr   -= d_mom_proj[iMom][0][2][1][3];
    *mesonsAddr   -= d_mom_proj[iMom][0][3][1][2];
    *mesonsAddr   -= d_mom_proj[iMom][1][0][0][1];
    *mesonsAddr   -= d_mom_proj[iMom][1][1][0][0];
    *mesonsAddr   -= d_mom_proj[iMom][1][2][0][3];
    *mesonsAddr   -= d_mom_proj[iMom][1][3][0][2];
    *mesonsAddr   -= d_mom_proj[iMom][2][0][3][1];
    *mesonsAddr   -= d_mom_proj[iMom][2][1][3][0];
    *mesonsAddr   -= d_mom_proj[iMom][2][2][3][3];
    *mesonsAddr   -= d_mom_proj[iMom][2][3][3][2];
    *mesonsAddr   -= d_mom_proj[iMom][3][0][2][1];
    *mesonsAddr   -= d_mom_proj[iMom][3][1][2][0];
    *mesonsAddr   -= d_mom_proj[iMom][3][2][2][3];
    *mesonsAddr++ -= d_mom_proj[iMom][3][3][2][2];     

    // Meson 10:    GAMMA = [gt^1][gz^0][gy^1][gx^0]  Catalin's 5
    //mesonsAddr++;
    //traceDirac(1,(IFloat *)WGamma2Gamma4Gamma5, (IFloat *)WGamma2Gamma4Gamma5,iMom, mesonsAddr);

    *mesonsAddr    = d_mom_proj[iMom][0][1][1][0];
    *mesonsAddr   -= d_mom_proj[iMom][0][0][1][1];
    *mesonsAddr   -= d_mom_proj[iMom][0][2][1][3];
    *mesonsAddr   += d_mom_proj[iMom][0][3][1][2];
    *mesonsAddr   += d_mom_proj[iMom][1][0][0][1];
    *mesonsAddr   -= d_mom_proj[iMom][1][1][0][0];
    *mesonsAddr   += d_mom_proj[iMom][1][2][0][3];
    *mesonsAddr   -= d_mom_proj[iMom][1][3][0][2];
    *mesonsAddr   -= d_mom_proj[iMom][2][0][3][1];
    *mesonsAddr   += d_mom_proj[iMom][2][1][3][0];
    *mesonsAddr   -= d_mom_proj[iMom][2][2][3][3];
    *mesonsAddr   += d_mom_proj[iMom][2][3][3][2];
    *mesonsAddr   += d_mom_proj[iMom][3][0][2][1];
    *mesonsAddr   -= d_mom_proj[iMom][3][1][2][0];
    *mesonsAddr   += d_mom_proj[iMom][3][2][2][3];
    *mesonsAddr++ -= d_mom_proj[iMom][3][3][2][2]; 

    // Meson 11:  GAMMA = [gt^1][gz^0][gy^1][gx^1]  Catalin's 13
    //mesonsAddr++;
    //traceDirac(-1,(IFloat *)WGamma3, (IFloat *)WGamma3, iMom, mesonsAddr);

    *mesonsAddr    = d_mom_proj[iMom][0][1][2][3];
    *mesonsAddr   -= d_mom_proj[iMom][0][0][2][2];
    *mesonsAddr   += d_mom_proj[iMom][0][2][2][0];
    *mesonsAddr   -= d_mom_proj[iMom][0][3][2][1];
    *mesonsAddr   += d_mom_proj[iMom][1][0][3][2];
    *mesonsAddr   -= d_mom_proj[iMom][1][1][3][3];
    *mesonsAddr   -= d_mom_proj[iMom][1][2][3][0];
    *mesonsAddr   += d_mom_proj[iMom][1][3][3][1];
    *mesonsAddr   += d_mom_proj[iMom][2][0][0][2];
    *mesonsAddr   -= d_mom_proj[iMom][2][1][0][3];
    *mesonsAddr   -= d_mom_proj[iMom][2][2][0][0];
    *mesonsAddr   += d_mom_proj[iMom][2][3][0][1];
    *mesonsAddr   -= d_mom_proj[iMom][3][0][1][2];
    *mesonsAddr   += d_mom_proj[iMom][3][1][1][3];
    *mesonsAddr   += d_mom_proj[iMom][3][2][1][0];
    *mesonsAddr++ -= d_mom_proj[iMom][3][3][1][1]; 

    // Meson 12:    GAMMA = [gt^1][gz^1][gy^0][gx^0] Catalin's 6
    //    mesonsAddr++;
    //traceDirac(1,(IFloat *)WGamma3Gamma4Gamma5, (IFloat *)WGamma3Gamma4Gamma5,iMom, mesonsAddr);

    *mesonsAddr    = d_mom_proj[iMom][0][1][0][1];
    *mesonsAddr   -= d_mom_proj[iMom][0][0][0][0];
    *mesonsAddr   -= d_mom_proj[iMom][0][2][0][2];
    *mesonsAddr   += d_mom_proj[iMom][0][3][0][3];
    *mesonsAddr   += d_mom_proj[iMom][1][0][1][0];
    *mesonsAddr   -= d_mom_proj[iMom][1][1][1][1];
    *mesonsAddr   += d_mom_proj[iMom][1][2][1][2];
    *mesonsAddr   -= d_mom_proj[iMom][1][3][1][3];
    *mesonsAddr   -= d_mom_proj[iMom][2][0][2][0];
    *mesonsAddr   += d_mom_proj[iMom][2][1][2][1];
    *mesonsAddr   -= d_mom_proj[iMom][2][2][2][2];
    *mesonsAddr   += d_mom_proj[iMom][2][3][2][3];
    *mesonsAddr   += d_mom_proj[iMom][3][0][3][0];
    *mesonsAddr   -= d_mom_proj[iMom][3][1][3][1];
    *mesonsAddr   += d_mom_proj[iMom][3][2][3][2];
    *mesonsAddr++ -= d_mom_proj[iMom][3][3][3][3];     

    // Meson 13:  GAMMA = [gt^1][gz^1][gy^0][gx^1]  Catalin's 14
    //mesonsAddr++;
    //traceDirac(-1,(IFloat *)WGamma2, (IFloat *)WGamma2, iMom, mesonsAddr);

    *mesonsAddr    = d_mom_proj[iMom][0][1][3][2];
    *mesonsAddr   -= d_mom_proj[iMom][0][0][3][3];
    *mesonsAddr   += d_mom_proj[iMom][0][2][3][1];
    *mesonsAddr   -= d_mom_proj[iMom][0][3][3][0];
    *mesonsAddr   += d_mom_proj[iMom][1][0][2][3];
    *mesonsAddr   -= d_mom_proj[iMom][1][1][2][2];
    *mesonsAddr   -= d_mom_proj[iMom][1][2][2][1];
    *mesonsAddr   += d_mom_proj[iMom][1][3][2][0];
    *mesonsAddr   += d_mom_proj[iMom][2][0][1][3];
    *mesonsAddr   -= d_mom_proj[iMom][2][1][1][2];
    *mesonsAddr   -= d_mom_proj[iMom][2][2][1][1];
    *mesonsAddr   += d_mom_proj[iMom][2][3][1][0];
    *mesonsAddr   -= d_mom_proj[iMom][3][0][0][3];
    *mesonsAddr   += d_mom_proj[iMom][3][1][0][2];
    *mesonsAddr   += d_mom_proj[iMom][3][2][0][1];
    *mesonsAddr++ -= d_mom_proj[iMom][3][3][0][0]; 

    // Meson 14:   GAMMA = [gt^1][gz^1][gy^1][gx^0] Catalin's 7
    // mesonsAddr++;
    //traceDirac(-1,(IFloat *)WGamma1, (IFloat *)WGamma1, iMom, mesonsAddr);

    *mesonsAddr    = d_mom_proj[iMom][0][2][3][1];
    *mesonsAddr   -= d_mom_proj[iMom][0][0][3][3];
    *mesonsAddr   -= d_mom_proj[iMom][0][1][3][2];
    *mesonsAddr   += d_mom_proj[iMom][0][3][3][0];
    *mesonsAddr   -= d_mom_proj[iMom][1][0][2][3];
    *mesonsAddr   -= d_mom_proj[iMom][1][1][2][2];
    *mesonsAddr   += d_mom_proj[iMom][1][2][2][1];
    *mesonsAddr   += d_mom_proj[iMom][1][3][2][0];
    *mesonsAddr   += d_mom_proj[iMom][2][0][1][3];
    *mesonsAddr   += d_mom_proj[iMom][2][1][1][2];
    *mesonsAddr   -= d_mom_proj[iMom][2][2][1][1];
    *mesonsAddr   -= d_mom_proj[iMom][2][3][1][0];
    *mesonsAddr   += d_mom_proj[iMom][3][0][0][3];
    *mesonsAddr   += d_mom_proj[iMom][3][1][0][2];
    *mesonsAddr   -= d_mom_proj[iMom][3][2][0][1];
    *mesonsAddr++ -= d_mom_proj[iMom][3][3][0][0]; 

    // Meson 15:  GAMMA = [gt^1][gz^1][gy^1][gx^1] Catalin's 15
    //mesonsAddr++;
    //traceDirac(1,(IFloat *)WUnitMat, (IFloat *)WUnitMat, iMom, mesonsAddr);

    *mesonsAddr    = d_mom_proj[iMom][0][0][0][0];
    *mesonsAddr   += d_mom_proj[iMom][0][1][0][1];
    *mesonsAddr   += d_mom_proj[iMom][0][2][0][2];
    *mesonsAddr   += d_mom_proj[iMom][0][3][0][3];
    *mesonsAddr   += d_mom_proj[iMom][1][0][1][0];
    *mesonsAddr   += d_mom_proj[iMom][1][1][1][1];
    *mesonsAddr   += d_mom_proj[iMom][1][2][1][2];
    *mesonsAddr   += d_mom_proj[iMom][1][3][1][3];
    *mesonsAddr   += d_mom_proj[iMom][2][0][2][0];
    *mesonsAddr   += d_mom_proj[iMom][2][1][2][1];
    *mesonsAddr   += d_mom_proj[iMom][2][2][2][2];
    *mesonsAddr   += d_mom_proj[iMom][2][3][2][3];
    *mesonsAddr   += d_mom_proj[iMom][3][0][3][0];
    *mesonsAddr   += d_mom_proj[iMom][3][1][3][1];
    *mesonsAddr   += d_mom_proj[iMom][3][2][3][2];
    *mesonsAddr   += d_mom_proj[iMom][3][3][3][3]; 

    // increment the pointer to the beginning of the next momentum
    mesonsAddr += numComplexToNextMom;    
  }
}


//---------------------------------------------------------------------------
// void WspectMesons::print(const CommonArg *common_arg) const
//--------------------------------------------------------------------------- 
void
WspectMesons::print(WspectOutput *w_spect_output) const
{

#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  char *fname = "print";
  VRB.Func(d_class_name, fname);

  // the first 8 enum's have to be in the same order as in WspectOutput
  enum {a0,           a0_prime,     a1,        b1,
        pion,         pion_prime,   rho,       rho_prime,
        NUM_DATAFILEs,
        a1_x,         a1_y,         a1_z,
        b1_x,         b1_y,         b1_z,
        rho_x,        rho_y,        rho_z,
        rho_x_prime,  rho_y_prime,  rho_z_prime,
        NUM_MAP
  };


  // map[x] : filenames x => meson data for this particular prop dir
  //-------------------------------------------------------------------------
  int map[NUM_MAP];
  {
    int time = d_whr.dir();
    int offset_time = (1 << time);

    int space = 0;    
    for (int l = 0; l < LORENTZs; ++l) {
      if (l != time) {
        int offset_space       = (1<<l);
        map[rho_x+space]       =              offset_space;
        map[rho_x_prime+space] =              offset_space + offset_time;
        map[b1_x+space]        = (MESONs-1) - offset_space - offset_time;
        map[a1_x+space]        = (MESONs-1) - offset_space;
        ++space;
      }
    }

    map[pion]       = (MESONs-1);
    map[pion_prime] = (MESONs-1) - offset_time;
    map[a0]         = 0;
    map[a0_prime]   =              offset_time;
    map[a1]         = map[a1_x];
    map[b1]         = map[b1_x];
    map[rho]        = map[rho_x];
    map[rho_prime]  = map[rho_x_prime];
  }

  char *meson_names[16];
  meson_names[0] = w_spect_output->meson_name00;
  meson_names[1] = w_spect_output->meson_name01;
  meson_names[2] = w_spect_output->meson_name02;
  meson_names[3] = w_spect_output->meson_name03;
  meson_names[4] = w_spect_output->meson_name04;
  meson_names[5] = w_spect_output->meson_name05;
  meson_names[6] = w_spect_output->meson_name06;
  meson_names[7] = w_spect_output->meson_name07;
  meson_names[8] = w_spect_output->meson_name08;
  meson_names[9] = w_spect_output->meson_name09;
  meson_names[10] = w_spect_output->meson_name10;
  meson_names[11] = w_spect_output->meson_name11;
  meson_names[12] = w_spect_output->meson_name12;
  meson_names[13] = w_spect_output->meson_name13;
  meson_names[14] = w_spect_output->meson_name14;
  meson_names[15] = w_spect_output->meson_name15;

  // -mcneile
  // this routine was only correct for the a0, a0_prime
  // I can not see a simple way to fix it in the previous style
  // The meson names were obtained by empirically matching 
  // against the regression tests.

#if 0
  char *meson_names[] = 
  {
#if 1
    "Scalar.dat" , 
    "Vector_x.dat" , 
    "Vector_y.dat"    ,
    "Tensor_xy.dat"     ,
    "Vector_z.dat"    ,
    "Tensor_xz.dat",
    "Tensor_yz.dat",
    "PV_t.dat",
    "Vector_t.dat",
    "Tensor_xt.dat",
    "Tensor_yt.dat",
    "PV_z.dat" ,
    "Tensor_zt.dat" , 
    "PV_y.dat" , 
    "PV_x.dat" , 
    "PS.dat" 
#else
    "a0.dat" , 
    "rho_x.dat" , 
    "rho_y.dat"    ,
    "b1_z.dat"     ,
    "rho_z.dat"    ,
    "b1_y.dat",
    "b1_x.dat",
    "pion_prime.dat",
    "a0_prime.dat",
    "rho_x_prime.dat",
    "rho_y_prime.dat",
    "a1_z.dat" ,
    "rho_z_prime.dat" , 
    "a1_y.dat" , 
    "a1_x.dat" , 
    "pion.dat" 
#endif

  };
#endif

  for(int imap= 0 ; imap < NUM_MAP ; ++imap)
    {  map[imap] = imap ; }


  //  --mcneile
  //  I am not sure what the lines below were meant to do 
  //  The attempt to average over x,y,z did not work

#ifdef GGGGGGGGGGGGGG
  // average over x,y,z for vector mesons
  //        Complex[d_num_mom][global_slices][MESONs]
  //-------------------------------------------------------------------------
  {
    int size_in_complex = d_size / COMPLEXs;
    int x1 = map[a1_x];       int y1 = map[a1_y];    int z1 = map[a1_z];
    int x2 = map[b1_x];       int y2 = map[b1_y];    int z2 = map[b1_z];
    int x3 = map[rho_x];      int y3 = map[rho_y];   int z3 = map[rho_z];
    int x4 = map[rho_x_prime];
    int y4 = map[rho_y_prime];
    int z4 = map[rho_z_prime];

    for (int i = 0; i < size_in_complex; i += MESONs) {
      d_data_p[i+x1] += d_data_p[i+y1] += d_data_p[i+z1]; d_data_p[i+x1] /= 3;
      d_data_p[i+x2] += d_data_p[i+y2] += d_data_p[i+z2]; d_data_p[i+x2] /= 3;
      d_data_p[i+x3] += d_data_p[i+y3] += d_data_p[i+z3]; d_data_p[i+x3] /= 3;
      d_data_p[i+x4] += d_data_p[i+y4] += d_data_p[i+z4]; d_data_p[i+x4] /= 3;
    }
  }
#endif 

  // Print out mesons
  //------------------------------------------------------------------
  FILE *fp;
  for (int meson = 0; meson < 16  ; ++meson) 
    {
      
      if ( !meson_names[meson] || !(fp = Fopen(meson_names[meson], "a")) ) {
	ERR.FileA(d_class_name,fname, meson_names[meson]);
      }

      printf("Writing data to file: %s\n",meson_names[meson]);
      for (int wall = 0; wall < d_glb_walls; ++wall) {
	for (int mom = 0; mom < d_num_mom; ++mom) {
	  int offset = mom * d_glb_walls * MESONs + map[meson];
	  Complex tmp = d_data_p[offset + MESONs*
				((d_whr.glbCoord() + wall)%d_glb_walls)];
	  // only print out the real parts
	  Fprintf(fp, "%d %d %d %.10e\n", 
		  AlgWspect::GetCounter(), wall, mom, tmp.real());
	}
      }
      Fclose(fp);

    }



}


//---------------------------------------------------------------------------
// void WspectMesons::print_mp(const CommonArg *common_arg) const
//--------------------------------------------------------------------------- 
void
WspectMesons::print_mp(char *filename) const
{

#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  char *fname = "print_mp";
  VRB.Func(d_class_name, fname);

  // Print out < \Delta J^5 \bar q \gamma^5 q> correlator
  // (mid-point sink)
  //------------------------------------------------------------------
  FILE *fp;

  if ( filename != 0 ) {
    if ( !(fp = Fopen(filename, "a")))
      ERR.FileA(d_class_name,fname, filename);
  
    for (int wall = 0; wall < d_glb_walls; ++wall) {
      for (int mom = 0; mom < d_num_mom; ++mom) {
        int offset = mom * d_glb_walls * MESONs + MESONs -1;
        Complex tmp = d_data_p[offset + MESONs*
                            ((d_whr.glbCoord() + wall)%d_glb_walls)];
        // only print out the real part
        Fprintf(fp, "%d %d %d %e\n",
                AlgWspect::GetCounter(), wall, mom, tmp.real());
      }
    }
    Fclose(fp);
  }
}

//---------------------------------------------------------------------------
// void WspectMesons::dumpData() const 
//--------------------------------------------------------------------------- 
void
WspectMesons::dumpData(char *filename) const {

#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif

  FILE *fp;
  
  if (filename && (fp = Fopen(filename, "a"))) {
    for (int i = 0; i < d_size/COMPLEXs; ++i) {
      Fprintf(fp, "%.10e %.10e\n", d_data_p[i].real(), d_data_p[i].imag());
    }
    Fclose(fp);    
  } else {
    ERR.FileA(d_class_name, "dumpData", filename);
  }
}





























//---------------------------------------------------------------------------
//  ALL THE FOLLOWING CODE ARE ONLY FOR DEBUGGING PURPOSE.
//  WE COMMENT THEM ALL OUT NOW.
//

//---------------------------------------------------------------------------
//  math lib for the calculation of non-zero spatial momenta
//---------------------------------------------------------------------------
// Float COS(Float x) 
//         calculate cos(2 PI x)
// Float SIN(Float x)
//         calculate sin(2 PI x)
// Note:
//  The Tartan library messes up the sign of cos and sin. It is correct
//  iff the angle falls [0, PI]. 
//  It is more accurate to comparing angle/(2 PI) against 0.5.
//---------------------------------------------------------------------------


//#define TEST_COS_SIN
#ifdef  TEST_COS_SIN
static void testCosSin()
{
  Float f, g;  
  f = MY_TWO_PI * 0.101;     
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
  f += MY_TWO_PI * 0.2;    
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)f, 
	 (IFloat)(COS(f)), (IFloat)f, (IFloat)(SIN(f)));
  g = -f;  
  printf("COS(%g) = %g SIN(%g) = %g\n", (IFloat)g, 
	 (IFloat)(COS(g)), (IFloat)g, (IFloat)(SIN(g)));
}
#endif   // #ifdef TEST_COS_SIN

CPS_END_NAMESPACE
