/*! \file

  $Id: w_ferm_vec.C,v 1.6 2004-08-17 03:33:11 chulwoo Exp $
*/

#include<config.h>
#include <alg/w_all.h>
#include <util/gjp.h>              // GJP
#include <util/error.h>            // ERR
#include <util/verbose.h>          // VRB
#include <util/lattice.h>          // Lattice::FixGaugePtr()
#include <util/qcdio.h>
#include <comms/glb.h>               // glb_sum(...)
#include <comms/scu.h>              //getMinusData, getPlusData

//Warning: do not use math64.h, seems to use some rsgisters that
//is in conflict with optimized code.
#include <math.h>

CPS_START_NAMESPACE

#define DEBUG_FERMIONVECTOR
//---------------------------------------------------------------------------
// static data members
//---------------------------------------------------------------------------
char * FermionVector::d_class_name = "FermionVector";
Vector FermionVector::m_tmp1;

//---------------------------------------------------------------------------
// FermionVector::FermionVector()
//---------------------------------------------------------------------------
FermionVector::FermionVector()
  : d_size(GJP.VolNodeSites() * (COLORs*DIRACs*COMPLEXs))
{
  VRB.Func(d_class_name, ctor_str);

  // allocate memory for d_fermion_p
  d_fermion_p = (IFloat *)smalloc(d_size * sizeof(IFloat));
  if (!d_fermion_p) 
    ERR.Pointer(d_class_name, ctor_str, empty_str);
  VRB.Smalloc(d_class_name, ctor_str, empty_str, 
	      d_fermion_p, d_size * sizeof(IFloat));

  // removed allocation of d_proj_fermion_p
  // last present in version v3, T&X March 29, 2000
  
#ifdef DEBUG_FERMIONVECTOR
  printf("allocate source vector - %i IFloat at %x\n",d_size,d_fermion_p);
#endif 

 }


//---------------------------------------------------------------------------
// FermionVector::~FermionVector()
//---------------------------------------------------------------------------
FermionVector::~FermionVector()
{
  VRB.Func(d_class_name, dtor_str);
  VRB.Sfree(d_class_name, dtor_str, empty_str, d_fermion_p);
  sfree(d_fermion_p);

  // removed de-allocation of d_proj_fermion_p
  // last present in version v3, T&X March 29, 2000
}


//---------------------------------------------------------------------------
// void FermionVector::zeroOut() const
//---------------------------------------------------------------------------
void
FermionVector::zeroOut() const
{
  for (int i = 0; i < d_size; ++i) 
    d_fermion_p[i] = 0.0;
}

//---------------------------------------------------------------------------
// FermionVector::print(const char *file)
//---------------------------------------------------------------------------
void FermionVector::print(char *file) const
{
  FILE *fp;
  if (!file || !(fp = Fopen(file, "a"))) 
    ERR.FileA(d_class_name,"print", file);


  int glb[LORENTZs];               // a lattice site on the whole machine
  int s;                           // spinor index = color * dirac * complex 
  IFloat spinor[SPINORs];           // a spinor  

  for (glb[0] = 0; glb[0] < glb_sites[0]; glb[0]++)   {
    for (glb[1] = 0; glb[1] < glb_sites[1]; glb[1]++)   {
      for (glb[2] = 0; glb[2] < glb_sites[2]; glb[2]++)   {
	for (glb[3] = 0; glb[3] < glb_sites[3]; glb[3]++)   {

	  printSite(fp, glb);

	  int lcl[LORENTZs];
	  int is_on_node = glb2lcl(lcl, glb);	  
	  int site_offset = siteOffset(lcl);
	  
	  for (s = 0; s < SPINORs; ++s) {
	    spinor[s] = is_on_node ? d_fermion_p[s+SPINORs*site_offset] : 0.0;
	    glb_sum((Float *)(spinor+s));
	  }

	  printSpinor(fp, spinor);
	}
      }
    }
  }

  Fclose(fp);  
}

//---------------------------------
//For debugging
//void FermionVector::printWaveFunc(char *filename)
//---------------------------------
void FermionVector::printWaveFunc(char *file) const{

#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
    FILE *fp;
  if (!file || !(fp = Fopen(file, "a"))) 
    ERR.FileA(d_class_name,"print", file);
  

  int glb[LORENTZs];               // a lattice site on the whole machine
  int s;                           // spinor index = color * dirac * complex 
  IFloat spinor[SPINORs];           // a spinor  
  
  for (glb[0] = 0; glb[0] < glb_sites[0]; glb[0]++)   {
    for (glb[1] = 0; glb[1] < glb_sites[1]; glb[1]++)   {
      for (glb[2] = 0; glb[2] < glb_sites[2]; glb[2]++)   {

	  glb[3]=0; //supose prop_dir==3 and source plane==0

	  int lcl[LORENTZs];
	  int is_on_node = glb2lcl(lcl, glb);	  
	  int site_offset = siteOffset(lcl);
	  
	  for (s = 0; s < SPINORs; ++s) {
	    spinor[s] = is_on_node ? d_fermion_p[s+SPINORs*site_offset] : 0.0;
	    glb_sum((Float *)(spinor+s));
	  }
	  
	  //calculate average amplitude
	  
	  IFloat sum;
	  IFloat radius=0;
	  IFloat time;
	  IFloat abs=0; 
	  s=0;
	  sum=0.0;
	  for(int D1=0;D1<DIRACs;D1++){
	    for(int C1=0; C1<COLORs; C1++){
	      //Note: Do not use rcomplex::abs which uses "math64.h"!
#ifdef _TARTAN
              abs=(sqrt(spinor[s]*spinor[s]+spinor[s+1]*spinor[s+1]));
#else
              abs=sqrt(spinor[s]*spinor[s]+spinor[s+1]*spinor[s+1]);
#endif   
	      sum=sum+abs;
	      s=s+COMPLEXs;
	    }
	  }

	  sum=sum/(DIRACs*COLORs);
#ifdef _TARTAN
	  radius=(sqrt(glb[0]*glb[0]+glb[1]*glb[1]+glb[2]*glb[2]));
#else
	  radius=sqrt(glb[0]*glb[0]+glb[1]*glb[1]+glb[2]*glb[2]);
#endif 
	  time=glb[3];
	  
	  Fprintf(fp,"%g %g %g\n", radius, time, sum);
	
      }
    }
  }

  Fclose(fp);  
}

//---------------------------------------------------------------------------
// void FermionVector::
// setPointSrc(int color, int spin, int src_global_coord[]) const
//---------------------------------------------------------------------------
void FermionVector::setPointSrc(int C, int D, const int *lcl_src) const
{
  char *fname = "setPointSrc";
  VRB.Func(d_class_name, fname);

  // zero the vector
  //-------------------------------------------------------------------------
  zeroOut();

  // set the source
  //-------------------------------------------------------------------------
  if (lcl_src) {
    d_fermion_p[COMPLEXs*(C+COLORs*(D+DIRACs*siteOffset(lcl_src)))] = 1.;
  }
}



//---------------------------------------------------------------------------
// removed overloaded setPointSrc for all colours and spins simultaneously
// last present in version v3. T&X; March 29 2000
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// void FermionVector::setSmearedSrc(...)
//---------------------------------------------------------------------------
void 
FermionVector::setSmearedSrc(int C, int D, const IFloat *src_mat,
			     int dir,
			     const int *lclMin,
			     const int *lclMax) const
{
  char *fname = "setSmearedSrc()";
  VRB.Func(d_class_name, fname);

  // Zero out the entire vector
  //-------------------------------------------------------------------------
  zeroOut();

  // Return if the wall does not go across this node.
  //-------------------------------------------------------------------------
  if (!src_mat)
    return;

  // Set the source  -- walk in 3D -- Note: lclMin[L] = lclMax[L] = lclW
  //-------------------------------------------------------------------------
  int lcl[LORENTZs];                             // local lattice site
  for (lcl[0] = lclMin[0]; lcl[0] <= lclMax[0]; lcl[0]++) {
    for (lcl[1] = lclMin[1]; lcl[1] <= lclMax[1]; lcl[1]++) {	
      for (lcl[2] = lclMin[2]; lcl[2] <= lclMax[2]; lcl[2]++) {
	for (lcl[3] = lclMin[3]; lcl[3] <= lclMax[3]; lcl[3]++) {
	  for (int c = 0; c < COLORs; c++) {
	    int frm_offset = COMPLEXs*(c+COLORs*(D+DIRACs*siteOffset(lcl)));
	    int gau_offset = COMPLEXs*(c+COLORs*(C+COLORs*
						 siteOffset(lcl, dir)));
	    d_fermion_p[frm_offset++] =  src_mat[gau_offset++];    // real
	    d_fermion_p[frm_offset]   = -src_mat[gau_offset];      // imag
	  }
	}
      }
    }
  }
}


//---------------------------------------------------------------------------
// removed overloaded function setSmearedSrc ((const IFloat *src_mat, 
//			     int dir,			     
//			     const int *lclMin,
//			     const int *lclMax) const
// last present in version v3. Thomas and Xiadong, March 29 2000
//---------------------------------------------------------------------------


// -----------------------------------------------------------------------
// added by Thomas and Xiadong: Jacobi-smeared source (1-epsi/n \Delta^2)
// not yet implemented
//---------------------------------------------------------------------------
void FermionVector::setJacobiSrc(const int *lcl_src, const Lattice &lat, Float epsi, int n_iter) const
{
  char *fname = "setJacobiSrc";
  VRB.Func(d_class_name, fname);
}


// removed symmetric derivative acting on source vector
// Thomas and Xiadong, March 29 2000.

// added by Thomas and Xiaodong:
// get Fermionvector at site defined by *site
// and with Dirac index defined by dirac
// result is a 3-vector with complex entries = 6 IFloats
// assumes that off-nodes are NEXT-NODES ?

const Vector *
FermionVector::GetFermion(const int *site, int dirac) const
{
//char *fname = "GetFermion()";
//VRB.Func(cname,fname);

//VRB.Debug(cname, fname, "fermion %3i %3i %3i %3i ; %i\n",
//  site[0], site[1], site[2], site[3], dirac) ;

  // offset out-of-range coordinates site[] into on_node_site[]
  // in order to locate the link
  //------------------------------------------------------------------------
  int on_node_site[4];
  int on_node = 1;
  const Vector *on_node_fermion;

  int node_sites[LORENTZs];
  node_sites[0]=GJP.XnodeSites();
  node_sites[1]=GJP.YnodeSites();
  node_sites[2]=GJP.ZnodeSites();
  node_sites[3]=GJP.TnodeSites();
  

  {
    for (int i = 0; i < 4; ++i) {
      on_node_site[i] = site[i] ;
      while (on_node_site[i] < 0) { // map negative values into smallest positive
        on_node_site[i] += node_sites[i] ;
      }
      on_node_site[i] %= node_sites[i];
      if (on_node_site[i] != site[i]) {  // 0%4=0, 1%4=1, ..., 3%4=3, 4%4=0 !=4 --> off-node
        on_node = 0;
      }
    }
    int FsiteOffset = (on_node_site[3]*node_sites[2]*node_sites[1]*node_sites[0]*DIRACs +
		       on_node_site[2]*node_sites[1]*node_sites[0]*DIRACs +
		       on_node_site[1]*node_sites[0]*DIRACs +
		       on_node_site[0]*DIRACs + 
		       dirac) *COLORs*COMPLEXs;
    
    // should be the same ! siteOffset is implemented in w_ginfo.C
    // int FsiteOffset = (siteOffset(on_node_site)*DIRACs + dirac) *COLORs*COMPLEXs;
    
    on_node_fermion = (Vector *)(d_fermion_p + FsiteOffset) ;
  }


#ifndef PARALLEL
//VRB.FuncEnd(cname, fname) ;
  return on_node_fermion;
#endif

  // send to the destination node if the site is off-node
  //------------------------------------------------------------------------
  if (on_node) {
//  VRB.FuncEnd(cname, fname) ;
    return on_node_fermion;
  } else {
    Vector send = *on_node_fermion;
    Vector &recv = m_tmp1 ;
    for (int i = 0; i < 4; ++i) {
      while (site[i] != on_node_site[i]) {
        if (site[i] < 0) {
          getMinusData((IFloat *)&recv, (IFloat *)&send, sizeof(recv)/sizeof(IFloat), i);
          on_node_site[i] -= node_sites[i];
        } else {
          getPlusData((IFloat *)&recv, (IFloat *)&send, sizeof(recv)/sizeof(IFloat), i);
          on_node_site[i] += node_sites[i];
        }
        send = recv;
      }
    }
//  VRB.FuncEnd(cname, fname) ;
    return &recv ;
  }
}


//---------------------------------------------------------------------------
// removed projectSource(int Cpick, int Dpick) const
// last present in version v3, T&X March 29, 2000
//---------------------------------------------------------------------------




//---------------------------------------------------------------------------
// void FermionVector::setSource(...)
// sets source based on a colour x colour source matrix src_mat(x)
// it projects out the column C 
// at the moment the spin part is assumed to be diagonal
// but this could be generalised by passing a bigger src_mat
// src_mat is defined for all sites, but it can be "boxed-out"
// if (lclMin,lclMax) is passed to this function
// Thomas & Xiaodong, April 3rd 2000.
//---------------------------------------------------------------------------
void 
FermionVector::setSource(int C, int D, IFloat rs_fac, const IFloat *src_mat,
			     int dir,
			     const int *lclMin,
			     const int *lclMax) const
{
  char *fname = "setSource()";
  VRB.Func(d_class_name, fname);

  // Zero out the entire vector
  //-------------------------------------------------------------------------
  zeroOut();

  // Return if the wall does not go across this node.
  //-------------------------------------------------------------------------
  if (!src_mat)
    return;
  //above code doesn't work at all! source_slice is always allocated in WspectQuark
  //waste of memory!
  // Set the source  -- walk in 3D -- Note: lclMin[L] = lclMax[L] = lclW
  //-------------------------------------------------------------------------
  int lcl[LORENTZs];                             // local lattice site
  for (lcl[0] = lclMin[0]; lcl[0] <= lclMax[0]; lcl[0]++) {
    for (lcl[1] = lclMin[1]; lcl[1] <= lclMax[1]; lcl[1]++) {	
      for (lcl[2] = lclMin[2]; lcl[2] <= lclMax[2]; lcl[2]++) {
	for (lcl[3] = lclMin[3]; lcl[3] <= lclMax[3]; lcl[3]++) {
	  for (int c = 0; c < COLORs; c++) {
	    int frm_offset = COMPLEXs*(c+COLORs*(D+DIRACs*siteOffset(lcl)));
	    int gau_offset = COMPLEXs*(C+COLORs*(c+COLORs*siteOffset(lcl, dir)));
	    d_fermion_p[frm_offset++] =  rs_fac * src_mat[gau_offset++];    // real
	    d_fermion_p[frm_offset]   =  rs_fac * src_mat[gau_offset];      // imag
	  }
	}
      }
    }
  }
}





CPS_END_NAMESPACE
