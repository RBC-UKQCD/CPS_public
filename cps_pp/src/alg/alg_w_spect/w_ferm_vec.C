/*! \file

  $Id: w_ferm_vec.C,v 1.11 2008-02-08 18:35:05 chulwoo Exp $
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
#include <util/time_cps.h>             //dclock(),print_flops()
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
int FermionVector::allocated = 0;

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

  zeroOut();
  // removed allocation of d_proj_fermion_p
  // last present in version v3, T&X March 29, 2000
  
#ifdef DEBUG_FERMIONVECTOR
  printf("allocate source vector - %i IFloat at %x\n",d_size,d_fermion_p);
#endif 
  allocated = 1;
 }

//---------------------------------------------------------------------------
// FermionVector::FermionVector(Float *)
// This is a dangerous hack. Since the array has to have the correct size.
// --- Meifeng Lin, 01/31/06
//---------------------------------------------------------------------------
FermionVector::FermionVector(Float *array)
  : d_size(GJP.VolNodeSites() * (COLORs*DIRACs*COMPLEXs))
{
  VRB.Func(d_class_name, ctor_str);
  if (array == NULL) ERR.General(d_class_name, ctor_str, "Array is not intialized\n");
  d_fermion_p = array;
  allocated = 0;
}



FermionVector& FermionVector::operator=(const FermionVector& p){
  for ( int n = 0; n < d_size; n++ )
    d_fermion_p[n] = p.d_fermion_p[n];
  return *this;     
}
  

FermionVector& FermionVector::operator+=(const FermionVector& f1){
  for ( int n = 0; n < d_size; n++ )
    d_fermion_p[n] += f1.d_fermion_p[n]; 
  return *this;
}

//---------------------------------------------------------------------------
// FermionVector::~FermionVector()
//---------------------------------------------------------------------------
FermionVector::~FermionVector()
{
  VRB.Func(d_class_name, dtor_str);
  
  if(allocated){
    VRB.Sfree(d_class_name, dtor_str, empty_str, d_fermion_p);
    sfree(d_fermion_p);
  }

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
	  radius=sqrt((Float)(glb[0]*glb[0]+glb[1]*glb[1]+glb[2]*glb[2]));
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

//----------------------------------------------------------------------------
// Coulomb gauge fix fermion solution vector in hyperplanes of direction dir.
// Gauge fixing matrices must exist before this routine is used.
// --- Meifeng Lin 01/31/06
//---------------------------------------------------------------------------- 
void FermionVector::gaugeFixSink(Lattice &lat, int dir) const
{
  char *fname = "gaugeFixSink()";
  VRB.Func(d_class_name, fname);
  
  Matrix **gm = lat.FixGaugePtr();
  //-------------------------------------------------------
  //If gauge fixing matrices do not exist, 
  //nothing has to be done.
  //-------------------------------------------------------
  if(gm==0){
    ERR.General(d_class_name,fname,"No gauge fixing matrices were found\n");
    return;
  }
  else{
    printf("Gauge Fixing Sink...\n");

  }

  IFloat *fv = d_fermion_p;

  int error = 0;
  switch(dir){
  case 0:
    if ( lat.FixGaugeKind() != FIX_GAUGE_COULOMB_X )
      error = 1;break;
  case 1:
    if ( lat.FixGaugeKind() != FIX_GAUGE_COULOMB_Y )
      error = 1;break;
  case 2:
    if ( lat.FixGaugeKind() != FIX_GAUGE_COULOMB_Z )
      error = 1;break;
  case 3:
    if ( lat.FixGaugeKind() != FIX_GAUGE_COULOMB_T )
      error = 1;break;
  }
  if (error)
    ERR.General(d_class_name,fname,"Mismatch gauge fixing kind\n");

  int lcl[LORENTZs];
  
  //-------------------------------------------------------------------------
  // Indices for all the directions
  // Here use "t" to indicate the propagation direction, but not necessarily 
  // the generic t direction. Other 3 directions are arbitrarily labelled x,y,z
  //-------------------------------------------------------------------------
  int t = dir;
  int z = (dir + 1) % LORENTZs;
  int y = (dir + 2) % LORENTZs;
  int x = (dir + 3) % LORENTZs;
  
  Vector tmp_vec;
  
  for (lcl[t]=0; lcl[t] < lcl_sites[t]; lcl[t]++) {
    
    Matrix *pM = gm[lcl[t]];
    
    if(pM != NULL ){
      for (lcl[z]  = 0; lcl[z] < lcl_sites[z]; lcl[z]++)
	for (lcl[y] = 0; lcl[y] < lcl_sites[y]; lcl[y]++)
	  for (lcl[x] = 0; lcl[x] < lcl_sites[x]; lcl[x]++)
	    {
	      // the matrix offset
	      int j =  siteOffset(lcl,dir);
	      for (int spin = 0; spin < DIRACs; spin++)
		{
		  // the vector offset
		  int i= 2 * COLORs * ( spin + 4 * siteOffset(lcl));
		  tmp_vec.CopyVec((Vector*)&fv[i], 6);
		  uDotXEqual((IFloat*)&fv[i],(const IFloat*)&pM[j],
			     (const IFloat*)&tmp_vec);
		}
	    }
    }
    else{
         ERR.General(d_class_name,fname,"Not gauge fixing matrices were found at slice = %d\n",lcl[t]); 
         Float *null_p = 0;
         *null_p = 1.;
    }
  }


}

//----------------------------------------------------------
//sum the fermion vector over all the hyperplanes bounded by 
//box_b[] and box_e[] 
//in the propagation direction.
//In fact, the result should be only of the length Nt*24.
//In order to use the existing meson correlator code, we preserve the
//structure of a full length fermion vector.
//All the sites on the same timeslice has the same spin-color vector, 
//as depicted in the following.
//|---------------|
//| sum at t = 0  | site (0,0,0,0) 
//|---------------|
//| sum at t = 0  | site (1,0,0,0)
//|---------------|
//|               | site (2,0,0,0)
//| ....          |
//|---------------|
//| sum at t = 1  | site (0,0,0,1)
//|---------------|
//| sum at t =1   | site (1,0,0,1)
//|---------------| 
//| ....          |
//|---------------|
//| sum at t=Nt-1 | site (0,0,0,Nt-1)
//|---------------|
//| sum at t=Nt-1 |
//|---------------|
//| ....          |
//|---------------|
//| sum at t=Nt-1 | site (Nx-1,Ny-1,Nz-1,Nt-1)
//|---------------|
//
//-----------------------------------------------------------------------------------------

void FermionVector::sumOverHyperPlane(int dir, int box_b[], int box_e[]){
  char *fname = "sumOverHyperPlane(int,int,int)";
  VRB.Func(d_class_name, fname);
#ifdef TIMING
  Float time = dclock();
#endif
  for( int lcl = 0; lcl < GJP.NodeSites(dir); lcl++){
    int hplane = lcl2glb_offset[dir] + lcl;
    sumOverHyperPlane(dir,hplane, box_b, box_e);      
  }
#ifdef TIMING
  time = dclock() - time;
  print_flops(d_class_name,fname,0,time);
#endif
}


#undef TIMING


//-----------------------------------------------------------------------------------------
//sum the fermion vector over a hyperplane with fixed t bounded by 
//box_b[] and box_e[] 
//in the propagation direction
//----------------------------------------------------------
void FermionVector::sumOverHyperPlane(int dir,int hplane, int box_begin[], int box_end[])
{
  char *fname = "sumOverHyperPlane(int,int,int,int)";
  VRB.Func(d_class_name,fname);
#ifdef TIMING
  Float time = dclock();
#endif
  //------------------------------------------
  //Do nothing is the timeslice is not 
  //on this node
  //------------------------------------------
  if ( hplane < lcl2glb_offset[dir] ||
	hplane > lcl2glb_offset[dir] + lcl_sites[dir] - 1 )
	return ;
 
  //--------------------------------------------------
  //Various index
  //--------------------------------------------------
  int lcl[LORENTZs],lclMin[LORENTZs],lclMax[LORENTZs];
  
  int glb[LORENTZs];
  
  int box_b[LORENTZs], box_e[LORENTZs];
  for(int n = 0; n < LORENTZs; n++){
    box_b[n] = box_begin[n];
    box_e[n] = box_end[n];
  }

  //the box doesn't extend in the propagation direction (except for volume source ,which I don't consider here)
  box_b[dir] = box_e[dir] = hplane;
  
  //fermion vector for one site
  int vec_len = COMPLEXs * COLORs * DIRACs;
  
  //------------------------------------
  //Temporary data containers
  //------------------------------------
  Vector *pv; //pointer to the value at each site
  Vector *sum = (Vector *)smalloc(d_class_name,fname,"sum",vec_len*sizeof(Float));

  for(int i = 0; i < LORENTZs; i++) {
    lcl[i] = 0;
    lclMin[i] = 0;
    lclMax[i] = lcl_sites[i]-1;
  }
  // no looping in the propagation direction
  lcl[dir] = lclMin[dir] = lclMax[dir] = hplane - lcl2glb_offset[dir]; 
 
  pv = (Vector *)(d_fermion_p + vec_len * siteOffset(lcl));
 
  // clears the summation container
  sum->VecZero(vec_len);
  for (lcl[0] = lclMin[0]; lcl[0] <= lclMax[0]; lcl[0]++) {
    for (lcl[1] = lclMin[1]; lcl[1] <= lclMax[1]; lcl[1]++) {	
      for (lcl[2] = lclMin[2]; lcl[2] <= lclMax[2]; lcl[2]++) {
	for (lcl[3] = lclMin[3]; lcl[3] <= lclMax[3]; lcl[3]++) {
	  //map local coordinates to global
	  lcl2glb(lcl,glb);
	
	  //-----------------------------------------
	  //flag to determine if the local sites 
	  //hit the hyperplane bounded by box_b and box_e
	  //-------------------------------------------
	  int hit = 1; 
	  for (int i = 0; i < LORENTZs; i++){
	    if ( i != dir) {
	      //check if the box is normal bound box
	      switch(box_b[i]<box_e[i]){
	      case 1:
		if ( glb[i] <= box_e[i] && glb[i] >= box_b[i] ) hit &= 1;
		else hit &= 0; break;
		//around the world box
	      case 0:
		if ( glb[i] <= box_e[i] || glb[i] >= box_b[i] ) hit &= 1;
		else hit &= 0; break;
		
	      }
	    }
	  }

	  pv = (Vector *)(d_fermion_p + vec_len * siteOffset(lcl));
	  //if fall in the box, do the sum
	  if ( hit ){
	    sum->VecAddEquVec(pv,vec_len);	  
       	  }
	  
	  //clean the vector to store the final results
	  pv->VecZero(vec_len);
	}// end lcl[3]
      }// end lcl[2]
    }//end lcl[1]
  }//end lcl[0]
	  
  //----------------------------------------------------------
  //Copy the summation into the timeslice
  //All the sites on the same timeslice have the same result.
  //----------------------------------------------------------
  slice_sum((Float *)sum,vec_len,dir);
  for (lcl[0] = lclMin[0]; lcl[0] <= lclMax[0]; lcl[0]++) 
    for (lcl[1] = lclMin[1]; lcl[1] <= lclMax[1]; lcl[1]++) 	
      for (lcl[2] = lclMin[2]; lcl[2] <= lclMax[2]; lcl[2]++) 
	for (lcl[3] = lclMin[3]; lcl[3] <= lclMax[3]; lcl[3]++) 
	  
	  {
	    pv = (Vector *)(d_fermion_p + vec_len * siteOffset(lcl));
	    
	    pv->CopyVec(sum,vec_len);
	  }
  

  sfree(sum);
#ifdef TIMING
  time = dclock() - time;
  print_flops(d_class_name,fname,0,time);
#endif

}

//#define TIMING

//----------------------------------------------------------
//sum the fermion vector over all the hyperplanes bounded by 
//box_b[] and box_e[] 
//in the propagation direction
//The result is stored in the first site of the timeslice. For example, if the sink is at
//(0,0,0) to (3,3,3), then the resulting fermion vector would look like
//|---------------|
//| sum at t = 0  | site (0,0,0,0) 
//|---------------|
//| 0             | site (1,0,0,0)
//|---------------|
//| 0             | site (2,0,0,0)
//| ....          |
//|---------------|
//| sum at t = 1  | site (0,0,0,1)
//|---------------|
//| 0             | site (1,0,0,1)
//|---------------| 
//| ....          |
//|---------------|
//| sum at t=Nt-1 | site (0,0,0,Nt-1)
//|---------------|
//| 0             |
//|---------------|
//| ....          |
//|---------------|
//| 0             | site (Nx-1,Ny-1,Nz-1,Nt-1)
//|---------------|
//
//-----------------------------------------------------------------------------------------

void FermionVector::sumOverHyperPlaneStride(int dir, int box_b[], int box_e[]){
  char *fname = "sumOverHyperPlane(int,int,int)";
  VRB.Func(d_class_name, fname);
#ifdef TIMING
  Float time = dclock();
#endif
  for( int lcl = 0; lcl < GJP.NodeSites(dir); lcl++){
    int hplane = lcl2glb_offset[dir] + lcl;
    sumOverHyperPlaneStride(dir,hplane, box_b, box_e);      
  }
#ifdef TIMING
  time = dclock() - time;
  print_flops(d_class_name,fname,0,time);
#endif
}


#undef TIMING


//-----------------------------------------------------------------------------------------
//sum the fermion vector over a hyperplane with fixed t bounded by 
//box_b[] and box_e[] 
//in the propagation direction
//----------------------------------------------------------
void FermionVector::sumOverHyperPlaneStride(int dir,int hplane, int box_begin[], int box_end[])
{
  char *fname = "sumOverHyperPlane(int,int,int,int)";
  VRB.Func(d_class_name,fname);
#ifdef TIMING
  Float time = dclock();
#endif
  //------------------------------------------
  //Do nothing is the timeslice is not 
  //on this node
  //------------------------------------------
  if ( hplane < lcl2glb_offset[dir] ||
	hplane > lcl2glb_offset[dir] + lcl_sites[dir] - 1 )
	return ;
 
  //--------------------------------------------------
  //Various index
  //--------------------------------------------------
  int lcl[LORENTZs],lclMin[LORENTZs],lclMax[LORENTZs];
  
  int glb[LORENTZs];
  
  int box_b[LORENTZs], box_e[LORENTZs];
  for(int n = 0; n < LORENTZs; n++){
    box_b[n] = box_begin[n];
    box_e[n] = box_end[n];
  }

  //the box doesn't extend in the propagation direction (except for volume source ,which I don't consider here)
  box_b[dir] = box_e[dir] = hplane;
  
  //fermion vector for one site
  int vec_len = COMPLEXs * COLORs * DIRACs;
  
  //------------------------------------
  //Temporary data containers
  //------------------------------------
  Vector *pv; //pointer to the value at each site
  Vector *sum = (Vector *)smalloc(d_class_name,fname,"sum",vec_len*sizeof(Float));

  for(int i = 0; i < LORENTZs; i++) {
    lcl[i] = 0;
    lclMin[i] = 0;
    lclMax[i] = lcl_sites[i]-1;
  }
  // no looping in the propagation direction
  lcl[dir] = lclMin[dir] = lclMax[dir] = hplane - lcl2glb_offset[dir]; 
 
  pv = (Vector *)(d_fermion_p + vec_len * siteOffset(lcl));
 
  // clears the summation container
  sum->VecZero(vec_len);
  for (lcl[0] = lclMin[0]; lcl[0] <= lclMax[0]; lcl[0]++) {
    for (lcl[1] = lclMin[1]; lcl[1] <= lclMax[1]; lcl[1]++) {	
      for (lcl[2] = lclMin[2]; lcl[2] <= lclMax[2]; lcl[2]++) {
	for (lcl[3] = lclMin[3]; lcl[3] <= lclMax[3]; lcl[3]++) {
	  //map local coordinates to global
	  lcl2glb(lcl,glb);
	
	  //-----------------------------------------
	  //flag to determine if the local sites 
	  //hit the hyperplane bounded by box_b and box_e
	  //-------------------------------------------
	  int hit = 1; 
	  for (int i = 0; i < LORENTZs; i++){
	    if ( i != dir) {
	      //check if the box is normal bound box
	      switch(box_b[i]<box_e[i]){
	      case 1:
		if ( glb[i] <= box_e[i] && glb[i] >= box_b[i] ) hit &= 1;
		else hit &= 0; break;
		//around the world box
	      case 0:
		if ( glb[i] <= box_e[i] || glb[i] >= box_b[i] ) hit &= 1;
		else hit &= 0; break;
		
	      }
	    }
	  }

	  pv = (Vector *)(d_fermion_p + vec_len * siteOffset(lcl));
	  //if fall in the box, do the sum
	  if ( hit ){
	    sum->VecAddEquVec(pv,vec_len);	  
       	  }
	  
	  //clean the vector to store the final results
	  pv->VecZero(vec_len);
	}// end lcl[3]
      }// end lcl[2]
    }//end lcl[1]
  }//end lcl[0]
	  
  
 
  slice_sum((Float *)sum,vec_len,dir);
  //-------------------------------------------------------------
  //check if the first site on the hyperplane is on the node
  //If yes, keep the slice sum on this node
  //-------------------------------------------------------------
  int is_on_node = glb2lcl(lcl,box_b);
  if ( is_on_node ) {
    pv = (Vector *)(d_fermion_p + vec_len * siteOffset(lcl));

    pv->CopyVec(sum,vec_len);
  }
  

  sfree(sum);
#ifdef TIMING
  time = dclock() - time;
  print_flops(d_class_name,fname,0,time);
#endif

}

#define TIMING

void FermionVector::sumOverHyperPlaneZeroMom(int dir, int box_begin[], int box_end[])
{
  char *fname = "sumOverHyperPlaneZeroMom(int,int,int)";
#ifdef TIMING
  Float time = dclock();
#endif
  FermionVector ferm_vec_cpy;
  FermionVector ferm_vec_sum;

  ferm_vec_cpy = *this;

  int glbMin[LORENTZs],glbMax[LORENTZs];
  int box_b[LORENTZs], box_e[LORENTZs],box_size[LORENTZs];
  
  for(int n = 0; n < LORENTZs; n++){
    glbMin[n] = 0;
    glbMax[n] = glb_sites[n] - 1;
    box_size[n] = box_end[n] - box_begin[n] + 1;
  }
  glbMin[dir] = glbMax[dir] = 0;
  
  for(box_b[0]=glbMin[0];box_b[0]<=glbMax[0];box_b[0]++)
    for(box_b[1]=glbMin[1];box_b[1]<=glbMax[1];box_b[1]++)
      for(box_b[2]=glbMin[2];box_b[2]<=glbMax[2];box_b[2]++)
	for(box_b[3]=glbMin[3];box_b[3]<=glbMax[3];box_b[3]++)
	  {
	    for(int n = 0; n < LORENTZs; n++)
	      {
		box_e[n] = (box_b[n] + box_size[n] - 1) % glb_sites[n];
	      }
	    ferm_vec_cpy.sumOverHyperPlaneStride(dir, box_b, box_e); 

	    ferm_vec_sum += ferm_vec_cpy;
	    ferm_vec_cpy = *this;
	  }

  *this = ferm_vec_sum;
#ifdef TIMING
  time = dclock() - time;
  print_flops(d_class_name,fname,0,time);
#endif
}


CPS_END_NAMESPACE
 
