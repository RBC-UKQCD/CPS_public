#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:40 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_w_spect/w_ginfo.C,v 1.7 2004-08-18 11:57:40 zs Exp $
//  $Id: w_ginfo.C,v 1.7 2004-08-18 11:57:40 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: w_ginfo.C,v $
//  $Revision: 1.7 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_w_spect/w_ginfo.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <alg/w_all.h>
#include <util/gjp.h>
#include <util/qcdio.h>
CPS_START_NAMESPACE

//---------------------------------------------------------------------------
// WspectGinfo:: static data members
//---------------------------------------------------------------------------
int WspectGinfo::  initialized = 0; 

int WspectGinfo::  lcl2glb_offset[WspectGinfo::LORENTZs];  
int WspectGinfo::  glb_sites     [WspectGinfo::LORENTZs];
int WspectGinfo::  lcl_sites     [WspectGinfo::LORENTZs]; 
int WspectGinfo::  lcl_node      [WspectGinfo::LORENTZs];  
int WspectGinfo::  bnd_cnd       [WspectGinfo::LORENTZs];  

char *WspectGinfo::ctor_str            = "CTOR";
char *WspectGinfo::dtor_str            = "DTOR";
char *WspectGinfo::empty_str           = "";
char *WspectGinfo::wrong_type_str      = "WrongType";
char *WspectGinfo::out_range_str       = "ArgOutOfRange";
char *WspectGinfo::inconsistent_str    = "Inconsistent";


//---------------------------------------------------------------------------
// WspectGinfo::WspectGinfo() 
//---------------------------------------------------------------------------
WspectGinfo::WspectGinfo() 
{
  if (initialized <= 0) {
    // global sites
    glb_sites[0] = GJP.Xnodes() * GJP.XnodeSites();
    glb_sites[1] = GJP.Ynodes() * GJP.YnodeSites();
    glb_sites[2] = GJP.Znodes() * GJP.ZnodeSites();
    glb_sites[3] = GJP.Tnodes() * GJP.TnodeSites();
  
    // number of sites on the local node
    lcl_sites[0] = GJP.XnodeSites();
    lcl_sites[1] = GJP.YnodeSites();
    lcl_sites[2] = GJP.ZnodeSites();
    lcl_sites[3] = GJP.TnodeSites();

    // the coordinates of the local node among the whole machine
    lcl_node[0] = GJP.XnodeCoor();
    lcl_node[1] = GJP.YnodeCoor();
    lcl_node[2] = GJP.ZnodeCoor();
    lcl_node[3] = GJP.TnodeCoor();
    
    // lcl2glb_offset
    for (int i = 0; i < LORENTZs; ++i) 
      lcl2glb_offset[i] = lcl_sites[i] * lcl_node[i];

    // the boundary conditions of fermions
    bnd_cnd[0] = GJP.Xbc() == BND_CND_PRD ? 1 : -1;
    bnd_cnd[1] = GJP.Ybc() == BND_CND_PRD ? 1 : -1;
    bnd_cnd[2] = GJP.Zbc() == BND_CND_PRD ? 1 : -1;
    bnd_cnd[3] = GJP.Tbc() == BND_CND_PRD ? 1 : -1;

    // set the flag
    initialized = 1;   
  }
}



//---------------------------------------------------------------------------
// WspectGinfo::isOutOfRange(const int point[], const int frame[]) const
//---------------------------------------------------------------------------
int 
WspectGinfo::isOutOfRange(const int point[], const int frame[]) const {
  for (int i = 0; i < LORENTZs; ++i) {
    if (point[i] < 0 || point[i] >= frame[i])
      return 1;
  }
  return 0;
}



//---------------------------------------------------------------------------
// WspectGinfo::siteOffset(const int lcl_site[], int exclude_dir = -1) const
//---------------------------------------------------------------------------
// The storage order [in the C convension of arrays]:
//      something[sites[n-1]] [.]..[sites[0]]
// If the point is specified as [0..3 as x..t], then the storage order is
//      something[t][z][y][x]
// More specifically, 
//    offset =         lcl[0] + lcl_sites[0] * (
//                     lcl[1] + lcl_sites[1] * (
//                     ....
//                  // lcl[e] + lcl_sites[e] * (
//                     ....
//                     lcl[n-2] + lcl_sites[n-2] * lcl[n-1] )../*)*/..)
//---------------------------------------------------------------------------
int 
WspectGinfo::siteOffset(const int lcl[], int e)  const {
  int l = (e == LORENTZs - 1 ? LORENTZs - 2 : LORENTZs - 1);
  int offset = lcl[l];
  while (l-- > 0) {
    if (l != e) {
      offset *= lcl_sites[l];
      offset += lcl[l];
    }
  }
  return offset;
}

//---------------------------------------------------------------------------
// int WspectGinfo::glb2lcl(int lcl[], const int glb[]) const
//---------------------------------------------------------------------------
int 
WspectGinfo::glb2lcl(int lcl[], const int glb[]) const
{
  int is_on_node = 1;  
  for (int l = 0; l < LORENTZs; ++l) {
    lcl[l] = glb[l] - lcl2glb_offset[l];
    if (lcl[l] < 0 || lcl[l] >= lcl_sites[l]) 
      is_on_node = 0;    
  }
  return is_on_node;
}

//---------------------------------------------------------------------------
// int WspectGinfo::lcl2glb(const int lcl[], int glb[]) const
//---------------------------------------------------------------------------
void WspectGinfo::lcl2glb(const int lcl[], int glb[]) const
{
  for (int l = 0; l < LORENTZs; ++l) {
    glb[l] = lcl[l] + lcl2glb_offset[l];
  }
}


//---------------------------------------------------------------------------
// void WspectGinfo::printSite(FILE *fp, const int site[]) const
//---------------------------------------------------------------------------
    
void 
WspectGinfo::printSite(FILE *fp, const int x[]) const
{  
  if (fp) {
    Fprintf(fp, "Site [x,y,z,t] = [%d, %d, %d, %d]\n",
            x[0],x[1],x[2],x[3]);
  } 
}


//---------------------------------------------------------------------------
// void WspectGinfo::printSpinor(FILE *fp, const IFloat *flt_p) const
//---------------------------------------------------------------------------
void 
WspectGinfo::printSpinor(FILE *fp, const IFloat *flt_p) const
{  
  if (fp) {
    for (int d = 0; d < DIRACs; ++d) {
      Fprintf(fp, " [[%g, %g], [%g, %g], [%g, %g]]\n",
	      flt_p[0], flt_p[1], flt_p[2], flt_p[3], flt_p[4], flt_p[5]);
      flt_p += COLORs * COMPLEXs;
    }
  }
}


CPS_END_NAMESPACE
