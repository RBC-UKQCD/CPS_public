#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_w_spect/w_axialcurr.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: w_axialcurr.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.7  2002/03/11 22:25:46  anj
//  This should now be the correct, fully merged code from our two versions. Anj
//
//  Revision 1.4.2.1  2002/03/08 16:35:05  anj
//  Checking in the Columbia code branch on tag Columbia4_1_1_test-branch, to be
//  merged with the UKQCD head branch shortly.  Anj
//
//  Revision 1.4  2001/08/16 10:49:42  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:11:34  anj
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
//  Revision 1.2  2001/05/25 06:16:00  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: w_axialcurr.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_w_spect/w_axialcurr.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include<alg/w_all.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <stdlib.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include<util/gjp.h>
#include<util/error.h>
#include<util/verbose.h>
#include<util/sproj_tr.h>
#include<util/lattice.h>
#include<util/data_types.h>
CPS_START_NAMESPACE


CPS_END_NAMESPACE
#include<comms/glb.h>
#include<comms/scu.h>
#include<alg/alg_w_spect.h>
CPS_START_NAMESPACE


//---------------------------------------------------------------------------
// static data members
//---------------------------------------------------------------------------
char *      WspectAxialCurrent::d_class_name = "WspectAxialCurrent";

//---------------------------------------------------------------------------
// WspectAxialCurrent::WspectAxialCurrent(...)
//--------------------------------------------------------------------------- 
WspectAxialCurrent::WspectAxialCurrent(Lattice &     lat, 
				       const WspectHyperRectangle & whr,
				       char * ap_corr_outfile )
  :  d_lat(lat), d_whr(whr)
{

  VRB.Func(d_class_name, ctor_str);

  if (lat.Fclass() == F_CLASS_DWF && ap_corr_outfile != 0) {
    
    ap_filename = ap_corr_outfile;

    // Initiate total Ls, propagation direction and its length, ...
    //-----------------------------------------------------------------------
    ls_glb = GJP.SnodeSites() * GJP.Snodes();
    prop_dir  = d_whr.dir();  
    glb_walls = glb_sites[prop_dir];
    {
      const int *low  = d_whr.lclMin();
      const int *high = d_whr.lclMax();
      for (int i = 0; i < LORENTZs; ++i) {
	lclMin[i] = low[i];
	lclMax[i] = high[i];
      }
    }
    
    // allocate and clear the space for the <A_0 P> results
    //-----------------------------------------------------------------------
    
    // since only real part of the trace is taken, save result as Float
    {
      int fsize = glb_walls * sizeof(Float);
      
      // For local axial current
      d_local_p = (Float *) smalloc(fsize);
      
      if (!d_local_p)
	ERR.Pointer(d_class_name, ctor_str, empty_str);
      
      VRB.Smalloc(d_class_name, ctor_str, empty_str, d_local_p, fsize);
      
      Float *flt_p = (Float *)d_local_p;
      {
	for ( int i = 0; i < glb_walls; i++)
	  *flt_p++ = 0.0;
      }
      
      // For conserved axial current
      d_conserved_p = (Float *) smalloc(fsize);
      
      if (!d_conserved_p)
	ERR.Pointer(d_class_name, ctor_str, empty_str);
      
      VRB.Smalloc(d_class_name, ctor_str, empty_str, d_conserved_p, fsize);
      
      flt_p = (Float *)d_conserved_p;
      {
	for ( int i = 0; i < glb_walls; i++)
	  *flt_p++ = 0.0;
      }
    }    
  
    // allocate space for two 4d field 
    //-----------------------------------------------------------------------
    
    {
      int d_size_4d = GJP.VolNodeSites() * SPINORs ;
      
      d_data_p1 = (IFloat *) smalloc(d_size_4d * sizeof(IFloat));
      if ( d_data_p1 == 0)
	ERR.Pointer(d_class_name,ctor_str, "d_data_p1");
      VRB.Smalloc(d_class_name,ctor_str, "d_data_p1", d_data_p1,
		  d_size_4d * sizeof(IFloat));
      
      d_data_p2 = (IFloat *) smalloc(d_size_4d * sizeof(IFloat));
      if ( d_data_p2 == 0)
	ERR.Pointer(d_class_name,ctor_str, "d_data_p2");
      VRB.Smalloc(d_class_name,ctor_str, "d_data_p2", d_data_p2,
		  d_size_4d * sizeof(IFloat));
    }
    
    //-----------------------------------------------------------------------
    // allocate space for two spinor fields for SCU transfer
    //-----------------------------------------------------------------------
    
    {
      tmp_p1 = (Float *)smalloc(SPINORs * sizeof(Float));
      if (tmp_p1 == 0)
	ERR.Pointer(d_class_name, ctor_str, "tmp_p1") ;
      VRB.Smalloc(d_class_name, ctor_str, "tmp_p1", tmp_p1,
		  SPINORs * sizeof(Float)) ;
      
      tmp_p2 = (Float *)smalloc(SPINORs * sizeof(Float));
      if (tmp_p2 == 0)
	ERR.Pointer(d_class_name, ctor_str, "tmp_p2") ;
      VRB.Smalloc(d_class_name, ctor_str, "tmp_p2", tmp_p2,
		  SPINORs * sizeof(Float)) ;
      
    }
    
    //-----------------------------------------------------------------------
    // allocate space for two spinor fields for gamma_5 multiplication 
    //-----------------------------------------------------------------------
    
    {
      v1_g5 = (Float *)smalloc(SPINORs * sizeof(Float));
      if (v1_g5 == 0)
	ERR.Pointer(d_class_name, ctor_str, "v1_g5") ;
      VRB.Smalloc(d_class_name, ctor_str, "v1_g5", v1_g5,
		  SPINORs * sizeof(Float)) ;
      
      v1_next_g5 = (Float *)smalloc(SPINORs * sizeof(Float));
      if (v1_next_g5 == 0)
	ERR.Pointer(d_class_name, ctor_str, "v1_next_g5") ;
      VRB.Smalloc(d_class_name, ctor_str, "v1_next_g5", v1_next_g5,
		  SPINORs * sizeof(Float)) ;
      
    }

    //-----------------------------------------------------------------------
    // Fill in the array of sproj_tr functions
    //-----------------------------------------------------------------------
    sproj_tr[SPROJ_XM] = sprojTrXm;
    sproj_tr[SPROJ_YM] = sprojTrYm;
    sproj_tr[SPROJ_ZM] = sprojTrZm;
    sproj_tr[SPROJ_TM] = sprojTrTm;
    sproj_tr[SPROJ_XP] = sprojTrXp;
    sproj_tr[SPROJ_YP] = sprojTrYp;
    sproj_tr[SPROJ_ZP] = sprojTrZp;
    sproj_tr[SPROJ_TP] = sprojTrTp;

  }    

}

//---------------------------------------------------------------------------
// WspectAxialCurrent::~WspectAxialCurrent()
//--------------------------------------------------------------------------- 
// Purpose:
//    Free all the memory on heap.
//--------------------------------------------------------------------------- 
WspectAxialCurrent::~WspectAxialCurrent()
{

  VRB.Func(d_class_name, dtor_str);

  if (d_lat.Fclass() == F_CLASS_DWF && ap_filename != 0) {

    VRB.Sfree(d_class_name, dtor_str, empty_str, v1_next_g5);
    sfree(v1_next_g5);
    
    VRB.Func(d_class_name, dtor_str);
    VRB.Sfree(d_class_name, dtor_str, empty_str, v1_g5);
    sfree(v1_g5);
    
    VRB.Func(d_class_name, dtor_str);
    VRB.Sfree(d_class_name, dtor_str, empty_str, tmp_p2);
    sfree(tmp_p2);
    
    VRB.Func(d_class_name, dtor_str);
    VRB.Sfree(d_class_name, dtor_str, empty_str, tmp_p1);
    sfree(tmp_p1);
    
    VRB.Func(d_class_name, dtor_str);
    VRB.Sfree(d_class_name, dtor_str, empty_str, d_data_p2);
    sfree(d_data_p2);
    
    VRB.Func(d_class_name, dtor_str);
    VRB.Sfree(d_class_name, dtor_str, empty_str, d_data_p1);
    sfree(d_data_p1);
    
    VRB.Func(d_class_name, dtor_str);
    VRB.Sfree(d_class_name, dtor_str, empty_str, d_conserved_p);
    sfree(d_conserved_p);
    
    VRB.Func(d_class_name, dtor_str);
    VRB.Sfree(d_class_name, dtor_str, empty_str, d_local_p);
    sfree(d_local_p);
  }
}

//---------------------------------------------------------------------------
// void WspectAxialCurrent::measureAll(...)
//--------------------------------------------------------------------------- 
// Purpose:
//    Does the <A_0 P> calculation for DWF lattice. 
//    Calls function measureConserved for conserved axial current A0,
//    and  function measureLocal for local axial current A0.
//---------------------------------------------------------------------------
void WspectAxialCurrent::measureAll(Vector * data_5d_p) {

    measureConserved(data_5d_p);
  
    measureLocal(data_5d_p);
}


//---------------------------------------------------------------------------
// void WspectAxialCurrent::measureConserved(...)
//--------------------------------------------------------------------------- 
// Purpose:
//    Does the <A_0 P> calculation for DWF lattice
//    where A_0 is the conserved axial current 
//---------------------------------------------------------------------------
void WspectAxialCurrent::measureConserved(Vector * data_5d_p) {

  Matrix *gauge_field = d_lat.GaugeField();  

  //  To avoid duplicate work, stop at the middle of the s direction
  for (int s = 0; s < ls_glb/2; s++) {
    d_lat.Ffive2four((Vector *) d_data_p1, data_5d_p, s, s);
    d_lat.Ffive2four((Vector *) d_data_p2, data_5d_p, ls_glb-1-s, ls_glb-1-s);

    {
      int lcl_walls = lcl_sites[prop_dir];
      
      for( int lclW = 0; lclW < lcl_walls; ++lclW) {
	int lcl[LORENTZs];
	int lcl_next[LORENTZs]; // Next site along propagation direction

	// Define hyperplane    
	lclMin[prop_dir]=lclMax[prop_dir] = lclW;

	for(lcl[0]=lclMin[0]; lcl[0]<=lclMax[0]; lcl[0]++) 
	for(lcl[1]=lclMin[1]; lcl[1]<=lclMax[1]; lcl[1]++) 
	for(lcl[2]=lclMin[2]; lcl[2]<=lclMax[2]; lcl[2]++) 
        for(lcl[3]=lclMin[3]; lcl[3]<=lclMax[3]; lcl[3]++) {

	  int lcl_offset = siteOffset(lcl) * SPINORs ;

	  // coordinates and offset for lcl_next
	  for (int i  = 0; i < LORENTZs; i++ ) 
	    lcl_next[i] = ( (i == prop_dir) ? (lcl[i]+1)%lcl_sites[i]
			    : lcl[i] );

	  int lcl_next_offset = siteOffset(lcl_next) * SPINORs;

	  // U_mu(x) where mu = prop_dir
	  Matrix * link = gauge_field + siteOffset(lcl) * 4 + prop_dir ;

	  { 
	    Float * v1, * v2;
	    Float * v1_next, * v2_next;
	    Float coeff = 1.0;
	    
	    // S_F(x, s)
	    v1 = (Float *)d_data_p1+lcl_offset ;
	    // S_F(x, ls_glb-1-s)
	    v2 = (Float *)d_data_p2+lcl_offset ;

	    // v1_next = S_F(x+prop_dir, s)
	    // v2_next = S_F(x+prop_dir, ls_glb-1-s)
	    if ((lcl[prop_dir]+1) == lcl_sites[prop_dir]) {
	      getPlusData( (IFloat *)tmp_p1,
			   (IFloat *)d_data_p1+lcl_next_offset, SPINORs, 
			   prop_dir) ;
	      getPlusData( (IFloat *)tmp_p2,
			   (IFloat *)d_data_p2+lcl_next_offset, SPINORs, 
			   prop_dir) ;
	      v1_next = tmp_p1;
	      v2_next = tmp_p2;

	      // fix boundary condition
	      switch( prop_dir ) {
	      case 0:
		if (GJP.XnodeBc()==BND_CND_APRD) coeff = -coeff ;
		break;
	      case 1:
		if (GJP.YnodeBc()==BND_CND_APRD) coeff = -coeff ;
		break;
	      case 2:
		if (GJP.ZnodeBc()==BND_CND_APRD) coeff = -coeff ;
		break;
	      case 3:
		if (GJP.TnodeBc()==BND_CND_APRD) coeff = -coeff ;
		break;
	      } // end switch

	    } else {
	      v1_next = (Float *)d_data_p1+lcl_next_offset ;
	      v2_next = (Float *)d_data_p2+lcl_next_offset ;
	    }
	    
	    {
	      // Gamma^5 S_F(x, s)
	      d_lat.Gamma5((Vector *) v1_g5, (Vector *) v1, 1 );

	      // Gamma^5 S_F(x+prop_dir, s)
	      d_lat.Gamma5((Vector *) v1_next_g5, (Vector *) v1_next, 1 );
	    
	      Matrix tmp1, tmp2, f;
	      Float result = 0.0 ;

	      // tmp1 = Tr_spin ( (1 + gamma_{prop_dir} ) 
	      //       gamma_5 S_F(x+prop_dir, s) S_F^dagger(x, ls_glb-1-s))
	      sproj_tr[prop_dir+4]((IFloat *)&tmp1,
				   (IFloat *)v1_next_g5,
				   (IFloat *)v2, 1, 0, 0);
	      
	      f.DotMEqual(*link, tmp1);

	      result = f.ReTr();
	      
	      // tmp2 = Tr_spin ( (1 - gamma_{prop_dir} ) 
	      //       gamma_5 S_F(x, s) S_F^dagger(x+prop_dir, ls_glb-1-s))
	      sproj_tr[prop_dir]( (IFloat *)&tmp2,
					(IFloat *)v1_g5,
					(IFloat *)v2_next, 1, 0, 0);
					
	      tmp1.Dagger(*link);
	      f.DotMEqual(tmp1, tmp2);
	      result -= f.ReTr();

	      result *= coeff;

	      *( d_conserved_p + lclW + lcl2glb_offset[prop_dir] ) += result;
	    }

	  }
	} // for(lcl[_]..)
      } // for (int lclW...)
    }
  } // for (int s ... )

}


//---------------------------------------------------------------------------
// void WspectAxialCurrent::measureLocal(...)
//--------------------------------------------------------------------------- 
// Purpose:
//    Does the <A_0 P> calculation for DWF lattice
//    where A_0 is the local axial current 
//---------------------------------------------------------------------------
void WspectAxialCurrent::measureLocal(Vector * data_5d_p) { 

  d_lat.Ffive2four((Vector *) d_data_p1, data_5d_p, ls_glb-1, 0);
  
  int lcl_walls = lcl_sites[prop_dir];
  
  for( int lclW = 0; lclW < lcl_walls; ++lclW) {
    int lcl[LORENTZs];
      
    // Define hyperplane    
    lclMin[prop_dir]=lclMax[prop_dir] = lclW;
      
    for(lcl[0]=lclMin[0]; lcl[0]<=lclMax[0]; lcl[0]++) 
    for(lcl[1]=lclMin[1]; lcl[1]<=lclMax[1]; lcl[1]++) 
    for(lcl[2]=lclMin[2]; lcl[2]<=lclMax[2]; lcl[2]++) 
    for(lcl[3]=lclMin[3]; lcl[3]<=lclMax[3]; lcl[3]++) {

      int lcl_offset = siteOffset(lcl) * SPINORs ;
      
      Complex * v1 = (Complex *) (d_data_p1 + lcl_offset) ;

      Complex d_proj[DIRACs][DIRACs];
      
      {
	IFloat *p = (IFloat *) d_proj;
	for(int i = 0; i < COMPLEXs*DIRACs*DIRACs; i++)
	  *p++ = 0;
      }

      // Color Algebra -- only necessary Dirac indexes are calculated

      {
	for(int D1x = 0; D1x < DIRACs/2; D1x ++ )
	  for(int D2x = DIRACs/2; D2x< DIRACs; D2x ++) {
	    
	    const Complex * v1_p = v1 + COLORs * D1x;
	    const Complex * v2_p = v1 + COLORs * D2x;
	    
	    d_proj[D1x][D2x] += v2_p[0]*conj(v1_p[0]);
	    d_proj[D1x][D2x] += v2_p[1]*conj(v1_p[1]);
	    d_proj[D1x][D2x] += v2_p[2]*conj(v1_p[2]);
	    
	    d_proj[D2x][D1x] += v1_p[0]*conj(v2_p[0]);
	    d_proj[D2x][D1x] += v1_p[1]*conj(v2_p[1]);
	    d_proj[D2x][D1x] += v1_p[2]*conj(v2_p[2]);
	  }

      }

      Complex result(0.0, 0.0), I(0, 1);

      // calculate {- Tr_{spin, color} (S_F^dagger gamma_{prop_dir} S_F)}

      switch(prop_dir) {
      case 0 :
	result -= d_proj[0][3];
	result -= d_proj[1][2];
	result += d_proj[2][1];
	result += d_proj[3][0];
	result *= I;
	break;
      case 1 :
	result += d_proj[0][3];
	result -= d_proj[1][2];
	result -= d_proj[2][1];
	result += d_proj[3][0];
	break;
      case 2 :
	result -= d_proj[0][2];
	result += d_proj[1][3];
	result += d_proj[2][0];
	result -= d_proj[3][1];
	result *= I;
	break;
      case 3 :
	result -= d_proj[0][2];
	result -= d_proj[1][3];
	result -= d_proj[2][0];
	result -= d_proj[3][1];
	break;
      }

      *( d_local_p + lclW + lcl2glb_offset[prop_dir] ) += result.real(); 

    } // for(lcl[_]..)
  } // for (int lclW...) 
  
}

//---------------------------------------------------------------------------
// void WspectAxialCurrent::doSum()
//--------------------------------------------------------------------------- 
void WspectAxialCurrent::doSum() 
{
  char *fname = "doSum";
  VRB.Func(d_class_name, fname);

  // Global sum over all data
  
  {
    Float *flt_p = d_conserved_p;
    for (int i=0; i < glb_walls; i++) 
      glb_sum(flt_p++);
  }
  
  // Global sum over all data
  
  {
    Float *flt_p = d_local_p;
    for (int i=0; i < glb_walls; i++) 
      glb_sum(flt_p++);
  }
}

//---------------------------------------------------------------------------
// void WspectAxialCurrent::print(char filename) const
//--------------------------------------------------------------------------- 
void WspectAxialCurrent::print() const
{
  char *fname = "print";
  VRB.Func(d_class_name, fname);

  // Print out correlator data 
  //------------------------------------------------------------------

  FILE *fp;

  if ( !ap_filename || !(fp = fopen(ap_filename, "a")) )
    ERR.FileA(d_class_name,fname, ap_filename);

  for (int wall = 0; wall < glb_walls; ++wall) {

    IFloat conserved_result, local_result;
    conserved_result = 
      d_conserved_p[(d_whr.glbCoord() + wall)%glb_walls];

    local_result = d_local_p[(d_whr.glbCoord() + wall)%glb_walls];
     
    fprintf(fp, "%d %d %e %e\n", 
	    AlgWspect::GetCounter(), wall, 
	    local_result,
	    conserved_result);
  }
  fclose(fp);
  
}

//---------------------------------------------------------------------------
// void WspectAxialCurrent::dumpData() const 
//--------------------------------------------------------------------------- 
void WspectAxialCurrent::dumpData(char *filename) const {
  FILE *fp;

  if (filename && (fp = fopen(filename, "a"))) {
    for (int i = 0; i < glb_walls; ++i) {
      fprintf(fp, "%e %e\n", d_local_p[i], d_conserved_p[i]);
    }
    fclose(fp);    
  } else {
    ERR.FileA(d_class_name, "dumpData", filename);
  }

}


CPS_END_NAMESPACE
