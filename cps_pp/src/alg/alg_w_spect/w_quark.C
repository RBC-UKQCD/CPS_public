/*! \file

  $Id: w_quark.C,v 1.16 2007/06/06 16:06:22 chulwoo Exp $
*/
#include<config.h>
#include <util/gjp.h>              // GJP
#include <util/error.h>            // ERR
#include <util/verbose.h>          // VRB
#include <util/lattice.h>          // Lattice::FixGaugePtr()
#include <util/qcdio.h> 
#include <comms/glb.h>               // glb_sum(...)
#include <comms/scu.h>               // getMinusData, getPlusData
#include <alg/w_all.h>
#include <alg/w_ginfo.h>

CPS_START_NAMESPACE

 
#undef DEBUG_W_QUARK
//#define DEBUG_CHECKSU3
//#define DEBUG_W_WAVEFUNC
#ifdef DEBUG_W_WAVEFUNC
static char *wavefilename="wavefunc.dat.*************"; //initial value
#endif

//---------------------------------------------------------------------------
// static data members
//---------------------------------------------------------------------------
char *   WspectQuark::d_class_name = "WspectQuark";
int      WspectQuark::d_weight_Dy;
int      WspectQuark::d_weight_Cy;
Matrix WspectQuark::m_tmp1; //buffer for GetSource
Vector WspectQuark::v_tmp1; //buffer for GetPropData
Vector WspectQuark::v_tmp2[WspectGinfo::DIRACs];


//---------------------------------------------------------------------------
// WspectQuark::WspectQuark(...)
//------------------------------------------------------
// src_op_kind has default value UNIT in function delcaration
// added arguments pbp_outfile, mid_point_outfile and a0_p_outfile as in phys_v4.0.0
//--------------------- ------------------------------------------------
WspectQuark::WspectQuark(Lattice &lat,
			 char *outfile,
                         char *        pbp_outfile,
                         char *        mid_point_outfile,
                         char *        a0_p_outfile,
			 WspectArg &warg,
			 CgArg &cg,
			 const WspectHyperRectangle & whr,
			 DEVOperatorKind source_operator_kind,
			 WspectFuzzing *src_fuzz_ptr,
			 WspectField *fld_ptr)
  : d_size(GJP.VolNodeSites() * (COLORs*COLORs*DIRACs*DIRACs*COMPLEXs)), prop_direction(whr.dir()), src_plane_position(whr.lclCoord()),d_lat(lat), src_op_kind(source_operator_kind), src_fuzz_p(src_fuzz_ptr), fld_p(fld_ptr), midplane(warg.midplane)
{
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif


  VRB.Func(d_class_name, ctor_str);
  
  // set d_weight_* to speed up the accesses to quark propagators
  //-------------------------------------------------------------------------
  {    
    d_weight_Cy = lcl_sites[0]*lcl_sites[1]*lcl_sites[2]*lcl_sites[3]
      *SPINORs;  
    d_weight_Dy = d_weight_Cy * COLORs;
  }

#ifdef DEBUG_W_QUARK
  //printf("constructing WspectQuark  \n");
#endif 
  
  // Open the file where the cg monitor parameters are written
  //-------------------------------------------------------------------------
  FILE *fp;
  if (outfile) {
    if (!(fp = Fopen(outfile, "a"))) {
      ERR.FileA(d_class_name, ctor_str, outfile);
    }
  }
  
  // Allocate space for the quark propagator
  //-------------------------------------------------------------------------
  d_data_p = (IFloat *) smalloc(d_size * sizeof(IFloat));
#ifdef DEBUG_W_QUARK
  printf("allocated quark propagator - %6i Floats at %x \n",d_size,d_data_p);
#endif 
  if (!d_data_p) 
    ERR.Pointer(d_class_name, ctor_str, "d_data_p");
  VRB.Smalloc(d_class_name, ctor_str, 
	      "d_data_p", d_data_p, d_size * sizeof(IFloat));
  
  // ======== added from phys_v4.0.0
  // In case of DWF, allocate space need for propagator for
  // calculating <\Delta J^5 \bar q \gamma^5 q> using mid_point sink
  if (lat.Fclass() == F_CLASS_DWF && mid_point_outfile !=0 ) {
     d_data_mid_point_1 = (IFloat *) smalloc(d_size * sizeof(IFloat));
    if (!d_data_mid_point_1)
      ERR.Pointer(d_class_name, ctor_str, "d_data_mid_point_1");
    VRB.Smalloc(d_class_name, ctor_str,
		"d_data_mid_point_1", d_data_mid_point_1,
		d_size * sizeof(Float));
    
    
    d_data_mid_point_2 = (IFloat *) smalloc(d_size * sizeof(IFloat));
    if (!d_data_mid_point_2)
      ERR.Pointer(d_class_name, ctor_str, "d_data_mid_point_2");
    VRB.Smalloc(d_class_name, ctor_str,
		"d_data_mid_point_2", d_data_mid_point_2,
		d_size * sizeof(Float));
    
  } else {
       d_data_mid_point_1 = 0;
       d_data_mid_point_2 = 0;  
  }
  // ======== end added from phys_v4.0.0 =========================

  //===============================================================
  // Allocate space for the source matrix S^{C1C2}(x)
  // added by T&X, March 30 2000
  //----------------------------------------------------------------
  source_matrix_size = 1;
  for (int l=0; l < LORENTZs; l++){
    if (l != prop_direction) source_matrix_size *=lcl_sites[l];
  }
  source_matrix_size *= COLORs*COLORs*COMPLEXs;
  
    
  source_matrix_p = (IFloat *) smalloc(source_matrix_size* sizeof(IFloat));
#ifdef DEBUG_W_QUARK
  printf("allocated source matrix - %6i Floats at %x\n",source_matrix_size, source_matrix_p);
#endif 
  if (!source_matrix_p) 
    ERR.Pointer(d_class_name, ctor_str, "source_matrix");
  VRB.Smalloc(d_class_name, ctor_str, 
	      "source_matrix", source_matrix_p, source_matrix_size * sizeof(Float));
  

  // initialisation of source_matrix
  {
    for (int i = 0; i < source_matrix_size; i++) source_matrix_p[i] = 0.0;
  }
  
  //===============================================================

  // added new to phys_v4.1.0
  // rescale factor for sources -- default=1.0
  Float rs_fac = 1.0;
 
  // now let me hold the throat of somebody who wrote line below!!!
  // if (!rs_fac) rs_fac = 1.0;
  //Float rs_fac = 1.0E+5;
  printf("rescaled source with rs_fac = %e \n",rs_fac);
  // end added
  
  // For POINT_W only, This is also the initial source for JACOBI source!
  //----------------------------------------------------------------
  const int *point = 0;  //0 means source point is off-node
  int src_point[LORENTZs];  //global coordinate of source point
  if (whr.onNode()) {
    src_point[0] = 0;
    src_point[1] = 0;
    src_point[2] = 0;
    src_point[3] = 0;
    src_point[prop_direction] = whr.glbCoord();

    if (glb2lcl(src_point, src_point)) {    
      point = src_point;   // if src_point on-node ->= local corrd
      {
	int i = siteOffset(point,prop_direction)*COMPLEXs*COLORs*COLORs;
       	source_matrix_p[i] = 1.0;
	source_matrix_p[i+COMPLEXs*(COLORs+1)] = 1.0;
	source_matrix_p[i+2*COMPLEXs*(COLORs+1)] = 1.0;
      }  
    }
  }


  //------------------------------------------------------------------------
  // Kuramashi source : assign unit source for every lattice site and every
  // spin and color without gauge fixing. 
  // added by Meifeng Lin, 05/2006
  //------------------------------------------------------------------------
  
  if ( warg.source_kind == KURAMASHI ) {
    int per_site = COLORs * COLORs * COMPLEXs;  
    for (int i = 0; i< source_matrix_size; i += per_site) {
      source_matrix_p[i] =  1.0;
      source_matrix_p[i+COMPLEXs*(COLORs+1)] = 1.0;
      source_matrix_p[i+2*COMPLEXs*(COLORs+1)] = 1.0;
    }
  }
  
  //------------------------------------------------------------------------
  // Z2 noise source. Every spin-color point has random noise of 1 or -1
  // added by Meifeng Lin, 05/2006
  //------------------------------------------------------------------------

  if ( warg.source_kind == Z2 || warg.source_kind == COMPLEX_Z2 ) {
//    LRGState rng_state;
    //store random number seeds
    
//    rng_state.GetStates();
    LRG.SetInterval(0.0,1.0);
    int per_site = COLORs * COLORs * COMPLEXs;  
    for (int i = 0; i< source_matrix_size; i += per_site) {
      
      source_matrix_p[i] =  ( LRG.Urand() > 0.5 )? 1.0 : -1.0;
      source_matrix_p[i+COMPLEXs*(COLORs+1)] = ( LRG.Urand() > 0.5 ) ? 1.0 : -1.0;
      source_matrix_p[i+2*COMPLEXs*(COLORs+1)] = ( LRG.Urand() > 0.5 ) ? 1.0 : -1.0;

      if ( warg.source_kind == COMPLEX_Z2 ) {
	source_matrix_p[i + 1] =  ( LRG.Urand() > 0.5 )? 1.0 : -1.0;
	source_matrix_p[i+COMPLEXs*(COLORs+1) + 1] = ( LRG.Urand() > 0.5 ) ? 1.0 : -1.0;
	source_matrix_p[i+2*COMPLEXs*(COLORs+1) + 1] = ( LRG.Urand() > 0.5 ) ? 1.0 : -1.0;
      }
    }
    
    //restore random number seeds
//    rng_state.SetStates();

  }

  //-----------------------------------------------------------------------
  //

  // For WALL_W and BOX_W only:  
  // Thomas & Xiaodong:
  // after change from setPointSrc and setSmearedSrc to --> setSource 
  // it is now important that src_lower_bound and source_upper_bound 
  // is defined for all sources
  // the projection is now unified
  //-------------------------------------------------------------------------
  const Float *src_matrix      = 0; //temporary src_matrix
  const int   *src_lower_bound = whr.lclMin(); 
  const int   *src_upper_bound = whr.lclMax(); 
  
  int          prop_dir        = prop_direction;  

  

  
  // For BOX_W only:             lower_bound, upper_bound
  //     shift origin and end of box source to aots time slice
  //-------------------------------------------------------------------------
  int src_low[LORENTZs]; 
  int src_upp[LORENTZs];

  if (whr.onNode() && warg.source_kind == BOX_W) {

    glb2lcl(src_low, warg.src_box_b);
    glb2lcl(src_upp, warg.src_box_e);

    src_low[prop_dir] = src_upp[prop_dir] = src_plane_position;

    int l;
    
    // check whether the source is valid [a point, at least]
    int src_dimension = 0;      
    for (l = 0; l < LORENTZs; ++l) {
      if (src_low[l] > src_upp[l])  
	ERR.General(d_class_name, ctor_str, "%s", inconsistent_str);
      if (src_low[l] < src_upp[l])  
	++src_dimension;
    }
    if (src_dimension < 1 || prop_dir != warg.prop_dir) {
      ERR.General(d_class_name, ctor_str, "%s", inconsistent_str);
    }

    // find out the on-node part of the source
    for (l = 0; l < LORENTZs; ++l) {
      // src_matrix = 0 if off-node in any one direction
      if (src_low[l] > src_upper_bound[l] || src_upp[l] < src_lower_bound[l])
	src_matrix = 0;
      
      if (src_low[l] < src_lower_bound[l])  src_low[l] = src_lower_bound[l];
      if (src_upp[l] > src_upper_bound[l])  src_upp[l] = src_upper_bound[l];
    }

    // redirect the pointers 
    src_lower_bound = src_low;
    src_upper_bound = src_upp;

      
     
  }
#ifdef DEBUG_W_QUARK
  //printf("src_lower_bound = %4i,%4i%4i,%4i\n",src_lower_bound[0],src_lower_bound[1],src_lower_bound[2],src_lower_bound[3]);
  //printf("src_upper_bound = %4i,%4i%4i,%4i\n",src_upper_bound[0],src_upper_bound[1],src_upper_bound[2],src_upper_bound[3]);
  //printf("source_kind %i  \n",warg.source_kind);
#endif 
  

  // d_source_center2[LORENTZs] for WspectMomenta
  //-------------------------------------------------------------------------
  //what about JACOBI_W ??
  switch (warg.source_kind) {
  case POINT_W:
    { 
      for (int l = 0; l < LORENTZs; ++l) 
	d_source_center2[l] = 0; 
      break; 
    }
  case BOX_W:
    { 
      for (int l = 0; l < LORENTZs; ++l) 
	d_source_center2[l] = warg.src_box_b[l] + warg.src_box_e[l]; 
      break; 
    }
  case WALL_W:
    { 
      for (int l = 0; l < LORENTZs; ++l) 
	d_source_center2[l] = glb_sites[l] - 1; 
      break;
    }
  case JACOBI_W:
    { 
      for (int l = 0; l < LORENTZs; ++l) 
	d_source_center2[l] = 0;
      break;
    }
  case Z2:
    { 
      for (int l = 0; l < LORENTZs; ++l) 
	d_source_center2[l] = glb_sites[l] - 1; 
      break; 
    }
  case COMPLEX_Z2:
    { 
      for (int l = 0; l < LORENTZs; ++l) 
	d_source_center2[l] = glb_sites[l] - 1; 
      break; 
    }
    
  case KURAMASHI:
    { 
      for (int l = 0; l < LORENTZs; ++l) 
	d_source_center2[l] = glb_sites[l] - 1; 
      break; 
    }
  }
  // ========= added from phys_v4.0.0
  Float pbp_norm;
  Float pbp = 0. ;
  Float pbg5p = 0. ;

  // Allocate memory for gamma5 * solution
  //-------------------------------------------------------------------------

  int f_size = GJP.VolNodeSites() * SPINORs;

  Vector *sol_g5 = 0;
  if ( (warg.source_kind==POINT_W) && pbp_outfile ) {
    sol_g5 = (Vector *) smalloc(f_size * sizeof(Float));

    if(sol_g5 == 0)
      ERR.Pointer(d_class_name, ctor_str, "sol_g5");
    VRB.Smalloc(d_class_name, ctor_str, "sol_g5", sol_g5,
                f_size * sizeof(Float));
  }

  // Calulate pbp normalization factor
  //-------------------------------------------------------------------------

  if ((lat.Fclass() == F_CLASS_WILSON) ||
      (lat.Fclass() == F_CLASS_CLOVER)) {
    pbp_norm = (4.0 + cg.mass) * COLORs * DIRACs;
  }
  else if (lat.Fclass() == F_CLASS_DWF) {
    pbp_norm = (4.0 + GJP.DwfA5Inv() - GJP.DwfHeight()) * COLORs * DIRACs;
  }
  else {
    ERR.General(d_class_name, ctor_str, wrong_type_str);
  }

  // ======= end added from phys_v4.0.0 ==================



  // Set the source and call CG  [COLORs x DIRACs] times
  // added from phys_v4.0.0: Caculate pbp and pbg5p at the same time
  //-------------------------------------------------------------------------
  FermionVector source;  

  // ======= added from phys_v4.0.0 ==================
  // For <A_0 P> calculation -- only used for DWF
  WspectAxialCurrent axialCurrent(lat, warg, whr, a0_p_outfile);
  // ======= end added from phys_v4.0.0 ==================
  

  // if Box source or Wall Source is set
  // then the source_slice should be the gaugefixed equivalent of
  // unit matrices everywhere
  // the upper and lower restricitions for boxes are implemented later
  // through setSource
  if(whr.onNode()){
    if (warg.source_kind == BOX_W || warg.source_kind == WALL_W){
      if ( !(lat.FixGaugePtr())[src_plane_position] ) 
	ERR.Pointer(d_class_name,ctor_str,"Gauge Fixing Pointer\n");
      
      src_matrix = (const Float *)((lat.FixGaugePtr())[src_plane_position]); // G(x)
      copyGaugefixToSource(src_matrix, whr);   // copy G^{\dag}(x) --> S(x) 
      //copied from src_matrix to source_matrix_p
#ifdef DEBUG_W_QUARK
      printf("copy gauge fixing matrix onto source done \n");
#endif
    }
  }

  // JACOBI SMEARING
  // modifies source_matrix_p = S(x) = complex 3x3 matrix at every lattice site
  // modifies source S(x) --> S(x) = (1-epsi/N*Delta^2)^N S(x)
  Float epsi=warg.g_epsi;
  // Jacobi-iterations from  w_spect_arg_Gaussian[epsi_count].g_n = GAUSSIAN_PARAM_N;
  int n_iter = warg.g_n;  
  if (warg.source_kind == JACOBI_W && n_iter>=1) {
#ifdef DEBUG_W_WAVEFUNC
    sprintf(wavefilename,"%sE%.2fN%d","wavefunc.dat",epsi,n_iter);
#endif
    Float r = epsi/n_iter;
#ifdef DEBUG_W_QUARK
    printf("JacobiSource, n_iter = %4i, epsilon = %18.8e, r = %18.8e  \n",n_iter,epsi,r);
#endif
    // loop over number of iterations
    for (int i_iter=0; i_iter < n_iter; i_iter++){
      JacobiSource(lat, whr, r);
    }
    
  }
  
  // modifies source S(x) --> S(x) = O S(x)
  if (source_operator_kind != UNIT) {  // modify source_matrix_p
    doSourceOperator(lat, whr, source_operator_kind);
#ifdef DEBUG_W_QUARK
    printf("finished source_operator_kind %i  \n",source_operator_kind);
#endif
  }
  
#ifdef DEBUG_W_QUARK
  // printf("\n");
  //printf("SOURCE MATRIX\n");
  //printf("=============\n");
  //DisplayAllSource(whr);
  //	printf("finished source set -- dump source \n");
  //	dumpSource("source.dat", source);
  //	int local_site[LORENTZs]={0,0,0,0};
  //	CheckSU3((const Float *)lat.GetLink(local_site, 1));
#endif 
  
  // now the source has been set 

  // ==================================================

  // initial guess for d_data_p = quark propagator [ N*COLORs*COLORs*DIRACs*DIRACs*COMPLEXs]
  {
    for (int i = 0; i < d_size;i++ )  d_data_p[i] = 0.;
    if(whr.onNode()){
      int site[4];
      
      for(site[0]=0;site[0]<lcl_sites[0];site[0]++){
	for(site[1]=0;site[1]<lcl_sites[1];site[1]++){
	  for(site[2]=0;site[2]<lcl_sites[2];site[2]++){
	    for(site[3]=0;site[3]<lcl_sites[3];site[3]++){
	      if(site[prop_direction]!=src_plane_position) continue;

	      for (int Cx = 0; Cx < COLORs; Cx++) {
		for (int Cy = 0; Cy < COLORs; Cy++) {
		  int src_offset = 
		    COMPLEXs*(Cx+COLORs*(Cy+COLORs*siteOffset(site, prop_dir)));

		  //set Dx=Dy  to source
		  for (int Dx = 0; Dx < DIRACs; Dx++) {
		    int prop_offset = 
		      COMPLEXs*Cx+
		      COMPLEXs*COLORS*Dx+
		      SPINORs*siteOffset(site)+
		      Cy*weightSrcColor()+
		      Dx*weightSrcDirac(); //Dy=Dx

		    d_data_p[prop_offset]  = rs_fac*source_matrix_p[src_offset];
		    d_data_p[prop_offset+1]= rs_fac*source_matrix_p[src_offset+1];

		  } // endfor (int Dx ...)
		}// endfor (int Cy ...)
	      } // endfor (int Cx ...)
	    }
	  }
	}
      } // endfor (int site[0] ...)
    } // endif (whr.onNode())
  } // end initial guess  
  // =============================================================

  // CHANGES COMPARED TO PHYS_V4.0.0
  // 1. initial guess has already been done above and is removed
  //    from the Cy, Dy loop below
  // 2. phys_v4.0.0 distinguishes between point, box and wall sources
  //    in the loop below, here we adopt phys_v3.11.4.xiaodong and
  //    set a general source defined above
  


  // start loop over source colours Cy and source Diracs Dy
  for (int Cy = 0; Cy < COLORs; ++Cy) {                 // source color
    for (int Dy = 0; Dy < DIRACs; ++Dy) {               // source dirac
      
      //Note: setSource must be passed zero src_matrix pointer if off-node!
      //Because src_lower/upper_bound[prop_dir] set by whr.lclMax/Min() 
      //is out of range if off-node
      {
	const IFloat *src_tmp=0;
	if(whr.onNode()){
	  src_tmp=source_matrix_p;
	}
	// source is a fermion vector [ N*COLORs*DIRACs*COMPLEXs ]
	// here it is initialised with the appropriate elements of
	// source_matrix_p
	source.setSource(Cy, Dy, rs_fac, src_tmp, prop_dir,
			 src_lower_bound, src_upper_bound);
      }
#ifdef DEBUG_W_WAVEFUNC
	if(src_op_kind==UNIT){
         source.printWaveFunc(wavefilename);
        }else{	
         source.printWaveFunc("wavefunc.dat.q2");
        }
	//source.print("source.dat");
#endif
      
#ifdef DEBUG_W_QUARK
	// printf("set and project sources for Cy, Dy = %4i,%4i\n",Cy,Dy);
	// dumpSource("source.dat", source);
#endif 
      
      // Do inversion seperately for each COLOR and SPIN
      //---------------------------------------------------------------------
      int f_size = GJP.VolNodeSites() * SPINORs;

      int          cg_iter;          // num of iter it took CG to converge
      Float        true_res;

      if (lat.Fclass() == F_CLASS_DWF) {
        int ls = GJP.SnodeSites();
        int ls_glb = GJP.SnodeSites() * GJP.Snodes();

	//If midplane is negative, use canonical location for the midpoint, i.e., Ls/2
        if ( midplane <= 0 ) midplane = ls_glb/2; 
  
        int f_size_5d = f_size * ls;
        Vector *src_4d = (Vector *)source.data();
        Vector *sol_4d = (Vector *)(d_data_p+f_size*(Cy+COLORs*Dy));

        Vector *src_5d = (Vector *) smalloc(f_size_5d * sizeof(Float));
        if(src_5d == 0)
          ERR.Pointer(d_class_name,ctor_str, "src_5d");
        VRB.Smalloc(d_class_name,ctor_str, "src_5d", src_5d,
                    f_size_5d * sizeof(Float));

        Vector *sol_5d = (Vector *) smalloc(f_size_5d * sizeof(Float));
        if(sol_5d == 0)
          ERR.Pointer(d_class_name,ctor_str, "sol_5d");
        VRB.Smalloc(d_class_name,ctor_str, "sol_5d", sol_5d,
                    f_size_5d * sizeof(Float));
	
        lat.Ffour2five(src_5d, src_4d, 0, ls_glb-1);
        lat.Ffour2five(sol_5d, sol_4d, ls_glb-1, 0);
	
        cg_iter = lat.FmatInv(sol_5d, src_5d, &(cg), &true_res);
	
	// ================= added from phys_v4.0.0 =================
        // Calculate <A_0 P> correlators for DWF lattice
        if (a0_p_outfile != 0)
          axialCurrent.measureAll(sol_5d);
	// ================= end added from phys_v4.0.0 ===================

        lat.Ffive2four(sol_4d, sol_5d, ls_glb-1, 0);

	//----------------------------------------------------------
	//Wall or Box sink -- added by mflin 03/16/06
	//---------------------------------------------------------
	if (warg.sink_kind == W_WALL || warg.sink_kind == W_BOX){
	  int snk_box_b[LORENTZs],snk_box_e[LORENTZs];
	  if (warg.sink_kind == W_WALL ){
	    for( int n = 0; n < LORENTZs; n++ ){
	      snk_box_b[n] = 0;
	      snk_box_e[n] = glb_sites[n] - 1;
	    }
	    snk_box_b[prop_dir] = snk_box_e[prop_dir] = 0;
	  }
	  else {
	    for( int n = 0; n < LORENTZs; n++ ){
	      snk_box_b[n] = warg.snk_box_b[n];
	      snk_box_e[n] = warg.snk_box_e[n];
	    }
	    snk_box_b[prop_dir] = snk_box_e[prop_dir] = 0;
	  }

	  FermionVector fv_sol_4d((Float *)sol_4d);
	  fv_sol_4d.gaugeFixSink(lat,prop_dir);

	  if(warg.sink_kind == W_BOX && warg.zero_mom_box_snk)
	    fv_sol_4d.sumOverHyperPlaneZeroMom(prop_dir,snk_box_b,snk_box_e);
	  else
	    fv_sol_4d.sumOverHyperPlane(prop_dir,snk_box_b,snk_box_e);
	}
	//-----------------------------
	//end added by mflin
	//-----------------------------

	//--------------------------------------------------------
	//Calculate midpoint correlator
	//midpoint is not necessarily at Ls/2
	//added by mflin 03/16/06
	//--------------------------------------------------------
        if (mid_point_outfile != 0){
          Vector *sol_4d_mid_point = 0;
	  sol_4d_mid_point =
	    (Vector *)(d_data_mid_point_1+f_size*(Cy+COLORs*Dy));
	  lat.Ffive2four(sol_4d_mid_point, sol_5d, midplane-1, midplane);
	  sol_4d_mid_point =
	    (Vector *)(d_data_mid_point_2+f_size*(Cy+COLORs*Dy)); 
	  lat.Ffive2four(sol_4d_mid_point, sol_5d, ls_glb-midplane-1, ls_glb-midplane);
	}
	// ==== end added by mflin 03/16/06 ===== 

        sfree(d_class_name,ctor_str, "sol_5d", sol_5d);

        sfree(d_class_name,ctor_str, "src_5d", src_5d);

	
      }
      else {
#ifdef DEBUG_W_QUARK
	//	printf("w_quark: enter FmatInv Cy=%4i, Dy=%4i \n",Cy,Dy);
#endif 
	
        
	cg_iter=0; 
	Vector *sol_4d = (Vector *)(d_data_p+f_size*(Cy+COLORs*Dy));
	
	cg_iter += lat.FmatInv(sol_4d, (Vector *)source.data(), 
			       &(cg), &true_res);
	
	
	//----------------------------------------------------------
	//Wall or Box sink -- added by mflin 03/16/06
	//---------------------------------------------------------
	if (warg.sink_kind == W_WALL || warg.sink_kind == W_BOX){
	  int snk_box_b[LORENTZs],snk_box_e[LORENTZs];
	  if (warg.sink_kind == W_WALL ){
	    for( int n = 0; n < LORENTZs; n++ ){
	      snk_box_b[n] = 0;
	      snk_box_e[n] = glb_sites[n] - 1;
	    }
	    snk_box_b[prop_dir] = snk_box_e[prop_dir] = 0;
	  }
	  else {
	    for( int n = 0; n < LORENTZs; n++ ){
	      snk_box_b[n] = warg.snk_box_b[n];
	      snk_box_e[n] = warg.snk_box_e[n];
	    }
	    snk_box_b[prop_dir] = snk_box_e[prop_dir] = 0;
	  }
	  
	  FermionVector fv_sol_4d((Float *)sol_4d);
	  fv_sol_4d.gaugeFixSink(lat,prop_dir);
	  
	  if(warg.sink_kind == W_BOX && warg.zero_mom_box_snk)
	    fv_sol_4d.sumOverHyperPlaneZeroMom(prop_dir,snk_box_b,snk_box_e);
	  else
	    fv_sol_4d.sumOverHyperPlane(prop_dir,snk_box_b,snk_box_e);
	}
	//-----------------------------
	//end added by mflin
	//-----------------------------
#ifdef DEBUG_W_QUARK
	//printf("w_quark: exit FmatInv %i %22.12e \n",cg_iter,true_res);
#endif 

      } // end else (if not DWF)
      
      // Print out number of iterations and true residual.
      //---------------------------------------------------------------------
      VRB.Result(d_class_name,ctor_str,"CG iterations = %d %.10e\n",
		 cg_iter, Float(true_res));
      Fprintf(fp, "%d %.10e\n", cg_iter, Float(true_res));

      // ============== added from phys_v4.0.0 ======================
      if ( (warg.source_kind==POINT_W) && pbp_outfile ) {
        // Calculate pbp
        pbp += ((Vector *)(d_data_p+f_size*(Cy+COLORs*Dy))) ->
          ReDotProductGlbSum((Vector *) source.data(), f_size);

        // Calculate pbg5p = Tr[ PsiBar * Gamma5 * Psi ]

        lat.Gamma5(sol_g5, (Vector *)(d_data_p+f_size*(Cy+COLORs*Dy)),
                   GJP.VolNodeSites());
        pbg5p += sol_g5->ReDotProductGlbSum((Vector *) source.data(), f_size);
      }
      // ============== end added from phys_v4.0.0 ==================


    } // end for (int Dy ...)
  } // end for (int Cy ... )

  //end of loop over colours and spins



  // ======= added from phys_v4.0.0
  if ( (warg.source_kind==POINT_W) && pbp_outfile ) {
    VRB.Sfree(d_class_name,ctor_str, "sol_g5", sol_g5);
    sfree(sol_g5);
  }
  // ======= end added from phys_v4.0.0
    
  
  // Close the file
  //-------------------------------------------------------------------------
  Fclose(fp);
  
  // === added from phys_v4.0.0 =================================
  if ( (warg.source_kind==POINT_W) && pbp_outfile) {
    if (!(fp = Fopen(pbp_outfile, "a"))) {
      ERR.FileA(d_class_name, ctor_str, pbp_outfile);
    }
    Fprintf(fp, "%e %e %e\n",
            IFloat(cg.mass),
            IFloat(pbp/pbp_norm),
            IFloat(pbg5p/pbp_norm));
    Fclose(fp);
  }

  if (lat.Fclass() == F_CLASS_DWF && a0_p_outfile != 0) {
    axialCurrent.doSum();
    axialCurrent.print();
  }


  // DEBUG -- PING
#ifdef DEBUG_W_QUARK
  //dumpData("quark_prop0.dat");
#endif
}



// END OF CONSTRUCTOR


// OTHER MEMBER FUNCTIONS


//---------------------------------------------------------------------------
// void WspectQuark::dumpData() const 
//--------------------------------------------------------------------------- 
void
WspectQuark::dumpData(char *filename) const {

#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif

  FILE *fp;
  
  if (filename && (fp = Fopen(filename, "a"))) {
    for (int i = 0; i < d_size; ++i)
      Fprintf(fp, "%g\n", d_data_p[i]);
  } else {
    ERR.FileA(d_class_name, "dumpData", filename);
  }

  Fclose(fp);  
}





//---------------------------------------------------------------------------
// WspectQuark::~WspectQuark()
//---------------------------------------------------------------------------
WspectQuark::~WspectQuark()
{
  VRB.Func(d_class_name, dtor_str);

  if(d_data_mid_point_2 != 0) {
    VRB.Sfree(d_class_name, dtor_str, empty_str, d_data_mid_point_2);
    sfree(d_data_mid_point_2);
  }
  if(d_data_mid_point_1 != 0) {
    VRB.Sfree(d_class_name, dtor_str, empty_str, d_data_mid_point_1);
    sfree(d_data_mid_point_1);
  }
  // end added from phys_v4.0.0

  VRB.Sfree(d_class_name, dtor_str, empty_str, d_data_p);
  sfree(d_data_p);

  // De-allocate source_matrix_p; added by T&X March 30.
  VRB.Func(d_class_name, dtor_str);
  VRB.Sfree(d_class_name, dtor_str, empty_str, source_matrix_p);
  sfree(source_matrix_p);
}




//---------------------------------------------------------------------------
// void WspectQuark::dumpSource() const 
//--------------------------------------------------------------------------- 
void
WspectQuark::dumpSource(char *filename, FermionVector &source) const {
  //print out non-zero source components

#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
   
  FILE *fp;
  
  Float *pp = source.data();
  //Float sum;
  //printf("dump source to file  %s \n",filename);

  if (filename && (fp = Fopen(filename, "a"))) {
    for (int i = 0; i < GJP.VolNodeSites()* COLORs*DIRACs*COMPLEXs; i++){
      //      sum=pp[i];
      //glb_sum((Float *)(&sum));
      if (pp[i] != 0.0) {
	Fprintf(fp,"i = %4i, source = %16.8e \n",i,pp[i]);
      }
    }
  } else {
    ERR.FileA(d_class_name, "dumpSource", filename);
  }

  Fclose(fp);  
}





//Added by T&X for debugging of GetLink
//---------------------------------------------------------------------------
// void WspectQuark::CheckSU3() const 
//--------------------------------------------------------------------------- 
void 
WspectQuark::CheckSU3(const Float *mat) const {

  int icol;
  Float fcol=0.;
	
  //DisplayMatrix(mat, "CHECK gauge for SU(3)");
  printf("UNITARITY CHECK: \n");
  
  for (icol = 0; icol < 2*COLORs; icol++)
    fcol+=mat[icol]*mat[icol];
  printf("first raw     =  1 = %22.12e \n", fcol);
  fcol = 0.;
  for (icol = 2*COLORs; icol < 4*COLORs; icol++)
    fcol+=mat[icol]*mat[icol];
  printf("second raw     = 1 = %22.12e \n", fcol);
  fcol = 0.;
  for (icol = 4*COLORs; icol < 6*COLORs; icol++)
    fcol+=mat[icol]*mat[icol];
  printf("third raw     =  1 = %22.12e \n", fcol);
  
  fcol = 0.;
  for (icol = 0; icol < COLORs; icol++)
    fcol+=mat[6*icol]*mat[6*icol]+mat[6*icol+1]*mat[6*icol+1];
  printf("first  column =  1 = %22.12e \n", fcol);
	
  fcol = 0.;
  for (icol = 0; icol < COLORs; icol++)
    fcol+=mat[6*icol+2]*mat[6*icol+2]+mat[6*icol+3]*mat[6*icol+3];
  printf("second column =  1 = %22.12e \n", fcol);

  fcol = 0.;
  for (icol = 0; icol < COLORs; icol++)
    fcol+=mat[6*icol+4]*mat[6*icol+4]+mat[6*icol+5]*mat[6*icol+5];
  printf("third column  =  1 = %22.12e \n", fcol);
  
}
/*




// added by T&X for debugging
//---------------------------------------------------------------------------
// void WspectQuark::DisplayMatrix(...) const 
//--------------------------------------------------------------------------- 
void 
WspectQuark::DisplayMatrix(const Float *matrix, const char *text) const {

  int display=0;
  for (int i = 0; i < COLORs*COLORs*COMPLEXs; i++){
    if (matrix[i] != 0.) display=1;
  }

  // display matrix only if some elements are non-zero  
  if (display != 0){
    printf("DISPLAY MATRIX -- %s \n",text);
    int icol;
    for (icol = 0; icol  < COLORs; icol++)
      printf("(%14.6e, %14.6e)", matrix[2*icol],matrix[2*icol+1]);
    printf("\n");
    for (icol = COLORs; icol  < 2*COLORs; icol++)
      printf("(%14.6e, %14.6e)", matrix[2*icol],matrix[2*icol+1]);
    printf("\n");
    for (icol = 2*COLORs; icol  < 3*COLORs; icol++)
      printf("(%14.6e, %14.6e)", matrix[2*icol],matrix[2*icol+1]);
    printf("\n");
  }else{
    printf("DISPLAY MATRIX: all elements are zero -- %s \n", text);
  }
}
*/
//--------------------------------------------------------------------------- 
// int WspectQuark::devDirToPolDir(...) const 
//--------------------------------------------------------------------------- 
/*DEV1, DEV2, DEV3 
 *Depending on proppagation direction, the actual direction of derivatives are:
 *   Note: prop_dir are 0(X),1(Y),2(Z),3(T)
 *      propdir-> 0  1  2  3
 *       DEV1     Y  X  X  X
 *       DEV2     Z  Z  Y  Y
 *       DEV3     T  T  T  Z
 */
int WspectQuark::devDirToPolDir(int dev_dir, int prop_dir) const{
  char *fname="devDirToPolDir";
  switch(dev_dir)
    {
    case 1:
      if (prop_dir > 0 ) return 0; // normal 0=x
      if (prop_dir == 0 ) return 1; 
      break;
    case 2:
      if (prop_dir >  1 ) return 1; // normal 1=y
      if (prop_dir <= 1 ) return 2;
      break;
    case 3:
      if (prop_dir >  2 ) return 2; // normal 2=z
      if (prop_dir <= 2 ) return 3;
      break;
    default:
      ERR.General(d_class_name, fname,"Invalid dev_dir");
    }
  return -1; //error
}

//Added by T&X
//--------------------------------------------------------------------------- 
// void WspectQuark::doSourceOperator(...) const 
//--------------------------------------------------------------------------- 
void WspectQuark::doSourceOperator(Lattice &lat, const WspectHyperRectangle &whr, DEVOperatorKind src_op) const

{
  char *fname = "doSourceOperator";
  VRB.Func(d_class_name, fname);
  
  // allocate temporary array to store output from symmetricDerivative
  // for src_op== SUM_F, SUM_S_ANTISYM, SUM_S_SYM, SUM_F_S_ANSYM, SUM_F_S_SYM
  Float *orig;
  Float *sum;
  int i;
  //if source operator is a SUM operator
  if(src_op>=SUM_F && src_op<END_SUM_OP){
    orig = (Float *) smalloc(source_matrix_size* sizeof(Float));
    if (!orig) {
      ERR.Pointer(d_class_name, ctor_str, "orig");
      VRB.Smalloc(d_class_name, ctor_str, 	   
		  "orig", orig, source_matrix_size * sizeof(Float));
    }
    sum = (Float *) smalloc(source_matrix_size* sizeof(Float));
    if (!sum) {
      ERR.Pointer(d_class_name, ctor_str, "sum");
      VRB.Smalloc(d_class_name, ctor_str, 	   
		  "sum", sum, source_matrix_size * sizeof(Float));
    }
    
    //initialize sum and orig
    for (i=0; i < source_matrix_size ; i++){
      sum[i]=0.0;
      orig[i] = source_matrix_p[i];
    }
  }
  
  int polarisation;
  Float weight=1.; //can be adjusted for each operator(the actual contribution may have extra -1 sign)
  
  
  // apply operator onto source_matrix and add result to dev with weight
  //-------------------------------------------------------------------------
  
  //-------------------------------------------------------
  // 1. Extended Meson Operators(Derivative Operators)
  // 2. ExtendedBE meson Operators(Local Field Operators)
  //-------------------------------------------------------

  // 1. Extended Meson Operators(Derivative Operators)

  weight=1.0; //weight of this operator
  if (src_op == DEV1 || src_op == SUM_F || 
      src_op == SUM_F_S_ANTISYM || src_op == SUM_UNIT_F_S_ANTISYM  ) {
    
#ifdef DEBUG_W_QUARK
    printf("calculate: DEV1\n");
#endif
    //  if (prop_direction >  0 ) polarisation = 0; // normal 1=x
    //  if (prop_direction == 0 ) polarisation = 1;
    //  To save code space replaced by:
    polarisation=devDirToPolDir(1,prop_direction);
    symmetricDerivative(lat, whr, polarisation);
    if (src_op != DEV1) {
      AddToSum(sum, weight);   // add weighted source_matrix (from symmetricDerivative) to sum
      Equal(orig);             // set source_matrix to original one
    }
  }

  weight=1.0; //weight of this operator
  if (src_op == DEV2 ||  src_op == SUM_F || 
      src_op == SUM_F_S_ANTISYM || src_op == SUM_UNIT_F_S_ANTISYM  ) {
#ifdef DEBUG_W_QUARK
    printf("calculate: DEV2\n");
#endif
  
    polarisation=devDirToPolDir(2,prop_direction);
    symmetricDerivative(lat, whr, polarisation);

    if (src_op != DEV2) {
      AddToSum(sum, weight);   
      Equal(orig);             // set source_matrix to original one
    }
  }

  weight=1.0; //weight of this operator
  if (src_op == DEV3 || src_op == SUM_F || 
      src_op == SUM_F_S_ANTISYM || src_op == SUM_UNIT_F_S_ANTISYM ) {
#ifdef DEBUG_W_QUARK
    printf("calculate: DEV3\n"); 
#endif
    
    polarisation=devDirToPolDir(3,prop_direction);
    symmetricDerivative(lat, whr, polarisation);

    if (src_op != DEV3) {
      AddToSum(sum, weight); 
      Equal(orig);            // set source_matrix to original one
    }
  }
  
  //commented out, not quite worth doing, may consider doing
  //this when using SUM operator 
  
  weight=1.0; //weight of this operator
  if (src_op == DEV1DEV1 || src_op == SUM_S_DIAG ||  
      src_op == SUM_S_SYM_DIAG ) {
#ifdef DEBUG_W_QUARK
    printf("calculate: DEV11\n");
#endif
    polarisation=devDirToPolDir(1,prop_direction);
    symmetricDerivative(lat, whr, polarisation);
    // DisplayAllSource(whr);
    polarisation=devDirToPolDir(1,prop_direction);
    symmetricDerivative(lat, whr, polarisation);
    // DisplayAllSource(whr);
    if (src_op != DEV1DEV1) {
      AddToSum(sum, weight);  // add weighted source_matrix (from symmetricDerivative) to sum
      Equal(orig);            // set source_matrix to original one
    }
  }
  
  if (src_op == DEV2DEV2 || src_op == SUM_S_DIAG ||  
      src_op == SUM_S_SYM_DIAG ) {
#ifdef DEBUG_W_QUARK
    printf("calculate: DEV22\n");
#endif
    polarisation=devDirToPolDir(2,prop_direction);
    symmetricDerivative(lat, whr, polarisation);
    // DisplayAllSource(whr);
    polarisation=devDirToPolDir(2,prop_direction);
    symmetricDerivative(lat, whr, polarisation);
    // DisplayAllSource(whr);
    if (src_op != DEV2DEV2) {
      AddToSum(sum, weight);  // add weighted source_matrix (from symmetricDerivative) to sum
      Equal(orig);            // set source_matrix to original one
    }
  }
    
    
  if (src_op == DEV3DEV3 || src_op == SUM_S_DIAG ||  
      src_op == SUM_S_SYM_DIAG) {
#ifdef DEBUG_W_QUARK
    printf("calculate: DEV33\n");
#endif
   polarisation=devDirToPolDir(3,prop_direction);
   symmetricDerivative(lat, whr, polarisation);
   // DisplayAllSource(whr);
   polarisation=devDirToPolDir(3,prop_direction);
   symmetricDerivative(lat, whr, polarisation);
   // DisplayAllSource(whr);
   if (src_op != DEV3DEV3) {
     AddToSum(sum, weight);  // add weighted source_matrix (from symmetricDerivative) to sum
     Equal(orig);            // set source_matrix to original one
   }
  }
  
  
  
  weight=1.0; //weight of this operator
  if (src_op == DEV1DEV2 || src_op == SUM_S_SYM || src_op== SUM_S_SYM_DIAG || 
      src_op == SUM_S_ANTISYM || src_op == SUM_F_S_ANTISYM || 
      src_op == SUM_UNIT_F_S_ANTISYM ){
#ifdef DEBUG_W_QUARK
    printf("calculate: DEV12\n");
#endif
    polarisation=devDirToPolDir(2,prop_direction);
    symmetricDerivative(lat, whr, polarisation);
    // DisplayAllSource(whr);
    polarisation=devDirToPolDir(1,prop_direction);
    symmetricDerivative(lat, whr, polarisation);
    // DisplayAllSource(whr);
    if (src_op != DEV1DEV2) {
      AddToSum(sum, weight);  
      Equal(orig);            // set source_matrix to original one
    }
  }

  weight=1.0; //weight of this operator
  if (src_op == DEV2DEV1 || src_op == SUM_S_SYM ||  src_op== SUM_S_SYM_DIAG ||
      src_op == SUM_S_ANTISYM || src_op == SUM_F_S_ANTISYM || 
      src_op == SUM_UNIT_F_S_ANTISYM ){
#ifdef DEBUG_W_QUARK
    printf("calculate: DEV21\n");
#endif
    polarisation=devDirToPolDir(1,prop_direction);
    symmetricDerivative(lat, whr, polarisation);
    // DisplayAllSource(whr);
    polarisation=devDirToPolDir(2,prop_direction);
    symmetricDerivative(lat, whr, polarisation);
    // DisplayAllSource(whr);

    if (src_op != DEV2DEV1) {
      if (src_op == SUM_S_SYM || src_op== SUM_S_SYM_DIAG) {
	AddToSum(sum, weight);  
	Equal(orig);            // set source_matrix to original one
      }else{
	AddToSum(sum, -weight);  
	Equal(orig);            // set source_matrix to original one
      }
    }
    
  }
  
  weight=1.0; //weight of this operator
  if (src_op == DEV2DEV3 || src_op == SUM_S_SYM ||  src_op== SUM_S_SYM_DIAG ||
      src_op == SUM_S_ANTISYM || src_op == SUM_F_S_ANTISYM || 
      src_op == SUM_UNIT_F_S_ANTISYM ){
#ifdef DEBUG_W_QUARK
    printf("calculate: DEV23\n"); 
#endif
    polarisation=devDirToPolDir(3,prop_direction);
    symmetricDerivative(lat, whr, polarisation);
    // DisplayAllSource(whr);
    polarisation=devDirToPolDir(2,prop_direction);
    symmetricDerivative(lat, whr, polarisation);
    // DisplayAllSource(whr);
    if (src_op != DEV2DEV3) {
      AddToSum(sum, weight);  
      Equal(orig);            // set source_matrix to original one
    }
  }

  weight=1.0; //weight of this operator
  if (src_op == DEV3DEV2 || src_op == SUM_S_SYM ||  src_op== SUM_S_SYM_DIAG ||
      src_op == SUM_S_ANTISYM || src_op == SUM_F_S_ANTISYM || 
      src_op == SUM_UNIT_F_S_ANTISYM ){
#ifdef DEBUG_W_QUARK
    printf("calculate: DEV32\n");
#endif
    polarisation=devDirToPolDir(2,prop_direction);
    symmetricDerivative(lat, whr, polarisation);
    // DisplayAllSource(whr);
    polarisation=devDirToPolDir(3,prop_direction);
    symmetricDerivative(lat, whr, polarisation);
    // DisplayAllSource(whr);
    if (src_op != DEV3DEV2 ) {
      if (src_op == SUM_S_SYM || src_op== SUM_S_SYM_DIAG) {
	AddToSum(sum, weight);  
	Equal(orig);            // set source_matrix to original one
      }else{
	AddToSum(sum, -weight);  
	Equal(orig);            // set source_matrix to original one
      }
    }
    
  }

  weight=1.0; //weight of this operator
  if (src_op == DEV1DEV3 || src_op == SUM_S_SYM ||  src_op== SUM_S_SYM_DIAG ||
      src_op == SUM_S_ANTISYM || src_op == SUM_F_S_ANTISYM || 
      src_op == SUM_UNIT_F_S_ANTISYM ){
#ifdef DEBUG_W_QUARK
    printf("calculate: DEV13\n"); 
#endif
    polarisation=devDirToPolDir(3,prop_direction);
    symmetricDerivative(lat, whr, polarisation);
    // DisplayAllSource(whr);
    polarisation=devDirToPolDir(1,prop_direction);
    symmetricDerivative(lat, whr, polarisation);
    // DisplayAllSource(whr);
    if (src_op != DEV1DEV3 ) {
      AddToSum(sum, weight); 
      Equal(orig);            // set source_matrix to original one
    }
  }

  weight=1.0; //weight of this operator
  if (src_op == DEV3DEV1 || src_op == SUM_S_SYM ||  src_op== SUM_S_SYM_DIAG ||
      src_op == SUM_S_ANTISYM || src_op == SUM_F_S_ANTISYM || 
      src_op == SUM_UNIT_F_S_ANTISYM ){
#ifdef DEBUG_W_QUARK
    printf("calculate: DEV31\n"); 
#endif
    polarisation=devDirToPolDir(1,prop_direction);
    symmetricDerivative(lat, whr, polarisation);
    // DisplayAllSource(whr);
    polarisation=devDirToPolDir(3,prop_direction);
    symmetricDerivative(lat, whr, polarisation);
    // DisplayAllSource(whr);
    if (src_op !=DEV3DEV1 ) {
      if (src_op == SUM_S_SYM || src_op== SUM_S_SYM_DIAG) {
	AddToSum(sum, weight);  
	Equal(orig);            // set source_matrix to original one
      }else{
	AddToSum(sum, -weight);  
	Equal(orig);            // set source_matrix to original one
      }
    }
  }


  // 2. ExtendedBE meson Operators(Local Field Operators)
  /*FB1_OP=0,FB2_OP,FB3_OP,
   *FE1_OP,FE2_OP,FE3_OP,
   *FUNIT_OP,  //for mixing
   *
   *group 2:
   *SUM_MAGN_OP,
   *SUM_ELEC_OP,
   *group 3
   *SUM_MAGN_ELEC_OP,
   */
  
  //multiply source by field tensor at every site
  if(src_op>BEGIN_BE_OP && src_op<END_BE_OP){
    doBEOperator(src_op,whr);
  }


  weight=1.0; //weight of this operator
  //copy sum to source matrix
  if (src_op >= SUM_F && src_op<END_SUM_OP) {
#ifdef DEBUG_W_QUARK
    printf("sum --> source_matrix\n");
#endif
    Equal(sum);            // set source_matrix to sum
    // moveMem((Float *)source_matrix_p, (Float *)sum, source_matrix_size**sizeof(Float));
    VRB.Func(d_class_name, dtor_str);
    VRB.Sfree(d_class_name, dtor_str, empty_str, orig);
    sfree(orig);
    
    
    VRB.Func(d_class_name, dtor_str);
    VRB.Sfree(d_class_name, dtor_str, empty_str, sum);
    sfree(sum);
  }

  
 
}



// added by Thomas and Xiaodong
// calculates the symmetric derivative on a given source matrxi
//---------------------------------------------------------------------------
// void WspectQuark::symmetricDerivative(...) const 
//--------------------------------------------------------------------------- 
void 
WspectQuark::symmetricDerivative(Lattice &lat, const WspectHyperRectangle &whr, int polarisation) const

{
  char *fname = "symmetricDerivative";
  VRB.Func(d_class_name, fname);
  
  const Float *source_matrix;    //  source matrix at a given site
  Float product[COLORs*COLORs*COMPLEXs];
  Matrix dag;
  const Matrix*  gauge;  // SU(3) matrix at a given site into given direction 

  int C1,C2;


  //should do nothing if the node is not on the source plane at all!
  if(!whr.onNode()) return;

  Float *dev = (Float *) smalloc(source_matrix_size* sizeof(Float));
  if (!dev) {
    ERR.Pointer(d_class_name, ctor_str, "dev");
  }
  VRB.Smalloc(d_class_name, ctor_str, 	   
	      "dev", dev, source_matrix_size * sizeof(Float));


  int lcl_shift[LORENTZs]; 
  int lcl[LORENTZs];
  int min[LORENTZs];
  int max[LORENTZs];
  int offset;
  int productOffset;
  int devOffset;

  if (polarisation == prop_direction) {
    ERR.General(d_class_name, fname, "Error:deriv_dir==prop_dir");
    //printf("prop_direction = %4i, polarisation = %4i \n",prop_direction,polarisation);
  }

  /*Wrong! in prop_direction min=max=src_plane_position!!
    or getLink will get the link on the wrong plane
    max[0] = lcl_sites[0]-1;
    max[1] = lcl_sites[1]-1;
    max[2] = lcl_sites[2]-1;
    max[3] = lcl_sites[3]-1;
    
    // do not sum over the direction of propagation, normally "time" = 3
    max[prop_direction] = 0;
  */
  
  //changed to:
  for(int i=0;i<LORENTZs;i++){
    max[i]=lcl_sites[i]-1;
    min[i]=0;
  }
  
  min[prop_direction]=max[prop_direction]=src_plane_position;
  if(min[prop_direction]<0) ERR.General(d_class_name, fname, "wrong src plane");
  

  for (lcl[3] = min[3]; lcl[3] <= max[3]; lcl[3]++) {
    for (lcl[2] = min[2]; lcl[2] <= max[2]; lcl[2]++) {   
      for (lcl[1] = min[1]; lcl[1] <= max[1] ; lcl[1]++) {
	for (lcl[0] = min[0]; lcl[0] <= max[0]; lcl[0]++) {
	  
	  // printf("symdev lcl = %i%i%i%i \n",lcl[0],lcl[1],lcl[2],lcl[3]);
	  offset = siteOffset(lcl,prop_direction);

	  lcl_shift[0] = lcl[0];
	  lcl_shift[1] = lcl[1];
	  lcl_shift[2] = lcl[2];
	  lcl_shift[3] = lcl[3];

	  lcl_shift[polarisation]=lcl[polarisation]+1;
	  //if (polarisation == 0) lcl_shift[0] = lcl[0] + 1;
	  //if (polarisation == 1) lcl_shift[1] = lcl[1] + 1;
	  //if (polarisation == 2) lcl_shift[2] = lcl[2] + 1;
	  //if (polarisation == 3) lcl_shift[3] = lcl[3] + 1;


	  if (lcl_shift[prop_direction] != min[prop_direction])
	    {
	      ERR.General(d_class_name, fname, "shifted source not on src_slice!");
	    }

	  // get gauge matrix U_i(x) at lcl from lat.
	  if(src_fuzz_p==0){
	    gauge  = lat.GetLink(lcl, polarisation);
	  }else{
	  gauge = src_fuzz_p->GetLink(lcl, polarisation);
	  }
	  //DisplayMatrix((const Float *) gauge, "U_i(x) in symmderiv");
#ifdef DEBUG_CHECKSU3
	  CheckSU3((const Float *)gauge);
#endif
	  // get 3x3 "source matrix" from shifted source_matrix (3x3xsites)
	  // which lives only on a 3-dim hyperplane orthogonal to 
	  // the propagation direction, GetSource does the job

	  source_matrix = (const Float *)GetSource(lcl_shift, whr);
	  // DisplayMatrix(source_matrix, "S(x+i) in symm");
	 
	  // multiply U_i(x) * S(x+i)
	  mDotMEqual((IFloat *)product, (const IFloat *)gauge,(const IFloat *)source_matrix); 
	  // DisplayMatrix(product, "forward product of U_i(x)*S(x+i)");

	  // copy product onto dev(x)
	  for (C2 = 0; C2 < COLORs; C2++) {
	    for (C1 = 0; C1 < COLORs; C1++) {
	      productOffset    = (C1 + C2*COLORs)*COMPLEXs;
	      devOffset        = productOffset + offset*COLORS*COLORs*COMPLEXs;

	      dev[devOffset]   = product[productOffset];     // real
	      dev[devOffset+1] = product[productOffset+1];   // imaginary
	    }
	  }
	  
	  lcl_shift[polarisation]=lcl[polarisation]-1;
	  //if (polarisation == 0) lcl_shift[0] = lcl[0] - 1;
	  //if (polarisation == 1) lcl_shift[1] = lcl[1] - 1;
	  //if (polarisation == 2) lcl_shift[2] = lcl[2] - 1;
	  //if (polarisation == 3) lcl_shift[3] = lcl[3] - 1;
	  
	  // get U_i(x-i)

	  if(src_fuzz_p==0){
	    gauge  = lat.GetLink(lcl_shift, polarisation);
	  }else{
	    gauge  = src_fuzz_p->GetLink(lcl_shift, polarisation);
	  }
 	  // DisplayMatrix((const Float *)gauge, "U_i(x-i) in symmderiv");
#ifdef DEBUG_CHECKSU3
	  CheckSU3((const Float *)gauge);
#endif
	  //  U_i(x-i) = U_i^{\dag}(x-i)
	  // couldn't get dag to work
	  dag.Dagger((const IFloat*)gauge); 
	  //Dagger(dag,(const Float *)gauge);
 	  // DisplayMatrix((const Float*)dag, "dag in symmderiv");

	  // get S(x-i)
	  source_matrix = (const Float *)GetSource(lcl_shift, whr);
 	  // DisplayMatrix(source_matrix, "S(x-i) in symmderiv");
	  
	  // multiply  U_i^{\dag}(x-i) S(x-i)
	  mDotMEqual((IFloat *)product, (const IFloat *)&dag,(const IFloat *)source_matrix);
	  
 	  // DisplayMatrix(product, "backward product of U_{-i}(x)*S(x-i)");
 
	  // subtract product from dev(x)
	  for (C2 = 0; C2 < COLORs; C2++) {
	    for (C1 = 0; C1 < COLORs; C1++) {
	      productOffset    = (C1 + C2*COLORs)*COMPLEXs;
	      devOffset        = productOffset + offset*COLORS*COLORs*COMPLEXs;
	      
	      dev[devOffset]   -= product[productOffset];     // real
	      dev[devOffset+1] -= product[productOffset+1];   // imaginary
	    }
	  }  
	}
      }
    }
  }

  //copy the result into source_matrix_p 
  moveMem((Float *)source_matrix_p, (Float *)dev, source_matrix_size*sizeof(Float));
  VRB.Func(d_class_name, dtor_str);
  VRB.Sfree(d_class_name, dtor_str, empty_str, dev);
  sfree(dev);
}


//multiply src matrix by field tensor
//-------------------------------------------
// doBEOperator(DEVOperatorKind src_op)
//-------------------------------------------
void WspectQuark::doBEOperator(DEVOperatorKind src_op, const WspectHyperRectangle &whr)const{
  

  const Matrix *src_mat_p;    //  source matrix at a given site
  
  Matrix field_mat, product;

  int lcl[LORENTZs];
  int min[LORENTZs];
  int max[LORENTZs];
  //loop over all sites
  for(int i=0;i<LORENTZs;i++){
    max[i]=lcl_sites[i]-1;
    min[i]=0;
  }
  
  min[prop_direction]=max[prop_direction]=src_plane_position;
  
  for (lcl[3] = min[3]; lcl[3] <= max[3]; lcl[3]++) {
    for (lcl[2] = min[2]; lcl[2] <= max[2]; lcl[2]++) {   
      for (lcl[1] = min[1]; lcl[1] <= max[1] ; lcl[1]++) {
	for (lcl[0] = min[0]; lcl[0] <= max[0]; lcl[0]++) {
	  
	  // printf("symdev lcl = %i%i%i%i \n",lcl[0],lcl[1],lcl[2],lcl[3]);
	  //original src_matrix element:
	  src_mat_p = GetSource(lcl, whr);
	  //multiply by field
	  switch(src_op){
	  case FB1_OP:
	    field_mat=fld_p->getField(lcl,FB1);
	    break;
	  case FB2_OP:
	    field_mat=fld_p->getField(lcl,FB2);
	    break;
	  case FB3_OP:
	    field_mat=fld_p->getField(lcl,FB2);
	    break;
	  case FE1_OP:
	    field_mat=fld_p->getField(lcl,FE1);
	    break;
	  case FE2_OP:
	    field_mat=fld_p->getField(lcl,FE2);
	    break;
	  case FE3_OP:
	    field_mat=fld_p->getField(lcl,FE2);
	    break;
	  case SUM_MAGN_OP:
	    field_mat=fld_p->getField(lcl,FB1);
	    field_mat+=fld_p->getField(lcl,FB2);
	    field_mat+=fld_p->getField(lcl,FB2);
	    break;
	  case SUM_ELEC_OP: 
	    field_mat=fld_p->getField(lcl,FE1);
	    field_mat+=fld_p->getField(lcl,FE2);
	    field_mat+=fld_p->getField(lcl,FE2);
	    break;
	  case SUM_MAGN_ELEC_OP:
	    field_mat=fld_p->getField(lcl,FB1);
	    field_mat+=fld_p->getField(lcl,FB2);
	    field_mat+=fld_p->getField(lcl,FB2);
	    field_mat+=fld_p->getField(lcl,FE1);
	    field_mat+=fld_p->getField(lcl,FE2);
	    field_mat+=fld_p->getField(lcl,FE2);
	    break;
	  default:
	    break;
	    //write back to src matrix
	  } //end switch
	  
	  
	  mDotMEqual((IFloat *)&product, (const IFloat *)&field_mat,(const IFloat *)src_mat_p);
	  moveMem((Float *)src_mat_p, (Float *)&product, sizeof(Matrix));	  

	}
      }
    }
  }

}


// added by Thomas and Xiaodong
// calculates S(x) = (1- r*Delta^2) S(x), where S is a given source matrix
// S(x) = S(x)+ \sum_i [  2*r S(x) - r*S(x+i) -r*S(x-i) ]
//      = (1+ 3*2*r) S(x) + \sum_i [ r*S(x+i) -r*S(x-i) ]
//Can be called several times
//---------------------------------------------------------------------------
// void WspectQuark::JacobiSource(...) const 
//--------------------------------------------------------------------------- 
void 
WspectQuark::JacobiSource(Lattice &lat, const WspectHyperRectangle &whr, Float r) const

{

  //should do nothing if the node is not on the source plane at all!
  //?? Might cause different memory map for different nodes??
  if(!whr.onNode()) return;

  char *fname = "JacobiSource";
  VRB.Func(d_class_name, fname);
  const Float *source_matrix;    //  source matrix at a given site
  Float product[COLORs*COLORs*COMPLEXs];
  
  Matrix dag;
  const Matrix*  gauge;  // SU(3) matrix at a given site into given direction 
 
  int C1,C2;

  Float *jacob = (Float *) smalloc(source_matrix_size* sizeof(Float));
  if (!jacob) {
    ERR.Pointer(d_class_name, ctor_str, "jacob");
    VRB.Smalloc(d_class_name, ctor_str, 	   
		"jacob", jacob, source_matrix_size * sizeof(Float));
  }

  int lcl_shift[LORENTZs]; 
  int lcl[LORENTZs];
  int min[LORENTZs];
  int max[LORENTZs];
  int offset; 
  int productOffset;
  int jacobOffset;
  int i;

  // initialise jacob before looping over directions
  for (i=0; i < source_matrix_size; i++){
    jacob[i]=0.;
  }

  for(i=0;i<LORENTZs;i++){
    max[i] = lcl_sites[i]-1;
    min[i]=0;
  }

  //for propagtion direction max=min=srcplane_local location
  //printf("prop_direction=%d, src_plane=%d\n",prop_direction,src_plane_position);
  min[prop_direction]=max[prop_direction]=src_plane_position;
  if(min[prop_direction]<0) ERR.General(d_class_name, fname, "src plane<0");
  
  
  // start loop over directions, ignore the prop_direction
  for (int polarisation=0; polarisation < LORENTZs; polarisation++){

    if (polarisation == prop_direction) {
      //      printf("JACOBI - error !\n");
      //      printf("prop. dir. and source polarisation must never be the same\n");
      //      printf("prop_direction = %4i, polarisation = %4i \n",prop_direction,polarisation);

    } else {

      
      // do not sum over the direction of propagation, normally "time" = 3
 
       

      for (lcl[0] = min[0]; lcl[0] <= max[0]; lcl[0]++) {
	for (lcl[1] = min[1]; lcl[1] <= max[1]; lcl[1]++) {   
	  for (lcl[2] = min[2]; lcl[2] <= max[2] ; lcl[2]++) {
	    for (lcl[3] = min[3]; lcl[3] <= max[3]; lcl[3]++) {

	      //	      printf("jacobi lcl = %i,%i,%i,%i \n",lcl[0],lcl[1],lcl[2],lcl[3]);
	      offset = siteOffset(lcl,prop_direction);

	      lcl_shift[0] = lcl[0];
	      lcl_shift[1] = lcl[1];
	      lcl_shift[2] = lcl[2];
	      lcl_shift[3] = lcl[3];

	      /*
	      if (polarisation == 0) lcl_shift[0] = lcl[0] + 1;
	      if (polarisation == 1) lcl_shift[1] = lcl[1] + 1;
	      if (polarisation == 2) lcl_shift[2] = lcl[2] + 1;
	      if (polarisation == 3) lcl_shift[3] = lcl[3] + 1;
	      */
	      lcl_shift[polarisation] = lcl[polarisation] + 1;

	      if (lcl_shift[prop_direction] != min[prop_direction])
		{
		  ERR.General(d_class_name, fname, "JACOBI -error !\n");
		}

	      // get gauge matrix U_i(x) at lcl from lat.
	      gauge  = lat.GetLink(lcl, polarisation);
	      
	      // get 3x3 "source matrix" from shifted source_matrix S(x+i) 
	      source_matrix = (const Float *)GetSource(lcl_shift, whr);
	      
	      // multiply U_i(x) * S(x+i)
	      mDotMEqual((IFloat *)product, (const IFloat *)gauge,(const IFloat *)source_matrix); 

	      // subtract product from jacob(x) -= r*U_i(x)S(x+i)
	      for (C1 = 0; C1 < COLORs; C1++) {
		for (C2 = 0; C2 < COLORs; C2++) {
		  productOffset    = (C1 + C2*COLORs)*COMPLEXs;
		  jacobOffset        = productOffset + offset*COLORS*COLORs*COMPLEXs;

		  jacob[jacobOffset]   += r*product[productOffset];     // real
		  jacob[jacobOffset+1] += r*product[productOffset+1];   // imaginary
		}
	      }
	      /*
	      if (polarisation == 0) lcl_shift[0] = lcl[0] - 1;
	      if (polarisation == 1) lcl_shift[1] = lcl[1] - 1;
	      if (polarisation == 2) lcl_shift[2] = lcl[2] - 1;
	      if (polarisation == 3) lcl_shift[3] = lcl[3] - 1;
	      */
	      lcl_shift[polarisation]=lcl[polarisation]-1;

	      // get U_i(x-i)
	      gauge  = lat.GetLink(lcl_shift, polarisation);
	      //  U_i(x-i) = U_i^{\dag}(x-i)
	      // dag.Dagger(gauge); 
	      //Dagger(dag,(const Float *)gauge);
	      dag.Dagger((const IFloat *)gauge);
	      // DisplayMatrix((const Float*)dag, "dag in symmderiv");

	      // get S(x-i)
	      source_matrix = (const Float *)GetSource(lcl_shift, whr);
	  
	      // multiply  U_i^{\dag}(x-i) S(x-i)
	      mDotMEqual((IFloat *)product, (const IFloat *)&dag,(const IFloat *)source_matrix); 	
 
	      // subtract product from jacob(x): jacob-= r*S(x-i)
	      for (C1 = 0; C1 < COLORs; C1++) {
		for (C2 = 0; C2 < COLORs; C2++) {
		  productOffset    = (C1 + C2*COLORs)*COMPLEXs;
		  jacobOffset        = productOffset + offset*COLORS*COLORs*COMPLEXs;
		  
		  jacob[jacobOffset]   += r*product[productOffset];     // real
		  jacob[jacobOffset+1] += r*product[productOffset+1];   // imaginary
		}
	      }
	    } // end of lcl[3]
	  } // end of lcl[2]
	} // end of lcl[1]
      } // end of lcl[0]
    } // end if polarisation == whr.dir
  } // end of loop over polarisations

  // add central value 3*(1 + 2*r) S(x), factor 3 for directions

  for (i=0; i < source_matrix_size; i++){
    jacob[i] += (1-6*r)*source_matrix_p[i];
    //jacob[i] += source_matrix_p[i];
  }

  // copy jacob --> S(x)
  moveMem((Float *)source_matrix_p, (Float *)jacob, source_matrix_size*sizeof(Float));
  
  // free jacob
  VRB.Func(d_class_name, dtor_str);
  VRB.Sfree(d_class_name, dtor_str, empty_str, jacob);
  sfree(jacob);
}






// added by Thomas and Xiaodong:
// get the source at site defined by *site
// Similar to Lattice::GetLink(...)
// result is a complex 3x3 matrix = 18 Floats
// assumes that off-nodes are NEXT-NODES ?
//---------------------------------------------------------------------------
// void WspectQuark::GetSource(...) const 
//--------------------------------------------------------------------------- 
const Matrix *
WspectQuark::GetSource(const int *site, const WspectHyperRectangle &whr) const
{
  // offset out-of-range coordinates site[] into on_node_site[]
  // in order to locate the link
  //-----------------------------------------------------

  int on_node_site[4];
  int on_node = 1;
  const Matrix *on_node_source;

  int node_sites[LORENTZs];

  node_sites[0]=lcl_sites[0];
  node_sites[1]=lcl_sites[1];
  node_sites[2]=lcl_sites[2];
  node_sites[3]=lcl_sites[3];
 

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

    int sourceOffset = siteOffset(on_node_site,prop_direction)*COMPLEXs*COLORs*COLORs;
    // checked projection to on_node_site and sourceOffset

    on_node_source = (Matrix *)(source_matrix_p + sourceOffset);
    // DisplayMatrix((Float *)on_node_source, "on_node_source");
  }


#ifndef PARALLEL
  //VRB.FuncEnd(cname, fname) ;
  return on_node_source;
#endif

  // send to the destination node if the site is off-node
  //------------------------------------------------------------------------ 

  if (on_node) {
    //  VRB.FuncEnd(cname, fname) ;
    return on_node_source;
  } else {
    Matrix send = *on_node_source;
    Matrix &recv = m_tmp1 ;
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
/*
// calculate dagger of complex matrix
// very primitive, but I couldn't get dag-> to work
// Thomas&Xiaodong, April 1st.

void
WspectQuark::Dagger(IFloat * dagger, const IFloat * matrix) const
{
  //   (0,1)     (2,3)     (4,5)
  //   (6,7)     (8,9)   (10,11)
  // (12,13)   (14,15)   (16,17)

  dagger[0] =  matrix[0];
  dagger[1] = -matrix[1];

  dagger[2] =  matrix[6];
  dagger[3] = -matrix[7];

  dagger[4] =  matrix[12];
  dagger[5] = -matrix[13];

  dagger[6] =  matrix[2];
  dagger[7] = -matrix[3];

  dagger[8] =  matrix[8];
  dagger[9] = -matrix[9];

  dagger[10] =  matrix[14];
  dagger[11] = -matrix[15];

  dagger[12] =  matrix[4];
  dagger[13] = -matrix[5];

  dagger[14] =  matrix[10];
  dagger[15] = -matrix[11];

  dagger[16] =  matrix[16];
  dagger[17] = -matrix[17];
}
*/  
//Added by T&X
//copy *src_mat to *WspectQuark::source_matrix_p (which is passed to FmatInv)
//---------------------------------------------------------------------------
// void WspectQuark::copyGaugefixToSource(...) const 
//--------------------------------------------------------------------------- 
void 
WspectQuark::copyGaugefixToSource(const Float *src_mat, const WspectHyperRectangle &whr) const

{
  char *fname = "copyGaugefixToSource()";
  VRB.Func(d_class_name, fname);
  int max[LORENTZs], min[LORENTZs];
  int prop_dir = prop_direction;

 
  max[0] = lcl_sites[0]-1;
  max[1] = lcl_sites[1]-1;
  max[2] = lcl_sites[2]-1;
  max[3] = lcl_sites[3]-1;

  min[0] = 0;
  min[1] = 0;
  min[2] = 0;
  min[3] = 0;

  // no loop over start slice, fix it to start slice coordinate
  min[prop_dir] = src_plane_position;
  max[prop_dir] = src_plane_position;


  // map Dagger(src_matrix) --> source_matrix_p 
  //-------------------------------------------------------------------------
  int lcl[LORENTZs];                             // local lattice site
  
  for (lcl[0] = min[0]; lcl[0] <= max[0]; lcl[0]++) {
    for (lcl[1] = min[0]; lcl[1] <= max[1]; lcl[1]++) {   
      for (lcl[2] = min[0]; lcl[2] <= max[2] ; lcl[2]++) {
	for (lcl[3] = min[0]; lcl[3] <= max[3]; lcl[3]++) {
	  for (int C1 = 0; C1 < COLORs; C1++) {
	    for (int C2 = 0; C2 < COLORs; C2++) {
	      int old_offset = COMPLEXs*(C1+COLORs*(C2+COLORs*siteOffset(lcl, prop_dir)));
	      int new_offset = COMPLEXs*(C2+COLORs*(C1+COLORs*siteOffset(lcl, prop_dir)));
	      source_matrix_p[new_offset]   =  src_mat[old_offset];        // real
	      //Negative sign for imaginary part is due to inverse Link Matrix
	      source_matrix_p[new_offset+1] = -src_mat[old_offset+1];      // imag
	    }
	  }
	}
      }
    }
  }
}




//-------------------------------------------------------
void 
WspectQuark::Equal(Float *orig) const 
{

  for (int i = 0; i < source_matrix_size; i++){
    source_matrix_p[i] = orig[i];
  }
  //can use moveMem
  //moveMem(src_matrix_p,orig,source_matrix_size*sizeof(IFloat);

}
//-------------------------------------------------------
void 
WspectQuark::AddToSum(Float *sum, Float weight) const

{
  for (int i = 0; i < source_matrix_size; i++){
    sum[i] += weight*source_matrix_p[i];
  }
}

/*

//-------------------------------------------------------
void
WspectQuark::DisplayAllSource(const WspectHyperRectangle &whr) const
{
  int lcl[LORENTZs];
  for (lcl[0] = 0; lcl[0] < GJP.XnodeSites() ; lcl[0]++) {
    for (lcl[1] = 0; lcl[1] < GJP.YnodeSites(); lcl[1]++) {   
      for (lcl[2] = 0; lcl[2] < GJP.ZnodeSites(); lcl[2]++) {
	for (lcl[3] = 0; lcl[3] < GJP.TnodeSites(); lcl[3]++) {
	  int i = siteOffset(lcl,prop_direction)*COMPLEXs*COLORs*COLORs;
	  printf("check source lcl = %i%i%i%i, i = %4i \n",lcl[0],lcl[1],lcl[2],lcl[3],i);
	  DisplayMatrix((Float*)source_matrix_p + i, "check source");
	}
      }
    }
  }
}

*/

CPS_END_NAMESPACE
