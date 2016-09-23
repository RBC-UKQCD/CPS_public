#include<config.h>
CPS_START_NAMESPACE
/*!\file

  $Id: w_ext_mesons.C,v 1.11 2008/02/08 18:35:05 chulwoo Exp $
*/
/* class WspectExtendedMesons
 * Thomas and Xiaodong. March 2000.
 * calculates all correlators from a given set of basis operators
 * according to a weighting table (correlator, basis operator)
 * basis operators are arbitrary combinations of GAMMAS and DERIVATIVES
 * Since these basis operators are defined both at source and
 * sink, the second index in the table runs to basis_operators*basis_operators
 * this allows also a mixing study between different operators
 * in class WspectMesons 
 * Including: magnetic hybrids, ....
 */

// #define MIXING_ON


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
#include <comms/cbuf.h> 
CPS_START_NAMESPACE

//---------------------------------------------------------------------------
// For the purpose of debugging or timing during code upgrade
//---------------------------------------------------------------------------
//#define DEBUG_W_EXT_MESON

//#define TIME_W_EXTENDED_MESON

#ifdef  TIME_W_EXTENDED_MESON
CPS_END_NAMESPACE
#include <time.h>                 // clock()
CPS_START_NAMESPACE
#endif




//---------------------------------------------------------------------------
// static data members
//---------------------------------------------------------------------------
char *WspectExtendedMesons::d_class_name = "WspectExtendedMesons";
int WspectExtendedMesons::tableInitialized=0;
struct WMesonOpInfo
        WspectExtendedMesons::WMesonOpTable[NUM_WMESON_OP_KIND];
struct WMesonStateInfo
        WspectExtendedMesons::WMesonStateTable[NUM_WMESON_STATE];

//---------------------------------------------------------------------------
// WspectExtendedMesons::WspectExtendedMesons(...)
//--------------------------------------------------------------------------- 
WspectExtendedMesons::WspectExtendedMesons(WspectArg *w_arg_p,const WspectHyperRectangle & whr, int fuzzing_index, int allocatemem)
  : arg_p(w_arg_p), fuzzing_c_index(fuzzing_index), d_whr(whr)
{
  
  //const members have been initialized in the initialization list
  d_prop_dir  = whr.dir();  
  d_glb_walls = glb_sites[d_prop_dir];  
  d_size   = d_glb_walls* NUM_WMESON_STATE * COMPLEXs; //correlator array in Float
  d_output_size   = d_glb_walls* NUM_WMESON_OUTPUT * COMPLEXs; //corrlator(polarisation average)in Float
  {
    const int *low  = whr.lclMin();
    const int *high = whr.lclMax();
    for (int i = 0; i < LORENTZs; ++i) {
      d_lclMin[i] = low[i];
      d_lclMax[i] = high[i];
    }
  }  
  
  //allocate memory for coor_data_p and zero out the array
  //dont't forget to release memory in destructor!
  if(allocatemem){
    coor_data_p = (Complex *) smalloc(d_size * sizeof(Float));
    if (!coor_data_p) 
      ERR.Pointer(d_class_name, ctor_str, "coor_data_p");
    VRB.Smalloc(d_class_name, ctor_str, 
		"coor_data_p", coor_data_p, d_size * sizeof(Float));
    {
      //zero out correlators
      Float *flt_p = (Float *)coor_data_p;
      for (int i = 0; i < d_size; ++i)
	*flt_p++ = 0.0;
    }
    
    //allocate memory for coor_output_p and zero out the array
  //dont't forget to release memory in destructor!
    coor_output_p = (Complex *) smalloc(d_output_size * sizeof(Float));
    if (!coor_output_p) 
      ERR.Pointer(d_class_name, ctor_str, "coor_output_p");
    VRB.Smalloc(d_class_name, ctor_str, 
		"coor_output_p", coor_output_p, d_output_size * sizeof(Float));
    {
      //zero out correlators
      Float *flt_p = (Float *)coor_output_p;
      for (int i = 0; i < d_output_size; ++i)
	*flt_p++ = 0.0;
    }
    
  }

  {
    //initialize map
    //!!!!!!!!!!!!!!!!!!!!!Need to be modified to make prop!=3 work!!!!!!!
    //for now assume prop_dir==3 only
    //should do manually mapping from coor_data_p (calculated aussuming prop==3)
    //to output data files later
    for(int i=0;i<NUM_WMESON_STATE;i++){
      map[i]=i; //if prop_dir=3(T)
    }
  }

  //initialize the stateTable according to spect argument
  if(!tableInitialized){
    initWMesonOpTable();
    initWMesonStateTable(arg_p);
    tableInitialized=1;
  }
#ifdef DEBUG_W_EXT_MESON
  //printf("ExtendedMeson constructor done.\n");
  //printf("prop_dir=%d,glb_walls=%d,d_size=%d,data_p=%x\n",d_prop_dir,d_glb_walls,d_size,coor_data_p);
#endif
  
}

//---------------------------------------------------------------------------
// WspectExtendedMesons::~WspectExtendedMesons(...)
//--------------------------------------------------------------------------- 
WspectExtendedMesons::~WspectExtendedMesons(){
  //free memory
  
  VRB.Func(d_class_name, dtor_str);
  if(coor_data_p){
    VRB.Sfree(d_class_name, dtor_str, empty_str, coor_data_p);
    sfree(coor_data_p);
  }
  if(coor_output_p){
    VRB.Sfree(d_class_name, dtor_str, empty_str, coor_output_p);
    sfree(coor_output_p);
  }
}


//---------------------------------------------------------------------------
// WspectExtendedMesons::collect(...)
//---------------------------------------------------------------------------
void WspectExtendedMesons::collect(const WspectQuark &q_l, WspectQuark &q_nl, WspectFuzzing *sink_fuzz_p){
  //collect correlators with given local quark propagator and non-local 
  //quark progagator,  
  
  char *fname="collect";
  VRB.Func(d_class_name, fname);
  
  // maximal conceivable number of states
  // 12 conventional mesons M_i (all polarisations)
  // 12 x 3  M_i x P_j
  // 12 x 3  M_i x B_j
  // 12 x 5  M_i x D_j
  // =================
  // 12 x 12 = 144
  // should be modified at program start, the following is for debugging:

  
  // determine the Source OperatorKind, eg. sour_op=DEV1
  // and the corresponding real operator, eg. sour_op_tmp=DEV2 (if prop_dir==0)
  
  DEVOperatorKind sour_op,sour_op_tmp;
  sour_op=q_nl.srcOpKind();
  sour_op_tmp=(DEVOperatorKind)sour_op;

#ifdef  DEBUG_W_EXT_MESON
  printf("colleting extmesons - source operator = %d %d \n", sour_op, sour_op_tmp);
#endif

  // int sink_op will later loop over all possible operators
  // sink_op_tmp will determine the real operator (according to the prop_dir)

  DEVOperatorKind sink_op_tmp;

  // size of non-local quark propagator on a timeslice 
  int sink_prop_size=q_nl.dataSize()/lcl_sites[d_prop_dir];

  Float *prop_sink_p; // pointer to result of sink operator, must allocate memory before usage
  const Float *ql_p=(Float *)q_l.Data(); //local propagator

  
  // in general: calculate the contraction with all possible sinks (operators)
  //             previously there was a selective choice depending on the
  //             source operator which enters this routine

  // in practice: now a table will determine the weighted contribution of 
  //              each particular combination for the states which are to be calculated
  //              If the table returns zero weight for all states of interests, 
  //              this combination will not be calculated

  // strategy:
  // loop over all sink_op
  //   step1: take into account the prop_dir: sink_op --> sink_op_tmp
  //   step2: check whether the particular combination of source and sink operators is needed
  //          test whether the SUM operator is compatible with the definition of the state
  //   step3: allocate memory for result of sink operator (timeslice)
  //   step4: loop over all timeslices 
  //          step5: apply sink_op_tmp
  //          step6: do all traces ("DiracAlgebra")
  //   step7: release memory
  //repeat above for other related sink operators

  //========================================================================

  
  //loop over all sink operators

  for(int sink_op=0; sink_op < DEV_OP_NUM; sink_op++){ 
  
    // step 1:  map sink_op to sink_op_tmp according to the propagation direction
    // not yet implemented !!!
    sink_op_tmp=(DEVOperatorKind)sink_op;
#ifdef  DEBUG_W_EXT_MESON
     printf("src_op=%d sink_op= %d\n", sour_op,sink_op);
#endif
    // step 2: check if (sink_op, sour_op) is relevant 
    // for any state to be calculated according to the weight table

    Float weight_test=0.;
    testCombination(sour_op,sink_op, weight_test);

#ifdef  DEBUG_W_EXT_MESON
     printf("weight = %e \n", weight_test);
#endif

    // if the combination (sour_op,sink_op) is not needed skip the remainder
    // and continue loop over the next sink operator

    if (weight_test == 0.) continue; 


    // for (weight_test != 0) there will be some state for which
    // the following is useful


    // step3: allocate memory to store output of sink operator
    prop_sink_p = (Float *) smalloc(sink_prop_size * sizeof(Float));
    
#ifdef DEBUG_W_EXT_MESON
    //printf("prop_sink_p allocated at %x(size=%x)\n",prop_sink_p,sink_prop_size);
#endif	  
    if (!prop_sink_p) ERR.Pointer(d_class_name, fname, "prop_sink_p");


#ifdef DEBUG_W_EXT_MESON
    printf("Apply sink operator %d\n",sink_op);
#endif
    
    //if sink_op is second derivative operator
    //allocate extra temprary buffer for it
    Float *prop_sink_tmp_p=0;
    if(sink_op_tmp>=DEV1DEV2){
      prop_sink_tmp_p = (Float *) smalloc(sink_prop_size * sizeof(Float));
      if (!prop_sink_tmp_p) ERR.Pointer(d_class_name, fname, "prop_sink_tmp_p");
    }
     
    // step4: loop over the different timeslices
    int lclWalls=lcl_sites[d_prop_dir];
    for(int lclW=0;lclW<lclWalls;lclW++){
      //run fuzzing here if do fuzzing slice by slice to save memory
      //There are a lot of redunacy, have to rerun fuzzing for this slice
      //for the next sink operator!!!
      
      // step5: apply sink_op_tmp
      
      q_nl.doSinkOperator(lclW, sink_op_tmp, prop_sink_p, prop_sink_tmp_p,sink_fuzz_p);

      // step6: do traces (color,momentum)
      // and calculate trace for all different Gamma combinations (where needed)
      // collect data to relevant state

      DiracAlgebra(ql_p,prop_sink_p,lclW, sour_op, sink_op);      
      
    } // endfor (lclW<lclWalls)
	  
    //release the memory
    if(prop_sink_tmp_p){
       VRB.Func(d_class_name, fname);
       VRB.Sfree(d_class_name, fname, empty_str, prop_sink_tmp_p);
       sfree(prop_sink_tmp_p);
    }

    VRB.Func(d_class_name, fname);
    VRB.Sfree(d_class_name, fname, empty_str, prop_sink_p);
    sfree(prop_sink_p);


  } // endfor sink_op
 
  return;
}
  



//===========================================================================

//---------------------------------------------------------------------------
// void WspectExtendedMesons::finish(...)
//---------------------------------------------------------------------------
void WspectExtendedMesons::finish(){
  //do global sum of Complex coor_data_p[glbwall][meson]
#ifdef DEBUG_W_EXT_MESON
  printf("finish().Do global sum.\n");
#endif
  Float *Flt_p = (Float *)coor_data_p;
  for (int i = 0; i < d_size; ++i)    glb_sum(Flt_p++);
  
}

//===========================================================================




//---------------------------------------------------------------------------
// void WspectExtendedMesons::DiracAlgebra(...)
//---------------------------------------------------------------------------
// qp1 is local prop, qp2 is non-local and after doSinkOperator()
void
WspectExtendedMesons::DiracAlgebra(const Float* qp1, const Float* qp2, int lclW, int sour_op, int sink_op){ 
  char *fname="DiracAlgebra";
  VRB.Func(d_class_name, fname);
  
  // zero out mom_proj result  -- more effecient this way
  Float *p = (Float *)d_zero_mom_proj;
  //debug
  for (int  i = 0; i < sizeof(d_zero_mom_proj)/sizeof(Float); ++i)
    *p++ = 0.0;
  
  // For fixed Dirac indices do Mom Projection and Color Trace
#ifdef DEBUG_W_EXT_MESON
  //  printf("MomProj..\n");
#endif
  int D1x, D2x, D1y, D2y;
  for (D1x = 0; D1x < DIRACs; D1x++) {    
    for (D1y = 0; D1y < DIRACs; D1y++) {      
      for (D2x = 0; D2x < DIRACs; D2x++) {	  
	for (D2y = 0; D2y < DIRACs; D2y++) {

	  // MomProject traces over colour and space
	  // returns d_zero_mom_proj[D1x][D1y][D2x][D2y], which is
	  // used in the subsequent traceDirac

	  MomProject(qp1, qp2, D1x, D2x, D1y, D2y, lclW);

	}                   
      }
    }
  } // endfor D1x

#ifdef DEBUG_W_EXT_MESON
  //  printf("MomProj done\n");
#endif

  for(int sour_gamma=0;sour_gamma< SOUR_GAMMAS;sour_gamma++){ 
    for(int sink_gamma=0;sink_gamma< SINK_GAMMAS;sink_gamma++){ 

      // this routine is only called if there is some relevant combination of (sour_op,sink_op)
      // now check also the (sour_gamma,sink_gamma)

      int state = 0;
      int traceDiracDone=0;
      Complex trace;
      int sign=1;
      Float weight = 0.;
      while (state < NUM_WMESON_STATE){
	weight = table(state, sour_gamma, sour_op, sink_gamma, sink_op);
	
	
	if(weight!=0){
	  if(!traceDiracDone){
	    //do traceDirac once(set flag traceDiracDone at end so that
	    //no traceDirac is called for other states

	    // modify the Gamma matrices at source and sink for the trace term
	    // 1. if prop_direction is not time --> not yet implemented
	    // 2. the gamma_5-algebra maps each G[i] into its bit complement
	    //    complement 0010 = G[4]  --> 1101 = G[15-4]
	    // 3. step 2 entail also an additional sign factor 
	    
	    // the binary values for sour_gamma and sink_gamma
	    int sour_n1,sour_n2,sour_n3,sour_n4;
	    int sink_n1,sink_n2,sink_n3,sink_n4;
	    
	    // the complementary gamma index (after gamma5 algebra)
	    int sour_index = 15 - sour_gamma;
	    int sink_index = 15 - sink_gamma;
	    
	    // Gamma(index) = g1^n1*g2^n2*g3^n3*g4^n4;
	    // binary index = n1 + n2*2 + n3*4 + n4*8
	    
	    getBinary(sour_gamma,sour_n1,sour_n2,sour_n3,sour_n4);
	    getBinary(sink_gamma,sink_n1,sink_n2,sink_n3,sink_n4);
	    
	    // the combinatorial signs are determined according to the formulae
	    
	    int sign_sour,sign_sink; // +1 or -1
	    
	    // sign at source: sign_sour = (-1)**(sour_n2+sour_n4); 
	    if ( ((sour_n2+sour_n4)%2) == 0){
	      sign_sour = 1;
	    }else{
	      sign_sour = -1;
	    }
	    
	    // sign at sink:
	    int help= sink_n1*sink_n2 + sink_n1*sink_n3 + sink_n1*sink_n4 + \
	      sink_n2*sink_n3 + sink_n2*sink_n4 + sink_n3*sink_n4 + sink_n2;
	    if ((help%2) == 0){
	      sign_sink = 1;
	    }else{
	      sign_sink = -1;
	    }
	    
#ifdef DEBUG_W_EXT_MESON
	    // printf("sour_gamma, sink_gamma %d, %d \n",sour_gamma, sink_gamma);
	    // printf("complementary index    %d, %d \n",sour_index, sink_index);
	    // printf("combinatorical sign    %d, %d \n",sign_sour, sign_sink);
#endif
	    
	    
	    // overall sign is product of sign at sink and sour
	    sign = sign_sink*sign_sour;
	    
	    // point to source and sink Gamma's
	    Float* sink_gamma_mat = (Float *)WGamma + sink_index*DIRACs*DIRACs*COMPLEXs;
	    Float* sour_gamma_mat = (Float *)WGamma + sour_index*DIRACs*DIRACs*COMPLEXs;
	    
#ifdef DEBUG_W_EXT_MESON
	    // printf("calculate trace: source-Gamma = %d,  sink-Gamma = %d \n",sour_index,sink_index); 
#endif
	
	    traceDirac(sink_gamma_mat, sour_gamma_mat, trace);
	    
	    traceDiracDone=1;//set flag to be true to 

	  }//if(!traceDiracDone)
	  
	  //add to related states
	  int coor_offset=(lclW + lcl2glb_offset[d_prop_dir])*NUM_WMESON_STATE;
	  Complex *mesonsAddr = coor_data_p + coor_offset+state;
	  *(mesonsAddr)  += sign*weight*trace;
	}//if(weight!=0)
	state++;
      }//while(state<..)
      
     traceDiracDone=0;
    } // endfor sink_gamma
  } // endfor sour_gamma

  //must call WspectExtendedMesons::finish() to do global sum before print!
}






// ==================================================================

//---------------------------------------------------------------------------
// void WspectExtendedMesons::ColorAlgebra(...)
//---------------------------------------------------------------------------
// The quark propagator data is stored as
//        Float[Dy][Cy][T][Z][Y][X][Dx][Cx][2]
// The same for all extended mesons
//---------------------------------------------------------------------------

/*#define USE_CBUF
#ifdef USE_CBUF
const unsigned CBUF_MODE2 = 0xcca52112;
const unsigned CBUF_MODE4 = 0xcca52112;
const unsigned MYBANK4_BASE=BANK4_BASE;
const unsigned MYBANK2_BASE=BANK2_BASE;
const unsigned MYBANK_SIZE=BANK_SIZE;
#else
const unsigned MYBANK4_BASE=0;
const unsigned MYBANK2_BASE=0;
const unsigned MYBANK_SIZE=0;
#endif
*/

void
WspectExtendedMesons::ColorAlgebra(const Float* q1_p, const Float* q2_p, 
				   int D1x, int D2x, int D1y, int D2y, 
				   const int lcl[LORENTZs], Complex &result) const
{
  /*
  // set cbuf registers
#ifdef USE_CBUF
  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(4, CBUF_MODE4);
#endif
  */
  int local_site_offset1 = siteOffset(lcl) * SPINORs;  
  int local_site_offset2 = siteOffset(lcl, d_prop_dir) * SPINORs;
  int lclWalls=lcl_sites[d_prop_dir]; //for sink slice propogator
 
  const Complex *q1c_p = (const Complex *)(q1_p +
					   WspectQuark::weightSrcDirac()*D1y + 
					   local_site_offset1 + 
					   (COMPLEXs*COLORs) * D1x);
  const Complex *q2c_p = (const Complex *)(q2_p + 
					   WspectQuark::weightSrcDirac()*D2y/lclWalls + 
					   local_site_offset2 + 
					   (COMPLEXs*COLORs) * D2x);
  
  int off_Cy_in_complex1 = WspectQuark::weightSrcColor() / COMPLEXs;
  int off_Cy_in_complex2 = WspectQuark::weightSrcColor() / (COMPLEXs*lclWalls);

  Complex temp;
  
  for (int Cy = 0; Cy < COLORs; Cy++) {
    //can get better performance by optimizing compDotProduct
    compDotProduct((IFloat *)&temp,((IFloat *)&temp)+1,(const IFloat *)q2c_p, (const IFloat *)q1c_p,6);
    result += temp;

    /*This is very inefficient. 
    result += q1c_p[0] * conj(q2c_p[0]);
    result += q1c_p[1] * conj(q2c_p[1]);
    result += q1c_p[2] * conj(q2c_p[2]);
    */
    q1c_p += off_Cy_in_complex1;    
    q2c_p += off_Cy_in_complex2;    
  }
}

//============================================================================

//---------------------------------------------------------------------------
// void WspectExtendedMesons::MomProject(...)
//---------------------------------------------------------------------------
// result stored in d_zero_mom_proj[D1x][D1y][D2x][D2y] 
//---------------------------------------------------------------------------
void
WspectExtendedMesons::MomProject(const Float* q1_p, const Float* q2_p, int D1x, int D2x, int D1y, int D2y, int lclW)
{
  char *fname="MomProject";
  VRB.Func(d_class_name, fname);
  
  // set the coordinate limit in the prop_dir to avoid looping in that dir.
  //-------------------------------------------------------------------------
  d_lclMin[d_prop_dir] = d_lclMax[d_prop_dir] = lclW; 
  
  // loop over local sites [first for the efficiency in calculation of
  // local_site_offset and non-zero momentum projection] and Dirac indexes.
  //-------------------------------------------------------------------------
  int lcl[LORENTZs];                         
  for (lcl[3] = d_lclMin[3]; lcl[3] <= d_lclMax[3]; lcl[3]++) {
    for (lcl[2] = d_lclMin[2]; lcl[2] <= d_lclMax[2]; lcl[2]++) {	
      for (lcl[1] = d_lclMin[1]; lcl[1] <= d_lclMax[1]; lcl[1]++) {
	for (lcl[0] = d_lclMin[0]; lcl[0] <= d_lclMax[0]; lcl[0]++) {
	  
	  Complex colorSum(0.0, 0.0);		
	  ColorAlgebra(q1_p, q2_p, D1x, D2x, D1y, D2y, lcl, colorSum);
	  
	  // zero momentum projection
	  d_zero_mom_proj[D1x][D1y][D2x][D2y] += colorSum;  
	  
	}       // for (lcl[3] = d_lclMin[3];
      }         // for (lcl[2] = d_lclMin[2]
    }           // for (lcl[1] = d_lclMin[1];
  }             // for (lcl[0] = d_lclMin[0];
}



// ===================================================================


void WspectExtendedMesons::traceDirac(Float* gam1, Float* gam2,Complex &result)
{
  char *fname="traceDirac";
  
  //check for errors
  if(DIRACs!=4) ERR.General(d_class_name,fname,"DIRACs!=4");
  if(gam1==0 || gam2==0) ERR.General(d_class_name,fname,"null ptr");
  
  result = Complex (0.,0.);

  for(int d1=0;d1<DIRACs;d1++){
    for(int d2=0;d2<DIRACs;d2++){
      for(int d3=0;d3<DIRACs;d3++){
	for(int d4=0;d4<DIRACs;d4++){
	  
	  Complex gam1_element(gam1[d1*4*2+d2*2],gam1[d1*4*2+d2*2+1]);
	  Complex gam2_element(gam2[d3*4*2+d4*2],gam2[d3*4*2+d4*2+1]);

	  
	  result+=gam1_element*gam2_element*d_zero_mom_proj[d2][d3][d1][d4];
	
	}//d1
      }//d2
    }//d3
  }//d4 
  //done
}


//---------------------------------------------------------------------------
// void WspectExtendedMesons::print(const CommonArg *common_arg) const
//--------------------------------------------------------------------------- 
void
WspectExtendedMesons::print() const {

#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  char *fname = "print";
  char outputfilename[100];//generated according to WpsectArg.filetail and stateName
  VRB.Func(d_class_name, fname);

  //caluculate average of x,y,z, polarizations
  //store average in coor_data_p's x polarization entry
  {
    for(int wall=0;wall< d_glb_walls;wall++){
      int offset_data_wall=wall*NUM_WMESON_STATE;
      int offset_output_wall=wall*NUM_WMESON_OUTPUT;
      int num_pol; //count number of polarization


      //meason is real name, must use map
      for(int mesonId=0;mesonId<NUM_WMESON_OUTPUT;mesonId++){
        //loop over all states to find all polarizations of the same
        //meson
        num_pol=0;
        int output_offset=offset_output_wall+mesonId;
        for(int mesonState=0;mesonState<NUM_WMESON_STATE;mesonState++){

          if(WMesonStateTable[mesonState].mesonId!=mesonId){
            continue;
          }

          num_pol++;
          //check whether the meson is measured
          if(WMesonStateTable[mesonState].measure==0) {
            //measured[mesonId]=0;
            break;
          }

          //printf("Averaging meson#%d...\n",meson);
          int data_offset=offset_data_wall+map[mesonState];
          coor_output_p[output_offset]+=coor_data_p[data_offset];
        }//mesonState

        //divide
        if(num_pol!=0) coor_output_p[output_offset]/=num_pol;
      }//mesonId
    }//wall
  }


	
  
  // Print out mesons
  //------------------------------------------------------------------
  FILE *fp;
  for (int meson = 0; meson < NUM_WMESON_OUTPUT; ++meson) {
    
    //    if(measured[meson]){
    int measure_flag=0;
    //generate the output filename
    //get statename from WMesonStateTable
    for(int state=0;state<NUM_WMESON_STATE;state++){
      if(WMesonStateTable[state].mesonId==meson){
	if(WMesonStateTable[state].measure) {
	  measure_flag=1;
	  if(!arg_p->fuzzing_on){
	    //    sprintf(outputfilename,"%s.%s",WMesonStateTable[state].stateName,arg_p->filetail);
  sprintf(outputfilename,"%s",WMesonStateTable[state].stateName);
	  }else{
	    //fuzzing is on
	    //   sprintf(outputfilename,"%s.%sF%.1f",WMesonStateTable[state].stateName,arg_p->filetail,arg_p->fuzzing_c[fuzzing_c_index]);
	    sprintf(outputfilename,"%sF%.1f",WMesonStateTable[state].stateName,arg_p->fuzzing_c[fuzzing_c_index]);
	  }
	  break;
	}
      }
    }

    if(measure_flag){
#ifdef DEBUG_W_EXT_MESON
      printf("Writing extended meson#%d data to file %s\n",meson,outputfilename);   
#endif
      
      if ( !(fp = Fopen(outputfilename, "a")) )
	ERR.FileA(d_class_name,fname, outputfilename);

      for (int wall = 0; wall <= d_glb_walls/2; ++wall) {
	
	Complex tmp = coor_output_p[meson + NUM_WMESON_OUTPUT*
				 ((d_whr.glbCoord() + wall)%d_glb_walls)];
	tmp += coor_output_p[meson +  NUM_WMESON_OUTPUT*
			  ((d_whr.glbCoord() + d_glb_walls - wall)%d_glb_walls)];
	tmp /= 2.0; //average of t and N-t
	// only print out the imaginary parts
	//mom=0
	Fprintf(fp, "%d %d %d %.10e\n", 
		AlgWspect::GetCounter(), wall, 0, tmp.real());
      }
      Fclose(fp);
    }//if(measur_flag)
  }//for meson
}

//---------------------------------------------------------------------------
// void WspectExtendedMesons::getBinary(i, &i1, &i2, &i3, &i4) const
// decompose an integer into its binary components i = i1 + i2*2 + i3*4 + i4*8
// is there any easier way to get binary numbers (i1,i2,i3,i4) ?
//--------------------------------------------------------------------------- 
void
WspectExtendedMesons::getBinary(const int i, int &i1, int &i2, int &i3, int &i4) const
{
  i1 = i%2;
  
  i2 = (i - i1)/2;
  i2 = i2%2;
  
  i3 = (i - i1 - 2*i2)/4;
  i3 = i3%2;
  
  i4 = (i - i1 - 2*i2 - 4*i3)/8;
}






//---------------------------------------------------------------------------
// void WspectExtendedMesons::testCombination(sour_op,sink_op, Float &weight_test) const
// tests whether or not to calculate a certain combination (sour_op,sink_op)
//
// an extra check is performed in case the source is one of the sum operators
// in this case the weighted sum must correspond to an state combination
// in table, e.g.  
// assume source SUM_F = (c1*DEV1 + c2*DEV2 + c3*DEV3)
// assume sink   DEV1
// assume state  w1*tr [DEV1 DEV1]  + w2*tr [DEV1 DEV2] + w3*tr [DEV1 DEV3]
// the weights are read from Table and they should be multiples of c_i
// If this is not the case a warning is printed and this combination [DEV1 SUM] 
// is not calculated
// 
// if the sour_op is a SUM (S = \sum c_l O_l) over basic operators
// check that this is compatible with the weights defined in table
// state += w(i,j,k) Tr[ G_i O_j ... G_k \sum_l c_l O_l ... ] 
//       ?= \sum_l w(i,j,k,l) Tr[ G_i O_j ... G_k  O_l ... ] 
// w(i,j,k,l) ?= w(i,j,k)*c_l 
// are all w(i,j,k,l) either zero or same multiples of c_l
// implementation ready only for c_l = 1.
//
// set sour_gamma=-1 and sink_gamma=-1 so this check is done
// independently of the DiracAlgebra 

//--------------------------------------------------------------------------- 
void
WspectExtendedMesons::testCombination(const int sour_op, const int sink_op, Float &weight_test) const
{
  int state = 0;
  // loop over all states as long as weight_test==0

  while (weight_test==0. && state < NUM_WMESON_STATE){
    weight_test = table(state, -1, sour_op, -1, sink_op);
    state += 1; // test next state
  } // endwhile (weight_test != 0)
  
}


Float WspectExtendedMesons::table(int stateId, 
				  int src_gamma,int src_devop, 
				  int sink_gamma, int sink_devop) const{

   if(WMesonStateTable[stateId].measure==0) return 0;

   //ignore the argument if the value is negative
   int sink_mesonOp=WMesonStateTable[stateId].sinkOp;
   int src_mesonOp=WMesonStateTable[stateId].srcOp;
   int cat=WMesonStateTable[stateId].category;
   //use WMesonOpTable to find the details
   int i;
   int sink_weight=0;
   int src_weight=0;
   int src_found=0;
   int sink_found=0;

   

   for(i=0;i<WMesonOpTable[sink_mesonOp].num_terms;i++){
     sink_weight=WMesonOpTable[sink_mesonOp].terms[i][0];
     if(sink_gamma>=0 && WMesonOpTable[sink_mesonOp].terms[i][1]!=sink_gamma) continue ;
     if(sink_devop>=0 && WMesonOpTable[sink_mesonOp].terms[i][2]!=sink_devop) continue ;
     sink_found=1;
     break;
   }


   if(sink_found==0) return 0;
   
   for(i=0;i<WMesonOpTable[src_mesonOp].num_terms;i++){
     src_weight=WMesonOpTable[src_mesonOp].terms[i][0];
     if(src_gamma>=0 && WMesonOpTable[src_mesonOp].terms[i][1]!=src_gamma) continue ;

     //check source operator. Be careful dealing with operators containing SUM_S_SYM
     //and SUM_S_ANTISYM
     
     int term_src_devop=WMesonOpTable[src_mesonOp].terms[i][2];
 
     if(src_devop>=SUM_F){
       //is SUM operator
       switch(src_devop)
	 {
	 case SUM_F:
	   if(matchSUMOp(term_src_devop,src_devop)) src_found=1;
 	   break;
	   
	 case SUM_S_SYM:
	   if(cat!=EXT_SECONDDEV_SYM_MESON) return 0;
	   if(matchSUMOp(term_src_devop,src_devop)) src_found=1;
	   break;
	 
	 case SUM_S_ANTISYM:
	   if(cat!=EXT_SECONDDEV_ANTISYM_MESON) return 0;
	   if(matchSUMOp(term_src_devop,src_devop)) src_found=1;
	   break;
	   
	 case SUM_S_DIAG:
	   if(cat!=EXT_SECONDDEV_DIAG_MESON) return 0;
	   if(matchSUMOp(term_src_devop,src_devop)) src_found=1;
	   break;
	
	 case SUM_F_S_ANTISYM:
	   if(cat!=EXT_FIRSTDEV_MESON &&
	      cat!=EXT_SECONDDEV_ANTISYM_MESON) return 0;
	   if(matchSUMOp(term_src_devop,src_devop)) src_found=1;
	   break;
	   
	 case SUM_S_SYM_DIAG:
	   if(cat!=EXT_SECONDDEV_SYM_MESON &&
	      cat!=EXT_SECONDDEV_DIAG_MESON) return 0;
	   if(matchSUMOp(term_src_devop,src_devop)) src_found=1;
	   break;
	   
	 case SUM_UNIT_F_S_ANTISYM:
	   if(cat!=NORMALMESON &&
	      cat!=EXT_FIRSTDEV_MESON &&
	      cat!=EXT_SECONDDEV_ANTISYM_MESON) return 0;
	   if(matchSUMOp(term_src_devop,src_devop)) src_found=1;
	   break;
	  
	 }//end of switch
       
       if(src_found==1) break;//out of for loop
       continue; //next term
     }else{
       //not sum operator, should exact match
       if(src_devop>=0 && term_src_devop !=src_devop) continue ;
       src_found=1;
       break;
     }
   }//for(src_terms)

   if(src_found*sink_found!=0){
     //found the term
     return (Float)(src_weight*sink_weight);
   }
   return 0;
} 
 
int WspectExtendedMesons::isInOpGroup(int op, int groupId) const{
  char *fname="isInOpGroup";
  int result=0;
  switch(groupId){
  case 0:
    if(op>=0 && op< DEV_OP_NUM) result=1;
    break;
  case 1:
    if(op==UNIT || op==SUM_F || op == SUM_S_ANTISYM || 
       op==SUM_S_SYM || op==SUM_S_DIAG) result=1;
    break;
  case 2:
    if(op==UNIT || op==SUM_F_S_ANTISYM || 
       op==SUM_S_SYM || op==SUM_S_DIAG) result= 1;
    break;
  case 3:
    if(op==UNIT || op==SUM_F_S_ANTISYM || 
       op==SUM_S_SYM_DIAG) result=1;
    break;
  case 4:
    if(op==SUM_UNIT_F_S_ANTISYM || 
       op==SUM_S_SYM_DIAG) result= 1;
    break;
    
  default:
    ERR.General(d_class_name, fname, "invalid grpId");
  }

  return result;
 
}
//check if sum operator sum_op contains op  
int WspectExtendedMesons::matchSUMOp(int op, int sum_op) const{
  
  int match=0;

  if(sum_op>=SUM_F){
    //is SUM operator
    switch(sum_op){
    case SUM_F:
      if(op==DEV1 || op==DEV2 || op==DEV3) match=1;
      break;
    case SUM_S_SYM:
      if(op==DEV1DEV2 || op==DEV1DEV3 || op==DEV2DEV3) match=1;
      break;
      
    case SUM_S_ANTISYM:
      if(op==DEV1DEV2 || op==DEV1DEV3 || op==DEV2DEV3) match=1;
      break;
      
    case SUM_S_DIAG:
      if(op==DEV1DEV1 || op==DEV2DEV2 || op==DEV3DEV3) match =1;
      break;
      
    case SUM_F_S_ANTISYM:
      if(op==DEV1 || op==DEV2 || op==DEV3 || 
	 op==DEV1DEV2 || op==DEV1DEV3 || op==DEV2DEV3) match =1;
      break;
      
    case SUM_S_SYM_DIAG:
      if(op==DEV1DEV2 || op==DEV1DEV3 || op==DEV2DEV3 ||
	 op==DEV1DEV1 || op==DEV2DEV2 || op==DEV3DEV3 ) match=1;
      break;
      
    case SUM_UNIT_F_S_ANTISYM:
      if(op==UNIT || op==DEV1 || op==DEV2 || op==DEV3 || 
	 op==DEV1DEV2 || op==DEV1DEV3 || op==DEV2DEV3) match =1;
      break;
    }//end of switch
  }else{
    //if sum_op is not a sum operator
    if(sum_op==op) match=1;
    
  }
  
  return match;
  
}

void WspectExtendedMesons::setWMesonOpTerm(int *term_p, int weight, WGammaMatrix gammaMat, DEVOperatorKind opKind){
  term_p[0]=weight;
  term_p[1]=gammaMat;
  term_p[2]=opKind;
}

//-------------------------------
// void initWMesonOpTable()
//-------------------------------
void WspectExtendedMesons::initWMesonOpTable(){
  //Normal MesonOperators
  /*
  //MO_a0,           
  WMesonOpTable[MO_a0].num_terms=1;
  setWMesonOpTerm(&(WMesonOpTable[MO_a0].terms[0][0]),1,WUNIT,UNIT);

  //MO_a0_prime,
  WMesonOpTable[MO_a0_prime].num_terms=1;
  setWMesonOpTerm(&(WMesonOpTable[MO_a0_prime].terms[0][0]),1,WGAM_4,UNIT);

  //MO_a1_x,         MO_a1_y,         MO_a1_z,
  WMesonOpTable[MO_a1_x].num_terms=1;
  setWMesonOpTerm(&(WMesonOpTable[MO_a1_x].terms[0][0]),-1,WGAM_1_5,UNIT);
  WMesonOpTable[MO_a1_y]=WMesonOpTable[MO_a1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a1_y].terms[0][0]),1,WGAM_5_2,UNIT);
  WMesonOpTable[MO_a1_z]=WMesonOpTable[MO_a1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a1_z].terms[0][0]),-1,WGAM_3_5,UNIT);

  //MO_b1_x,         MO_b1_y,         MO_b1_z,
  WMesonOpTable[MO_b1_x].num_terms=1;
  setWMesonOpTerm(&(WMesonOpTable[MO_b1_x].terms[0][0]),-1,WGAM_2_3,UNIT);
  WMesonOpTable[MO_b1_y]=WMesonOpTable[MO_b1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_b1_y].terms[0][0]),1,WGAM_1_3,UNIT);
  WMesonOpTable[MO_b1_z]=WMesonOpTable[MO_b1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_b1_z].terms[0][0]),-1,WGAM_1_2,UNIT);
 
  //MO_rho_x,        MO_rho_y,        MO_rho_z,
  WMesonOpTable[MO_rho_x].num_terms=1;
  setWMesonOpTerm(&(WMesonOpTable[MO_rho_x].terms[0][0]),1,WGAM_1,UNIT);
  WMesonOpTable[MO_rho_y]=WMesonOpTable[MO_rho_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_rho_y].terms[0][0]),1,WGAM_2,UNIT);
  WMesonOpTable[MO_rho_z]=WMesonOpTable[MO_rho_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_rho_z].terms[0][0]),1,WGAM_3,UNIT);

  //MO_rho_x_prime,  MO_rho_y_prime,  MO_rho_z_prime,
  WMesonOpTable[MO_rho_prime_x].num_terms=1;
  setWMesonOpTerm(&(WMesonOpTable[MO_rho_prime_x].terms[0][0]),-1,WGAM_1_4,UNIT);
  WMesonOpTable[MO_rho_prime_y]=WMesonOpTable[MO_rho_prime_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_rho_prime_y].terms[0][0]),-1,WGAM_2_4,UNIT);
  WMesonOpTable[MO_rho_prime_z]=WMesonOpTable[MO_rho_prime_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_rho_prime_z].terms[0][0]),-1,WGAM_3_4,UNIT);

  //MO_pion,         
  WMesonOpTable[MO_pion].num_terms=1;
  setWMesonOpTerm(&(WMesonOpTable[MO_pion].terms[0][0]),1,WGAM_5,UNIT);
  
  //MO_pion_prime,
  WMesonOpTable[MO_pion_prime].num_terms=1;
  setWMesonOpTerm(&(WMesonOpTable[MO_pion_prime].terms[0][0]),-1,WGAM_5_4,UNIT);
  */

  //-----------------------
  //ExtendedMeson Operators
  //-----------------------
  //MO_a0xP_x, MO_a0xP_y, MO_a0xP_z,
  WMesonOpTable[MO_a0xP_x].num_terms=1;
  setWMesonOpTerm(&(WMesonOpTable[MO_a0xP_x].terms[0][0]),1,WUNIT,DEV1);

  WMesonOpTable[MO_a0xP_y]=WMesonOpTable[MO_a0xP_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a0xP_y].terms[0][0]),1,WUNIT,DEV2);

  WMesonOpTable[MO_a0xP_z]=WMesonOpTable[MO_a0xP_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a0xP_z].terms[0][0]),1,WUNIT,DEV3);

  //MO_pionxP_x, MO_pionxP_y, MO_pionxP_z,
  WMesonOpTable[MO_pionxP_x].num_terms=1;
  setWMesonOpTerm(&(WMesonOpTable[MO_pionxP_x].terms[0][0]),1,WGAM_5,DEV1);

  WMesonOpTable[MO_pionxP_y]=WMesonOpTable[MO_pionxP_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_pionxP_y].terms[0][0]),1,WGAM_5,DEV2);

  WMesonOpTable[MO_pionxP_z]=WMesonOpTable[MO_pionxP_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_pionxP_z].terms[0][0]),1,WGAM_5,DEV3);

  //MO_a0_primexP_x, MO_a0_primexP_y, MO_a0_primexP_z,
  WMesonOpTable[MO_a0_primexP_x].num_terms=1;
  setWMesonOpTerm(&(WMesonOpTable[MO_a0_primexP_x].terms[0][0]),1,WGAM_4,DEV1);

  WMesonOpTable[MO_a0_primexP_y]=WMesonOpTable[MO_a0_primexP_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a0_primexP_y].terms[0][0]),1,WGAM_4,DEV2);

  WMesonOpTable[MO_a0_primexP_z]=WMesonOpTable[MO_a0_primexP_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a0_primexP_z].terms[0][0]),1,WGAM_4,DEV3);

  // =================================================================
  // rho x P

  //MO_rhoxP_A1
  WMesonOpTable[MO_rhoxP_A1].num_terms=3;
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxP_A1].terms[0][0]),1,WGAM_1,DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxP_A1].terms[1][0]),1,WGAM_2,DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxP_A1].terms[2][0]),1,WGAM_3,DEV3);

  //MO_rhoxP_T1_x, MO_rhoxP_T1_y, MO_rhoxP_T1_z,
  WMesonOpTable[MO_rhoxP_T1_x].num_terms=2;
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxP_T1_x].terms[0][0]), 1,WGAM_2,DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxP_T1_x].terms[1][0]),-1,WGAM_3,DEV2);

  WMesonOpTable[MO_rhoxP_T1_y]=WMesonOpTable[MO_rhoxP_T1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxP_T1_y].terms[0][0]), 1,WGAM_3,DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxP_T1_y].terms[1][0]),-1,WGAM_1,DEV3);

  WMesonOpTable[MO_rhoxP_T1_z]=WMesonOpTable[MO_rhoxP_T1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxP_T1_z].terms[0][0]), 1,WGAM_1,DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxP_T1_z].terms[1][0]),-1,WGAM_2,DEV1);

  //MO_rhoxP_T2_x, MO_rhoxP_T2_y, MO_rhoxP_T2_z,
  WMesonOpTable[MO_rhoxP_T2_x].num_terms=2;
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxP_T2_x].terms[0][0]),1,WGAM_2,DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxP_T2_x].terms[1][0]),1,WGAM_3,DEV2);

  WMesonOpTable[MO_rhoxP_T2_y]=WMesonOpTable[MO_rhoxP_T2_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxP_T2_y].terms[0][0]),1,WGAM_3,DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxP_T2_y].terms[1][0]),1,WGAM_1,DEV3);

  WMesonOpTable[MO_rhoxP_T2_z]=WMesonOpTable[MO_rhoxP_T2_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxP_T2_z].terms[0][0]),1,WGAM_1,DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxP_T2_z].terms[1][0]),1,WGAM_2,DEV1);

  // =============================================================================
  // a1 x P = A1 + T1 + T2 + E

  //MO_a1xP_A1
  WMesonOpTable[MO_a1xP_A1].num_terms=3;
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xP_A1].terms[0][0]),-1,WGAM_1_5,DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xP_A1].terms[1][0]), 1,WGAM_5_2,DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xP_A1].terms[2][0]),-1,WGAM_3_5,DEV3);

  // a1xP (T2)
  //MO_a1xP_T2_x, MO_a1xP_T2_y, MO_a1xP_T2_z,
  WMesonOpTable[MO_a1xP_T2_x].num_terms=2;
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xP_T2_x].terms[0][0]), 1,WGAM_5_2,DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xP_T2_x].terms[1][0]),-1,WGAM_3_5,DEV2);

  WMesonOpTable[MO_a1xP_T2_y]=WMesonOpTable[MO_a1xP_T2_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xP_T2_y].terms[0][0]),-1,WGAM_1_5,DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xP_T2_y].terms[1][0]),-1,WGAM_3_5,DEV1);
  
  WMesonOpTable[MO_a1xP_T2_z]=WMesonOpTable[MO_a1xP_T2_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xP_T2_z].terms[0][0]),-1,WGAM_1_5,DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xP_T2_z].terms[1][0]), 1,WGAM_5_2,DEV1);

  // a1 x P (E) =   Sajk g5gj Dk

  WMesonOpTable[MO_a1xP_E_1].num_terms=2;
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xP_E_1].terms[0][0]),-1,WGAM_1_5,DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xP_E_1].terms[1][0]),-1,WGAM_5_2,DEV2);

  WMesonOpTable[MO_a1xP_E_2].num_terms=2;
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xP_E_2].terms[0][0]), 1,WGAM_5_2,DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xP_E_2].terms[1][0]), 1,WGAM_3_5,DEV3);

  // =============================================================================

  //MO_b1xP_T1_x, MO_b1xP_T1_y, MO_b1xP_T1_z,
  WMesonOpTable[MO_b1xP_T1_x].num_terms=2;
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xP_T1_x].terms[0][0]),1,WGAM_1_3,DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xP_T1_x].terms[1][0]),1,WGAM_1_2,DEV2);

  WMesonOpTable[MO_b1xP_T1_y]=WMesonOpTable[MO_b1xP_T1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xP_T1_y].terms[0][0]),-1,WGAM_1_2,DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xP_T1_y].terms[1][0]), 1,WGAM_2_3,DEV3);

  WMesonOpTable[MO_b1xP_T1_z]=WMesonOpTable[MO_b1xP_T1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xP_T1_z].terms[0][0]),-1,WGAM_2_3,DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xP_T1_z].terms[1][0]),-1,WGAM_1_3,DEV1);


  // =============================================================================
  // a0_prime x D:  T2 = 2+-

  //MO_a0_primexD_x, MO_a0_primexD_y, MO_a0_primexD_z,
  WMesonOpTable[MO_a0_primexD_x].num_terms=2;
  setWMesonOpTerm(&(WMesonOpTable[MO_a0_primexD_x].terms[0][0]),1,WGAM_4,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a0_primexD_x].terms[1][0]),1,WGAM_4,DEV3DEV2);

  WMesonOpTable[MO_a0_primexD_y]=WMesonOpTable[MO_a0_primexD_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a0_primexD_y].terms[0][0]),1,WGAM_4,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a0_primexD_y].terms[1][0]),1,WGAM_4,DEV1DEV3);

  WMesonOpTable[MO_a0_primexD_z]=WMesonOpTable[MO_a0_primexD_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a0_primexD_z].terms[0][0]),1,WGAM_4,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a0_primexD_z].terms[1][0]),1,WGAM_4,DEV2DEV1);
  
  // =============================================================================
  // rho x B  = 

  // rho x B (T1)   1-+ 
  // MO_rhoxB_T1_x, MO_rhoxB_T1_y, MO_rhoxB_T1_z,
  WMesonOpTable[MO_rhoxB_T1_x].num_terms=4;
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T1_x].terms[0][0]), 1,WGAM_2,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T1_x].terms[1][0]),-1,WGAM_2,DEV2DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T1_x].terms[2][0]),-1,WGAM_3,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T1_x].terms[3][0]), 1,WGAM_3,DEV1DEV3);

  WMesonOpTable[MO_rhoxB_T1_y]=WMesonOpTable[MO_rhoxB_T1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T1_y].terms[0][0]), 1,WGAM_3,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T1_y].terms[1][0]),-1,WGAM_3,DEV3DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T1_y].terms[2][0]),-1,WGAM_1,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T1_y].terms[3][0]), 1,WGAM_1,DEV2DEV1);
  
  WMesonOpTable[MO_rhoxB_T1_z]=WMesonOpTable[MO_rhoxB_T1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T1_z].terms[0][0]), 1,WGAM_1,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T1_z].terms[1][0]),-1,WGAM_1,DEV1DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T1_z].terms[2][0]),-1,WGAM_2,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T1_z].terms[3][0]), 1,WGAM_2,DEV3DEV2);

  // rho x B (T2)   2-+ 
  // MO_rhoxB_T2_x, MO_rhoxB_T2_y, MO_rhoxB_T2_z,
  WMesonOpTable[MO_rhoxB_T2_x].num_terms=4;
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T2_x].terms[0][0]), 1,WGAM_2,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T2_x].terms[1][0]),-1,WGAM_2,DEV2DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T2_x].terms[2][0]), 1,WGAM_3,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T2_x].terms[3][0]),-1,WGAM_3,DEV1DEV3);

  WMesonOpTable[MO_rhoxB_T2_y]=WMesonOpTable[MO_rhoxB_T2_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T2_y].terms[0][0]), 1,WGAM_3,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T2_y].terms[1][0]),-1,WGAM_3,DEV3DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T2_y].terms[2][0]), 1,WGAM_1,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T2_y].terms[3][0]),-1,WGAM_1,DEV2DEV1);
  
  WMesonOpTable[MO_rhoxB_T2_z]=WMesonOpTable[MO_rhoxB_T2_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T2_z].terms[0][0]), 1,WGAM_1,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T2_z].terms[1][0]),-1,WGAM_1,DEV1DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T2_z].terms[2][0]), 1,WGAM_2,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxB_T2_z].terms[3][0]),-1,WGAM_2,DEV3DEV2);


  //=============================================================================
  // a1 x B
  
  // a1 x B (A1) =   gig5 Bi

  WMesonOpTable[MO_a1xB_A1].num_terms=6;
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_A1].terms[0][0]), 1,WGAM_1_5,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_A1].terms[1][0]),-1,WGAM_1_5,DEV3DEV2);

  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_A1].terms[2][0]),-1,WGAM_5_2,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_A1].terms[3][0]), 1,WGAM_5_2,DEV1DEV3);

  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_A1].terms[4][0]), 1,WGAM_3_5,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_A1].terms[5][0]),-1,WGAM_3_5,DEV2DEV1);
 
  // a1xB_T1 = MO_a1xB_T1_x, MO_a1xB_T1_y, MO_a1xB_T1_z,
  WMesonOpTable[MO_a1xB_T1_x].num_terms=4;
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T1_x].terms[0][0]),1,WGAM_5_2,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T1_x].terms[1][0]),-1,WGAM_5_2,DEV2DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T1_x].terms[2][0]),1,WGAM_3_5,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T1_x].terms[3][0]),-1,WGAM_3_5,DEV1DEV3);

  WMesonOpTable[MO_a1xB_T1_y]=WMesonOpTable[MO_a1xB_T1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T1_y].terms[0][0]),-1,WGAM_3_5,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T1_y].terms[1][0]),1,WGAM_3_5,DEV3DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T1_y].terms[2][0]),1,WGAM_1_5,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T1_y].terms[3][0]),-1,WGAM_1_5,DEV2DEV1);
  
  WMesonOpTable[MO_a1xB_T1_z]=WMesonOpTable[MO_a1xB_T1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T1_z].terms[0][0]),-1,WGAM_1_5,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T1_z].terms[1][0]),1,WGAM_1_5,DEV1DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T1_z].terms[2][0]),-1,WGAM_5_2,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T1_z].terms[3][0]),1,WGAM_5_2,DEV3DEV2);

  //MO_a1xB_T2_x, MO_a1xB_T2_y, MO_a1xB_T2_z, 
  WMesonOpTable[MO_a1xB_T2_x].num_terms=4;
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T2_x].terms[0][0]), 1,WGAM_5_2,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T2_x].terms[1][0]),-1,WGAM_5_2,DEV2DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T2_x].terms[2][0]),-1,WGAM_3_5,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T2_x].terms[3][0]), 1,WGAM_3_5,DEV1DEV3);

  WMesonOpTable[MO_a1xB_T2_y]=WMesonOpTable[MO_a1xB_T2_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T2_y].terms[0][0]),-1,WGAM_3_5,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T2_y].terms[1][0]),1,WGAM_3_5,DEV3DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T2_y].terms[2][0]),-1,WGAM_1_5,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T2_y].terms[3][0]),1,WGAM_1_5,DEV2DEV1);
  
  WMesonOpTable[MO_a1xB_T2_z]=WMesonOpTable[MO_a1xB_T2_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T2_z].terms[0][0]),-1,WGAM_1_5,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T2_z].terms[1][0]),1,WGAM_1_5,DEV1DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T2_z].terms[2][0]),1,WGAM_5_2,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xB_T2_z].terms[3][0]),-1,WGAM_5_2,DEV3DEV2);

  // ==========================================================================
  // a1 x D

  // a1 x D (A2) =  g5 g_i s_ijk D_j D_k
  WMesonOpTable[MO_a1xD_A2].num_terms=6;
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_A2].terms[0][0]),-1,WGAM_1_5,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_A2].terms[1][0]),-1,WGAM_1_5,DEV3DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_A2].terms[2][0]), 1,WGAM_5_2,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_A2].terms[3][0]), 1,WGAM_5_2,DEV1DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_A2].terms[4][0]),-1,WGAM_3_5,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_A2].terms[5][0]),-1,WGAM_3_5,DEV2DEV1);

  //MO_a1xD_T1_x, MO_a1xD_T1_y, MO_a1xD_T1_z = s_{ijk} g5gj s_{klm} Dl Dm
  WMesonOpTable[MO_a1xD_T1_x].num_terms=4;
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T1_x].terms[0][0]), 1,WGAM_5_2,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T1_x].terms[1][0]), 1,WGAM_5_2,DEV2DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T1_x].terms[2][0]),-1,WGAM_3_5,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T1_x].terms[3][0]),-1,WGAM_3_5,DEV1DEV3);

  WMesonOpTable[MO_a1xD_T1_y]=WMesonOpTable[MO_a1xD_T1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T1_y].terms[0][0]),-1,WGAM_3_5,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T1_y].terms[1][0]),-1,WGAM_3_5,DEV3DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T1_y].terms[2][0]),-1,WGAM_1_5,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T1_y].terms[3][0]),-1,WGAM_1_5,DEV2DEV1);

  WMesonOpTable[MO_a1xD_T1_z]=WMesonOpTable[MO_a1xD_T1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T1_z].terms[0][0]),-1,WGAM_1_5,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T1_z].terms[1][0]),-1,WGAM_1_5,DEV1DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T1_z].terms[2][0]), 1,WGAM_5_2,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T1_z].terms[3][0]), 1,WGAM_5_2,DEV3DEV2);


  //MO_a1xD_T2_x, MO_a1xD_T2_y, MO_a1xD_T2_z = e_{ijk} g5gj s_{klm} Dl Dm
  WMesonOpTable[MO_a1xD_T2_x].num_terms=4;
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T2_x].terms[0][0]), 1,WGAM_5_2,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T2_x].terms[1][0]), 1,WGAM_5_2,DEV2DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T2_x].terms[2][0]),1,WGAM_3_5,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T2_x].terms[3][0]),1,WGAM_3_5,DEV1DEV3);

  WMesonOpTable[MO_a1xD_T2_y]=WMesonOpTable[MO_a1xD_T2_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T2_y].terms[0][0]),-1,WGAM_3_5,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T2_y].terms[1][0]),-1,WGAM_3_5,DEV3DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T2_y].terms[2][0]),1,WGAM_1_5,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T2_y].terms[3][0]),1,WGAM_1_5,DEV2DEV1);

  WMesonOpTable[MO_a1xD_T2_z]=WMesonOpTable[MO_a1xD_T2_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T2_z].terms[0][0]),-1,WGAM_1_5,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T2_z].terms[1][0]),-1,WGAM_1_5,DEV1DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T2_z].terms[2][0]), -1,WGAM_5_2,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_T2_z].terms[3][0]), -1,WGAM_5_2,DEV3DEV2);

  //MO_a1xD_E_1, MO_a1xD_E_2;  S_{ajk} gjg5 s_{klm} Dl Dm
  WMesonOpTable[MO_a1xD_E_1].num_terms=4;
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_E_1].terms[0][0]), 1,WGAM_1_5,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_E_1].terms[1][0]), 1,WGAM_1_5,DEV3DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_E_1].terms[2][0]), 1,WGAM_5_2,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_E_1].terms[3][0]), 1,WGAM_5_2,DEV1DEV3);

  WMesonOpTable[MO_a1xD_E_2]=WMesonOpTable[MO_a1xD_E_1];
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_E_2].terms[0][0]),-1,WGAM_5_2,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_E_2].terms[1][0]),-1,WGAM_5_2,DEV1DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_E_2].terms[2][0]),-1,WGAM_3_5,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_a1xD_E_2].terms[3][0]),-1,WGAM_3_5,DEV2DEV1);
  

  // b1 x D

  // b1 x D (A2) =  g4 g5 g_i s_ijk D_j D_k

  WMesonOpTable[MO_b1xD_A2].num_terms=6;
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_A2].terms[0][0]),-1,WGAM_2_3,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_A2].terms[1][0]),-1,WGAM_2_3,DEV3DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_A2].terms[2][0]), 1,WGAM_1_3,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_A2].terms[3][0]), 1,WGAM_1_3,DEV1DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_A2].terms[4][0]),-1,WGAM_1_2,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_A2].terms[5][0]),-1,WGAM_1_2,DEV2DEV1);

  //MO_b1xD_T1_x, MO_b1xD_T1_y, MO_b1xD_T1_z = s_{ijk} g4g5gj s_{klm} Dl Dm
  WMesonOpTable[MO_b1xD_T1_x].num_terms=4;
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T1_x].terms[0][0]), 1,WGAM_1_3,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T1_x].terms[1][0]), 1,WGAM_1_3,DEV2DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T1_x].terms[2][0]), -1,WGAM_1_2,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T1_x].terms[3][0]), -1,WGAM_1_2,DEV1DEV3);

  WMesonOpTable[MO_b1xD_T1_y]=WMesonOpTable[MO_b1xD_T1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T1_y].terms[0][0]),-1,WGAM_1_2,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T1_y].terms[1][0]),-1,WGAM_1_2,DEV3DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T1_y].terms[2][0]), -1,WGAM_2_3,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T1_y].terms[3][0]), -1,WGAM_2_3,DEV2DEV1);
  
  WMesonOpTable[MO_b1xD_T1_z]=WMesonOpTable[MO_b1xD_T1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T1_z].terms[0][0]),-1,WGAM_2_3,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T1_z].terms[1][0]),-1,WGAM_2_3,DEV1DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T1_z].terms[2][0]), 1,WGAM_1_3,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T1_z].terms[3][0]), 1,WGAM_1_3,DEV3DEV2);

  //MO_b1xD_T2_x, MO_b1xD_T2_y, MO_b1xD_T2_z = e_{ijk} g4g5gj s_{klm} Dl Dm
   WMesonOpTable[MO_b1xD_T2_x].num_terms=4;
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T2_x].terms[0][0]), 1,WGAM_1_3,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T2_x].terms[1][0]), 1,WGAM_1_3,DEV2DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T2_x].terms[2][0]), 1,WGAM_1_2,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T2_x].terms[3][0]), 1,WGAM_1_2,DEV1DEV3);

  WMesonOpTable[MO_b1xD_T2_y]=WMesonOpTable[MO_b1xD_T2_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T2_y].terms[0][0]),-1,WGAM_1_2,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T2_y].terms[1][0]),-1,WGAM_1_2,DEV3DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T2_y].terms[2][0]), 1,WGAM_2_3,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T2_y].terms[3][0]), 1,WGAM_2_3,DEV2DEV1);
  
  WMesonOpTable[MO_b1xD_T2_z]=WMesonOpTable[MO_b1xD_T2_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T2_z].terms[0][0]),-1,WGAM_2_3,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T2_z].terms[1][0]),-1,WGAM_2_3,DEV1DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T2_z].terms[2][0]), -1,WGAM_1_3,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_T2_z].terms[3][0]), -1,WGAM_1_3,DEV3DEV2);

  //MO_b1xD_E_1, MO_b1xD_E_2;  S_{ajk} g4g5gj s_{klm} Dl Dm
  WMesonOpTable[MO_b1xD_E_1].num_terms=4;
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_E_1].terms[0][0]), -1,WGAM_2_3,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_E_1].terms[1][0]), -1,WGAM_2_3,DEV3DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_E_1].terms[2][0]), -1,WGAM_1_3,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_E_1].terms[3][0]), -1,WGAM_1_3,DEV1DEV3);

  WMesonOpTable[MO_b1xD_E_2]=WMesonOpTable[MO_b1xD_E_1];
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_E_2].terms[0][0]),1,WGAM_1_3,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_E_2].terms[1][0]),1,WGAM_1_3,DEV1DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_E_2].terms[2][0]),1,WGAM_1_2,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_b1xD_E_2].terms[3][0]),1,WGAM_1_2,DEV2DEV1);


  // =============================================================================
  // rho x D

  // rho x D (A2) =  g_i s_ijk D_j D_k

  WMesonOpTable[MO_rhoxD_A2].num_terms=6;
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_A2].terms[0][0]), 1,WGAM_1, DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_A2].terms[1][0]), 1,WGAM_1, DEV3DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_A2].terms[2][0]), 1,WGAM_2,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_A2].terms[3][0]), 1,WGAM_2,DEV1DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_A2].terms[4][0]), 1,WGAM_3,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_A2].terms[5][0]), 1,WGAM_3,DEV2DEV1);


  // rho x D (T1) = s_{ijk}  g_j  s_{klm}D_lD_m
  WMesonOpTable[MO_rhoxD_T1_x].num_terms=4;
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T1_x].terms[0][0]), 1,WGAM_2,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T1_x].terms[1][0]), 1,WGAM_2,DEV2DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T1_x].terms[2][0]), 1,WGAM_3,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T1_x].terms[3][0]), 1,WGAM_3,DEV1DEV3);

  WMesonOpTable[MO_rhoxD_T1_y]=WMesonOpTable[MO_rhoxD_T1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T1_y].terms[0][0]), 1,WGAM_3,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T1_y].terms[1][0]), 1,WGAM_3,DEV3DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T1_y].terms[2][0]), 1,WGAM_1,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T1_y].terms[3][0]), 1,WGAM_1,DEV2DEV1);

  WMesonOpTable[MO_rhoxD_T1_z]=WMesonOpTable[MO_rhoxD_T1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T1_z].terms[0][0]), 1,WGAM_1,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T1_z].terms[1][0]), 1,WGAM_1,DEV1DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T1_z].terms[2][0]), 1,WGAM_2,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T1_z].terms[3][0]), 1,WGAM_2,DEV3DEV2);


  // rho x D (T2) =   e_{ijk}  g_j  s_{klm}D_lD_m
  WMesonOpTable[MO_rhoxD_T2_x].num_terms=4;
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T2_x].terms[0][0]), 1,WGAM_2,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T2_x].terms[1][0]), 1,WGAM_2,DEV2DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T2_x].terms[2][0]),-1,WGAM_3,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T2_x].terms[3][0]),-1,WGAM_3,DEV1DEV3);

  WMesonOpTable[MO_rhoxD_T2_y]=WMesonOpTable[MO_rhoxD_T2_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T2_y].terms[0][0]), 1,WGAM_3,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T2_y].terms[1][0]), 1,WGAM_3,DEV3DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T2_y].terms[2][0]),-1,WGAM_1,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T2_y].terms[3][0]),-1,WGAM_1,DEV2DEV1);

  WMesonOpTable[MO_rhoxD_T2_z]=WMesonOpTable[MO_rhoxD_T2_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T2_z].terms[0][0]), 1,WGAM_1,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T2_z].terms[1][0]), 1,WGAM_1,DEV1DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T2_z].terms[2][0]),-1,WGAM_2,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_rhoxD_T2_z].terms[3][0]),-1,WGAM_2,DEV3DEV2);


  

  //added pionxB
  WMesonOpTable[MO_pionxB_T1_x].num_terms=2;
  setWMesonOpTerm(&(WMesonOpTable[MO_pionxB_T1_x].terms[0][0]), 1,WGAM_5,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_pionxB_T1_x].terms[1][0]), -1,WGAM_5,DEV3DEV2);

  WMesonOpTable[MO_pionxB_T1_y]=WMesonOpTable[MO_pionxB_T1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_pionxB_T1_y].terms[0][0]), 1,WGAM_5,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_pionxB_T1_y].terms[1][0]), -1,WGAM_5,DEV1DEV3);

  WMesonOpTable[MO_pionxB_T1_z]=WMesonOpTable[MO_pionxB_T1_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_pionxB_T1_z].terms[0][0]), 1,WGAM_5,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_pionxB_T1_z].terms[1][0]), -1,WGAM_5,DEV2DEV1);
  //pionxD
   WMesonOpTable[MO_pionxD_T2_x].num_terms=2;
  setWMesonOpTerm(&(WMesonOpTable[MO_pionxD_T2_x].terms[0][0]), 1,WGAM_5,DEV2DEV3);
  setWMesonOpTerm(&(WMesonOpTable[MO_pionxD_T2_x].terms[1][0]), 1,WGAM_5,DEV3DEV2);

  WMesonOpTable[MO_pionxD_T2_y]=WMesonOpTable[MO_pionxD_T2_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_pionxD_T2_y].terms[0][0]), 1,WGAM_5,DEV3DEV1);
  setWMesonOpTerm(&(WMesonOpTable[MO_pionxD_T2_y].terms[1][0]), 1,WGAM_5,DEV1DEV3);

  WMesonOpTable[MO_pionxD_T2_z]=WMesonOpTable[MO_pionxD_T2_x];
  setWMesonOpTerm(&(WMesonOpTable[MO_pionxD_T2_z].terms[0][0]), 1,WGAM_5,DEV1DEV2);
  setWMesonOpTerm(&(WMesonOpTable[MO_pionxD_T2_z].terms[1][0]), 1,WGAM_5,DEV2DEV1);

}

//------------------------------
//initWMesonStateTable()
//------------------------------
void WspectExtendedMesons::initWMesonStateTable(WspectArg *arg){
  //---------------
  //Normal Mesons
  //--------------
  /* MS_a0,           MS_a0_prime,
   * MS_a1_x,         MS_a1_y,         MS_a1_z,
   * MS_b1_x,         MS_b1_y,         MS_b1_z,
   * MS_rho_x,        MS_rho_y,        MS_rho_z,
   * MS_rho_x_prime,  MS_rho_y_prime,  MS_rho_z_prime,
   * MS_pion,         MS_pion_prime,
   */
  /*
  // MS_a0,   
  WMesonStateTable[MS_a0].stateName="a0";
  WMesonStateTable[MS_a0].mesonId=a0;
  WMesonStateTable[MS_a0].category=NORMALMESON;
  WMesonStateTable[MS_a0].polarization=0;
  WMesonStateTable[MS_a0].measure=0;
  WMesonStateTable[MS_a0].sinkOp=MO_a0;
  WMesonStateTable[MS_a0].srcOp=MO_a0;

  //        MS_a0_prime,
  WMesonStateTable[MS_a0_prime].stateName="a0_prime";
  WMesonStateTable[MS_a0_prime].mesonId=a0_prime;
  WMesonStateTable[MS_a0_prime].category=NORMALMESON;
  WMesonStateTable[MS_a0_prime].polarization=0;
  WMesonStateTable[MS_a0_prime].measure=0;
  WMesonStateTable[MS_a0_prime].sinkOp=MO_a0_prime;
  WMesonStateTable[MS_a0_prime].srcOp=MO_a0_prime;

  // MS_a1_x,         MS_a1_y,         MS_a1_z,
  WMesonStateTable[MS_a1_x].stateName="a1";
  WMesonStateTable[MS_a1_x].mesonId=a1;
  WMesonStateTable[MS_a1_x].category=NORMALMESON;
  WMesonStateTable[MS_a1_x].polarization=0;
  WMesonStateTable[MS_a1_x].measure=0;
  WMesonStateTable[MS_a1_x].sinkOp=MO_a1_x;
  WMesonStateTable[MS_a1_x].srcOp=MO_a1_x;
  
  WMesonStateTable[MS_a1_y]=WMesonStateTable[MS_a1_x];
  WMesonStateTable[MS_a1_y].polarization=1;
  WMesonStateTable[MS_a1_y].sinkOp=MO_a1_y;
  WMesonStateTable[MS_a1_y].srcOp=MO_a1_y;
  
  WMesonStateTable[MS_a1_z]=WMesonStateTable[MS_a1_x];
  WMesonStateTable[MS_a1_z].polarization=2;
  WMesonStateTable[MS_a1_z].sinkOp=MO_a1_z;
  WMesonStateTable[MS_a1_z].srcOp=MO_a1_z;

  // MS_b1_x,         MS_b1_y,         MS_b1_z,
  WMesonStateTable[MS_b1_x].stateName="b1";
  WMesonStateTable[MS_b1_x].mesonId=b1;
  WMesonStateTable[MS_b1_x].category=NORMALMESON;
  WMesonStateTable[MS_b1_x].polarization=0;
  WMesonStateTable[MS_b1_x].measure=0;
  WMesonStateTable[MS_b1_x].sinkOp=MO_b1_x;
  WMesonStateTable[MS_b1_x].srcOp=MO_b1_x;
  
  WMesonStateTable[MS_b1_y]=WMesonStateTable[MS_b1_x];
  WMesonStateTable[MS_b1_y].polarization=1;
  WMesonStateTable[MS_b1_y].sinkOp=MO_b1_y;
  WMesonStateTable[MS_b1_y].srcOp=MO_b1_y;
  
  WMesonStateTable[MS_b1_z]=WMesonStateTable[MS_b1_x];
  WMesonStateTable[MS_b1_z].polarization=2;
  WMesonStateTable[MS_b1_z].sinkOp=MO_b1_z;
  WMesonStateTable[MS_b1_z].srcOp=MO_b1_z;

  // MS_rho_x,        MS_rho_y,        MS_rho_z,
  WMesonStateTable[MS_rho_x].stateName="rho";
  WMesonStateTable[MS_rho_x].mesonId=rho;
  WMesonStateTable[MS_rho_x].category=NORMALMESON;
  WMesonStateTable[MS_rho_x].polarization=0;
  WMesonStateTable[MS_rho_x].measure=0;
  WMesonStateTable[MS_rho_x].sinkOp=MO_rho_x;
  WMesonStateTable[MS_rho_x].srcOp=MO_rho_x;
  
  WMesonStateTable[MS_rho_y]=WMesonStateTable[MS_rho_x];
  WMesonStateTable[MS_rho_y].polarization=1;
  WMesonStateTable[MS_rho_y].sinkOp=MO_rho_y;
  WMesonStateTable[MS_rho_y].srcOp=MO_rho_y;
  
  WMesonStateTable[MS_rho_z]=WMesonStateTable[MS_rho_x];
  WMesonStateTable[MS_rho_z].polarization=2;
  WMesonStateTable[MS_rho_z].sinkOp=MO_rho_z;
  WMesonStateTable[MS_rho_z].srcOp=MO_rho_z;
  // MS_rho_x_prime,  MS_rho_y_prime,  MS_rho_z_prime,
  WMesonStateTable[MS_rho_prime_x].stateName="rho_prime";
  WMesonStateTable[MS_rho_prime_x].mesonId=rho_prime;
  WMesonStateTable[MS_rho_prime_x].category=NORMALMESON;
  WMesonStateTable[MS_rho_prime_x].polarization=0;
  WMesonStateTable[MS_rho_prime_x].measure=0;
  WMesonStateTable[MS_rho_prime_x].sinkOp=MO_rho_prime_x;
  WMesonStateTable[MS_rho_prime_x].srcOp=MO_rho_prime_x;
  
  WMesonStateTable[MS_rho_prime_y]=WMesonStateTable[MS_rho_prime_x];
  WMesonStateTable[MS_rho_prime_y].polarization=1;
  WMesonStateTable[MS_rho_prime_y].sinkOp=MO_rho_prime_y;
  WMesonStateTable[MS_rho_prime_y].srcOp=MO_rho_prime_y;
  
  WMesonStateTable[MS_rho_prime_z]=WMesonStateTable[MS_rho_prime_x];
  WMesonStateTable[MS_rho_prime_z].polarization=2;
  WMesonStateTable[MS_rho_prime_z].sinkOp=MO_rho_prime_z;
  WMesonStateTable[MS_rho_prime_z].srcOp=MO_rho_prime_z;

  // MS_pion,         
  WMesonStateTable[MS_pion].stateName="pion";
  WMesonStateTable[MS_pion].mesonId=pion;
  WMesonStateTable[MS_pion].category=NORMALMESON;
  WMesonStateTable[MS_pion].polarization=0;
  WMesonStateTable[MS_pion].measure=0;
  WMesonStateTable[MS_pion].sinkOp=MO_pion;
  WMesonStateTable[MS_pion].srcOp=MO_pion;

  //MS_pion_prime,
  WMesonStateTable[MS_pion_prime].stateName="pion_prime";
  WMesonStateTable[MS_pion_prime].mesonId=pion_prime;
  WMesonStateTable[MS_pion_prime].category=NORMALMESON;
  WMesonStateTable[MS_pion_prime].polarization=0;
  WMesonStateTable[MS_pion_prime].measure=0;
  WMesonStateTable[MS_pion_prime].sinkOp=MO_pion_prime;
  WMesonStateTable[MS_pion_prime].srcOp=MO_pion_prime;
  */

  //---------------
  //Extended Mesons
  //--------------
  //MS_a0xP_x, MS_a0xP_y, MS_a0xP_z
  WMesonStateTable[MS_a0xP_x].stateName="a0xP";
  WMesonStateTable[MS_a0xP_x].mesonId=a0xP;
  WMesonStateTable[MS_a0xP_x].category=EXT_FIRSTDEV_MESON;
  WMesonStateTable[MS_a0xP_x].polarization=0;
  WMesonStateTable[MS_a0xP_x].measure=0;
  WMesonStateTable[MS_a0xP_x].sinkOp=MO_a0xP_x;
  WMesonStateTable[MS_a0xP_x].srcOp=MO_a0xP_x;

  WMesonStateTable[MS_a0xP_y]=WMesonStateTable[MS_a0xP_x];
  WMesonStateTable[MS_a0xP_y].polarization=1;
  WMesonStateTable[MS_a0xP_y].sinkOp=MO_a0xP_y;
  WMesonStateTable[MS_a0xP_y].srcOp=MO_a0xP_y;

  WMesonStateTable[MS_a0xP_z]=WMesonStateTable[MS_a0xP_x];
  WMesonStateTable[MS_a0xP_z].polarization=2;
  WMesonStateTable[MS_a0xP_z].sinkOp=MO_a0xP_z;
  WMesonStateTable[MS_a0xP_z].srcOp=MO_a0xP_z;

  //MS_pionxP_x, MS_pionxP_y, MS_pionxP_z,
  WMesonStateTable[MS_pionxP_x].stateName="pionxP";
  WMesonStateTable[MS_pionxP_x].mesonId=pionxP;
  WMesonStateTable[MS_pionxP_x].category=EXT_FIRSTDEV_MESON;
  WMesonStateTable[MS_pionxP_x].polarization=0;
  WMesonStateTable[MS_pionxP_x].measure=0;
  WMesonStateTable[MS_pionxP_x].sinkOp=MO_pionxP_x;
  WMesonStateTable[MS_pionxP_x].srcOp=MO_pionxP_x;
  
  WMesonStateTable[MS_pionxP_y]=WMesonStateTable[MS_pionxP_x];
  WMesonStateTable[MS_pionxP_y].polarization=1;
  WMesonStateTable[MS_pionxP_y].sinkOp=MO_pionxP_y;
  WMesonStateTable[MS_pionxP_y].srcOp=MO_pionxP_y;
  
  WMesonStateTable[MS_pionxP_z]=WMesonStateTable[MS_pionxP_x];
  WMesonStateTable[MS_pionxP_z].polarization=2;
  WMesonStateTable[MS_pionxP_z].sinkOp=MO_pionxP_z;
  WMesonStateTable[MS_pionxP_z].srcOp=MO_pionxP_z;
  
  
  //MS_a0_primexP_x, MS_a0_primexP_y, MS_a0_primexP_z,
  WMesonStateTable[MS_a0_primexP_x].stateName="a0_primexP";
  WMesonStateTable[MS_a0_primexP_x].mesonId=a0_primexP;
  WMesonStateTable[MS_a0_primexP_x].category=EXT_FIRSTDEV_MESON;
  WMesonStateTable[MS_a0_primexP_x].polarization=0;
  WMesonStateTable[MS_a0_primexP_x].measure=0;
  WMesonStateTable[MS_a0_primexP_x].sinkOp=MO_a0_primexP_x;
  WMesonStateTable[MS_a0_primexP_x].srcOp=MO_a0_primexP_x;
  
  WMesonStateTable[MS_a0_primexP_y]=WMesonStateTable[MS_a0_primexP_x];
  WMesonStateTable[MS_a0_primexP_y].polarization=1;
  WMesonStateTable[MS_a0_primexP_y].sinkOp=MO_a0_primexP_y;
  WMesonStateTable[MS_a0_primexP_y].srcOp=MO_a0_primexP_y;
  
  WMesonStateTable[MS_a0_primexP_z]=WMesonStateTable[MS_a0_primexP_x];
  WMesonStateTable[MS_a0_primexP_z].polarization=2;
  WMesonStateTable[MS_a0_primexP_z].sinkOp=MO_a0_primexP_z;
  WMesonStateTable[MS_a0_primexP_z].srcOp=MO_a0_primexP_z;

  //MS_rhoxP_A1
  WMesonStateTable[MS_rhoxP_A1_1].stateName="rhoxP_A1";
  WMesonStateTable[MS_rhoxP_A1_1].mesonId=rhoxP_A1;
  WMesonStateTable[MS_rhoxP_A1_1].category=EXT_FIRSTDEV_MESON;
  WMesonStateTable[MS_rhoxP_A1_1].polarization=0;
  WMesonStateTable[MS_rhoxP_A1_1].measure=0;
  WMesonStateTable[MS_rhoxP_A1_1].sinkOp=MO_rhoxP_A1;
  WMesonStateTable[MS_rhoxP_A1_1].srcOp=MO_rhoxP_A1;

  //MS_rhoxP_T1_x, MS_rhoxP_T1_y, MS_rhoxP_T1_z,
  WMesonStateTable[MS_rhoxP_T1_x].stateName="rhoxP_T1";
  WMesonStateTable[MS_rhoxP_T1_x].mesonId=rhoxP_T1;
  WMesonStateTable[MS_rhoxP_T1_x].category=EXT_FIRSTDEV_MESON;
  WMesonStateTable[MS_rhoxP_T1_x].polarization=0;
  WMesonStateTable[MS_rhoxP_T1_x].measure=0;
  WMesonStateTable[MS_rhoxP_T1_x].sinkOp=MO_rhoxP_T1_x;
  WMesonStateTable[MS_rhoxP_T1_x].srcOp=MO_rhoxP_T1_x;
  
  WMesonStateTable[MS_rhoxP_T1_y]=WMesonStateTable[MS_rhoxP_T1_x];
  WMesonStateTable[MS_rhoxP_T1_y].polarization=1;
  WMesonStateTable[MS_rhoxP_T1_y].sinkOp=MO_rhoxP_T1_y;
  WMesonStateTable[MS_rhoxP_T1_y].srcOp=MO_rhoxP_T1_y;
  
  WMesonStateTable[MS_rhoxP_T1_z]=WMesonStateTable[MS_rhoxP_T1_x];
  WMesonStateTable[MS_rhoxP_T1_z].polarization=2;
  WMesonStateTable[MS_rhoxP_T1_z].sinkOp=MO_rhoxP_T1_z;
  WMesonStateTable[MS_rhoxP_T1_z].srcOp=MO_rhoxP_T1_z;
  
  //MS_rhoxP_T2_x, MS_rhoxP_T2_y, MS_rhoxP_T2_z,
  WMesonStateTable[MS_rhoxP_T2_x].stateName="rhoxP_T2";
  WMesonStateTable[MS_rhoxP_T2_x].mesonId=rhoxP_T2;
  WMesonStateTable[MS_rhoxP_T2_x].category=EXT_FIRSTDEV_MESON;
  WMesonStateTable[MS_rhoxP_T2_x].polarization=0;
  WMesonStateTable[MS_rhoxP_T2_x].measure=0;
  WMesonStateTable[MS_rhoxP_T2_x].sinkOp=MO_rhoxP_T2_x;
  WMesonStateTable[MS_rhoxP_T2_x].srcOp=MO_rhoxP_T2_x;
  
  WMesonStateTable[MS_rhoxP_T2_y]=WMesonStateTable[MS_rhoxP_T2_x];
  WMesonStateTable[MS_rhoxP_T2_y].polarization=1;
  WMesonStateTable[MS_rhoxP_T2_y].sinkOp=MO_rhoxP_T2_y;
  WMesonStateTable[MS_rhoxP_T2_y].srcOp=MO_rhoxP_T2_y;
  
  WMesonStateTable[MS_rhoxP_T2_z]=WMesonStateTable[MS_rhoxP_T2_x];
  WMesonStateTable[MS_rhoxP_T2_z].polarization=2;
  WMesonStateTable[MS_rhoxP_T2_z].sinkOp=MO_rhoxP_T2_z;
  WMesonStateTable[MS_rhoxP_T2_z].srcOp=MO_rhoxP_T2_z;

  //  MS_a1xP_A1 = 0--
  WMesonStateTable[MS_a1xP_A1_1].stateName="a1xP_A1";
  WMesonStateTable[MS_a1xP_A1_1].mesonId=a1xP_A1;
  WMesonStateTable[MS_a1xP_A1_1].category=EXT_FIRSTDEV_MESON;
  WMesonStateTable[MS_a1xP_A1_1].polarization=0;
  WMesonStateTable[MS_a1xP_A1_1].measure=0;
  WMesonStateTable[MS_a1xP_A1_1].sinkOp=MO_a1xP_A1;
  WMesonStateTable[MS_a1xP_A1_1].srcOp=MO_a1xP_A1;
  
  //MS_a1xP_T2_x, MS_a1xP_T2_y, MS_a1xP_T2_z,
  WMesonStateTable[MS_a1xP_T2_x].stateName="a1xP_T2";
  WMesonStateTable[MS_a1xP_T2_x].mesonId=a1xP_T2;
  WMesonStateTable[MS_a1xP_T2_x].category=EXT_FIRSTDEV_MESON;
  WMesonStateTable[MS_a1xP_T2_x].polarization=0;
  WMesonStateTable[MS_a1xP_T2_x].measure=0;
  WMesonStateTable[MS_a1xP_T2_x].sinkOp=MO_a1xP_T2_x;
  WMesonStateTable[MS_a1xP_T2_x].srcOp=MO_a1xP_T2_x;
  
  WMesonStateTable[MS_a1xP_T2_y]=WMesonStateTable[MS_a1xP_T2_x];
  WMesonStateTable[MS_a1xP_T2_y].polarization=1;
  WMesonStateTable[MS_a1xP_T2_y].sinkOp=MO_a1xP_T2_y;
  WMesonStateTable[MS_a1xP_T2_y].srcOp=MO_a1xP_T2_y;
  
  WMesonStateTable[MS_a1xP_T2_z]=WMesonStateTable[MS_a1xP_T2_x];
  WMesonStateTable[MS_a1xP_T2_z].polarization=2;
  WMesonStateTable[MS_a1xP_T2_z].sinkOp=MO_a1xP_T2_z;
  WMesonStateTable[MS_a1xP_T2_z].srcOp=MO_a1xP_T2_z;

  //MS_a1xP_E_1, MS_a1xP_E_2
  WMesonStateTable[MS_a1xP_E_1].stateName="a1xP_E";
  WMesonStateTable[MS_a1xP_E_1].mesonId=a1xP_E;
  WMesonStateTable[MS_a1xP_E_1].category=EXT_FIRSTDEV_MESON;
  WMesonStateTable[MS_a1xP_E_1].polarization=0;
  WMesonStateTable[MS_a1xP_E_1].measure=0;
  WMesonStateTable[MS_a1xP_E_1].sinkOp=MO_a1xP_E_1;
  WMesonStateTable[MS_a1xP_E_1].srcOp=MO_a1xP_E_1;
  
  WMesonStateTable[MS_a1xP_E_2]=WMesonStateTable[MS_a1xP_E_1];
  WMesonStateTable[MS_a1xP_E_2].polarization=1;
  WMesonStateTable[MS_a1xP_E_2].sinkOp=MO_a1xP_E_2;
  WMesonStateTable[MS_a1xP_E_2].srcOp=MO_a1xP_E_2;
  
  //MS_b1xP_T1_x, MS_b1xP_T1_y, MS_b1xP_T1_z,
  WMesonStateTable[MS_b1xP_T1_x].stateName="b1xP_T1";
  WMesonStateTable[MS_b1xP_T1_x].mesonId=b1xP_T1;
  WMesonStateTable[MS_b1xP_T1_x].category=EXT_FIRSTDEV_MESON;
  WMesonStateTable[MS_b1xP_T1_x].polarization=0;
  WMesonStateTable[MS_b1xP_T1_x].measure=0;
  WMesonStateTable[MS_b1xP_T1_x].sinkOp=MO_b1xP_T1_x;
  WMesonStateTable[MS_b1xP_T1_x].srcOp=MO_b1xP_T1_x;
  
  WMesonStateTable[MS_b1xP_T1_y]=WMesonStateTable[MS_b1xP_T1_x];
  WMesonStateTable[MS_b1xP_T1_y].polarization=1;
  WMesonStateTable[MS_b1xP_T1_y].sinkOp=MO_b1xP_T1_y;
  WMesonStateTable[MS_b1xP_T1_y].srcOp=MO_b1xP_T1_y;
  
  WMesonStateTable[MS_b1xP_T1_z]=WMesonStateTable[MS_b1xP_T1_x];
  WMesonStateTable[MS_b1xP_T1_z].polarization=2;
  WMesonStateTable[MS_b1xP_T1_z].sinkOp=MO_b1xP_T1_z;
  WMesonStateTable[MS_b1xP_T1_z].srcOp=MO_b1xP_T1_z;

  // MS_b1xD_A2 
  WMesonStateTable[MS_b1xD_A2_1].stateName="b1xD_A2";
  WMesonStateTable[MS_b1xD_A2_1].mesonId=b1xD_A2;
  WMesonStateTable[MS_b1xD_A2_1].category=EXT_SECONDDEV_SYM_MESON;
  WMesonStateTable[MS_b1xD_A2_1].polarization=0;
  WMesonStateTable[MS_b1xD_A2_1].measure=0;
  WMesonStateTable[MS_b1xD_A2_1].srcOp=MO_b1xD_A2;
  WMesonStateTable[MS_b1xD_A2_1].sinkOp=MO_b1xD_A2;

  // MS_b1xD_T1 

  WMesonStateTable[MS_b1xD_T1_x].stateName="b1xD_T1";
  WMesonStateTable[MS_b1xD_T1_x].mesonId=b1xD_T1;
  WMesonStateTable[MS_b1xD_T1_x].category=EXT_SECONDDEV_SYM_MESON;
  WMesonStateTable[MS_b1xD_T1_x].measure=0;
  WMesonStateTable[MS_b1xD_T1_x].polarization=0;
  WMesonStateTable[MS_b1xD_T1_x].srcOp=MO_b1xD_T1_x;
  WMesonStateTable[MS_b1xD_T1_x].sinkOp=MO_b1xD_T1_x;

  WMesonStateTable[MS_b1xD_T1_y]=WMesonStateTable[MS_b1xD_T1_x];
  WMesonStateTable[MS_b1xD_T1_y].polarization=1;
  WMesonStateTable[MS_b1xD_T1_y].srcOp=MO_b1xD_T1_y;
  WMesonStateTable[MS_b1xD_T1_y].sinkOp=MO_b1xD_T1_y;

  WMesonStateTable[MS_b1xD_T1_z]=WMesonStateTable[MS_b1xD_T1_x];
  WMesonStateTable[MS_b1xD_T1_z].polarization=2;
  WMesonStateTable[MS_b1xD_T1_z].srcOp=MO_b1xD_T1_z;
  WMesonStateTable[MS_b1xD_T1_z].sinkOp=MO_b1xD_T1_z;


  // MS_b1xD_T2 

  WMesonStateTable[MS_b1xD_T2_x].stateName="b1xD_T2";
  WMesonStateTable[MS_b1xD_T2_x].mesonId=b1xD_T2;
  WMesonStateTable[MS_b1xD_T2_x].category=EXT_SECONDDEV_SYM_MESON;
  WMesonStateTable[MS_b1xD_T2_x].measure=0;
  WMesonStateTable[MS_b1xD_T2_x].polarization=0;
  WMesonStateTable[MS_b1xD_T2_x].srcOp=MO_b1xD_T2_x;
  WMesonStateTable[MS_b1xD_T2_x].sinkOp=MO_b1xD_T2_x;

  WMesonStateTable[MS_b1xD_T2_y]=WMesonStateTable[MS_b1xD_T2_x];
  WMesonStateTable[MS_b1xD_T2_y].polarization=1;
  WMesonStateTable[MS_b1xD_T2_y].srcOp=MO_b1xD_T2_y;
  WMesonStateTable[MS_b1xD_T2_y].sinkOp=MO_b1xD_T2_y;

  WMesonStateTable[MS_b1xD_T2_z]=WMesonStateTable[MS_b1xD_T2_x];
  WMesonStateTable[MS_b1xD_T2_z].polarization=2;
  WMesonStateTable[MS_b1xD_T2_z].srcOp=MO_b1xD_T2_z;
  WMesonStateTable[MS_b1xD_T2_z].sinkOp=MO_b1xD_T2_z;

   //MS_b1xD_E_1, MS_b1xD_E_2
  WMesonStateTable[MS_b1xD_E_1].stateName="b1xD_E";
  WMesonStateTable[MS_b1xD_E_1].mesonId=b1xD_E;
  WMesonStateTable[MS_b1xD_E_1].category=EXT_SECONDDEV_SYM_MESON;
  WMesonStateTable[MS_b1xD_E_1].polarization=0;
  WMesonStateTable[MS_b1xD_E_1].measure=0;
  WMesonStateTable[MS_b1xD_E_1].sinkOp=MO_b1xD_E_1;
  WMesonStateTable[MS_b1xD_E_1].srcOp=MO_b1xD_E_1;
  
  WMesonStateTable[MS_b1xD_E_2]=WMesonStateTable[MS_b1xD_E_1];
  WMesonStateTable[MS_b1xD_E_2].polarization=1;
  WMesonStateTable[MS_b1xD_E_2].sinkOp=MO_b1xD_E_2;
  WMesonStateTable[MS_b1xD_E_2].srcOp=MO_b1xD_E_2;


  //MS_a0_primexD_x, MS_a0_primexD_y, MS_a0_primexD_z,
  WMesonStateTable[MS_a0_primexD_x].stateName="a0_primexD";
  WMesonStateTable[MS_a0_primexD_x].mesonId=a0_primexD;
  WMesonStateTable[MS_a0_primexD_x].category=EXT_SECONDDEV_SYM_MESON; 
  WMesonStateTable[MS_a0_primexD_x].polarization=0;
  WMesonStateTable[MS_a0_primexD_x].measure=0;
  WMesonStateTable[MS_a0_primexD_x].sinkOp=MO_a0_primexD_x;
  WMesonStateTable[MS_a0_primexD_x].srcOp=MO_a0_primexD_x;
  
  WMesonStateTable[MS_a0_primexD_y]=WMesonStateTable[MS_a0_primexD_x];
  WMesonStateTable[MS_a0_primexD_y].polarization=1;
  WMesonStateTable[MS_a0_primexD_y].sinkOp=MO_a0_primexD_y;
  WMesonStateTable[MS_a0_primexD_y].srcOp=MO_a0_primexD_y;
  
  WMesonStateTable[MS_a0_primexD_z]=WMesonStateTable[MS_a0_primexD_x];
  WMesonStateTable[MS_a0_primexD_z].polarization=2;
  WMesonStateTable[MS_a0_primexD_z].sinkOp=MO_a0_primexD_z;
  WMesonStateTable[MS_a0_primexD_z].srcOp=MO_a0_primexD_z;

  //MS_rhoxB_T1_x, MS_rhoxB_T1_y, MS_rhoxB_T1_z,
  WMesonStateTable[MS_rhoxB_T1_x].stateName="rhoxB_T1";
  WMesonStateTable[MS_rhoxB_T1_x].mesonId=rhoxB_T1;
  WMesonStateTable[MS_rhoxB_T1_x].category=EXT_SECONDDEV_ANTISYM_MESON; 
  WMesonStateTable[MS_rhoxB_T1_x].polarization=0;
  WMesonStateTable[MS_rhoxB_T1_x].measure=0;
  WMesonStateTable[MS_rhoxB_T1_x].sinkOp=MO_rhoxB_T1_x;
  WMesonStateTable[MS_rhoxB_T1_x].srcOp=MO_rhoxB_T1_x;
  
  WMesonStateTable[MS_rhoxB_T1_y]=WMesonStateTable[MS_rhoxB_T1_x];
  WMesonStateTable[MS_rhoxB_T1_y].polarization=1;
  WMesonStateTable[MS_rhoxB_T1_y].sinkOp=MO_rhoxB_T1_y;
  WMesonStateTable[MS_rhoxB_T1_y].srcOp=MO_rhoxB_T1_y;
  
  WMesonStateTable[MS_rhoxB_T1_z]=WMesonStateTable[MS_rhoxB_T1_x];
  WMesonStateTable[MS_rhoxB_T1_z].polarization=2;
  WMesonStateTable[MS_rhoxB_T1_z].sinkOp=MO_rhoxB_T1_z;
  WMesonStateTable[MS_rhoxB_T1_z].srcOp=MO_rhoxB_T1_z;


  //MS_rhoxB_T2_x, MS_rhoxB_T2_y, MS_rhoxB_T2_z,
  WMesonStateTable[MS_rhoxB_T2_x].stateName="rhoxB_T2";
  WMesonStateTable[MS_rhoxB_T2_x].mesonId=rhoxB_T2;
  WMesonStateTable[MS_rhoxB_T2_x].category=EXT_SECONDDEV_ANTISYM_MESON; 
  WMesonStateTable[MS_rhoxB_T2_x].polarization=0;
  WMesonStateTable[MS_rhoxB_T2_x].measure=0;
  WMesonStateTable[MS_rhoxB_T2_x].sinkOp=MO_rhoxB_T2_x;
  WMesonStateTable[MS_rhoxB_T2_x].srcOp=MO_rhoxB_T2_x;
  
  WMesonStateTable[MS_rhoxB_T2_y]=WMesonStateTable[MS_rhoxB_T2_x];
  WMesonStateTable[MS_rhoxB_T2_y].polarization=1;
  WMesonStateTable[MS_rhoxB_T2_y].sinkOp=MO_rhoxB_T2_y;
  WMesonStateTable[MS_rhoxB_T2_y].srcOp=MO_rhoxB_T2_y;
  
  WMesonStateTable[MS_rhoxB_T2_z]=WMesonStateTable[MS_rhoxB_T2_x];
  WMesonStateTable[MS_rhoxB_T2_z].polarization=2;
  WMesonStateTable[MS_rhoxB_T2_z].sinkOp=MO_rhoxB_T2_z;
  WMesonStateTable[MS_rhoxB_T2_z].srcOp=MO_rhoxB_T2_z;



  // MS_a1xB_A1

  WMesonStateTable[MS_a1xB_A1_1].stateName="a1xB_A1";
  WMesonStateTable[MS_a1xB_A1_1].mesonId=a1xB_A1;
  WMesonStateTable[MS_a1xB_A1_1].category=EXT_SECONDDEV_ANTISYM_MESON;
  WMesonStateTable[MS_a1xB_A1_1].polarization=0;
  WMesonStateTable[MS_a1xB_A1_1].measure=0;
  WMesonStateTable[MS_a1xB_A1_1].srcOp=MO_a1xB_A1;
  WMesonStateTable[MS_a1xB_A1_1].sinkOp=MO_a1xB_A1;


  //MS_a1xB_T1_x, MS_a1xB_T1_y, MS_a1xB_T1_z,
  WMesonStateTable[MS_a1xB_T1_x].stateName="a1xB_T1";
  WMesonStateTable[MS_a1xB_T1_x].mesonId=a1xB_T1;
  WMesonStateTable[MS_a1xB_T1_x].category=EXT_SECONDDEV_ANTISYM_MESON;
  WMesonStateTable[MS_a1xB_T1_x].polarization=0;
  WMesonStateTable[MS_a1xB_T1_x].measure=0;
  WMesonStateTable[MS_a1xB_T1_x].sinkOp=MO_a1xB_T1_x;
  WMesonStateTable[MS_a1xB_T1_x].srcOp=MO_a1xB_T1_x;
  
  WMesonStateTable[MS_a1xB_T1_y]=WMesonStateTable[MS_a1xB_T1_x];
  WMesonStateTable[MS_a1xB_T1_y].polarization=1;
  WMesonStateTable[MS_a1xB_T1_y].sinkOp=MO_a1xB_T1_y;
  WMesonStateTable[MS_a1xB_T1_y].srcOp=MO_a1xB_T1_y;
  
  WMesonStateTable[MS_a1xB_T1_z]=WMesonStateTable[MS_a1xB_T1_x];
  WMesonStateTable[MS_a1xB_T1_z].polarization=2;
  WMesonStateTable[MS_a1xB_T1_z].sinkOp=MO_a1xB_T1_z;
  WMesonStateTable[MS_a1xB_T1_z].srcOp=MO_a1xB_T1_z;


  //  MS_a1xB_T2_x, MS_a1xB_T2_y, MS_a1xB_T2_z, 
  WMesonStateTable[MS_a1xB_T2_x].stateName="a1xB_T2";
  WMesonStateTable[MS_a1xB_T2_x].mesonId=a1xB_T2;
  WMesonStateTable[MS_a1xB_T2_x].category=EXT_SECONDDEV_ANTISYM_MESON;
  WMesonStateTable[MS_a1xB_T2_x].polarization=0;
  WMesonStateTable[MS_a1xB_T2_x].measure=0;
  WMesonStateTable[MS_a1xB_T2_x].sinkOp=MO_a1xB_T2_x;
  WMesonStateTable[MS_a1xB_T2_x].srcOp=MO_a1xB_T2_x;
  
  WMesonStateTable[MS_a1xB_T2_y]=WMesonStateTable[MS_a1xB_T2_x];
  WMesonStateTable[MS_a1xB_T2_y].polarization=1;
  WMesonStateTable[MS_a1xB_T2_y].sinkOp=MO_a1xB_T2_y;
  WMesonStateTable[MS_a1xB_T2_y].srcOp=MO_a1xB_T2_y;
  
  WMesonStateTable[MS_a1xB_T2_z]=WMesonStateTable[MS_a1xB_T2_x];
  WMesonStateTable[MS_a1xB_T2_z].polarization=2;
  WMesonStateTable[MS_a1xB_T2_z].sinkOp=MO_a1xB_T2_z;
  WMesonStateTable[MS_a1xB_T2_z].srcOp=MO_a1xB_T2_z;

  // MS_a1xD_A2 
  WMesonStateTable[MS_a1xD_A2_1].stateName="a1xD_A2";
  WMesonStateTable[MS_a1xD_A2_1].mesonId=a1xD_A2;
  WMesonStateTable[MS_a1xD_A2_1].category=EXT_SECONDDEV_SYM_MESON;
  WMesonStateTable[MS_a1xD_A2_1].polarization=0;
  WMesonStateTable[MS_a1xD_A2_1].measure=0;
  WMesonStateTable[MS_a1xD_A2_1].srcOp=MO_a1xD_A2;
  WMesonStateTable[MS_a1xD_A2_1].sinkOp=MO_a1xD_A2;

  // MS_a1xD_T1 

  WMesonStateTable[MS_a1xD_T1_x].stateName="a1xD_T1";
  WMesonStateTable[MS_a1xD_T1_x].mesonId=a1xD_T1;
  WMesonStateTable[MS_a1xD_T1_x].category=EXT_SECONDDEV_SYM_MESON;
  WMesonStateTable[MS_a1xD_T1_x].measure=0;
  WMesonStateTable[MS_a1xD_T1_x].polarization=0;
  WMesonStateTable[MS_a1xD_T1_x].srcOp=MO_a1xD_T1_x;
  WMesonStateTable[MS_a1xD_T1_x].sinkOp=MO_a1xD_T1_x;

  WMesonStateTable[MS_a1xD_T1_y]=WMesonStateTable[MS_a1xD_T1_x];
  WMesonStateTable[MS_a1xD_T1_y].polarization=1;
  WMesonStateTable[MS_a1xD_T1_y].srcOp=MO_a1xD_T1_y;
  WMesonStateTable[MS_a1xD_T1_y].sinkOp=MO_a1xD_T1_y;

  WMesonStateTable[MS_a1xD_T1_z]=WMesonStateTable[MS_a1xD_T1_x];
  WMesonStateTable[MS_a1xD_T1_z].polarization=2;
  WMesonStateTable[MS_a1xD_T1_z].srcOp=MO_a1xD_T1_z;
  WMesonStateTable[MS_a1xD_T1_z].sinkOp=MO_a1xD_T1_z;


  // MS_a1xD_T2 

  WMesonStateTable[MS_a1xD_T2_x].stateName="a1xD_T2";
  WMesonStateTable[MS_a1xD_T2_x].mesonId=a1xD_T2;
  WMesonStateTable[MS_a1xD_T2_x].category=EXT_SECONDDEV_SYM_MESON;
  WMesonStateTable[MS_a1xD_T2_x].measure=0;
  WMesonStateTable[MS_a1xD_T2_x].polarization=0;
  WMesonStateTable[MS_a1xD_T2_x].srcOp=MO_a1xD_T2_x;
  WMesonStateTable[MS_a1xD_T2_x].sinkOp=MO_a1xD_T2_x;

  WMesonStateTable[MS_a1xD_T2_y]=WMesonStateTable[MS_a1xD_T2_x];
  WMesonStateTable[MS_a1xD_T2_y].polarization=1;
  WMesonStateTable[MS_a1xD_T2_y].srcOp=MO_a1xD_T2_y;
  WMesonStateTable[MS_a1xD_T2_y].sinkOp=MO_a1xD_T2_y;

  WMesonStateTable[MS_a1xD_T2_z]=WMesonStateTable[MS_a1xD_T2_x];
  WMesonStateTable[MS_a1xD_T2_z].polarization=2;
  WMesonStateTable[MS_a1xD_T2_z].srcOp=MO_a1xD_T2_z;
  WMesonStateTable[MS_a1xD_T2_z].sinkOp=MO_a1xD_T2_z;

   //MS_a1xD_E_1, MS_a1xD_E_2
  WMesonStateTable[MS_a1xD_E_1].stateName="a1xD_E";
  WMesonStateTable[MS_a1xD_E_1].mesonId=a1xD_E;
  WMesonStateTable[MS_a1xD_E_1].category=EXT_SECONDDEV_SYM_MESON;
  WMesonStateTable[MS_a1xD_E_1].polarization=0;
  WMesonStateTable[MS_a1xD_E_1].measure=0;
  WMesonStateTable[MS_a1xD_E_1].sinkOp=MO_a1xD_E_1;
  WMesonStateTable[MS_a1xD_E_1].srcOp=MO_a1xD_E_1;
  
  WMesonStateTable[MS_a1xD_E_2]=WMesonStateTable[MS_a1xD_E_1];
  WMesonStateTable[MS_a1xD_E_2].polarization=1;
  WMesonStateTable[MS_a1xD_E_2].sinkOp=MO_a1xD_E_2;
  WMesonStateTable[MS_a1xD_E_2].srcOp=MO_a1xD_E_2;

  // MS_rhoxD_A2 

  WMesonStateTable[MS_rhoxD_A2_1].stateName="rhoxD_A2";
  WMesonStateTable[MS_rhoxD_A2_1].mesonId=rhoxD_A2;
  WMesonStateTable[MS_rhoxD_A2_1].category=EXT_SECONDDEV_SYM_MESON;
  WMesonStateTable[MS_rhoxD_A2_1].polarization=0;
  WMesonStateTable[MS_rhoxD_A2_1].measure=0;
  WMesonStateTable[MS_rhoxD_A2_1].srcOp=MO_rhoxD_A2;
  WMesonStateTable[MS_rhoxD_A2_1].sinkOp=MO_rhoxD_A2;

  // MS_rhoxD_T1 

  WMesonStateTable[MS_rhoxD_T1_x].stateName="rhoxD_T1";
  WMesonStateTable[MS_rhoxD_T1_x].mesonId=rhoxD_T1;
  WMesonStateTable[MS_rhoxD_T1_x].category=EXT_SECONDDEV_SYM_MESON;
  WMesonStateTable[MS_rhoxD_T1_x].measure=0;
  WMesonStateTable[MS_rhoxD_T1_x].polarization=0;
  WMesonStateTable[MS_rhoxD_T1_x].srcOp=MO_rhoxD_T1_x;
  WMesonStateTable[MS_rhoxD_T1_x].sinkOp=MO_rhoxD_T1_x;

  WMesonStateTable[MS_rhoxD_T1_y]=WMesonStateTable[MS_rhoxD_T1_x];
  WMesonStateTable[MS_rhoxD_T1_y].polarization=1;
  WMesonStateTable[MS_rhoxD_T1_y].srcOp=MO_rhoxD_T1_y;
  WMesonStateTable[MS_rhoxD_T1_y].sinkOp=MO_rhoxD_T1_y;

  WMesonStateTable[MS_rhoxD_T1_z]=WMesonStateTable[MS_rhoxD_T1_x];
  WMesonStateTable[MS_rhoxD_T1_z].polarization=2;
  WMesonStateTable[MS_rhoxD_T1_z].srcOp=MO_rhoxD_T1_z;
  WMesonStateTable[MS_rhoxD_T1_z].sinkOp=MO_rhoxD_T1_z;


  // MS_rhoxD_T2 

  WMesonStateTable[MS_rhoxD_T2_x].stateName="rhoxD_T2";
  WMesonStateTable[MS_rhoxD_T2_x].mesonId=rhoxD_T2;
  WMesonStateTable[MS_rhoxD_T2_x].category=EXT_SECONDDEV_SYM_MESON;
  WMesonStateTable[MS_rhoxD_T2_x].measure=0;
  WMesonStateTable[MS_rhoxD_T2_x].polarization=0;
  WMesonStateTable[MS_rhoxD_T2_x].srcOp=MO_rhoxD_T2_x;
  WMesonStateTable[MS_rhoxD_T2_x].sinkOp=MO_rhoxD_T2_x;

  WMesonStateTable[MS_rhoxD_T2_y]=WMesonStateTable[MS_rhoxD_T2_x];
  WMesonStateTable[MS_rhoxD_T2_y].polarization=1;
  WMesonStateTable[MS_rhoxD_T2_y].srcOp=MO_rhoxD_T2_y;
  WMesonStateTable[MS_rhoxD_T2_y].sinkOp=MO_rhoxD_T2_y;

  WMesonStateTable[MS_rhoxD_T2_z]=WMesonStateTable[MS_rhoxD_T2_x];
  WMesonStateTable[MS_rhoxD_T2_z].polarization=2;
  WMesonStateTable[MS_rhoxD_T2_z].srcOp=MO_rhoxD_T2_z;
  WMesonStateTable[MS_rhoxD_T2_z].sinkOp=MO_rhoxD_T2_z;


  //added pionxB
  WMesonStateTable[MS_pionxB_T1_x].stateName="pionxB_T1";
  WMesonStateTable[MS_pionxB_T1_x].mesonId=pionxB_T1;
  WMesonStateTable[MS_pionxB_T1_x].category=EXT_SECONDDEV_ANTISYM_MESON;
  WMesonStateTable[MS_pionxB_T1_x].measure=0;
  WMesonStateTable[MS_pionxB_T1_x].polarization=0;
  WMesonStateTable[MS_pionxB_T1_x].srcOp=MO_pionxB_T1_x;
  WMesonStateTable[MS_pionxB_T1_x].sinkOp=MO_pionxB_T1_x;

  WMesonStateTable[MS_pionxB_T1_y]=WMesonStateTable[MS_pionxB_T1_x];
  WMesonStateTable[MS_pionxB_T1_y].polarization=1;
  WMesonStateTable[MS_pionxB_T1_y].srcOp=MO_pionxB_T1_y;
  WMesonStateTable[MS_pionxB_T1_y].sinkOp=MO_pionxB_T1_y;

  WMesonStateTable[MS_pionxB_T1_z]=WMesonStateTable[MS_pionxB_T1_x];
  WMesonStateTable[MS_pionxB_T1_z].polarization=2;
  WMesonStateTable[MS_pionxB_T1_z].srcOp=MO_pionxB_T1_z;
  WMesonStateTable[MS_pionxB_T1_z].sinkOp=MO_pionxB_T1_z;

  //pionxD_T2
  //MS_pionxD_T2_x, MS_pionxD_T2_y, MS_pionxD_T2_z,
  WMesonStateTable[MS_pionxD_T2_x].stateName="pionxD_T2";
  WMesonStateTable[MS_pionxD_T2_x].mesonId=pionxD_T2;
  WMesonStateTable[MS_pionxD_T2_x].category=EXT_SECONDDEV_SYM_MESON;
  WMesonStateTable[MS_pionxD_T2_x].polarization=0;
  WMesonStateTable[MS_pionxD_T2_x].measure=0;
  WMesonStateTable[MS_pionxD_T2_x].sinkOp=MO_pionxD_T2_x;
  WMesonStateTable[MS_pionxD_T2_x].srcOp=MO_pionxD_T2_x;
 
  WMesonStateTable[MS_pionxD_T2_y]=WMesonStateTable[MS_pionxD_T2_x];
  WMesonStateTable[MS_pionxD_T2_y].polarization=1;
  WMesonStateTable[MS_pionxD_T2_y].sinkOp=MO_pionxD_T2_y;
  WMesonStateTable[MS_pionxD_T2_y].srcOp=MO_pionxD_T2_y;
 
  WMesonStateTable[MS_pionxD_T2_z]=WMesonStateTable[MS_pionxD_T2_x];
  WMesonStateTable[MS_pionxD_T2_z].polarization=2;
  WMesonStateTable[MS_pionxD_T2_z].sinkOp=MO_pionxD_T2_z;
  WMesonStateTable[MS_pionxD_T2_z].srcOp=MO_pionxD_T2_z;
  
  // Mixing Analysis
  // add terms from mixing_terms.C here


  // end of mixing states

  //measure control
  //!!!!!!!!!!!!!Turn On/Off manually now!!!!!
  //Only do mixing at present

  for(int stateId=0;stateId<NUM_WMESON_STATE;stateId++){
    int cat=WMesonStateTable[stateId].category;
    if(cat==NORMALMESON){
      //WMesonStateTable[stateId].measure=measure_normal_mesons;
    }
    else if(cat==EXT_FIRSTDEV_MESON){
#ifdef DEBUG_W_EXT_MESON
     printf("InitTable: first dev on=%d\n",arg->extended_mesons_first_dev_on);
#endif	
      WMesonStateTable[stateId].measure=arg->extended_mesons_first_dev_on;
    }
    else if(cat==EXT_SECONDDEV_SYM_MESON){
#ifdef DEBUG_W_EXT_MESON
     printf("InitTable: second sym dev on=%d\n",arg->extended_mesons_second_sym_dev_on);
#endif
     WMesonStateTable[stateId].measure=arg->extended_mesons_second_sym_dev_on;
    }
    else if(cat==EXT_SECONDDEV_ANTISYM_MESON){
#ifdef DEBUG_W_EXT_MESON
     printf("InitTable: second antisym dev on=%d\n",arg->extended_mesons_second_antisym_dev_on);
#endif
      WMesonStateTable[stateId].measure=arg->extended_mesons_second_antisym_dev_on;
    }
    else if(cat==EXT_SECONDDEV_DIAG_MESON){
#ifdef DEBUG_W_EXT_MESON
      printf("InitTable: second diag dev on=%d\n",arg->extended_mesons_second_diag_dev_on);
#endif
      WMesonStateTable[stateId].measure=arg->extended_mesons_second_diag_dev_on;
    }
    else if(cat==MIXING){
#ifdef MIXING_ON
      WMesonStateTable[stateId].measure=1;
#endif
    }
    //if(WMesonStateTable[stateId].measure) printf("%d State measure==1\n");
  }
  //can have other categorization criterio, for example:
  //charge conj. number
  //final thing, still have chance to flag out some states explicitly for testing purpose
  
}










CPS_END_NAMESPACE
