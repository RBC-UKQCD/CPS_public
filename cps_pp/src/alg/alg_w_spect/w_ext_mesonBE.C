#include<config.h>
CPS_START_NAMESPACE
/*! \file

  $Id: w_ext_mesonBE.C,v 1.11 2008/02/08 18:35:05 chulwoo Exp $
*/  

/*w_ext_mesonBE.C
 *
 */

CPS_END_NAMESPACE
#include <alg/w_all.h>
#include <alg/w_gamma_mat.h>
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
#include <util/error.h>                // ERR
#include <util/qcdio.h>                // ERR
#include <util/verbose.h>              // VRB
#include <util/vector.h>               // dotProduct
#include <comms/glb.h>                   // glb_sum
#include <alg/alg_w_spect.h>          // AlgWspect::GetCounter()
#include <comms/cbuf.h> 

CPS_START_NAMESPACE

//for debugging
#define DEBUG_W_EXT_MESON_BE

#define MATRIX_SIZE 18



//+++++++++++++++++++++++++++Major Functions++++++++++++++++++++++

//======================================================================
//Static Member Declarations
//======================================================================
char *WspectExtendedMesonsBE::d_class_name="WspectExtendedMesonsBE";   
int WspectExtendedMesonsBE::tableInitialized=0;
struct WMesonBEOpInfo
      WspectExtendedMesonsBE::WMesonOpTable[NUM_WEXTMESON_BE_OPS];
struct WMesonBEStateInfo
      WspectExtendedMesonsBE::WMesonStateTable[NUM_WEXTMESON_BE_STATES];

//-------------------
// CTOR
//-------------------
WspectExtendedMesonsBE::WspectExtendedMesonsBE(WspectArg *w_arg_p,
					       const WspectHyperRectangle & whr, 
					       const int fuzzing_id,
					       const WspectField *field_ptr)
  : WspectExtendedMesons(w_arg_p,whr,fuzzing_id),field_p(field_ptr){

  //-----------------------------------
  //allocate memory for correlators
  //-----------------------------------
  d_size   = d_glb_walls* NUM_WEXTMESON_BE_STATES * COMPLEXs; //correlator array in Float
  d_output_size   = d_glb_walls* NUM_WEXTMESON_BE_OUTPUT * COMPLEXs; //corrlator(polarisation average)in Float
  
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

 
  

  {
    //initialize map
    //!!!!!!!!!!!!!!!!!!!!!Need to be modified to make prop!=3 work!!!!!!!
    //for now assume prop_dir==3 only
    //should do manually mapping from coor_data_p (calculated aussuming prop==3)
    //to output data files later
    for(int i=0;i<NUM_WEXTMESON_BE_STATES;i++){
      map[i]=i; //if prop_dir=3(T)
    }
  }
  
  //initialize the stateTable according to spect argument
  if(!tableInitialized){
    initWMesonOpTable();
    initWMesonStateTable(arg_p);
    tableInitialized=1;
  }
  
}

//--------------------
// DTOR
//--------------------
WspectExtendedMesonsBE::~WspectExtendedMesonsBE(){
//deallocate memory
  VRB.Func(d_class_name, dtor_str);

  VRB.Sfree(d_class_name, dtor_str, empty_str, coor_data_p);
  sfree(coor_data_p);
  
  
  VRB.Sfree(d_class_name, dtor_str, empty_str, coor_output_p);
  sfree(coor_output_p);
  

  
}



//-----------------
// collect
//-----------------
void WspectExtendedMesonsBE::collect(const WspectQuark &q_l, WspectQuark &q_nl){
    
  char *fname="collect";
  VRB.Func(d_class_name, fname);
  
  DEVOperatorKind sour_op;
  sour_op=q_nl.srcOpKind();
 
#ifdef  DEBUG_W_EXT_MESON_BE
  printf("colleting BE extmesons - source operator = %d \n", sour_op);
#endif
 
  const Float *ql_p=(Float *)q_l.Data(); //local propagator
  const Float *qnl_p=(Float *)q_nl.Data(); //non_local propagator

  for(int sink_op=BEGIN_BE_OP+1; sink_op < END_BE_OP; sink_op++){ 
    //check 
    Float weight_test=0.;
    int state=0;

    while (weight_test==0. && state <NUM_WEXTMESON_BE_STATES ){
      weight_test = table(state, -1, sour_op, -1, sink_op);
      state += 1; // test next state
    } // endwhile (weight_test != 0)
    
    if (weight_test == 0.) continue; 

    int lclWalls=lcl_sites[d_prop_dir];
    for(int lclW=0;lclW<lclWalls;lclW++){
      DiracAlgebra(ql_p,qnl_p,lclW, sour_op, sink_op); 
    }
  }
  return;
  
}

//-------------------
// doAllAlgebra
//--------------------


//-------------------
// DiracAlgebra
//--------------------


//---------------------
// finish
//--------------------
void WspectExtendedMesonsBE::finish(){
  //do global sum of Complex coor_data_p[glbwall][meson]
#ifdef DEBUG_W_EXT_MESON_BE
  printf("finish().Do global sum.\n");
#endif
  Float *Flt_p = (Float *)coor_data_p;
  for (int i = 0; i < d_size; ++i)    glb_sum(Flt_p++);
  
}

//--------------------
// print
//--------------------
void WspectExtendedMesonsBE::print()const{

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
      int offset_data_wall=wall*NUM_WEXTMESON_BE_STATES;
      int offset_output_wall=wall*NUM_WEXTMESON_BE_OUTPUT;
      int num_pol; //count number of polarization


      //meason is real name, must use map
      for(int mesonId=0;mesonId<NUM_WEXTMESON_BE_OUTPUT;mesonId++){
        //loop over all states to find all polarizations of the same
        //meson
        num_pol=0;
        int output_offset=offset_output_wall+mesonId;
        for(int mesonState=0;mesonState<NUM_WEXTMESON_BE_STATES;mesonState++){

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
  for (int meson = 0; meson < NUM_WEXTMESON_BE_OUTPUT; ++meson) {
    
    //    if(measured[meson]){
    int measure_flag=0;
    //generate the output filename
    //get statename from WMesonStateTable
    for(int state=0;state<NUM_WEXTMESON_BE_STATES;state++){
      if(WMesonStateTable[state].mesonId==meson){
	if(WMesonStateTable[state].measure) {
	  measure_flag=1;
	  if(!arg_p->fuzzing_on){
	    //	    sprintf(outputfilename,"%s.%s",WMesonStateTable[state].stateName,arg_p->filetail);
      	    sprintf(outputfilename,"%s",WMesonStateTable[state].stateName);
	  }else{
	    //fuzzing is on
	    // sprintf(outputfilename,"%s.%sF%.1f",WMesonStateTable[state].stateName,arg_p->filetail,arg_p->fuzzing_c[fuzzing_c_index]);
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
	
	Complex tmp = coor_output_p[meson + NUM_WEXTMESON_BE_OUTPUT*
				 ((d_whr.glbCoord() + wall)%d_glb_walls)];
	tmp += coor_output_p[meson +  NUM_WEXTMESON_BE_OUTPUT*
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

//++++++++++++++++++++End of Major Functions+++++++++++++++++

//++++++++++++++++++++Beginning of Utility Functions+++++++++

Float WspectExtendedMesonsBE::table(int state, int sour_gamma, int sour_op, int sink_gamma, int sink_int) const{
  return 0;

}

int WspectExtendedMesonsBE::isInOpGroup(int op, int groupId) const{
   char *fname="isInOpGroup";
   int result=0;
   switch(groupId){
   case 0:
     if(op>=FB1_OP && op<=FUNIT_OP) result=1;
     break;
   case 1:
     if(op==SUM_MAGN_OP || op==SUM_ELEC_OP || op == FUNIT) result=1;
     break;
   case 2:
     if(op==SUM_MAGN_ELEC_OP ||
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

int WspectExtendedMesonsBE::matchSUMOp(int op, int sum_op) const{

  return 0;
}


//--------------------
// initStateTable
//--------------------
void WspectExtendedMesonsBE::initWMesonStateTable(WspectArg *arg){
  //pionxB
  WMesonStateTable[BE_MS_pionxB_x].stateName="BEpionxB";
  WMesonStateTable[BE_MS_pionxB_x].mesonId=BE_pionxB;
  WMesonStateTable[BE_MS_pionxB_x].category=MAG_HYBRID_BE;
  WMesonStateTable[BE_MS_pionxB_x].polarization=0;
  WMesonStateTable[BE_MS_pionxB_x].measure=0;
  WMesonStateTable[BE_MS_pionxB_x].srcOp=BE_MO_pionxB_x;
  WMesonStateTable[BE_MS_pionxB_x].sinkOp=BE_MO_pionxB_x;
  
  WMesonStateTable[BE_MS_pionxB_y]=WMesonStateTable[BE_MS_pionxB_x];
  WMesonStateTable[BE_MS_pionxB_y].polarization=1;
  WMesonStateTable[BE_MS_pionxB_y].srcOp=BE_MO_pionxB_y;
  WMesonStateTable[BE_MS_pionxB_y].sinkOp=BE_MO_pionxB_y;
  
  WMesonStateTable[BE_MS_pionxB_z]=WMesonStateTable[BE_MS_pionxB_x];
  WMesonStateTable[BE_MS_pionxB_z].polarization=1;
  WMesonStateTable[BE_MS_pionxB_z].srcOp=BE_MO_pionxB_z;
  WMesonStateTable[BE_MS_pionxB_z].sinkOp=BE_MO_pionxB_z;

  //rhoxB_T1
  WMesonStateTable[BE_MS_rhoxB_T1_x].stateName="BErhoxB_T1";
  WMesonStateTable[BE_MS_rhoxB_T1_x].mesonId=BE_rhoxB_T1;
  WMesonStateTable[BE_MS_rhoxB_T1_x].category=MAG_HYBRID_BE;
  WMesonStateTable[BE_MS_rhoxB_T1_x].polarization=0;
  WMesonStateTable[BE_MS_rhoxB_T1_x].measure=0;
  WMesonStateTable[BE_MS_rhoxB_T1_x].srcOp=BE_MO_rhoxB_T1_x;
  WMesonStateTable[BE_MS_rhoxB_T1_x].sinkOp=BE_MO_rhoxB_T1_x;
  
  WMesonStateTable[BE_MS_rhoxB_T1_y]=WMesonStateTable[BE_MS_rhoxB_T1_x];
  WMesonStateTable[BE_MS_rhoxB_T1_y].polarization=1;
  WMesonStateTable[BE_MS_rhoxB_T1_y].srcOp=BE_MO_rhoxB_T1_y;
  WMesonStateTable[BE_MS_rhoxB_T1_y].sinkOp=BE_MO_rhoxB_T1_y;
  
  WMesonStateTable[BE_MS_rhoxB_T1_z]=WMesonStateTable[BE_MS_rhoxB_T1_x];
  WMesonStateTable[BE_MS_rhoxB_T1_z].polarization=1;
  WMesonStateTable[BE_MS_rhoxB_T1_z].srcOp=BE_MO_rhoxB_T1_z;
  WMesonStateTable[BE_MS_rhoxB_T1_z].sinkOp=BE_MO_rhoxB_T1_z;

  //change measure property according to argument
  //!!!!!!Turn everything on for now
  for(int stateId=0;stateId<NUM_WEXTMESON_BE_STATES;stateId++){
    int cat=WMesonStateTable[stateId].category;
    WMesonStateTable[stateId].measure=1;
  }
}

//--------------------
// initOpTable
//--------------------
void WspectExtendedMesonsBE::initWMesonOpTable(){

  //pionxB
  WMesonOpTable[BE_MO_pionxB_x].num_terms=1;
  setWMesonOpTerm(&(WMesonOpTable[BE_MO_pionxB_x].terms[0][0]),1,WGAM_5,FB1);
  
  WMesonOpTable[BE_MO_pionxB_y].num_terms=1;
  setWMesonOpTerm(&(WMesonOpTable[BE_MO_pionxB_y].terms[0][0]),1,WGAM_5,FB2);
  
  WMesonOpTable[BE_MO_pionxB_z].num_terms=1;
  setWMesonOpTerm(&(WMesonOpTable[BE_MO_pionxB_z].terms[0][0]),1,WGAM_5,FB3);
  
  //rhoxB_T1
  WMesonOpTable[BE_MO_rhoxB_T1_x].num_terms=2;
  setWMesonOpTerm(&(WMesonOpTable[BE_MO_rhoxB_T1_x].terms[0][0]),1,WGAM_2,FB3);
  setWMesonOpTerm(&(WMesonOpTable[BE_MO_rhoxB_T1_x].terms[1][0]),-1,WGAM_3,FB2);

  WMesonOpTable[BE_MO_rhoxB_T1_y].num_terms=2;
  setWMesonOpTerm(&(WMesonOpTable[BE_MO_rhoxB_T1_y].terms[0][0]),1,WGAM_3,FB1);
  setWMesonOpTerm(&(WMesonOpTable[BE_MO_rhoxB_T1_y].terms[1][0]),-1,WGAM_1,FB3);
  
  WMesonOpTable[BE_MO_rhoxB_T1_z].num_terms=2;
  setWMesonOpTerm(&(WMesonOpTable[BE_MO_rhoxB_T1_z].terms[0][0]),1,WGAM_1,FB2);
  setWMesonOpTerm(&(WMesonOpTable[BE_MO_rhoxB_T1_z].terms[1][0]),-1,WGAM_2,FB1);
  
 
}

//---------------------
//setWMesonOpTerm
//---------------------
void WspectExtendedMesonsBE::setWMesonOpTerm(int *term_p, int weight, WGammaMatrix gammaMat, FieldTensorId fieldId){
  term_p[0]=weight;
  term_p[1]=gammaMat;
  term_p[2]=fieldId;
}



CPS_END_NAMESPACE
