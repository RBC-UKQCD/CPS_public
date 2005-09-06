#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------------
// alg_w_spc_ex.C
//
// AlgWspectExtMeson is derived from AlgWspect and is relevant to  
// meson spectroscopy with Wilson type fermions (i.e. Wilson, Clover, Dwf). 
// The type of fermion is determined by the argument to the 
// constructor but it should not be a non Wilson type lattice.
//--------------------------------------------------------------------------


CPS_END_NAMESPACE
#include <alg/alg_w_spect.h>         // class AlgWspect
#include <alg/common_arg.h>          // class CommonArg
#include <alg/w_spect_arg.h>         // class WspectArg
#include <util/lattice.h>            // class Lattice

#include <util/verbose.h>            // VRB
#include <util/error.h>              // ERR
#include <alg/w_all.h>                                 // class Wspect*
#include <stdlib.h>      // exit
CPS_START_NAMESPACE

#define DEBUG_ALG_W_SPECT_EXT_MESONS

//---------------------------------------------------------------------------
// For the purpose of debugging or timing during code upgrade
//---------------------------------------------------------------------------

//#define TIMING_ALG_W_SPECT

#ifdef  TIMING_ALG_W_SPECT
CPS_END_NAMESPACE
#include <time.h>                 // clock()
CPS_START_NAMESPACE
#endif


//--------------------------------------------------------------------------
// Static data members
//--------------------------------------------------------------------------
char * AlgWspectExtMeson::d_class_name = "AlgWspectExtMeson";

//---------------------------------------------------------------------------
// CTOR
//---------------------------------------------------------------------------
AlgWspectExtMeson::AlgWspectExtMeson(Lattice & latt, 
		  CommonArg *c_arg, 
		  WspectArg *arg,
				     CgArg *cg,
		  int n_quark_masses) : AlgWspect(latt,c_arg,arg,cg,n_quark_masses){

}

AlgWspectExtMeson::~AlgWspectExtMeson(){
}

//--------------------------------------------------------------------------
// run() override AlgWspect::run()
//--------------------------------------------------------------------------
void AlgWspectExtMeson::run()
{
#ifdef  TIMING_ALG_W_SPECT
  int quark_b, quark_e;
  int meson_b, meson_m, meson_e;
  int nucleon_b, nucleon_m, nucleon_e;
  int total_b = clock();
  int total_e;  
#endif

  char *fname = "run()";
  VRB.Func(d_class_name,fname);

  //common_arg->results is initialized in main program
  //for each individual source type/param
  WspectOutput * output = (WspectOutput *)common_arg->results;
  

  // Set the Lattice pointer
  //------------------------------------------------------------------------
  Lattice& lat = AlgLattice();

  int src_slice = d_arg_p->aots_start;
  int src_slice_step = d_arg_p->aots_step;
  int src_slice_end  = src_slice_step * d_arg_p->aots_num;
  CgArg cg;
  for ( ; src_slice < src_slice_end; src_slice += src_slice_step) {


    // Calculate quark propagator
    //----------------------------------------------------------------------
    // Ping:  certainly more work here to be done about the desired
    //        combinations of non-degenerate quarks.
    //        Presumably, more arguments will have to be passed in.
    //        One way: for three flavors, [100] means use only q1
    //                 to caculate spectrum.
    //        Also some care needed to get the scope (CTOR and DTOR) 
    //        of each quark propagator right.
    //    const WspectQuark & q2 = q;
    //    const WspectQuark & q3 = q;
    //Xiaodong&Thomas:
    //Modified to calculate extended mesons
    //q1 is the usual propagator(no source operator), which can be 
    //
    
    //used to pass propagation direction and src_slice infomation
    //to spectrum class
    WspectHyperRectangle hyperRect(d_arg_p->prop_dir, src_slice);    


    //create local quark propagator
    // phys_v3.11.4.xiaodong
    // WspectQuark q1(lat, output->cg, d_arg_p[0], hyperRect);
    // as from phys_v4.0.0
    WspectQuark q1(lat, output->cg, output->pbp, 
		   output->mid_point, output->a0_p, d_arg_p[0], cg,hyperRect);
    
    
    //Note: for ExtendedMesons, do only zero momentum projection
    WspectMomenta  mom(hyperRect, q1.SourceCenter2(), d_arg_p->num_mom - 1);
       
    
    //---------------------------------------------------------------------
    // Normal mesons
    //----------------------------------------------------------------------
    //added control by Thomas and Xiaodong
    //
    {
      if(d_arg_p->normal_mesons_on) {
	WspectMesons mes(q1, q1, hyperRect, mom);
	//write data to files
	mes.print(output);	
      }
    }//end of normal mesons 
    
    
    //----------------------------------------------------------------------
    //calculate extended mesons from Derivative operators
    //----------------------------------------------------------------------
    {
      //calculate additional quark propagator with source_operator
      
      if(d_arg_p->extended_mesons_on){
#ifdef DEBUG_ALG_W_EXT_SPECT
	printf("Starting extended meson...\n");
	printf("Operator groupID = %d\n",d_arg_p->extended_mesons_op_groupId);
#endif
	//extended mesons here

	//fuzzing object pointers, initialized to be zero
	WspectFuzzing *src_fuzz_p=0;
	WspectFuzzing *sink_fuzz_p=0;

	WspectExtendedMesons *ext_meson_p[MAX_FUZZING_C_NUM]; //for multiple sink_only fuzzing

	int src_fuzz_count, sink_fuzz_count;
	int src_fuzz_num=1;//default 1, even with out fuzzing(equivalent to fuzzing_c=Infi)
	int sink_fuzz_num=0;

	//set up src_sink_num and sink_src_num
	if(d_arg_p->fuzzing_on){
	  if(d_arg_p->sink_fuzzing_only){
	    src_fuzz_num=1;
	    sink_fuzz_num=d_arg_p->fuzzing_c_num;
	    //else keep default value 1
	  }else{
	    src_fuzz_num=d_arg_p->fuzzing_c_num;
	    sink_fuzz_num=1;
	  }
       	}else{
	  //no fuzzing
	  src_fuzz_num=1; 
	  sink_fuzz_num=1;
	}

	//create extended meson objects
	                                                            

	//loop over all source fuzzing parameters
	//Only one if no fuzzing or sink_fuzzing_only
	for(src_fuzz_count=0; src_fuzz_count<src_fuzz_num;src_fuzz_count++){

	  if(d_arg_p->fuzzing_on && d_arg_p->sink_fuzzing_only){
	    for(sink_fuzz_count=0;sink_fuzz_count<sink_fuzz_num;sink_fuzz_count++){
	      ext_meson_p[sink_fuzz_count]=new WspectExtendedMesons(d_arg_p,hyperRect,sink_fuzz_count);
	    }
	  }else{
	    ext_meson_p[0]=new WspectExtendedMesons(d_arg_p,hyperRect,src_fuzz_count);
	  }

	  //RUN fuzzing
	  if(d_arg_p->fuzzing_on){
#ifdef DEBUG_ALG_W_EXT_SPECT
	    printf("Start fuzzing at source...\n");
#endif
	    src_fuzz_p=new WspectFuzzing(lat,d_arg_p[0],0,src_fuzz_count);
	    src_fuzz_p->run();
	    sink_fuzz_p=src_fuzz_p;
#ifdef DEBUG_ALG_W_EXT_SPECT
	    printf("Fuzzing done\n");
#endif
	  }

	  //loop over all possible source_operators
	  int src_op;
	  
	  for(src_op=0;src_op<END_SUM_OP;src_op++){
	    if(ext_meson_p[0]->isInOpGroup(src_op,d_arg_p->extended_mesons_op_groupId)){
	      int measure=0;
	      //check if exists related states
	      
	      for(int state=0;state<NUM_WMESON_STATE;state++){
		IFloat weight=ext_meson_p[0]->table(state,-1, src_op, -1, -1);
		if(weight!=0) measure=1;
	      }
	      
	      WspectQuark *q2=0;

	      if(measure){

		if(d_arg_p->fuzzing_on && d_arg_p->sink_fuzzing_only){
		  
		  for(sink_fuzz_count=0;sink_fuzz_count<sink_fuzz_num; sink_fuzz_count++){
#ifdef DEBUG_ALG_W_EXT_SPECT
		    printf("sink fuzz only mode, sink_fuzz index #%d\n",sink_fuzz_count);
#endif
		    //if sink_fuzz_only, need to rerun fuzzing for each fuzzing parameter
		    //run fuzzing
		    if(!sink_fuzz_p){
		      sink_fuzz_p=new WspectFuzzing(lat,d_arg_p[0],0,sink_fuzz_count);
		      sink_fuzz_p->run();
		    }
		    
		    //collect
		     if(src_op!=UNIT){

		       // from phys_v3.11.4.xiaodong
		       // if(!q2) q2=new WspectQuark(lat, output->cg2, d_arg_p[0], hyperRect, (DEVOperatorKind)src_op,src_fuzz_p);//run once is enough!
		       
		       if(!q2) q2=new WspectQuark(lat, output->cg2, output->pbp, 
                  output->mid_point, output->a0_p, d_arg_p[0], cg,hyperRect, (DEVOperatorKind)src_op,src_fuzz_p );

		       ext_meson_p[sink_fuzz_count]->collect(q1,*q2,sink_fuzz_p);
		     }else{
		       ext_meson_p[sink_fuzz_count]->collect(q1,q1,sink_fuzz_p);
		     }
		     
		     delete(sink_fuzz_p);
		     sink_fuzz_p=0;
		     
		  }//for each sink_fuzzing param
		  
		}else{
		  //no fuzzing or src+sink fuzzing
		  if(src_op!=UNIT){
		    // as in phys_v3.11.4.xiaodong
		    // if(!q2) q2=new WspectQuark(lat, output->cg2, d_arg_p[0], hyperRect, (DEVOperatorKind)src_op,src_fuzz_p);
		    if(!q2) q2=new WspectQuark(lat, output->cg2, output->pbp, 
                  output->mid_point, output->a0_p, d_arg_p[0], cg,hyperRect, (DEVOperatorKind)src_op, src_fuzz_p );
		    //then call collect
		    ext_meson_p[0]->collect(q1,*q2,sink_fuzz_p);
		  }else{
		    ext_meson_p[0]->collect(q1,q1,sink_fuzz_p);
		  }

		}
	      	      
		
	      }//if(measure)
	      
	      if(q2) delete(q2);

	    }//if(srcOp is in operator Group)
	  }//for(src_op,...)  
	    
	  //call finish to do global sum and print out data, then delete extmeson objects
	  for(sink_fuzz_count=0;sink_fuzz_count<sink_fuzz_num;sink_fuzz_count++){
#ifdef DEBUG_ALG_W_EXT_SPECT
	    printf("Finish and print out extended mesons(fuzzing index #%d)\n",sink_fuzz_count);
#endif
	    ext_meson_p[sink_fuzz_count]->finish();
	    //print out all correlators(if non-zero) to files
	    ext_meson_p[sink_fuzz_count]->print();
	    delete(ext_meson_p[sink_fuzz_count]);
	  }
	  
	  if(src_fuzz_p) delete src_fuzz_p;
	  src_fuzz_p=0;
	  if(sink_fuzz_p && d_arg_p->sink_fuzzing_only) delete sink_fuzz_p;
	  sink_fuzz_p=0;
	  
	}//loop over source fuzzing parameters
#ifdef DEBUG_ALG_W_EXT_SPECT
	printf("Extended mesons done...\n");
#endif
      }//if ext_meson_on
	 
    }//end of extended mesons
    

    //----------------------------------------------------------
    // Extended Mesons from local B/E construction
    //----------------------------------------------------------
    {
      //loop over fuzzing param
      //for each fuzzing_id, run fuzzing, calculate B/E fields
      //create ext_mesonsBE object
      //loop over src_field operator, create second quark prop
      //  call ext_meson.collect
      //print out the result, delete ext_meson object
      //
       if(d_arg_p->extended_mesonsBE_on){
#ifdef DEBUG_ALG_W_EXT_SPECT
	printf("Starting extended meson BE ...\n");
	printf("Operator groupID = %d\n",d_arg_p->extended_mesonsBE_op_groupId);
#endif
	//extended mesons here

	//fuzzing object pointers, initialized to be zero
	WspectFuzzing *fuzz_p=0;
	WspectField *field_p=0;
	
	WspectExtendedMesonsBE *ext_mesonBE_p=0; 

	int fuzz_count;
	int fuzz_num=1; //no fuzzing or 1 fuzzing

	if(d_arg_p->BEfuzzing_on){
	  fuzz_num=d_arg_p->BEfuzzing_c_num;
	}
	
	
	//loop over all source fuzzing parameters
	for(fuzz_count=0; fuzz_count<fuzz_num;fuzz_count++){

	

	  //RUN fuzzing
	  if(d_arg_p->BEfuzzing_on){
#ifdef DEBUG_ALG_W_EXT_SPECT
	    printf("Start fuzzing.....\n");
#endif
	    fuzz_p=new WspectFuzzing(lat,d_arg_p[0],0,fuzz_count);
	    fuzz_p->run();
#ifdef DEBUG_ALG_W_EXT_SPECT
	    printf("Fuzzing done\n");
#endif
	  }

	  //Calculate Fields
	  field_p=new WspectField(lat,fuzz_p);
	  field_p->run();

	  //release fuzzed links
	  if(fuzz_p) delete(fuzz_p);

	  ext_mesonBE_p=new WspectExtendedMesonsBE(d_arg_p,hyperRect,fuzz_count,field_p);
	  

	  //loop over all possible source_operators
	  int src_op;
	  
	  for(src_op=BEGIN_BE_OP;src_op<END_BE_OP;src_op++){
	    if(ext_mesonBE_p->isInOpGroup(src_op,d_arg_p->extended_mesonsBE_op_groupId)){
	      int measure=0;
	      //check if exists related states
	      
	      for(int state=0;state<NUM_WEXTMESON_BE_STATES;state++){
		IFloat weight=ext_mesonBE_p->table(state,-1, src_op, -1, -1);
		if(weight!=0) measure=1;
	      }
	      
	      WspectQuark *q2=0;

	      if(measure){

		//collect
		if(src_op!=UNIT){
		  // as in phys_v3.11.4.xiaodong
		  // q2=new WspectQuark(lat, output->cg2, d_arg_p[0], hyperRect, (DEVOperatorKind)src_op,0,field_p);
		  q2=new WspectQuark(lat, output->cg2, output->pbp, 
                  output->mid_point, output->a0_p, d_arg_p[0], cg, hyperRect, (DEVOperatorKind)src_op, 0, field_p );
		       
		  ext_mesonBE_p->collect(q1,*q2);
		}else{
		  ext_mesonBE_p->collect(q1,q1);
		}
		
	      }//if(measure)
	      
	      delete(q2);

	    }//if(srcOp is in operator Group)
	  }//for(src_op,...)  
	    
	  //call finish to do global sum and print out data, then delete extmeson objects
#ifdef DEBUG_ALG_W_EXT_SPECT
	  printf("Finish and print out BE extended mesons(fuzzing index #%d)\n",fuzz_count);
#endif
	  ext_mesonBE_p->finish();
	  //print out all correlators(if non-zero) to files
	  ext_mesonBE_p->print();
	  delete(ext_mesonBE_p);
       
	  if(fuzz_p) delete fuzz_p;
	  
	}//loop over source fuzzing parameters
#ifdef DEBUG_ALG_W_EXT_SPECT
	printf("Extended BE mesons done...\n");
#endif
       }//if ext_mesonBE_on
    }
   
    // Increment the counter
    d_counter += d_count_step;
  }//end of for(sc_slice,..)
  

}


CPS_END_NAMESPACE
