#include<config.h>
CPS_START_NAMESPACE
/*w_quark_d.C
* new derivative functions for quark propagators
* Xiaodong & Thomas
*/

CPS_END_NAMESPACE
#include <alg/w_all.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <stdlib.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <util/gjp.h>              // GJP
#include <util/error.h>            // ERR
#include <util/verbose.h>          // VRB
#include <util/vector.h>
#include <util/lattice.h>          // Lattice::GetLink()
#include <comms/scu.h>               //getMinusData, getPlusData
#include <comms/nga_reg.h>
#include <comms/cbuf.h>
CPS_START_NAMESPACE

//for debugging
//#define DEBUG_WQUARK_D
//#define DEBUG_CHECKSU3

//Used to turn on/off Circular buffer access
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

//Note: The actual direction of derivative depends on propagation direction
//      See end of w_quark.h for detail
//---------------------------------------------------------------------------
// void WspectQuark::doSinkOperator(...)
//---------------------------------------------------------------------------
//assume memory has been allocated for prop_out

void WspectQuark::doSinkOperator(int lclW, DEVOperatorKind sink_op_kind, Float*prop_out, Float *prop_tmp,WspectFuzzing *sink_fuzz_ptr){
  char *fname="doSinkOperator";
  VRB.Func(d_class_name, fname);
  
  if(!prop_out){
    //error: must allocate memory for prop_out before calling this function
    ERR.Pointer(d_class_name, fname, "prop_out");
  }

  //set pointer to fuzzing links object

  sink_fuzz_p = sink_fuzz_ptr;

  /*Note: for second derivatives,need extra memory to store result from first derivative!
   *may run out of memory! Release immediately after usage*/
  int pol_dir; //actual derivative direction of DEVI(I=1,2,3)
  switch(sink_op_kind)
    {
    case UNIT:
      //extract propagator on a single
      //time slice from (Float *)d_data_p, and store in prop_out
      pol_dir=-1;//UNIT in sink_dev
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,1,(Float *)d_data_p, prop_out,pol_dir);
      break;
      
    case DEV1:
      //pol_dir is determined from propagation direction
      //At present, pol_dir can't be the same as prop_direction
      pol_dir=devDirToPolDir(1,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,1,(Float *)d_data_p, prop_out,pol_dir);
      break;
    case DEV2:
      pol_dir=devDirToPolDir(2,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,1,(Float *)d_data_p, prop_out,pol_dir);
      break;
    case DEV3:
      pol_dir=devDirToPolDir(3,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,1,(Float *)d_data_p, prop_out,pol_dir);
      break;
    case DEV1DEV2:
      pol_dir=devDirToPolDir(2,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,1,(Float *)d_data_p, prop_tmp,pol_dir);
      pol_dir=devDirToPolDir(1,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,0,prop_tmp, prop_out,pol_dir); 
      break;
    case DEV2DEV1:
      pol_dir=devDirToPolDir(1,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,1,(Float *)d_data_p, prop_tmp,pol_dir);
      pol_dir=devDirToPolDir(2,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,0,prop_tmp, prop_out,pol_dir); 
      break;
    case DEV2DEV3:
      pol_dir=devDirToPolDir(3,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,1,(Float *)d_data_p, prop_tmp,pol_dir);
      pol_dir=devDirToPolDir(2,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,0,prop_tmp, prop_out,pol_dir); 
      break;
    case DEV3DEV2:
      pol_dir=devDirToPolDir(2,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,1,(Float *)d_data_p, prop_tmp,pol_dir);
      pol_dir=devDirToPolDir(3,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,0,prop_tmp, prop_out,pol_dir); 
      break;
    case DEV1DEV3:
      pol_dir=devDirToPolDir(3,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,1,(Float *)d_data_p, prop_tmp,pol_dir);
      pol_dir=devDirToPolDir(1,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,0,prop_tmp, prop_out,pol_dir); 
      break;
    case DEV3DEV1:
      pol_dir=devDirToPolDir(1,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,1,(Float *)d_data_p, prop_tmp,pol_dir);
      pol_dir=devDirToPolDir(3,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,0,prop_tmp, prop_out,pol_dir); 
      break;
    case DEV1DEV1:
      pol_dir=devDirToPolDir(1,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,1,(Float *)d_data_p, prop_tmp,pol_dir);
      pol_dir=devDirToPolDir(1,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,0,prop_tmp, prop_out,pol_dir);
      break;

    case DEV2DEV2:
      pol_dir=devDirToPolDir(2,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,1,(Float *)d_data_p, prop_tmp,pol_dir);
      pol_dir=devDirToPolDir(2,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,0,prop_tmp, prop_out,pol_dir);
      break;

    case DEV3DEV3:
      pol_dir=devDirToPolDir(3,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,1,(Float *)d_data_p, prop_tmp,pol_dir);
      pol_dir=devDirToPolDir(3,prop_direction);
#ifdef DEBUG_WQUARK_D
      printf("call sink_dev for direction %d\n",pol_dir);
#endif
      sink_deriv(lclW,0,prop_tmp, prop_out,pol_dir);
      break;
    
    default:
      ERR.General(d_class_name, fname, "Unknown Sink Operator!");
      //error!
    }
}

//---------------------------------------------------------------------------
// void WspectQuark::sink_deriv(...)
//---------------------------------------------------------------------------
//dir=0(X),1(Y),2(Z),3(T)
//if dir<0, then do slice extraction only
void WspectQuark::sink_deriv(int lclW, int isFullProp,const Float *prop_in, Float*prop_out, int dir) const{
//Out=Deriv1 prop_in
/*Note: quark propagator's storage order is:
   *[y_Dirac][y_Color][x_T][x_Z][x_Y][x_X][x_Dirac][x_Color][Complex]
   *Deriv1 operates in sink indices(x) only
   *Deriv1 prop_in[Cx=C1] = SumOverC2{ULink(x,1)[C1][C2]*prop[Cx=C2]-
   *                        SumOverC2{Dag(ULink(x-1,1)[C1,C2]*prop[Cx=C2]
   */
  char *fname="sink_deriv"; 


  /* set cbuf registers
#ifdef USE_CBUF
  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(4, CBUF_MODE4);
#endif
  */

  //local dimensions are in WspectGinfo::lcl_sites[LORENTZs]
  int offset;

  //loop variables
  int Dy,Cy; //source Dirac and Color index
  int site[LORENTZs]; //sink site loop index
  int siteMin[LORENTZs];
  int siteMax[LORENTZs];
  
  int site_plus[LORENTZs]; //positive shift
  int site_minus[LORENTZs]; //negative shift
  int Dx; //sink Dirac
  int lclWalls=lcl_sites[prop_direction];
  int sink_prop_size=d_size/lcl_sites[prop_direction];
  
  if(dir==prop_direction) ERR.General(d_class_name,fname,"dir can't be prop_dir!");

  for(int i=0;i<LORENTZs;i++){
    siteMin[i]=0;
    siteMax[i]=lcl_sites[i]-1;
    if(i==prop_direction) siteMin[i]=siteMax[i]=lclW;
  }

  //calculate derivative
  IFloat *out; //complex 3D vector
  const IFloat *in;
  Vector in1[DIRACs]; //complex 3D vector
  Vector in2[DIRACs]; 

  const Matrix *glink;
  Matrix glink1;
  Matrix glink2;
  Matrix glink_dag; //store Dagger of glink
  Matrix glink_trans; //transpose of glink
  
  //loop over prop_out indices:  
  for(site[3]=siteMin[3];site[3]<=siteMax[3];site[3]++){
    for(site[2]=siteMin[2];site[2]<=siteMax[2];site[2]++){
      for(site[1]=siteMin[1];site[1]<=siteMax[1];site[1]++){
	for(site[0]=siteMin[0];site[0]<=siteMax[0];site[0]++){
	  
	  if(dir<0){
	    //UNIT operator, copy in to out
	    //careful if in is a full propagator, must do conversion
	    //SINK_UNIT, use sink_dev(....,dir=-1)
	    if(!isFullProp){
	      moveMem(out,in,sink_prop_size);
	    }else{
	      //convert full prop to single slice prop
	      for(Dy=0; Dy<DIRACs; Dy++){
		for(Cy=0;Cy<COLORs;Cy++){
		  for(Dx=0; Dx<DIRACs; Dx++){
		    for(int Cx=0;Cx<COLORs;Cx++){
		      int outoffset=Cx*COMPLEXs+COMPLEXs*COLORS*Dx+SPINORs*siteOffset(site,prop_direction)+Cy*weightSrcColor()/lclWalls+Dy*weightSrcDirac()/lclWalls;
		      int inoffset=Cx*COMPLEXs+COMPLEXs*COLORS*Dx+SPINORs*siteOffset(site)+Cy*weightSrcColor()+Dy*weightSrcDirac();
		      *(out+outoffset)=*(in+inoffset);//real
		      *(out+outoffset+1)=*(in+inoffset+1);//imaginary
		    }
		  }
		}
	      }
	    }
	    //done with UNIT operator
	  }else{
	    //calculate site_plus and site_minus
	    //get the gauge links
	    for(int i=0;i<LORENTZs;i++){
	      if(i==dir){
		site_plus[i]=site[i]+1; //??? if > maxsite? periodic BC?
		site_minus[i]=site[i]-1;
	      }else{
		site_plus[i]=site[i];
		site_minus[i]=site[i];
	      }
	    }//i
	    
	    //sink_op != SINK_UNIT
	    //get U_i(x)
	    if(sink_fuzz_p==0){
	      glink=d_lat.GetLink(site, dir);
	    }else{
	      glink=sink_fuzz_p->GetLink(site, dir);
	    }
#ifdef DEBUG_CHECKSU3
	    CheckSU3((const Float *)glink);
#endif
	  
	    glink1=*glink;
	    
	    //get U_i(x-i)
	    if(sink_fuzz_p==0){
	      glink=d_lat.GetLink(site_minus, dir);
	    }else{
	      glink=sink_fuzz_p->GetLink(site_minus, dir);
	    }
	  
	    glink2=*glink;
	    	    
	    //apply the derivative
	    for(Dy=0; Dy<DIRACs; Dy++){
	      for(Cy=0;Cy<COLORs;Cy++){
		//get prop[x+i]
		in=(IFloat *)getPropData(isFullProp,(IFloat *)prop_in,site_plus,Dy,Cy,-1);
		moveMem((IFloat *)in1,in,sizeof(Vector)*DIRACs);
		in=(IFloat *)getPropData(isFullProp,(IFloat *)prop_in,site_minus,Dy,Cy,-1);
		moveMem((IFloat *)in2,in,sizeof(Vector)*DIRACs);
		for(Dx=0; Dx<DIRACs; Dx++){
		  
		  //set offset in output propagator
		  offset=COMPLEXs*COLORS*Dx+SPINORs*siteOffset(site,prop_direction)+Cy*weightSrcColor()/lclWalls+Dy*weightSrcDirac()/lclWalls;
		  out=(IFloat *)(prop_out+offset); //point to output 3 complex vector
		  
		  //out=glink*in=U_i(x)*prop(x+i)
		  uDotXEqual(out,(const IFloat *)&glink1,(const IFloat *)&in1[Dx]);
		  
		  //out=out-Dagger(glink)*in==U_i(x)*prop(x+i)-U_i~(x-i)*prop(x-i)
		  glink_dag.Dagger((const IFloat *)&glink2);
		  uDotXMinus(out,(const IFloat *)&glink_dag,(const IFloat *)&in2[Dx]);
		  //dummay read
		  //moveMem((IFloat *)&glink_dag, (const IFloat*)glink+MYBANK4_BASE+MYBANK_SIZE,sizeof(glink_dag));
		 
		}//Dx
	      }//Cy
	    }//Dy
	  }//else(dir>==0)
	}//site[0]
      }//site[1]
    }//site[2]
  }//site[3]
}//finish derivative calculation

//---------------------------------------------------------
//const Vector *WspectQuark::getPropData(int isFullProp,const IFloat *prop_data_p, 
//                          const int *site, int Dy, int Cy, int Dx) const
//
//If Dx>=0 return Vector [Cx][Re,Im] (Color vector, 3x2) 
//If Dx<= return Vector[Dx][Cx][Re,Im]     (a whole spinor,4x3x2) to reduce internode communication
//---------------------------------------------------------

const Vector *WspectQuark::getPropData(int isFullProp,const IFloat *prop_data_p, const int *site, int Dy, int Cy, int Dx) const{  

  // offset out-of-range coordinates site[] into on_node_site[]
  // in order to locate the link
  //-----------------------------------------------------

  int on_node_site[4];
  int on_node = 1;
  const IFloat *on_node_data;

  int node_sites[LORENTZs];
  int lclWalls=lcl_sites[prop_direction];

  node_sites[0]=lcl_sites[0];
  node_sites[1]=lcl_sites[1];
  node_sites[2]=lcl_sites[2];
  node_sites[3]=lcl_sites[3];
 
  int Dx1=Dx;
  if(Dx<0) Dx1=0;
  int buf_byte_size=DIRACs*COLORs*COMPLEXs*sizeof(IFloat);
  
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

    //offset if no-node
    int propDataOffset;

    if(isFullProp){
      propDataOffset=COMPLEXs*COLORS*Dx1+SPINORs*siteOffset(on_node_site)+Cy*d_weight_Cy+Dy*d_weight_Dy;
    }else{
      propDataOffset=COMPLEXs*COLORS*Dx1+SPINORs*siteOffset(on_node_site,prop_direction)+Cy*d_weight_Cy/lclWalls+Dy*d_weight_Dy/lclWalls;
    }
    // checked projection to on_node_site and sourceOffset

    on_node_data = (IFloat *)(prop_data_p + propDataOffset);
    
  }


#ifndef PARALLEL
  //VRB.FuncEnd(cname, fname) ;
  return (Vector *)on_node_data;
#endif

  // send to the destination node if the site is off-node
  //------------------------------------------------------------------------ 
  Vector send[DIRACs];
  Vector *recv;

  if (on_node) {
    //  VRB.FuncEnd(cname, fname) ;
    return (Vector *)on_node_data;
  } else {
    moveMem((IFloat *)send, (const IFloat *)on_node_data,buf_byte_size);
    recv = &v_tmp2[0] ;
    for (int i = 0; i < 4; ++i) {
      while (site[i] != on_node_site[i]) {
        if (site[i] < 0) {
          getMinusData((IFloat *)recv, (IFloat *)&send[0], buf_byte_size/sizeof(IFloat), i);
          on_node_site[i] -= node_sites[i];
        } else {
          getPlusData((IFloat *)recv, (IFloat *)&send[0], buf_byte_size/sizeof(IFloat), i);
          on_node_site[i] += node_sites[i];
        }
        moveMem((IFloat *)send, (const IFloat *)recv,buf_byte_size);
      }
    }
    //  VRB.FuncEnd(cname, fname) ;
    return recv ;
  }
}


CPS_END_NAMESPACE
