#include<config.h>
CPS_START_NAMESPACE
/*
 *class WspectField
 */

CPS_END_NAMESPACE
#include <alg/w_all.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <stdlib.h>
CPS_START_NAMESPACE

#ifdef PARALLEL
CPS_END_NAMESPACE
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE
#endif

CPS_END_NAMESPACE
#include <util/error.h>                // ERR
#include <util/verbose.h>              // VRB
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <util/vector.h>               // dotProduct
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <comms/glb.h>                   // glb_sum
#include <alg/common_arg.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <alg/alg_w_spect.h>          // AlgWspect::GetCounter()
#include <comms/cbuf.h> 
CPS_START_NAMESPACE

#define MATRIX_SIZE 18

char *WspectField::d_class_name="WspectField"; 



//----------------------------------------------------------------
// CTOR
//----------------------------------------------------------------
WspectField::WspectField(const Lattice &lat, const WspectFuzzing* fuzzing_p)
  :d_lat(lat),fuzz_p(fuzzing_p){

  //buffer for fields
  ft_size=lcl_sites[0]*lcl_sites[1]*lcl_sites[2]*lcl_sites[3]*NUM_FLDS*MATRIX_SIZE;

  ft_p = (Matrix *) smalloc(ft_size * sizeof(Float));
  if (!ft_p) 
    ERR.Pointer(d_class_name, ctor_str, "ft_p");
  VRB.Smalloc(d_class_name, ctor_str, 
	      "ft_p", ft_p, ft_size * sizeof(Float));

  {
    IFloat *fp=(IFloat *)ft_p;
    for(int i=0;i<ft_size;i++) *(fp+i)=0;
  }
}

WspectField::~WspectField(){
  if(ft_p){
    VRB.Sfree(d_class_name, dtor_str, empty_str, ft_p);
    sfree(ft_p);
  }
  
}

void WspectField::run(){
  calcField();
}

//--------------------
// calcField()
// calculate B and E fields from links(original or fuzzed)
//--------------------
void WspectField::calcField(){
  int site[4];
  int siteMin[4];
  int siteMax[4];
  int siteoffset; //site offset in field storage buffer
  Matrix fld;

  for(int i=0;i<4;i++){
    siteMin[i]=0;
    siteMax[i]=lcl_sites[i]-1;
  }
  
  for(site[3]=siteMin[3];site[3]<=siteMax[3];site[3]++){
    for(site[2]=siteMin[2];site[2]<=siteMax[2];site[2]++){
      for(site[1]=siteMin[1];site[1]<=siteMax[1];site[1]++){
	for(site[0]=siteMin[0];site[0]<=siteMax[0];site[0]++){
	  
	  siteoffset=siteOffset(site,-1)*NUM_FLDS; 
	  //B1=F23
	  fld=calcLclField(site,2,3);
	  moveMem((IFloat *)(ft_p+siteoffset+FB1),(const IFloat *)&fld,sizeof(Matrix));
	  //B2=F31
	  calcLclField(site,3,1);
	  moveMem((IFloat *)(ft_p+siteoffset+FB2),(const IFloat *)&fld,sizeof(Matrix));
	  //B3=F12
	  calcLclField(site,1,2);
	  moveMem((IFloat *)(ft_p+siteoffset+FB3),(const IFloat *)&fld,sizeof(Matrix));
	  //E1=F01
	  calcLclField(site,0,1);
	  moveMem((IFloat *)(ft_p+siteoffset+FE1),(const IFloat *)&fld,sizeof(Matrix));
	  //E2=F02
	  calcLclField(site,0,2);
	  moveMem((IFloat *)(ft_p+siteoffset+FE2),(const IFloat *)&fld,sizeof(Matrix));
	  //E3=F03
	  calcLclField(site,0,3);
	  moveMem((IFloat *)(ft_p+siteoffset+FE3),(const IFloat *)&fld,sizeof(Matrix));
	}
      }
    }
  }
}

//-----------------------------
//Matrix calcLclField(int site[], int mu, int nu, WspectFuzzing *fuzz_p){
//-----------------------------
Matrix WspectField::calcLclField(int site[], int mu, int nu){
  Matrix fld;
  fld=calcPlaq(site,1,mu,1,nu);
  fld-=calcPlaq(site,1,nu,1,mu);
  fld+=calcPlaq(site,1,nu,-1,mu);
  fld-=calcPlaq(site,-1,mu,1,nu);
  fld+=calcPlaq(site,-1,mu,-1,nu);
  fld-=calcPlaq(site,-1,nu,-1,mu);
  fld+=calcPlaq(site,-1,nu,1,mu);
  fld-=calcPlaq(site,1,mu,-1,nu);
  return fld;
}

//-----------------------------
//Matrix getField(int site[],FieldTensorId ftId);
//-----------------------------
Matrix WspectField::getField(int site[],FieldTensorId ftId){
  Matrix fld;
  moveMem((IFloat *)&fld,(const IFloat *)(ft_p+siteOffset(site,-1)*NUM_FLDS+ftId),sizeof(Matrix));
  return fld;
}


//-------------------
// Matrix calcPlaq(x,signu,u,signv,v)
//-------------------
/*
 *if fuzz_p!=0 use fuzzed links
 *   else      use original links
 *
 *P_uv=U_u(x)U_v(x+u)U_u(x+v)~U_v(x)~
 *    =U_u(x)U_v(x+u)U_{-u}(x+u+v)U_{-v){x+v);
 *Note: numsign and nusign indicate positive or negative direction
 */
Matrix WspectField::calcPlaq(int x[], int musign, int mu, int nusign, int nu){
  int site[4];
  int i;
  Matrix plaq;
  Matrix tmp1;
  Matrix tmp2;
  
  for(i=0;i<4;i++) site[i]=x[i];
  //tmp1=U_u(x)
  tmp1=getLink(site,musign,mu);
  //tmp2=U_v(x+u)
  site[mu]+=musign;
  tmp2=getLink(site,nusign,nu);
  //tmp1=plaq=U_u(x)U_v(x+u)
  mDotMEqual((IFloat *)&plaq, (const IFloat *)&tmp1, (const IFloat *)&tmp2);
  moveMem((IFloat *)&tmp1,(const IFloat*)&plaq,sizeof(Matrix));
  //tmp1=plaq=U_u(x)U_v(x+u)U_{-u}(x+u+v)
  site[nu]+=nusign;
  tmp2=getLink(site,-musign,mu);
  mDotMEqual((IFloat *)&plaq, (const IFloat *)&tmp1, (const IFloat *)&tmp2);
  moveMem((IFloat *)&tmp1,(const IFloat*)&plaq,sizeof(Matrix));
  //last one
  site[mu]-=musign;
  tmp2=getLink(site,-nusign,nu);
  mDotMEqual((IFloat *)&plaq, (const IFloat *)&tmp1, (const IFloat *)&tmp2);

  return plaq;
}

//-------------------
// Matrix getLink(int x[], int musign, int mu, WspectFuzzing *fuzz_p)
//-------------------
//wrapper function which hides details of getting links
//use Lattice::GetLink
//    or WspectFuzzing::GetLink
//Also add capability of dealing with negative mu value!
//U_{-u}(x)=dag(U_u(x-u))

Matrix WspectField::getLink(int x[],int musign,int mu){
  Matrix link;
  Matrix dag;
  int i;
  int site[4];
  for(i=0;i<4;i++) site[i]=x[i];

  if(musign>0){
    if(fuzz_p==0) return *(d_lat.GetLink(x,mu));
    return *(fuzz_p->GetLink(x,mu));
  }else{
    site[mu]-=1;
    if(fuzz_p==0){
      link=*(d_lat.GetLink(site,mu));
    }else{
      link=*(fuzz_p->GetLink(site,mu));
    }
    dag.Dagger((IFloat *)(&link));
    return dag;
  }
}



CPS_END_NAMESPACE
