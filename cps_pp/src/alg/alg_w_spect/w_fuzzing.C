#include<config.h>
CPS_START_NAMESPACE
/*! \file

  $Id :$
*/
  
/*w_fuzzing.C
 *class WspectFuzzing:public WspectGinfo
 *
 *Construct fuzzed links (APE smearing)
 *
 *Usage: after creating WspectFuzzing object, call run() method
 *       to construct fuzzed links, then use GetLink(site,dir) 
 *       to get fuzzed links
 *
 *
 */

CPS_END_NAMESPACE
#include <alg/w_all.h>

#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif

#include <util/error.h>                // ERR
#include <util/verbose.h>              // VRB
#include <util/vector.h>               // dotProduct
#include <alg/alg_w_spect.h>          // AlgWspect::GetCounter()

#include <comms/glb.h>               // glb_sum(...)
#include <comms/scu.h>              //getMinusData, getPlusData
#include <comms/cbuf.h>              

#include <math.h>
CPS_START_NAMESPACE



//#define DEBUG_FUZZING

#define DIM 4 //fuzzedlink_p's dirac index dimesion
enum { MATRIX_SIZE=18};

//static buffer for offnode communication in GetLink();
static Matrix scuBufM; 


//------------------------------------------------------------
//class WspectFuzzing definition
//------------------------------------------------------------
char *WspectFuzzing::d_class_name="WspectFuzzing";

WspectFuzzing::WspectFuzzing(Lattice &lat, WspectArg & warg, int singleslice,int fuzzing_c_index)
  :d_lat(lat),prop_dir(warg.prop_dir),sliceonly(singleslice){
  //set up local reference to gauge field, prop_dir,...
  //allocate memory to store fuzzed links
  if(singleslice){
    d_size=lcl_sites[0]*lcl_sites[1]*lcl_sites[2]*lcl_sites[3]/lcl_sites[prop_dir]*DIM*MATRIX_SIZE;  //lclVolume*3*3*3*2 IFloats
    slice=-1;
  }else{
    d_size=lcl_sites[0]*lcl_sites[1]*lcl_sites[2]*lcl_sites[3]*DIM*MATRIX_SIZE;
  }
  slice=-1;
  
  fuzzedlink_p = (Matrix *)smalloc(d_size*sizeof(Float));
  if (!fuzzedlink_p)
    ERR.Pointer(d_class_name, ctor_str, empty_str);
  VRB.Smalloc(d_class_name, ctor_str, empty_str, fuzzedlink_p, d_size*sizeof(Float));
  
  
  
#ifdef DEBUG_FUZZING
   printf("Fuzzedlink_p allocated at %x(size=%x)\n",fuzzedlink_p,d_size);
#endif

   //zero out
   IFloat *fp=(IFloat *) fuzzedlink_p;
   for(int i=0;i<d_size;i++) *(fp+i)=0.0;
			       
   //fuzzing algorithm parameters
   fuzzing_on=warg.fuzzing_on;
   fuzzing_level=warg.fuzzing_level;
   if(warg.fuzzing_c_num<=fuzzing_c_index) ERR.General(d_class_name,ctor_str,inconsistent_str);
   fuzzing_c=warg.fuzzing_c[fuzzing_c_index];
   fuzzing_hits=warg.fuzzing_hits;
   printf("Fuzzing constructor: level=%d, c=%e, hits=%d\n", fuzzing_level,fuzzing_c,fuzzing_hits);
   
}

WspectFuzzing::~WspectFuzzing(){
  
  VRB.Sfree(d_class_name, dtor_str, "fuzzedlink_p", fuzzedlink_p);
  sfree(fuzzedlink_p);
}

void WspectFuzzing::run(int wall){
  //construct fuzzed links on single time slice
  char *fname="run";
  //if already done, return immediately
  if(slice==wall) return;

  slice=wall;

  fuzzedlink_tmp_p = (Matrix *)smalloc(d_size*sizeof(Float));
  if (!fuzzedlink_tmp_p)
    ERR.Pointer(d_class_name, ctor_str, empty_str);
   VRB.Smalloc(d_class_name, ctor_str, empty_str, fuzzedlink_tmp_p, d_size*sizeof(Float));
   IFloat *fp=(IFloat *) fuzzedlink_tmp_p;
   for(int i=0;i<d_size;i++) *(fp+i)=0.0;
   
   if(!fuzzedlink_p) ERR.Pointer(d_class_name, fname, empty_str);
   
   const Matrix *gf_p=d_lat.GaugeField();
   
   if(fuzzing_on==0) fuzzing_level=0;
   
   printf("Fuzzing::run(), c=%f, level=%d, hits=%d\n",fuzzing_c, fuzzing_level, fuzzing_hits);
   
#ifdef DEBUG_FUZZING
     printf("Call fuzz() wall=%d, level=%d\n",wall,fuzzing_level);
#endif
    fuzz(fuzzedlink_p, gf_p, wall,fuzzing_c,fuzzing_level,fuzzing_hits);
  
#ifdef DEBUG_FUZZING
  //  display();
#endif

  //deallocate memory of fuzzedlink_tmp_p
  VRB.Sfree(d_class_name, fname, "fuzzedlink_tmp_p", fuzzedlink_tmp_p);
  sfree(fuzzedlink_tmp_p);

}

void WspectFuzzing::run(){
  char *fname="run";
  if(sliceonly) ERR.General(d_class_name, fname, "sliceonly mode");
  
  fuzzedlink_tmp_p = (Matrix *)smalloc(d_size*sizeof(Float));
  if (!fuzzedlink_tmp_p)
    ERR.Pointer(d_class_name, ctor_str, empty_str);
  VRB.Smalloc(d_class_name, ctor_str, empty_str, fuzzedlink_tmp_p, d_size*sizeof(Float));
  IFloat *fp=(IFloat *) fuzzedlink_tmp_p;
  for(int i=0;i<d_size;i++) *(fp+i)=0.0;
   
  int lclWalls=lcl_sites[prop_dir];
  
  if(!fuzzedlink_p) ERR.Pointer(d_class_name, fname, empty_str);
  
  const Matrix *gf_p=d_lat.GaugeField();
  
  if(fuzzing_on==0) fuzzing_level=0;

  printf("Fuzzing::run(), c=%f, level=%d, hits=%d\n",fuzzing_c, fuzzing_level, fuzzing_hits);
  for(int wall=0;wall<lclWalls;wall++){
#ifdef DEBUG_FUZZING
    printf("Call fuzz() wall=%d, level=%d\n",wall,fuzzing_level);
#endif
    fuzz(fuzzedlink_p, gf_p, wall,fuzzing_c,fuzzing_level,fuzzing_hits);
  }
#ifdef DEBUG_FUZZING
  //  display();
#endif

  //deallocate memory of fuzzedlink_tmp_p
  VRB.Sfree(d_class_name, fname, "fuzzedlink_tmp_p", fuzzedlink_tmp_p);
  sfree(fuzzedlink_tmp_p);
}




//----------------------------------------------------------------------------------
//fuzz(Matrix *FU, const Matrix *U, int lclWall, Float c, int f_iter, int hit_iter)
//
//Generic fuzzing fuction
//do spatial links fuzzing on single time slice
// 
//----------------------------------------------------------------------------------
//U[site][u][c1][c2], site[4],u=0,1,2,3
//FU[site][i][c1][c2], site[4],i=0,1,2

void  WspectFuzzing::fuzz(Matrix *FU, const Matrix *U, int lclWall, Float c, int flevel, int hit_iter){
  char *fname="fuzz";
  //printf("Inside WpsectFuzzing::fuzz\n");
  //assume memory has been allocated
  if(!FU || !U) ERR.General(d_class_name,fname,"Null pointer");

  int site[4], siteMin[4], siteMax[4];
  int i,dir,dir1;
  /*Step1:
   *FU(x,t,dir)=c*U(x,t,dir)+Sum of Staples(spatial staples only)
   *
   *The sum of staples is similar to Lattice::staple
   */

  //set site boundary
  for(i=0;i<4;i++){
    if(i==prop_dir) {
      siteMin[i]=siteMax[i]=lclWall;
    }else{
      siteMin[i]=0;
      siteMax[i]=lcl_sites[i]-1;
    }
  }

  // printf("Start copying\n");
  //copy U to FU
  for(site[3]=siteMin[3];site[3]<=siteMax[3];site[3]++){
    for(site[2]=siteMin[2];site[2]<=siteMax[2];site[2]++){
      for(site[1]=siteMin[1];site[1]<=siteMax[1];site[1]++){
	for(site[0]=siteMin[0];site[0]<=siteMax[0];site[0]++){
	  for(dir=0;dir<=3;dir++){
	    //printf("dir=%d\n",dir);
	    if(sliceonly){
	      *(FU+siteOffset(site,prop_dir)*DIM+dir) = *(U+siteOffset(site)*4+dir);
	    }else{
	      *(FU+siteOffset(site)*DIM+dir) = *(U+siteOffset(site)*4+dir);
	    }
	    //printf("After copying.........\n");
	    //displaySU3(FU+siteOffset(site)*4+dir);
	    //displaySU3((Matrix *)U+siteOffset(site)*4+dir);
	  }//dir
	}
      }
    }
  }
  
  //if no fuzzing, done 
  if(fuzzing_on==0) return;

  //do fuzzing
  //temporary buffers(time slice)
#ifdef DEBUG_FUZZING
  //printf("WspectFuzzing::fuzz() with fuzzing_on==1\n");
#endif

 
  Matrix FUZZM;
  Matrix FUZZOld,dagM,multiM;
  Matrix *FUM_p; //pointer to fuzzedlink_p
  Matrix *FUTMPM_p; //pointer to fuzzedlink_tmp_p buffer
  Matrix staplev; //sum over spatial staples
  int siteoffset;

  for(int iter=0;iter<flevel;iter++){
#ifdef DEBUG_FUZZING
    printf("Fuzzing iteration #%d\n",iter);
#endif
    for(dir=0;dir<=2;dir++){
      //convert dir to one of the three directions!=prop_direction
      dir1=(dir+1+prop_dir)%4;
#ifdef DEBUG_FUZZING
      printf("Dir #%d\n",dir1);
#endif
      //for all sites on the wall
      for(site[3]=siteMin[3];site[3]<=siteMax[3];site[3]++){
	for(site[2]=siteMin[2];site[2]<=siteMax[2];site[2]++){
	  for(site[1]=siteMin[1];site[1]<=siteMax[1];site[1]++){
	    for(site[0]=siteMin[0];site[0]<=siteMax[0];site[0]++){
	      //calculate stape(spatial), NOT implemented yet
	      staplev=staple_s((Matrix *)FU,site,dir1);
	      //FUZZ=c*FU+staple
	      if(sliceonly) {
		siteoffset=siteOffset(site,prop_dir);
	      }else{
		siteoffset=siteOffset(site);
	      }
	      FUM_p=FU+siteoffset*DIM+dir1;
	      FUZZM=*FUM_p;
	      (FUZZM)*=c;
	      FUZZM+=staplev;
	      
	      FUZZOld=FUZZM;
	      
	      FUZZM.Unitarize();
#ifdef DEBUG_FUZZING
	      printf("---------After projSU3: check=%f\n",FUZZM.ErrorSU3());
	      //displaySU3(&FUZZM);
#endif
	      //multiM=Dag(FUZZOLD)*FUZZM;
	      dagM.Dagger((IFloat *)&FUZZOld);
	      mDotMEqual((IFloat *)&multiM,(const IFloat *)&dagM, (const IFloat *)&FUZZM);
	    	      
	      //cool
	      for(int hit=1;hit<=hit_iter;hit++){
		cabbibo(hit,hit_iter,multiM,FUZZM);
	      }//hit
	      //check SU3 of fuzzed link
#ifdef DEBUG_FUZZING
	      printf("---------After cooling: SU3 check=%f\n",FUZZM.ErrorSU3());
	      // displaySU3(&FUZZM);
#endif
	      //copy FUZZ to temp buffer
	      FUTMPM_p=fuzzedlink_tmp_p+siteoffset*DIM+dir1;
	      *(FUTMPM_p)=FUZZM;
	    }
	  }
	} 
      }//site[0]
      
    }//dir loop
    //copy fuzzedlink_tmp_p to fuzzedlink_p
    for(dir=0;dir<=2;dir++){
      //convert dir to one of the three directions!=prop_direction
      dir1=(dir+1+prop_dir)%4;
#ifdef DEBUG_FUZZING
      printf("Dir #%d\n",dir1);
#endif
      //for all sites on the wall
      for(site[3]=siteMin[3];site[3]<=siteMax[3];site[3]++){
	for(site[2]=siteMin[2];site[2]<=siteMax[2];site[2]++){
	  for(site[1]=siteMin[1];site[1]<=siteMax[1];site[1]++){
	    for(site[0]=siteMin[0];site[0]<=siteMax[0];site[0]++){
	      if(sliceonly) {
		siteoffset=siteOffset(site,prop_dir);
	      }else{
		siteoffset=siteOffset(site);
	      }
	      *(FU+siteoffset*DIM+dir1)=*(fuzzedlink_tmp_p+siteoffset*DIM+dir1);
	    }
	  }
	}
      }
    }

  }//fuzzing level loop

}


//----------------------------------------------------------
//
//----------------------------------------------------------
//Actually the GetLinkOld is also OK

const Matrix *
WspectFuzzing::GetLink(const int *site, int dir) const
{
  char *fname = "GetLink";
  
  //VRB.Debug(cname, fname, "link %3i %3i %3i %3i ; %i\n",
  //  site[0], site[1], site[2], site[3], dir) ;
#ifdef DEBUG_FUZZING
  //printf("Getlink %3i %3i %3i %3i ; %i\n",
  //	    site[0], site[1], site[2], site[3], dir) ;
#endif

  // offset out-of-range coordinates site[] into on_node_site[]
  // in order to locate the link
  //------------------------------------------------------------------------

  if(dir==prop_dir) ERR.General(d_class_name,fname,"Invalid dir");
  //check if the buffer is holding the right slice
  //must call run(wall) right before using fuzzed links on this slice!!
  if(sliceonly && site[prop_dir]!=slice) ERR.General(d_class_name,fname,"Invalid slice");
  
  int on_node_site[4];
  int on_node = 1;
  
 
  const Matrix *on_node_link;
  {
    for (int i = 0; i < 4; ++i) {
      on_node_site[i] = site[i] ;
      while (on_node_site[i] < 0) {
        on_node_site[i] += lcl_sites[i] ;
      }
      on_node_site[i] %= lcl_sites[i];
      if (on_node_site[i] != site[i]) {
        on_node = 0;
      }
    }
    if(sliceonly){
      on_node_link =fuzzedlink_p + siteOffset(on_node_site,prop_dir)*DIM + dir ;
    }else{ 
      on_node_link =fuzzedlink_p + siteOffset(on_node_site)*DIM + dir ;
    }
    
  }
  
#ifndef PARALLEL
  //VRB.FuncEnd(cname, fname) ;
  return on_node_link;
#endif
  
  // send to the destination node if the site is off-node
  //------------------------------------------------------------------------
  if (on_node) {
    //  VRB.FuncEnd(cname, fname) ;
    //printf("On node data\n");
    return on_node_link;
  } else {
    Matrix send = *on_node_link;
    Matrix &recv = scuBufM ;
    for (int i = 0; i < 4; ++i) {
      while (site[i] != on_node_site[i]) {
        if (site[i] < 0) {
          getMinusData((IFloat *)&recv, (IFloat *)&send, sizeof(recv)/sizeof(IFloat), i);
          on_node_site[i] -= lcl_sites[i];
        } else {
          getPlusData((IFloat *)&recv, (IFloat *)&send, sizeof(recv)/sizeof(IFloat), i);
          on_node_site[i] += lcl_sites[i];
        }
        send = recv;
      }
    }
    //  VRB.FuncEnd(cname, fname) ;
    return &recv ;
  }
}

//------------------------------------------------------------------
// Matrix staple_s(Complex *gf_p,int site[4],int mu)
// It calculates the staple field(spatial only) at x, mu.
// Note: input gf_p is fuzzed gauge field (spatial links only)
// The staple field is:
//
//      V_u(x) = \sum_v(!=u !=T) {
//              U_v(x) U_u(x+v) U_v(x+u)~ 
//           +  U_v(x-v)~ U_u(x-v) U_v(x+u-v) }
//
//------------------------------------------------------------------
//Used to turn on/off Circular buffer access
//#define USE_CBUF
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



Matrix WspectFuzzing::staple_s(Matrix *gf_p,int x[4],int mu){
  
  //char *fname = "Staple";
 
  // set cbuf registers
#ifdef USE_CBUF
  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(4, CBUF_MODE4);
#endif

  int site[4];
  site[0]=x[0]; site[1]=x[1]; site[2]=x[2]; site[3]=x[3];
  
  const Matrix *p1;
  Matrix stap; //result
  stap.ZeroMatrix();

  //temporary buffers
  Matrix mp2,mp3;
  //int offset_x = siteOffset(site)*3;
  //Matrix *g_offset = gf_p+offset_x;

  for(int nu = 0; nu < 4; ++nu) {
    if(nu != mu && nu!=prop_dir) {

      // calculate U_v(x) U_u(x+v) U_v(x+u)~ 

      //----------------------------------------------------------
      // mp3 = U_v(x)
      //----------------------------------------------------------
      //should use WspectFuzzing::getLink()!!!
      p1 = GetLink(site, nu);
      moveMem((IFloat *)&mp3, (IFloat *)p1+MYBANK4_BASE+MYBANK_SIZE,
              MATRIX_SIZE * sizeof(IFloat));
      
      
      
      //----------------------------------------------------------
      // p1 = &U_u(x+v)
      //----------------------------------------------------------
      site[nu]++;
      p1 = GetLink(site, mu);
      site[nu]--;

      //----------------------------------------------------------
      // mp2 = U_v(x) U_u(x+v)
      //----------------------------------------------------------
      mDotMEqual((IFloat *)&mp2, (const IFloat *)&mp3, 
		 (const IFloat *)p1+MYBANK2_BASE);
                 


      //----------------------------------------------------------
      //  mp3 = U_v(x+u)~
      //----------------------------------------------------------
      site[mu]++;
      mp3.Dagger((IFloat *)(GetLink(site,nu))+MYBANK4_BASE);
      site[mu]--;

      //----------------------------------------------------------
      // calculate  U_v(x) U_u(x+v) U_v(x+u)~= mp2 * mp3
      //----------------------------------------------------------
    
      mDotMPlus((IFloat *)&stap, (const IFloat *)&mp2,
                  (const IFloat *)&mp3);




      //----------------------------------------------------------
      //  calculate U_v(x-v)~ U_u(x-v) U_v(x+u-v)
      //----------------------------------------------------------
      //----------------------------------------------------------
      // p1 = U_v(x-v)
      // mp3 = U_v(x-v)~
      //---------------------------------------------------------
      site[nu]--;
      p1 = GetLink(site, nu);
      site[nu]++;
      mp3.Dagger((IFloat *)p1+MYBANK2_BASE);
      
      
      //----------------------------------------------------------
      // mp2 = U_v(x-v)~ U_u(x-v) 
      //----------------------------------------------------------
      site[nu]--;
      mDotMEqual((IFloat *)&mp2, (const IFloat *)&mp3,  (const IFloat *)(GetLink(site,mu))+MYBANK4_BASE);
      site[nu]++;
      
      
      //----------------------------------------------------------
      // mp3 = U_v(x+u-v)
      //----------------------------------------------------------
      site[nu]--;
      site[mu]++;
      moveMem((IFloat *)&mp3, (const IFloat *)(GetLink(site,nu))+MYBANK2_BASE,
	      MATRIX_SIZE * sizeof(IFloat));
      site[mu]--;
      site[nu]++;

      
      //----------------------------------------------------------
      // stap += mp2 * mp3
      //----------------------------------------------------------
      mDotMPlus((IFloat *)&stap, (const IFloat *)&mp2,
		(const IFloat *)&mp3);
      // dummy read
      mp2 = *((IFloat *)&mp3+MYBANK4_BASE);
     
    }
  }
  return stap;
}


//---------------------------------------------------------------
//void cabbibo(int hit, int maxhit, Matrix X, Matrix &T);
//---------------------------------------------------------------
void WspectFuzzing::cabbibo(int hit, int maxhit, Matrix &X, Matrix &T){
  /*
   *  + + o
   *  + + o
   *  o o o
   */
  
  Complex h1,h2,h;
  Matrix SU2,S;

  h1=(X(0,0)+conj(X(1,1)))*0.5;
  h2=(X(1,0)-conj(X(0,1)))*0.5;
  h=sqrt(norm(h1)+norm(h2));
  h=1./h;
  
  h1=h1*h;
  h2=h2*h;
    
  SU2.UnitMatrix();
  SU2(0,0)=conj(h1);
  SU2(0,1)=conj(h2);
  SU2(1,0)=-h2;
  SU2(1,1)= h1;
    
  S.DotMEqual(T,SU2);
    
  T.DotMEqual(X,SU2);
  
  X=T;
  T=S;

    // o o o
    // o + +
    // o + +
    
    h1=(X(1,1)+conj(X(2,2)))*0.5;
    h2=(X(2,1)-conj(X(1,2)))*0.5;
    h=sqrt(norm(h1)+norm(h2));
    h=1./h;
    
    h1=h1*h;
    h2=h2*h;
    //--------------START
    SU2.UnitMatrix();
    SU2(1,1)=conj(h1);
    SU2(1,2)=conj(h2);
    SU2(2,1)=-h2;
    SU2(2,2)= h1;
    
    S.DotMEqual(T,SU2);

    T.DotMEqual(X,SU2);

    X=T;
    T=S;

    // + o +
    // o o o
    // + o +
    h1=(X(0,0)+conj(X(2,2)))*0.5;
    h2=(X(2,0)-conj(X(0,2)))*0.5;
    h=sqrt(norm(h1)+norm(h2));
    h=1./h;
    
    h1=h1*h;
    h2=h2*h;
    
    SU2.UnitMatrix();
    SU2(0,0)=conj(h1);
    SU2(0,2)=conj(h2);
    SU2(2,0)=-h2;
    SU2(2,2)= h1 ;             
  
    S.DotMEqual(T,SU2);

    if(hit<maxhit){
    //    return S=T*SU2 as new trial T
    //    return T=X*SU2 as new X
      T.DotMEqual(X,SU2);
      
      X=T;
      T=S;
  }//if(ht<maxhit)

  if(hit==maxhit){
    //return S (solution) as final matrix as T
    T=S;
  }

}

#ifdef DEBUG_FUZZING
/*
void WspectFuzzing::display(){

  //dislay fuzzedlinks and original links
  int site[4], siteMin[4], siteMax[4];
  int i,dir;
  //Step1:
  //FU(x,t,dir)=c*U(x,t,dir)+Sum of Staples(spatial staples only)
  //
  //The sum of staples is similar to Lattice::staple
  //
  printf("=========Links on slice #%d========\n",slice);

  File *fp=Fopen("fuzzedlink.dat","a");
  

  //set site boundary
  for(i=0;i<4;i++){
    siteMin[i]=0;
    siteMax[i]=lcl_sites[i]-1;
    if(i==prop_dir) siteMin[i]=siteMax[i]=slice;
  }

  int prop_dir1=prop_dir;
  if(!sliceonly) prop_dir1=-1;

  Matrix *FU=fuzzedlink_p;
  Matrix *U=d_lat.GaugeField();
  //copy U to FU
  for(site[3]=siteMin[3];site[3]<=siteMax[3];site[3]++){
    for(site[2]=siteMin[2];site[2]<=siteMax[2];site[2]++){
      for(site[1]=siteMin[1];site[1]<=siteMax[1];site[1]++){
	for(site[0]=siteMin[0];site[0]<=siteMax[0];site[0]++){
	  for(dir=0;dir<=2;dir++){
	    //dir1=(dir+1+prop_dir)%4;//directions different from prop_dir
	    //printf("dir=%d\n",dir);
	    printf("SU3 check, fuzzed=%f, org=%f\n", (FU+siteOffset(site,prop_dir1)*DIM+dir)->ErrorSU3(), (U+siteOffset(site)*4+dir)->ErrorSU3());
	    displaySU3(FU+siteOffset(site,prop_dir1)*DIM+dir);
	    displaySU3(U+siteOffset(site)*4+dir);
	    printf("Getlink check, fuzzed=%f, org=%f\n",(*GetLink(site,dir)).ErrorSU3(), (*d_lat.GetLink(site,dir)).ErrorSU3());
	    if((*GetLink(site,dir)).ErrorSU3()>0.0001) ERR.General(d_class_name,"fuzz","SU3 error"); 

	    //write to file
	    int icol;
	    Matrix *matrix=GetLink(site,dir);
	    for (icol = 0; icol  < COLORs; icol++)
	      Fprintf(fp,"(%14.6e, %14.6e)", matrix[2*icol],matrix[2*icol+1]);
	    Fprintf(fp,"\n");
	    for (icol = COLORs; icol  < 2*COLORs; icol++)
	      Fprintf(fp"(%14.6e, %14.6e)", matrix[2*icol],matrix[2*icol+1]);
	    Fprintf(fp,"\n");
	    for (icol = 2*COLORs; icol  < 3*COLORs; icol++)
	      Fprintf(fp,"(%14.6e, %14.6e)", matrix[2*icol],matrix[2*icol+1]);
	    Fprintf(fp,"\n");
	  }//dir
	}
      }
    }
  }

  Fclose(fp);

}

void WspectFuzzing::displaySU3(Matrix *m){

  Float *matrix=(Float *)m;
   int display=0;
  for (int i = 0; i < COLORs*COLORs*COMPLEXs; i++){
    if (matrix[i] != 0.) display=1;
  }
  int icol;
  Float fcol;
  // display matrix only if some elements are non-zero  
  if (display != 0){
    printf("DISPLAY MATRIX ----------\n");
    
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
    printf("DISPLAY MATRIX: all elements are zero -- \n");
  }

  //SU3  col/row check

  for (icol = 0; icol < 2*COLORs; icol++)
    fcol+=matrix[icol]*matrix[icol];
  printf("first raw     =  1 = %22.12e \n", fcol);
  fcol = 0.;
  for (icol = 2*COLORs; icol < 4*COLORs; icol++)
    fcol+=matrix[icol]*matrix[icol];
  printf("second raw     = 1 = %22.12e \n", fcol);
  fcol = 0.;
  for (icol = 4*COLORs; icol < 6*COLORs; icol++)
    fcol+=matrix[icol]*matrix[icol];
  printf("third raw     =  1 = %22.12e \n", fcol);
  
  fcol = 0.;
  for (icol = 0; icol < COLORs; icol++)
    fcol+=matrix[6*icol]*matrix[6*icol]+matrix[6*icol+1]*matrix[6*icol+1];
  printf("first  column =  1 = %22.12e \n", fcol);
	
  fcol = 0.;
  for (icol = 0; icol < COLORs; icol++)
    fcol+=matrix[6*icol+2]*matrix[6*icol+2]+matrix[6*icol+3]*matrix[6*icol+3];
  printf("second column =  1 = %22.12e \n", fcol);

  fcol = 0.;
  for (icol = 0; icol < COLORs; icol++)
    fcol+=matrix[6*icol+4]*matrix[6*icol+4]+matrix[6*icol+5]*matrix[6*icol+5];
  printf("third column  =  1 = %22.12e \n", fcol);

}
*/
#endif //ifdef DEBUG_FUZZING





CPS_END_NAMESPACE
