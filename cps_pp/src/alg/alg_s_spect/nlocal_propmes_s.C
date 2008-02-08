//**********************************************************************
//Non-local Staggered Meson Propagator
//
//Originally written by chateau, annotated by michaelc
//
//This class calculates the N*N-1 = 15 pi mesons associated with
//a 4-flavor staggered fermions.  Because of spin-flavor mixing, these
//pions are not degenerate, with the degneracy broken by O(a^2) terms
//in the staggered action.  13 of the 15 pions involve quark propagators
//that are shifted within a unit cube on a constant time slice
//
//Note, to calculate these nonlocal pions, one must calculate eight
//different quark propogators, 
//
//For more information on staggered meson operators see:
//
//Golterman, Nuc. Phys. B273(1986) 663-676
//Golterman and Smit, Nuc. Phys. B245(1984) 61-88
//Ishizuka et al., Nuc. Phys. B411(1994) 875-901
//**********************************************************************

#include<config.h>
#include <alg/nlocal_propmes_s.h>
#include <alg/quark_prop_s.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <comms/scu.h>	
#include <alg/alg_s_spect.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <alg/enum.h>

#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
CPS_START_NAMESPACE
//enum { VECT_LEN = 6, MATRIX_SIZE = 18 };


char NLSMesonPropS::cname[] = "NLSMesonPropS";
//-----------------------------------------------------------------------
// local to this file
//-----------------------------------------------------------------------
static int nx[4];

void NLSMesonPropS::getNeighbors(int t)
{
  for(int src_i=0; src_i<8 ; src_i++)
    transfer(qp0[src_i],buffer0[src_i],t);
}

void NLSMesonPropS::transfer(Float **qp, Float **buffer, int t)
{

  //------------------------------------------------------------
  // * 26 kinds off node data transfer
  //------------------------------------------------------------

  int x[4],color,mu,mu_1,mu_2;
  Float *rcv, *snd;

  // first transfer +x,+y,+z
  // get data from (0,0,-1) (0,-1,0) (-1,0,0)
  x[3] = t;
  for(mu = 0 ; mu < 3 ; mu++)  {
    mu_1 = (mu + 1) % 3;
    mu_2 = (mu + 2) % 3;
    for(x[mu_1] = 0 ; x[mu_1] < nx[mu_1] ; x[mu_1]++) {
      for(x[mu_2] = 0 ; x[mu_2] < nx[mu_2] ; x[mu_2]++) {  
        for(color = 0 ; color < 3 ; color++)  {
          x[mu] = -1;
          rcv = buffer[color] + x_map(x);
          x[mu] = nx[mu]-1;
          snd = qp[color] + X_OFFSET(x);
          getMinusData((IFloat *)rcv, (IFloat *)snd, VECT_LEN, mu );
//          printf(" %f %f \n\n",*rcv,*snd);  it was here before the error happens
          x[mu] = nx[mu];
          rcv = buffer[color] + x_map(x);
          x[mu] = 0; 
          snd = qp[color] + X_OFFSET(x);
          getPlusData((IFloat *)rcv, (IFloat *)snd, VECT_LEN, mu );
//          printf(" %f %f \n\n",*rcv,*snd);
        }
      }
    } 
  }

  // second transfer x,y,z
  // get data
  for(mu = 0 ; mu < 3 ; mu++)  {
    mu_1 = (mu + 1) % 3;
    mu_2 = (mu + 2) % 3;
      for(x[mu_2] = 0 ; x[mu_2] < nx[mu_2] ; x[mu_2]++)  
        for(color = 0 ; color < 3 ; color++)  {
          x[mu] = -1;
          x[mu_1] = -1;
          rcv = buffer[color] + x_map(x);
          x[mu] = nx[mu] - 1;
          snd = buffer[color] + x_map(x);
          getMinusData((IFloat *)rcv, (IFloat *)snd, VECT_LEN, mu );
//          printf(" %f %f \n\n",*rcv,*snd);

          x[mu] = -1;
          x[mu_1] = nx[mu_1];
          rcv = buffer[color] + x_map(x);
          x[mu] = nx[mu]-1;
          snd = buffer[color] + x_map(x);
          getMinusData((IFloat *)rcv, (IFloat *)snd, VECT_LEN, mu );
//          printf(" %f %f \n\n",*rcv,*snd);

          x[mu] = nx[mu];
          x[mu_1] = -1;
          rcv = buffer[color] + x_map(x);
          x[mu] = 0; 
          snd = buffer[color] + x_map(x);
          getPlusData((IFloat *)rcv, (IFloat *)snd, VECT_LEN, mu );
//          printf(" %f %f \n\n",*rcv,*snd);

          x[mu] = nx[mu];
          x[mu_1] = nx[mu_1];
          rcv = buffer[color] + x_map(x);
          x[mu] = 0; 
          snd = buffer[color] + x_map(x);
          getPlusData((IFloat *)rcv, (IFloat *)snd, VECT_LEN, mu );
//          printf(" %f %f \n\n",*rcv,*snd);
          }
  }

  int pos[4][4]={{-1,-1,-1,t},{-1,-1,nx[2],t},{-1,nx[1],-1,t},{-1,nx[1],nx[2],t}};
// mu is for different purpose.
  for(mu = 0 ; mu < 4 ; mu++ )  
    for(color=0;color<3;color++)  {
    pos[mu][0] = -1;
    rcv = buffer[color] + x_map(pos[mu]);
    pos[mu][0] = nx[0] - 1;
    snd = buffer[color] + x_map(pos[mu]);
    getMinusData((IFloat *)rcv, (IFloat *)snd, VECT_LEN, 0);
//          printf(" %f %f \n\n",*rcv,*snd);
  
    pos[mu][0] = nx[0];
    rcv = buffer[color] + x_map(pos[mu]);
    pos[mu][0] =  0;
    snd = buffer[color] + x_map(pos[mu]);
    getPlusData((IFloat *)rcv, (IFloat *)snd, VECT_LEN, 0);
//          printf(" %f %f \n\n",*rcv,*snd);
  }
}

int NLSMesonPropS::x_map(int *x)
{
  int offset1=0;
  int y[4]={ x[0]+1, x[1]+1, x[2]+1,0};  
  int offset0 = y[0] + ( nx[0]+2 ) * y[1] + ( nx[0]+2 ) * ( nx[1]+2 ) * y[2] ;
  if (x[2]== -1) offset1=0;
  else if (x[2] == nx[2] ) offset1 = nx[0] * nx[1] * nx[2] ;
    else if( x[1] == -1 ) offset1 = nx[0] * nx[1] * x[2] ;
      else if ( x[1] == nx[1] ) offset1 = nx[0] * nx[1] * ( x[2] + 1 );
        else if ( x[0] == -1 ) offset1 = nx[0] * nx[1] * x[2] + nx[0] * x[1]; 
          else if ( x[0] == nx[0] ) offset1 = nx[0] * nx[1] * x[2] + nx[0] * ( x[1] + 1); 
//  printf("%d %d %d %d %d off0=%d 0ff1=%d",IsOffNode(x),x[0],x[1],x[2],x[3],offset0,offset1);
  return (offset0 - offset1)*VECT_LEN;
}

/* int NLSMesonPropS::X_OFFSET(int *x)
// {
//  return x[0] + nx[0] * x[1] + nx[0] * nx[1] * x[2];
// }
*/

int NLSMesonPropS::IsOffNode(int *y)
{
  if ( ( (y[0] + 1) % (nx[0]+1) ) * ( (y[1] + 1) % (nx[1]+1) ) * ( (y[2] + 1) % (nx[2]+1) ) ) return 0;
else return 1;
}    

Complex NLSMesonPropS::traceG1DagG2(Float *G1[], Float *G2[], int offset1,int offset2)
{
  Complex res;
  for(int color = 0; color < 3; ++color) {
    Complex *a = (Complex *) (G1[color]+offset1);
    Complex *b = (Complex *) (G2[color]+offset2);
    for (int i = 0; i < 3; ++i) {
      res += conj(a[i]) * b[i];
    }
  }
  return res;
}


// Here  ( x y z t ) -> ( 0 1 2 3 ) 
// But, in Golterman paper they used ( x y z t ) -> ( 1 2 3 4 )
// Thus, the sign factor convention is different

void NLSMesonPropS::localVal(Complex *currp, int *s)
{
  int x[4],y[4],i,j,k,l,m,B[4],A[4],C[4];
  Complex temp;
  A[3]=B[3]=C[3]=0;
  y[3]=x[3]=s[3];
  int offset1 = X_OFFSET(s);

//---------------------------------------------------------------------
// type I class 1++
// type II class 1+-
//---------------------------------------------------------------------

  for( i = 0 ; i < 8 ; i ++ )  //for each source 
    { 
    B[0]=int( i/4 );
    B[1]=int( (i%4) / 2 );  
    B[2]=i%2;


    currp[0] -= epsilon(s) * epsilon(B) * traceG1DagG2(qp0[i],qp0[i] , offset1,offset1 ); 
    currp[1] -= eta(s,3)*eta(B,3) * epsilon(s)*epsilon(B) * traceG1DagG2(qp0[i] , qp0[i] , offset1, offset1); 

    }
 
//---------------------------------------------------------------------
// type VII class 3"--  3-pions 
// type VIII class 3"-+  3-pions 
//---------------------------------------------------------------------

  for( l = 0 ; l < 3 ; l++ )
    {
    y[0]= s[0] ;
    y[1]= s[1] ;
    y[2]= s[2] ;   // y is the shift of anti-quark in operator
    for( i = 0 ; i < 8 ; i++ )  //for each source 
      { 
      B[0]=int(i/4);
      B[1]=int( (i%4) / 2 );  
      B[2]=i%2;
  
      B[l] = ( B[l] + 1 ) % 2;    // B+/delta for source
      j = B[0] * 4 + B[1] * 2 + B[2];  // j the shifted source number
      B[l] = ( B[l] + 1 ) % 2;    // modify B for sign factor 

      for( y[l]=s[l]-1 ; y[l]<=s[l]+1 ; y[l]+=2 )
        {
//        if ( IsOffNode(y) )  
//          temp = traceG1DagG2(buffer0[i] , qp0[j] , x_map(y) , offset1);
//          else temp = traceG1DagG2( qp0[i] , qp0[j] , X_OFFSET(y) , offset1);
        if ( IsOffNode(y) )  
          temp = traceG1DagG2(qp0[j] , buffer0[i] , offset1 , x_map(y) );
          else temp = traceG1DagG2( qp0[j] , qp0[i] , offset1 , X_OFFSET(y) );
  
        currp[l+2] -= zeta(s,l)*zeta(B,l)  * temp; 
        currp[l+5] -= eta(s,3)*eta(B,3) * zeta(s,l)*zeta(B,l) * temp; 
        }
      }
    } 

//---------------------------------------------------------------------
// type XIII class 3"++  3-pions 
// type XIV  class 3"+-  3-pions 
// epsilon_klm zeta_k D_k ( zeta_l D_l chi ) 
// epsilon_klm zeta_k zeta_l ( x + k ) ( D_k D_l chi )   
// 2 zeta_0 zeta_1 ( D_k D_l chi )       
// 2 zeta_1 zeta_2 ( D_k D_l chi )       
// 2 zeta_2 zeta_0 ( D_k D_l chi )       
//---------------------------------------------------------------------

  for( m = 0 ; m < 3 ; m++ )
    {

    y[0] = s[0] ;
    y[1] = s[1] ;
    y[2] = s[2] ;   

    k = ( m + 1 ) % 3;
    l = ( m + 2 ) % 3;

    for( i = 0 ; i < 8 ; i++ )  //for each source 
      { 
      B[0]=C[0]=int(i/4);
      B[1]=C[1]=int( (i%4) / 2 );  
      B[2]=C[2]=i%2;
  
      B[k] = ( B[k] + 1 ) % 2;    
      B[l] = ( B[l] + 1 ) % 2;    
      j = B[0] * 4 + B[1] * 2 + B[2];  //source number

      for( y[k]=s[k]-1 ; y[k]<=s[k]+1 ; y[k]+=2 )
		  for( y[l]=s[l]-1 ; y[l]<=s[l]+1 ; y[l]+=2 )
          {
//          if ( IsOffNode(y) )  
//            temp = traceG1DagG2(buffer0[i] , qp0[j] , x_map(y) , offset1);
//            else temp = traceG1DagG2( qp0[i] , qp0[j] , X_OFFSET(y) , offset1);
        if ( IsOffNode(y) )  
          temp = traceG1DagG2(qp0[j] , buffer0[i] , offset1 , x_map(y) );
          else temp = traceG1DagG2( qp0[j] , qp0[i] , offset1 , X_OFFSET(y) );
    
          currp[m+8] +=  2*zeta(s,k)*zeta(C,k) * zeta(s,l)*zeta(C,l) * epsilon(s)*epsilon(C) * temp; 
          currp[m+11] += 2*eta(s,3)*eta(C,3) * zeta(s,k)*zeta(C,k) * zeta(s,l)*zeta(C,l) * epsilon(s)*epsilon(C) * temp; 
          }
        
      }
    }                   
//---------------------------------------------------------------------
// Class 3 type XVIII representation 1--  1-pion 
//---------------------------------------------------------------------
  for( i = 0 ; i < 8 ; i++ )  //for each source 
    { 
    A[0]=B[0]=C[0]=int(i/4);
    A[1]=B[1]=C[1]=int( (i%4) / 2 );  
    A[2]=B[2]=C[2]=i%2;         // C non shift
                         
    B[0] = ( B[0] + 1 ) % 2;    // B+/delta for source
    B[1] = ( B[1] + 1 ) % 2;    // B+/delta for source
    B[2] = ( B[2] + 1 ) % 2;    // B+/delta for source
    j = B[0] * 4 + B[1] * 2 + B[2];  //source number

    B[2] = C[2];                // B double shift 
    A[0] = B[0];                // A sigle shift
    for(y[0] = s[0] -1 ; y[0] <= s[0] + 1 ; y[0]+=2)
      for(y[1] = s[1] -1 ; y[1] <= s[1] + 1 ; y[1]+=2)
        for(y[2] = s[2] -1 ; y[2] <= s[2] + 1 ; y[2]+=2)
          {
//          if ( IsOffNode(y) )  
//            temp = traceG1DagG2(buffer0[i] , qp0[j] , x_map(y) , offset1);
//            else temp = traceG1DagG2(qp0[i] , qp0[j] , X_OFFSET(y) , offset1);

        if ( IsOffNode(y) )  
          {
			 temp = traceG1DagG2(qp0[j] , buffer0[i] , offset1 , x_map(y) );
//          printf("%d,%d,%d,%d,%d\n",y[0],y[1],y[2],y[3],IsOffNode(y) );
          } 
          else temp = traceG1DagG2( qp0[j] , qp0[i] , offset1 , X_OFFSET(y) );
      
          x[0]=s[0]+1; x[1]=s[1]+1; x[2]=s[2];  // double shift
          temp *= eta(x,2) * eta(B,2);

          x[0]=s[0]+1; x[1]=s[1]; x[2]=s[2];   // single shift
          temp *= eta(x,1) * eta(A,1);

          temp *= eta(s,3)*eta(C,3) * epsilon(s)*epsilon(C);
          currp[14] -= temp; 
          } 
    }

  return;
}

/*
int NLSMesonPropS::IsOffNode(int *y)
{
if ( ( (y[0] + 1) % nx[0] ) * ( (y[1] + 1) % nx[1] ) * ( (y[2] + 1) % nx[2] ) ) return 0;
else return 1;
}    
*/

//eta_l(y) = (-1)^(y[0]+y[1]+...y[l-1])
int NLSMesonPropS::eta(int *y, int l)
{
  int tot=0;
//  l--;  // to match the notation
  for(int i=0;i<l;i++) tot+=y[i];
  if (((tot + 100)%2) == 0) return 1;
  else return -1;
}

//zeta_l(y) = (-1)^(y[3]+y[2]+...y[l+1])
int NLSMesonPropS::zeta(int *y, int l)
{
  int tot=0;
//  l--;
  for(int i=3;i>l;i--) tot+=y[i];
  if (((tot + 100)%2) == 0) return 1;
  else return -1;
}

//epsilon(y) = (-1)^(y[0]+y[1]+y[2]+y[3])
int NLSMesonPropS::epsilon(int *y)
{
  int tot=0;
  for(int i=0;i<4;i++) tot+=y[i];
  if (((tot + 100)%2) == 0) return 1;
  else return -1;
}

//------------------------------------------------------------------
// CTOR
//------------------------------------------------------------------


NLSMesonPropS::NLSMesonPropS(Lattice& lattice, NLStagMesonArg& narg)
: HadronPropS(lattice, 15, narg.dir, QuarkPropSMng::srcSlice(narg.qid0[0]),1) ,
  qp00(QuarkPropSMng::prop(narg.qid0[0])),
  qp01(QuarkPropSMng::prop(narg.qid0[1])),
  qp02(QuarkPropSMng::prop(narg.qid0[2])),
  qp03(QuarkPropSMng::prop(narg.qid0[3])),
  qp04(QuarkPropSMng::prop(narg.qid0[4])),
  qp05(QuarkPropSMng::prop(narg.qid0[5])),
  qp06(QuarkPropSMng::prop(narg.qid0[6])),
  qp07(QuarkPropSMng::prop(narg.qid0[7]))
{
  char *fname = "NLSMesonPropS(Lattice&, NLSMesonArg&)";
  nx[0] = GJP.XnodeSites();
  nx[1] = GJP.YnodeSites();
  nx[2] = GJP.ZnodeSites();
  nx[3] = GJP.TnodeSites();

  int i = 0; int j = 1; int k = 2;
  qp0[0] = qp00;
  qp0[1] = qp01;
  qp0[2] = qp02;
  qp0[3] = qp03;
  qp0[4] = qp04;
  qp0[5] = qp05;
  qp0[6] = qp06;
  qp0[7] = qp07;

  int bufferSize = VECT_LEN * (8 + ( nx[i] + nx[j] + nx[k] ) * 4 + 2 * (nx[i]*nx[j] + nx[i]*nx[k] + nx[j]*nx[k]) ); 

/*
// Only for the test purpose!!
  Float *testbuffer[3];
  Float *testpropagator[3];  
  FILE *testfile;
  testfile = fopen("test.out","w");
  //set up the buffer for test
  for(int color = 0; color < 3; color++) 
        testbuffer[color] = (Float *)smalloc(bufferSize*sizeof(Float));
  //set up the propagator for test
  int buffersize = nx[0] * nx[1] * nx[2] * nx[3] * VECT_LEN;
  for(color = 0; color < 3; color++) 
        testpropagator[color] = (Float *)smalloc(buffersize*sizeof(Float));
  //fill the test propagator
*/
/*  int x[4];
  int node_origin[4] = {  GJP.XnodeCoor() * nx[0], GJP.YnodeCoor() * nx[1], GJP.ZnodeCoor() * nx[2], GJP.TnodeCoor() * nx[3] };
  for(x[0]=0;x[0]<nx[0];x[0]++)
    for(x[1]=0;x[1]<nx[1];x[1]++)
      for(x[2]=0;x[2]<nx[2];x[2]++)    
        for(x[3]=0;x[3]<nx[3];x[3]++)    
          for( color=0;color<3;color++)  
            for( int l=0; l<VECT_LEN ; l++)   
             *(testpropagator[color]+X_OFFSET(x)+l) = x[(l%4)] + node_origin[(l%4)] + 0.1*color;
  for(x[3]=0;x[3]<nx[3];x[3]++) {
    transfer(testpropagator,testbuffer,x[3]);   
    for(x[0]=-1;x[0]<=nx[0];x[0]++)
      for(x[1]=-1;x[1]<=nx[1];x[1]++)
        for(x[2]=-1;x[2]<=nx[2];x[2]++)  {
          for( color=0;color<3;color++)  { 
            for( int l=0; l<VECT_LEN ; l++)   
              if (IsOffNode(x)) fprintf(testfile,"%d:%1.1f ",x[(l%4)] + node_origin[(l%4)],*(testbuffer[color]+x_map(x)+l));    
              else fprintf(testfile,"%d:%1.1f ",x[(l%4)] + node_origin[(l%4)],*(testpropagator[color]+X_OFFSET(x)+l));    
            fprintf(testfile," next color ");
            }
          fprintf(testfile,"\n");
          }
  }
  for(color=0;color<3;color++) {
    sfree(testpropagator[color]);
    sfree(testbuffer[color]);
  }          
  bufferSize = VECT_LEN * (8 + ( nx[i] + nx[j] + nx[k] ) * 4 + 2 * (nx[i]*nx[j] + nx[i]*nx[k] + nx[j]*nx[k]) ); 
  fclose(testfile);
//end of test  
*/

  for(int src_i = 0 ; src_i < 8 ; src_i++ )
    for(int color = 0; color < 3; color++) {
      if(isDegQuarks()) {
        buffer0[src_i][color] = (Float *)smalloc(bufferSize*sizeof(Float));
        if(buffer0[src_i][color] == 0)
          ERR.Pointer(cname,fname, "buffer0[color]");
          VRB.Smalloc(cname,fname, "buffer0[color]", buffer0[color], 
    		 bufferSize*sizeof(Float));
        }
        else {
          ERR.Pointer(cname,fname, "Not Degerate Quark!!");
        }
     } 


}

//------------------------------------------------------------------
// DTOR
//------------------------------------------------------------------
NLSMesonPropS::~NLSMesonPropS()
{
  char *fname = "~NLSMesonPropS()";
  if(isDegQuarks()) {
    for(int src_i = 0; src_i < 8 ; src_i++)
      for(int color = 0; color < 3; color++) {
        VRB.Sfree(cname,fname, "buffer0[color]",buffer0[color]);
        sfree(buffer0[src_i][color]);
      }
  }
  else {
        ERR.Pointer(cname,fname, "Not Degerate Quark!!");
    }
}
CPS_END_NAMESPACE



