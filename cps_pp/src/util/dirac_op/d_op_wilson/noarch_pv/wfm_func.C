#include <config.h>
CPS_START_NAMESPACE
/***************************************************************************/
/*                                                                         */
/* Various optimized functions for the wilson operator                     */
/* The "in", "out" vectors are half-spinnors (12 components).              */
/* The "inf", "outf" are full-spinors (24 components).                     */
/*                                                                         */
/***************************************************************************/

CPS_END_NAMESPACE
#include <util/data_types.h>
#include <util/vector.h>
#include <util/wilson.h>
#include <util/error.h>
#include <sys/bgl/bgl_sys.h>
#include <comms/scu.h>
CPS_START_NAMESPACE

#define U(r,row,col,d)  *(u+(r+2*(row+3*(col+3*d))))

#define IN(r,c,s)       *(in   +(r+2*(c+3*s)))
#define INF(r,c,s)      *(inf  +(r+2*(c+3*s)))

#define OUT(r,c,s)      *(out  +(r+2*(c+3*s)))
#define OUTF(r,c,s)     *(outf +(r+2*(c+3*s)))

#define IN0(r,c,s)      *(in0  +(r+2*(c+3*s)))
#define IN1(r,c,s)      *(in1  +(r+2*(c+3*s)))
#define IN2(r,c,s)      *(in2  +(r+2*(c+3*s)))
#define IN3(r,c,s)      *(in3  +(r+2*(c+3*s)))

#define OUT0(r,c,s)     *(out0 +(r+2*(c+3*s)))
#define OUT1(r,c,s)     *(out1 +(r+2*(c+3*s)))
#define OUT2(r,c,s)     *(out2 +(r+2*(c+3*s)))
#define OUT3(r,c,s)     *(out3 +(r+2*(c+3*s)))



/***********************************************************************/
/* Color multiply the half spinors "in" with U_mu. Result in "out".    */ 
/***********************************************************************/
void wilson_dslash_csmat(double *out, double *u, double *in, int mu)
{
  int s;
  int c;

  for(s=0;s<2;s++){
    for(c=0;c<3;c++){

      OUT(0,c,s) = (  U(0,c,0,mu) * IN(0,0,s)
		    + U(0,c,1,mu) * IN(0,1,s)
		    + U(0,c,2,mu) * IN(0,2,s) 
		    - U(1,c,0,mu) * IN(1,0,s)
		    - U(1,c,1,mu) * IN(1,1,s)
		    - U(1,c,2,mu) * IN(1,2,s) );
      OUT(1,c,s) = (  U(0,c,0,mu) * IN(1,0,s)
		    + U(0,c,1,mu) * IN(1,1,s)
		    + U(0,c,2,mu) * IN(1,2,s) 
		    + U(1,c,0,mu) * IN(0,0,s)
		    + U(1,c,1,mu) * IN(0,1,s)
		    + U(1,c,2,mu) * IN(0,2,s) );
    }
  }

}


/***********************************************************************/
/* Color multiply the half spinors "in" with U_mu^dag. Result in "out".*/ 
/***********************************************************************/
void wilson_dslash_csmatdag(double *out, double *u, double *in, int mu)
{
  int s;
  int c;

  for(s=0;s<2;s++){
    for(c=0;c<3;c++){

      OUT(0,c,s) = (  U(0,0,c,mu) * IN(0,0,s)
		    + U(0,1,c,mu) * IN(0,1,s)
		    + U(0,2,c,mu) * IN(0,2,s) 
		    + U(1,0,c,mu) * IN(1,0,s)
		    + U(1,1,c,mu) * IN(1,1,s)
		    + U(1,2,c,mu) * IN(1,2,s) );
      OUT(1,c,s) = (  U(0,0,c,mu) * IN(1,0,s)
		    + U(0,1,c,mu) * IN(1,1,s)
		    + U(0,2,c,mu) * IN(1,2,s) 
		    - U(1,0,c,mu) * IN(0,0,s)
		    - U(1,1,c,mu) * IN(0,1,s)
		    - U(1,2,c,mu) * IN(0,2,s) );
    }
  }

}


/***********************************************************************/
/* Spin project with [1 + sign * gamma] the full spinor "inf" onto the */
/* four half spinors "out0", "out1", "out2", "out3".                   */
/***********************************************************************/
void wilson_dslash_spproj(double *out0, double *out1, double *out2, double *out3, 
			  double *inf, double sign)
{
  int c;

  for(c=0;c<3;c++){

    /* spin project with [1 + sign * gamma_0] */
    OUT0(0,c,0) = INF(0,c,0) + sign * ( -INF(1,c,3) ); 
    OUT0(1,c,0) = INF(1,c,0) + sign * (  INF(0,c,3) ); 
    
    OUT0(0,c,1) = INF(0,c,1) + sign * ( -INF(1,c,2) ); 
    OUT0(1,c,1) = INF(1,c,1) + sign * (  INF(0,c,2) ); 
    
    /* spin project with [1 + sign * gamma_1] */
    OUT1(0,c,0) = INF(0,c,0) + sign * ( -INF(0,c,3) ); 
    OUT1(1,c,0) = INF(1,c,0) + sign * ( -INF(1,c,3) ); 
    
    OUT1(0,c,1) = INF(0,c,1) + sign * (  INF(0,c,2) ); 
    OUT1(1,c,1) = INF(1,c,1) + sign * (  INF(1,c,2) ); 
    
    /* spin project with [1 + sign * gamma_2] */
    OUT2(0,c,0) = INF(0,c,0) + sign * ( -INF(1,c,2) ); 
    OUT2(1,c,0) = INF(1,c,0) + sign * (  INF(0,c,2) ); 
    
    OUT2(0,c,1) = INF(0,c,1) + sign * (  INF(1,c,3) ); 
    OUT2(1,c,1) = INF(1,c,1) + sign * ( -INF(0,c,3) ); 
    
    /* spin project with [1 + sign * gamma_3] */
    OUT3(0,c,0) = INF(0,c,0) + sign * (  INF(0,c,2) ); 
    OUT3(1,c,0) = INF(1,c,0) + sign * (  INF(1,c,2) ); 
    
    OUT3(0,c,1) = INF(0,c,1) + sign * (  INF(0,c,3) ); 
    OUT3(1,c,1) = INF(1,c,1) + sign * (  INF(1,c,3) ); 
  }

}


/***********************************************************************/
/* Expand the half spinors "in" to a full spinor "outf" (trick), where */
/* the spin projection was done with [1+sign*gamma].                   */
/* If acc=0 store the result into the full spinor "outf".              */
/* If acc=1 accumulate the result into the full spinor "outf".         */
/***********************************************************************/
void wilson_dslash_trick(double *outf, double *in0, double *in1, double *in2, double *in3, 
			 double sign, int accum)
{
  int c;

  if(accum == 0) {
    for(c=0;c<3;c++){
      /* trick into spin component 0*/
      OUTF(0,c,0) =          IN0(0,c,0) + IN1(0,c,0) + IN2(0,c,0) + IN3(0,c,0);
      OUTF(1,c,0) =          IN0(1,c,0) + IN1(1,c,0) + IN2(1,c,0) + IN3(1,c,0); 
      
      /* trick into spin component 1*/
      OUTF(0,c,1) =          IN0(0,c,1) + IN1(0,c,1) + IN2(0,c,1) + IN3(0,c,1); 
      OUTF(1,c,1) =          IN0(1,c,1) + IN1(1,c,1) + IN2(1,c,1) + IN3(1,c,1); 
      
      /* trick into spin component 2*/
      OUTF(0,c,2) = sign * ( IN0(1,c,1) + IN1(0,c,1) + IN2(1,c,0) + IN3(0,c,0) );
      OUTF(1,c,2) = sign * (-IN0(0,c,1) + IN1(1,c,1) - IN2(0,c,0) + IN3(1,c,0) );
      
      /* trick into spin component 3*/
      OUTF(0,c,3) = sign * ( IN0(1,c,0) - IN1(0,c,0) - IN2(1,c,1) + IN3(0,c,1) );
      OUTF(1,c,3) = sign * (-IN0(0,c,0) - IN1(1,c,0) + IN2(0,c,1) + IN3(1,c,1) );
    }

  } else {
    for(c=0;c<3;c++){
      /* trick and accumulate into spin component 0*/
      OUTF(0,c,0) = OUTF(0,c,0) +          IN0(0,c,0) + IN1(0,c,0) + IN2(0,c,0) + IN3(0,c,0);
      OUTF(1,c,0) = OUTF(1,c,0) +          IN0(1,c,0) + IN1(1,c,0) + IN2(1,c,0) + IN3(1,c,0); 
      
      /* trick and accumulate into spin component 1*/
      OUTF(0,c,1) = OUTF(0,c,1) +          IN0(0,c,1) + IN1(0,c,1) + IN2(0,c,1) + IN3(0,c,1); 
      OUTF(1,c,1) = OUTF(1,c,1) +          IN0(1,c,1) + IN1(1,c,1) + IN2(1,c,1) + IN3(1,c,1); 
      
      /* trick and accumulate into spin component 2*/
      OUTF(0,c,2) = OUTF(0,c,2) + sign * ( IN0(1,c,1) + IN1(0,c,1) + IN2(1,c,0) + IN3(0,c,0) );
      OUTF(1,c,2) = OUTF(1,c,2) + sign * (-IN0(0,c,1) + IN1(1,c,1) - IN2(0,c,0) + IN3(1,c,0) );
      
      /* trick and accumulate into spin component 3*/
      OUTF(0,c,3) = OUTF(0,c,3) + sign * ( IN0(1,c,0) - IN1(0,c,0) - IN2(1,c,1) + IN3(0,c,1) );
      OUTF(1,c,3) = OUTF(1,c,3) + sign * (-IN0(0,c,0) - IN1(1,c,0) + IN2(0,c,1) + IN3(1,c,1) );
    }

  }

}


CPS_END_NAMESPACE
