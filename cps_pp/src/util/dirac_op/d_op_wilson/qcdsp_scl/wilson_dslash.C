#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_scl/wilson_dslash.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/***************************************************************************/
/*                                                                         */
/* wilson_dslash: It calculates chi = Dslash * psi, or                     */
/*                                chi = DslashDag * psi, where Dslassh is  */
/* the Wilson fermion Dslash matrix. Dslash is a function of the gauge     */
/* fields u.                                                               */
/* cb = 0/1 denotes input on even/odd checkerboards i.e.                   */
/* cb = 1 --> Dslash_EO, acts on an odd column vector and produces an      */
/* even column vector                                                      */
/* cb = 0 --> Dslash_OE, acts on an even column vector and produces an     */
/* odd column vector,                                                      */
/*                                                                         */
/* dag = 0/1  results into calculating Dslash or DslashDagger.             */
/* lx,ly,lz,lt is the lattice size.                                        */
/*                                                                         */
/* This routine is to be used with scalar machines.                        */
/*                                                                         */
/* WARNING:                                                                 */
/*                                                                          */
/* This set of routines will work only if the node sublattices have         */
/* even number of sites in each direction.                                  */
/*                                                                          */
/***************************************************************************/

CPS_END_NAMESPACE
#include<util/data_types.h>
#include<util/wilson.h>
#include<util/error.h>
CPS_START_NAMESPACE


#define U(r,row,col,d,n,cb) *(u+(r+2*(row+3*(col+3*(d+4*(n+vol*(cb)))))))
#define PSI(r,c,s,n)     *(psi+(r+2*(c+3*(s+4*(n)))))
#define CHI(r,c,s,n)     *(chi+(r+2*(c+3*(s+4*(n)))))
#define TMP(r,c,s)     *(tmp+(r+2*(c+3*(s))))
#define TMP1(r,c,s)     *(tmp1+(r+2*(c+3*(s))))
#define TMP2(r,c,s)     *(tmp2+(r+2*(c+3*(s))))
#define TMP3(r,c,s)     *(tmp3+(r+2*(c+3*(s))))
#define TMP4(r,c,s)     *(tmp4+(r+2*(c+3*(s))))
#define TMP5(r,c,s)     *(tmp5+(r+2*(c+3*(s))))
#define TMP6(r,c,s)     *(tmp6+(r+2*(c+3*(s))))
#define TMP7(r,c,s)     *(tmp7+(r+2*(c+3*(s))))
#define TMP8(r,c,s)     *(tmp8+(r+2*(c+3*(s))))



/*--------------------------------------------------------------------------*/
/* wilson_dslash:                                                           */
/*--------------------------------------------------------------------------*/
void wilson_dslash(IFloat *chi_f, 
		   IFloat *u_f, 
		   IFloat *psi_f, 
		   int cb,
		   int dag,
		   Wilson *wilson_p)
{
  char *fname = "wilson_dslash";
  int lx, ly, lz, lt;
  int x, y, z, t;
  int xp, yp, zp, tp;
  int xm, ym, zm, tm;
  int xyzt;
  int xpyzt, xypzt, xyzpt, xyztp;
  int xmyzt, xymzt, xyzmt, xyztm;
  int cbn;
  int parity;
  int sdag;                 /* = +/-1 if dag = 0/1 */
  int r, c, s, mu;
  int vol;
  Float tmp[SPINOR_SIZE];
  Float tmp1[SPINOR_SIZE];
  Float tmp2[SPINOR_SIZE];
  Float tmp3[SPINOR_SIZE];
  Float tmp4[SPINOR_SIZE];
  Float tmp5[SPINOR_SIZE];
  Float tmp6[SPINOR_SIZE];
  Float tmp7[SPINOR_SIZE];
  Float tmp8[SPINOR_SIZE];
  
/*--------------------------------------------------------------------------*/
/* Initializations                                                          */
/*--------------------------------------------------------------------------*/
  Float *chi = (Float *) chi_f;
  Float *u = (Float *) u_f;
  Float *psi = (Float *) psi_f;

  lx = wilson_p->ptr[0];
  ly = wilson_p->ptr[1];
  lz = wilson_p->ptr[2];
  lt = wilson_p->ptr[3];
  vol = wilson_p->vol[0];

  if(dag == 0)
    sdag = 1;
  else if(dag == 1)
    sdag = -1;
  else{
    ERR.General(" ",fname,"dag must be 0 or 1");
  }

  if(cb == 0)
    cbn = 1;
  else if(cb == 1)
    cbn = 0;
  else{
    ERR.General(" ",fname,"cb must be 0 or 1");
  }

/*--------------------------------------------------------------------------*/
/* Loop over sites                                                          */
/*--------------------------------------------------------------------------*/
  for(x=0; x<lx; x++){
    for(y=0; y<ly; y++){
      for(z=0; z<lz; z++){
	for(t=0; t<lt; t++){
	 parity = x+y+z+t;
	 parity = parity % 2;
	 if(parity == cbn){

	   /* x,y,z,t addressing of cbn checkerboard */
	   xyzt = (x/2)+(lx/2)*(y+ly*(z+lz*t));

	   /* x,y,z,t addressing of cb checkerboard */
	   xp = (x+1) % lx;
	   yp = (y+1) % ly;
	   zp = (z+1) % lz;
	   tp = (t+1) % lt;
	   xm = x-1 + ( (lx-x)/lx ) * lx;
	   ym = y-1 + ( (ly-y)/ly ) * ly;
	   zm = z-1 + ( (lz-z)/lz ) * lz;
	   tm = t-1 + ( (lt-t)/lt ) * lt;
	   xpyzt = (xp/2)+(lx/2)*(y+ly*(z+lz*t));
	   xmyzt = (xm/2)+(lx/2)*(y+ly*(z+lz*t));
	   xypzt = (x/2)+(lx/2)*(yp+ly*(z+lz*t));
	   xymzt = (x/2)+(lx/2)*(ym+ly*(z+lz*t));
	   xyzpt = (x/2)+(lx/2)*(y+ly*(zp+lz*t));
	   xyzmt = (x/2)+(lx/2)*(y+ly*(zm+lz*t));
	   xyztp = (x/2)+(lx/2)*(y+ly*(z+lz*tp));
	   xyztm = (x/2)+(lx/2)*(y+ly*(z+lz*tm));


	   /* 1-gamma_0 */
	   /*-----------*/
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = PSI(0,c,0,xpyzt) - sdag * ( -PSI(1,c,3,xpyzt) ); 
	     TMP(1,c,0) = PSI(1,c,0,xpyzt) - sdag * (  PSI(0,c,3,xpyzt) ); 

	     TMP(0,c,1) = PSI(0,c,1,xpyzt) - sdag * ( -PSI(1,c,2,xpyzt) ); 
	     TMP(1,c,1) = PSI(1,c,1,xpyzt) - sdag * (  PSI(0,c,2,xpyzt) ); 

	     TMP(0,c,2) = PSI(0,c,2,xpyzt) - sdag * (  PSI(1,c,1,xpyzt) ); 
	     TMP(1,c,2) = PSI(1,c,2,xpyzt) - sdag * ( -PSI(0,c,1,xpyzt) ); 

	     TMP(0,c,3) = PSI(0,c,3,xpyzt) - sdag * (  PSI(1,c,0,xpyzt) ); 
	     TMP(1,c,3) = PSI(1,c,3,xpyzt) - sdag * ( -PSI(0,c,0,xpyzt) ); 
	   }
	   /* multiply by U_mu */
	   mu = 0;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       TMP1(0,c,s) = (  U(0,c,0,mu,xyzt,cbn) * TMP(0,0,s)
			      + U(0,c,1,mu,xyzt,cbn) * TMP(0,1,s)
			      + U(0,c,2,mu,xyzt,cbn) * TMP(0,2,s) 
			      - U(1,c,0,mu,xyzt,cbn) * TMP(1,0,s)
			      - U(1,c,1,mu,xyzt,cbn) * TMP(1,1,s)
			      - U(1,c,2,mu,xyzt,cbn) * TMP(1,2,s) );
	       TMP1(1,c,s) = (  U(0,c,0,mu,xyzt,cbn) * TMP(1,0,s)
			      + U(0,c,1,mu,xyzt,cbn) * TMP(1,1,s)
			      + U(0,c,2,mu,xyzt,cbn) * TMP(1,2,s) 
			      + U(1,c,0,mu,xyzt,cbn) * TMP(0,0,s)
			      + U(1,c,1,mu,xyzt,cbn) * TMP(0,1,s)
			      + U(1,c,2,mu,xyzt,cbn) * TMP(0,2,s) );
	     }
	   }


	   /* 1-gamma_1 */
	   /*-----------*/
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = PSI(0,c,0,xypzt) - sdag * ( -PSI(0,c,3,xypzt) ); 
	     TMP(1,c,0) = PSI(1,c,0,xypzt) - sdag * ( -PSI(1,c,3,xypzt) ); 

	     TMP(0,c,1) = PSI(0,c,1,xypzt) - sdag * (  PSI(0,c,2,xypzt) ); 
	     TMP(1,c,1) = PSI(1,c,1,xypzt) - sdag * (  PSI(1,c,2,xypzt) ); 

	     TMP(0,c,2) = PSI(0,c,2,xypzt) - sdag * (  PSI(0,c,1,xypzt) ); 
	     TMP(1,c,2) = PSI(1,c,2,xypzt) - sdag * (  PSI(1,c,1,xypzt) ); 

	     TMP(0,c,3) = PSI(0,c,3,xypzt) - sdag * ( -PSI(0,c,0,xypzt) ); 
	     TMP(1,c,3) = PSI(1,c,3,xypzt) - sdag * ( -PSI(1,c,0,xypzt) ); 
	   }
	   /* multiply by U_mu */
	   mu = 1;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       TMP2(0,c,s) = (  U(0,c,0,mu,xyzt,cbn) * TMP(0,0,s)
			      + U(0,c,1,mu,xyzt,cbn) * TMP(0,1,s)
			      + U(0,c,2,mu,xyzt,cbn) * TMP(0,2,s) 
			      - U(1,c,0,mu,xyzt,cbn) * TMP(1,0,s)
			      - U(1,c,1,mu,xyzt,cbn) * TMP(1,1,s)
			      - U(1,c,2,mu,xyzt,cbn) * TMP(1,2,s) );
	       TMP2(1,c,s) = (  U(0,c,0,mu,xyzt,cbn) * TMP(1,0,s)
			      + U(0,c,1,mu,xyzt,cbn) * TMP(1,1,s)
			      + U(0,c,2,mu,xyzt,cbn) * TMP(1,2,s) 
			      + U(1,c,0,mu,xyzt,cbn) * TMP(0,0,s)
			      + U(1,c,1,mu,xyzt,cbn) * TMP(0,1,s)
			      + U(1,c,2,mu,xyzt,cbn) * TMP(0,2,s) );
	     }
	   }


	   /* 1-gamma_2 */
	   /*-----------*/
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = PSI(0,c,0,xyzpt) - sdag * ( -PSI(1,c,2,xyzpt) ); 
	     TMP(1,c,0) = PSI(1,c,0,xyzpt) - sdag * (  PSI(0,c,2,xyzpt) ); 

	     TMP(0,c,1) = PSI(0,c,1,xyzpt) - sdag * (  PSI(1,c,3,xyzpt) ); 
	     TMP(1,c,1) = PSI(1,c,1,xyzpt) - sdag * ( -PSI(0,c,3,xyzpt) ); 

	     TMP(0,c,2) = PSI(0,c,2,xyzpt) - sdag * (  PSI(1,c,0,xyzpt) ); 
	     TMP(1,c,2) = PSI(1,c,2,xyzpt) - sdag * ( -PSI(0,c,0,xyzpt) ); 

	     TMP(0,c,3) = PSI(0,c,3,xyzpt) - sdag * ( -PSI(1,c,1,xyzpt) ); 
	     TMP(1,c,3) = PSI(1,c,3,xyzpt) - sdag * (  PSI(0,c,1,xyzpt) ); 
	   }
	   /* multiply by U_mu */
	   mu = 2;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       TMP3(0,c,s) = (  U(0,c,0,mu,xyzt,cbn) * TMP(0,0,s)
			      + U(0,c,1,mu,xyzt,cbn) * TMP(0,1,s)
			      + U(0,c,2,mu,xyzt,cbn) * TMP(0,2,s) 
			      - U(1,c,0,mu,xyzt,cbn) * TMP(1,0,s)
			      - U(1,c,1,mu,xyzt,cbn) * TMP(1,1,s)
			      - U(1,c,2,mu,xyzt,cbn) * TMP(1,2,s) );
	       TMP3(1,c,s) = (  U(0,c,0,mu,xyzt,cbn) * TMP(1,0,s)
			      + U(0,c,1,mu,xyzt,cbn) * TMP(1,1,s)
			      + U(0,c,2,mu,xyzt,cbn) * TMP(1,2,s) 
			      + U(1,c,0,mu,xyzt,cbn) * TMP(0,0,s)
			      + U(1,c,1,mu,xyzt,cbn) * TMP(0,1,s)
			      + U(1,c,2,mu,xyzt,cbn) * TMP(0,2,s) );
	     }
	   }


	   /* 1-gamma_3 */
	   /*-----------*/
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = PSI(0,c,0,xyztp) - sdag * (  PSI(0,c,2,xyztp) ); 
	     TMP(1,c,0) = PSI(1,c,0,xyztp) - sdag * (  PSI(1,c,2,xyztp) ); 

	     TMP(0,c,1) = PSI(0,c,1,xyztp) - sdag * (  PSI(0,c,3,xyztp) ); 
	     TMP(1,c,1) = PSI(1,c,1,xyztp) - sdag * (  PSI(1,c,3,xyztp) ); 

	     TMP(0,c,2) = PSI(0,c,2,xyztp) - sdag * (  PSI(0,c,0,xyztp) ); 
	     TMP(1,c,2) = PSI(1,c,2,xyztp) - sdag * (  PSI(1,c,0,xyztp) ); 

	     TMP(0,c,3) = PSI(0,c,3,xyztp) - sdag * (  PSI(0,c,1,xyztp) ); 
	     TMP(1,c,3) = PSI(1,c,3,xyztp) - sdag * (  PSI(1,c,1,xyztp) ); 
	   }
	   /* multiply by U_mu */
	   mu = 3;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       TMP4(0,c,s) = (  U(0,c,0,mu,xyzt,cbn) * TMP(0,0,s)
			      + U(0,c,1,mu,xyzt,cbn) * TMP(0,1,s)
			      + U(0,c,2,mu,xyzt,cbn) * TMP(0,2,s) 
			      - U(1,c,0,mu,xyzt,cbn) * TMP(1,0,s)
			      - U(1,c,1,mu,xyzt,cbn) * TMP(1,1,s)
			      - U(1,c,2,mu,xyzt,cbn) * TMP(1,2,s) );
	       TMP4(1,c,s) = (  U(0,c,0,mu,xyzt,cbn) * TMP(1,0,s)
			      + U(0,c,1,mu,xyzt,cbn) * TMP(1,1,s)
			      + U(0,c,2,mu,xyzt,cbn) * TMP(1,2,s) 
			      + U(1,c,0,mu,xyzt,cbn) * TMP(0,0,s)
			      + U(1,c,1,mu,xyzt,cbn) * TMP(0,1,s)
			      + U(1,c,2,mu,xyzt,cbn) * TMP(0,2,s) );
	     }
	   }


	   /* 1+gamma_0 */
	   /*-----------*/
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = PSI(0,c,0,xmyzt) + sdag * ( -PSI(1,c,3,xmyzt) ); 
	     TMP(1,c,0) = PSI(1,c,0,xmyzt) + sdag * (  PSI(0,c,3,xmyzt) ); 

	     TMP(0,c,1) = PSI(0,c,1,xmyzt) + sdag * ( -PSI(1,c,2,xmyzt) ); 
	     TMP(1,c,1) = PSI(1,c,1,xmyzt) + sdag * (  PSI(0,c,2,xmyzt) ); 

	     TMP(0,c,2) = PSI(0,c,2,xmyzt) + sdag * (  PSI(1,c,1,xmyzt) ); 
	     TMP(1,c,2) = PSI(1,c,2,xmyzt) + sdag * ( -PSI(0,c,1,xmyzt) ); 

	     TMP(0,c,3) = PSI(0,c,3,xmyzt) + sdag * (  PSI(1,c,0,xmyzt) ); 
	     TMP(1,c,3) = PSI(1,c,3,xmyzt) + sdag * ( -PSI(0,c,0,xmyzt) ); 
	   }
	   /* multiply by U_mu */
	   mu = 0;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       TMP5(0,c,s) = (  U(0,0,c,mu,xmyzt,cb) * TMP(0,0,s)
			      + U(0,1,c,mu,xmyzt,cb) * TMP(0,1,s)
			      + U(0,2,c,mu,xmyzt,cb) * TMP(0,2,s) 
			      + U(1,0,c,mu,xmyzt,cb) * TMP(1,0,s)
			      + U(1,1,c,mu,xmyzt,cb) * TMP(1,1,s)
			      + U(1,2,c,mu,xmyzt,cb) * TMP(1,2,s) );
	       TMP5(1,c,s) = (  U(0,0,c,mu,xmyzt,cb) * TMP(1,0,s)
			      + U(0,1,c,mu,xmyzt,cb) * TMP(1,1,s)
			      + U(0,2,c,mu,xmyzt,cb) * TMP(1,2,s) 
			      - U(1,0,c,mu,xmyzt,cb) * TMP(0,0,s)
			      - U(1,1,c,mu,xmyzt,cb) * TMP(0,1,s)
			      - U(1,2,c,mu,xmyzt,cb) * TMP(0,2,s) );
	     }
	   }


	   /* 1+gamma_1 */
	   /*-----------*/
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = PSI(0,c,0,xymzt) + sdag * ( -PSI(0,c,3,xymzt) ); 
	     TMP(1,c,0) = PSI(1,c,0,xymzt) + sdag * ( -PSI(1,c,3,xymzt) ); 

	     TMP(0,c,1) = PSI(0,c,1,xymzt) + sdag * (  PSI(0,c,2,xymzt) ); 
	     TMP(1,c,1) = PSI(1,c,1,xymzt) + sdag * (  PSI(1,c,2,xymzt) ); 

	     TMP(0,c,2) = PSI(0,c,2,xymzt) + sdag * (  PSI(0,c,1,xymzt) ); 
	     TMP(1,c,2) = PSI(1,c,2,xymzt) + sdag * (  PSI(1,c,1,xymzt) ); 

	     TMP(0,c,3) = PSI(0,c,3,xymzt) + sdag * ( -PSI(0,c,0,xymzt) ); 
	     TMP(1,c,3) = PSI(1,c,3,xymzt) + sdag * ( -PSI(1,c,0,xymzt) ); 
	   }
	   /* multiply by U_mu */
	   mu = 1;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       TMP6(0,c,s) = (  U(0,0,c,mu,xymzt,cb) * TMP(0,0,s)
			      + U(0,1,c,mu,xymzt,cb) * TMP(0,1,s)
			      + U(0,2,c,mu,xymzt,cb) * TMP(0,2,s) 
			      + U(1,0,c,mu,xymzt,cb) * TMP(1,0,s)
			      + U(1,1,c,mu,xymzt,cb) * TMP(1,1,s)
			      + U(1,2,c,mu,xymzt,cb) * TMP(1,2,s) );
	       TMP6(1,c,s) = (  U(0,0,c,mu,xymzt,cb) * TMP(1,0,s)
			      + U(0,1,c,mu,xymzt,cb) * TMP(1,1,s)
			      + U(0,2,c,mu,xymzt,cb) * TMP(1,2,s) 
			      - U(1,0,c,mu,xymzt,cb) * TMP(0,0,s)
			      - U(1,1,c,mu,xymzt,cb) * TMP(0,1,s)
			      - U(1,2,c,mu,xymzt,cb) * TMP(0,2,s) );
	     }
	   }


	   /* 1+gamma_2 */
	   /*-----------*/
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = PSI(0,c,0,xyzmt) + sdag * ( -PSI(1,c,2,xyzmt) ); 
	     TMP(1,c,0) = PSI(1,c,0,xyzmt) + sdag * (  PSI(0,c,2,xyzmt) ); 

	     TMP(0,c,1) = PSI(0,c,1,xyzmt) + sdag * (  PSI(1,c,3,xyzmt) ); 
	     TMP(1,c,1) = PSI(1,c,1,xyzmt) + sdag * ( -PSI(0,c,3,xyzmt) ); 

	     TMP(0,c,2) = PSI(0,c,2,xyzmt) + sdag * (  PSI(1,c,0,xyzmt) ); 
	     TMP(1,c,2) = PSI(1,c,2,xyzmt) + sdag * ( -PSI(0,c,0,xyzmt) ); 

	     TMP(0,c,3) = PSI(0,c,3,xyzmt) + sdag * ( -PSI(1,c,1,xyzmt) ); 
	     TMP(1,c,3) = PSI(1,c,3,xyzmt) + sdag * (  PSI(0,c,1,xyzmt) ); 
	   }
	   /* multiply by U_mu */
	   mu = 2;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       TMP7(0,c,s) = (  U(0,0,c,mu,xyzmt,cb) * TMP(0,0,s)
			      + U(0,1,c,mu,xyzmt,cb) * TMP(0,1,s)
			      + U(0,2,c,mu,xyzmt,cb) * TMP(0,2,s) 
			      + U(1,0,c,mu,xyzmt,cb) * TMP(1,0,s)
			      + U(1,1,c,mu,xyzmt,cb) * TMP(1,1,s)
			      + U(1,2,c,mu,xyzmt,cb) * TMP(1,2,s) );
	       TMP7(1,c,s) = (  U(0,0,c,mu,xyzmt,cb) * TMP(1,0,s)
			      + U(0,1,c,mu,xyzmt,cb) * TMP(1,1,s)
			      + U(0,2,c,mu,xyzmt,cb) * TMP(1,2,s) 
			      - U(1,0,c,mu,xyzmt,cb) * TMP(0,0,s)
			      - U(1,1,c,mu,xyzmt,cb) * TMP(0,1,s)
			      - U(1,2,c,mu,xyzmt,cb) * TMP(0,2,s) );
	     }
	   }


	   /* 1+gamma_3 */
	   /*-----------*/
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = PSI(0,c,0,xyztm) + sdag * (  PSI(0,c,2,xyztm) ); 
	     TMP(1,c,0) = PSI(1,c,0,xyztm) + sdag * (  PSI(1,c,2,xyztm) ); 

	     TMP(0,c,1) = PSI(0,c,1,xyztm) + sdag * (  PSI(0,c,3,xyztm) ); 
	     TMP(1,c,1) = PSI(1,c,1,xyztm) + sdag * (  PSI(1,c,3,xyztm) ); 

	     TMP(0,c,2) = PSI(0,c,2,xyztm) + sdag * (  PSI(0,c,0,xyztm) ); 
	     TMP(1,c,2) = PSI(1,c,2,xyztm) + sdag * (  PSI(1,c,0,xyztm) ); 

	     TMP(0,c,3) = PSI(0,c,3,xyztm) + sdag * (  PSI(0,c,1,xyztm) ); 
	     TMP(1,c,3) = PSI(1,c,3,xyztm) + sdag * (  PSI(1,c,1,xyztm) ); 
	   }
	   /* multiply by U_mu */
	   mu = 3;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       TMP8(0,c,s) = (  U(0,0,c,mu,xyztm,cb) * TMP(0,0,s)
			      + U(0,1,c,mu,xyztm,cb) * TMP(0,1,s)
			      + U(0,2,c,mu,xyztm,cb) * TMP(0,2,s) 
			      + U(1,0,c,mu,xyztm,cb) * TMP(1,0,s)
			      + U(1,1,c,mu,xyztm,cb) * TMP(1,1,s)
			      + U(1,2,c,mu,xyztm,cb) * TMP(1,2,s) );
	       TMP8(1,c,s) = (  U(0,0,c,mu,xyztm,cb) * TMP(1,0,s)
			      + U(0,1,c,mu,xyztm,cb) * TMP(1,1,s)
			      + U(0,2,c,mu,xyztm,cb) * TMP(1,2,s) 
			      - U(1,0,c,mu,xyztm,cb) * TMP(0,0,s)
			      - U(1,1,c,mu,xyztm,cb) * TMP(0,1,s)
			      - U(1,2,c,mu,xyztm,cb) * TMP(0,2,s) );
	     }
	   }



	   /* Add all contributions */

	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       for(r=0;r<2;r++){
		 CHI(r,c,s,xyzt) = (  TMP1(r,c,s)
				    + TMP2(r,c,s)
				    + TMP3(r,c,s)
				    + TMP4(r,c,s)
				    + TMP5(r,c,s)
				    + TMP6(r,c,s)
				    + TMP7(r,c,s)
				    + TMP8(r,c,s) );
	       }
	     }
	   }




	 }
       }
      }
    }
  }
  
  



}

CPS_END_NAMESPACE
