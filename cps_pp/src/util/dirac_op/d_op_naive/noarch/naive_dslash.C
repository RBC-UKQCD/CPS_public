#include <config.h>
#include <stdio.h>
#include <math.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.

  $Id: naive_dslash.C,v 1.2 2013-05-28 14:52:58 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-05-28 14:52:58 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_naive/noarch/naive_dslash.C,v 1.2 2013-05-28 14:52:58 chulwoo Exp $
//  $Id: naive_dslash.C,v 1.2 2013-05-28 14:52:58 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_naive/noarch/naive_dslash.C,v $
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
#include <util/data_types.h>
#include <util/vector.h>
#include <util/naive.h>
#include <util/error.h>
#include <comms/scu.h>
CPS_START_NAMESPACE


//! Access to the elements of the \e SU(3) matrix
/*!
  Gets the element of the \e SU(3) matrix \e u with row \e row,
  column \e col and complex component \e d
*/
#define U(r,row,col,d)  *(u+(r+2*(row+3*(col+3*d))))
//! Access to the elements of a spinor vector.
/*!
  Gets the element of the spinor \e psi with spin \e s,
  colour \e c and complex component \e r
*/
#define PSI(r,c,s)      *(psi +(r+2*(c+3*s)))

//! As above, but the vector is called chi
#define CHI(r,c,s)      *(chi +(r+2*(c+3*s)))
#define TMP(r,c,s)      *(tmp +(r+2*(c+3*s))) 
#define TMP1(r,c,s)     *(tmp1+(r+2*(c+3*s))) 
#define TMP2(r,c,s)     *(tmp2+(r+2*(c+3*s))) 
#define TMP3(r,c,s)     *(tmp3+(r+2*(c+3*s))) 
#define TMP4(r,c,s)     *(tmp4+(r+2*(c+3*s))) 
#define TMP5(r,c,s)     *(tmp5+(r+2*(c+3*s))) 
#define TMP6(r,c,s)     *(tmp6+(r+2*(c+3*s))) 
#define TMP7(r,c,s)     *(tmp7+(r+2*(c+3*s))) 
#define TMP8(r,c,s)     *(tmp8+(r+2*(c+3*s))) 
#define FBUF(r,c,s)     *(fbuf+(r+2*(c+3*s))) 






void naive_dslash(IFloat *chi_p_f, 
		   IFloat *u_p_f, 
		   IFloat *psi_p_f, 
		   int cb,
		   int dag,
		   Wilson *wilson_p)
{
  char *fname = "naive_dslash";
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
  Float fbuf[SPINOR_SIZE];
  
/*--------------------------------------------------------------------------*/
/* Initializations                                                          */
/*--------------------------------------------------------------------------*/
  Float *chi_p = (Float *) chi_p_f;
  Float *u_p = (Float *) u_p_f;
  Float *psi_p = (Float *) psi_p_f;
  Float *chi;
  Float *u;
  Float *psi;


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
	   u   = u_p + GAUGE_SIZE * ( xyzt  + vol * cbn);
	   psi = psi_p + SPINOR_SIZE * xpyzt;
	   if(x == lx-1){
	     getPlusData((IFloat *)fbuf, (IFloat *) psi, SPINOR_SIZE, 0);
	     psi = fbuf;
	   }
#if 0
 	     for(int l = 0; l<2;l++)
 	     for(int m = 0; m<3;m++)
 	     for(int n = 0; n<3;n++)
	     if( fabs(PSI(l,m,n)) > 1e-10) {
		printf("PSI_xp(%d %d %d %d) (%d %d %d) = %e\n",
		x,y,z,t,l,m,n,PSI(l,m,n));
 	     }
#endif
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = 0.*(PSI(0,c,0)) - sdag * ( -PSI(1,c,3) ); 
	     TMP(1,c,0) = 0.*(PSI(1,c,0)) - sdag * (  PSI(0,c,3) ); 

	     TMP(0,c,1) = 0.*(PSI(0,c,1)) - sdag * ( -PSI(1,c,2) ); 
	     TMP(1,c,1) = 0.*(PSI(1,c,1)) - sdag * (  PSI(0,c,2) ); 

	     TMP(0,c,2) = 0.*(PSI(0,c,2)) - sdag * (  PSI(1,c,1) ); 
	     TMP(1,c,2) = 0.*(PSI(1,c,2)) - sdag * ( -PSI(0,c,1) ); 

	     TMP(0,c,3) = 0.*(PSI(0,c,3)) - sdag * (  PSI(1,c,0) ); 
	     TMP(1,c,3) = 0.*(PSI(1,c,3)) - sdag * ( -PSI(0,c,0) ); 
	   }
	   /* multiply by U_mu */
	   mu = 0;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       TMP1(0,c,s) = (  U(0,c,0,mu) * TMP(0,0,s)
			      + U(0,c,1,mu) * TMP(0,1,s)
			      + U(0,c,2,mu) * TMP(0,2,s) 
			      - U(1,c,0,mu) * TMP(1,0,s)
			      - U(1,c,1,mu) * TMP(1,1,s)
			      - U(1,c,2,mu) * TMP(1,2,s) );
	       TMP1(1,c,s) = (  U(0,c,0,mu) * TMP(1,0,s)
			      + U(0,c,1,mu) * TMP(1,1,s)
			      + U(0,c,2,mu) * TMP(1,2,s) 
			      + U(1,c,0,mu) * TMP(0,0,s)
			      + U(1,c,1,mu) * TMP(0,1,s)
			      + U(1,c,2,mu) * TMP(0,2,s) );
	     }
	   }


	   /* 1-gamma_1 */
	   /*-----------*/
	   u   = u_p + GAUGE_SIZE * ( xyzt  + vol * cbn);
	   psi = psi_p + SPINOR_SIZE * xypzt;
	   if(y == ly-1){
	     getPlusData((IFloat *)fbuf, (IFloat *) psi, SPINOR_SIZE, 1);
	     psi = fbuf;
#if 0
 	     for(int l = 0; l<2;l++)
 	     for(int m = 0; m<3;m++)
 	     for(int n = 0; n<3;n++)
	     if( fabs(PSI(l,m,n)) > 1e-10) {
		printf("PSI_yp(%d %d %d %d) (%d %d %d) = %e\n",
		x,y,z,t,l,m,n,PSI(l,m,n));
 	     }
#endif
	   }
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = 0.*(PSI(0,c,0)) - sdag * ( -PSI(0,c,3) ); 
	     TMP(1,c,0) = 0.*(PSI(1,c,0)) - sdag * ( -PSI(1,c,3) ); 

	     TMP(0,c,1) = 0.*(PSI(0,c,1)) - sdag * (  PSI(0,c,2) ); 
	     TMP(1,c,1) = 0.*(PSI(1,c,1)) - sdag * (  PSI(1,c,2) ); 

	     TMP(0,c,2) = 0.*(PSI(0,c,2)) - sdag * (  PSI(0,c,1) ); 
	     TMP(1,c,2) = 0.*(PSI(1,c,2)) - sdag * (  PSI(1,c,1) ); 

	     TMP(0,c,3) = 0.*(PSI(0,c,3)) - sdag * ( -PSI(0,c,0) ); 
	     TMP(1,c,3) = 0.*(PSI(1,c,3)) - sdag * ( -PSI(1,c,0) ); 
	   }
	   /* multiply by U_mu */
	   mu = 1;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       TMP2(0,c,s) = (  U(0,c,0,mu) * TMP(0,0,s)
			      + U(0,c,1,mu) * TMP(0,1,s)
			      + U(0,c,2,mu) * TMP(0,2,s) 
			      - U(1,c,0,mu) * TMP(1,0,s)
			      - U(1,c,1,mu) * TMP(1,1,s)
			      - U(1,c,2,mu) * TMP(1,2,s) );
	       TMP2(1,c,s) = (  U(0,c,0,mu) * TMP(1,0,s)
			      + U(0,c,1,mu) * TMP(1,1,s)
			      + U(0,c,2,mu) * TMP(1,2,s) 
			      + U(1,c,0,mu) * TMP(0,0,s)
			      + U(1,c,1,mu) * TMP(0,1,s)
			      + U(1,c,2,mu) * TMP(0,2,s) );
	     }
	   }


	   /* 1-gamma_2 */
	   /*-----------*/
	   u   = u_p + GAUGE_SIZE * ( xyzt  + vol * cbn);
	   psi = psi_p + SPINOR_SIZE * xyzpt;
	   if(z == lz-1){
	     getPlusData((IFloat *)fbuf, (IFloat *) psi, SPINOR_SIZE, 2);
	     psi = fbuf;
#if 0
 	     for(int l = 0; l<2;l++)
 	     for(int m = 0; m<3;m++)
 	     for(int n = 0; n<3;n++)
	     if( fabs(PSI(l,m,n)) > 1e-10) {
		printf("PSI_zp(%d %d %d %d) (%d %d %d) = %e\n",
		x,y,z,t,l,m,n,PSI(l,m,n));
 	     }
#endif
	   }
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = 0.*(PSI(0,c,0)) - sdag * ( -PSI(1,c,2) ); 
	     TMP(1,c,0) = 0.*(PSI(1,c,0)) - sdag * (  PSI(0,c,2) ); 

	     TMP(0,c,1) = 0.*(PSI(0,c,1)) - sdag * (  PSI(1,c,3) ); 
	     TMP(1,c,1) = 0.*(PSI(1,c,1)) - sdag * ( -PSI(0,c,3) ); 

	     TMP(0,c,2) = 0.*(PSI(0,c,2)) - sdag * (  PSI(1,c,0) ); 
	     TMP(1,c,2) = 0.*(PSI(1,c,2)) - sdag * ( -PSI(0,c,0) ); 

	     TMP(0,c,3) = 0.*(PSI(0,c,3)) - sdag * ( -PSI(1,c,1) ); 
	     TMP(1,c,3) = 0.*(PSI(1,c,3)) - sdag * (  PSI(0,c,1) ); 
	   }
	   /* multiply by U_mu */
	   mu = 2;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       TMP3(0,c,s) = (  U(0,c,0,mu) * TMP(0,0,s)
			      + U(0,c,1,mu) * TMP(0,1,s)
			      + U(0,c,2,mu) * TMP(0,2,s) 
			      - U(1,c,0,mu) * TMP(1,0,s)
			      - U(1,c,1,mu) * TMP(1,1,s)
			      - U(1,c,2,mu) * TMP(1,2,s) );
	       TMP3(1,c,s) = (  U(0,c,0,mu) * TMP(1,0,s)
			      + U(0,c,1,mu) * TMP(1,1,s)
			      + U(0,c,2,mu) * TMP(1,2,s) 
			      + U(1,c,0,mu) * TMP(0,0,s)
			      + U(1,c,1,mu) * TMP(0,1,s)
			      + U(1,c,2,mu) * TMP(0,2,s) );
	     }
	   }


	   /* 1-gamma_3 */
	   /*-----------*/
	   u   = u_p + GAUGE_SIZE * ( xyzt  + vol * cbn);
	   psi = psi_p + SPINOR_SIZE * xyztp;
	   if(t == lt-1){
	     getPlusData((IFloat *)fbuf, (IFloat *) psi, SPINOR_SIZE, 3);
	     psi = fbuf;
#if 0
 	     for(int l = 0; l<2;l++)
 	     for(int m = 0; m<3;m++)
 	     for(int n = 0; n<3;n++)
	     if( fabs(PSI(l,m,n)) > 1e-10) {
		printf("PSI_tp(%d %d %d %d) (%d %d %d) = %e\n",
		x,y,z,t,l,m,n,PSI(l,m,n));
 	     }
#endif
	   }
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = 0.*(PSI(0,c,0)) - sdag * (  PSI(0,c,2) ); 
	     TMP(1,c,0) = 0.*(PSI(1,c,0)) - sdag * (  PSI(1,c,2) ); 

	     TMP(0,c,1) = 0.*(PSI(0,c,1)) - sdag * (  PSI(0,c,3) ); 
	     TMP(1,c,1) = 0.*(PSI(1,c,1)) - sdag * (  PSI(1,c,3) ); 

	     TMP(0,c,2) = 0.*(PSI(0,c,2)) - sdag * (  PSI(0,c,0) ); 
	     TMP(1,c,2) = 0.*(PSI(1,c,2)) - sdag * (  PSI(1,c,0) ); 

	     TMP(0,c,3) = 0.*(PSI(0,c,3)) - sdag * (  PSI(0,c,1) ); 
	     TMP(1,c,3) = 0.*(PSI(1,c,3)) - sdag * (  PSI(1,c,1) ); 
	   }
	   /* multiply by U_mu */
	   mu = 3;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       TMP4(0,c,s) = (  U(0,c,0,mu) * TMP(0,0,s)
			      + U(0,c,1,mu) * TMP(0,1,s)
			      + U(0,c,2,mu) * TMP(0,2,s) 
			      - U(1,c,0,mu) * TMP(1,0,s)
			      - U(1,c,1,mu) * TMP(1,1,s)
			      - U(1,c,2,mu) * TMP(1,2,s) );
	       TMP4(1,c,s) = (  U(0,c,0,mu) * TMP(1,0,s)
			      + U(0,c,1,mu) * TMP(1,1,s)
			      + U(0,c,2,mu) * TMP(1,2,s) 
			      + U(1,c,0,mu) * TMP(0,0,s)
			      + U(1,c,1,mu) * TMP(0,1,s)
			      + U(1,c,2,mu) * TMP(0,2,s) );
	     }
	   }


	   /* 1+gamma_0 */
	   /*-----------*/
	   u   = u_p + GAUGE_SIZE * ( xmyzt  + vol * cb);
	   psi = psi_p + SPINOR_SIZE * xmyzt;
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = 0.*(PSI(0,c,0)) + sdag * ( -PSI(1,c,3) ); 
	     TMP(1,c,0) = 0.*(PSI(1,c,0)) + sdag * (  PSI(0,c,3) ); 

	     TMP(0,c,1) = 0.*(PSI(0,c,1)) + sdag * ( -PSI(1,c,2) ); 
	     TMP(1,c,1) = 0.*(PSI(1,c,1)) + sdag * (  PSI(0,c,2) ); 

	     TMP(0,c,2) = 0.*(PSI(0,c,2)) + sdag * (  PSI(1,c,1) ); 
	     TMP(1,c,2) = 0.*(PSI(1,c,2)) + sdag * ( -PSI(0,c,1) ); 

	     TMP(0,c,3) = 0.*(PSI(0,c,3)) + sdag * (  PSI(1,c,0) ); 
	     TMP(1,c,3) = 0.*(PSI(1,c,3)) + sdag * ( -PSI(0,c,0) ); 
	   }
	   /* multiply by U_mu */
	   mu = 0;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       TMP5(0,c,s) = (  U(0,0,c,mu) * TMP(0,0,s)
			      + U(0,1,c,mu) * TMP(0,1,s)
			      + U(0,2,c,mu) * TMP(0,2,s) 
			      + U(1,0,c,mu) * TMP(1,0,s)
			      + U(1,1,c,mu) * TMP(1,1,s)
			      + U(1,2,c,mu) * TMP(1,2,s) );
	       TMP5(1,c,s) = (  U(0,0,c,mu) * TMP(1,0,s)
			      + U(0,1,c,mu) * TMP(1,1,s)
			      + U(0,2,c,mu) * TMP(1,2,s) 
			      - U(1,0,c,mu) * TMP(0,0,s)
			      - U(1,1,c,mu) * TMP(0,1,s)
			      - U(1,2,c,mu) * TMP(0,2,s) );
	     }
	   }
	   if(x == 0){
	     getMinusData((IFloat *)fbuf, (IFloat *)tmp5, SPINOR_SIZE, 0);
	     moveMem((IFloat *)tmp5, (IFloat *)fbuf, 
			SPINOR_SIZE*sizeof(Float) / sizeof(char));
	   }



	   /* 1+gamma_1 */
	   /*-----------*/
	   u   = u_p + GAUGE_SIZE * ( xymzt  + vol * cb);
	   psi = psi_p + SPINOR_SIZE * xymzt;
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = 0.*(PSI(0,c,0)) + sdag * ( -PSI(0,c,3) ); 
	     TMP(1,c,0) = 0.*(PSI(1,c,0)) + sdag * ( -PSI(1,c,3) ); 

	     TMP(0,c,1) = 0.*(PSI(0,c,1)) + sdag * (  PSI(0,c,2) ); 
	     TMP(1,c,1) = 0.*(PSI(1,c,1)) + sdag * (  PSI(1,c,2) ); 

	     TMP(0,c,2) = 0.*(PSI(0,c,2)) + sdag * (  PSI(0,c,1) ); 
	     TMP(1,c,2) = 0.*(PSI(1,c,2)) + sdag * (  PSI(1,c,1) ); 

	     TMP(0,c,3) = 0.*(PSI(0,c,3)) + sdag * ( -PSI(0,c,0) ); 
	     TMP(1,c,3) = 0.*(PSI(1,c,3)) + sdag * ( -PSI(1,c,0) ); 
	   }
	   /* multiply by U_mu */
	   mu = 1;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       TMP6(0,c,s) = (  U(0,0,c,mu) * TMP(0,0,s)
			      + U(0,1,c,mu) * TMP(0,1,s)
			      + U(0,2,c,mu) * TMP(0,2,s) 
			      + U(1,0,c,mu) * TMP(1,0,s)
			      + U(1,1,c,mu) * TMP(1,1,s)
			      + U(1,2,c,mu) * TMP(1,2,s) );
	       TMP6(1,c,s) = (  U(0,0,c,mu) * TMP(1,0,s)
			      + U(0,1,c,mu) * TMP(1,1,s)
			      + U(0,2,c,mu) * TMP(1,2,s) 
			      - U(1,0,c,mu) * TMP(0,0,s)
			      - U(1,1,c,mu) * TMP(0,1,s)
			      - U(1,2,c,mu) * TMP(0,2,s) );
	     }
	   }
	   if(y == 0){
	     getMinusData((IFloat *)fbuf, (IFloat *)tmp6, SPINOR_SIZE, 1);
	     moveMem((IFloat *)tmp6, (IFloat *)fbuf, 
		SPINOR_SIZE*sizeof(Float) / sizeof(char));
	   }


	   /* 1+gamma_2 */
	   /*-----------*/
	   u   = u_p + GAUGE_SIZE * ( xyzmt  + vol * cb);
	   psi = psi_p + SPINOR_SIZE * xyzmt;
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = 0.*(PSI(0,c,0)) + sdag * ( -PSI(1,c,2) ); 
	     TMP(1,c,0) = 0.*(PSI(1,c,0)) + sdag * (  PSI(0,c,2) ); 

	     TMP(0,c,1) = 0.*(PSI(0,c,1)) + sdag * (  PSI(1,c,3) ); 
	     TMP(1,c,1) = 0.*(PSI(1,c,1)) + sdag * ( -PSI(0,c,3) ); 

	     TMP(0,c,2) = 0.*(PSI(0,c,2)) + sdag * (  PSI(1,c,0) ); 
	     TMP(1,c,2) = 0.*(PSI(1,c,2)) + sdag * ( -PSI(0,c,0) ); 

	     TMP(0,c,3) = 0.*(PSI(0,c,3)) + sdag * ( -PSI(1,c,1) ); 
	     TMP(1,c,3) = 0.*(PSI(1,c,3)) + sdag * (  PSI(0,c,1) ); 
	   }
	   /* multiply by U_mu */
	   mu = 2;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       TMP7(0,c,s) = (  U(0,0,c,mu) * TMP(0,0,s)
			      + U(0,1,c,mu) * TMP(0,1,s)
			      + U(0,2,c,mu) * TMP(0,2,s) 
			      + U(1,0,c,mu) * TMP(1,0,s)
			      + U(1,1,c,mu) * TMP(1,1,s)
			      + U(1,2,c,mu) * TMP(1,2,s) );
	       TMP7(1,c,s) = (  U(0,0,c,mu) * TMP(1,0,s)
			      + U(0,1,c,mu) * TMP(1,1,s)
			      + U(0,2,c,mu) * TMP(1,2,s) 
			      - U(1,0,c,mu) * TMP(0,0,s)
			      - U(1,1,c,mu) * TMP(0,1,s)
			      - U(1,2,c,mu) * TMP(0,2,s) );
	     }
	   }
	   if(z == 0){
	     getMinusData((IFloat *)fbuf, (IFloat *)tmp7, SPINOR_SIZE, 2);
	     moveMem((IFloat *)tmp7, (IFloat *)fbuf, 
		SPINOR_SIZE*sizeof(Float) / sizeof(char));
	   }


	   /* 1+gamma_3 */
	   /*-----------*/
	   u   = u_p + GAUGE_SIZE * ( xyztm  + vol * cb);
	   psi = psi_p + SPINOR_SIZE * xyztm;
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = 0.*(PSI(0,c,0)) + sdag * (  PSI(0,c,2) ); 
	     TMP(1,c,0) = 0.*(PSI(1,c,0)) + sdag * (  PSI(1,c,2) ); 

	     TMP(0,c,1) = 0.*(PSI(0,c,1)) + sdag * (  PSI(0,c,3) ); 
	     TMP(1,c,1) = 0.*(PSI(1,c,1)) + sdag * (  PSI(1,c,3) ); 

	     TMP(0,c,2) = 0.*(PSI(0,c,2)) + sdag * (  PSI(0,c,0) ); 
	     TMP(1,c,2) = 0.*(PSI(1,c,2)) + sdag * (  PSI(1,c,0) ); 

	     TMP(0,c,3) = 0.*(PSI(0,c,3)) + sdag * (  PSI(0,c,1) ); 
	     TMP(1,c,3) = 0.*(PSI(1,c,3)) + sdag * (  PSI(1,c,1) ); 
	   }
	   /* multiply by U_mu */
	   mu = 3;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       TMP8(0,c,s) = (  U(0,0,c,mu) * TMP(0,0,s)
			      + U(0,1,c,mu) * TMP(0,1,s)
			      + U(0,2,c,mu) * TMP(0,2,s) 
			      + U(1,0,c,mu) * TMP(1,0,s)
			      + U(1,1,c,mu) * TMP(1,1,s)
			      + U(1,2,c,mu) * TMP(1,2,s) );
	       TMP8(1,c,s) = (  U(0,0,c,mu) * TMP(1,0,s)
			      + U(0,1,c,mu) * TMP(1,1,s)
			      + U(0,2,c,mu) * TMP(1,2,s) 
			      - U(1,0,c,mu) * TMP(0,0,s)
			      - U(1,1,c,mu) * TMP(0,1,s)
			      - U(1,2,c,mu) * TMP(0,2,s) );
	     }
	   }
	   if(t == 0){
	     getMinusData((IFloat *)fbuf, (IFloat *)tmp8, SPINOR_SIZE, 3);
	     moveMem((IFloat *)tmp8, (IFloat *)fbuf, 
		SPINOR_SIZE*sizeof(Float) / sizeof(char));
	   }



	   /* Add all contributions */

	   chi = chi_p + SPINOR_SIZE * xyzt;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       for(r=0;r<2;r++){
		 CHI(r,c,s) = (  TMP1(r,c,s)
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

void naive_dslash_two(Float *chi0, Float *chi1,
                   Float *u,
                   Float *psi0, Float *psi1,
                   int cb0, int cb1,
                   int dag,
                   Wilson *wp)
{
  naive_dslash(chi0,u,psi0,cb0,dag,wp);
  naive_dslash(chi1,u,psi1,cb1,dag,wp);
}

CPS_END_NAMESPACE
