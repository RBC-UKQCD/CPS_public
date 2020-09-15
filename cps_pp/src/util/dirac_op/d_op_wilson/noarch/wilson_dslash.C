#include <config.h>
#if (!defined USE_SSE) && (!defined SSE_TO_C)
#include <stdio.h>
#include <math.h>
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.

*/
//--------------------------------------------------------------------
//
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

#include <util/data_types.h>
#include <util/dirac_op.h>
#include <util/vector.h>
#include <util/wilson.h>
#include <util/error.h>
#include <comms/scu.h>

#include <util/gjp.h>

#ifndef USE_QMP //if USE_QMP then it compiles the vectorised wilson dslash instead

CPS_START_NAMESPACE
//! Access to the elements of the \e SU(3) matrix
/*!
  Gets the element of the \e SU(3) matrix \e u with row \e row,
  column \e col and complex component \e d
*/
#define UELEM(u,r,row,col,d) *(u+(r+2*(row+3*(col+3*d))))
//! Access to the elements of a spinor vector.
/*!
  Gets the element of the spinor \e f with spin \e s,
  colour \e c and complex component \e r
*/
#define FERM(f,r,c,s) *(f +(r+2*(c+3*s)))


//pm is minus (0) or plus (1)
static void OneMinusPlusGammaTimesPsi(Float *into, Float *psi, const int &sdag, const int &mu, const int &mp){
  Float sign = -1.0;
  if(mp==1) sign = 1.0;
  
  int c;

  if(mu==0){
    //1+-GAMMA_X
    for(c=0;c<3;c++){
      FERM(into,0,c,0) = FERM(psi,0,c,0) + sign*sdag * ( -FERM(psi,1,c,3) ); 
      FERM(into,1,c,0) = FERM(psi,1,c,0) + sign*sdag * (  FERM(psi,0,c,3) ); 

      FERM(into,0,c,1) = FERM(psi,0,c,1) + sign*sdag * ( -FERM(psi,1,c,2) ); 
      FERM(into,1,c,1) = FERM(psi,1,c,1) + sign*sdag * (  FERM(psi,0,c,2) ); 

      FERM(into,0,c,2) = FERM(psi,0,c,2) + sign*sdag * (  FERM(psi,1,c,1) ); 
      FERM(into,1,c,2) = FERM(psi,1,c,2) + sign*sdag * ( -FERM(psi,0,c,1) ); 

      FERM(into,0,c,3) = FERM(psi,0,c,3) + sign*sdag * (  FERM(psi,1,c,0) ); 
      FERM(into,1,c,3) = FERM(psi,1,c,3) + sign*sdag * ( -FERM(psi,0,c,0) ); 
    }
  }else if(mu==1){
    for(c=0;c<3;c++){
      FERM(into,0,c,0) = FERM(psi,0,c,0) + sign*sdag * ( -FERM(psi,0,c,3) ); 
      FERM(into,1,c,0) = FERM(psi,1,c,0) + sign*sdag * ( -FERM(psi,1,c,3) ); 

      FERM(into,0,c,1) = FERM(psi,0,c,1) + sign*sdag * (  FERM(psi,0,c,2) ); 
      FERM(into,1,c,1) = FERM(psi,1,c,1) + sign*sdag * (  FERM(psi,1,c,2) ); 

      FERM(into,0,c,2) = FERM(psi,0,c,2) + sign*sdag * (  FERM(psi,0,c,1) ); 
      FERM(into,1,c,2) = FERM(psi,1,c,2) + sign*sdag * (  FERM(psi,1,c,1) ); 

      FERM(into,0,c,3) = FERM(psi,0,c,3) + sign*sdag * ( -FERM(psi,0,c,0) ); 
      FERM(into,1,c,3) = FERM(psi,1,c,3) + sign*sdag * ( -FERM(psi,1,c,0) ); 
    }
  }else if(mu==2){
    for(c=0;c<3;c++){
      FERM(into,0,c,0) = FERM(psi,0,c,0) + sign*sdag * ( -FERM(psi,1,c,2) ); 
      FERM(into,1,c,0) = FERM(psi,1,c,0) + sign*sdag * (  FERM(psi,0,c,2) ); 

      FERM(into,0,c,1) = FERM(psi,0,c,1) + sign*sdag * (  FERM(psi,1,c,3) ); 
      FERM(into,1,c,1) = FERM(psi,1,c,1) + sign*sdag * ( -FERM(psi,0,c,3) ); 

      FERM(into,0,c,2) = FERM(psi,0,c,2) + sign*sdag * (  FERM(psi,1,c,0) ); 
      FERM(into,1,c,2) = FERM(psi,1,c,2) + sign*sdag * ( -FERM(psi,0,c,0) ); 

      FERM(into,0,c,3) = FERM(psi,0,c,3) + sign*sdag * ( -FERM(psi,1,c,1) ); 
      FERM(into,1,c,3) = FERM(psi,1,c,3) + sign*sdag * (  FERM(psi,0,c,1) ); 
    }
  }else if(mu==3){
    for(c=0;c<3;c++){
      FERM(into,0,c,0) = FERM(psi,0,c,0) + sign*sdag * (  FERM(psi,0,c,2) ); 
      FERM(into,1,c,0) = FERM(psi,1,c,0) + sign*sdag * (  FERM(psi,1,c,2) ); 

      FERM(into,0,c,1) = FERM(psi,0,c,1) + sign*sdag * (  FERM(psi,0,c,3) ); 
      FERM(into,1,c,1) = FERM(psi,1,c,1) + sign*sdag * (  FERM(psi,1,c,3) ); 

      FERM(into,0,c,2) = FERM(psi,0,c,2) + sign*sdag * (  FERM(psi,0,c,0) ); 
      FERM(into,1,c,2) = FERM(psi,1,c,2) + sign*sdag * (  FERM(psi,1,c,0) ); 

      FERM(into,0,c,3) = FERM(psi,0,c,3) + sign*sdag * (  FERM(psi,0,c,1) ); 
      FERM(into,1,c,3) = FERM(psi,1,c,3) + sign*sdag * (  FERM(psi,1,c,1) ); 
    }
  }
}

//rowcol is whether the multiplication is on the row or column index of the gauge field
static void UmuTimes(Float *into, Float *u, Float *psi, const int &mu, const int &rowcol){
  if(rowcol==0){
    for(int s=0;s<4;s++){
      for(int c=0;c<3;c++){
	FERM(into,0,c,s) = (  UELEM(u,0,c,0,mu) * FERM(psi,0,0,s)
			      + UELEM(u,0,c,1,mu) * FERM(psi,0,1,s)
			      + UELEM(u,0,c,2,mu) * FERM(psi,0,2,s) 
			      - UELEM(u,1,c,0,mu) * FERM(psi,1,0,s)
			      - UELEM(u,1,c,1,mu) * FERM(psi,1,1,s)
			      - UELEM(u,1,c,2,mu) * FERM(psi,1,2,s) );
	FERM(into,1,c,s) = (  UELEM(u,0,c,0,mu) * FERM(psi,1,0,s)
			      + UELEM(u,0,c,1,mu) * FERM(psi,1,1,s)
			      + UELEM(u,0,c,2,mu) * FERM(psi,1,2,s) 
			      + UELEM(u,1,c,0,mu) * FERM(psi,0,0,s)
			      + UELEM(u,1,c,1,mu) * FERM(psi,0,1,s)
			      + UELEM(u,1,c,2,mu) * FERM(psi,0,2,s) );
      }
    }
  }else{
    for(int s=0;s<4;s++){
      for(int c=0;c<3;c++){
	FERM(into,0,c,s) = (  UELEM(u,0,0,c,mu) * FERM(psi,0,0,s)
			 + UELEM(u,0,1,c,mu) * FERM(psi,0,1,s)
			 + UELEM(u,0,2,c,mu) * FERM(psi,0,2,s) 
			 + UELEM(u,1,0,c,mu) * FERM(psi,1,0,s)
			 + UELEM(u,1,1,c,mu) * FERM(psi,1,1,s)
			 + UELEM(u,1,2,c,mu) * FERM(psi,1,2,s) );
	FERM(into,1,c,s) = (  UELEM(u,0,0,c,mu) * FERM(psi,1,0,s)
			 + UELEM(u,0,1,c,mu) * FERM(psi,1,1,s)
			 + UELEM(u,0,2,c,mu) * FERM(psi,1,2,s) 
			 - UELEM(u,1,0,c,mu) * FERM(psi,0,0,s)
			 - UELEM(u,1,1,c,mu) * FERM(psi,0,1,s)
			 - UELEM(u,1,2,c,mu) * FERM(psi,0,2,s) );
      }
    }
  }
}

//CK: vastly improved wilson code + added G-parity functionality
//Calculates the wilson Dslash for the part of the fermion field on checkerboard cb (not the mass term)
//Note that due to the fact that the massless Dslash at a site x uses only the fields at x +- mu, not at x,
//the actual sites at which the Dlash is evaluated have the opposite parity to cb.
	void wilson_dslash(IFloat *chi_p_f, 
		   IFloat *u_p_f, 
		   IFloat *psi_p_f, 
		   int cb,
		   int dag,
		   Wilson *wilson_p)
{
  char *fname = "wilson_dslash";
  int lx, ly, lz, lt;
  int x, y, z, t;
  int xyzt;
  int cbn;
  int parity;
  int sdag;                 /* = +/-1 if dag = 0/1 */
  int r, c, s, mu;
  int vol;

  int buf_size = SPINOR_SIZE;
  if(GJP.Gparity()) buf_size*=2;
  int flav_fbuf_offset = SPINOR_SIZE;

  Float tmp[buf_size];
  Float tmp1[buf_size];
  Float tmp2[buf_size];
  Float tmp3[buf_size];
  Float tmp4[buf_size];
  Float tmp5[buf_size];
  Float tmp6[buf_size];
  Float tmp7[buf_size];
  Float tmp8[buf_size];
  Float fbuf[buf_size];

/*--------------------------------------------------------------------------*/
/* Initializations                                                          */
/*--------------------------------------------------------------------------*/
  Float *chi_p = (Float *) chi_p_f;
  Float *u_p = (Float *) u_p_f;
  Float *psi_p = (Float *) psi_p_f;
  Float *chi;
  Float *u;
  Float *psi;

  //for G-parity
  Float *chi_f1;
  Float *u_f1;
  Float *psi_f1;

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

  //DEBUG
  // for(x=0; x<lx; x++){
  //   for(y=0; y<ly; y++){
  //     for(z=0; z<lz; z++){
  // 	for(t=0; t<lt; t++){
  // 	  xyzt = (x/2)+(lx/2)*(y+ly*(z+lz*t));
  // 	  psi = psi_p + SPINOR_SIZE * xyzt;
  // 	  if(GJP.Gparity()){
  // 	    printf("x=(%d,%d,%d,%d)  psi0_x = (%f,%f)\n",x,y,z,t,*psi,*(psi+1));
  // 	  }else{
  // 	    int gpx = x; 
  // 	    int gpf = 0;
  // 	    if(gpx>=lx/2){ gpx-=lx/2; gpf = 1; }
  // 	    printf("x=(%d,%d,%d,%d)  psi%d_x = (%f,%f)\n",gpx,y,z,t,gpf,*psi,*(psi+1));
  // 	    fflush(stdout);
  // 	  }
  // 	}
  //     }
  //   }
  // }
  // if(GJP.Gparity()){
  //   for(x=0; x<lx; x++){
  //     for(y=0; y<ly; y++){
  // 	for(z=0; z<lz; z++){
  // 	  for(t=0; t<lt; t++){
  // 	    xyzt = (x/2)+(lx/2)*(y+ly*(z+lz*t));
  // 	    psi_f1 = psi_p + SPINOR_SIZE * (xyzt + vol);
  // 	    printf("x=(%d,%d,%d,%d)  psi1_x = (%f,%f)\n",x,y,z,t,*psi_f1,*(psi_f1+1));
  // 	  }
  // 	}
  //     }
  //   }
  // }
  //DEBUG

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

	   int pos[4] = {x,y,z,t};
	   int lattsz[4] = {lx,ly,lz,lt};
	   Float *temps[4][2] = { {tmp1,tmp5}, {tmp2,tmp6}, {tmp3,tmp7}, {tmp4,tmp8} };

	   int u_cboff = vol; //vol is the checkerboard volume, i.e. half the 4d volume
	   if(GJP.Gparity()) u_cboff*=2; //[(f0 odd)(f1 odd)(f0 even)(f1 even)]  each bracket is one cbvolume

	   for(int mu=0;mu<4;mu++){
	     /* x,y,z,t addressing of cb checkerboard */
	     int posp[4] = {x,y,z,t};
	     posp[mu] = (posp[mu]+1)%lattsz[mu];

	     int posp_xyzt = (posp[0]/2)+(lx/2)*(posp[1]+ly*(posp[2]+lz*posp[3]));

	     /* 1-gamma_mu */
	     /*-----------*/

	     u   = u_p + GAUGE_SIZE * ( xyzt  + u_cboff * cbn);
	     psi = psi_p + SPINOR_SIZE * posp_xyzt;

	     if(GJP.Gparity()){
	       u_f1 = u_p + GAUGE_SIZE * ( xyzt  + u_cboff * cbn + vol);
	       psi_f1 = psi_p + SPINOR_SIZE * (posp_xyzt + vol);
	     }

	     if(pos[mu] == lattsz[mu]-1){
	       if(!GJP.Gparity()){
		 getPlusData((IFloat *)fbuf, (IFloat *) psi, SPINOR_SIZE, mu);
		 psi = fbuf;
	       }else{
		 IFloat* psi_f0_recv_to = (IFloat*)fbuf; //f0 receives f0
		 IFloat* psi_f1_recv_to = (IFloat*)(fbuf + flav_fbuf_offset); //f1 receives f1
		 if(GJP.Bc(mu) == BND_CND_GPARITY && GJP.NodeCoor(mu)==GJP.Nodes(mu)-1){
		   //do the G-parity flavour twist. (assumes minus sign is on the gauge links on the boundary)
		   psi_f0_recv_to = (IFloat*)(fbuf + flav_fbuf_offset); //f0 receives f1
		   psi_f1_recv_to = (IFloat*)fbuf; //f1 receives f0
		 }
		 getPlusData(psi_f0_recv_to, (IFloat *) psi, SPINOR_SIZE, mu);
		 getPlusData(psi_f1_recv_to, (IFloat *) psi_f1, SPINOR_SIZE, mu);
		 psi = fbuf;//is either f0 or f1 depending on how the buffers were set
		 psi_f1 = fbuf+flav_fbuf_offset;
	       }
	     }

	     OneMinusPlusGammaTimesPsi(tmp,psi,sdag,mu,0);
	     if(GJP.Gparity()) OneMinusPlusGammaTimesPsi(tmp+SPINOR_SIZE,psi_f1,sdag,mu,0);
	     
	     /* multiply by U_mu */
	     UmuTimes(temps[mu][0],u,tmp,mu,0);
	     if(GJP.Gparity()) UmuTimes(temps[mu][0]+SPINOR_SIZE,u_f1,tmp+SPINOR_SIZE,mu,0);
	   }
	   for(int mu=0;mu<4;mu++){
	     /* x,y,z,t addressing of cb checkerboard */
	     int posm[4] = {x,y,z,t};
	     posm[mu] = posm[mu]-1 + ( (lattsz[mu]-posm[mu])/lattsz[mu] ) * lattsz[mu];
	     
	     //looks like checkerboard is divided across the x-direction
	     int posm_xyzt = (posm[0]/2)+(lx/2)*(posm[1]+ly*(posm[2]+lz*posm[3]));

	     /* 1+gamma_mu */
	     /*-----------*/

	     u   = u_p + GAUGE_SIZE * ( posm_xyzt  + u_cboff * cb);
	     psi = psi_p + SPINOR_SIZE * posm_xyzt;
	   
	     if(GJP.Gparity()){
	       u_f1 = u_p + GAUGE_SIZE * ( posm_xyzt  + u_cboff * cb + vol);
	       psi_f1 = psi_p + SPINOR_SIZE * (posm_xyzt + vol);
	     }

	     OneMinusPlusGammaTimesPsi(tmp,psi,sdag,mu,1);
	     if(GJP.Gparity()) OneMinusPlusGammaTimesPsi(tmp+SPINOR_SIZE,psi_f1,sdag,mu,1);

	     /* multiply by U_mu */
	     UmuTimes(temps[mu][1],u,tmp,mu,1);
	     if(GJP.Gparity()) UmuTimes(temps[mu][1]+SPINOR_SIZE,u_f1,tmp+SPINOR_SIZE,mu,1);

	     if(pos[mu] == 0){
	       if(!GJP.Gparity()){
		 getMinusData((IFloat *)fbuf, (IFloat *)temps[mu][1], SPINOR_SIZE, mu);
		 moveMem((IFloat *)temps[mu][1], (IFloat *)fbuf, 
			 SPINOR_SIZE*sizeof(Float) / sizeof(char));
	       }else{
		 IFloat* psi_f0_recv_to = (IFloat*)fbuf;
		 IFloat* psi_f1_recv_to = (IFloat*)(fbuf + flav_fbuf_offset);
		 if(GJP.Bc(mu) == BND_CND_GPARITY && GJP.NodeCoor(mu)==0){
		   //do the G-parity flavour twist. (assumes minus sign is on the gauge links on the boundary)
		   psi_f0_recv_to = (IFloat*)(fbuf + flav_fbuf_offset);
		   psi_f1_recv_to = (IFloat*)fbuf;
		 }
		 getMinusData((IFloat *)psi_f0_recv_to, (IFloat *)temps[mu][1], SPINOR_SIZE, mu);
		 getMinusData((IFloat *)psi_f1_recv_to, (IFloat *)temps[mu][1]+SPINOR_SIZE, SPINOR_SIZE, mu);

		 moveMem((IFloat *)temps[mu][1], (IFloat *)fbuf, 
			 SPINOR_SIZE*sizeof(Float) / sizeof(char));
		 moveMem((IFloat *)(temps[mu][1]+SPINOR_SIZE), (IFloat *)(fbuf+ flav_fbuf_offset), 
			 SPINOR_SIZE*sizeof(Float) / sizeof(char));

	       }
	     }
	   }


	   /* Add all contributions */
	   //DEBUG
	   // if(GJP.Gparity() && x==0 && y==0 && z==1 && t==1){
	   //   printf("Temp contributions to chi0 (%d %d %d %d)\n",x,y,z,t);
	   //   for(int mu=0;mu<4;mu++){
	   //     for(int mp=0; mp<2; mp++){
	   // 	 Float sum(0.0);
	   // 	 for(int i=0;i<SPINOR_SIZE;i++){
	   // 	   sum += temps[mu][mp][i];
	   // 	 }
	   // 	 printf("tmp[%d][%d] = %f\n",mu,mp,sum);
	   //     }
	   //   }
	   //   printf("Temp contributions to chi1 (%d %d %d %d)\n",x,y,z,t);
	   //   for(int mu=0;mu<4;mu++){
	   //     for(int mp=0; mp<2; mp++){
	   // 	 Float sum(0.0);
	   // 	 for(int i=0;i<SPINOR_SIZE;i++){
	   // 	   sum += temps[mu][mp][i+SPINOR_SIZE];
	   // 	 }
	   // 	 printf("tmp[%d][%d] = %f\n",mu,mp,sum);
	   //     }
	   //   }
	   // }else if(!GJP.Gparity() && (x==0 || x==2) && y==0 && z==1 && t==1){
	   //   if(x<lx/2){
	   //     printf("Temp contributions to chi0 (%d %d %d %d)\n",x,y,z,t);
	   //   }else{
	   //     printf("Temp contributions to chi1 (%d %d %d %d)\n",x-lx/2,y,z,t);
	   //   }

	   //   for(int mu=0;mu<4;mu++){
	   //     for(int mp=0; mp<2; mp++){
	   // 	 Float sum(0.0);
	   // 	 for(int i=0;i<SPINOR_SIZE;i++){
	   // 	   sum += temps[mu][mp][i];
	   // 	 }
	   // 	 printf("tmp[%d][%d] = %f\n",mu,mp,sum);
	   //     }
	   //   }
	   // }
	   //DEBUG
	 

	   chi = chi_p + SPINOR_SIZE * xyzt;
	   for(s=0;s<4;s++){
	     for(c=0;c<3;c++){
	       for(r=0;r<2;r++){
		 FERM(chi,r,c,s) = (  FERM(tmp1,r,c,s)
				      + FERM(tmp2,r,c,s)
				      + FERM(tmp3,r,c,s)
				      + FERM(tmp4,r,c,s)
				      + FERM(tmp5,r,c,s)
				      + FERM(tmp6,r,c,s)
				      + FERM(tmp7,r,c,s)
				      + FERM(tmp8,r,c,s) );
	       }
	     }
	   }
	   

	   if(GJP.Gparity()){
	     chi_f1 = chi_p + SPINOR_SIZE * (xyzt + vol);
	     for(s=0;s<4;s++){
	       for(c=0;c<3;c++){
		 for(r=0;r<2;r++){
		   FERM(chi_f1,r,c,s) = (  FERM(tmp1+SPINOR_SIZE,r,c,s)
					   + FERM(tmp2+SPINOR_SIZE,r,c,s)
					   + FERM(tmp3+SPINOR_SIZE,r,c,s)
					   + FERM(tmp4+SPINOR_SIZE,r,c,s)
					   + FERM(tmp5+SPINOR_SIZE,r,c,s)
					   + FERM(tmp6+SPINOR_SIZE,r,c,s)
					   + FERM(tmp7+SPINOR_SIZE,r,c,s)
					   + FERM(tmp8+SPINOR_SIZE,r,c,s) );
		 }
	       }
	     }	     
	   }
	   


	 }
       }
      }
    }
  }

  //DEBUG
  // for(x=0; x<lx; x++){
  //   for(y=0; y<ly; y++){
  //     for(z=0; z<lz; z++){
  // 	for(t=0; t<lt; t++){
  // 	  xyzt = (x/2)+(lx/2)*(y+ly*(z+lz*t));
  // 	  chi = chi_p + SPINOR_SIZE * xyzt;
  // 	  if(GJP.Gparity()){
  // 	    printf("x=(%d,%d,%d,%d)  chi0_x = (%f,%f)\n",x,y,z,t,*chi,*(chi+1));
  // 	  }else{
  // 	    int gpx = x; 
  // 	    int gpf = 0;
  // 	    if(gpx>=lx/2){ gpx-=lx/2; gpf = 1; }
  // 	    printf("x=(%d,%d,%d,%d)  chi%d_x = (%f,%f)\n",gpx,y,z,t,gpf,*chi,*(chi+1));
  // 	    fflush(stdout);
  // 	  }
  // 	}
  //     }
  //   }
  // }
  // if(GJP.Gparity()){
  //   for(x=0; x<lx; x++){
  //     for(y=0; y<ly; y++){
  // 	for(z=0; z<lz; z++){
  // 	  for(t=0; t<lt; t++){
  // 	    xyzt = (x/2)+(lx/2)*(y+ly*(z+lz*t));
  // 	    chi_f1 = chi_p + SPINOR_SIZE * (xyzt + vol);
  // 	    printf("x=(%d,%d,%d,%d)  chi1_x = (%f,%f)\n",x,y,z,t,*chi_f1,*(chi_f1+1));
  // 	  }
  // 	}
  //     }
  //   }
  // }
  //DEBUG

//not correct because of non spin projection, but relevant for comparison
	DiracOp::CGflops += 1320*vol;


CPS_END_NAMESPACE
#else //USE_QMP
#ifdef USE_SSE
#include "sse-wilson_dslash.C"
#else
#include "../qmp/wilson_dslash_vec.C"
CPS_START_NAMESPACE
void wilson_dslash(IFloat *chi_p_f, 
			IFloat *u_p_f, 
			IFloat *psi_p_f, 
			int cb,
			int dag,
			Wilson *wilson_p)
{
	wilson_dslash_vec(chi_p_f,u_p_f,psi_p_f,cb,dag,wilson_p,1,0);
}
CPS_END_NAMESPACE
#endif //USE_SSE
#endif //USE_QMP

CPS_START_NAMESPACE
void wilson_dslash_two(Float *chi0, Float *chi1,
                   Float *u,
                   Float *psi0, Float *psi1,
                   int cb0, int cb1,
                   int dag,
                   Wilson *wp)
{
  wilson_dslash(chi0,u,psi0,cb0,dag,wp);
  wilson_dslash(chi1,u,psi1,cb1,dag,wp);
}

CPS_END_NAMESPACE
#endif
