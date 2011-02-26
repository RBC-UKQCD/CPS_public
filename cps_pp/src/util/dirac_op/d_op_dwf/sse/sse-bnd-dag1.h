#include "fake_omp.h"
void wilson_dslash_bnd_dag1(
                   IFloat *chi_p_f, 
		   IFloat *u_p_f, 
		   IFloat *psi_p_f, 
		   int cb,
		   Wilson *wilson_p)
{

  const int dag=1;
  
  int lx, ly, lz, lt;
  int sdag;                 /* = +/-1 if dag = 0/1 */
  int cbn;
  int vol;

  char* fname ="";
  
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


  lx = wilson_p->ptr[0];
  ly = wilson_p->ptr[1];
  lz = wilson_p->ptr[2];
  lt = wilson_p->ptr[3];
  vol = wilson_p->vol[0];
  

  Float* const recv_buf1 = wilson_p->recv_buf[0];
  Float* const recv_buf2 = wilson_p->recv_buf[1];
  Float* const recv_buf3 = wilson_p->recv_buf[2];
  Float* const recv_buf4 = wilson_p->recv_buf[3];
  Float* const recv_buf5 = wilson_p->recv_buf[4];
  Float* const recv_buf6 = wilson_p->recv_buf[5];
  Float* const recv_buf7 = wilson_p->recv_buf[6];
  Float* const recv_buf8 = wilson_p->recv_buf[7];




  
  const int block0=HALF_SPINOR_SIZE*ly*lz*lt/2;
  const int block1=HALF_SPINOR_SIZE*lx*lz*lt/2;
  const int block2=HALF_SPINOR_SIZE*lx*ly*lt/2;
  const int block3=HALF_SPINOR_SIZE*lx*ly*lz/2;

/*  mpi comminucation start 2006.8.3 S.AOKI */
{
int shft; 
int shft6;
int ixp=0;	
int ixm=0;
int iyp=0;	
int iym=0;
int izp=0;	
int izm=0;
int itp=0;	
int itm=0;
Float* const send_buf0 = wilson_p->send_buf[0];
Float* const send_buf1 = wilson_p->send_buf[1];
Float* const send_buf2 = wilson_p->send_buf[2];
Float* const send_buf3 = wilson_p->send_buf[3];
Float* const send_buf4 = wilson_p->send_buf[4];
Float* const send_buf5 = wilson_p->send_buf[5];
Float* const send_buf6 = wilson_p->send_buf[6];
Float* const send_buf7 = wilson_p->send_buf[7];

  int x, y, z, t;
  int xp, yp, zp, tp;
  int xm, ym, zm, tm;
  int xyzt;
  int xpyzt, xypzt, xyzpt, xyztp;
  int xmyzt, xymzt, xyzmt, xyztm;
  int parity;
  int r, c, s, mu;

  Float *chi_p = (Float *) chi_p_f;
  Float *u_p = (Float *) u_p_f;
  Float *psi_p = (Float *) psi_p_f;
  Float *chi;
  Float *u;
  Float *psi;
  //mu=0
  //SA 2007.5.22
  if(GJP.Xnodes() != 1) {
    //send_buf0=(Float *) malloc(block2*sizeof(Float));
    //recv_buf1=(Float *) malloc(block2*sizeof(Float));
    //recv_buf5=(Float *) malloc(block2*sizeof(Float));
    //plus direction
    x=lx-1; xp=0;
    for(t=0; t<lt; t++){
      for(z=0; z<lz; z++){
	for(y=0; y<ly; y++){
	 parity = x+y+z+t;
	 parity = parity % 2;
	 if(parity == cbn){
	   xpyzt = (xp/2)+(lx/2)*(y+ly*(z+lz*t));
	   psi = psi_p + SPINOR_SIZE * xpyzt;
	   shft=ixp*HALF_SPINOR_SIZE;
	   shft6=shft+6;
	   for(c=0;c<3;c++){
	     *(send_buf0+shft)   = PSI(0,c,0) - sdag * ( -PSI(1,c,3) ); 
	     *(send_buf0+shft+1) = PSI(1,c,0) - sdag * (  PSI(0,c,3) ); 

	     //	     int my_xyzt=x+lx*(y+ly*(z+lz*t));
	     //	     ::printf("send_buf (%d,%d,%d,%d) %d %d %d %e %e\n",lx*GJP.XnodeCoor()+x,y,z,t, my_xyzt, my_xyzt/lx/2, ixp,*(send_buf0+shft), *(send_buf0+shft+1));
   
	     *(send_buf0+shft6)   = PSI(0,c,1) - sdag * ( -PSI(1,c,2) ); 
	     *(send_buf0+shft6+1) = PSI(1,c,1) - sdag * (  PSI(0,c,2) ); 
	   shft+=2;
	   shft6+=2;
	   }
	   ++ixp;
	 }
	}
      }
    }
    //getPlusData((IFloat *)recv_buf1, (IFloat *)send_buf0, block0, 0);
    QMP_start(wilson_p->multiple[0]);    

    //minus direction
    x=0; xm=lx-1;
    for(t=0; t<lt; t++){
      for(z=0; z<lz; z++){
	for(y=0; y<ly; y++){
	 parity = x+y+z+t;
	 parity = parity % 2;
	 if(parity == cbn){
	   xmyzt = (xm/2)+(lx/2)*(y+ly*(z+lz*t));
	   u   = u_p + GAUGE_SIZE * ( xmyzt  + vol * cb);
	   psi = psi_p + SPINOR_SIZE * xmyzt;
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = PSI(0,c,0) + sdag * ( -PSI(1,c,3) ); 
	     TMP(1,c,0) = PSI(1,c,0) + sdag * (  PSI(0,c,3) ); 

	     TMP(0,c,1) = PSI(0,c,1) + sdag * ( -PSI(1,c,2) ); 
	     TMP(1,c,1) = PSI(1,c,1) + sdag * (  PSI(0,c,2) ); 
	   }
	   /* multiply by U_mu */
	   mu = 0;
	   shft=ixm*HALF_SPINOR_SIZE;
	   for(s=0;s<2;s++){
	     for(c=0;c<3;c++){
	       *(send_buf4+shft) = ( U(0,0,c,mu) * TMP(0,0,s)
				   + U(0,1,c,mu) * TMP(0,1,s)
				   + U(0,2,c,mu) * TMP(0,2,s) 
			           + U(1,0,c,mu) * TMP(1,0,s)
			           + U(1,1,c,mu) * TMP(1,1,s)
			           + U(1,2,c,mu) * TMP(1,2,s) );
	       ++shft;
	       *(send_buf4+shft) = ( U(0,0,c,mu) * TMP(1,0,s)
				   + U(0,1,c,mu) * TMP(1,1,s)
			           + U(0,2,c,mu) * TMP(1,2,s) 
			           - U(1,0,c,mu) * TMP(0,0,s)
			           - U(1,1,c,mu) * TMP(0,1,s)
			           - U(1,2,c,mu) * TMP(0,2,s) );
	       ++shft;
	     }
	   }
	   ++ixm;
	 }
	}
      }
    }
    //getMinusData((IFloat *)recv_buf5, (IFloat *)send_buf0, block0, 0);
    QMP_start(wilson_p->multiple[4]);    
	     //	     free(send_buf0);
  }
  //mu=1
  if(GJP.Ynodes() != 1) {
  //plus direction
    //send_buf1=(Float *) malloc(block2*sizeof(Float));
    //recv_buf2=(Float *) malloc(block2*sizeof(Float));
    //recv_buf6=(Float *) malloc(block2*sizeof(Float));
    y=ly-1; yp=0;
    for(t=0; t<lt; t++){
      for(z=0; z<lz; z++){
	for(x=0; x<lx; x++){
	 parity = x+y+z+t;
	 parity = parity % 2;
	 if(parity == cbn){
	   xypzt = (x/2)+(lx/2)*(yp+ly*(z+lz*t));
	   psi = psi_p + SPINOR_SIZE * xypzt;
	   shft=iyp*HALF_SPINOR_SIZE;
	   shft6=shft+6;
	   for(c=0;c<3;c++){
	     *(send_buf1+shft)   = PSI(0,c,0) - sdag * ( -PSI(0,c,3) ); 
	     *(send_buf1+shft+1) = PSI(1,c,0) - sdag * ( -PSI(1,c,3) ); 

	     *(send_buf1+shft6)   = PSI(0,c,1) - sdag * (  PSI(0,c,2) ); 
	     *(send_buf1+shft6+1) = PSI(1,c,1) - sdag * (  PSI(1,c,2) ); 
	   shft+=2;
	   shft6+=2;
	   }
	   ++iyp;
	 }
	}
      }
    }
    // getPlusData((IFloat *)recv_buf2, (IFloat *)send_buf1, block1, 1);
    QMP_start(wilson_p->multiple[1]);    
  //minus direction
    y=0; ym=ly-1;
    for(t=0; t<lt; t++){
      for(z=0; z<lz; z++){
	for(x=0; x<lx; x++){
	 parity = x+y+z+t;
	 parity = parity % 2;
	 if(parity == cbn){
	   xymzt = (x/2)+(lx/2)*(ym+ly*(z+lz*t));
	   u   = u_p + GAUGE_SIZE * ( xymzt  + vol * cb);
	   psi = psi_p + SPINOR_SIZE * xymzt;
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = PSI(0,c,0) + sdag * ( -PSI(0,c,3) ); 
	     TMP(1,c,0) = PSI(1,c,0) + sdag * ( -PSI(1,c,3) ); 

	     TMP(0,c,1) = PSI(0,c,1) + sdag * (  PSI(0,c,2) ); 
	     TMP(1,c,1) = PSI(1,c,1) + sdag * (  PSI(1,c,2) ); 
	   }
	   /* multiply by U_mu */
	   mu = 1;
	   shft=iym*HALF_SPINOR_SIZE;
	   for(s=0;s<2;s++){
	     for(c=0;c<3;c++){
	       *(send_buf5+shft) = ( U(0,0,c,mu) * TMP(0,0,s)
				   + U(0,1,c,mu) * TMP(0,1,s)
				   + U(0,2,c,mu) * TMP(0,2,s) 
			           + U(1,0,c,mu) * TMP(1,0,s)
			           + U(1,1,c,mu) * TMP(1,1,s)
			           + U(1,2,c,mu) * TMP(1,2,s) );
	       ++shft;
	       *(send_buf5+shft) = ( U(0,0,c,mu) * TMP(1,0,s)
				   + U(0,1,c,mu) * TMP(1,1,s)
			           + U(0,2,c,mu) * TMP(1,2,s) 
			           - U(1,0,c,mu) * TMP(0,0,s)
			           - U(1,1,c,mu) * TMP(0,1,s)
			           - U(1,2,c,mu) * TMP(0,2,s) );
	       ++shft;
	     }
	   }
	   ++iym;
	 }
	}
      }
    }
    //getMinusData((IFloat *)recv_buf6, (IFloat *)send_buf1, block1, 1);
    QMP_start(wilson_p->multiple[5]);    
    ///free(send_buf1);
  }
  //mu=2
  if(GJP.Znodes() != 1) {
    //send_buf2=(Float *) malloc(block2*sizeof(Float));
    //recv_buf3=(Float *) malloc(block2*sizeof(Float));
    //recv_buf7=(Float *) malloc(block2*sizeof(Float));
  //plus direction
    z=lz-1; zp=0;
    for(t=0; t<lt; t++){
      for(y=0; y<ly; y++){
	for(x=0; x<lx; x++){
	 parity = x+y+z+t;
	 parity = parity % 2;
	 if(parity == cbn){
	   xyzpt = (x/2)+(lx/2)*(y+ly*(zp+lz*t));
	   psi = psi_p + SPINOR_SIZE * xyzpt;
	   shft=izp*HALF_SPINOR_SIZE;
	   shft6=shft+6;
	   for(c=0;c<3;c++){
	     *(send_buf2+shft)   = PSI(0,c,0) - sdag * ( -PSI(1,c,2) ); 
	     *(send_buf2+shft+1) = PSI(1,c,0) - sdag * (  PSI(0,c,2) ); 

	     *(send_buf2+shft6)   = PSI(0,c,1) - sdag * (  PSI(1,c,3) ); 
	     *(send_buf2+shft6+1) = PSI(1,c,1) - sdag * ( -PSI(0,c,3) ); 
	   shft+=2;
	   shft6+=2;
	   }
	   ++izp;
	 }
	}
      }
    }
    //getPlusData((IFloat *)recv_buf3, (IFloat *)send_buf2, block2, 2);
    QMP_start(wilson_p->multiple[2]);    
  //minus direction
    z=0; zm=lz-1;
    for(t=0; t<lt; t++){
      for(y=0; y<ly; y++){
	for(x=0; x<lx; x++){
	 parity = x+y+z+t;
	 parity = parity % 2;
	 if(parity == cbn){
	   xyzmt = (x/2)+(lx/2)*(y+ly*(zm+lz*t));
	   u   = u_p + GAUGE_SIZE * ( xyzmt  + vol * cb);
	   psi = psi_p + SPINOR_SIZE * xyzmt;
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = PSI(0,c,0) + sdag * ( -PSI(1,c,2) ); 
	     TMP(1,c,0) = PSI(1,c,0) + sdag * (  PSI(0,c,2) ); 

	     TMP(0,c,1) = PSI(0,c,1) + sdag * (  PSI(1,c,3) ); 
	     TMP(1,c,1) = PSI(1,c,1) + sdag * ( -PSI(0,c,3) ); 
	   }
	   /* multiply by U_mu */
	   mu = 2;
	   shft=izm*HALF_SPINOR_SIZE;
	   for(s=0;s<2;s++){
	     for(c=0;c<3;c++){
	       *(send_buf6+shft) = ( U(0,0,c,mu) * TMP(0,0,s)
				   + U(0,1,c,mu) * TMP(0,1,s)
				   + U(0,2,c,mu) * TMP(0,2,s) 
			           + U(1,0,c,mu) * TMP(1,0,s)
			           + U(1,1,c,mu) * TMP(1,1,s)
			           + U(1,2,c,mu) * TMP(1,2,s) );
	       ++shft;
	       *(send_buf6+shft) = ( U(0,0,c,mu) * TMP(1,0,s)
				   + U(0,1,c,mu) * TMP(1,1,s)
			           + U(0,2,c,mu) * TMP(1,2,s) 
			           - U(1,0,c,mu) * TMP(0,0,s)
			           - U(1,1,c,mu) * TMP(0,1,s)
			           - U(1,2,c,mu) * TMP(0,2,s) );
	       ++shft;
	     }
	   }
	   ++izm;
	 }
	}
      }
    }
    //getMinusData((IFloat *)recv_buf7, (IFloat *)send_buf2, block2, 2);
    QMP_start(wilson_p->multiple[6]);    
    //free(send_buf2);
  }
  //mu=3
  if(GJP.Tnodes() != 1) {
  //send_buf3=(Float *) malloc(block3*sizeof(Float));
  //recv_buf4=(Float *) malloc(block3*sizeof(Float));
  //recv_buf8=(Float *) malloc(block3*sizeof(Float));
  //plus direction
    t=lt-1; tp=0;
    for(z=0; z<lz; z++){
      for(y=0; y<ly; y++){
	for(x=0; x<lx; x++){
	 parity = x+y+z+t;
	 parity = parity % 2;
	 if(parity == cbn){
	   xyztp = (x/2)+(lx/2)*(y+ly*(z+lz*tp));
	   psi = psi_p + SPINOR_SIZE * xyztp;
	   shft=itp*HALF_SPINOR_SIZE;
	   shft6=shft+6;
	   for(c=0;c<3;c++){
	     *(send_buf3+shft)   = PSI(0,c,0) - sdag * (  PSI(0,c,2) ); 
	     *(send_buf3+shft+1) = PSI(1,c,0) - sdag * (  PSI(1,c,2) ); 

	     *(send_buf3+shft6)   = PSI(0,c,1) - sdag * (  PSI(0,c,3) ); 
	     *(send_buf3+shft6+1) = PSI(1,c,1) - sdag * (  PSI(1,c,3) ); 
	   shft+=2;
	   shft6+=2;
	   }
	   ++itp;
	 }
	}
      }
    }
    //	     getPlusData((IFloat *)recv_buf4, (IFloat *)send_buf3, block3, 3);
    QMP_start(wilson_p->multiple[3]);    
	     
  //minus direction
    t=0; tm=lt-1;
    for(z=0; z<lz; z++){
      for(y=0; y<ly; y++){
	for(x=0; x<lx; x++){
	 parity = x+y+z+t;
	 parity = parity % 2;
	 if(parity == cbn){
	   xyztm = (x/2)+(lx/2)*(y+ly*(z+lz*tm));
	   u   = u_p + GAUGE_SIZE * ( xyztm  + vol * cb);
	   psi = psi_p + SPINOR_SIZE * xyztm;
	   for(c=0;c<3;c++){
	     TMP(0,c,0) = PSI(0,c,0) + sdag * (  PSI(0,c,2) ); 
	     TMP(1,c,0) = PSI(1,c,0) + sdag * (  PSI(1,c,2) ); 

	     TMP(0,c,1) = PSI(0,c,1) + sdag * (  PSI(0,c,3) ); 
	     TMP(1,c,1) = PSI(1,c,1) + sdag * (  PSI(1,c,3) ); 
	   }
	   /* multiply by U_mu */
	   mu = 3;
	   shft=itm*HALF_SPINOR_SIZE;
	   for(s=0;s<2;s++){
	     for(c=0;c<3;c++){
	       *(send_buf7+shft) = ( U(0,0,c,mu) * TMP(0,0,s)
				   + U(0,1,c,mu) * TMP(0,1,s)
				   + U(0,2,c,mu) * TMP(0,2,s) 
			           + U(1,0,c,mu) * TMP(1,0,s)
			           + U(1,1,c,mu) * TMP(1,1,s)
			           + U(1,2,c,mu) * TMP(1,2,s) );
	       ++shft;
	       *(send_buf7+shft) = ( U(0,0,c,mu) * TMP(1,0,s)
				   + U(0,1,c,mu) * TMP(1,1,s)
			           + U(0,2,c,mu) * TMP(1,2,s) 
			           - U(1,0,c,mu) * TMP(0,0,s)
			           - U(1,1,c,mu) * TMP(0,1,s)
			           - U(1,2,c,mu) * TMP(0,2,s) );
	       ++shft;
	     }
	   }
	   ++itm;
	 }
	}
      }
    }
    //getMinusData((IFloat *)recv_buf8, (IFloat *)send_buf3, block3, 3);
    QMP_start(wilson_p->multiple[7]);    
	     //	     free(send_buf3);
  }
}

}

  
