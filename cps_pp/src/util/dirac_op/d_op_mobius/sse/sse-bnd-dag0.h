//#include "fake_omp.h"
void wilson_dslash_bnd_dag0(
                   IFloat *chi_p_f, 
		   IFloat *u_p_f, 
		   IFloat *psi_p_f, 
		   int cb,
		   Wilson *wilson_p)
{
  int lx, ly, lz, lt;
  int cbn;
  int vol;

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

  __m128d* send_buf0 = (__m128d*)(wilson_p->send_buf[0]);
  __m128d* send_buf1 = (__m128d*)(wilson_p->send_buf[1]);
  __m128d* send_buf2 = (__m128d*)(wilson_p->send_buf[2]);
  __m128d* send_buf3 = (__m128d*)(wilson_p->send_buf[3]);
  __m128d* send_buf4 = (__m128d*)(wilson_p->send_buf[4]); 
  __m128d* send_buf5 = (__m128d*)(wilson_p->send_buf[5]);
  __m128d* send_buf6 = (__m128d*)(wilson_p->send_buf[6]);
  __m128d* send_buf7 = (__m128d*)(wilson_p->send_buf[7]);

  // fixme: do it in wilson_p later
  const int block0=HALF_SPINOR_SIZE*ly*lz*lt/2;
  const int block1=HALF_SPINOR_SIZE*lx*lz*lt/2;
  const int block2=HALF_SPINOR_SIZE*lx*ly*lt/2;
  const int block3=HALF_SPINOR_SIZE*lx*ly*lz/2;

  if(cb == 0)  cbn = 1;
  else   cbn = 0;

  Float *u_p = (Float *) u_p_f;
  Float *psi_p = (Float *) psi_p_f;

  int tt;
  if(GJP.Xnodes() != 1) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1)
#endif
    for(tt = 0; tt < lt; tt++){

      int x, y, z, t, s;
      int i;
      int by, bz;
      int bt1, by1, bz1;
      int xp, yp, zp, tp;
      int xm, ym, zm, tm;
      int xyzt, _xyzt, __xyzt;
      int xpyzt, xypzt, xyzpt, xyztp;
      int xmyzt, xymzt, xyzmt, xyztm;
      int _xpyzt, _xypzt, _xyzpt, _xyztp;
      int _xmyzt, _xymzt, _xyzmt, _xyztm;
      Float __RESTRICT *u;
      Float __RESTRICT *chi;
      Float __RESTRICT *psi;
      int shft,_shft;      
      int ip,im,iu;
      ip=0; im=0;iu=0;
    

      x=lx-1; xp=0; xm=lx-1;
#ifdef _OPENMP
      if (omp_get_num_threads() == 4) {
	if (tt & 2)
	  t = (lt >> 1) + ((tt >> 2) << 1) + (tt & 1);
	else 
	  t = (lt >> 1) - ((tt >> 2) << 1) - (tt & 1) - 1;
      } else
	t = (tt + 1) % lt;
#else
	t = (tt + 1) % lt;
#endif
      //      tp = (t + 1)%lt;
      //tm = (t+lt-1)%lt; 
      for(z = 0; z < lz; z++){
	//    zp = (z + 1) % lz;
	//zm = (z - 1 +lz)%lz;
	_xpyzt = (xp>>1)+(lx>>1)*(ly*(z+lz*t));
	_xmyzt = (xm>>1)-(xp>>1)+(lx>>1);
	_shft=  (SPINOR_SIZE>>2)* ((ly*(z+lz*t))>>1);
	const int parity= (cbn^((x+z+t)&1)); // by is even
	for(y = parity; y < ly; y+=2){
	  {
	    register __m128d _a,_b,_c,_d;
	    
	    //		  xp = (x + 1) & ((x + 1 - lx) >> 31);
	    //		  xm = x - 1 + (((x - 1) >> 31) & lx);
	    xpyzt = _xpyzt  + (lx>>1)*y;
	    psi = psi_p + SPINOR_SIZE * xpyzt;
	    //shft=  (SPINOR_SIZE/4)* ((y+ly*(z+lz*t))/2);

	    shft=  (SPINOR_SIZE>>2)*(y>>1) + _shft;

	    
#if 0
	    *(send_buf0+shft)   = P_PROJ_X_03(0);
	    *(send_buf0+shft+1) = P_PROJ_X_03(1);
	    *(send_buf0+shft+2) = P_PROJ_X_03(2);
	    
	    *(send_buf0+shft+3) = P_PROJ_X_12(0);
	    *(send_buf0+shft+4) = P_PROJ_X_12(1);
	    *(send_buf0+shft+5) = P_PROJ_X_12(2);
#else
	    _mm_store_pd((double*)(send_buf0+shft), P_PROJ_X_03(0));
	    _mm_store_pd((double*)(send_buf0+shft+1) , P_PROJ_X_03(1));
	    _mm_store_pd((double*)(send_buf0+shft+2) , P_PROJ_X_03(2));
	    
	    _mm_store_pd((double*)(send_buf0+shft+3) , P_PROJ_X_12(0));
	    _mm_store_pd((double*)(send_buf0+shft+4) , P_PROJ_X_12(1));
	    _mm_store_pd((double*)(send_buf0+shft+5) , P_PROJ_X_12(2));
#endif
	    BND_PREFETCH_PSI;
	    
	    xmyzt = xpyzt + _xmyzt  - (lx>>1)*(parity<<1);
	    
	    psi = psi_p+ SPINOR_SIZE * xmyzt;
	    u   = u_p + GAUGE_SIZE * ( xmyzt  + vol * cb);

	    
	    //stay same: shft=  (SPINOR_SIZE/4)* ((yD+ly*(z+lz*t))/2);
#if 0
	    __m128d* wxm = send_buf4+shft;
	    P_KERN_XM_noadd;
	    BND_PREFETCH_U0; 
	    BND_PREFETCH_PSI;
#else
	    __m128d wxm[6];
	    P_KERN_XM_noadd;
	    BND_PREFETCH_U0; 
	    BND_PREFETCH_PSI;
	    
	    _mm_store_pd( (double*)(send_buf4+shft), wxm[0]);
	    _mm_store_pd( (double*)(send_buf4+shft+1), wxm[1]);
	    _mm_store_pd( (double*)(send_buf4+shft+2), wxm[2]);
	    _mm_store_pd( (double*)(send_buf4+shft+3), wxm[3]);
	    _mm_store_pd( (double*)(send_buf4+shft+4), wxm[4]);
	    _mm_store_pd( (double*)(send_buf4+shft+5), wxm[5]);
#endif
	  }//

	}
      }
    }
#ifdef BND_COMM
    //getPlusData((IFloat *)recv_buf1, (IFloat *)send_buf0, block0, 0);
    //getMinusData((IFloat *)recv_buf5, (IFloat *)send_buf4, block0, 0);
    QMP_start(wilson_p->multiple[0]);
    QMP_start(wilson_p->multiple[4]);
#endif
  }
    //----------------------------------------------------------//
  if(GJP.Ynodes() != 1) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1)
#endif
    for(tt = 0; tt < lt; tt++){

      int x, y, z, t, s;
      int i;
      int by, bz;
      int bt1, by1, bz1;
      int xp, yp, zp, tp;
      int xm, ym, zm, tm;
      int xyzt, _xyzt, __xyzt;
      int xpyzt, xypzt, xyzpt, xyztp;
      int xmyzt, xymzt, xyzmt, xyztm;
      int _xpyzt, _xypzt, _xyzpt, _xyztp;
      int _xmyzt, _xymzt, _xyzmt, _xyztm;
      Float __RESTRICT *u;
      Float __RESTRICT *chi;
      Float __RESTRICT *psi;
      int shft,_shft;      
      
      y=ly-1; yp=0; ym=ly-1;
#ifdef _OPENMP
      if (omp_get_num_threads() == 4) {
	if (tt & 2)
	  t = (lt >> 1) + ((tt >> 2) << 1) + (tt & 1);
	else 
	  t = (lt >> 1) - ((tt >> 2) << 1) - (tt & 1) - 1;
      } else
	t = (tt + 1) % lt;
#else
	t = (tt + 1) % lt;
#endif
      //tp = (t + 1)%lt;
      //tm = (t+lt-1)%lt; 
      for (bz = 0; bz < lz; bz += Z_BLOCK) {
	bz1 = bz + Z_BLOCK;
	if (bz1 >= lz) bz1 = lz;
	for(z = bz; z < bz1; z++){
	  //zp = (z + 1) % lz;
	  //zm = (z - 1 +lz)%lz;
	  int parity= (cbn^((y+z+t)&1)); 
	  _xypzt = (lx>>1)*(yp+ly*(z+lz*t));
	  _shft=  (SPINOR_SIZE>>2)* ((lx*(z+lz*t))>>1);
	  _xymzt =(lx>>1)*(ym+ly*(z+lz*t));
	  for(x = parity; x < lx; x+=2)
	    {
		  register __m128d _a,_b,_c,_d;

		  //		  xp = (x + 1) & ((x + 1 - lx) >> 31);
		  //		  xm = x - 1 + (((x - 1) >> 31) & lx);
 
		  xypzt = (x>>1) + _xypzt;
		  psi = psi_p + SPINOR_SIZE * xypzt;
		  //shft = (SPINOR_SIZE/4)* ((x+lx*(z+lz*t))/2);
		  shft= _shft +  (SPINOR_SIZE>>2)* (x>>1);
		  
		  *(send_buf1+shft)   = P_PROJ_Y_03(0);
		  *(send_buf1+shft+1) = P_PROJ_Y_03(1);
		  *(send_buf1+shft+2) = P_PROJ_Y_03(2);

		  *(send_buf1+shft+3) = P_PROJ_Y_12(0);
		  *(send_buf1+shft+4) = P_PROJ_Y_12(1);
		  *(send_buf1+shft+5) = P_PROJ_Y_12(2);
		  BND_PREFETCH_PSI;

		  //int xD = x+1-(parity<<1) ;// x+1 or x-1; 
		  //xymzt = (xD/2)+(lx/2)*(ym+ly*(z+lz*t));


		  xymzt = ((x+1-(parity<<1))>>1) + _xymzt;

		  psi = psi_p + SPINOR_SIZE * xymzt;
		  //shft=  (SPINOR_SIZE/4)* ((xD+lx*(z+lz*t))/2);
		  u   = u_p + GAUGE_SIZE * ( xymzt  + vol * cb);

		  __m128d* wym = send_buf5+shft;
		  P_KERN_YM_noadd;
		  BND_PREFETCH_U1;
		  BND_PREFETCH_PSI;


		  //		  printf("deb %d %d %d %d %d  %e %e %e\n",
		  //			 xD,ym,z,t,shft, *psi, *u, *(double*)&(wym[0]));

	  }//loop(x)
	}//loop(z)
      }//loop(bz)
    }//loop(t)
#ifdef BND_COMM
    //getPlusData((IFloat *)recv_buf2, (IFloat *)send_buf1, block1, 1);
    //getMinusData((IFloat *)recv_buf6, (IFloat *)send_buf5, block1, 1);
    QMP_start(wilson_p->multiple[1]);
    QMP_start(wilson_p->multiple[5]);
#endif
  }    
    //----------------------------------------------------------//

  if(GJP.Znodes() != 1) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1)
#endif
    for(tt = 0; tt < lt; tt++){

      int x, y, z, t, s;
      int i;
      int by, bz;
      int bt1, by1, bz1;
      int xp, yp, zp, tp;
      int xm, ym, zm, tm;
      int xyzt, _xyzt, __xyzt;
      int xpyzt, xypzt, xyzpt, xyztp;
      int xmyzt, xymzt, xyzmt, xyztm;
      int _xpyzt, _xypzt, _xyzpt, _xyztp;
      int _xmyzt, _xymzt, _xyzmt, _xyztm;
      Float __RESTRICT *u;
      Float __RESTRICT *chi;
      Float __RESTRICT *psi;
      int shft,_shft;      
      
      z=lz-1; zp=0; zm=lz-1;
#ifdef _OPENMP
      if (omp_get_num_threads() == 4) {
	if (tt & 2)
	  t = (lt >> 1) + ((tt >> 2) << 1) + (tt & 1);
	else 
	  t = (lt >> 1) - ((tt >> 2) << 1) - (tt & 1) - 1;
      } else
	t = (tt + 1) % lt;
#else
	t = (tt + 1) % lt;
#endif
      //      tp = (t + 1)%lt;
      //      tm = (t+lt-1)%lt; 

      for (by = 0; by < ly; by += Y_BLOCK) {
	by1 = by + Y_BLOCK;
	if (by1 >= ly) by1 = ly;
	for(y = by; y < by1; y++)   {

	  _xyzpt= (lx>>1)*(y+ly*(zp+lz*t));
	  _shft=(SPINOR_SIZE>>2) * ((lx*(y+ly*t))>>1);
	  _xyzmt=(lx>>1)*(y+ly*(zm+lz*t));
	  int parity= (cbn^((y+z+t)&1)); 
	  for(x = parity; x < lx; x+=2)
	    {
	      register __m128d _a,_b,_c,_d;
	      
	      //		  xp = (x + 1) & ((x + 1 - lx) >> 31);
	      //		  xm = x - 1 + (((x - 1) >> 31) & lx);
	      
	      //xyzpt = (x/2)+(lx/2)*(y+ly*(zp+lz*t));

	      xyzpt = (x>>1) + _xyzpt;
	      psi = psi_p + SPINOR_SIZE * xyzpt;
	      //shft=  (SPINOR_SIZE/4)* ((x+lx*(y+ly*t))/2);
	      shft=(SPINOR_SIZE>>2)*(x>>1)+_shft;

	      *(send_buf2+shft)   = P_PROJ_Z_02(0);
	      *(send_buf2+shft+1) = P_PROJ_Z_02(1);
	      *(send_buf2+shft+2) = P_PROJ_Z_02(2);
	      
	      *(send_buf2+shft+3) = P_PROJ_Z_13(0);
	      *(send_buf2+shft+4) = P_PROJ_Z_13(1);
	      *(send_buf2+shft+5) = P_PROJ_Z_13(2);
	      
	      //int xD = x+1-(parity<<1) ;// x+1 or x-1; 
	      //xyzmt = (xD/2)+(lx/2)*(y+ly*(zm+lz*t));

	      xyzmt = ((x+1-(parity<<1))>>1) + _xyzmt;

	      psi = psi_p + SPINOR_SIZE * xyzmt;
	      //shft=  (SPINOR_SIZE/4)* ((xD+lx*(y+ly*t))/2);
	      u   = u_p + GAUGE_SIZE * ( xyzmt  + vol * cb);
	      
	      __m128d* wzm = send_buf6+shft;
	      P_KERN_ZM_noadd;

		  
	      //printf("deb %d %d %d %d %d  %e %e %e\n",
	      //			 xD,ym,z,t,shft, *psi, *u, *(double*)&(wym[0]));
	      
	  }//loop(x)
	}//loop(z)
      }//loop(bz)
    }//loop(t)
#ifdef BND_COMM
    //getPlusData((IFloat *)recv_buf3, (IFloat *)send_buf2, block2, 2);
    //getMinusData((IFloat *)recv_buf7, (IFloat *)send_buf6, block2, 2);
    QMP_start(wilson_p->multiple[2]);
    QMP_start(wilson_p->multiple[6]);
#endif
  }    
    //----------------------------------------------------------//
  int zz;

  if(GJP.Tnodes() != 1) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1)
#endif
      for(zz = 0; zz < lz; zz++){
      int x, y,  t, z;
      int i;
      int by, bz;
      int bt1, by1, bz1;
      int xp, yp, zp, tp;
      int xm, ym, zm, tm;
      int xyzt, _xyzt, __xyzt;
      int xpyzt, xypzt, xyzpt, xyztp;
      int xmyzt, xymzt, xyzmt, xyztm;
      int _xpyzt, _xypzt, _xyzpt, _xyztp;
      int _xmyzt, _xymzt, _xyzmt, _xyztm;
      Float __RESTRICT *u;
      Float __RESTRICT *chi;
      Float __RESTRICT *psi;
      int shft,_shft;      

#ifdef _OPENMP
      if (omp_get_num_threads() == 4) {
	if (zz & 2)
	  z = (lz >> 1) + ((zz >> 2) << 1) + (zz & 1);
	else 
	  z = (lz >> 1) - ((zz >> 2) << 1) - (zz & 1) - 1;
      } else
	z = (zz + 1) % lz;
#else
	z = (zz + 1) % lz;
#endif

      t=lt-1; tp=0; tm=lt-1;
      for (by = 0; by < ly; by += Y_BLOCK) {
	by1 = by + Y_BLOCK;
	if (by1 >= ly) by1 = ly;
	for(y = by; y < by1; y++)   {
	  _xyztp=(lx>>1)*(y+ly*(z+lz*tp));
	  _xyztm=(lx/2)*(y+ly*(z+lz*tm));
	  _shft=  (SPINOR_SIZE>>2)* ((lx*(y+ly*z))>>1);
	  int parity= (cbn^((y+z+t)&1)); 
	  for(x = parity; x < lx; x+=2)
	    {
		  register __m128d _a,_b,_c,_d;

		  //		  xp = (x + 1) & ((x + 1 - lx) >> 31);
		  //		  xm = x - 1 + (((x - 1) >> 31) & lx);
 
		  xyztp = (x>>1)+_xyztp;
		  
		  psi = psi_p + SPINOR_SIZE * xyztp;
		  //shft=  (SPINOR_SIZE/4)* ((x+lx*(y+ly*z))/2);
		  shft=  (SPINOR_SIZE>>2)* (x>>1)+_shft;
		  

		  *(send_buf3+shft)   = P_PROJ_T_02(0);
		  *(send_buf3+shft+1) = P_PROJ_T_02(1);
		  *(send_buf3+shft+2) = P_PROJ_T_02(2);

		  *(send_buf3+shft+3) = P_PROJ_T_13(0);
		  *(send_buf3+shft+4) = P_PROJ_T_13(1);
		  *(send_buf3+shft+5) = P_PROJ_T_13(2);
		  BND_PREFETCH_PSI;

		  //printf("deb %d %d %d %d %d  %e %e %e\n",
		  //x,y,z,t,shft, *psi, *u, *(double*)(send_buf3+shft));

		  //int xD = x+1-(parity<<1) ;// x+1 or x-1; 
		  //xyztm = (xD/2)+(lx/2)*(y+ly*(z+lz*tm));

		  //int xD = x+1-(parity<<1) ;// x+1 or x-1; 
		  xyztm = ((x+1-(parity<<1))>>1)+_xyztm;

		  psi = psi_p + SPINOR_SIZE * xyztm;
		  //shft=  (SPINOR_SIZE/4)* ((xD+lx*(y+ly*z))/2);
		  u   = u_p + GAUGE_SIZE * ( xyztm  + vol * cb);

		  __m128d* wtm = send_buf7+shft;
		  P_KERN_TM_noadd;
		  BND_PREFETCH_U2;
		  BND_PREFETCH_PSI;



		  //		  printf("deb %d %d %d %d %d  %e %e %e\n",
		  //			 xD,y,z,t,shft, *psi, *u, *(double*)&(wtm[0]));

	    }//loop(x)
	}//loop(y)
	}//loop(z)
      }//loop(zz)

#ifdef BND_COMM
      //getPlusData((IFloat *)recv_buf4, (IFloat *)send_buf3, block3, 3);
      //getMinusData((IFloat *)recv_buf8, (IFloat *)send_buf7, block3, 3);
      QMP_start(wilson_p->multiple[3]);
      QMP_start(wilson_p->multiple[7]);
#endif
     
  }    
    //----------------------------------------------------------//

    //////////////////////////////////////////////////////////////////
  
#if 0
    /*  mpi comminucation start 2006.8.3 S.AOKI */
{


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

  


  int dag=0;
  int sdag;
  
  if(dag == 0)
    sdag = 1;
  else if(dag == 1)
    sdag = -1;
  else{
  }
  
		
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

#if 0
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
    //    getPlusData((IFloat *)recv_buf1, (IFloat *)send_buf0, block0, 0);

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
	       *(send_buf0+shft) = ( U(0,0,c,mu) * TMP(0,0,s)
				   + U(0,1,c,mu) * TMP(0,1,s)
				   + U(0,2,c,mu) * TMP(0,2,s) 
			           + U(1,0,c,mu) * TMP(1,0,s)
			           + U(1,1,c,mu) * TMP(1,1,s)
			           + U(1,2,c,mu) * TMP(1,2,s) );
	       ++shft;
	       *(send_buf0+shft) = ( U(0,0,c,mu) * TMP(1,0,s)
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
    //	     getMinusData((IFloat *)recv_buf5, (IFloat *)send_buf0, block0, 0);
	     //	     free(send_buf0);
  }
#endif
#define CHECK(A,B)					\
  {							\
    double reldiff = fabs( ((A)-(B))/((A)+(B)) ) *0.5;	\
    if(  reldiff > 1e-12 )				\
      {							\
	int gx=GJP.XnodeCoor()*lx+x;			\
	int gy=GJP.YnodeCoor()*ly+y;			\
	int gz=GJP.ZnodeCoor()*lz+z;			\
	int gt=GJP.TnodeCoor()*lt+t;			\
	printf("check failed (%d,%d,%d,%d): %e %e\n",	\
	       gx,gy,gz,gt, A,B);			\
      }}						\


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
#if 0
	     *(send_buf1+shft)   = PSI(0,c,0) - sdag * ( -PSI(0,c,3) ); 
	     *(send_buf1+shft+1) = PSI(1,c,0) - sdag * ( -PSI(1,c,3) ); 
	     *(send_buf1+shft6)   = PSI(0,c,1) - sdag * (  PSI(0,c,2) ); 
	     *(send_buf1+shft6+1) = PSI(1,c,1) - sdag * (  PSI(1,c,2) ); 
#else
	     CHECK(*(send_buf1+shft)  , PSI(0,c,0) - sdag * ( -PSI(0,c,3) )); 
	     CHECK(*(send_buf1+shft+1) , PSI(1,c,0) - sdag * ( -PSI(1,c,3) )); 
	     CHECK(*(send_buf1+shft6)   , PSI(0,c,1) - sdag * (  PSI(0,c,2) )); 
	     CHECK(*(send_buf1+shft6+1) , PSI(1,c,1) - sdag * (  PSI(1,c,2) )); 
#endif

	     shft+=2;
	     shft6+=2;
	   }
	   ++iyp;
	 }
	}
      }
    }
    //    getPlusData((IFloat *)recv_buf2, (IFloat *)send_buf1, block1, 1);
#define CHECK2(sft, A,B)				\
  if( (A) != (B) )				\
  {						\
    int gx=GJP.XnodeCoor()*lx+x;		\
    int gy=GJP.YnodeCoor()*ly+y;		\
    int gz=GJP.ZnodeCoor()*lz+z;		\
    int gt=GJP.TnodeCoor()*lt+t;			\
    printf("check failed (%d,%d,%d,%d) %d: %e %e\n",	\
  	   gx,gy,gz,gt, sft, A,B);				\
  }							\

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
	   //	   printf("PSI %d %d %d %d %e %e\n", x,y,z,t,PSI(0,0,0),*u);
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
#if 1
	       *(send_buf1+shft) = ( U(0,0,c,mu) * TMP(0,0,s)
				   + U(0,1,c,mu) * TMP(0,1,s)
				   + U(0,2,c,mu) * TMP(0,2,s) 
			           + U(1,0,c,mu) * TMP(1,0,s)
			           + U(1,1,c,mu) * TMP(1,1,s)
			           + U(1,2,c,mu) * TMP(1,2,s) );
	       ++shft;
	       *(send_buf1+shft) = ( U(0,0,c,mu) * TMP(1,0,s)
				   + U(0,1,c,mu) * TMP(1,1,s)
			           + U(0,2,c,mu) * TMP(1,2,s) 
			           - U(1,0,c,mu) * TMP(0,0,s)
			           - U(1,1,c,mu) * TMP(0,1,s)
			           - U(1,2,c,mu) * TMP(0,2,s) );
	       ++shft;
#else
	       CHECK2(shft,*(send_buf5+shft) , ( U(0,0,c,mu) * TMP(0,0,s)
				   + U(0,1,c,mu) * TMP(0,1,s)
				   + U(0,2,c,mu) * TMP(0,2,s) 
			           + U(1,0,c,mu) * TMP(1,0,s)
			           + U(1,1,c,mu) * TMP(1,1,s)
					   + U(1,2,c,mu) * TMP(1,2,s) ));
	       ++shft;
	       CHECK2(shft,*(send_buf5+shft) , ( U(0,0,c,mu) * TMP(1,0,s)
				   + U(0,1,c,mu) * TMP(1,1,s)
			           + U(0,2,c,mu) * TMP(1,2,s) 
			           - U(1,0,c,mu) * TMP(0,0,s)
			           - U(1,1,c,mu) * TMP(0,1,s)
					   - U(1,2,c,mu) * TMP(0,2,s) ));
	       ++shft;
#endif
	     }
	   }
	   ++iym;
	 }
	}
      }
    }
    //    	     getMinusData((IFloat *)recv_buf6, (IFloat *)send_buf1, block1, 1);
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
#if 0
	     *(send_buf2+shft)   = PSI(0,c,0) - sdag * ( -PSI(1,c,2) ); 
	     *(send_buf2+shft+1) = PSI(1,c,0) - sdag * (  PSI(0,c,2) ); 

	     *(send_buf2+shft6)   = PSI(0,c,1) - sdag * (  PSI(1,c,3) ); 
	     *(send_buf2+shft6+1) = PSI(1,c,1) - sdag * ( -PSI(0,c,3) ); 
#else
	     CHECK(*(send_buf2+shft)   , PSI(0,c,0) - sdag * ( -PSI(1,c,2) )); 
	     CHECK(*(send_buf2+shft+1) , PSI(1,c,0) - sdag * (  PSI(0,c,2) )); 

	     CHECK(*(send_buf2+shft6)   , PSI(0,c,1) - sdag * (  PSI(1,c,3) )); 
	     CHECK(*(send_buf2+shft6+1) , PSI(1,c,1) - sdag * ( -PSI(0,c,3) )); 

#endif
	   shft+=2;
	   shft6+=2;
	   }
	   ++izp;
	 }
	}
      }
    }
    //    getPlusData((IFloat *)recv_buf3, (IFloat *)send_buf2, block2, 2);

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
#if 0
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
#else
	       CHECK(*(send_buf6+shft) , ( U(0,0,c,mu) * TMP(0,0,s)
				   + U(0,1,c,mu) * TMP(0,1,s)
				   + U(0,2,c,mu) * TMP(0,2,s) 
			           + U(1,0,c,mu) * TMP(1,0,s)
			           + U(1,1,c,mu) * TMP(1,1,s)
					   + U(1,2,c,mu) * TMP(1,2,s) ));
	       ++shft;
	       CHECK(*(send_buf6+shft) , ( U(0,0,c,mu) * TMP(1,0,s)
				   + U(0,1,c,mu) * TMP(1,1,s)
			           + U(0,2,c,mu) * TMP(1,2,s) 
			           - U(1,0,c,mu) * TMP(0,0,s)
			           - U(1,1,c,mu) * TMP(0,1,s)
					   - U(1,2,c,mu) * TMP(0,2,s) ));
	       ++shft;
#endif
	     }
	   }
	   ++izm;
	 }
	}
      }
    }
    //    getMinusData((IFloat *)recv_buf7, (IFloat *)send_buf6, block2, 2);
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
	   //	   printf("PSI %d %d %d %d %e %e\n", x,y,z,t,PSI(0,0,0),*u);
       	   for(c=0;c<3;c++){
#if 1
	     *(send_buf3+shft)   = PSI(0,c,0) - sdag * (  PSI(0,c,2) ); 
	     *(send_buf3+shft+1) = PSI(1,c,0) - sdag * (  PSI(1,c,2) ); 

	     *(send_buf3+shft6)   = PSI(0,c,1) - sdag * (  PSI(0,c,3) ); 
	     *(send_buf3+shft6+1) = PSI(1,c,1) - sdag * (  PSI(1,c,3) ); 
#else
	     CHECK(*(send_buf3+shft)   , PSI(0,c,0) - sdag * (  PSI(0,c,2) )); 
	     CHECK(*(send_buf3+shft+1) , PSI(1,c,0) - sdag * (  PSI(1,c,2) )); 

	     CHECK(*(send_buf3+shft6)   , PSI(0,c,1) - sdag * (  PSI(0,c,3) )); 
	     CHECK(*(send_buf3+shft6+1) , PSI(1,c,1) - sdag * (  PSI(1,c,3) )); 
#endif
	   shft+=2;
	   shft6+=2;
	   }
	   ++itp;
	 }
	}
      }
    }
    //    getPlusData((IFloat *)recv_buf4, (IFloat *)send_buf3, block3, 3);

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
	   //	   printf("PSI %d %d %d %d %e %e\n", x,y,z,t,PSI(0,0,0),*u);
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
#if 0
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
#else
	       CHECK(*(send_buf7+shft) , ( U(0,0,c,mu) * TMP(0,0,s)
				   + U(0,1,c,mu) * TMP(0,1,s)
				   + U(0,2,c,mu) * TMP(0,2,s) 
			           + U(1,0,c,mu) * TMP(1,0,s)
			           + U(1,1,c,mu) * TMP(1,1,s)
					   + U(1,2,c,mu) * TMP(1,2,s) ));
	       ++shft;
	       CHECK(*(send_buf7+shft) , ( U(0,0,c,mu) * TMP(1,0,s)
				   + U(0,1,c,mu) * TMP(1,1,s)
			           + U(0,2,c,mu) * TMP(1,2,s) 
			           - U(1,0,c,mu) * TMP(0,0,s)
			           - U(1,1,c,mu) * TMP(0,1,s)
					   - U(1,2,c,mu) * TMP(0,2,s) ));
	       ++shft;
#endif
	     }
	   }
	   ++itm;
	 }
	}
      }
    }
    //       getMinusData((IFloat *)recv_buf8, (IFloat *)send_buf7, block3, 3);
	     //	     free(send_buf3);
  }
 }
#endif
}


