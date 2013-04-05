//#include "fake_omp.h"
void wilson_dslash_blk_dag1(IFloat *chi_p_f, 
			    IFloat *u_p_f, 
			    IFloat *psi_p_f, 
			    int cb,
			    Wilson *wilson_p)
{

  int tt;
  int cbn;
  
  
  Float *chi_p = (Float *) chi_p_f;
  Float *u_p = (Float *) u_p_f;
  Float *psi_p = (Float *) psi_p_f;

#if defined(_OPENMP) && defined(__linux__)
  static int init = 0;

  //  if (init == 0) { cpubind(); init = 1; }
#endif

  const int lx = wilson_p->ptr[0];
  const int ly = wilson_p->ptr[1];
  const int lz = wilson_p->ptr[2];
  const int lt = wilson_p->ptr[3];
  const int vol = wilson_p->vol[0];

  const int x_loc  =   GJP.Xnodes() == 1 ? 1 : 0;
  const int x_nloc =   !x_loc;
  const int y_nloc =  GJP.Ynodes() == 1 ? 0 : 1;
  const int z_nloc =  GJP.Znodes() == 1 ? 0 : 1;
  const int t_nloc =  GJP.Tnodes() == 1 ? 0 : 1;
    

  const Float* const recv_buf1 =  wilson_p->recv_buf[0];
  const Float* const recv_buf2 =  wilson_p->recv_buf[1];
  const Float* const recv_buf3 =  wilson_p->recv_buf[2];
  const Float* const recv_buf4 =  wilson_p->recv_buf[3];
  const Float* const recv_buf5 =  wilson_p->recv_buf[4];
  const Float* const recv_buf6 =  wilson_p->recv_buf[5];
  const Float* const recv_buf7 =  wilson_p->recv_buf[6];
  const Float* const recv_buf8 =  wilson_p->recv_buf[7];


  //for(int itr=0;itr<NITR;++itr){
  //printf("%d x %d x %d x %d = %d\n",lx,ly,lz,lt,vol);


  if(cb == 0)   cbn = 1;
  else    cbn = 0;


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

    //printf("tt=%d\n",tt);
    
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

    tp = (t + 1) % lt;
    tm = t - 1 + ((lt - t) / lt) * lt;

    
    for (bz = 0; bz < lz; bz += Z_BLOCK) {
      bz1 = bz + Z_BLOCK;

      if (bz1 >= lz) bz1 = lz;
      
      for (by = 0; by < ly; by += Y_BLOCK) {
	
	by1 = by + Y_BLOCK;
	if (by1 >= ly) by1 = ly;
	
	for(z = bz; z < bz1; z++){
	  
	  zp = (z + 1) % lz;
	  zm = z - 1 + ((lz - z) / lz) * lz;
	  
	  
	  _xyzt = (lx >> 1) * (y + ly * (z + lz * t));
	  _xpyzt = (lx >> 1) * (y + ly * (z + lz * t));
	  _xmyzt = (lx >> 1) * (y + ly * (z + lz * t));
	  _xypzt = (lx >> 1) * (yp + ly * (z + lz * t));
	  _xymzt = (lx >> 1) * (ym + ly * (z + lz * t));
	  _xyzpt = (lx >> 1) * (y + ly * (zp + lz * t));
	  _xyzmt = (lx >> 1) * (y + ly * (zm + lz * t));
	  _xyztp = (lx >> 1) * (y + ly * (z + lz * tp));
	  _xyztm = (lx >> 1) * (y + ly * (z + lz * tm));
	  
	  for(y = by; y < by1; y++){
	    int x_check, yzt_edge;
	    
	    yp = (y + 1) % ly;
	    ym = y - 1 + ((ly - y) / ly ) * ly;
	    
	    _xyzt = (lx >> 1) * (y + ly * (z + lz * t));
	    
	    _xpyzt = (lx >> 1) * (y + ly * (z + lz * t));
	    _xmyzt = (lx >> 1) * (y + ly * (z + lz * t));
	    _xypzt = (lx >> 1) * (yp + ly * (z + lz * t));
	    _xymzt = (lx >> 1) * (ym + ly * (z + lz * t));
	    _xyzpt = (lx >> 1) * (y + ly * (zp + lz * t));
	    _xyzmt = (lx >> 1) * (y + ly * (zm + lz * t));
	    _xyztp = (lx >> 1) * (y + ly * (z + lz * tp));
	    _xyztm = (lx >> 1) * (y + ly * (z + lz * tm));
	      
	    if ( ((z == 0 || z == lz - 1)&& z_nloc) ||
		 ((t == 0 || t == lt - 1)&& t_nloc) ||
		 ((y == 0 || y == ly - 1)&& y_nloc) ) yzt_edge = 0;
	    else yzt_edge = -1;
	
	    if ((y + z + t + cbn) & 1) x_check = lx - 1;
	    else x_check = 0;


	    for (x = cbn ^ ((y + z + t) & 1); x < lx; x += 2) 
	      {
		//::printf("%d %d %d %d\n",x,y,z,t);
		
		const int x_edge = x==x_check ? x_loc : 1;

		
		//		DECLARE;
		
#ifdef ADD2REG
		//		__m128d __RESTRICT wxp[6],wyp[6], wzp[6], wtp[6];	
		//		__m128d __RESTRICT wxm[6],wym[6], wzm[6], wtm[6];	
		register __m128d t00, t01, t02, t03, t04, t05;	
		register __m128d t06, t07, t08, t09, t10, t11; 
#else

		__m128d __RESTRICT wxp[6];
		__m128d __RESTRICT wyp[6];
		__m128d __RESTRICT wzp[6];
		__m128d __RESTRICT wtp[6];

		__m128d __RESTRICT wxm[6];
		__m128d __RESTRICT wym[6];
		__m128d __RESTRICT wzm[6];
		__m128d __RESTRICT wtm[6];
#endif

#ifndef USE_HERN
		register __m128d _a, _b, _c, _d;		
#endif

		xyzt = (x >> 1) + _xyzt;
	  
		xp = (x + 1) & ((x + 1 - lx) >> 31);
		xm = x - 1 + (((x - 1) >> 31) & lx);
	  
		xpyzt = (xp >> 1) + _xpyzt;
		xmyzt = (xm >> 1) + _xmyzt;
		xypzt = (x >> 1) + _xypzt;
		xymzt = (x >> 1) + _xymzt;
		xyzpt = (x >> 1) + _xyzpt;
		xyzmt = (x >> 1) + _xyzmt;
		xyztp = (x >> 1) + _xyztp;
		xyztm = (x >> 1) + _xyztm;
	  
#ifndef  DEBUG_NOEDGE
		if ((yzt_edge & x_edge ) == 0) {	    

		  ZERO;

		  u = u_p + GAUGE_SIZE * (xyzt + vol * cbn);
		  psi = psi_p + SPINOR_SIZE * xpyzt;
		  if(x_nloc&&(x == lx - 1)){
		    const size_t shft  = (SPINOR_SIZE/2)* ((y+ly*(z+lz*t))/2);
		    N_KERN_XP_EDG( recv_buf1+shft, recv_buf1+shft+6 );
		  }else  {N_KERN_XP;}
		  PREFETCH_U0;
		  PREFETCH_PSI;

		  u = u_p + GAUGE_SIZE * (xyzt + vol * cbn);
		  psi = psi_p + SPINOR_SIZE * xypzt;
		  if( y_nloc && (y == ly-1) ) {
		    const size_t shft  = (SPINOR_SIZE/2)* ((x+lx*(z+lz*t))/2);

#if 0
		    {//debug

		      int dag=1;
		      int gx=GJP.XnodeCoor()*lx+x;
		      int gy=GJP.YnodeCoor()*ly+y;
		      int gz=GJP.ZnodeCoor()*lz+z;
		      int gt=GJP.TnodeCoor()*lt+t;

		      if( gx==0 && gy==3 && gz==6 && gt==0){
			printf("debug4 %d %d (%d %d %d %d) a:",
			       cb,dag, gx,gy,gz,gt );

			for(int i=0;i<6;++i)
			  printf("%e ", *(i+recv_bf2+shft));
			printf(" | b: ");
			for(int i=0;i<6;++i)
			  printf("%e ", *(i+recv_bf2+6+shft));
		      }
		  }//debug
#endif

		    N_KERN_YP_EDG( recv_buf2+shft, recv_buf2+shft+6 );


#if 0
		    {//debug

		      int dag=1;
		      int gx=GJP.XnodeCoor()*lx+x;
		      int gy=GJP.YnodeCoor()*ly+y;
		      int gz=GJP.ZnodeCoor()*lz+z;
		      int gt=GJP.TnodeCoor()*lt+t;

		      if( gx==0 && gy==3 && gz==6 && gt==0){
			printf("debug4 %d %d (%d %d %d %d) wyp ",
			       cb,dag, gx,gy,gz,gt );

			for(int i=0;i<24;++i)
			  printf("%e ", *((double*)wyp+i));
		      }

		  }//debug
#endif

		  }else  {N_KERN_YP;}
		  PREFETCH_U1;
		  PREFETCH_PSI;

		  u = u_p + GAUGE_SIZE * (xyzt + vol * cbn);
		  psi = psi_p + SPINOR_SIZE * xyzpt;
		  if( z_nloc && (z == lz-1) ) {
		    const size_t shft  = (SPINOR_SIZE/2)* ((x+lx*(y+ly*t))/2);
		    N_KERN_ZP_EDG( recv_buf3+shft, recv_buf3+shft+6 );
		  }else  {N_KERN_ZP;}
		  PREFETCH_U2;
		  PREFETCH_PSI;

		  u = u_p + GAUGE_SIZE * (xyzt + vol * cbn);
		  psi = psi_p + SPINOR_SIZE * xyztp;
		  if(t_nloc && (t == lt - 1)){
		    const size_t shft  = (SPINOR_SIZE/2)* ((x+lx*(y+ly*z))/2);
		    N_KERN_TP_EDG( recv_buf4+shft, recv_buf4+shft+6 );
		  }else  {N_KERN_TP;}
		  PREFETCH_U3;
		  PREFETCH_PSI;

		  u = u_p + GAUGE_SIZE * (xmyzt + vol * cb);
		  psi = psi_p + SPINOR_SIZE * xmyzt;
		  if(x == 0 && x_nloc)  {
		    const size_t shft  = (SPINOR_SIZE/2)* ((y+ly*(z+lz*t))/2);
		    wxm[0] = _mm_load_pd( recv_buf5+shft );
		    wxm[1] = _mm_load_pd( recv_buf5+shft + 2);
		    wxm[2] = _mm_load_pd( recv_buf5+shft + 4);
		    wxm[3] = _mm_load_pd( recv_buf5+shft + 6);
		    wxm[4] = _mm_load_pd( recv_buf5+shft + 8);
		    wxm[5] = _mm_load_pd( recv_buf5+shft + 10);
		  } else  {N_KERN_XM;}
		  PREFETCH_U0;
		  PREFETCH_PSI;

		  u = u_p + GAUGE_SIZE * (xymzt + vol * cb);
		  psi = psi_p + SPINOR_SIZE * xymzt;

		  if( y_nloc && y==0)  {
		    const size_t shft  = (SPINOR_SIZE/2)* ((x+lx*(z+lz*t))/2);
		    wym[0] = _mm_load_pd( recv_buf6+shft );
		    wym[1] = _mm_load_pd( recv_buf6+shft + 2);
		    wym[2] = _mm_load_pd( recv_buf6+shft + 4);
		    wym[3] = _mm_load_pd( recv_buf6+shft + 6);
		    wym[4] = _mm_load_pd( recv_buf6+shft + 8);
		    wym[5] = _mm_load_pd( recv_buf6+shft + 10);

		  }
		  else  {N_KERN_YM;}
		  PREFETCH_U1;
		  PREFETCH_PSI;

		  u = u_p + GAUGE_SIZE * (xyzmt + vol * cb);
		  psi = psi_p + SPINOR_SIZE * xyzmt;
		  if ( z_nloc && z == 0 )  {
		    const size_t shft  = (SPINOR_SIZE/2)* ((x+lx*(y+ly*t))/2);

		      wzm[0] = _mm_load_pd( recv_buf7+shft );
		      wzm[1] = _mm_load_pd( recv_buf7+shft + 2);
		      wzm[2] = _mm_load_pd( recv_buf7+shft + 4);
		      wzm[3] = _mm_load_pd( recv_buf7+shft + 6);
		      wzm[4] = _mm_load_pd( recv_buf7+shft + 8);
		      wzm[5] = _mm_load_pd( recv_buf7+shft + 10);

		  } else    {N_KERN_ZM;}
		  PREFETCH_U2;
		  PREFETCH_PSI;

		  u = u_p + GAUGE_SIZE * (xyztm + vol * cb);
		  psi = psi_p + SPINOR_SIZE * xyztm;
		  if ( t_nloc && t == 0 ) {
	    const size_t shft  = (SPINOR_SIZE/2)* ((x+lx*(y+ly*z))/2);
		      wtm[0] = _mm_load_pd( recv_buf8+shft );
		      wtm[1] = _mm_load_pd( recv_buf8+shft + 2);
		      wtm[2] = _mm_load_pd( recv_buf8+shft + 4);
		      wtm[3] = _mm_load_pd( recv_buf8+shft + 6);
		      wtm[4] = _mm_load_pd( recv_buf8+shft + 8);
		      wtm[5] = _mm_load_pd( recv_buf8+shft + 10);
		      
		  } else   {N_KERN_TM;}
		  PREFETCH_U3;
		  PREFETCH_PSI;

		  chi = chi_p + SPINOR_SIZE * xyzt;

#if 0		  
		  {//deb
		    int gx=GJP.XnodeCoor()*lx+x;
		    int gy=GJP.YnodeCoor()*ly+y;
		    int gz=GJP.ZnodeCoor()*lz+z;
		    int gt=GJP.TnodeCoor()*lt+t;

		    if(gx==4&&gy==0&&gz==0&&gt==1)
		      {
			::printf("edge (%d %d %d %d) %4.3e %4.3e\n",
			     gx,gy,gz,gt,
			     *(0+(double*)&(wxm[0])),
			     *(1+(double*)&(wxm[0])));
		    }
		  }
#endif

		  N_STORE_XP_noadd;
		  N_STORE_YP;
		  N_STORE_ZP;
		  N_STORE_TP;

		  N_STORE_XM;
		  N_STORE_YM;
		  N_STORE_ZM;
		  N_STORE_TM;

		  PREFETCH_CHI;

		} else
#endif  //DEBUG_NOEDGE
		  {// if (yzt_edge & (x - x_check) == 0) 

		  ZERO;
#ifdef POS_DIR_ON	    
#ifdef X_DIR_ON
		  u = u_p + GAUGE_SIZE * (xyzt + vol * cbn);
		  psi = psi_p + SPINOR_SIZE *xpyzt;
		  N_KERN_XP;
		  PREFETCH_U0;
		  PREFETCH_PSI;
#endif
#ifdef Y_DIR_ON
		  u = u_p + GAUGE_SIZE * (xyzt + vol * cbn);
		  psi = psi_p +SPINOR_SIZE * xypzt;
		  N_KERN_YP;
		  PREFETCH_U1;
		  PREFETCH_PSI;
#endif
#ifdef Z_DIR_ON

		  u = u_p + GAUGE_SIZE * (xyzt + vol * cbn);
		  psi = psi_p +SPINOR_SIZE * xyzpt;
		  N_KERN_ZP;
		  PREFETCH_U2;
		  PREFETCH_PSI;
#endif
#ifdef T_DIR_ON
		  u = u_p + GAUGE_SIZE * (xyzt + vol * cbn);
		  psi = psi_p +SPINOR_SIZE * xyztp;
		  N_KERN_TP;
		  PREFETCH_U3;
		  PREFETCH_PSI;
#endif
#endif //POS_DIR_ON


#ifdef NEG_DIR_ON

#ifdef X_DIR_ON
		  u = u_p + GAUGE_SIZE * (xmyzt + vol * cb);
		  psi = psi_p +SPINOR_SIZE * xmyzt;
		  N_KERN_XM;
		  PREFETCH_U0;
		  PREFETCH_PSI;
#endif
#ifdef Y_DIR_ON
		  u = u_p + GAUGE_SIZE * (xymzt + vol * cb);
		  psi = psi_p +SPINOR_SIZE * xymzt;
		  N_KERN_YM;
		  PREFETCH_U1;
		  PREFETCH_PSI;
#endif
#ifdef Z_DIR_ON

		  u = u_p + GAUGE_SIZE * (xyzmt + vol * cb);
		  psi = psi_p +SPINOR_SIZE * xyzmt;
		  N_KERN_ZM;
		  PREFETCH_U2;
		  PREFETCH_PSI;
#endif
#ifdef T_DIR_ON

		  u = u_p + GAUGE_SIZE * (xyztm + vol * cb);
		  psi = psi_p +SPINOR_SIZE * xyztm;
		  N_KERN_TM;
		  PREFETCH_U3;
		  PREFETCH_PSI;
#endif
#endif //NEG_DIR_ON
		  chi = chi_p + SPINOR_SIZE * xyzt;

		  N_STORE_XP_noadd;
		  N_STORE_YP;
		  N_STORE_ZP;
		  N_STORE_TP;


		  N_STORE_XM;
		  N_STORE_YM;
		  N_STORE_ZM;
		  N_STORE_TM;

		  //STORE(chi);
		  PREFETCH_CHI;
		  }// else of if (yzt_edge & (x - x_check) == 0) 
	      }//  for (x = cbn ^ ((y + z + t) & 1); x < lx; x += 2) 
	    }
	  }
	}
      }
    }  



#ifdef DEBUG_PRINT_DAG1
  int dag=1;
  int mpisize=NumNodes();
  int mpirank=UniqueID();
  //  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);    
  //  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  for(int irank=0; irank< mpisize;++irank)
    {
      int cbn=!cb;
      
      if(mpirank == irank)
	{
	  int idx=0;
	  for(int t=0; t<lt; t++){
	    for(int z=0; z<lz; z++){
	      for(int y=0; y<ly; y++){
		for(int x=0; x<lx; x++){
		  int parity = x+y+z+t;
		  parity = parity % 2;
		  if(parity == cbn){
		    ::printf("debug3 %d %d (%d %d %d %d)  ",
			     cb,dag,
			     GJP.XnodeCoor()*lx+x,
			     GJP.YnodeCoor()*ly+y,
			     GJP.ZnodeCoor()*lz+z,
			     GJP.TnodeCoor()*lt+t);
		    for(int i=0;i<24;++i)
		      //::printf(" %4.3e", *(chi_p+idx+i));
		      ::printf(" %e", *(chi_p+idx+i));
		    ::printf("\n");
		    idx +=24;
		  }
		}
	      }
	    }
	  }

	}

      //MPI_Barrier( MPI_COMM_WORLD);
      cps::sync();
    }
  
  //exit(1);
#endif
}
