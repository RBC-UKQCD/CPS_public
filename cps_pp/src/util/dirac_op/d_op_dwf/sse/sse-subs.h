/* wilson dlash kernel pieces */
#if 1



// original opteron
//
//#define DIST_U 4
//#define DIST_PSI 2
//#define DIST_CHI 4


#define DIST_U 6
#define DIST_PSI 4
#define DIST_CHI 4


//#define PREFETCH_U						\
//  _mm_prefetch((char *)(u + GAUGE_SIZE * 4), _MM_HINT_T0);

#define PREFETCH_U0      \
  _mm_prefetch((char *)(u + GAUGE_SIZE * DIST_U), _MM_HINT_T0);              \
  _mm_prefetch((char *)(u + GAUGE_SIZE * DIST_U + 8), _MM_HINT_T0);	\

#define PREFETCH_U1     \
  _mm_prefetch((char *)(u + GAUGE_SIZE * DIST_U + 18), _MM_HINT_T0);         \
   _mm_prefetch((char *)(u + GAUGE_SIZE * DIST_U + 18 + 8), _MM_HINT_T0); \

#define PREFETCH_U2     \
  _mm_prefetch((char *)(u + GAUGE_SIZE * DIST_U + 18 * 2), _MM_HINT_T0);     \
   _mm_prefetch((char *)(u + GAUGE_SIZE * DIST_U + 18 * 2 + 8), _MM_HINT_T0); \

#define PREFETCH_U3     \
  _mm_prefetch((char *)(u + GAUGE_SIZE * DIST_U + 18 * 3), _MM_HINT_T0);     \
   _mm_prefetch((char *)(u + GAUGE_SIZE * DIST_U + 18 * 3 + 8), _MM_HINT_T0); \

#define PREFETCH_PSI                                                    \
  _mm_prefetch((char *)(psi + SPINOR_SIZE * DIST_PSI), _MM_HINT_T0); \

#define PREFETCH_CHI    \
      _mm_prefetch((char *)(chi + SPINOR_SIZE * DIST_CHI), _MM_HINT_T0);
#endif

#if 0
#define PREFETCH_U      \
  _mm_prefetch((char *)(u + GAUGE_SIZE * 4), _MM_HINT_T0);
#define PREFETCH_U0      \
  _mm_prefetch((char *)(u + GAUGE_SIZE * 4), _MM_HINT_T0);              \
  _mm_prefetch((char *)(u + GAUGE_SIZE * 4 + 8), _MM_HINT_T0);          \
  _mm_prefetch((char *)(u + GAUGE_SIZE * 4 + 16), _MM_HINT_T0);          \

#define PREFETCH_U1     \
  _mm_prefetch((char *)(u + GAUGE_SIZE * 4 + 18), _MM_HINT_T0);         \
  _mm_prefetch((char *)(u + GAUGE_SIZE * 4 + 18 + 8), _MM_HINT_T0);     \
  _mm_prefetch((char *)(u + GAUGE_SIZE * 4 + 18 + 16), _MM_HINT_T0);     \

#define PREFETCH_U2     \
  _mm_prefetch((char *)(u + GAUGE_SIZE * 4 + 18 * 2), _MM_HINT_T0);     \
  _mm_prefetch((char *)(u + GAUGE_SIZE * 4 + 18 * 2 + 8), _MM_HINT_T0); \
  _mm_prefetch((char *)(u + GAUGE_SIZE * 4 + 18 * 2 + 16), _MM_HINT_T0); \

#define PREFETCH_U3     \
  _mm_prefetch((char *)(u + GAUGE_SIZE * 4 + 18 * 3), _MM_HINT_T0);     \
  _mm_prefetch((char *)(u + GAUGE_SIZE * 4 + 18 * 3 + 8), _MM_HINT_T0); \
  _mm_prefetch((char *)(u + GAUGE_SIZE * 4 + 18 * 3 + 16), _MM_HINT_T0); \

#define PREFETCH_PSI                                                    \
  _mm_prefetch((char *)(psi + SPINOR_SIZE * 2), _MM_HINT_T0); \

#define PREFETCH_CHI    \
  _mm_prefetch((char *)(chi + SPINOR_SIZE * 4), _MM_HINT_T0);
#endif


#if 0

#define PREFETCH_U0 TOUCH(u+72,72);
#define PREFETCH_U1
#define PREFETCH_U2
#define PREFETCH_U3

#define PREFETCH_PSI  TOUCH( psi+24,24);
#define PREFETCH_CHI

#endif

//-----------------------------------------

#define BDIST_U 4
#define BDIST_PSI 4
#define BDIST_CHI 4

#define BND_PREFETCH_U0      \
  _mm_prefetch((char *)(u + GAUGE_SIZE * BDIST_U), _MM_HINT_T0);              \
  _mm_prefetch((char *)(u + GAUGE_SIZE * BDIST_U + 8), _MM_HINT_T0);	\

#define BND_PREFETCH_U1     \
  _mm_prefetch((char *)(u + GAUGE_SIZE * BDIST_U + 18), _MM_HINT_T0);         \
   _mm_prefetch((char *)(u + GAUGE_SIZE * BDIST_U + 18 + 8), _MM_HINT_T0); \

#define BND_PREFETCH_U2     \
  _mm_prefetch((char *)(u + GAUGE_SIZE * BDIST_U + 18 * 2), _MM_HINT_T0);     \
   _mm_prefetch((char *)(u + GAUGE_SIZE * BDIST_U + 18 * 2 + 8), _MM_HINT_T0); \

#define BND_PREFETCH_U3     \
  _mm_prefetch((char *)(u + GAUGE_SIZE * BDIST_U + 18 * 3), _MM_HINT_T0);     \
   _mm_prefetch((char *)(u + GAUGE_SIZE * BDIST_U + 18 * 3 + 8), _MM_HINT_T0); \

#define BND_PREFETCH_PSI                                                    \
  _mm_prefetch((char *)(psi + SPINOR_SIZE * BDIST_PSI), _MM_HINT_T0); \

#define BND_PREFETCH_CHI    \
      _mm_prefetch((char *)(chi + SPINOR_SIZE * BDIST_CHI), _MM_HINT_T0);


//---------------------------------------------------------------




#ifdef ADD2REG
#define ZERO                                            \
   t00 = _mm_xor_pd(t00, t00);				\
   t01 = _mm_xor_pd(t01, t01);                          \
   t02 = _mm_xor_pd(t02, t02);                          \
   t03 = _mm_xor_pd(t03, t03);                          \
   t04 = _mm_xor_pd(t04, t04);                          \
   t05 = _mm_xor_pd(t05, t05);                          \
   t06 = _mm_xor_pd(t06, t06);                          \
   t07 = _mm_xor_pd(t07, t07);                          \
   t08 = _mm_xor_pd(t08, t08);                          \
   t09 = _mm_xor_pd(t09, t09);                          \
   t10 = _mm_xor_pd(t10, t10);                          \
   t11 = _mm_xor_pd(t11, t11);				\

#else
#define ZERO_(i)					   \
  wxp[i]=_mm_xor_pd(wxp[i],wxp[i]);			   \
  wyp[i]=_mm_xor_pd(wyp[i],wyp[i]);				   \
  wzp[i]=_mm_xor_pd(wzp[i],wzp[i]);				   \
  wtp[i]=_mm_xor_pd(wtp[i],wtp[i]);				   \
								   \
  wxm[i]=_mm_xor_pd(wxm[i],wxm[i]);				   \
  wym[i]=_mm_xor_pd(wym[i],wym[i]);				   \
  wzm[i]=_mm_xor_pd(wzm[i],wzm[i]);				   \
  wtm[i]=_mm_xor_pd(wtm[i],wtm[i]);				   \

#define ZERO					\
  ZERO_(0); ZERO_(1); ZERO_(2);			\
  ZERO_(3); ZERO_(4); ZERO_(5);			\

#endif



#ifdef ADD2REG
#ifdef USE_HERN
#define DECLARE						\
  register __m128d t00, t01, t02, t03, t04, t05;	\
  register __m128d t06, t07, t08, t09, t10, t11;	\
  register __m128d _a, _b, _c, _d;			\

#else
#define DECLARE						\
  register __m128d t00, t01, t02, t03, t04, t05;	\
  register __m128d t06, t07, t08, t09, t10, t11;	\

#endif
#else
#ifdef USE_HERN
#define DECLARE						\
  __m128d __RESTRICT wxp[6];				\
  __m128d __RESTRICT wyp[6];				\
  __m128d __RESTRICT wzp[6];				\
  __m128d __RESTRICT wtp[6];				\
							\
  __m128d __RESTRICT wxm[6];				\
  __m128d __RESTRICT wym[6];				\
  __m128d __RESTRICT wzm[6];				\
  __m128d __RESTRICT wtm[6];				\
  register __m128d _a, _b, _c, _d;			\
  
#else
#define DECLARE						\
  __m128d __RESTRICT wxp[6];				\
  __m128d __RESTRICT wyp[6];				\
  __m128d __RESTRICT wzp[6];				\
  __m128d __RESTRICT wtp[6];				\
							\
  __m128d __RESTRICT wxm[6];				\
  __m128d __RESTRICT wym[6];				\
  __m128d __RESTRICT wzm[6];				\
  __m128d __RESTRICT wtm[6];				\

#endif
#endif





#define STORE(addr)                                     \
   _mm_store_pd(addr + 0, t00);                         \
   _mm_store_pd(addr + 2, t01);                         \
   _mm_store_pd(addr + 4, t02);                         \
   _mm_store_pd(addr + 6, t03);                         \
   _mm_store_pd(addr + 8, t04);                         \
   _mm_store_pd(addr + 10, t05);                        \
   _mm_store_pd(addr + 12, t06);                        \
   _mm_store_pd(addr + 14, t07);                        \
   _mm_store_pd(addr + 16, t08);                        \
   _mm_store_pd(addr + 18, t09);                        \
   _mm_store_pd(addr + 20, t10);                        \
   _mm_store_pd(addr + 22, t11);                        \

/*  sdag = +1,  a = TMP(*,c,0), c = TMP(*,c,3)  */
//    c  = I a
// a = [psi(0).r + psi(3).i,  psi(0).i - psi(3).re]
#define P_KERN_XP_03(C)                                 \
   _a = _mm_load_pd(psi + 0 + C * 2);                   \
   _c = _mm_load_pd(psi + 18 + C * 2);                  \
   _a = _mm_shuffle_pd(_a, _a, 1);			\
   _a = _mm_addsub_pd(_a, _c);                          \
   _a = _mm_shuffle_pd(_a, _a, 1);

#define PPX0_A0(C) _mm_load_pd(psi + 0 + C * 2) 
#define PPX0_C0(C) _mm_load_pd(psi + 18 + C * 2) 
#define PPX0_A1(C) _mm_shuffle_pd(PPX0_A0(C), PPX0_A0(C), 1) 
#define PPX0_A2(C) _mm_addsub_pd(PPX0_A1(C), PPX0_C0(C)) 
#define PPX0_A3(C) _mm_shuffle_pd(PPX0_A2(C),PPX0_A2(C), 1) 

#define P_PROJ_X_03(C) PPX0_A3(C)


/*  sdag = +1,  a = TMP(*,c,1), c = TMP(*,c,2)  */
//    c = I a
#define P_KERN_XP_12(C)                                 \
   _a = _mm_load_pd(psi + 6 + C * 2);                   \
   _c = _mm_load_pd(psi + 12 + C * 2);                  \
   _a = _mm_shuffle_pd(_a, _a, 1);			\
   _a = _mm_addsub_pd(_a, _c);                          \
   _a = _mm_shuffle_pd(_a, _a, 1);

#define PPX1_A0(C) _mm_load_pd(psi + 6 + C * 2)
#define PPX1_C0(C) _mm_load_pd(psi + 12 + C * 2)
#define PPX1_A1(C) _mm_shuffle_pd(PPX1_A0(C), PPX1_A0(C), 1)
#define PPX1_A2(C) _mm_addsub_pd(PPX1_A1(C), PPX1_C0(C))
#define PPX1_A3(C) _mm_shuffle_pd(PPX1_A2(C), PPX1_A2(C), 1)
#define P_PROJ_X_12(C) PPX1_A3(C)

/*  sdag = -1,  a = TMP(*,c,0), c = TMP(*,c,3)  */
//   c = -I a
#define N_KERN_XP_03(C)                                 \
   _a = _mm_load_pd(psi + 0 + C * 2);                   \
   _c = _mm_load_pd(psi + 18 + C * 2);                  \
   _c = _mm_shuffle_pd(_c, _c, 1);			\
   _a = _mm_addsub_pd(_a, _c);                          \

#define NPX0_A0(C) _mm_load_pd(psi + 0 + C * 2)
#define NPX0_C0(C) _mm_load_pd(psi + 18 + C * 2)
#define NPX0_C1(C) _mm_shuffle_pd(NPX0_C0(C), NPX0_C0(C), 1)
#define NPX0_A1(C)  _mm_addsub_pd(NPX0_A0(C), NPX0_C1(C))
#define N_PROJ_X_03(C)  NPX0_A1(C)


/*  sdag = -1,  a = TMP(*,c,1), c = TMP(*,c,2)  */
//   c = -I a
#define N_KERN_XP_12(C)                                 \
   _a = _mm_load_pd(psi + 6 + C * 2);                   \
   _c = _mm_load_pd(psi + 12 + C * 2);                  \
   _c = _mm_shuffle_pd(_c, _c, 1);                 \
   _a = _mm_addsub_pd(_a, _c);                          \

#define NPX1_A0(C) _mm_load_pd(psi + 6 + C * 2)
#define NPX1_C0(C) _mm_load_pd(psi + 12 + C * 2)
#define NPX1_C1(C) _mm_shuffle_pd(NPX1_C0(C), NPX1_C0(C), 1)
#define NPX1_A1(C) _mm_addsub_pd(NPX1_A0(C), NPX1_C1(C))
#define N_PROJ_X_12(C) NPX1_A1(C)

/*  sdag = +1, a = TMP(*,c,0), c = TMP(*,c,3)  */
// c = a
#define P_KERN_YP_03(C)                                 \
   _a = _mm_load_pd(psi + 0 + C * 2);                   \
   _c = _mm_load_pd(psi + 18 + C * 2);                  \
   _a = _mm_add_pd(_a, _c);


#define PPY0_A0(C) _mm_load_pd(psi + 0 + C * 2)
#define PPY0_C(C) _mm_load_pd(psi + 18 + C * 2)
#define PPY0_A1(C) _mm_add_pd(PPY0_A0(C), PPY0_C(C))
#define P_PROJ_Y_03(C) PPY0_A1(C)

/*  sdag = +1, a = TMP(*,c,1), c = TMP(*,c,2)  */
// c = -a
#define P_KERN_YP_12(C)                                 \
   _a = _mm_load_pd(psi + 6 + C * 2);              \
   _c = _mm_load_pd(psi + 12 + C * 2);                  \
   _a = _mm_sub_pd(_a, _c);                             \

#define PPY1_A0(C) _mm_load_pd(psi + 6 + C * 2)
#define PPY1_C0(C) _mm_load_pd(psi + 12 + C * 2)
#define PPY1_A1(C)  _mm_sub_pd(PPY1_A0(C),PPY1_C0(C))
#define P_PROJ_Y_12(C) PPY1_A1(C)

/*  sdag = -1, a = TMP(*,c,0), c = TMP(*,c,3)  */
// c = -a
#define N_KERN_YP_03(C)                                 \
   _a = _mm_load_pd(psi + 0 + C * 2);              \
   _c = _mm_load_pd(psi + 18 + C * 2);                  \
   _a = _mm_sub_pd(_a, _c);                             \

#define NPY0_A0(C) _mm_load_pd(psi + 0 + C * 2)
#define NPY0_C0(C) _mm_load_pd(psi + 18 + C * 2)
#define NPY0_A1(C) _mm_sub_pd(NPY0_A0(C), NPY0_C0(C))
#define N_PROJ_Y_03(C) NPY0_A1(C)

/*  sdag = -1, a = TMP(*,c,1), c = TMP(*,c,2)  */
// c = a
#define N_KERN_YP_12(C)                                 \
  _a = _mm_load_pd(psi + 6 + C * 2);			\
  _c = _mm_load_pd(psi + 12 + C * 2);			\
  _a = _mm_add_pd(_a, _c);

#define NPY1_A0(C) _mm_load_pd(psi + 6 + C * 2)
#define NPY1_C0(C) _mm_load_pd(psi + 12 + C * 2)
#define NPY1_A1(C) _mm_add_pd(NPY1_A0(C), NPY1_C0(C))
#define N_PROJ_Y_12(C) NPY1_A1(C)

/*  sdag = +1,  a = TMP(*,c,0), c = TMP(*,c,2)  */
// c = I a
#define P_KERN_ZP_02(C)                                 \
   _a = _mm_load_pd(psi + 0 + C * 2);                   \
   _c = _mm_load_pd(psi + 12 + C * 2);                  \
   _a = _mm_shuffle_pd(_a, _a, 1);			\
   _a = _mm_addsub_pd(_a, _c);                          \
   _a = _mm_shuffle_pd(_a, _a, 1);

#define PPZ0_A0(C) _mm_load_pd(psi + 0 + C * 2)
#define PPZ0_C0(C) _mm_load_pd(psi + 12 + C * 2)
#define PPZ0_A1(C) _mm_shuffle_pd(PPZ0_A0(C), PPZ0_A0(C), 1)
#define PPZ0_A2(C) _mm_addsub_pd(PPZ0_A1(C), PPZ0_C0(C))
#define PPZ0_A3(C) _mm_shuffle_pd(PPZ0_A2(C), PPZ0_A2(C), 1)
#define P_PROJ_Z_02(C) PPZ0_A3(C)

/*  sdag = +1,  a = TMP(*,c,1), c = TMP(*,c,3)  */
// c = -I a
#define P_KERN_ZP_13(C)                                 \
   _a = _mm_load_pd(psi + 6 + C * 2);                   \
   _c = _mm_load_pd(psi + 18 + C * 2);                  \
   _c = _mm_shuffle_pd(_c, _c, 1);			\
   _a = _mm_addsub_pd(_a, _c);                          \

#define PPZ1_A0(C) _mm_load_pd(psi + 6 + C * 2)
#define PPZ1_C0(C) _mm_load_pd(psi + 18 + C * 2)
#define PPZ1_C1(C) _mm_shuffle_pd(PPZ1_C0(C),PPZ1_C0(C), 1)
#define PPZ1_A1(C) _mm_addsub_pd(PPZ1_A0(C),PPZ1_C1(C))
#define P_PROJ_Z_13(C) PPZ1_A1(C)

/*  sdag = -1,  a = TMP(*,c,0), c = TMP(*,c,2)  */
// c = -I a
#define N_KERN_ZP_02(C)                                 \
   _a = _mm_load_pd(psi + 0 + C * 2);                   \
   _c = _mm_load_pd(psi + 12 + C * 2);                  \
   _c = _mm_shuffle_pd(_c, _c, 1);                 \
   _a = _mm_addsub_pd(_a, _c);                          \

#define NPZ0_A0(C) _mm_load_pd(psi + 0 + C * 2)
#define NPZ0_C0(C) _mm_load_pd(psi + 12 + C * 2)
#define NPZ0_C1(C) _mm_shuffle_pd(NPZ0_C0(C), NPZ0_C0(C), 1)
#define NPZ0_A1(C) _mm_addsub_pd(NPZ0_A0(C), NPZ0_C1(C))
#define N_PROJ_Z_02(C) NPZ0_A1(C)

/*  sdag = -1,  a = TMP(*,c,1), c = TMP(*,c,3)  */
// c = I a
#define N_KERN_ZP_13(C)                                 \
   _a = _mm_load_pd(psi + 6 + C * 2);                   \
   _c = _mm_load_pd(psi + 18 + C * 2);                  \
   _a = _mm_shuffle_pd(_a, _a, 1);                 \
   _a = _mm_addsub_pd(_a, _c);                          \
   _a = _mm_shuffle_pd(_a, _a, 1);

#define NPZ1_A0(C) _mm_load_pd(psi + 6 + C * 2)
#define NPZ1_C0(C) _mm_load_pd(psi + 18 + C * 2)
#define NPZ1_A1(C) _mm_shuffle_pd(NPZ1_A0(C), NPZ1_A0(C), 1)
#define NPZ1_A2(C) _mm_addsub_pd(NPZ1_A1(C), NPZ1_C0(C))
#define NPZ1_A3(C) _mm_shuffle_pd(NPZ1_A2(C), NPZ1_A2(C), 1)
#define N_PROJ_Z_13(C) NPZ1_A3(C)

/*  sdag = +1, a = TMP(*,c,0), c = TMP(*,c,2)  */
// c = -a
#define P_KERN_TP_02(C)                                 \
   _a = _mm_load_pd(psi + 0 + C * 2);              \
   _c = _mm_load_pd(psi + 12 + C * 2);                  \
   _a = _mm_sub_pd(_a, _c);                             \

#define PPT0_A0(C) _mm_load_pd(psi + 0 + C * 2)
#define PPT0_C0(C) _mm_load_pd(psi + 12 + C * 2)
#define PPT0_A1(C) _mm_sub_pd(PPT0_A0(C), PPT0_C0(C))
#define P_PROJ_T_02(C) PPT0_A1(C)

/*  sdag = +1, a = TMP(*,c,1), c = TMP(*,c,3)  */
// c = -a
#define P_KERN_TP_13(C)                                 \
   _a = _mm_load_pd(psi + 6 + C * 2);              \
   _c = _mm_load_pd(psi + 18 + C * 2);                  \
   _a = _mm_sub_pd(_a, _c);                             \

#define PPT1_A0(C) _mm_load_pd(psi + 6 + C * 2)
#define PPT1_C0(C) _mm_load_pd(psi + 18 + C * 2)
#define PPT1_A1(C) _mm_sub_pd(PPT1_A0(C), PPT1_C0(C))
#define P_PROJ_T_13(C) PPT1_A1(C)

/*  sdag = -1, a = TMP(*,c,0), c = TMP(*,c,2)  */
// c = a
#define N_KERN_TP_02(C)                                 \
   _a = _mm_load_pd(psi + 0 + C * 2);                   \
   _c = _mm_load_pd(psi + 12 + C * 2);                  \
   _a = _mm_add_pd(_a, _c);

#define NPT0_A0(C) _mm_load_pd(psi + 0 + C * 2)
#define NPT0_C0(C) _mm_load_pd(psi + 12 + C * 2)
#define NPT0_A1(C) _mm_add_pd(NPT0_A0(C), NPT0_C0(C))
#define N_PROJ_T_02(C) NPT0_A1(C) 

/*  sdag = -1, a = TMP(*,c,1), c = TMP(*,c,3)  */
// c = a
#define N_KERN_TP_13(C)                                 \
   _a = _mm_load_pd(psi + 6 + C * 2);                   \
   _c = _mm_load_pd(psi + 18 + C * 2);                  \
   _a = _mm_add_pd(_a, _c);

#define NPT1_A0(C) _mm_load_pd(psi + 6 + C * 2)
#define NPT1_C0(C) _mm_load_pd(psi + 18 + C * 2)
#define NPT1_A1(C) _mm_add_pd(NPT1_A0(C), NPT1_C0(C))
#define N_PROJ_T_13(C) NPT1_A1(C)


#define P_KERN(t0,t1,t2,t3,t4,t5,C,MU)                  \
   _d = _b = _mm_loaddup_pd(u + 0 + C * 6 + 18 * MU);   \
   _b = _mm_mul_pd(_b, _a);                             \
   t0 = _mm_add_pd(t0, _b);                             \
   _d = _mm_mul_pd(_d, _c);                             \
   t1 = _mm_add_pd(t1, _d);                             \
   _d = _b = _mm_loaddup_pd(u + 2 + C * 6 + 18 * MU);   \
   _b = _mm_mul_pd(_b, _a);                             \
   t2 = _mm_add_pd(t2, _b);                             \
   _d = _mm_mul_pd(_d, _c);                             \
   t3 = _mm_add_pd(t3, _d);                             \
   _d = _b = _mm_loaddup_pd(u + 4 + C * 6 + 18 * MU);   \
   _b = _mm_mul_pd(_b, _a);                             \
   t4 = _mm_add_pd(t4, _b);                             \
   _d = _mm_mul_pd(_d, _c);                             \
   t5 = _mm_add_pd(t5, _d);                             \
   _a = _mm_shuffle_pd(_a, _a, 1);                      \
   _c = _mm_shuffle_pd(_c, _c, 1);                      \
   _d = _b = _mm_loaddup_pd(u + 1 + C * 6 + 18 * MU);   \
   _b = _mm_mul_pd(_b, _a);                             \
   t0 = _mm_addsub_pd(t0, _b);                          \
   _d = _mm_mul_pd(_d, _c);                             \
   t1 = _mm_addsub_pd(t1, _d);                          \
   _d = _b = _mm_loaddup_pd(u + 3 + C * 6 + 18 * MU);   \
   _b = _mm_mul_pd(_b, _a);                             \
   t2 = _mm_addsub_pd(t2, _b);                          \
   _d = _mm_mul_pd(_d, _c);                             \
   t3 = _mm_addsub_pd(t3, _d);                          \
   _d = _b = _mm_loaddup_pd(u + 5 + C * 6 + 18 * MU);   \
   _b = _mm_mul_pd(_b, _a);                             \
   t4 = _mm_addsub_pd(t4, _b);                          \
   _d = _mm_mul_pd(_d, _c);                             \
   t5 = _mm_addsub_pd(t5, _d);                          \


// c = I a
//
//   (u_r a)     = [ u_r a_r, u_r a_i ]
//   (CI u_i a)  = [ -u_i a_i, u_i a_r ]
//   (u_r c)     = (u_r I a)
//               = [ -u_r a_i, u_r a_r ]
//   (CI u_i c)  = -(u_i  a)
//               = [-u_i a_r, -u_i a_i ]
//
#define P_KERN_IEQ(t0,t1,t2,t3,t4,t5,C,MU)		\
  _b = _mm_loaddup_pd(u + 0 + C * 6 + 18 * MU);		\
  _b = _mm_mul_pd(_b, _a);				\
  t0 = _mm_add_pd(t0, _b);				\
  _c = _mm_shuffle_pd(_b,_b,1);				\
  t1 = _mm_addsub_pd(t1, _c);				\
  _d = _mm_loaddup_pd(u + 2 + C * 6 + 18 * MU);		\
  _d = _mm_mul_pd(_d, _a);				\
  t2 = _mm_add_pd(t2, _d);				\
  _c = _mm_shuffle_pd(_d,_d,1);				\
  t3 = _mm_addsub_pd(t3, _c);				\
  _b = _mm_loaddup_pd(u + 4 + C * 6 + 18 * MU);		\
  _b = _mm_mul_pd(_b, _a);				\
  t4 = _mm_add_pd(t4, _b);				\
  _c = _mm_shuffle_pd(_b,_b,1);				\
  t5 = _mm_addsub_pd(t5, _c);				\
  \
  _a = _mm_shuffle_pd(_a,_a,1);                         \
  _d = _mm_loaddup_pd(u + 1 + C * 6 + 18 * MU);	        \
  _d = _mm_mul_pd(_d,_a);                               \
  t0 = _mm_addsub_pd(t0, _d);				\
  _c = _mm_shuffle_pd(_d,_d,1);         		\
  t1 = _mm_sub_pd(t1, _c);				\
  _b = _mm_loaddup_pd(u + 3 + C * 6 + 18 * MU);	        \
  _b = _mm_mul_pd(_b, _a);				\
  t2 = _mm_addsub_pd(t2, _b);				\
  _c = _mm_shuffle_pd(_b,_b,1);	         		\
  t3 = _mm_sub_pd(t3, _c);				\
  _d = _mm_loaddup_pd(u + 5 + C * 6 + 18 * MU);	\
  _d = _mm_mul_pd(_d, _a);				\
  t4 = _mm_addsub_pd(t4, _d);				\
  _c = _mm_shuffle_pd(_d,_d,1);				\
  t5 = _mm_sub_pd(t5, _c);				\


#if 1
// This was OK, NAME orig
// a bit faster on Nehalem
#define P_HERN(t0,t1,t2,C,MU)				\
  _b = _mm_loaddup_pd(u + 0 + C * 6 + 18 * MU);		\
  t0 = _mm_add_pd(t0,_mm_mul_pd(_b, _a));		\
  _d = _mm_loaddup_pd(u + 2 + C * 6 + 18 * MU);		\
  t1 = _mm_add_pd(t1,_mm_mul_pd(_d, _a));		\
  _c = _mm_loaddup_pd(u + 4 + C * 6 + 18 * MU);		\
  t2 = _mm_add_pd(t2,_mm_mul_pd(_c, _a));		\
  							\
  _a = _mm_shuffle_pd(_a,_a,1);                         \
  								\
  _b = _mm_mul_pd(_mm_loaddup_pd(u + 1 + C * 6 + 18 * MU),_a);	\
  t0 = _mm_addsub_pd(t0,_b);					\
  _d = _mm_mul_pd(_mm_loaddup_pd(u + 3 + C * 6 + 18 * MU),_a);	\
  t1 = _mm_addsub_pd(t1,_d);					\
  _c = _mm_mul_pd(_mm_loaddup_pd(u + 5 + C * 6 + 18 * MU),_a);	\
  t2 = _mm_addsub_pd(t2,_c);					\

// this was slower than P_HERN_NOADD even after removing the initial ZERO
//
#define P_HERN_noadd(t0,t1,t2,C,MU)			\
  _b = _mm_loaddup_pd(u + 0 + C * 6 + 18 * MU);		\
  t0 = _mm_mul_pd(_b, _a);				\
  _d = _mm_loaddup_pd(u + 2 + C * 6 + 18 * MU);		\
  t1 = _mm_mul_pd(_d, _a);				\
  _c = _mm_loaddup_pd(u + 4 + C * 6 + 18 * MU);		\
  t2 = _mm_mul_pd(_c, _a);				\
  							\
  _a = _mm_shuffle_pd(_a,_a,1);                         \
  								\
  _b = _mm_mul_pd(_mm_loaddup_pd(u + 1 + C * 6 + 18 * MU),_a);	\
  t0 = _mm_addsub_pd(t0,_b);					\
  _d = _mm_mul_pd(_mm_loaddup_pd(u + 3 + C * 6 + 18 * MU),_a);	\
  t1 = _mm_addsub_pd(t1,_d);					\
  _c = _mm_mul_pd(_mm_loaddup_pd(u + 5 + C * 6 + 18 * MU),_a);	\
  t2 = _mm_addsub_pd(t2,_c);					\

#endif

#if 0
// Reading either of real/imaginary part of Psi rather than those of Link
#define P_HERN(t0,t1,t2,C,MU)						\
  _b = _mm_unpacklo_pd( _a, _a );					\
  t0 = _mm_add_pd(t0,_mm_mul_pd(_b, *(__m128d*)(u + 0 + C * 6 + 18 * MU))); \
  t1 = _mm_add_pd(t1,_mm_mul_pd(_b, *(__m128d*)(u + 2 + C * 6 + 18 * MU))); \
  t2 = _mm_add_pd(t2,_mm_mul_pd(_b, *(__m128d*)(u + 4 + C * 6 + 18 * MU))); \
  _c = _mm_unpackhi_pd( _a, _a );					\
  t0 = _mm_addsub_pd(t0,_mm_mul_pd(_c, _mm_shuffle_pd(*(__m128d*)(u + 0 + C * 6 + 18 * MU), *(__m128d*)(u + 0 + C * 6 + 18 * MU),1))); \
  t1 = _mm_addsub_pd(t1,_mm_mul_pd(_c, _mm_shuffle_pd(*(__m128d*)(u + 2 + C * 6 + 18 * MU), *(__m128d*)(u + 2 + C * 6 + 18 * MU),1))); \
  t2 = _mm_addsub_pd(t2,_mm_mul_pd(_c, _mm_shuffle_pd(*(__m128d*)(u + 4 + C * 6 + 18 * MU), *(__m128d*)(u + 4 + C * 6 + 18 * MU),1))); \
  
#endif


#define PH_B(A) _mm_unpacklo_pd( A, A )
#define PH_C(A)  _mm_unpackhi_pd( A, A )
#define P_HERN2(t0,t1,t2,C,MU,A)					\
  t0 = _mm_add_pd(t0,_mm_mul_pd(PH_B(A), *(__m128d*)(u + 0 + C * 6 + 18 * MU))); \
  t1 = _mm_add_pd(t1,_mm_mul_pd(PH_B(A), *(__m128d*)(u + 2 + C * 6 + 18 * MU))); \
  t2 = _mm_add_pd(t2,_mm_mul_pd(PH_B(A), *(__m128d*)(u + 4 + C * 6 + 18 * MU))); \
  t0 = _mm_addsub_pd(t0,_mm_mul_pd(PH_C(A), _mm_shuffle_pd(*(__m128d*)(u + 0 + C * 6 + 18 * MU), *(__m128d*)(u + 0 + C * 6 + 18 * MU),1))); \
  t1 = _mm_addsub_pd(t1,_mm_mul_pd(PH_C(A), _mm_shuffle_pd(*(__m128d*)(u + 2 + C * 6 + 18 * MU), *(__m128d*)(u + 2 + C * 6 + 18 * MU),1))); \
  t2 = _mm_addsub_pd(t2,_mm_mul_pd(PH_C(A), _mm_shuffle_pd(*(__m128d*)(u + 4 + C * 6 + 18 * MU), *(__m128d*)(u + 4 + C * 6 + 18 * MU),1))); \
  


#if 0
// !  This part is wrong, showing the idea only  ! fixme
// Use the horizontal operation HADDPD, HSUBPD
// Reading either of real/imaginary part of Psi rather than those of Link
// a bit slower than orig on Nehalem
#define P_HERN(t0,t1,t2,C,MU)						\
  _b = _mm_shuffle_pd(_a,_a,1);						\
  									\
  t0 = _mm_add_pd(t0,_mm_hadd_pd(_mm_mul_pd(_a,*(__m128d*)(u + 0 + C * 6 + 18 * MU)),_mm_mul_pd(_b,*(__m128d*)(u + 0 + C * 6 + 18 * MU)))); \
  t1 = _mm_add_pd(t1,_mm_hadd_pd(_mm_mul_pd(_a,*(__m128d*)(u + 2 + C * 6 + 18 * MU)),_mm_mul_pd(_b,*(__m128d*)(u + 2 + C * 6 + 18 * MU)))); \
  t2 = _mm_add_pd(t2,_mm_hadd_pd(_mm_mul_pd(_a,*(__m128d*)(u + 4 + C * 6 + 18 * MU)),_mm_mul_pd(_b,*(__m128d*)(u + 4 + C * 6 + 18 * MU)))); \

#endif


#if 0
// Using unpack{lo,hi}  good as original
#define P_HERN(t0,t1,t2,C,MU)						\
  _b = _mm_unpacklo_pd(*((__m128d*)(u + 0 + C * 6 + 18 * MU)),*((__m128d*)(u + 0 + C * 6 + 18 * MU))); \
  t0 = _mm_add_pd(t0,_mm_mul_pd(_b, _a));		\
 _d = _mm_unpacklo_pd(*((__m128d*)(u + 2 + C * 6 + 18 * MU)),*((__m128d*)(u + 2 + C * 6 + 18 * MU))); \
  t1 = _mm_add_pd(t1,_mm_mul_pd(_d, _a));		\
  _c = _mm_unpacklo_pd(*((__m128d*)(u + 4 + C * 6 + 18 * MU)),*((__m128d*)(u + 4 + C * 6 + 18 * MU))); \
  t2 = _mm_add_pd(t2,_mm_mul_pd(_c, _a));				\
  									\
  _a = _mm_shuffle_pd(_a,_a,1);						\
  									\
  _b = _mm_unpackhi_pd(*((__m128d*)(u + 0 + C * 6 + 18 * MU)),*((__m128d*)(u + 0 + C * 6 + 18 * MU))); \
  t0 = _mm_addsub_pd(t0,_mm_mul_pd(_b,_a));				\
  _d = _mm_unpackhi_pd(*((__m128d*)(u + 2 + C * 6 + 18 * MU)),*((__m128d*)(u + 2 + C * 6 + 18 * MU))); \
  t1 = _mm_addsub_pd(t1,_mm_mul_pd(_d,_a));				\
  _c = _mm_unpackhi_pd(*((__m128d*)(u + 4 + C * 6 + 18 * MU)),*((__m128d*)(u + 4 + C * 6 + 18 * MU))); \
  t2 = _mm_addsub_pd(t2,_mm_mul_pd(_c,_a));				\

#endif

#if 0
// this is not different from ORIG
#define P_HERN(t0,t1,t2,C,MU)				\
  t0 = _mm_add_pd(t0,_mm_mul_pd(_mm_loaddup_pd(u + 0 + C * 6 + 18 * MU), _a)); \
  t1 = _mm_add_pd(t1,_mm_mul_pd(_mm_loaddup_pd(u + 2 + C * 6 + 18 * MU), _a)); \
  t2 = _mm_add_pd(t2,_mm_mul_pd(_mm_loaddup_pd(u + 4 + C * 6 + 18 * MU), _a)); \
  									\
  _a = _mm_shuffle_pd(_a,_a,1);						\
  									\
  t0 = _mm_addsub_pd(t0,_mm_mul_pd(_mm_loaddup_pd(u + 1 + C * 6 + 18 * MU),_a)); \
  t1 = _mm_addsub_pd(t1,_mm_mul_pd(_mm_loaddup_pd(u + 3 + C * 6 + 18 * MU),_a)); \
  t2 = _mm_addsub_pd(t2,_mm_mul_pd(_mm_loaddup_pd(u + 5 + C * 6 + 18 * MU),_a)); \

#endif

#define save2_P_KERN_IEQ(t0,t1,t2,t3,t4,t5,C,MU)		\
  _c = _mm_shuffle_pd(_a,_a,1);					\
  _c =  _mm_xor_pd(_d, _mm_set_pd(0.0,-0.0));			\
  P_KERN(t0,t1,t2,t3,t4,t5,C,MU);  

#define op_IEQ_U(i, __X,__Y, __tA, __tB,  C,MU)				\
  __Y = __X  = _mm_mul_pd(_mm_loaddup_pd(u + i + C * 6 + 18 * MU), _a);	\
  __tB = _mm_addsub_pd(__tB, __Y);					\
  __tA = _mm_add_pd(__tA, __X);						\
  
//  __tB = _mm_addsub_pd(__tB, _mm_shuffle_pd(__X,__X,1));	\

#define op_IEQ_D(i, __X,__Y, __tA, __tB,  C,MU)				\
  __X = _mm_mul_pd(_mm_loaddup_pd(u + i + C * 6 + 18 * MU), _a);	\
  __tA = _mm_addsub_pd(__tA, __X);					\
  __tB = _mm_addsub_pd(__tB, __X);					\
  //  __tB = _mm_sub_pd(__tB, _mm_shuffle_pd(__X,__X,1)) ;		\


#define save_3_P_KERN_IEQ(t0,t1,t2,t3,t4,t5,C,MU)			\
  op_IEQ_U(0, _b, _c, t0, t1, C,MU);				\
  op_IEQ_U(2, _d, _b, t2, t3, C,MU);				\
  op_IEQ_U(4, _c, _d, t4, t5, C,MU);				\
								\
  _a = _mm_shuffle_pd(_a,_a,1);					\
								\
  op_IEQ_D(1, _b, _c, t0, t1, C,MU);				\
  op_IEQ_D(3, _d, _b, t2, t3, C,MU);				\
  op_IEQ_D(5, _c, _d, t4, t5, C,MU);				\

// c = I a
//
//   (u_r a)     = [ u_r a_r, u_r a_i ]
//   (u_r c)     = (u_r I a)
//               = [ -u_r a_i, u_r a_r ]
//   (CI u_i* a)  = [ u_i a_i, -u_i a_r ]
//   (CI u_i* c)  = (u_i  a)
//               = [u_i a_r, u_i a_i ]
//
#define N_KERN_IEQ(t0,t1,t2,t3,t4,t5,C,MU)		\
  _b = _mm_loaddup_pd(u + 0 + C * 2 + 18 * MU);		\
  _b = _mm_mul_pd(_b, _a);				\
  t0 = _mm_add_pd(t0, _b);				\
  _c = _mm_shuffle_pd(_b,_b,1);				\
  t1 = _mm_addsub_pd(t1, _c);				\
  _d = _mm_loaddup_pd(u + 6 + C * 2 + 18 * MU);		\
  _d = _mm_mul_pd(_d, _a);				\
  t2 = _mm_add_pd(t2, _d);				\
  _c = _mm_shuffle_pd(_d,_d,1);				\
  t3 = _mm_addsub_pd(t3, _c);				\
  _b = _mm_loaddup_pd(u + 12 + C * 2 + 18 * MU);	\
  _b = _mm_mul_pd(_b, _a);				\
  t4 = _mm_add_pd(t4, _b);				\
  _c = _mm_shuffle_pd(_b,_b,1);				\
  t5 = _mm_addsub_pd(t5, _c);				\
  							\
  _a = _mm_shuffle_pd(_a,_a,1);				   \
  _d = _mm_loaddup_pd(u + 1 + C * 2 + 18 * MU);		   \
  _c=_d = _mm_mul_pd(_d,_a);				   \
  _d= _mm_xor_pd(_d, _mm_set_pd(-0.0,0.0));		   \
  t0 = _mm_add_pd(t0, _d);				   \
  _c = _mm_shuffle_pd(_c,_c,1);				   \
  t1 = _mm_add_pd(t1, _c);				   \
  _b = _mm_loaddup_pd(u + 7 + C * 2 + 18 * MU);		   \
  _c=_b = _mm_mul_pd(_b, _a);				   \
  _b= _mm_xor_pd(_b, _mm_set_pd(-0.0,0.0));		   \
  t2 = _mm_add_pd(t2, _b);				   \
  _c = _mm_shuffle_pd(_c,_c,1);				   \
  t3 = _mm_add_pd(t3, _c);				   \
  _d = _mm_loaddup_pd(u + 13 + C * 2 + 18 * MU);	   \
  _c=_d = _mm_mul_pd(_d,_a);				   \
  _d= _mm_xor_pd(_d, _mm_set_pd(-0.0,0.0));		   \
  t4 = _mm_add_pd(t4, _d);				   \
  _c = _mm_shuffle_pd(_c,_c,1);				   \
  t5 = _mm_add_pd(t5, _c);				   \

// c = - I a
//
//   (u_r a)     = [ u_r a_r, u_r a_i ]
//   (CI u_i a)  = [ -u_i a_i, u_i a_r ]
//   (u_r c)     = -(u_r I a)
//               = [ u_r a_i,- u_r a_r ]
//   (CI u_i c)  = (u_i  a)
//               = [u_i a_r, u_i a_i ]
//
#define P_KERN_INEG(t0,t1,t2,t3,t4,t5,C,MU)		\
  _b = _mm_loaddup_pd(u + 0 + C * 6 + 18 * MU);		\
  _b = _mm_mul_pd(_b, _a);				\
  t0 = _mm_add_pd(t0, _b);				\
  _c = _mm_shuffle_pd(_b,_b,1);				\
  _c = _mm_xor_pd(_c,_mm_set_pd(-0.0,0.0));             \
  t1 = _mm_add_pd(t1, _c);				\
  _d = _mm_loaddup_pd(u + 2 + C * 6 + 18 * MU);		\
  _d = _mm_mul_pd(_d, _a);				\
  t2 = _mm_add_pd(t2, _d);				\
  _c = _mm_shuffle_pd(_d,_d,1);				\
  _c = _mm_xor_pd(_c,_mm_set_pd(-0.0,0.0));             \
  t3 = _mm_add_pd(t3, _c);				\
  _b = _mm_loaddup_pd(u + 4 + C * 6 + 18 * MU);		\
  _b = _mm_mul_pd(_b, _a);				\
  t4 = _mm_add_pd(t4, _b);				\
  _c = _mm_shuffle_pd(_b,_b,1);				\
  _c = _mm_xor_pd(_c,_mm_set_pd(-0.0,0.0));             \
  t5 = _mm_add_pd(t5, _c);				\
  							\
  _a = _mm_shuffle_pd(_a,_a,1);                         \
  _d = _mm_loaddup_pd(u + 1 + C * 6 + 18 * MU);	        \
  _d = _mm_mul_pd(_d,_a);                               \
  t0 = _mm_addsub_pd(t0, _d);				\
  _c = _mm_shuffle_pd(_d,_d,1);         		\
  t1 = _mm_add_pd(t1, _c);				\
  _b = _mm_loaddup_pd(u + 3 + C * 6 + 18 * MU);	        \
  _b = _mm_mul_pd(_b, _a);				\
  t2 = _mm_addsub_pd(t2, _b);				\
  _c = _mm_shuffle_pd(_b,_b,1);	         		\
  t3 = _mm_add_pd(t3, _c);				\
  _d = _mm_loaddup_pd(u + 5 + C * 6 + 18 * MU);	\
  _d = _mm_mul_pd(_d, _a);				\
  t4 = _mm_addsub_pd(t4, _d);				\
  _c = _mm_shuffle_pd(_d,_d,1);				\
  t5 = _mm_add_pd(t5, _c);				\

// c = - I a
//
//   (u_r a)     = [ u_r a_r, u_r a_i ]
//   (u_r c)     = -(u_r I a)
//               = [ u_r a_i,- u_r a_r ]
//   (CI u_i* a) = [ u_i a_i, -u_i a_r ]
//   (CI u_i* c)  = (-u_i  a)
//               = [-u_i a_r, -u_i a_i ]
//
#define N_KERN_INEG(t0,t1,t2,t3,t4,t5,C,MU)		\
  _b = _mm_loaddup_pd(u + 0 + C * 2 + 18 * MU);		\
  _b = _mm_mul_pd(_b, _a);				\
  t0 = _mm_add_pd(t0, _b);				\
  _c = _mm_shuffle_pd(_b,_b,1);				\
  _c = _mm_xor_pd(_c,_mm_set_pd(-0.0,0.0));             \
  t1 = _mm_add_pd(t1, _c);				\
  _d = _mm_loaddup_pd(u + 6 + C * 2 + 18 * MU);		\
  _d = _mm_mul_pd(_d, _a);				\
  t2 = _mm_add_pd(t2, _d);				\
  _c = _mm_shuffle_pd(_d,_d,1);				\
  _c = _mm_xor_pd(_c,_mm_set_pd(-0.0,0.0));             \
  t3 = _mm_add_pd(t3, _c);				\
  _b = _mm_loaddup_pd(u + 12 + C * 2 + 18 * MU);	\
  _b = _mm_mul_pd(_b, _a);				\
  t4 = _mm_add_pd(t4, _b);				\
  _c = _mm_shuffle_pd(_b,_b,1);				\
  _c = _mm_xor_pd(_c,_mm_set_pd(-0.0,0.0));             \
  t5 = _mm_add_pd(t5, _c);				\
  							\
  _a = _mm_shuffle_pd(_a,_a,1);				\
  _d = _mm_loaddup_pd(u + 1 + C * 2 + 18 * MU);	        \
  _c=_d = _mm_mul_pd(_d, _a);				\
  _d = _mm_xor_pd(_d,_mm_set_pd(-0.0,0.0));		\
  t0 = _mm_add_pd(t0, _d);				\
  _c = _mm_shuffle_pd(_c,_c,1);				\
  t1 = _mm_sub_pd(t1, _c);				\
  _b = _mm_loaddup_pd(u + 7 + C * 2 + 18 * MU);	        \
  _c=_b = _mm_mul_pd(_b, _a);				\
  _b = _mm_xor_pd(_b,_mm_set_pd(-0.0,0.0));		\
  t2 = _mm_add_pd(t2, _b);				\
  _c = _mm_shuffle_pd(_c,_c,1);				\
  t3 = _mm_sub_pd(t3, _c);				\
  _d = _mm_loaddup_pd(u + 13 + C * 2 + 18 * MU);        \
  _c=_d = _mm_mul_pd(_d, _a);				\
  _d = _mm_xor_pd(_d,_mm_set_pd(-0.0,0.0));		\
  t4 = _mm_add_pd(t4, _d);				\
  _c = _mm_shuffle_pd(_c,_c,1);				\
  t5 = _mm_sub_pd(t5, _c);				\


/*  xyztp : _a = _c */
#define P_KERN_EQ(t0,t1,t2,t3,t4,t5,C,MU)               \
   _b = _mm_loaddup_pd(u + 0 + C * 6 + 18 * MU);        \
   _b = _mm_mul_pd(_b, _a);                             \
   t0 = _mm_add_pd(t0, _b);                             \
   t1 = _mm_add_pd(t1, _b);                             \
   _c = _mm_loaddup_pd(u + 2 + C * 6 + 18 * MU);        \
   _c = _mm_mul_pd(_c, _a);                             \
   t2 = _mm_add_pd(t2, _c);                             \
   t3 = _mm_add_pd(t3, _c);                             \
   _d = _mm_loaddup_pd(u + 4 + C * 6 + 18 * MU);        \
   _d = _mm_mul_pd(_d, _a);                             \
   _a = _mm_shuffle_pd(_a, _a, 1);                      \
   t4 = _mm_add_pd(t4, _d);                             \
   t5 = _mm_add_pd(t5, _d);                             \
   _b = _mm_loaddup_pd(u + 1 + C * 6 + 18 * MU);        \
   _b = _mm_mul_pd(_b, _a);                             \
   t0 = _mm_addsub_pd(t0, _b);                          \
   t1 = _mm_addsub_pd(t1, _b);                          \
   _c = _mm_loaddup_pd(u + 3 + C * 6 + 18 * MU);        \
   _c = _mm_mul_pd(_c, _a);                             \
   t2 = _mm_addsub_pd(t2, _c);                          \
   t3 = _mm_addsub_pd(t3, _c);                          \
   _d = _mm_loaddup_pd(u + 5 + C * 6 + 18 * MU);        \
   _d = _mm_mul_pd(_d, _a);                             \
   t4 = _mm_addsub_pd(t4, _d);                          \
   t5 = _mm_addsub_pd(t5, _d);                          \

// c = - a
//
//   (u_r a)     = [ u_r a_r, u_r a_i ]
//   (u_r c)     = -(u_r a)
//               = [- u_r a_r,- u_r a_i ]
//   (CI u_i a) = [- u_i a_i, u_i a_r ]
//   (CI u_i c)  = (-CI u_i  a)
//               = [ u_i a_i, -u_i a_r ]
//
/*  xyztp : _a = - _c */
#define P_KERN_NEG(t0,t1,t2,t3,t4,t5,C,MU)              \
   _b = _mm_loaddup_pd(u + 0 + C * 6 + 18 * MU);        \
   _b = _mm_mul_pd(_b, _a);                             \
   t0 = _mm_add_pd(t0, _b);                             \
   t1 = _mm_sub_pd(t1, _b);                             \
   _b = _mm_loaddup_pd(u + 2 + C * 6 + 18 * MU);        \
   _b = _mm_mul_pd(_b, _a);                             \
   t2 = _mm_add_pd(t2, _b);                             \
   t3 = _mm_sub_pd(t3, _b);                             \
   _d = _mm_loaddup_pd(u + 4 + C * 6 + 18 * MU);        \
   _d = _mm_mul_pd(_d, _a);                             \
   _a = _mm_shuffle_pd(_a, _a, 1);                      \
   t4 = _mm_add_pd(t4, _d);                             \
   t5 = _mm_sub_pd(t5, _d);                             \
							\
   _b = _mm_loaddup_pd(u + 1 + C * 6 + 18 * MU);	\
   _b = _mm_mul_pd(_b, _a);                             \
   t0 = _mm_addsub_pd(t0, _b);                          \
   _b = _mm_xor_pd(_b, _mm_set_pd(-0.0,0.0));		\
   t1 = _mm_add_pd(t1, _b);				\
   _b = _mm_loaddup_pd(u + 3 + C * 6 + 18 * MU);        \
   _b = _mm_mul_pd(_b, _a);                             \
   t2 = _mm_addsub_pd(t2, _b);                          \
   _b = _mm_xor_pd(_b, _mm_set_pd(-0.0,0.0));		\
   t3 = _mm_add_pd(t3, _b);				\
   _b = _mm_loaddup_pd(u + 5 + C * 6 + 18 * MU);        \
   _b = _mm_mul_pd(_b, _a);                             \
   t4 = _mm_addsub_pd(t4, _b);                          \
   _b = _mm_xor_pd(_b, _mm_set_pd(-0.0,0.0));		\
   t5 = _mm_add_pd(t5, _b);				\


#if 1
// ORIGINAL
#define N_HERN(t0,t1,t2,C,MU)				\
  _b  = _mm_loaddup_pd(u + 0 + C * 2 + 18 * MU);	\
  _b = _mm_mul_pd(_b, _a);				\
  t0 = _mm_add_pd(t0, _b);				\
  _c = _mm_loaddup_pd(u + 6 + C * 2 + 18 * MU);		\
  _c = _mm_mul_pd(_c, _a);				\
  t1 = _mm_add_pd(t1, _c);				\
  _d = _mm_loaddup_pd(u + 12 + C * 2 + 18 * MU);	\
  _d = _mm_mul_pd(_d, _a);				\
  t2 = _mm_add_pd(t2, _d);				\
 							\
  _b = _mm_xor_pd(_b, _b);				\
  _a = _mm_sub_pd(_b, _a);				\
  _a = _mm_shuffle_pd(_a, _a, 1);			\
  _c = _mm_loaddup_pd(u + 1 + C * 2 + 18 * MU);		\
  _c = _mm_mul_pd(_c, _a);				\
  t0 = _mm_addsub_pd(t0, _c);				\
  _d = _mm_loaddup_pd(u + 7 + C * 2 + 18 * MU);		\
  _d = _mm_mul_pd(_d, _a);				\
  t1 = _mm_addsub_pd(t1, _d);				\
  _b = _mm_loaddup_pd(u + 13 + C * 2 + 18 * MU);	\
  _b = _mm_mul_pd(_b, _a);				\
  t2 = _mm_addsub_pd(t2, _b);				\

#define N_HERN_noadd(t0,t1,t2,C,MU)				\
  _b = _mm_loaddup_pd(u + 0 + C * 2 + 18 * MU);	\
  t0 = _mm_mul_pd(_b, _a);				\
  _c = _mm_loaddup_pd(u + 6 + C * 2 + 18 * MU);		\
  t1 = _mm_mul_pd(_c, _a);				\
  _d = _mm_loaddup_pd(u + 12 + C * 2 + 18 * MU);	\
  t2= _mm_mul_pd(_d, _a);				\
 							\
  _b = _mm_xor_pd(_b, _b);				\
  _a = _mm_sub_pd(_b, _a);				\
  _a = _mm_shuffle_pd(_a, _a, 1);			\
  _c = _mm_loaddup_pd(u + 1 + C * 2 + 18 * MU);		\
  _c = _mm_mul_pd(_c, _a);				\
  t0 = _mm_addsub_pd(t0, _c);				\
  _d = _mm_loaddup_pd(u + 7 + C * 2 + 18 * MU);		\
  _d = _mm_mul_pd(_d, _a);				\
  t1 = _mm_addsub_pd(t1, _d);				\
  _b = _mm_loaddup_pd(u + 13 + C * 2 + 18 * MU);	\
  _b = _mm_mul_pd(_b, _a);				\
  t2 = _mm_addsub_pd(t2, _b);				\

#endif


#define NH_B0(C,MU,A) _mm_loaddup_pd(u + 0 + C * 2 + 18 * MU)
#define NH_B1(C,MU,A) _mm_mul_pd(NH_B0(C,MU,A),A) 
#define NH_C0(C,MU,A) _mm_loaddup_pd(u + 6 + C * 2 + 18 * MU)
#define NH_C1(C,MU,A) _mm_mul_pd(NH_C0(C,MU,A), A)
#define NH_D0(C,MU,A) _mm_loaddup_pd(u + 12 + C * 2 + 18 * MU)
#define NH_D1(C,MU,A) _mm_mul_pd(NH_D0(C,MU,A), A)
#define NH_B2(C,MU,A) _mm_xor_pd(NH_B1(C,MU,A),NH_B1(C,MU,A))
#define NH_A0(C,MU,A) _mm_sub_pd(NH_B2(C,MU,A), A)
#define NH_A1(C,MU,A) _mm_shuffle_pd(NH_A0(C,MU,A),NH_A0(C,MU,A), 1)
#define NH_C2(C,MU,A) _mm_loaddup_pd(u + 1 + C * 2 + 18 * MU)
#define NH_C3(C,MU,A) _mm_mul_pd(NH_C2(C,MU,A), NH_A1(C,MU,A))
#define NH_D2(C,MU,A) _mm_loaddup_pd(u + 7 + C * 2 + 18 * MU)
#define NH_D3(C,MU,A) _mm_mul_pd(NH_D2(C,MU,A), NH_A1(C,MU,A))
#define NH_B3(C,MU,A) _mm_loaddup_pd(u + 13 + C * 2 + 18 * MU)
#define NH_B4(C,MU,A) _mm_mul_pd(NH_B3(C,MU,A), NH_A1(C,MU,A))

#define N_HERN2(t0,t1,t2,C,MU,A)			\
  t0 = _mm_add_pd(t0, NH_B1(C,MU,A)) ;			\
  t1 = _mm_add_pd(t1, NH_C1(C,MU,A));			\
  t2 = _mm_add_pd(t2, NH_D1(C,MU,A));			\
							\
  t0 = _mm_addsub_pd(t0, NH_C3(C,MU,A));		\
  t1 = _mm_addsub_pd(t1, NH_D3(C,MU,A));			\
  t2 = _mm_addsub_pd(t2, NH_B4(C,MU,A));			\





/*  xyztm */
#define N_KERN(t0,t1,t2,t3,t4,t5,C,MU)                  \
   _d = _b = _mm_loaddup_pd(u + 0 + C * 2 + 18 * MU);   \
   _b = _mm_mul_pd(_b, _a);                             \
   t0 = _mm_add_pd(t0, _b);                             \
   _d = _mm_mul_pd(_d, _c);                             \
   t1 = _mm_add_pd(t1, _d);                             \
   _d = _b = _mm_loaddup_pd(u + 6 + C * 2 + 18 * MU);   \
   _b = _mm_mul_pd(_b, _a);                             \
   t2 = _mm_add_pd(t2, _b);                             \
   _d = _mm_mul_pd(_d, _c);                             \
   t3 = _mm_add_pd(t3, _d);                             \
   _d = _b = _mm_loaddup_pd(u + 12 + C * 2 + 18 * MU);  \
   _b = _mm_mul_pd(_b, _a);                             \
   t4 = _mm_add_pd(t4, _b);                             \
   _d = _mm_mul_pd(_d, _c);                             \
   t5 = _mm_add_pd(t5, _d);                             \
   _b = _mm_xor_pd(_b, _b);                             \
   _a = _b = _mm_sub_pd(_b, _a);                        \
   _b = _mm_xor_pd(_b, _b);                             \
   _c = _b = _mm_sub_pd(_b, _c);                        \
   _a = _mm_shuffle_pd(_a, _a, 1);                      \
   _c = _mm_shuffle_pd(_c, _c, 1);                      \
   _d = _b = _mm_loaddup_pd(u + 1 + C * 2 + 18 * MU);   \
   _b = _mm_mul_pd(_b, _a);                             \
   t0 = _mm_addsub_pd(t0, _b);                          \
   _d = _mm_mul_pd(_d, _c);                             \
   t1 = _mm_addsub_pd(t1, _d);                          \
   _d = _b = _mm_loaddup_pd(u + 7 + C * 2 + 18 * MU);   \
   _b = _mm_mul_pd(_b, _a);                             \
   t2 = _mm_addsub_pd(t2, _b);                          \
   _d = _mm_mul_pd(_d, _c);                             \
   t3 = _mm_addsub_pd(t3, _d);                          \
   _d = _b = _mm_loaddup_pd(u + 13 + C * 2 + 18 * MU);  \
   _b = _mm_mul_pd(_b, _a);                             \
   t4 = _mm_addsub_pd(t4, _b);                          \
   _d = _mm_mul_pd(_d, _c);                             \
   t5 = _mm_addsub_pd(t5, _d);                          \

/*  xyztm : _a = _c */
#define N_KERN_EQ(t0,t1,t2,t3,t4,t5,C,MU)               \
   _b = _mm_loaddup_pd(u + 0 + C * 2 + 18 * MU);        \
   _b = _mm_mul_pd(_b, _a);                             \
   t0 = _mm_add_pd(t0, _b);                             \
   t1 = _mm_add_pd(t1, _b);                             \
   _c = _mm_loaddup_pd(u + 6 + C * 2 + 18 * MU);        \
   _c = _mm_mul_pd(_c, _a);                             \
   t2 = _mm_add_pd(t2, _c);                             \
   t3 = _mm_add_pd(t3, _c);                             \
   _b = _a;                                             \
   _b = _mm_shuffle_pd(_b, _b, 1);                      \
   _d = _mm_loaddup_pd(u + 12 + C * 2 + 18 * MU);       \
   _d = _mm_mul_pd(_d, _a);                             \
   t4 = _mm_add_pd(t4, _d);                             \
   t5 = _mm_add_pd(t5, _d);                             \
   _a = _mm_xor_pd(_a, _a);                             \
   _a = _mm_sub_pd(_a, _b);                             \
   _b = _mm_loaddup_pd(u + 1 + C * 2 + 18 * MU);        \
   _b = _mm_mul_pd(_b, _a);                             \
   t0 = _mm_addsub_pd(t0, _b);                          \
   t1 = _mm_addsub_pd(t1, _b);                          \
   _c = _mm_loaddup_pd(u + 7 + C * 2 + 18 * MU);        \
   _c = _mm_mul_pd(_c, _a);                             \
   t2 = _mm_addsub_pd(t2, _c);                          \
   t3 = _mm_addsub_pd(t3, _c);                          \
   _d = _mm_loaddup_pd(u + 13 + C * 2 + 18 * MU);       \
   _d = _mm_mul_pd(_d, _a);                             \
   t4 = _mm_addsub_pd(t4, _d);                          \
   t5 = _mm_addsub_pd(t5, _d);                          \

// c = - a
//
//   (u_r a)     = [ u_r a_r, u_r a_i ]
//   (u_r c)     = -(u_r a)
//               = [- u_r a_r,- u_r a_i ]
//   (CI u_i* a) = [ u_i a_i, -u_i a_r ]
//   (CI u_i* c)  = (CI u_i  a)
//               = [ -u_i a_i, u_i a_r ]
//
/*  xyztm : _a = - _c */
#define N_KERN_NEG(t0,t1,t2,t3,t4,t5,C,MU)              \
   _b = _mm_loaddup_pd(u + 0 + C * 2 + 18 * MU);        \
   _b = _mm_mul_pd(_b, _a);                             \
   t0 = _mm_add_pd(t0, _b);                             \
   t1 = _mm_sub_pd(t1, _b);                             \
   _b = _mm_loaddup_pd(u + 6 + C * 2 + 18 * MU);        \
   _b = _mm_mul_pd(_b, _a);                             \
   t2 = _mm_add_pd(t2, _b);                             \
   t3 = _mm_sub_pd(t3, _b);                             \
   _d = _mm_loaddup_pd(u + 12 + C * 2 + 18 * MU);       \
   _d = _mm_mul_pd(_d, _a);                             \
   _a = _mm_shuffle_pd(_a, _a, 1);                      \
   t4 = _mm_add_pd(t4, _d);                             \
   t5 = _mm_sub_pd(t5, _d);                             \
							\
   _b = _mm_loaddup_pd(u + 1 + C * 2 + 18 * MU);	\
   _d = _b = _mm_mul_pd(_b, _a);			\
   _b = _mm_xor_pd(_b, _mm_set_pd(-0.0,0.0));		\
   t0 = _mm_add_pd(t0, _b);				\
   t1 = _mm_addsub_pd(t1, _d);                          \
   _b = _mm_loaddup_pd(u + 7 + C * 2 + 18 * MU);	\
   _d =_b = _mm_mul_pd(_b, _a);				\
   _b = _mm_xor_pd(_b, _mm_set_pd(-0.0,0.0));		\
   t2 = _mm_add_pd(t2, _b);				\
   t3 = _mm_addsub_pd(t3, _d);                          \
   _d = _mm_loaddup_pd(u + 13 + C * 2 + 18 * MU);       \
   _d = _b = _mm_mul_pd(_d, _a);			\
   _b = _mm_xor_pd(_b, _mm_set_pd(-0.0,0.0));		\
   t4 = _mm_add_pd(t4, _b);				\
   t5 = _mm_addsub_pd(t5, _d);                          \


#ifdef USE_MACROS
#include "sse-macros-dag1.h"
#include "sse-macros-dag0.h"
#else
#include "sse-inlines.h"
#endif


