#if 1
#define N_KERN_XP				    \
  N_KERN_XP_03(0);				    \
  P_HERN(wxp[0],wxp[1],wxp[2],0,0);		    \
  N_KERN_XP_03(1);				    \
  P_HERN(wxp[0],wxp[1],wxp[2],1,0);		    \
  N_KERN_XP_03(2);				    \
  P_HERN(wxp[0],wxp[1],wxp[2],2,0);		    \
  N_KERN_XP_12(0);				    \
  P_HERN(wxp[3],wxp[4],wxp[5],0,0);		    \
  N_KERN_XP_12(1);				    \
  P_HERN(wxp[3],wxp[4],wxp[5],1,0);		    \
  N_KERN_XP_12(2);				    \
  P_HERN(wxp[3],wxp[4],wxp[5],2,0);		    \
  
#define N_KERN_XP_EDG(D_03, D_12)		    \
  _a = _mm_load_pd( (D_03) );			    \
  P_HERN(wxp[0],wxp[1],wxp[2],0,0);		    \
  _a = _mm_load_pd( (D_03)+2 );			    \
  P_HERN(wxp[0],wxp[1],wxp[2],1,0);		    \
  _a = _mm_load_pd( (D_03)+4 );			    \
  P_HERN(wxp[0],wxp[1],wxp[2],2,0);		    \
  _a = _mm_load_pd( (D_12) );			    \
  P_HERN(wxp[3],wxp[4],wxp[5],0,0);		    \
  _a = _mm_load_pd( (D_12)+2 );			    \
  P_HERN(wxp[3],wxp[4],wxp[5],1,0);		    \
  _a = _mm_load_pd( (D_12)+4 );			    \
  P_HERN(wxp[3],wxp[4],wxp[5],2,0);		    \

#else
#define N_KERN_XP				    \
  N_KERN_XP_03(0);				    \
  P_KERN_INEG(t00,t09,t01,t10,t02,t11,0,0);	    \
  N_KERN_XP_03(1);				    \
  P_KERN_INEG(t00,t09,t01,t10,t02,t11,1,0);	    \
  N_KERN_XP_03(2);				    \
  P_KERN_INEG(t00,t09,t01,t10,t02,t11,2,0);	    \
  N_KERN_XP_12(0);				    \
  P_KERN_INEG(t03,t06,t04,t07,t05,t08,0,0);	    \
  N_KERN_XP_12(1);				    \
  P_KERN_INEG(t03,t06,t04,t07,t05,t08,1,0);	    \
  N_KERN_XP_12(2);				    \
  P_KERN_INEG(t03,t06,t04,t07,t05,t08,2,0);	    \

#define N_KERN_XP_EDG(D_03, D_12)		    \
  _a = _mm_load_pd( (D_03) );			    \
  P_KERN_INEG(t00,t09,t01,t10,t02,t11,0,0);	    \
  _a = _mm_load_pd( (D_03)+2 );			    \
  P_KERN_INEG(t00,t09,t01,t10,t02,t11,1,0);	    \
  _a = _mm_load_pd( (D_03)+4 );			    \
  P_KERN_INEG(t00,t09,t01,t10,t02,t11,2,0);	    \
  _a = _mm_load_pd( (D_12) );			    \
  P_KERN_INEG(t03,t06,t04,t07,t05,t08,0,0);	    \
  _a = _mm_load_pd( (D_12)+2 );			    \
  P_KERN_INEG(t03,t06,t04,t07,t05,t08,1,0);	    \
  _a = _mm_load_pd( (D_12)+4 );			    \
  P_KERN_INEG(t03,t06,t04,t07,t05,t08,2,0);	    \

#endif

#define N_STORE_XP \
  _mm_store_pd(chi    ,  _mm_add_pd(*(__m128d*)(chi+0),wxp[0]));	\
  _mm_store_pd(chi+2  ,  _mm_add_pd(*(__m128d*)(chi+2),wxp[1]));	\
  _mm_store_pd(chi+4  ,  _mm_add_pd(*(__m128d*)(chi+4),wxp[2]));	\
  _mm_store_pd(chi+6  ,  _mm_add_pd(*(__m128d*)(chi+6),wxp[3]));	\
  _mm_store_pd(chi+8  ,  _mm_add_pd(*(__m128d*)(chi+8),wxp[4]));	\
  _mm_store_pd(chi+10 ,  _mm_add_pd(*(__m128d*)(chi+10),wxp[5]));	\
  									\
  _mm_store_pd(chi+12, _mm_add_pd(*(__m128d*)(chi+12),			\
   _mm_xor_pd(_mm_shuffle_pd(wxp[3],wxp[3],1), _mm_set_pd(-0.0,0.0)))); \
  _mm_store_pd(chi+14, _mm_add_pd(*(__m128d*)(chi+14),			\
   _mm_xor_pd(_mm_shuffle_pd(wxp[4],wxp[5],1), _mm_set_pd(-0.0,0.0)))); \
  _mm_store_pd(chi+16, _mm_add_pd(*(__m128d*)(chi+16),			\
   _mm_xor_pd(_mm_shuffle_pd(wxp[5],wxp[5],1), _mm_set_pd(-0.0,0.0)))); \
  _mm_store_pd(chi+18, _mm_add_pd(*(__m128d*)(chi+18),			\
   _mm_xor_pd(_mm_shuffle_pd(wxp[0],wxp[3],1), _mm_set_pd(-0.0,0.0)))); \
  _mm_store_pd(chi+20, _mm_add_pd(*(__m128d*)(chi+20),			\
   _mm_xor_pd(_mm_shuffle_pd(wxp[1],wxp[5],1), _mm_set_pd(-0.0,0.0)))); \
  _mm_store_pd(chi+22, _mm_add_pd(*(__m128d*)(chi+22),			\
   _mm_xor_pd(_mm_shuffle_pd(wxp[2],wxp[5],1), _mm_set_pd(-0.0,0.0)))); \

#define N_STORE_XP_noadd \
  _mm_store_pd(chi    ,  wxp[0]);				\
  _mm_store_pd(chi+2  ,  wxp[1]);				\
  _mm_store_pd(chi+4  ,  wxp[2]);				\
  _mm_store_pd(chi+6  ,  wxp[3]);				\
  _mm_store_pd(chi+8  ,  wxp[4]);				\
  _mm_store_pd(chi+10 ,  wxp[5]);				\
  									\
  _mm_store_pd(chi+12,							\
   _mm_xor_pd(_mm_shuffle_pd(wxp[3],wxp[3],1), _mm_set_pd(-0.0,0.0)));  \
  _mm_store_pd(chi+14,							\
   _mm_xor_pd(_mm_shuffle_pd(wxp[4],wxp[4],1), _mm_set_pd(-0.0,0.0)));  \
  _mm_store_pd(chi+16,							\
   _mm_xor_pd(_mm_shuffle_pd(wxp[5],wxp[5],1), _mm_set_pd(-0.0,0.0))) ; \
  _mm_store_pd(chi+18,							\
   _mm_xor_pd(_mm_shuffle_pd(wxp[0],wxp[0],1), _mm_set_pd(-0.0,0.0)));  \
  _mm_store_pd(chi+20,							\
   _mm_xor_pd(_mm_shuffle_pd(wxp[1],wxp[1],1), _mm_set_pd(-0.0,0.0)));  \
  _mm_store_pd(chi+22,							\
   _mm_xor_pd(_mm_shuffle_pd(wxp[2],wxp[2],1), _mm_set_pd(-0.0,0.0))) ; \

#if 1
#define N_KERN_YP                              \
  N_KERN_YP_03(0);			       \
  P_HERN(wyp[0],wyp[1],wyp[2],0,1);	       \
  N_KERN_YP_03(1);			       \
  P_HERN(wyp[0],wyp[1],wyp[2],1,1);	       \
  N_KERN_YP_03(2);			       \
  P_HERN(wyp[0],wyp[1],wyp[2],2,1);	       \
  N_KERN_YP_12(0);			       \
  P_HERN(wyp[3],wyp[4],wyp[5],0,1);	       \
  N_KERN_YP_12(1);			       \
  P_HERN(wyp[3],wyp[4],wyp[5],1,1);	       \
  N_KERN_YP_12(2);			       \
  P_HERN(wyp[3],wyp[4],wyp[5],2,1);	       \

#define N_KERN_YP_EDG(D_03,D_12)				    \
  _a = _mm_load_pd((D_03));					    \
  P_HERN(wyp[0],wyp[1],wyp[2],0,1);				    \
  _a = _mm_load_pd((D_03)+2);					    \
  P_HERN(wyp[0],wyp[1],wyp[2],1,1);				    \
  _a = _mm_load_pd((D_03)+4);					    \
  P_HERN(wyp[0],wyp[1],wyp[2],2,1);				    \
  _a = _mm_load_pd((D_12));					    \
  P_HERN(wyp[3],wyp[4],wyp[5],0,1);				    \
  _a = _mm_load_pd((D_12)+2);					    \
  P_HERN(wyp[3],wyp[4],wyp[5],1,1);				    \
  _a = _mm_load_pd((D_12)+4);					    \
  P_HERN(wyp[3],wyp[4],wyp[5],2,1);				    \
  
//
//  2010-04-07 T.I
// 
// The above two lines has a bug 1,2 at the fourth argument were 0,0
//  I should be hung....
//


#else
#define N_KERN_YP                              \
   N_KERN_YP_03(0);                            \
   P_KERN_NEG(t00,t09,t01,t10,t02,t11,0,1);    \
   N_KERN_YP_03(1);                            \
   P_KERN_NEG(t00,t09,t01,t10,t02,t11,1,1);    \
   N_KERN_YP_03(2);                            \
   P_KERN_NEG(t00,t09,t01,t10,t02,t11,2,1);    \
   N_KERN_YP_12(0);                            \
   P_KERN_EQ(t03,t06,t04,t07,t05,t08,0,1);     \
   N_KERN_YP_12(1);                            \
   P_KERN_EQ(t03,t06,t04,t07,t05,t08,1,1);     \
   N_KERN_YP_12(2);                            \
   P_KERN_EQ(t03,t06,t04,t07,t05,t08,2,1);     \

#define N_KERN_YP_EDG(D_03,D_12)	       \
  _a = _mm_load_pd((D_03));			\
  P_KERN_NEG(t00,t09,t01,t10,t02,t11,0,1);	\
  _a = _mm_load_pd((D_03)+2);			\
  P_KERN_NEG(t00,t09,t01,t10,t02,t11,1,1);	\
  _a = _mm_load_pd((D_03)+4);			\
  P_KERN_NEG(t00,t09,t01,t10,t02,t11,2,1);	\
  _a = _mm_load_pd((D_12));			\
   P_KERN_EQ(t03,t06,t04,t07,t05,t08,0,1);     \
  _a = _mm_load_pd((D_12)+2);			\
   P_KERN_EQ(t03,t06,t04,t07,t05,t08,1,1);     \
  _a = _mm_load_pd((D_12)+4);			\
   P_KERN_EQ(t03,t06,t04,t07,t05,t08,2,1);     \

#endif

#define N_STORE_YP							\
  _mm_store_pd(chi    ,  _mm_add_pd(*(__m128d*)(chi+0),wyp[0]));	\
  _mm_store_pd(chi+2  ,  _mm_add_pd(*(__m128d*)(chi+2),wyp[1]));	\
  _mm_store_pd(chi+4  ,  _mm_add_pd(*(__m128d*)(chi+4),wyp[2]));	\
  _mm_store_pd(chi+6  ,  _mm_add_pd(*(__m128d*)(chi+6),wyp[3]));	\
  _mm_store_pd(chi+8  ,  _mm_add_pd(*(__m128d*)(chi+8),wyp[4]));	\
  _mm_store_pd(chi+10 ,  _mm_add_pd(*(__m128d*)(chi+10),wyp[5]));	\
  									\
  _mm_store_pd(chi+12, _mm_add_pd(*(__m128d*)(chi+12),	wyp[3]));	\
  _mm_store_pd(chi+14, _mm_add_pd(*(__m128d*)(chi+14),	wyp[4]));	\
  _mm_store_pd(chi+16, _mm_add_pd(*(__m128d*)(chi+16),	wyp[5]));	\
  _mm_store_pd(chi+18, _mm_sub_pd(*(__m128d*)(chi+18),	wyp[0]));	\
  _mm_store_pd(chi+20, _mm_sub_pd(*(__m128d*)(chi+20),	wyp[1]));	\
  _mm_store_pd(chi+22, _mm_sub_pd(*(__m128d*)(chi+22),	wyp[2]));	\



#if 1
#define N_KERN_ZP							\
  N_KERN_ZP_02(0);							\
  P_HERN(wzp[0],wzp[1],wzp[2],0,2);					\
  N_KERN_ZP_02(1);							\
  P_HERN(wzp[0],wzp[1],wzp[2],1,2);					\
  N_KERN_ZP_02(2);							\
  P_HERN(wzp[0],wzp[1],wzp[2],2,2);					\
  N_KERN_ZP_13(0);							\
  P_HERN(wzp[3],wzp[4],wzp[5],0,2);					\
  N_KERN_ZP_13(1);							\
  P_HERN(wzp[3],wzp[4],wzp[5],1,2);					\
  N_KERN_ZP_13(2);							\
  P_HERN(wzp[3],wzp[4],wzp[5],2,2);					\
  
#define N_KERN_ZP_EDG( DADDR_02, DADDR_13)	       \
  _a = _mm_load_pd((DADDR_02));			       \
  P_HERN(wzp[0],wzp[1],wzp[2],0,2);		       \
  _a = _mm_load_pd((DADDR_02)+2);		       \
  P_HERN(wzp[0],wzp[1],wzp[2],1,2);		       \
  _a = _mm_load_pd((DADDR_02)+4);		       \
  P_HERN(wzp[0],wzp[1],wzp[2],2,2);		       \
  _a = _mm_load_pd((DADDR_13));			       \
  P_HERN(wzp[3],wzp[4],wzp[5],0,2);		       \
  _a = _mm_load_pd((DADDR_13)+2);		       \
  P_HERN(wzp[3],wzp[4],wzp[5],1,2);		       \
  _a = _mm_load_pd((DADDR_13)+4);		       \
  P_HERN(wzp[3],wzp[4],wzp[5],2,2);		       \
  
#else
#define N_KERN_ZP                              \
   N_KERN_ZP_02(0);                            \
   P_KERN_INEG(t00,t06,t01,t07,t02,t08,0,2);				\
   N_KERN_ZP_02(1);                            \
   P_KERN_INEG(t00,t06,t01,t07,t02,t08,1,2);        \
   N_KERN_ZP_02(2);                            \
   P_KERN_INEG(t00,t06,t01,t07,t02,t08,2,2);        \
   N_KERN_ZP_13(0);                            \
   P_KERN_IEQ(t03,t09,t04,t10,t05,t11,0,2);        \
   N_KERN_ZP_13(1);                            \
   P_KERN_IEQ(t03,t09,t04,t10,t05,t11,1,2);        \
   N_KERN_ZP_13(2);                            \
   P_KERN_IEQ(t03,t09,t04,t10,t05,t11,2,2);        \

#define N_KERN_ZP_EDG(D_02,D_13)	       \
  _a = _mm_load_pd((D_02));			    \
 P_KERN_INEG(t00,t06,t01,t07,t02,t08,0,2);	    \
  _a = _mm_load_pd((D_02)+2);			    \
  P_KERN_INEG(t00,t06,t01,t07,t02,t08,1,2);	    \
  _a = _mm_load_pd((D_02)+4);			    \
  P_KERN_INEG(t00,t06,t01,t07,t02,t08,2,2);	    \
  _a = _mm_load_pd((D_13));			    \
  P_KERN_IEQ(t03,t09,t04,t10,t05,t11,0,2);	    \
  _a = _mm_load_pd((D_13)+2);			    \
  P_KERN_IEQ(t03,t09,t04,t10,t05,t11,1,2);	    \
  _a = _mm_load_pd((D_13)+4);			    \
  P_KERN_IEQ(t03,t09,t04,t10,t05,t11,2,2);	    \

#endif
#define N_STORE_ZP							\
  _mm_store_pd(chi    ,  _mm_add_pd(*(__m128d*)(chi+0),wzp[0]));	\
  _mm_store_pd(chi+2  ,  _mm_add_pd(*(__m128d*)(chi+2),wzp[1]));	\
  _mm_store_pd(chi+4  ,  _mm_add_pd(*(__m128d*)(chi+4),wzp[2]));	\
  _mm_store_pd(chi+6  ,  _mm_add_pd(*(__m128d*)(chi+6),wzp[3]));	\
  _mm_store_pd(chi+8  ,  _mm_add_pd(*(__m128d*)(chi+8),wzp[4]));	\
  _mm_store_pd(chi+10 ,  _mm_add_pd(*(__m128d*)(chi+10),wzp[5]));	\
  									\
  _mm_store_pd(chi+12, _mm_add_pd(*(__m128d*)(chi+12),_mm_xor_pd(_mm_shuffle_pd(wzp[0],wzp[0],1),_mm_set_pd(-0.0,0.)))); \
  _mm_store_pd(chi+14, _mm_add_pd(*(__m128d*)(chi+14),_mm_xor_pd(_mm_shuffle_pd(wzp[1],wzp[1],1),_mm_set_pd(-0.0,0.)))); \
  _mm_store_pd(chi+16, _mm_add_pd(*(__m128d*)(chi+16),_mm_xor_pd(_mm_shuffle_pd(wzp[2],wzp[2],1),_mm_set_pd(-0.0,0.)))); \
  _mm_store_pd(chi+18, _mm_addsub_pd(*(__m128d*)(chi+18),_mm_shuffle_pd(wzp[3],wzp[3],1))); \
  _mm_store_pd(chi+20, _mm_addsub_pd(*(__m128d*)(chi+20),_mm_shuffle_pd(wzp[4],wzp[4],1))); \
  _mm_store_pd(chi+22, _mm_addsub_pd(*(__m128d*)(chi+22),_mm_shuffle_pd(wzp[5],wzp[5],1))); \


#if 1
#define N_KERN_TP                              \
   N_KERN_TP_02(0);                            \
   P_HERN(wtp[0],wtp[1],wtp[2],0,3);	       \
   N_KERN_TP_02(1);                            \
   P_HERN(wtp[0],wtp[1],wtp[2],1,3);	       \
   N_KERN_TP_02(2);                            \
   P_HERN(wtp[0],wtp[1],wtp[2],2,3);	       \
   N_KERN_TP_13(0);                            \
   P_HERN(wtp[3],wtp[4],wtp[5],0,3);	       \
   N_KERN_TP_13(1);                            \
   P_HERN(wtp[3],wtp[4],wtp[5],1,3);	       \
   N_KERN_TP_13(2);                            \
   P_HERN(wtp[3],wtp[4],wtp[5],2,3);	       \


#define N_KERN_TP_EDG( DADDR_02, DADDR_13)     \
  _a = _mm_load_pd( (DADDR_02) );	       \
  P_HERN(wtp[0],wtp[1],wtp[2],0,3);	       \
  _a = _mm_load_pd( (DADDR_02) +2);	       \
  P_HERN(wtp[0],wtp[1],wtp[2],1,3);	       \
  _a = _mm_load_pd( (DADDR_02)+4 );	       \
  P_HERN(wtp[0],wtp[1],wtp[2],2,3);	       \
  _a = _mm_load_pd( (DADDR_13) );	       \
  P_HERN(wtp[3],wtp[4],wtp[5],0,3);	       \
  _a = _mm_load_pd( (DADDR_13)+2 );	       \
  P_HERN(wtp[3],wtp[4],wtp[5],1,3);	       \
  _a = _mm_load_pd( (DADDR_13)+4 );	       \
  P_HERN(wtp[3],wtp[4],wtp[5],2,3);	       \
   
#else
#define N_KERN_TP                              \
   N_KERN_TP_02(0);                            \
   P_KERN_EQ(t00,t06,t01,t07,t02,t08,0,3);     \
   N_KERN_TP_02(1);                            \
   P_KERN_EQ(t00,t06,t01,t07,t02,t08,1,3);     \
   N_KERN_TP_02(2);                            \
   P_KERN_EQ(t00,t06,t01,t07,t02,t08,2,3);     \
   N_KERN_TP_13(0);                            \
   P_KERN_EQ(t03,t09,t04,t10,t05,t11,0,3);     \
   N_KERN_TP_13(1);                            \
   P_KERN_EQ(t03,t09,t04,t10,t05,t11,1,3);     \
   N_KERN_TP_13(2);                            \
   P_KERN_EQ(t03,t09,t04,t10,t05,t11,2,3);     \

#define N_KERN_TP_EDG(D_02, D_13)		\
  _a = _mm_load_pd( (D_02) );			\
   P_KERN_EQ(t00,t06,t01,t07,t02,t08,0,3);     \
  _a = _mm_load_pd( (D_02)+2 );			\
   P_KERN_EQ(t00,t06,t01,t07,t02,t08,1,3);     \
  _a = _mm_load_pd( (D_02)+4 );			\
   P_KERN_EQ(t00,t06,t01,t07,t02,t08,2,3);     \
  _a = _mm_load_pd( (D_13) );			\
   P_KERN_EQ(t03,t09,t04,t10,t05,t11,0,3);     \
  _a = _mm_load_pd( (D_13)+2 );			\
   P_KERN_EQ(t03,t09,t04,t10,t05,t11,1,3);     \
   _a = _mm_load_pd( (D_13)+4 );		       \
   P_KERN_EQ(t03,t09,t04,t10,t05,t11,2,3);     \

#endif
#define N_STORE_TP							\
  _mm_store_pd(chi    ,  _mm_add_pd(*(__m128d*)(chi+0),wtp[0]));	\
  _mm_store_pd(chi+2  ,  _mm_add_pd(*(__m128d*)(chi+2),wtp[1]));	\
  _mm_store_pd(chi+4  ,  _mm_add_pd(*(__m128d*)(chi+4),wtp[2]));	\
  _mm_store_pd(chi+6  ,  _mm_add_pd(*(__m128d*)(chi+6),wtp[3]));	\
  _mm_store_pd(chi+8  ,  _mm_add_pd(*(__m128d*)(chi+8),wtp[4]));	\
  _mm_store_pd(chi+10 ,  _mm_add_pd(*(__m128d*)(chi+10),wtp[5]));	\
  									\
  _mm_store_pd(chi+12, _mm_add_pd(*(__m128d*)(chi+12),wtp[0]));		\
  _mm_store_pd(chi+14, _mm_add_pd(*(__m128d*)(chi+14),wtp[1]));		\
  _mm_store_pd(chi+16, _mm_add_pd(*(__m128d*)(chi+16),wtp[2]));		\
  _mm_store_pd(chi+18, _mm_add_pd(*(__m128d*)(chi+18),wtp[3]));		\
  _mm_store_pd(chi+20, _mm_add_pd(*(__m128d*)(chi+20),wtp[4]));		\
  _mm_store_pd(chi+22, _mm_add_pd(*(__m128d*)(chi+22),wtp[5]));		\

#if 1
#define N_KERN_XM                              \
   P_KERN_XP_03(0);                            \
   N_HERN(wxm[0],wxm[1],wxm[2],0,0);	       \
   P_KERN_XP_03(1);                            \
   N_HERN(wxm[0],wxm[1],wxm[2],1,0);	       \
   P_KERN_XP_03(2);                            \
   N_HERN(wxm[0],wxm[1],wxm[2],2,0);	       \
   P_KERN_XP_12(0);                            \
   N_HERN(wxm[3],wxm[4],wxm[5],0,0);	       \
   P_KERN_XP_12(1);                            \
   N_HERN(wxm[3],wxm[4],wxm[5],1,0);	       \
   P_KERN_XP_12(2);                            \
   N_HERN(wxm[3],wxm[4],wxm[5],2,0);	       \

#else
#define N_KERN_XM				   \
  P_KERN_XP_03(0);				   \
  N_KERN_IEQ(t00,t09,t01,t10,t02,t11,0,0);	   \
  P_KERN_XP_03(1);				   \
  N_KERN_IEQ(t00,t09,t01,t10,t02,t11,1,0);	   \
  P_KERN_XP_03(2);				   \
  N_KERN_IEQ(t00,t09,t01,t10,t02,t11,2,0);	   \
  P_KERN_XP_12(0);				   \
  N_KERN_IEQ(t03,t06,t04,t07,t05,t08,0,0);	   \
  P_KERN_XP_12(1);				   \
  N_KERN_IEQ(t03,t06,t04,t07,t05,t08,1,0);	   \
  P_KERN_XP_12(2);				   \
  N_KERN_IEQ(t03,t06,t04,t07,t05,t08,2,0);	   \

#endif

#define N_STORE_XM							\
  _mm_store_pd(chi    ,  _mm_add_pd(*(__m128d*)(chi+0),wxm[0]));	\
  _mm_store_pd(chi+2  ,  _mm_add_pd(*(__m128d*)(chi+2),wxm[1]));	\
  _mm_store_pd(chi+4  ,  _mm_add_pd(*(__m128d*)(chi+4),wxm[2]));	\
  _mm_store_pd(chi+6  ,  _mm_add_pd(*(__m128d*)(chi+6),wxm[3]));	\
  _mm_store_pd(chi+8  ,  _mm_add_pd(*(__m128d*)(chi+8),wxm[4]));	\
  _mm_store_pd(chi+10 ,  _mm_add_pd(*(__m128d*)(chi+10),wxm[5]));	\
  									\
  _mm_store_pd(chi+12, _mm_addsub_pd(*(__m128d*)(chi+12),_mm_shuffle_pd(wxm[3],wxm[3],1))); \
  _mm_store_pd(chi+14, _mm_addsub_pd(*(__m128d*)(chi+14),_mm_shuffle_pd(wxm[4],wxm[4],1))); \
  _mm_store_pd(chi+16, _mm_addsub_pd(*(__m128d*)(chi+16),_mm_shuffle_pd(wxm[5],wxm[5],1))); \
  _mm_store_pd(chi+18, _mm_addsub_pd(*(__m128d*)(chi+18),_mm_shuffle_pd(wxm[0],wxm[0],1))); \
  _mm_store_pd(chi+20, _mm_addsub_pd(*(__m128d*)(chi+20),_mm_shuffle_pd(wxm[1],wxm[1],1))); \
  _mm_store_pd(chi+22, _mm_addsub_pd(*(__m128d*)(chi+22),_mm_shuffle_pd(wxm[2],wxm[2],1))); \
  


#if 1
#define N_KERN_YM                              \
   P_KERN_YP_03(0);                            \
   N_HERN(wym[0],wym[1],wym[2],0,1);	       \
   P_KERN_YP_03(1);                            \
   N_HERN(wym[0],wym[1],wym[2],1,1);	       \
   P_KERN_YP_03(2);                            \
   N_HERN(wym[0],wym[1],wym[2],2,1);	       \
   P_KERN_YP_12(0);                            \
   N_HERN(wym[3],wym[4],wym[5],0,1);	       \
   P_KERN_YP_12(1);                            \
   N_HERN(wym[3],wym[4],wym[5],1,1);	       \
   P_KERN_YP_12(2);                            \
   N_HERN(wym[3],wym[4],wym[5],2,1);	       \


#else
#define N_KERN_YM                              \
   P_KERN_YP_03(0);                            \
   N_KERN_EQ(t00,t09,t01,t10,t02,t11,0,1);     \
   P_KERN_YP_03(1);                            \
   N_KERN_EQ(t00,t09,t01,t10,t02,t11,1,1);     \
   P_KERN_YP_03(2);                            \
   N_KERN_EQ(t00,t09,t01,t10,t02,t11,2,1);     \
   P_KERN_YP_12(0);                            \
   N_KERN_NEG(t03,t06,t04,t07,t05,t08,0,1);    \
   P_KERN_YP_12(1);                            \
   N_KERN_NEG(t03,t06,t04,t07,t05,t08,1,1);    \
   P_KERN_YP_12(2);                            \
   N_KERN_NEG(t03,t06,t04,t07,t05,t08,2,1);    \

#endif
#define N_STORE_YM \
  _mm_store_pd(chi    ,  _mm_add_pd(*(__m128d*)(chi+0),wym[0]));	\
  _mm_store_pd(chi+2  ,  _mm_add_pd(*(__m128d*)(chi+2),wym[1]));	\
  _mm_store_pd(chi+4  ,  _mm_add_pd(*(__m128d*)(chi+4),wym[2]));	\
  _mm_store_pd(chi+6  ,  _mm_add_pd(*(__m128d*)(chi+6),wym[3]));	\
  _mm_store_pd(chi+8  ,  _mm_add_pd(*(__m128d*)(chi+8),wym[4]));	\
  _mm_store_pd(chi+10 ,  _mm_add_pd(*(__m128d*)(chi+10),wym[5]));	\
  									\
  _mm_store_pd(chi+12, _mm_sub_pd(*(__m128d*)(chi+12),	wym[3]));	\
  _mm_store_pd(chi+14, _mm_sub_pd(*(__m128d*)(chi+14),	wym[4]));	\
  _mm_store_pd(chi+16, _mm_sub_pd(*(__m128d*)(chi+16),	wym[5]));	\
  _mm_store_pd(chi+18, _mm_add_pd(*(__m128d*)(chi+18),	wym[0]));	\
  _mm_store_pd(chi+20, _mm_add_pd(*(__m128d*)(chi+20),	wym[1]));	\
  _mm_store_pd(chi+22, _mm_add_pd(*(__m128d*)(chi+22),	wym[2]));	\

#if 1
#define N_KERN_ZM                              \
   P_KERN_ZP_02(0);                            \
   N_HERN(wzm[0],wzm[1],wzm[2],0,2);	       \
   P_KERN_ZP_02(1);                            \
   N_HERN(wzm[0],wzm[1],wzm[2],1,2);	       \
   P_KERN_ZP_02(2);                            \
   N_HERN(wzm[0],wzm[1],wzm[2],2,2);	       \
   					       \
   P_KERN_ZP_13(0);				    \
   N_HERN(wzm[3],wzm[4],wzm[5],0,2);		    \
   P_KERN_ZP_13(1);				    \
   N_HERN(wzm[3],wzm[4],wzm[5],1,2);		    \
   P_KERN_ZP_13(2);				    \
   N_HERN(wzm[3],wzm[4],wzm[5],2,2);		    \

#else
#define N_KERN_ZM                              \
   P_KERN_ZP_02(0);                            \
   N_KERN_IEQ(t00,t06,t01,t07,t02,t08,0,2);        \
   P_KERN_ZP_02(1);                            \
   N_KERN_IEQ(t00,t06,t01,t07,t02,t08,1,2);        \
   P_KERN_ZP_02(2);                            \
   N_KERN_IEQ(t00,t06,t01,t07,t02,t08,2,2);        \
   P_KERN_ZP_13(0);                            \
   N_KERN_INEG(t03,t09,t04,t10,t05,t11,0,2);        \
   P_KERN_ZP_13(1);                            \
   N_KERN_INEG(t03,t09,t04,t10,t05,t11,1,2);        \
   P_KERN_ZP_13(2);                            \
   N_KERN_INEG(t03,t09,t04,t10,t05,t11,2,2);        \

#endif
#define N_STORE_ZM \
  _mm_store_pd(chi    ,  _mm_add_pd(*(__m128d*)(chi+0),wzm[0]));	\
  _mm_store_pd(chi+2  ,  _mm_add_pd(*(__m128d*)(chi+2),wzm[1]));	\
  _mm_store_pd(chi+4  ,  _mm_add_pd(*(__m128d*)(chi+4),wzm[2]));	\
  _mm_store_pd(chi+6  ,  _mm_add_pd(*(__m128d*)(chi+6),wzm[3]));	\
  _mm_store_pd(chi+8  ,  _mm_add_pd(*(__m128d*)(chi+8),wzm[4]));	\
  _mm_store_pd(chi+10 ,  _mm_add_pd(*(__m128d*)(chi+10),wzm[5]));	\
  									\
  _mm_store_pd(chi+12, _mm_addsub_pd(*(__m128d*)(chi+12),_mm_shuffle_pd(wzm[0],wzm[0],1))); \
  _mm_store_pd(chi+14, _mm_addsub_pd(*(__m128d*)(chi+14),_mm_shuffle_pd(wzm[1],wzm[1],1))); \
  _mm_store_pd(chi+16, _mm_addsub_pd(*(__m128d*)(chi+16),_mm_shuffle_pd(wzm[2],wzm[2],1))); \
  _mm_store_pd(chi+18, _mm_add_pd(*(__m128d*)(chi+18),	_mm_xor_pd(_mm_shuffle_pd(wzm[3],wzm[3],1), _mm_set_pd(-0.0,0.0)))); \
  _mm_store_pd(chi+20, _mm_add_pd(*(__m128d*)(chi+20),	_mm_xor_pd(_mm_shuffle_pd(wzm[4],wzm[4],1), _mm_set_pd(-0.0,0.0)))); \
  _mm_store_pd(chi+22, _mm_add_pd(*(__m128d*)(chi+22),	_mm_xor_pd(_mm_shuffle_pd(wzm[5],wzm[5],1), _mm_set_pd(-0.0,0.0)))); \

#if 1
#define N_KERN_TM                              \
   P_KERN_TP_02(0);                            \
   N_HERN(wtm[0],wtm[1],wtm[2],0,3);	       \
   P_KERN_TP_02(1);                            \
   N_HERN(wtm[0],wtm[1],wtm[2],1,3);	       \
   P_KERN_TP_02(2);                            \
   N_HERN(wtm[0],wtm[1],wtm[2],2,3);	       \
   P_KERN_TP_13(0);                            \
   N_HERN(wtm[3],wtm[4],wtm[5],0,3);	       \
   P_KERN_TP_13(1);                            \
   N_HERN(wtm[3],wtm[4],wtm[5],1,3);	       \
   P_KERN_TP_13(2);                            \
   N_HERN(wtm[3],wtm[4],wtm[5],2,3);	       \

#else
#define N_KERN_TM                              \
   P_KERN_TP_02(0);                            \
   N_KERN_NEG(t00,t06,t01,t07,t02,t08,0,3);    \
   P_KERN_TP_02(1);                            \
   N_KERN_NEG(t00,t06,t01,t07,t02,t08,1,3);    \
   P_KERN_TP_02(2);                            \
   N_KERN_NEG(t00,t06,t01,t07,t02,t08,2,3);    \
   P_KERN_TP_13(0);                            \
   N_KERN_NEG(t03,t09,t04,t10,t05,t11,0,3);    \
   P_KERN_TP_13(1);                            \
   N_KERN_NEG(t03,t09,t04,t10,t05,t11,1,3);    \
   P_KERN_TP_13(2);                            \
   N_KERN_NEG(t03,t09,t04,t10,t05,t11,2,3);    \

#endif

#define N_STORE_TM \
  _mm_store_pd(chi    ,  _mm_add_pd(*(__m128d*)(chi+0),wtm[0]));	\
  _mm_store_pd(chi+2  ,  _mm_add_pd(*(__m128d*)(chi+2),wtm[1]));	\
  _mm_store_pd(chi+4  ,  _mm_add_pd(*(__m128d*)(chi+4),wtm[2]));	\
  _mm_store_pd(chi+6  ,  _mm_add_pd(*(__m128d*)(chi+6),wtm[3]));	\
  _mm_store_pd(chi+8  ,  _mm_add_pd(*(__m128d*)(chi+8),wtm[4]));	\
  _mm_store_pd(chi+10 ,  _mm_add_pd(*(__m128d*)(chi+10),wtm[5]));	\
									\
  _mm_store_pd(chi+12, _mm_sub_pd(*(__m128d*)(chi+12),	wtm[0]));	\
  _mm_store_pd(chi+14, _mm_sub_pd(*(__m128d*)(chi+14),	wtm[1]));	\
  _mm_store_pd(chi+16, _mm_sub_pd(*(__m128d*)(chi+16),	wtm[2]));	\
  _mm_store_pd(chi+18, _mm_sub_pd(*(__m128d*)(chi+18),	wtm[3]));	\
  _mm_store_pd(chi+20, _mm_sub_pd(*(__m128d*)(chi+20),	wtm[4]));	\
  _mm_store_pd(chi+22, _mm_sub_pd(*(__m128d*)(chi+22),	wtm[5]));	\

