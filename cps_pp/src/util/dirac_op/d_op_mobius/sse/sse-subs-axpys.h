




#define DECLARE_AXPY_TMPVARS			\
  register __m128d x0, x1, x2, x3, x4, x5;	\
  register __m128d y0, y1, y2, y3, y4, y5;	\



#define SSE_DAXPY( _A, _X, _Y, _x, _y )		\
    _x = _mm_load_pd( _X );			\
    _y = _mm_load_pd( _Y );			\
    _x = _mm_add_pd( _y, _mm_mul_pd(_A, _x) );	\
    _mm_store_pd( _Y, _x );			\




#define DIST_XY 6
#define DAXPY12( _A, _X, _Y )			\
  SSE_DAXPY(_A, _X+0,  _Y+0, x0, y0);		\
    SSE_DAXPY(_A, _X+2,  _Y+2, x1, y1);		\
    SSE_DAXPY(_A, _X+4,  _Y+4, x2, y2);		\
    SSE_DAXPY(_A, _X+6,  _Y+6, x3, y3);		\
    SSE_DAXPY(_A, _X+8,  _Y+8, x4, y4);		\
    SSE_DAXPY(_A, _X+10, _Y+10, x5, y5);		\
    \
    _mm_prefetch((char *)(_X + DIST_XY*24), _MM_HINT_T0);	\
    _mm_prefetch((char *)(_Y + DIST_XY*24), _MM_HINT_T0);	\





#if 0
#define DECLARE_AXPY_TMPVARS			\
  register __m128d x0, x1, x2, x3, x4, x5;	\
  register __m128d y0, y1, y2, y3, y4, y5;	\
  register __m128d z0, z1;			\



#define DAXPYW( _A, _X, _Y, _W, _x, _y, _z )	\
    _x = _mm_load_pd( _X );			\
    _y = _mm_load_pd( _Y );			\
    _z= _mm_mul_pd(_A, _x);			\
    _z = _mm_add_pd( _y, _z);			\
    _mm_store_pd(_Y, _z);			\

#define DIST_XYW 6
#define DAXPYW12( _A, _X, _Y)				\
    DAXPYW(_A, _X+0,  _Y+0, w0, x0, y0, z0);		\
    DAXPYW(_A, _X+2,  _Y+2, w1, x1, y1, z1);		\
    DAXPYW(_A, _X+4,  _Y+4, w2, x2, y2, z0);		\
    DAXPYW(_A, _X+6,  _Y+6, w2, x3, y3, z1);		\
    DAXPYW(_A, _X+8,  _Y+8, w4, x4, y4, z0);		\
    DAXPYW(_A, _X+10, _Y+10, w5, x5, y5,z1);		\
									\
    _mm_prefetch((char *)(_X + DIST_XYW*24), _MM_HINT_T0);		\
    _mm_prefetch((char *)(_Y + DIST_XYW*24), _MM_HINT_T0);		\
    _mm_prefetch((char *)(_X + DIST_XYW*24+8), _MM_HINT_T0);		\
    _mm_prefetch((char *)(_Y + DIST_XYW*24+8), _MM_HINT_T0);		\


#define DO24(_al, _x, _y)				\
    DAXPYW12(al, f_in  + _shift,			\
	     f_out + _shift + first12_target);		\
							\
    DAXPYW12(al, f_in  + _shift + 12,			\
	     f_out + _shift + second12_target);		\


#endif

