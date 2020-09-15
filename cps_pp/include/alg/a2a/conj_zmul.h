//Complex multiply with complex conjugate of first and/or second argument (or neither)

#ifndef _CONJ_ZMUL
#define _CONJ_ZMUL

CPS_START_NAMESPACE

//For std::complex types
template<typename mf_Complex, bool conj_left, bool conj_right>
struct Mconj{};

// (re*re - im*im, re*im + im*re )
template<typename mf_Complex>
struct Mconj<mf_Complex,false,false>{
  static inline std::complex<double> doit(const mf_Complex *const l, const mf_Complex *const r){
    return std::complex<double>(l->real()*r->real() - l->imag()*r->imag(), l->real()*r->imag() + l->imag()*r->real());
  }
};
template<typename mf_Complex>
struct Mconj<mf_Complex,false,true>{
  static inline std::complex<double> doit(const mf_Complex *const l, const mf_Complex *const r){
    return std::complex<double>(l->real()*r->real() + l->imag()*r->imag(), l->imag()*r->real() - l->real()*r->imag());
  }
};
template<typename mf_Complex>
struct Mconj<mf_Complex,true,false>{
  static inline std::complex<double> doit(const mf_Complex *const l, const mf_Complex *const r){
    return std::complex<double>(l->real()*r->real() + l->imag()*r->imag(), l->real()*r->imag() - l->imag()*r->real());
  }
};
template<typename mf_Complex>
struct Mconj<mf_Complex,true,true>{
  static inline std::complex<double> doit(const mf_Complex *const l, const mf_Complex *const r){
    return std::complex<double>(l->real()*r->real() - l->imag()*r->imag(), -l->real()*r->imag() - l->imag()*r->real());
  }
};


//For Grid SIMD complex types
#ifdef USE_GRID

template<typename vComplexType, bool conj_left, bool conj_right>
struct MconjGrid{};

template<typename vComplexType>
struct MconjGrid<vComplexType,false,false>{
  static inline vComplexType doit(const vComplexType *const l, const vComplexType *const r){
    return (*l) * (*r);
  }
};
template<typename vComplexType>
struct MconjGrid<vComplexType,false,true>{
  static inline vComplexType doit(const vComplexType *const l, const vComplexType *const r){
    return (*l) * conjugate(*r);
  }
};
template<typename vComplexType>
struct MconjGrid<vComplexType,true,true>{
  static inline vComplexType doit(const vComplexType *const l, const vComplexType *const r){
    return conjugate(*l)*conjugate(*r);
  }
};

//Most common version heavily optimized
#if defined (AVX2)
inline Grid::vComplexD conjMult(const Grid::vComplexD *const l, const Grid::vComplexD *const r){
  Grid::vComplexD out;
  __m256d a_real = _mm256_movedup_pd( l->v ); // Ar Ar
  __m256d a_imag = _mm256_shuffle_pd(l->v,l->v,0xF);//aiai
  a_imag = _mm256_mul_pd( a_imag, _mm256_permute_pd( r->v, 0x5 ) );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
  out.v = _mm256_fmsubadd_pd( a_real, r->v, a_imag ); // Ar Br , Ar Bi   +- Ai Bi             = ArBr+AiBi , ArBi-AiBr
  return out;
}
inline Grid::vComplexF conjMult(const Grid::vComplexF *const l, const Grid::vComplexF *const r){
  Grid::vComplexF out;
  __m256 a_real = _mm256_moveldup_ps( l->v ); // Ar Ar
  __m256 a_imag = _mm256_movehdup_ps( l->v ); // Ai Ai
  a_imag = _mm256_mul_ps( a_imag, _mm256_shuffle_ps( r->v,r->v, _MM_SELECT_FOUR_FOUR(2,3,0,1) ));  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
  out.v = _mm256_fmsubadd_ps( a_real, r->v, a_imag ); // Ar Br , Ar Bi   +- Ai Bi             = ArBr+AiBi , ArBi-AiBr
  return out;
}
#endif

#if defined (AVX512)
inline Grid::vComplexD conjMult(const Grid::vComplexD *const l, const Grid::vComplexD *const r){
  Grid::vComplexD out;
  __m512d a_real = _mm512_shuffle_pd( l->v, l->v, 0x00 );
  __m512d a_imag = _mm512_shuffle_pd( l->v, l->v, 0xFF );
  a_imag = _mm512_mul_pd( a_imag, _mm512_permute_pd( r->v, 0x55 ) );
  out.v = _mm512_fmsubadd_pd( a_real, r->v, a_imag );
  return out;
}

inline Grid::vComplexF conjMult(const Grid::vComplexF *const l, const Grid::vComplexF *const r){
  Grid::vComplexF out;
  __m512 a_real = _mm512_moveldup_ps( l->v ); // Ar Ar
  __m512 a_imag = _mm512_movehdup_ps( l->v ); // Ai Ai
  a_imag = _mm512_mul_ps( a_imag, _mm512_permute_ps( r->v, 0xB1 ) );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
  out.v = _mm512_fmsubadd_ps( a_real, r->v, a_imag ); // Ar Br , Ar Bi   +- Ai Bi             = ArBr+AiBi , ArBi-AiBr 
  return out;
}
#endif


template<typename vComplexType>
struct MconjGrid<vComplexType,true,false>{
  static inline vComplexType doit(const vComplexType *const l, const vComplexType *const r){
#if defined (AVX2) || defined (AVX512)
    return conjMult(l,r);
#else    
    return conjugate(*l) * (*r);
#endif
  }
};


#endif




CPS_END_NAMESPACE


#endif
