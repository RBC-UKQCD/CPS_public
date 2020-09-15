#ifndef _SPIN_COLOR_MATRICES_H
#define _SPIN_COLOR_MATRICES_H

//These are alternatives to Matrix, WilsonMatrix, SpinColorFlavorMatrix
#include<util/flavormatrix.h>

CPS_START_NAMESPACE

template<typename T>
void CPSsetZero(T &what){
  what = 0.;
}
#ifdef USE_GRID
template<>
void CPSsetZero(Grid::vComplexD &what){
  zeroit(what);
}
template<>
void CPSsetZero(Grid::vComplexF &what){
  zeroit(what);
}
#endif

template<typename T>
void CPSprintT(std::ostream &into, const T &what){
  into << what;
}
template<>
void CPSprintT(std::ostream &into, const std::complex<double> &what){
  into << "(" << what.real() << "," << what.imag() << ")";
}
template<>
void CPSprintT(std::ostream &into, const std::complex<float> &what){
  into << "(" << what.real() << "," << what.imag() << ")";
}

template<typename T>
class isCPSsquareMatrix{

  template<typename U, U> struct Check;

  template<typename U>
  static char test(Check<int, U::isDerivedFromCPSsquareMatrix> *);

  template<typename U>
  static double test(...);
  
public:

  enum {value = sizeof(test<T>(0)) == sizeof(char) };
};


struct cps_square_matrix_mark;

template<typename T>
struct _MatrixClassify{
  typedef typename TestElem< isCPSsquareMatrix<T>::value, cps_square_matrix_mark,LastElem >::type type;
};

template<typename T, typename Tclass>
struct _timespmI{};


template<typename T>
struct _timespmI<T, no_mark>{
  static inline void timesMinusOne(T &out, const T &in){
    out = -in;
  }  
  static inline void timesI(T &out, const T &in){
    out = cps::timesI(in);
  }
  static inline void timesMinusI(T &out, const T &in){
    out = cps::timesMinusI(in);
  }
};
template<typename T>
struct _timespmI<T,cps_square_matrix_mark>{
  static inline void timesMinusOne(T &out, const T &in){
    for(int i=0;i<T::Size;i++)
      for(int j=0;j<T::Size;j++)
	_timespmI<typename T::value_type, typename _MatrixClassify<typename T::value_type>::type>::timesMinusOne(out(i,j), in(i,j));
  }    
  static inline void timesI(T &out, const T &in){    
    for(int i=0;i<T::Size;i++)
      for(int j=0;j<T::Size;j++)
	_timespmI<typename T::value_type, typename _MatrixClassify<typename T::value_type>::type>::timesI(out(i,j), in(i,j));
  }
  static inline void timesMinusI(T &out, const T &in){    
    for(int i=0;i<T::Size;i++)
      for(int j=0;j<T::Size;j++)
	_timespmI<typename T::value_type, typename _MatrixClassify<typename T::value_type>::type>::timesMinusI(out(i,j), in(i,j));
  }
};

template<typename T, typename TypeClass>
struct _CPSsetZeroOne{};

template<typename T>
struct _CPSsetZeroOne<T,no_mark>{
  inline static void setone(T &what){
    what = 1.0;
  }
  inline static void setzero(T &what){
    CPSsetZero(what);
  }
};
template<typename T>
struct _CPSsetZeroOne<T, cps_square_matrix_mark>{
  inline static void setone(T &what){
    what.unit();
  }
  inline static void setzero(T &what){
    what.zero();
  }
};

template<typename T, typename TypeClass>
struct _RecursiveTraceFindScalarType{};

template<typename T>
struct _RecursiveTraceFindScalarType<T, cps_square_matrix_mark>{
  typedef typename _RecursiveTraceFindScalarType<typename T::value_type,typename _MatrixClassify<typename T::value_type>::type>::scalar_type scalar_type;
};
template<typename T>
struct _RecursiveTraceFindScalarType<T, no_mark>{
  typedef T scalar_type;
};
 
template<typename scalar_type, typename T, typename TypeClass>
struct _RecursiveTraceImpl{};

template<typename scalar_type, typename T>
struct _RecursiveTraceImpl<scalar_type, T, cps_square_matrix_mark>{
  static inline void doit(scalar_type &into, const T &what){
    for(int i=0;i<T::Size;i++)
      _RecursiveTraceImpl<scalar_type, typename T::value_type, typename _MatrixClassify<typename T::value_type>::type>::doit(into,what(i,i));    
  }
  static inline void trace_prod(scalar_type &into, const T &a, const T&b){
    for(int i=0;i<T::Size;i++)
      for(int j=0;j<T::Size;j++)
	_RecursiveTraceImpl<scalar_type, typename T::value_type, typename _MatrixClassify<typename T::value_type>::type>::trace_prod(into, a(i,j), b(j,i));
  }
  
};
template<typename scalar_type>
struct _RecursiveTraceImpl<scalar_type, scalar_type, no_mark>{
  static inline void doit(scalar_type &into, const scalar_type &what){
    into = into + what;
  }
  static inline void trace_prod(scalar_type &into, const scalar_type &a, const scalar_type &b){
    into = into + a*b;
  }
};

template<typename T, int RemoveDepth>
struct _PartialTraceFindReducedType{
  typedef typename T::template Rebase< typename _PartialTraceFindReducedType<typename T::value_type,RemoveDepth-1>::type >::type type;
};
template<typename T>
struct _PartialTraceFindReducedType<T,0>{
  typedef typename T::value_type type;
};

template<typename U,typename T, int RemoveDepth>
struct _PartialTraceImpl{
  static inline void doit(U &into, const T&from){
    for(int i=0;i<T::Size;i++)
      for(int j=0;j<T::Size;j++)
	_PartialTraceImpl<typename U::value_type, typename T::value_type, RemoveDepth-1>::doit(into(i,j), from(i,j));
  }
};
template<typename U, typename T>
struct _PartialTraceImpl<U,T,0>{
  static inline void doit(U &into, const T&from){
    for(int i=0;i<T::Size;i++)
      into += from(i,i);
  }
};



template<typename T, int RemoveDepth1, int RemoveDepth2> 
struct _PartialDoubleTraceFindReducedType{
  typedef typename my_enable_if<RemoveDepth1 < RemoveDepth2,int>::type test;
  typedef typename T::template Rebase< typename _PartialDoubleTraceFindReducedType<typename T::value_type,RemoveDepth1-1,RemoveDepth2-1>::type >::type type;
};
template<typename T,int RemoveDepth2>
struct _PartialDoubleTraceFindReducedType<T,0,RemoveDepth2>{
  typedef typename _PartialTraceFindReducedType<typename T::value_type,RemoveDepth2-1>::type type;
};



template<typename U,typename T, int RemoveDepth1,int RemoveDepth2>
struct _PartialDoubleTraceImpl{
  static inline void doit(U &into, const T&from){
    for(int i=0;i<T::Size;i++)
      for(int j=0;j<T::Size;j++)
	_PartialDoubleTraceImpl<typename U::value_type, typename T::value_type, RemoveDepth1-1,RemoveDepth2-1>::doit(into(i,j), from(i,j));
  }
};
template<typename U, typename T, int RemoveDepth2>
struct _PartialDoubleTraceImpl<U,T,0,RemoveDepth2>{
  static inline void doit(U &into, const T&from){
    for(int i=0;i<T::Size;i++)
      _PartialTraceImpl<U,typename T::value_type, RemoveDepth2-1>::doit(into, from(i,i));
  }
};

template<typename T, int TransposeDepth>
struct _IndexTransposeImpl{
  static inline void doit(T &into, const T&from){
    for(int i=0;i<T::Size;i++)
      for(int j=0;j<T::Size;j++)
	_IndexTransposeImpl<typename T::value_type, TransposeDepth-1>::doit(into(i,j), from(i,j));
  }
};
template<typename T>
struct _IndexTransposeImpl<T,0>{
  static inline void doit(T &into, const T&from){
    for(int i=0;i<T::Size;i++)
      for(int j=0;j<T::Size;j++)
	into(i,j) = from(j,i);
  }
};




template<typename T, int N>
class CPSsquareMatrix{
protected:
  T v[N][N];
  template<typename U>  //, typename my_enable_if< my_is_base_of<U, CPSsquareMatrix<T,N> >::value, int>::type = 0
  static inline void mult(U &into, const U &a, const U &b){
    into.zero();
    for(int i=0;i<N;i++)
      for(int k=0;k<N;k++)
	for(int j=0;j<N;j++)
	  into(i,k) = into(i,k) + a(i,j)*b(j,k);      
  }
public:
  enum { isDerivedFromCPSsquareMatrix };
  enum { Size = N };
  typedef T value_type;
  typedef typename _RecursiveTraceFindScalarType<CPSsquareMatrix<T,N>,cps_square_matrix_mark>::scalar_type scalar_type;
  
  template<typename U>
  struct Rebase{
    typedef CPSsquareMatrix<U,N> type;
  };
  
  T & operator()(const int i, const int j){ return v[i][j]; }
  const T & operator()(const int i, const int j) const{ return v[i][j]; }

  inline void equalsTranspose(const CPSsquareMatrix<T,N> &r){
    for(int i=0;i<N;i++)
      for(int j=0;j<N;j++)
	v[i][j] = r.v[j][i];
  }
  inline CPSsquareMatrix<T,N> Transpose() const{
    CPSsquareMatrix<T,N> out; out.equalsTranspose(*this);
    return out;
  }				
  
  //Trace on just the index associated with this matrix
  // T traceIndex() const{
  //   T out; CPSsetZero(out);
  //   for(int i=0;i<N;i++) out = out + v[i][i];
  //   return out;
  // }
  //Trace on all indices recursively
  scalar_type Trace() const{
    scalar_type ret; CPSsetZero(ret);
    _RecursiveTraceImpl<scalar_type, CPSsquareMatrix<T,N>, cps_square_matrix_mark>::doit(ret, *this);
    return ret;
  }

  template<int RemoveDepth>
  typename _PartialTraceFindReducedType<CPSsquareMatrix<T,N>, RemoveDepth>::type TraceIndex() const{
    typedef typename _PartialTraceFindReducedType<CPSsquareMatrix<T,N>, RemoveDepth>::type ReducedType;
    ReducedType into; _CPSsetZeroOne<ReducedType, typename _MatrixClassify<ReducedType>::type>::setzero(into);//  into.zero();
    _PartialTraceImpl<ReducedType, CPSsquareMatrix<T,N>, RemoveDepth>::doit(into, *this);
    return into;
  }

  template<int RemoveDepth1, int RemoveDepth2>
  typename _PartialDoubleTraceFindReducedType<CPSsquareMatrix<T,N>,RemoveDepth1,RemoveDepth2>::type TraceTwoIndices() const{
    typedef typename _PartialDoubleTraceFindReducedType<CPSsquareMatrix<T,N>,RemoveDepth1,RemoveDepth2>::type ReducedType;
    ReducedType into; _CPSsetZeroOne<ReducedType, typename _MatrixClassify<ReducedType>::type>::setzero(into);//  into.zero();
    _PartialDoubleTraceImpl<ReducedType, CPSsquareMatrix<T,N>, RemoveDepth1,RemoveDepth2>::doit(into, *this);
    return into;
  }

  //Set this matrix equal to the transpose of r on its compound tensor index TransposeDepth. If it is not a compound matrix then TranposeDepth=0 does the same thing as equalsTranspose above
  template<int TransposeDepth>
  inline void equalsTransposeOnIndex(const CPSsquareMatrix<T,N> &r){
    assert(&r != this);
    _IndexTransposeImpl<CPSsquareMatrix<T,N>, TransposeDepth>::doit(*this,r);
  }
  template<int TransposeDepth>
  inline CPSsquareMatrix<T,N> & TransposeOnIndex(){
    CPSsquareMatrix<T,N> cp(*this);
    this->equalsTransposeOnIndex<TransposeDepth>(cp);
    return *this;
  }

  CPSsquareMatrix<T,N> & zero(){
    for(int i=0;i<N;i++)
      for(int j=0;j<N;j++)
	_CPSsetZeroOne<T,  typename _MatrixClassify<T>::type>::setzero(v[i][j]);
    return *this;    
  }
  CPSsquareMatrix<T,N>& unit(){
    zero();
    for(int i=0;i<N;i++)
      _CPSsetZeroOne<T,  typename _MatrixClassify<T>::type>::setone(v[i][i]);
    return *this;
  }
  //this = -this
  CPSsquareMatrix<T,N> & timesMinusOne(){
    _timespmI<CPSsquareMatrix<T,N>,cps_square_matrix_mark>::timesMinusOne(*this,*this);
    return *this;
  }
    
  CPSsquareMatrix<T,N> & operator+=(const CPSsquareMatrix<T,N> &r){
    for(int i=0;i<N;i++)
      for(int j=0;j<N;j++)
	v[i][j] = v[i][j] + r.v[i][j];
    return *this;
  }
  CPSsquareMatrix<T,N> & operator*=(const T &r){
    for(int i=0;i<N;i++)
      for(int j=0;j<N;j++)
	v[i][j] = v[i][j] * r;
    return *this;
  }
  bool operator==(const CPSsquareMatrix<T,N> &r) const{
    for(int i=0;i<N;i++)
      for(int j=0;j<N;j++)
	if(v[i][j] != r.v[i][j]) return false;
    return true;
  }
  friend std::ostream & operator<<(std::ostream &os, const CPSsquareMatrix<T,N> &m){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	CPSprintT(os, m.v[i][j]);
	if(j<N-1) os << " ";
      }
      os << '\n';
    }
    return os;
  }
};

template<typename U> 
typename my_enable_if<isCPSsquareMatrix<U>::value, U>::type operator*(const U &a, const U &b){
  U into;
  into.zero();
  for(int i=0;i<U::Size;i++)
    for(int k=0;k<U::Size;k++)
      for(int j=0;j<U::Size;j++)
	into(i,k) = into(i,k) + a(i,j)*b(j,k);
  return into;
}
template<typename U> 
typename my_enable_if<isCPSsquareMatrix<U>::value, U>::type operator+(const U &a, const U &b){
  U into;
  for(int i=0;i<U::Size;i++)
    for(int j=0;j<U::Size;j++)
      into(i,j) = a(i,j) + b(i,j);
  return into;
}

template<typename U>
typename my_enable_if<isCPSsquareMatrix<U>::value, typename U::scalar_type>::type Trace(const U &a, const U &b){
  typename U::scalar_type out;
  CPSsetZero(out);
  _RecursiveTraceImpl<typename U::scalar_type, U, cps_square_matrix_mark>::trace_prod(out, a,b);
  return out;
}

template<typename T>
T Transpose(const T& r){
  T out;
  out.equalsTranspose(r);
  return out;
}



//Old annoying CPS conventions with output on RHS (Fortran user or something??)
#define TIMESPLUSONE(a,b) { b=a; }
#define TIMESMINUSONE(a,b) { _timespmI<T, typename _MatrixClassify<T>::type>::timesMinusOne(b,a); }
#define TIMESPLUSI(a,b) { _timespmI<T, typename _MatrixClassify<T>::type>::timesI(b,a); }
#define TIMESMINUSI(a,b) { _timespmI<T, typename _MatrixClassify<T>::type>::timesMinusI(b,a); }
#define SETZERO(a){ _CPSsetZeroOne<T, typename _MatrixClassify<T>::type>::setzero(a); }

template<typename T>
class CPSflavorMatrix: public CPSsquareMatrix<T,2>{
public:
  typedef typename CPSsquareMatrix<T,2>::value_type value_type;
  typedef typename CPSsquareMatrix<T,2>::scalar_type scalar_type;

  inline operator CPSsquareMatrix<T,2>(){ return static_cast<CPSsquareMatrix<T,2> &>(*this); }

  CPSflavorMatrix(const CPSsquareMatrix<T,2> &r): CPSsquareMatrix<T,2>(r){}
  CPSflavorMatrix(): CPSsquareMatrix<T,2>(){}  
  
  template<typename U>
  struct Rebase{
    typedef CPSflavorMatrix<U> type;
  };
  
  //multiply on left by a flavor matrix
  CPSflavorMatrix<T> & pl(const FlavorMatrixType &type){
    T tmp1, tmp2;
    T (&v)[2][2] = this->v;
    
    switch( type ){
    case F0:
      SETZERO(v[1][0]);
      SETZERO(v[1][1]);
      break;
    case F1:
      SETZERO(v[0][0]);
      SETZERO(v[0][1]);
      break;
    case Fud:
      tmp1 = v[0][0];
      tmp2 = v[0][1];
      v[0][0] = v[1][0];
      v[0][1] = v[1][1];
      v[1][0] = tmp1;
      v[1][1] = tmp2;
      break;
    case sigma0:
      break;
    case sigma1:
      tmp1 = v[0][0];
      tmp2 = v[0][1];
      v[0][0] = v[1][0];
      v[0][1] = v[1][1];
      v[1][0] = tmp1;
      v[1][1] = tmp2;
      break;      
    case sigma2:
      TIMESPLUSI(v[0][0], tmp1);
      TIMESPLUSI(v[0][1], tmp2);
      TIMESMINUSI(v[1][0], v[0][0]);
      TIMESMINUSI(v[1][1], v[0][1]);
      v[1][0] = tmp1;
      v[1][1] = tmp2;
      break;
    case sigma3:
      TIMESMINUSONE(v[1][0],v[1][0]);
      TIMESMINUSONE(v[1][1],v[1][1]);
      break;
    default:
      ERR.General("FlavorMatrixGeneral","pl(const FlavorMatrixGeneralType &type)","Unknown FlavorMatrixGeneralType");
      break;
    }
    return *this;

  }

  //multiply on right by a flavor matrix
  CPSflavorMatrix<T> & pr(const FlavorMatrixType &type){
    T tmp1, tmp2;
    T (&v)[2][2] = this->v;
    
    switch(type){
    case F0:     
      SETZERO(v[0][1]);
      SETZERO(v[1][1]);
      break;
    case F1:
      SETZERO(v[0][0]);
      SETZERO(v[1][0]);
      break;
    case Fud:
      tmp1 = v[0][0];
      tmp2 = v[1][0];
      v[0][0] = v[0][1];
      v[1][0] = v[1][1];
      v[0][1] = tmp1;
      v[1][1] = tmp2;
      break;
    case sigma0:
      break;
    case sigma1:
      tmp1 = v[0][0];
      tmp2 = v[1][0];
      v[0][0] = v[0][1];
      v[1][0] = v[1][1];
      v[0][1] = tmp1;
      v[1][1] = tmp2;
      break;      
    case sigma2:
      TIMESMINUSI(v[0][0], tmp1);
      TIMESMINUSI(v[1][0], tmp2);
      TIMESPLUSI(v[0][1], v[0][0]);
      TIMESPLUSI(v[1][1], v[1][0]);
      v[0][1] = tmp1;
      v[1][1] = tmp2;
      break;
    case sigma3:
      TIMESMINUSONE(v[0][1],v[0][1]);
      TIMESMINUSONE(v[1][1],v[1][1]);
      break;
    default:
      ERR.General("FlavorMatrixGeneral","pr(const FlavorMatrixGeneralType &type)","Unknown FlavorMatrixGeneralType");
      break;
    }
    return *this;
  }

};


  
  
//  Chiral basis
//  gamma(XUP)    gamma(YUP)    gamma(ZUP)    gamma(TUP)    gamma(FIVE)
//  0  0  0  i    0  0  0 -1    0  0  i  0    0  0  1  0    1  0  0  0
//  0  0  i  0    0  0  1  0    0  0  0 -i    0  0  0  1    0  1  0  0
//  0 -i  0  0    0  1  0  0   -i  0  0  0    1  0  0  0    0  0 -1  0
// -i  0  0  0   -1  0  0  0    0  i  0  0    0  1  0  0    0  0  0 -1

template<typename T>
class CPSspinMatrix: public CPSsquareMatrix<T,4>{
public:
  typedef typename CPSsquareMatrix<T,4>::value_type value_type;
  typedef typename CPSsquareMatrix<T,4>::scalar_type scalar_type;

  inline operator CPSsquareMatrix<T,4>(){ return static_cast<CPSsquareMatrix<T,4> &>(*this); }

  CPSspinMatrix(const CPSsquareMatrix<T,4> &r): CPSsquareMatrix<T,4>(r){}
  CPSspinMatrix(): CPSsquareMatrix<T,4>(){}  

  
  template<typename U>
  struct Rebase{
    typedef CPSspinMatrix<U> type;
  };
  
  //Left Multiplication by Dirac gamma's
  CPSspinMatrix<T> & gl(int dir){
    int s2;
    CPSspinMatrix<T> cp(*this);
    const T (&src)[4][4] = cp.v;
    T (&p)[4][4] = this->v;
    
    switch(dir){
    case 0:
      for(s2=0;s2<4;s2++){
	TIMESPLUSI(  src[3][s2], p[0][s2] );
	TIMESPLUSI(  src[2][s2], p[1][s2] );
	TIMESMINUSI( src[1][s2], p[2][s2] );
	TIMESMINUSI( src[0][s2], p[3][s2] );
      }
      break;
    case 1:
      for(s2=0;s2<4;s2++){
	TIMESMINUSONE( src[3][s2], p[0][s2] );
	TIMESPLUSONE(  src[2][s2], p[1][s2] );
	TIMESPLUSONE(  src[1][s2], p[2][s2] );
	TIMESMINUSONE( src[0][s2], p[3][s2] );
      }
      break;
    case 2:
      for(s2=0;s2<4;s2++){
	TIMESPLUSI(  src[2][s2], p[0][s2] );
	TIMESMINUSI( src[3][s2], p[1][s2] );
	TIMESMINUSI( src[0][s2], p[2][s2] );
	TIMESPLUSI(  src[1][s2], p[3][s2] );
      }
      break;
    case 3:
      for(s2=0;s2<4;s2++){
	TIMESPLUSONE( src[2][s2], p[0][s2] );
	TIMESPLUSONE( src[3][s2], p[1][s2] );
	TIMESPLUSONE( src[0][s2], p[2][s2] );
	TIMESPLUSONE( src[1][s2], p[3][s2] );
      }
      break;
    case -5:
      for(s2=0;s2<4;s2++){
	TIMESPLUSONE(  src[0][s2], p[0][s2] );
	TIMESPLUSONE(  src[1][s2], p[1][s2] );
	TIMESMINUSONE( src[2][s2], p[2][s2] );
	TIMESMINUSONE( src[3][s2], p[3][s2] );
      }
      break;
    default:
      assert(0);
      break;
    }
    return *this;
  }


  //Right Multiplication by Dirac gamma's

  CPSspinMatrix<T>& gr(int dir)
  {
    int s1;
    CPSspinMatrix<T> cp(*this);
    const T (&src)[4][4] = cp.v;
    T (&p)[4][4] = this->v;

    switch(dir){
    case 0:
      for(s1=0;s1<4;s1++){
	TIMESMINUSI( src[s1][3], p[s1][0] );
	TIMESMINUSI( src[s1][2], p[s1][1] );
	TIMESPLUSI(  src[s1][1], p[s1][2] );
	TIMESPLUSI(  src[s1][0], p[s1][3] );
      }
      break;
    case 1:
      for(s1=0;s1<4;s1++){
	TIMESMINUSONE( src[s1][3], p[s1][0] );
	TIMESPLUSONE(  src[s1][2], p[s1][1] );
	TIMESPLUSONE(  src[s1][1], p[s1][2] );
	TIMESMINUSONE( src[s1][0], p[s1][3] );
      }
      break;
    case 2:
      for(s1=0;s1<4;s1++){
	TIMESMINUSI( src[s1][2], p[s1][0] );
	TIMESPLUSI(  src[s1][3], p[s1][1] );
	TIMESPLUSI(  src[s1][0], p[s1][2] );
	TIMESMINUSI( src[s1][1], p[s1][3] );
      }
      break;
    case 3:
      for(s1=0;s1<4;s1++){
	TIMESPLUSONE( src[s1][2], p[s1][0] );
	TIMESPLUSONE( src[s1][3], p[s1][1] );
	TIMESPLUSONE( src[s1][0], p[s1][2] );
	TIMESPLUSONE( src[s1][1], p[s1][3] );
      }
      break;
    case -5:
      for(s1=0;s1<4;s1++){
	TIMESPLUSONE(  src[s1][0], p[s1][0] );
	TIMESPLUSONE(  src[s1][1], p[s1][1] );
	TIMESMINUSONE( src[s1][2], p[s1][2] );
	TIMESMINUSONE( src[s1][3], p[s1][3] );
      }
      break;
    default:
      assert(0);
      break;
    }
    return *this;
  }

  //multiply gamma(i)gamma(5) on the left: result = gamma(i)*gamma(5)*from
  CPSspinMatrix<T>& glAx(const int dir){
    int s2;
    CPSspinMatrix<T> cp(*this);
    const T (&from_mat)[4][4] = cp.v;
    T (&p)[4][4] = this->v;
    
    switch(dir){
    case 0:
      for(s2=0;s2<4;s2++){
            TIMESMINUSI( from_mat[3][s2], p[0][s2] );
            TIMESMINUSI( from_mat[2][s2], p[1][s2] );
            TIMESMINUSI( from_mat[1][s2], p[2][s2] );
            TIMESMINUSI( from_mat[0][s2], p[3][s2] );
        }
        break;
    case 1:
        for(s2=0;s2<4;s2++){
            TIMESPLUSONE(  from_mat[3][s2], p[0][s2] );
            TIMESMINUSONE( from_mat[2][s2], p[1][s2] );
            TIMESPLUSONE(  from_mat[1][s2], p[2][s2] );
            TIMESMINUSONE( from_mat[0][s2], p[3][s2] );
        }
        break;
    case 2:
        for(s2=0;s2<4;s2++){
            TIMESMINUSI( from_mat[2][s2], p[0][s2] );
            TIMESPLUSI(  from_mat[3][s2], p[1][s2] );
            TIMESMINUSI( from_mat[0][s2], p[2][s2] );
            TIMESPLUSI(  from_mat[1][s2], p[3][s2] );
        }
	break;
    case 3:
        for(s2=0;s2<4;s2++){
            TIMESMINUSONE( from_mat[2][s2], p[0][s2] );
            TIMESMINUSONE( from_mat[3][s2], p[1][s2] );
            TIMESPLUSONE(  from_mat[0][s2], p[2][s2] );
            TIMESPLUSONE(  from_mat[1][s2], p[3][s2] );
        }
        break;
    default:
      assert(0);
      break;
    }
    return *this;
  }

  
};



template<typename T>
class CPScolorMatrix: public CPSsquareMatrix<T,3>{
public:
  typedef typename CPSsquareMatrix<T,3>::value_type value_type;
  typedef typename CPSsquareMatrix<T,3>::scalar_type scalar_type;

  template<typename U>
  struct Rebase{
    typedef CPScolorMatrix<U> type;
  };
  inline operator CPSsquareMatrix<T,3>(){ return static_cast<CPSsquareMatrix<T,3> &>(*this); }

  CPScolorMatrix(const CPSsquareMatrix<T,3> &r): CPSsquareMatrix<T,3>(r){}
  CPScolorMatrix(): CPSsquareMatrix<T,3>(){}  
};

template<typename ComplexType>
class CPSspinColorFlavorMatrix: public CPSspinMatrix<CPScolorMatrix<CPSflavorMatrix<ComplexType> > >{
public:
  typedef CPSspinMatrix<CPScolorMatrix<CPSflavorMatrix<ComplexType> > > SCFmat;
  typedef typename SCFmat::value_type value_type;
  typedef typename SCFmat::scalar_type scalar_type;

  CPSspinColorFlavorMatrix(): SCFmat(){}
  CPSspinColorFlavorMatrix(const SCFmat &r): SCFmat(r){}
  inline operator SCFmat(){ return static_cast<SCFmat&>(*this); }
  
  template<typename U>
  struct Rebase{
    typedef CPSspinMatrix<U> type;
  };

  inline CPSspinColorFlavorMatrix<ComplexType>& unit(){
    return static_cast<CPSspinColorFlavorMatrix<ComplexType>& >(this->CPSsquareMatrix<value_type,4>::unit());
  }
  inline CPSspinColorFlavorMatrix<ComplexType>& zero(){
    return static_cast<CPSspinColorFlavorMatrix<ComplexType>& >(this->CPSsquareMatrix<value_type,4>::zero());
  }
  inline CPSspinColorFlavorMatrix<ComplexType>& timesMinusOne(){
    return static_cast<CPSspinColorFlavorMatrix<ComplexType>& >(this->CPSsquareMatrix<value_type,4>::timesMinusOne());
  }

  inline CPScolorMatrix<ComplexType> SpinFlavorTrace() const{
    return this->CPSsquareMatrix<value_type,4>::template TraceTwoIndices<0,2>();
  }
  inline CPSspinMatrix<CPSflavorMatrix<ComplexType> > ColorTrace() const{
    return this->CPSsquareMatrix<value_type,4>::template TraceIndex<1>();
  }
  inline CPSspinColorFlavorMatrix<ComplexType>& TransposeColor(){
    return static_cast<CPSspinColorFlavorMatrix<ComplexType> &>(this->CPSsquareMatrix<value_type,4>::template TransposeOnIndex<1>());
  }
  inline void equalsColorTranspose(const CPSspinColorFlavorMatrix<ComplexType> &r){
    this->CPSsquareMatrix<value_type,4>::template equalsTransposeOnIndex<1>(r);
  }
  
  //multiply on left by a flavor matrix
  inline SCFmat & pl(const FlavorMatrixType type){
    for(int s1=0;s1<4;s1++)
      for(int s2=0;s2<4;s2++)
	for(int c1=0;c1<3;c1++)
	  for(int c2=0;c2<3;c2++)
	    this->operator()(s1,s2)(c1,c2).pl(type);
    return *this;
  }
  //multiply on left by a flavor matrix
  inline SCFmat & pr(const FlavorMatrixType type){
    for(int s1=0;s1<4;s1++)
      for(int s2=0;s2<4;s2++)
	for(int c1=0;c1<3;c1++)
	  for(int c2=0;c2<3;c2++)
	    this->operator()(s1,s2)(c1,c2).pr(type);
    return *this;
  }

  
};





#undef TIMESPLUSONE
#undef TIMESMINUSONE
#undef TIMESPLUSI
#undef TIMESMINUSI
#undef SETZERO



CPS_END_NAMESPACE

#endif
