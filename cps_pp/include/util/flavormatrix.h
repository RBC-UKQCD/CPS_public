#ifndef _FLAVOR_MATRIX_H
#define _FLAVOR_MATRIX_H
CPS_START_NAMESPACE

template<typename T>
inline T timesI(const T &v){
  return T( -v.imag(), v.real() );
}
template<typename T>
inline T timesMinusI(const T &v){
  return T( v.imag(), -v.real() );
}

#ifdef USE_GRID
CPS_END_NAMESPACE
#include<Grid/Grid.h>
CPS_START_NAMESPACE

template<>
inline Grid::vComplexD timesI(const Grid::vComplexD &v){
  return Grid::timesI(v);
}
template<>
inline Grid::vComplexD timesMinusI(const Grid::vComplexD &v){
  return Grid::timesMinusI(v);
}
template<>
inline Grid::vComplexF timesI(const Grid::vComplexF &v){
  return Grid::timesI(v);
}
template<>
inline Grid::vComplexF timesMinusI(const Grid::vComplexF &v){
  return Grid::timesMinusI(v);
}
#endif

//Note: F0 = 1/2(1+sigma3)  and  F1 = 1/2(1-sigma3). These are called F11 and F22 in the paper, respectively
enum FlavorMatrixType {F0, F1, Fud, sigma0, sigma1, sigma2, sigma3};

//Added by CK, optimized from Daiqian's equivalent
template<typename T = Complex>
class FlavorMatrixGeneral{
protected:
  T fmat[2][2];
public:
  FlavorMatrixGeneral(){
  }
  FlavorMatrixGeneral(const FlavorMatrixGeneral<T> &from){
    for(int i = 0 ; i < 2; i++)
      for(int j = 0 ; j < 2; j++)
	fmat[i][j] = from.fmat[i][j];
  }

  template<typename U>
  FlavorMatrixGeneral(const U rhs){
    for(int i = 0 ; i < 2; i++)
      for(int j = 0 ; j < 2; j++)
	fmat[i][j] = rhs;
  }

  void zero(){
    for(int i = 0 ; i < 2; i++)
      for(int j = 0 ; j < 2; j++)
	fmat[i][j] = 0.0;
  }

  FlavorMatrixGeneral<T>& operator=(const FlavorMatrixGeneral<T> &from){
    for(int i = 0 ; i < 2; i++)
      for(int j = 0 ; j < 2; j++)
	fmat[i][j] = from.fmat[i][j];
    return *this;
  }

  inline T &operator()(int f_row,int f_col){
    return fmat[f_row][f_col];
  }
  inline const T &operator()(int f_row,int f_col) const{
    return fmat[f_row][f_col];
  }

  FlavorMatrixGeneral<T> & operator+=(const FlavorMatrixGeneral<T> &r){
    for(int i = 0 ; i < 2; i++)
      for(int j = 0 ; j < 2; j++)
	fmat[i][j] += r.fmat[i][j];
    return *this;
  }
  template<typename U>
  FlavorMatrixGeneral<T> & operator*=(const U rhs){ 
    for(int i = 0 ; i < 2; i++)
      for(int j = 0 ; j < 2; j++)
	fmat[i][j] *= rhs;
    return *this;
  }

  FlavorMatrixGeneral<T> operator*(const FlavorMatrixGeneral<T>& rhs) const{
    FlavorMatrixGeneral<T> out;
    out.fmat[0][0] = fmat[0][0]*rhs.fmat[0][0] + fmat[0][1]*rhs.fmat[1][0];
    out.fmat[0][1] = fmat[0][0]*rhs.fmat[0][1] + fmat[0][1]*rhs.fmat[1][1];
    out.fmat[1][0] = fmat[1][0]*rhs.fmat[0][0] + fmat[1][1]*rhs.fmat[1][0];
    out.fmat[1][1] = fmat[1][0]*rhs.fmat[0][1] + fmat[1][1]*rhs.fmat[1][1];
    return out;
  }

  FlavorMatrixGeneral<T> transpose() const{
    FlavorMatrixGeneral<T> out;
    out.fmat[0][0] = fmat[0][0];
    out.fmat[0][1] = fmat[1][0];
    out.fmat[1][0] = fmat[0][1];
    out.fmat[1][1] = fmat[1][1];
    return out;
  }

  T Trace() const{
    return fmat[0][0] + fmat[1][1];
  }

  //Added by CK
  //multiply on left by a flavor matrix
  FlavorMatrixGeneral<T> & pl(const FlavorMatrixType &type){
    T tmp1, tmp2;

    switch( type ){
    case F0:
      (*this)(1,0) = 0.0;
      (*this)(1,1) = 0.0;
      break;
    case F1:
      (*this)(0,0) = 0.0;
      (*this)(0,1) = 0.0;
      break;
    case Fud:
      tmp1 = (*this)(0,0);
      tmp2 = (*this)(0,1);
      (*this)(0,0) = (*this)(1,0);
      (*this)(0,1) = (*this)(1,1);
      (*this)(1,0) = tmp1;
      (*this)(1,1) = tmp2;
      break;
    case sigma0:
      break;
    case sigma1:
      tmp1 = (*this)(0,0);
      tmp2 = (*this)(0,1);
      (*this)(0,0) = (*this)(1,0);
      (*this)(0,1) = (*this)(1,1);
      (*this)(1,0) = tmp1;
      (*this)(1,1) = tmp2;
      break;      
    case sigma2:
      tmp1 = timesI( this->operator()(0,0) ); //(*this)(0,0)*T(0.0,1.0);
      tmp2 = timesI( this->operator()(0,1) ); //(*this)(0,1)*T(0.0,1.0);
      (*this)(0,0) = timesMinusI( this->operator()(1,0) );   //(*this)(1,0)*T(0.0,-1.0);
      (*this)(0,1) = timesMinusI( this->operator()(1,1) );   //(*this)(1,1)*T(0.0,-1.0);
      (*this)(1,0) = tmp1;
      (*this)(1,1) = tmp2;
      break;
    case sigma3:
      (*this)(1,0)*=-1.0;
      (*this)(1,1)*=-1.0;
      break;
    default:
      ERR.General("FlavorMatrixGeneral","pl(const FlavorMatrixGeneralType &type)","Unknown FlavorMatrixGeneralType");
      break;
    }
    return *this;

  }

  //multiply on right by a flavor matrix
  FlavorMatrixGeneral<T> & pr(const FlavorMatrixType &type){
    T tmp1, tmp2;

    switch(type){
    case F0:
      (*this)(0,1) = 0.0;
      (*this)(1,1) = 0.0;
      break;
    case F1:
      (*this)(0,0) = 0.0;
      (*this)(1,0) = 0.0;
      break;
    case Fud:
      tmp1 = (*this)(0,0);
      tmp2 = (*this)(1,0);
      (*this)(0,0) = (*this)(0,1);
      (*this)(1,0) = (*this)(1,1);
      (*this)(0,1) = tmp1;
      (*this)(1,1) = tmp2;
      break;
    case sigma0:
      break;
    case sigma1:
      tmp1 = (*this)(0,0);
      tmp2 = (*this)(1,0);
      (*this)(0,0) = (*this)(0,1);
      (*this)(1,0) = (*this)(1,1);
      (*this)(0,1) = tmp1;
      (*this)(1,1) = tmp2;
      break;      
    case sigma2:
      tmp1 = timesMinusI(this->operator()(0,0)); //(*this)(0,0) *  T(0.0,-1.0);
      tmp2 = timesMinusI(this->operator()(1,0)); //(*this)(1,0) *  T(0.0,-1.0);
      (*this)(0,0) = timesI(this->operator()(0,1)); //(*this)(0,1)* T(0.0,1.0); 
      (*this)(1,0) = timesI(this->operator()(1,1)); //(*this)(1,1)* T(0.0,1.0);
      (*this)(0,1) = tmp1;
      (*this)(1,1) = tmp2;
      break;
    case sigma3:
      (*this)(0,1)*=-1.0;
      (*this)(1,1)*=-1.0;
      break;
    default:
      ERR.General("FlavorMatrixGeneral","pr(const FlavorMatrixGeneralType &type)","Unknown FlavorMatrixGeneralType");
      break;
    }
    return *this;
  }
};
//Trace(A * B);
template<typename T>
inline static T Trace(const FlavorMatrixGeneral<T> &a, const FlavorMatrixGeneral<T> &b){
  return a(0,0)*b(0,0) + a(0,1)*b(1,0) + a(1,0)*b(0,1) + a(1,1)*b(1,1);
}
//Trace(A^T * B)
template<typename T>
inline static T TransLeftTrace(const FlavorMatrixGeneral<T> &a, const FlavorMatrixGeneral<T> &b){
  return a(0,0)*b(0,0) + a(1,0)*b(1,0) + a(0,1)*b(0,1) + a(1,1)*b(1,1);
}

typedef FlavorMatrixGeneral<Complex> FlavorMatrix;

CPS_END_NAMESPACE
#endif
