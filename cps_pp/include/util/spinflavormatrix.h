#ifndef SPIN_FLAVOR_MATRIX_H
#define SPIN_FLAVOR_MATRIX_H

#include<config.h>
#include <alg/alg_base.h>
#include <alg/wilson_matrix.h>
CPS_START_NAMESPACE

class SpinFlavorMatrix
{
  static const int fsize = 2*4*4*2*2;
  Float u[fsize];	// The matrix


  public:
    SpinFlavorMatrix();
    SpinFlavorMatrix(Float c);
    SpinFlavorMatrix(const Complex& c);
    SpinFlavorMatrix(const SpinFlavorMatrix& m);

    //Set all elements equal to float or complex
    SpinFlavorMatrix& operator=(Float c);
    SpinFlavorMatrix& operator=(const Complex& c);
    SpinFlavorMatrix& operator=(const SpinFlavorMatrix& c);

    Complex& operator()(const int &s1, const int &f1, const int &s2, const int &f2);
    const Complex& operator()(const int &s1, const int &f1, const int &s2, const int &f2) const;
    Complex& operator[](int i) { return ((Complex*)u)[i]; }
    const Complex& operator[](int i) const { return ((Complex*)u)[i]; }

    SpinFlavorMatrix operator*(const Complex &c)const {
        SpinFlavorMatrix ret;
	for(int i=0;i<4*4*2*2;i++)
	  ret[i] = (*this)[i] * c;	
        return ret;
    }

    const SpinFlavorMatrix &operator*=(const Complex &c) {
	for(int i=0;i<4*4*2*2;i++)
	  (*this)[i] *= c;
        return *this;
    }

    SpinFlavorMatrix operator+(const SpinFlavorMatrix &s)const {
      SpinFlavorMatrix ret;
      for(int i=0;i<4*4*2*2;i++)
	ret[i] = (*this)[i] + s[i];	
      return ret;
    }

    const SpinFlavorMatrix &operator+=(const SpinFlavorMatrix &s) {
      for(int i=0;i<4*4*2*2;i++)
	(*this)[i] += s[i];
      return *this;
    }

    SpinFlavorMatrix operator-(const SpinFlavorMatrix &s)const {
      SpinFlavorMatrix ret;
      for(int i=0;i<4*4*2*2;i++)
	ret[i] = (*this)[i] - s[i];	
      return ret;
    }

    const SpinFlavorMatrix &operator-=(const SpinFlavorMatrix &s) {
      for(int i=0;i<4*4*2*2;i++)
	(*this)[i] -= s[i];
      return *this;
    }


    SpinFlavorMatrix operator*(const SpinFlavorMatrix &r)const {
        SpinFlavorMatrix ret(0.0);
	for(int s1=0;s1<4;s1++)
	  for(int f1=0;f1<2;f1++)
	    for(int s2=0;s2<4;s2++)
	      for(int f2=0;f2<2;f2++)
		for(int s3=0;s3<4;s3++)
		  for(int f3=0;f3<2;f3++)
		    ret(s1,f1,s2,f2) += (*this)(s1,f1,s3,f3) * r(s3,f3,s2,f2);
	return ret;
    }

    Complex Trace() const;
};

static inline Rcomplex Trace(const SpinFlavorMatrix &a, const SpinFlavorMatrix &b) {
    Rcomplex ret = 0;
    for(int s1=0;s1<4;s1++)
      for(int f1=0;f1<2;f1++)
	for(int s2=0;s2<4;s2++)
	  for(int f2=0;f2<2;f2++)
	    ret += a(s1,f1,s2,f2) * b(s2,f2,s1,f1);
    return ret;
}






//CK: Ported and renamed from Daiqian's A2A code 'SpinFlavorMatrix'. I also replaced the heap-allocated SpinMatrix with stack-allocated
class FlavorSpinMatrix{
protected:
  SpinMatrix smat[2][2];
  const char *cname;
public:
  FlavorSpinMatrix(): cname("FlavorSpinMatrix"){
  }

  FlavorSpinMatrix(const Float &val): cname("FlavorSpinMatrix"){
    for(int i=0;i<2;i++){
      smat[0][i] = val;
      smat[1][i] = val;
    }
  }

  FlavorSpinMatrix(const FlavorSpinMatrix &from): cname("FlavorSpinMatrix"){
    smat[0][0] = from.smat[0][0];
    smat[0][1] = from.smat[0][1];
    smat[1][0] = from.smat[1][0];
    smat[1][1] = from.smat[1][1];
  }

  FlavorSpinMatrix& operator=(const FlavorSpinMatrix &from){
    smat[0][0] = from.smat[0][0];
    smat[0][1] = from.smat[0][1];
    smat[1][0] = from.smat[1][0];
    smat[1][1] = from.smat[1][1];
    return *this;
  }

  inline SpinMatrix &operator()(int f1,int f2){
    return smat[f1][f2];
  }

  inline const SpinMatrix &operator()(int f1,int f2) const{
    return smat[f1][f2];
  }

  FlavorSpinMatrix operator*(const FlavorSpinMatrix& rhs){
    FlavorSpinMatrix out;
    out.smat[0][0] = smat[0][0]*rhs.smat[0][0] + smat[0][1]*rhs.smat[1][0];
    out.smat[1][0] = smat[1][0]*rhs.smat[0][0] + smat[1][1]*rhs.smat[1][0];
    out.smat[0][1] = smat[0][0]*rhs.smat[0][1] + smat[0][1]*rhs.smat[1][1];
    out.smat[1][1] = smat[1][0]*rhs.smat[0][1] + smat[1][1]*rhs.smat[1][1];
    return out;
  }

  inline Complex Trace(){
    return smat[0][0].Tr() + smat[1][1].Tr();
  }
};

static inline Rcomplex Trace(const FlavorSpinMatrix &a, const FlavorSpinMatrix &b) {
    Rcomplex ret = 0;
    for(int s1=0;s1<4;s1++)
      for(int f1=0;f1<2;f1++)
	for(int s2=0;s2<4;s2++)
	  for(int f2=0;f2<2;f2++)
	    ret += a(f1,f2).operator()(s1,s2) * b(f2,f1).operator()(s2,s1);
    return ret;
}




CPS_END_NAMESPACE

#endif
