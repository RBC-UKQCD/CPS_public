#include <stdio.h>
#include <util/spinflavormatrix.h>


CPS_START_NAMESPACE

SpinFlavorMatrix::SpinFlavorMatrix() {}
SpinFlavorMatrix::SpinFlavorMatrix(IFloat c) { *this = c; }
SpinFlavorMatrix::SpinFlavorMatrix(const Complex& c) { *this = c; }
SpinFlavorMatrix::SpinFlavorMatrix(const SpinFlavorMatrix& m) { *this = m; }


//------------------------------------------------------------------
SpinFlavorMatrix& SpinFlavorMatrix::operator=(IFloat c) {
  for(int i=0;i<fsize;i++) u[i] = c;
  return *this;
}

//------------------------------------------------------------------
SpinFlavorMatrix& SpinFlavorMatrix::operator=(const Complex& c) {
  for(int i=0;i<fsize/2;i++) ((Complex*)u)[i] = c;
  return *this;
}

SpinFlavorMatrix& SpinFlavorMatrix::operator=(const SpinFlavorMatrix& c){
  for(int i=0;i<fsize;i++) u[i] = c.u[i];
  return *this;
}

Complex& SpinFlavorMatrix::operator()(const int &s1, const int &f1, const int &s2, const int &f2){
  return ((Complex*)u)[f2 + 2*(s2 + 4*( f1 + 2*s1 ))];
}
const Complex& SpinFlavorMatrix::operator()(const int &s1, const int &f1, const int &s2, const int &f2) const{
  return ((Complex*)u)[f2 + 2*(s2 + 4*( f1 + 2*s1 ))];
}
Complex SpinFlavorMatrix::Trace() const{ 
  Complex ret(0.0);
  for(int s1=0;s1<4;s1++)
    for(int f1=0;f1<2;f1++)
      ret += (*this)(s1,f1,s1,f1);
  return ret;
}

CPS_END_NAMESPACE
