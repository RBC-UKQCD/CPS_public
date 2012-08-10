#ifndef __FOURMOM__CD
#define __FOURMOM__CD

#include <config.h>
#include <vector>

CPS_START_NAMESPACE

/*!
  Simple wrapper class around an array of four
  integers used to pass the momenta in units
  of 2\pi/L to alg_fourier_prop. If the boundary 
  conditions are anti-periodic then the values 
  should be shifted *up* by 1/2 when used.
*/
class FourMom
{
private:
  
  int val[4];

public:
  
  FourMom();
  FourMom(const int[]);
  FourMom(int,int,int,int);
  FourMom(const FourMom&);
  
  FourMom& operator=(const FourMom& a);
  
  inline int x() const { return val[0]; }
  inline int y() const { return val[1]; }
  inline int z() const { return val[2]; }
  inline int t() const { return val[3]; }

  const int* as_array() const { return val; }

  FourMom& operator+=(const FourMom&);
  FourMom& operator-=(const FourMom&);
  FourMom& operator*=(int);
};


/*!
  class to hold an array with elements of
  type FourMom
*/

class MomentaList
{
private:

  FourMom* vals;
  int      _size;

  char* cname;


  /*! Only pass by referance once created */
  MomentaList(MomentaList&){};
  
  /*! free memory if pointer is not null */
  void dealloc();
  
  /*! allocate (or re-allocate) array of the passed size. */
  void alloc(int);
  
public:

  MomentaList();
  MomentaList(int size);
  
  ~MomentaList() { dealloc(); }
  
  FourMom& operator[](int i) const { return vals[i]; }
 
  int   size()  const { return _size; }
  void  size(int a)   { alloc(a); }
  
};

CPS_END_NAMESPACE

#endif 
