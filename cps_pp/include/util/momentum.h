#ifndef INCLUDED_MOMENTUM_H
#define INCLUDED_MOMENTUM_H

#include "site.h"
//---------------------------------------------------------------------------
//  math lib for the calculation of non-zero spatial momenta
//---------------------------------------------------------------------------
#include <math.h>

class ThreeMom
{
  int p[3] ; // holds the momentum
  Float pp[3] ; /* holds the minimum lattice momentum in each direction  *
		 * x--> 0   y--> 1   z--> 2                              */

  bool ZeroMom ;

  void CalcLatMom(void) ;

  void setZeroMomFlag(void) ;

 public:
  ThreeMom() ;

  ThreeMom(int *q) ;

  ThreeMom(int q0,int q1,int q2) ;

  ThreeMom(const ThreeMom& rhs) ;

  ~ThreeMom(){} ;

  Complex Fact(Site& s) ; // exp(-i p*x) 
  
  Complex Fact(Site& s, int *sx) ; // exp(-i p*(x+sx)) 

  int cmp(int i){return p[i];}
  
  void conj(){    
    p[0]=-p[0] ;
    p[1]=-p[1] ;
    p[2]=-p[2] ;
  }
} ;
#endif // !INCLUDED_MOMENTUM_H  
