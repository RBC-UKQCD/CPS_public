#include <config.h>
#ifndef INCLUDED_MOMENTUM_H
#define INCLUDED_MOMENTUM_H

#include <util/site.h>
//---------------------------------------------------------------------------
//  math lib for the calculation of non-zero spatial momenta
//---------------------------------------------------------------------------
#include <math.h>

CPS_START_NAMESPACE

class ThreeMom
{
  int p[3] ; // holds the momentum

  bool ZeroMom ;

  void CalcLatMom(void) ;

  void setZeroMomFlag(void) ;

 protected:
  Float pp[3] ; /* holds the minimum lattice momentum in each direction  *
		 * x--> 0   y--> 1   z--> 2                              */

 public:
  ThreeMom() ;

  ThreeMom(const int *q) ;

  ThreeMom(int q0,int q1,int q2) ;

  ThreeMom(const ThreeMom& rhs) ;

  ~ThreeMom(){} ;

  Complex Fact(Site& s) ; // exp(-i p*x) 

  Complex FactCos(Site& s) ; // cos(p_x x)*cos(p_y y)*cos(p_z z)

  Complex Fact(Site& s, int *sx) ; // exp(-i p*(x+sx)) 

  int cmp(int i){return p[i];}
  
  void conj(){    
    p[0]=-p[0] ;
    p[1]=-p[1] ;
    p[2]=-p[2] ;
  }
} ;


// Momentum for twisted boundary conditions, in units
// of Pi/L rather than 2*Pi/L
class ThreeMomTwist : public ThreeMom
{
 public:
  ThreeMomTwist() ;
  
  ThreeMomTwist(const int *q) ;
  
  ThreeMomTwist(int q0,int q1,int q2) ;
  
  ThreeMomTwist(const ThreeMomTwist& rhs) ;
} ;


CPS_END_NAMESPACE
#endif // !INCLUDED_MOMENTUM_H  
