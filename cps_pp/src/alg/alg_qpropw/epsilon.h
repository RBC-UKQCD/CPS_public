#ifndef INCLUDED_EPSILON_H
#define INCLUDED_EPSILON_H

// The antisymetric Epsilon class.
//------------------------------------------------------------------

class datum {
  public:
    int a_ ;
    int b_ ;
    int c_ ;
    float s_ ;
   
    datum(int a,int b, int c, float s): a_(a),b_(b),c_(c),s_(s){}
    //datum(){}
  } ;

class Epsilon
{
  
  int perm ; 
  static datum p[6];

public:
  float sign() const {return p[perm].s_ ;}
  int a() const {return p[perm].a_ ;}
  int b() const {return p[perm].b_ ;}
  int c() const {return p[perm].c_ ;}

  void Begin(){
    perm = 0 ;
  }
  void Begin(int a){
    perm = 0 ;
    while((p[perm].a_ != a) && (perm<6))
      perm++ ;
  }
  int Contracting()
    {
      return (perm<6) ;
    }
  void Next(){
    perm++ ;
  }
  void Next(int a){// loops only over the permutations with the first index = a
    perm++ ;
    while((p[perm].a_ != a) && (perm<6))
      perm++ ;
  }
  Epsilon(): perm(0){}
};
#ifdef COMPILE_THE_EPSILON
datum Epsilon::p[6] ={ 
     datum(0,1,2,+1.0) ,
     datum(1,2,0,+1.0) ,
     datum(2,0,1,+1.0) ,
     datum(0,2,1,-1.0) ,
     datum(1,0,2,-1.0) ,
     datum(2,1,0,-1.0)
};
#endif
#endif

