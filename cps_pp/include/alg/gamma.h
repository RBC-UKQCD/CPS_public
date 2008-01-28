#include <stdio.h>
#include <util/error.h>
#include <util/qcdio.h>
#ifndef INCLUDED_GAMMA_H
#define INCLUDED_GAMMA_H

CPS_START_NAMESPACE

/*!------------------------------------------------------------------*
//
// gamma.h
//
// The  Gamma class. Holds the idices for a string of up to 3 gamma
// matrices.  Used in the construction of three point functions
//
// Used to represent lists of idices of operators
//
// April 2002
//
// Kostas Orginos
//
//------------------------------------------------------------------*/
class Gamma
{
  int _n ;
  int *_mu ;
  char *tag ;

  /*! Computes the operator Tag, a string to be used in 3pt IO.
    Also checks all the pointers in the class for properly being allocated.
  */
  void Tag() ;

 public:
  Gamma(                   ) ;
  Gamma(int m              ) ;
  Gamma(int m, int n       ) ;
  Gamma(int m, int n, int r) ; 

  /*! Copy constructor */
  Gamma(const Gamma& tt) ;
  
  ~Gamma() ; 

  int operator[](int i) const {return _mu[i] ;} 

  void setIndices(int *I)
  {
    for(int i=0; i < _n ; i++)
      _mu[i] = I[i] ; 
  }

  /*! prints the tag */
  void printTag(FILE *fp){Fprintf(fp,"%s",tag);}

  /*! Returns the the number of indices */
  int N(){return _n;}

} ;

CPS_END_NAMESPACE

#endif // !INCLUDED_GAMMA_H  





