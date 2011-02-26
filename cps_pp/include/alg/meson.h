//------------------------------------------------------------------
//
// meson.h
//
// Header file for the Meson class 
// 
// February 2002
//
// Kostas Orginos
//
//------------------------------------------------------------------


#ifndef INCLUDED_MESON_H
#define INCLUDED_MESON_H
#include <stdio.h>
#include <string.h>
#include <util/data_types.h>
#include <alg/qpropw.h>
#include <comms/glb.h>
#include "corrfunc.h"

CPS_START_NAMESPACE

class Meson {
  CorrFunc func ;
  char src[9];
  char snk[8];
  int gamma[2] ;
  Float m[2] ;


 public:
  Meson(int mu, int nu, char *src_i) ;
  Meson(int mu, char *src_i) ;
  Meson(char *src_i) ;
  Meson(){} ;

  void setMass(Float mm){  m[0]=mm ; } 
  void setMass(Float mm, Float mmm){  m[0]=mm ; m[1] = mmm ; } 
  void setGamma()
    {  
      gamma[0]=gamma[1]=1969 ; 
      sprintf(snk,"GAM_0"); 
    } 
  void setGamma(int g)
    {  
      gamma[0]=g ; 
      gamma[1] =1969 ; 
      if(g<0)
	sprintf(snk,"GAM_5");
      else
	sprintf(snk,"GAM_%i",g+1);
    } 
  void setGamma(int g, int gg)
    {  
      gamma[0]=g ; 
      gamma[1] =gg ;
      if((g<0)&&(gg>-1))
	sprintf(snk,"GAM_5%i",gg+1);
      else if((g>-1)&&(gg<0))
	sprintf(snk,"GAM_%i5",g+1);
      else
	sprintf(snk,"GAM_%i%i",g+1,gg+1);
    } 
  void setSrc(char *src_i) ;

  void calcMeson(QPropW& q1, QPropW& q2) ;
  void calcMidPointPion(QPropW& q1, QPropW& q2) ;

  void Zero(){func.Zero();}

  void Print(FILE *fp)
  {
    Fprintf(fp,"STARTPROP\nMASSES:  %e   %e\nSOURCE: %s\nSINKS: %s\n",
	    //(float)m[0],(float)m[1],src,snk);
	    (Float)m[0],(Float)m[1],src,snk);
    func.print(fp) ;
    Fprintf(fp,"ENDPROP\n");
  } 

  Meson& operator*=(const Float& r){func *=r  ; return *this ; }
 
  ~Meson(){} ;
} ;

CPS_END_NAMESPACE

#endif //!INCLUDED_MESON_H
