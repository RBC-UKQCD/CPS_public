//------------------------------------------------------------------
//
// nuc3pt_gammar.h
//
// Header file for the Nucl3pt_GammaR  class 
//
//
// September 2003
//
// Federico Berruto
//
//------------------------------------------------------------------


#ifndef INCLUDED_GAMMAR_H
#define INCLUDED_GAMMAR_H
#include "nuc3pt.h"

CPS_START_NAMESPACE

/*!
  Class for construncting the nucleon three point function
  of operators containing one gamma matrices time r

  \f[{\cal O} =\gamma_\mu\,\times r_\nu \f]

  The indices can be 0 1 2 3 : (\f$\gamma_x\f$,\f$\gamma_y\f$,\f$\gamma_z\f$)
  and -5 : \f$ \gamma_5 \f$
*/
class Nuc3ptGammaR : public Nuc3pt
{
  int dir ;
  Gamma G ;
  
 public:  
  Nuc3ptGammaR(Gamma op, int dd) ; 
  
  void InsertOp(CorrFunc& tmp,QPropW& seqQ, QPropW& Quark) ;
  
  void setDir(int d){dir=d;}

  void PrintTheTag(FILE *fp) ;
  
  ~Nuc3ptGammaR(){} ;
} ;

CPS_END_NAMESPACE

#endif //!INCLUDED_GAMMAR_H
