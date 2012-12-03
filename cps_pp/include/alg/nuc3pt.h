//------------------------------------------------------------------
//
// nuc3pt.h
//
// Header file for the Nucl3pt  class 
//
// WARNING: I need to check the convetions of the projectors and the
//          operator definitions. There are factors of i missing
//          in several places....
//
// April 2002
//
// Kostas Orginos
//
//------------------------------------------------------------------


#ifndef INCLUDED_NUC3PT_H
#define INCLUDED_NUC3PT_H
#include <string.h>
#include <util/site.h>
#include <alg/qpropw.h>
#include "corrfunc.h"
#include "gamma.h"
#include "derivative.h"
#include "dterms.h"

CPS_START_NAMESPACE

class Nuc3pt
{

  CorrFunc u_cf ;
  CorrFunc d_cf ;

  Complex factor ; // A complex number to be used as multiplicative factor
  Float quark_mass ;

  char *source ;
  char *sink ;
  char *project ;

  int srcdesc[3] ; // Source descriptor (box size / point location etc.)
  int src_t ; // Source timeslice
  int snk_t ; // Sink timeslice
  ThreeMom snk_mom ; // The momentum at the sink

protected:
  ThreeMom mom ; // Momentum injection at the operator
  char *cname ;

public:
  CorrFunc* u_cfm[56] ;
  CorrFunc* d_cfm[56] ;

  Nuc3pt() ; 
  Nuc3pt(ThreeMom m) ; 
  Nuc3pt(Complex cc) ;   
  Nuc3pt(ThreeMom m, Complex cc) ; 
  
  virtual void InsertOp(CorrFunc& tmp,QPropW& seqQ, QPropW& Quark)=0 ;
  virtual void InsertOp(CorrFunc* tmp,QPropW& seqQ, QPropW& Quark, int Nmom, ThreeMom* mom)=0 ;
  virtual void PrintTheTag(FILE *fp)=0 ;

  void Calc3pt(QPropWSeqBar& seqQ, QPropW& Quark);
  void Calc3pt(QPropWSeqBar& seqQ, QPropW& Quark, int Nmom, ThreeMom* mom);
  void Calc3pt(QPropWMultSeqBar& seqQ, QPropW& Quark);
  void Calc3pt(QPropWMultSeqBar& seqQ, QPropW& Quark, int Nmom, ThreeMom* mom);

  void Zero()
  {
    u_cf.Zero();
    d_cf.Zero();
  }
  
  virtual ~Nuc3pt(){}
 
  void Print(FILE *fp) ;
  void Print(FILE *fp, int Nmom, ThreeMom*) ;
  
  Nuc3pt& operator+=(Nuc3pt& rhs) ;
  Nuc3pt& operator-=(Nuc3pt& rhs) ;

} ;

/*!
  Class for construncting the nucleon three point function
  of operators containing only gamma matrices

  \f[{\cal O} =\gamma_\mu\,\gamma_\nu\,\gamma_\rho\f]

  The indices can be 0 1 2 3 : (\f$\gamma_x\f$,\f$\gamma_y\f$,\f$\gamma_z\f$)
  and -5 : \f$ \gamma_5 \f$
*/
class Nuc3ptGamma : public Nuc3pt
{
  Gamma G ;
 
 public:  
  Nuc3ptGamma(Gamma op) ; 
  Nuc3ptGamma(ThreeMom m, Gamma op) ; 

  Nuc3ptGamma(Complex cc, Gamma op) ; 

  Nuc3ptGamma(ThreeMom m, Complex cc, Gamma op) ; 
  
  void InsertOp(CorrFunc& tmp,QPropW& seqQ, QPropW& Quark) ;
  void InsertOp(CorrFunc* tmp,QPropW& seqQ, QPropW& Quark, int Nmom, ThreeMom* mom) ;

  void PrintTheTag(FILE *fp){G.printTag(fp);} 

  ~Nuc3ptGamma(){} ;
} ;


/*!
  Class for construncting the nucleon three point function
  of operators containing  gamma matrices and derivatives

  \f[{\cal O} =\gamma_\mu\,\gamma_\nu\,\gamma_\rho 
      \stackrel{\displaystyle \leftrightarrow}{D}_\lambda 
      \stackrel{\displaystyle \leftrightarrow}{D}_\sigma .... \f]

  The indices can be 0 1 2 3 : (x y z t)
  and -5 : \f$ \gamma_5 \f$
*/
class Nuc3ptStru : public Nuc3pt
{
  Gamma G ;
  Derivative D ;

  char tag[6] ; // holds a characteristic string for the operator

 public:

  Nuc3ptStru(Gamma gg, Derivative dd); 

  Nuc3ptStru(Complex cc, Gamma gg, Derivative dd); 

  Nuc3ptStru(ThreeMom m, Gamma gg, Derivative dd); 

  Nuc3ptStru(ThreeMom m, Complex cc, Gamma gg, Derivative dd) ;

  void setTag(const char *tt){if(strlen(tt)<5)strcpy(tag,tt) ;}
    
  void InsertOp(CorrFunc& tmp,QPropW& seqQ, QPropW& Quark) ;
  void InsertOp(CorrFunc* tmp,QPropW& seqQ, QPropW& Quark, int Nmom, ThreeMom* mom) ;

  void PrintTheTag(FILE *fp) ;

  ~Nuc3ptStru(){} ;
} ;

/*!
  Class for construncting the nucleon three point function
  of conserved current operators containing only gamma matrices

  \f[{\cal O} =\gamma_\mu\,\gamma_\nu\,\gamma_\rho\f]

  The indices can be 0 1 2 3 : (\f$\gamma_x\f$,\f$\gamma_y\f$,\f$\gamma_z\f$)
  and -5 : \f$ \gamma_5 \f$
*/
class Nuc3ptCons : public Nuc3pt
{
  Gamma G ;
 
 public:  
  Nuc3ptCons(Gamma op) ; 
  Nuc3ptCons(ThreeMom m, Gamma op) ; 

  Nuc3ptCons(Complex cc, Gamma op) ; 

  Nuc3ptCons(ThreeMom m, Complex cc, Gamma op) ; 
  
  Nuc3ptCons(int Nmom, Gamma op) ; 
  Nuc3ptCons(int Nmom, Complex cc, Gamma op) ; 

  void InsertOp(CorrFunc& tmp,QPropW& seqQ, QPropW& Quark) ;
  void InsertOp(CorrFunc* tmp,QPropW& seqQ, QPropW& Quark, int Nmom, ThreeMom* mom) ;

  void PrintTheTag(FILE *fp){G.printTag(fp);} 

  int siteOffset(const int lcl_site[], const int lcl_sites[]) const; 

  ~Nuc3ptCons(){} ;
} ;

CPS_END_NAMESPACE

#endif //!INCLUDED_NUC3PT_H


